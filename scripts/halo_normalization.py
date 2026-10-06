"""Explicit halo-normalization diagnostics; no installed function is replaced.

The caller supplies this adapter as the covariance backend. It changes only
returned halo moments, by algebra exactly equivalent to rescaling their mass
weights. Spectra, profiles, mass limits and projection kernels remain native.
This is a sensitivity experiment, not a proposed production implementation.
"""

import numpy as np
from scipy.special import gamma


def mass_integral(interface, a):
    """Integrate CoCoA's Tinker f(nu,a) over 0 < nu < infinity analytically.

    For a term nu**p exp(-g*nu**2/2), its integral is
    Gamma((p+1)/2) / [2*(g/2)**((p+1)/2)]. The two powers in the
    fitted shape are 2*eta and 2*eta-2*phi. Read the actual installed
    amplitude through fnu(1,a), so this includes its table interpolation.
    This is the full extrapolated fit, not a finite-mass normalization.
    """
    scale = max(float(a), 0.25)
    beta = 0.589*scale**(-0.2)
    width = 0.864*scale**0.01
    phi = -0.729*scale**0.08
    eta = -0.243*scale**(-0.27)
    shape_at_one = (1+beta**(-2*phi))*np.exp(-width/2)
    amplitude = interface.fnu(1.0, float(a))/shape_at_one

    powers = np.array([2*eta, 2*eta-2*phi])
    exponents = (powers+1)/2
    terms = gamma(exponents)/(2*(width/2)**exponents)
    return float(amplitude*(terms[0]+beta**(-2*phi)*terms[1]))


class NormalizedMoments:
    """Forward the native API except for explicitly rescaled halo moments.

    mode = mass_only: f -> f/S, b unchanged, with S = integral f dnu.
    mode = mass_and_bias: f -> f/S, b -> S*b at the same redshift.
    The latter preserves every biased weight n*b, while reducing the
    unbiased weights n. It imposes both full-fit constraints to the
    accuracy of the original bias-weighted normalization.
    """

    def __init__(self, interface, mode):
        if mode not in ("mass_only", "mass_and_bias"):
            raise ValueError("choose mass_only or mass_and_bias")
        self.interface = interface
        self.native = interface.covariance
        self.mode = mode

    def __getattr__(self, name):
        return getattr(self.native, name)

    def covariance_halo_moments(self, *, a, k, lnm_edges, nquad,
                                pair_moments=True):
        """Return I11[na,nk] and moments[5,na,npair] in core units.

        Roles are I02, I12, I13(K,Q,Q), I13(K,K,Q), I04. The native
        additive completion must be recomputed when n*b changes:
        I11 = resolved + (1-resolved_at_zero)*u_min. Therefore dividing
        all resolved biased weights by S gives
        I11_new = u_min + (I11_old-u_min)/S, not I11_old/S.
        """
        single, moments = self.native.covariance_halo_moments(
            a=a, k=k, lnm_edges=lnm_edges, nquad=nquad,
            pair_moments=pair_moments,
        )
        minimum = float(np.exp(lnm_edges[0]))
        for row, scale in enumerate(a):
            total_mass = mass_integral(interface=self.interface, a=scale)
            if self.mode == "mass_only":
                concentration = self.interface.conc(minimum, float(scale))
                profile = np.empty(k.shape[1])
                for node, wave in enumerate(k[row]):
                    # The scalar profile formula is singular at k=0;
                    # the normalized profile's analytic limit is one.
                    profile[node] = 1.0
                    if wave != 0.0:
                        profile[node] = self.interface.u_nfw_c(
                            concentration, float(wave), minimum, float(scale),
                        )
                single[row] = profile+(single[row]-profile)/total_mass
                if pair_moments:
                    moments[:, row] /= total_mass
            elif pair_moments:
                # n*b is unchanged in this experiment. Only moments
                # without halo bias receive the mass normalization.
                moments[0, row] /= total_mass
                moments[4, row] /= total_mass
        return single, moments


def check_ingredients(interface, edges, cosmology):
    """Check full-fit sum rules and finite-mass moments at several redshifts.

    Numerical integrals of the public fitted functions check the analytic
    mass factor. Small direct mass sums check the adapter, including its
    additive completion, against the stated modified integrands. These are
    tests of the diagnostic transformation, not another physical model.
    """
    from scipy.integrate import simpson

    nu = np.exp(np.linspace(-90, 3.5, 8193))
    # The public fnu reader requires a < 1, so stay inside its domain.
    redshifts = [0.01, 0.1, 0.5, 1.0, 2.0, 3.0]
    records = []
    for redshift in redshifts:
        scale = 1/(1+redshift)
        mass = mass_integral(interface=interface, a=scale)
        fitted = np.array([interface.fnu(float(n), scale) for n in nu])
        bias = np.array([interface.hb1nu(float(n), scale) for n in nu])
        numerical = float(simpson(fitted*nu, x=np.log(nu)))
        biased = float(simpson(fitted*bias*nu, x=np.log(nu)))
        np.testing.assert_allclose(mass, numerical, rtol=1e-10)
        records.append(dict(
            redshift=redshift,
            native_mass=mass,
            native_biased_mass=biased,
            mass_only_biased_mass=biased/mass,
            both_mass=numerical/mass,
            both_biased_mass=biased,
            abundance_multiplier=1/mass,
            paired_bias_multiplier=mass,
        ))

    # Compare a few physical scales, including k=0 and off-diagonal
    # pairs. Use the same precomputed GSL rule as the C mass integral.
    nodes, weights = interface.covariance.covariance_integration_rule(nquad=96)
    masses = []
    measures = []
    for lower, upper in zip(edges[:-1], edges[1:]):
        half = (upper-lower)/2
        masses.extend(np.exp((upper+lower)/2+half*nodes))
        measures.extend(half*weights)
    masses = np.asarray(masses)
    measures = np.asarray(measures)
    density = 7.4775e21*cosmology["omegam"]
    volume = masses/density
    wave = np.array([0.0, 0.01, 1.0, 100.0])*2997.92458
    first, second = np.triu_indices(len(wave))
    minimum = float(np.exp(edges[0]))
    finite = []
    for redshift in (0.1, 0.5, 1.0):
        scale = 1/(1+redshift)
        factor = mass_integral(interface=interface, a=scale)
        number = np.empty(len(masses))
        bias = np.empty(len(masses))
        profile = np.empty((len(wave), len(masses)))
        for node, mass in enumerate(masses):
            peak = 1.686/np.sqrt(interface.sigma2(float(mass), scale, 1))
            slope = interface.dlognudlogm(float(mass), scale)
            number[node] = (measures[node]*interface.fnu(peak, scale)
                            *peak*slope/volume[node])
            bias[node] = interface.hb1nu(peak, scale)
            concentration = interface.conc(float(mass), scale)
            for index, value in enumerate(wave):
                profile[index, node] = 1.0
                if value != 0.0:
                    profile[index, node] = interface.u_nfw_c(
                        concentration, float(value), float(mass), scale,
                    )
        resolved_mass = float(np.sum(number*volume))
        resolved_response = float(np.sum(number*volume*bias))
        finite.append(dict(
            redshift=redshift,
            resolved_mass=resolved_mass,
            resolved_biased_mass=resolved_response,
            missing_mass=1-resolved_mass,
            missing_response=1-resolved_response,
            effective_unresolved_bias=(1-resolved_response)/(1-resolved_mass),
        ))
        concentration = interface.conc(minimum, scale)
        floor_profile = np.ones(len(wave))
        for index, value in enumerate(wave):
            if value != 0.0:
                floor_profile[index] = interface.u_nfw_c(
                    concentration, float(value), minimum, scale,
                )
        for mode in ("mass_only", "mass_and_bias"):
            backend = NormalizedMoments(interface=interface, mode=mode)
            single, moments = backend.covariance_halo_moments(
                a=np.array([scale]), k=wave[None, :],
                lnm_edges=edges, nquad=96,
            )
            new_number = number/factor
            new_bias = bias
            if mode == "mass_and_bias":
                new_bias = bias*factor
            biased_weight = new_number*new_bias*volume
            expected_single = profile@biased_weight
            expected_single += (1-sum(biased_weight))*floor_profile
            np.testing.assert_allclose(single[0], expected_single, rtol=2e-10)
            np.testing.assert_allclose(single[0, 0], 1.0, atol=2e-12)

            # Profile products match the five named moment roles above.
            left, right = profile[first], profile[second]
            products = [left*right, left*right, left*right**2,
                        left**2*right, left**2*right**2]
            for role, power in enumerate((2, 2, 3, 3, 4)):
                weight = new_number*volume**power
                if role in (1, 2, 3):
                    weight = weight*new_bias
                expected = products[role]@weight
                np.testing.assert_allclose(moments[role, 0], expected,
                                           rtol=2e-10, atol=0)
    return dict(full_fit=records, finite_mass=finite,
                adapter_direct_quadrature_checks="passed")
