"""Change only mass-node placement in an isolated covariance build.

The upper integral retains its requested production quadrature. Below
1e4 Msun/h, a diagnostic environment variable chooses a cheaper rule.
All halo weights, profiles, SIMD sums and additive completion stay native.
This helper never edits the installed halo_cov.c.
"""


def private_source(source):
    """Return the native source with a separate low-mass quadrature rule."""
    replacements = {
        "  const int nmass = npanel*nquad;": """
  // Diagnostic only: keep the ordinary rule above 1e4, and count the
  // samples assigned to each lower panel separately. A single lower
  // panel tests 32 nodes across the whole tail; several panels test
  // the successive partial sums with the same rule on each interval.
  const char* tail_option = getenv("COCOA_DIAGNOSTIC_TAIL_NQUAD");
  const int tail_nquad = tail_option == NULL ? nquad : atoi(tail_option);
  const double split_mass = log(1.e4);

  if (tail_nquad != 32
      && tail_nquad != 64
      && tail_nquad != 96
      && tail_nquad != 128
      && tail_nquad != 256) {
    log_fatal("diagnostic low-mass rule is unsupported");
    exit(1);
  }

  int nmass = 0; // total number of mass samples in both pieces

  for (int panel=0; panel<npanel; panel++) {
    nmass += lnm_edges[panel] < split_mass ? tail_nquad : nquad;
  }
""",
        "  gsl_integration_glfixed_table* rule = malloc_gslint_glfixed(nquad);":
        """  gsl_integration_glfixed_table* rule = malloc_gslint_glfixed(nquad);
  gsl_integration_glfixed_table* tail_rule = malloc_gslint_glfixed(tail_nquad);
  int offset = 0; // first sample of the current panel in the shared arrays
""",
        """    for (int node=0; node<nquad; node++) {
      const int index = panel*nquad+node;""":
        """    const int low_mass = lnm_edges[panel] < split_mass;
    const int panel_nodes = low_mass ? tail_nquad : nquad;
    gsl_integration_glfixed_table* panel_rule = low_mass ? tail_rule : rule;

    for (int node=0; node<panel_nodes; node++) {
      const int index = offset+node;""",
        "node, &lnm, &measure, rule);": "node, &lnm, &measure, panel_rule);",
        """      mass[2][index] = measure/mass[1][index];
    }
  }

  gsl_integration_glfixed_table_free(rule);""":
        """      mass[2][index] = measure/mass[1][index];
    }
    offset += panel_nodes;
  }

  gsl_integration_glfixed_table_free(rule);
  gsl_integration_glfixed_table_free(tail_rule);""",
    }
    for before, after in replacements.items():
        if source.count(before) != 1:
            raise ValueError("the native mass-rule construction has changed")
        source = source.replace(before, after)
    return source
