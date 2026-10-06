"""Read selected intermediate arrays without changing OneCov's methods."""

import sys


def capture_return(function, names, call):
    """Call an unchanged Python function and save named locals at its return.

    Some OneCov ingredients have no public export. A temporary profiler
    reads them when the named function returns; it does not replace the
    function or change its arguments, arithmetic or output. Profiled times
    must not be presented as an uninstrumented performance measurement.
    """
    captured = {}
    previous = sys.getprofile()

    def observe(frame, event, result):
        if event == "return" and frame.f_code is function.__code__:
            for name in names:
                value = frame.f_locals[name]
                captured[name] = value.copy() if hasattr(value, "copy") else value

    sys.setprofile(observe)
    try:
        result = call()
    finally:
        sys.setprofile(previous)
    if set(captured) != set(names):
        raise RuntimeError("native function did not expose the expected inputs")
    return result, captured
