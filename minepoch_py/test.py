#!/usr/bin/env python3

def run_tests():
    import numpy as np
    from minepoch_py.energy import energy_check
    from minepoch_py.two_stream import two_stream_analysis

    energy_error = energy_check(check_error=True, plot=False)

    rate, err = two_stream_analysis(check_growth=True, plot=False)

    # Allow 5% error on growth rate
    rate_tolerance = 0.05
    # Analytic rate
    rate_analytic = 0.7
    rate_error = np.abs(rate - rate_analytic) / rate_analytic

    # Return 1 or 0, depending on error > tolerance
    if (energy_error > 0 or rate_error > rate_tolerance):
        if (energy_error > 0):
            print('Energy conservation check failed!')
        if (rate_error > rate_tolerance):
            print('Two stream growth check failed!')
            print('Observed growth rate = %.4E (expected %.2F)'
                  % (rate, rate_analytic))
        return 1
    else:
        return 0


if __name__ == "__main__":
    # Run error analysis
    run_tests()
