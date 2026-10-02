import unittest
import numpy as np
import pandas as pd
from orphans.correlations import (
    minimal_magnitude,
    peak_time,
    duration_between_first_and_peak,
    rate,
    color,
    heatmap
)

def make_mock_lc(
    n_obs: int = 10,
    filters: list = ["b", "c", "g", "orange", "r", "m"],
    seed: int | None = None,
    config: dict | None = None,
):
    """Create a mock light curve for testing purposes."""
    rng = np.random.default_rng(seed)
    time = np.sort(rng.uniform(0, 100, size=n_obs))
    filt = rng.choice(filters, size=n_obs)

    # Create realistic magnitude data
    true_mag = 20 - 5 * np.exp(-0.05 * time) + 0.5 * np.sin(0.2 * time)
    mags_lim = true_mag + rng.uniform(1.5, 3.0, size=n_obs)

    # make ~30% observable
    mask = rng.random(n_obs) < 0.3
    mags_lim[mask] -= rng.uniform(0.5, 1.5, size=mask.sum())

    if config is None:
        config = {
            "E0": 1e53,
            "n0": 1e-2,
            "z": 0.5,
            "thetaObs": rng.uniform(0, np.pi / 2),
            "thetaCore": rng.uniform(0, np.pi / 4),
            "thetaWing": rng.uniform(0, np.pi / 3),
        }

    return {
        "time": time,
        "mags": true_mag,
        "mags_lim": mags_lim,
        "filt": filt,
        "config": config,
    }

class TestCorrelationsFixed(unittest.TestCase):

    def setUp(self):
        """Set up test fixtures before each test method."""
        self.mock_lc_simu = make_mock_lc(seed=42)
        self.mock_lc_nonobs = {
            "mags": np.array([100.0]),
            "mags_lim": np.array([99.0]),
            "time": np.array([0.0]),
            "filt": np.array(["g"]),
            "config": {
                "E0": np.nan, "n0": np.nan, "z": np.nan,
                "thetaObs": np.nan, "thetaCore": np.nan, "thetaWing": np.nan,
            },
        }

    def test_minimal_magnitude_basic(self):
        """Test basic functionality of minimal_magnitude function."""
        lc_open = [self.mock_lc_simu]
        result = minimal_magnitude(lc_open)
        
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        # Should return a valid magnitude or np.nan if no observable points
        self.assertIsInstance(result[0], (float, np.float64))
        
    def test_peak_time_basic(self):
        """Test basic functionality of peak_time function."""
        lc_open = [self.mock_lc_simu]
        result = peak_time(lc_open, data="simu")
        
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        # Should return a valid time or np.nan if no observable points
        self.assertIsInstance(result[0], (float, np.float64, np.nan))
        
    def test_duration_between_first_and_peak_basic(self):
        """Test basic functionality of duration_between_first_and_peak function."""
        lc_open = [self.mock_lc_simu]
        dt, t0 = duration_between_first_and_peak(lc_open)
        
        # Should return two structures - dt (list of lists) and t0 (list)
        self.assertIsInstance(dt, list)
        self.assertIsInstance(t0, list)
        
        # Check that dt is structured as expected (6 filters, each with a list)
        self.assertEqual(len(dt), 6)
        for sublist in dt:
            self.assertIsInstance(sublist, list)
            
    def test_rate_basic(self):
        """Test basic functionality of rate function."""
        lc_open = [self.mock_lc_simu]
        rate_inc, rate_dec_1, rate_dec_3 = rate(lc_open, data="simu")
        
        # Should return three lists, each with 6 elements (one per filter)
        self.assertIsInstance(rate_inc, list)
        self.assertIsInstance(rate_dec_1, list)
        self.assertIsInstance(rate_dec_3, list)
        self.assertEqual(len(rate_inc), 6)
        self.assertEqual(len(rate_dec_1), 6)
        self.assertEqual(len(rate_dec_3), 6)
        
    def test_color_basic(self):
        """Test basic functionality of color function."""
        lc_open = [self.mock_lc_simu]
        result = color(lc_open, data="simu")
        
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        # Should return a valid color value or np.nan
        self.assertIsInstance(result[0], (float, np.float64, np.nan))
        
    def test_edge_cases(self):
        """Test edge cases for all functions."""
        # Test with empty or invalid input
        result = minimal_magnitude([0])
        self.assertEqual(len(result), 1)
        self.assertTrue(np.isnan(result[0]))
        
        result = peak_time([0], data="simu")
        self.assertEqual(len(result), 1)
        self.assertTrue(np.isnan(result[0]))
        
        # Test with a non-observable light curve
        result = minimal_magnitude([self.mock_lc_nonobs])
        self.assertEqual(len(result), 1)
        self.assertTrue(np.isnan(result[0]))
        
        result = peak_time([self.mock_lc_nonobs], data="simu")
        self.assertEqual(len(result), 1)
        self.assertTrue(np.isnan(result[0]))

    def test_function_integration(self):
        """Test that functions can work together without crashing."""
        # Test that we can call multiple functions in sequence
        lc_open = [self.mock_lc_simu]
        
        # These should all complete without raising exceptions
        mag_min = minimal_magnitude(lc_open)
        time_max = peak_time(lc_open, data="simu")
        dt, t0 = duration_between_first_and_peak(lc_open)
        rate_inc, rate_dec_1, rate_dec_3 = rate(lc_open, data="simu")
        color_val = color(lc_open, data="simu")
        
        # Check all return types are correct
        self.assertIsInstance(mag_min, list)
        self.assertIsInstance(time_max, list)
        self.assertIsInstance(dt, list)
        self.assertIsInstance(t0, list)
        self.assertIsInstance(rate_inc, list)
        self.assertIsInstance(rate_dec_1, list)
        self.assertIsInstance(rate_dec_3, list)
        self.assertIsInstance(color_val, list)

if __name__ == '__main__':
    unittest.main()