""" Unit test for the tools module
"""
import unittest
import numpy as np

from orphans.tools import flux_to_mag


class ToolsTestCase(unittest.TestCase):
    """ Test class for the ghosts.tools module"""
    def test_flux_to_mag(self):
        """ Verify we get the correct magnitude"""
        mag = flux_to_mag(3631000)
        self.assertAlmostEqual(mag, 0, delta=0.001)  # add assertion here


if __name__ == '__main__':
    unittest.main()
