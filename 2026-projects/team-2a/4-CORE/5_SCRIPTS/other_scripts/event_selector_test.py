import unittest
import numpy

import event_selector


class MyTestCase(unittest.TestCase):
    def test_known_energies(self):

        energies = numpy.array([1.17, 1.33])
        E1 = numpy.array([1.0, 1.33,  0.96, 1.11, 1.15, 1.12, 0.96])
        E2 = numpy.array([0.17, 0.17, 0.17, 0.17, 0.20, 0.22, 0.23])

        Ek_truth = numpy.array([1.33, -100.0, 1.17, 1.33, -100.0, -100.0, 1.17])

        Ek = event_selector.get_known_energies(E1, E2, energies)

        for ek, et in zip(Ek, Ek_truth):
            self.assertEqual(ek, et)

    def test_get_dca(self):

        dca_0 = energy_matcher.get_DCA(np.array([0.0, 0.0, 0.0]), 6.13, 4.7718739682, np.array([0.0, 100.0, 0.0]), np.array([10.0, 110.0, 0.0]))
        dca_10 = energy_matcher.get_DCA(np.array([0.0, 0.0, 0.0]), 6.13, 4.7718739682, np.array([10.0, 100.0, 0.0]),np.array([20.0, 110.0, 0.0]))
        dca_20 = energy_matcher.get_DCA(np.array([10.0, 0.0, 0.0]), 6.13, 4.7718739682, np.array([10.0, 100.0, 20.0]), np.array([10.0, 120.0, 40.0]))

        self.assertAlmostEqual(dca_0, 0.0)
        self.assertAlmostEqual(dca_10, 10.0)
        self.assertAlmostEqual(dca_20, 20.0)


        energies = numpy.array([1.17, 1.33])
        E1 = numpy.array([1.0, 1.33,  0.96, 1.11, 1.15, 1.12, 0.96])
        E2 = numpy.array([0.17, 0.17, 0.17, 0.17, 0.20, 0.22, 0.23])

        Ek_truth = numpy.array([1.33, -100.0, 1.17, 1.33, -100.0, -100.0, 1.17])

        Ek = event_selector.get_known_energies(E1, E2, energies)

        for ek, et in zip(Ek, Ek_truth):
            self.assertEqual(ek, et)

if __name__ == '__main__':
    unittest.main()
