import unittest
import math
import numpy as np
import pyximport
pyximport.install()
import energy_matcher


class MyTestCase(unittest.TestCase):
    def test_dca(self):
        dca_0 = energy_matcher.get_DCA(np.array([0.0, 0.0, 0.0]), 6.13, 4.7718739682, np.array([0.0, 100.0, 0.0]), np.array([10.0, 110.0, 0.0]))
        dca_10 = energy_matcher.get_DCA(np.array([0.0, 0.0, 0.0]), 6.13, 4.7718739682, np.array([10.0, 100.0, 0.0]),np.array([20.0, 110.0, 0.0]))
        dca_20 = energy_matcher.get_DCA(np.array([10.0, 0.0, 0.0]), 6.13, 4.7718739682, np.array([10.0, 100.0, 20.0]), np.array([10.0, 120.0, 40.0]))

        self.assertAlmostEqual(dca_0, 0.0)
        self.assertAlmostEqual(dca_10, 10.0)
        self.assertAlmostEqual(dca_20, 20.0)


    def test_get_scattering_angle(self):
        theta = (180.0/math.pi)*energy_matcher.get_scattering_angle(np.array([0.0, 0.0, 0.0]), np.array([0.0, 100.0, 0.0]), np.array([10.0, 110.0, 0.0]))
        self.assertAlmostEqual(theta, 45.0)

        theta = (180.0 / math.pi) * energy_matcher.get_scattering_angle(np.array([50.0, 0.0, 50.0]),
                                                                        np.array([0.0, 50 * math.sqrt(2), 0.0]),
                                                                        np.array([50.0, 100.0 * math.sqrt(2), 50.0]))
        self.assertAlmostEqual(theta, 90.0)

        theta = (180.0 / math.pi) * energy_matcher.get_scattering_angle(np.array([0.0, 0.0, 0.0]),
                                                                        np.array([0.0, 100.0, 0.0]),
                                                                        np.array([10.0, 110.0, 0.0]))
        self.assertAlmostEqual(theta, 45.0)

        theta = (180.0 / math.pi) * energy_matcher.get_scattering_angle(np.array([0.0, 0.0, 0.0]),
                                                                        np.array([0.0, 100.0, 0.0]),
                                                                        np.array([0.0, 120.0, 10.0]))
        self.assertAlmostEqual(theta, (180.0 / math.pi) * math.atan(10.0/20.0))

        theta = (180.0 / math.pi) * energy_matcher.get_scattering_angle(np.array([0.0, 0.0, 0.0]),
                                                                        np.array([0.0, 100.0, 0.0]),
                                                                        np.array([78.0, 295.0, 52.0]))
        self.assertAlmostEqual(theta, (180.0 / math.pi) * math.atan(math.sqrt(78.0*78.0 + 52.0*52.0)/195.0))

        theta = (180.0 / math.pi) * energy_matcher.get_scattering_angle(np.array([0.0, 0.0, 0.0]),
                                                                        np.array([0.0, 150.0, 0.0]),
                                                                        np.array([0.0, 100.0, 50.0]))
        self.assertAlmostEqual(theta, 135.0)

        theta = (180.0 / math.pi) * energy_matcher.get_scattering_angle(np.array([0.0, 0.0, 0.0]),
                                                                        np.array([0.0, 150.0, 0.0]),
                                                                        np.array([-50.0, 100.0, 50.0]))
        self.assertAlmostEqual(theta, 180.0 - (180.0 / math.pi) * math.atan(math.sqrt(50*50+50*50)/50.0))


    def test_get_scattering_angle_measured(self):
        theta = (180.0 / math.pi) * energy_matcher.get_scattering_angle_measured(0.4696597647, 1.17 - 0.4696597647)
        self.assertAlmostEqual(theta, 45.0)

        theta = (180.0 / math.pi) * energy_matcher.get_scattering_angle_measured(0.1358118506, 1.33 - 0.1358118506)
        self.assertAlmostEqual(theta, 17.0)

        theta = (180.0 / math.pi) * energy_matcher.get_scattering_angle_measured(1.2, 1.33 - 1.2)
        self.assertNotEqual(theta, theta)


    def test_get_E1_known(self):
        e1 = energy_matcher.get_E1_known(np.array([0.0, 0.0, 0.0]), 6.13,
                                         np.array([0.0, 100.0, 0.0]),
                                         np.array([10.0, 110.0, 0.0]))
        self.assertAlmostEqual(e1, 4.7718739682)

        e1 = energy_matcher.get_E1_known(np.array([50.0, 0.0, 50.0]), 1.33,
                                         np.array([0.0, 50 * math.sqrt(2), 0.0]),
                                         np.array([50.0, 100.0 * math.sqrt(2), 50.0]))
        self.assertAlmostEqual(e1, 0.9608365019)


    def test_get_KN_xs(self):

        '''
        Using Compton scattering of 662 keV gamma rays proposed by Klein-Nishina formula by Hossain et al
        as a reference. Also, Klein-Nishina Formula for Compton Effect from Wolfram, a CDF workbook.
        To compare the two sources, you need to input the energy into Wolfram as E/E_electron. So for
        0.662 MeV use 1.29.
        '''

        # x117 = np.array([energy_matcher.get_KN_xs(1.17, float(i)*math.pi/180.0) for i in range(0,360)])
        # x133 = np.array([energy_matcher.get_KN_xs(1.33, float(i)*math.pi/180.0) for i in range(0,360)])
        # for pair in zip(x117, x133):
        #     print "%.3e, %.3e, %.3e" % (pair[0], pair[1], pair[0] - pair[1])

        x117 = np.array([energy_matcher.get_KN_xs(1.17, energy_matcher.get_scattering_angle_measured(float(i)/100.0, 1.17 - float(i)/100.0)) for i in range(117)])
        x133 = np.array([energy_matcher.get_KN_xs(1.33, energy_matcher.get_scattering_angle_measured(float(i)/100.0, 1.33 - float(i)/100.0)) for i in range(117)])
        # for pair in zip(range(150), x117, x133):
        #     print "%d, %.3e, %.3e, %.3e" % (pair[0], pair[1], pair[2], pair[1] - pair[2])

        xs = energy_matcher.get_KN_xs(1.17, 0.175)
        print "1.17", xs, 0.920 * 7.94E-26 * 1.0E24
        self.assertAlmostEqual(xs, 0.920 * 7.94E-26 * 1.0E24, places=2)

        xs = energy_matcher.get_KN_xs(1.33, 0.698)
        self.assertAlmostEqual(xs, 0.027873, places=2)


    def test_probability_E1(self):

        xs117 = energy_matcher.probability_E1(1.17, 0.1, 0.01, rectangles_per_radian = 10)
        xs133 = energy_matcher.probability_E1(1.33, 0.1, 0.01, rectangles_per_radian = 10)
        self.assertGreater(xs117, xs133, '%.3e is not > than %.3e . . .' % (xs117, xs133))

        xs117 = energy_matcher.probability_E1(1.17, 0.2, 0.01, rectangles_per_radian = 1000)
        xs133 = energy_matcher.probability_E1(1.33, 0.2, 0.01, rectangles_per_radian = 1000)
        self.assertGreater(xs117, xs133, '%.3e is not > than %.3e . . .' % (xs117, xs133))

    # def test_get_energy_label(self):
    #
    #     e1_array =         np.array([0.98, 0.900, 1.13, 1.05, 0.24, 0.95, 0.12])
    #     e2_array =         np.array([0.19, 0.269, 0.20, 0.28, 0.54, 0.20, 1.12])
    #     correct_energies = np.array([1.33, 1.170, 0.00, 1.33, 1.17, 1.17, 1.33])
    #     # e1_array = np.array([0.90])
    #     # e2_array = np.array([0.269])
    #     # correct_energies = np.array([1.17])
    #
    #     energy_spectrum = np.array([1.17, 1.33])
    #     relative_cross_sections = np.array([0.5, 0.5])
    #     selected_energies = energy_matcher.get_energy_labels(e1_array, e2_array, energy_spectrum, relative_cross_sections)
    #
    #     for i in range(len(e1_array)):
    #         print "Testing: ", i, e1_array[i], e2_array[i], correct_energies[i], selected_energies[i]
    #         self.assertAlmostEqual(correct_energies[i], selected_energies[i], places=2)


if __name__ == '__main__':
    unittest.main()
