import unittest
import numpy as np
import numpy.testing as npt

from transformations import cartesian_to_kepler, kepler_to_cartesian, eccentric_anomaly_from_mean_anomaly


class TestTransformations(unittest.TestCase):
    def setUp(self) -> None:
        self.example_data = [
            {
                "cartesian": ([-2700816.14, -3314092.80, 5266346.42], [5168.606550, -5597.546618, -868.878445]),
                "kepler": (6787746.891, 0.000731104,
                           *np.radians((51.68714486, 127.5486706, 74.21987137, 24.10027677, 24.08317766, 24.06608426)))
            },
            {
                "cartesian": ([3126974.99, -6374445.74, 28673.59], [-254.91197, -83.30107, 7485.70674]),
                "kepler": (7096137.00, 0.0011219,
                           *np.radians((92.0316, 296.1384, 120.6878, 239.5437, 239.5991, 239.6546)))
            }
        ]

    def test_cartesian_to_kepler(self):
        for example_data in self.example_data:
            kepler = cartesian_to_kepler(*example_data['cartesian'])
            npt.assert_allclose(kepler, example_data['kepler'], rtol=1e-6)

    def test_kepler_to_cartesian__true_anomaly(self):
        # Test by passing only the true anomaly.
        for example_data in self.example_data:
            cartesian = kepler_to_cartesian(*(example_data['kepler'][:6]))
            npt.assert_allclose(cartesian, example_data['cartesian'], rtol=1e-2)

    def test_kepler_to_cartesian__eccentric_anomaly(self):
        # Test by passing only the eccentric anomaly.
        for example_data in self.example_data:
            cartesian = kepler_to_cartesian(*(example_data['kepler'][:5]), e_ano=example_data['kepler'][6])
            npt.assert_allclose(cartesian, example_data['cartesian'], rtol=1e-2)

    def test_kepler_to_cartesian__mean_anomaly(self):
        # Test by passing only the mean anomaly
        for example_data in self.example_data:
            cartesian = kepler_to_cartesian(*(example_data['kepler'][:5]), m=example_data['kepler'][7])
            npt.assert_allclose(cartesian, example_data['cartesian'], rtol=1e-2)
            
    def test_eccentric_anomaly_from_mean_anomaly(self):
        example_data = self.example_data[0]
        e = example_data["kepler"][1]
        m = example_data["kepler"][7]
        e_ano = eccentric_anomaly_from_mean_anomaly(e, m)
        npt.assert_allclose(e_ano, example_data['kepler'][6], rtol=1e-2)
        
    def test_true_anomaly_from_eccentric_anomaly(self):
        example_data = self.example_data[0]
        e = example_data["kepler"][1]
        e_ano = example_data["kepler"][6]
        theta = eccentric_anomaly_from_mean_anomaly(e, e_ano)
        npt.assert_allclose(theta, example_data['kepler'][5], rtol=1e-2)
        
