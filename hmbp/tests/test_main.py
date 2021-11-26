import pytest
from pytest import approx

import astropy.units as u
import hmbp


class TestInVegaSpectrum:
    @pytest.mark.parametrize("filter_name, ph_exp",
                             [("J", 2.56e9), ("H", 2.76e9), ("Ks", 1.27e9)])
    def test_returns_expected_photon_count_for_hawki(self, filter_name, ph_exp):
        vega_photons = hmbp.in_vega_spectrum(filter_name, "HAWKI", "Paranal")
        assert vega_photons == approx(ph_exp, rel=0.02)


class TestForMagnitudeInFilter:
    def test_returns_same_photon_count_for_V_for_AB_Vega_mags(self):
        vega_phs = hmbp.for_magnitude_in_filter("V", 0*u.mag)
        ab_phs = hmbp.for_magnitude_in_filter("V", 0*u.ABmag)

        assert vega_phs == approx(ab_phs, rel=0.02)

    def test_1_85mag_difference_in_photon_count_for_Ks_for_AB_Vega_mags(self):
        vega_phs = hmbp.for_magnitude_in_filter("Ks", 30 * u.mag)
        ab_phs = hmbp.for_magnitude_in_filter("Ks", 30 * u.ABmag)
        scale_factor = 10 ** (-0.4 * 1.85)

        assert vega_phs == approx(ab_phs * scale_factor, rel=0.02)
