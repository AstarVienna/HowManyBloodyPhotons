import pytest
from pytest import approx

import astropy.units as u
import hmbp


class TestForMagnitudeInFilter:
    def test_returns_same_photon_count_for_V_for_AB_Vega_mags(self):
        vega_phs = hmbp.for_flux_in_filter("V", 0 * u.mag)
        ab_phs = hmbp.for_flux_in_filter("V", 0 * u.ABmag)

        assert vega_phs.value == approx(ab_phs.value, rel=0.02)

    def test_1_85mag_difference_in_photon_count_for_Ks_for_AB_Vega_mags(self):
        vega_phs = hmbp.for_flux_in_filter("Ks", 30 * u.mag)
        ab_phs = hmbp.for_flux_in_filter("Ks", 30 * u.ABmag)
        scale_factor = 10 ** (-0.4 * 1.85)

        assert vega_phs.value == approx(ab_phs.value * scale_factor, rel=0.02)

    def test_ABmag_and_Jansky_are_compatible(self):
        ab_phs = hmbp.for_flux_in_filter("J", 0 * u.ABmag)
        jy_phs = hmbp.for_flux_in_filter("J", 3631 * u.Jy)

        assert ab_phs.value == approx(jy_phs.value, rel=0.02)


class TestInZeroVegaMags:
    @pytest.mark.parametrize("filter_name, ph_exp",
                             [("J", 2.56e9), ("H", 2.76e9), ("Ks", 1.27e9)])
    def test_returns_expected_photon_count_for_hawki(self, filter_name, ph_exp):
        vega_phs = hmbp.in_zero_vega_mags(filter_name, "HAWKI", "Paranal")

        assert vega_phs.value == approx(ph_exp, rel=0.02)


class TestInOneJansky:
    @pytest.mark.parametrize("filter_name, vega_ab",
                             [("I", 0.45), ("H", 1.39)])
    def test_returns_scaled_counts_to_vega_mags(self, filter_name, vega_ab):
        scale_factor = 10 ** (-0.4 * vega_ab)
        vega_phs = hmbp.in_zero_vega_mags(filter_name) / 3631 / scale_factor
        jy_phs = hmbp.in_one_jansky(filter_name)

        assert vega_phs.value == approx(jy_phs.value, rel=0.02)


class TestConvert:
    def test_0_vega_mag_is_3631_jansky_in_V(self):
        n_jy = hmbp.convert(0*u.mag, u.Jy, "V")

        assert n_jy.value == approx(3631, rel=0.02)

    def test_0_vega_mag_is_0_91_AB_mag_in_J(self):
        ab_mag = hmbp.convert(0 * u.mag, u.ABmag, "J")

        assert ab_mag.value == approx(0.91, rel=0.02)

    @pytest.mark.parametrize("filter_name, diff",
                             [("J", 0.93), ("H", 1.34), ("Ks", 1.85)])
    def test_returns_correct_difference_from_vega_to_AB(self, filter_name, diff):
        """HAWKI Vega-AB colours differ slightly from standard Bessel filters"""
        mag_diff = hmbp.convert(from_quantity=0*u.mag, to_unit=u.ABmag,
                                filter_name=filter_name, instrument="HAWKI",
                                observatory="Paranal")

        assert mag_diff.value == approx(diff, rel=0.02)
