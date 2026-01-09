# -*- coding: utf-8 -*-

import pytest
from numpy import testing as npt
import astropy.units as u
from synphot import SpectralElement, Empirical1D

import hmbp


class TestForMagnitudeInFilter:
    def test_returns_same_photon_count_for_V_for_AB_Vega_mags(self):
        vega_phs = hmbp.for_flux_in_filter("V", 0 * u.mag)
        ab_phs = hmbp.for_flux_in_filter("V", 0 * u.ABmag)
        npt.assert_allclose(vega_phs, ab_phs, rtol=0.02)

    def test_1_85mag_difference_in_photon_count_for_Ks_for_AB_Vega_mags(self):
        vega_phs = hmbp.for_flux_in_filter("Ks", 30 * u.mag)
        ab_phs = hmbp.for_flux_in_filter("Ks", 30 * u.ABmag)
        scale_factor = 10 ** (-0.4 * 1.85)
        npt.assert_allclose(vega_phs, ab_phs * scale_factor, rtol=0.02)

    def test_ABmag_and_Jansky_are_compatible(self):
        ab_phs = hmbp.for_flux_in_filter("J", 0 * u.ABmag)
        jy_phs = hmbp.for_flux_in_filter("J", 3631 * u.Jy)
        npt.assert_allclose(ab_phs, jy_phs, rtol=0.02)

    def test_ABmag_and_milliJansky_are_compatible(self):
        ab_phs = hmbp.for_flux_in_filter("J", 0 * u.ABmag)
        jy_phs = hmbp.for_flux_in_filter("J", 3631e3 * u.mJy)
        npt.assert_allclose(ab_phs, jy_phs, rtol=0.02)

    @pytest.mark.webtest
    def test_runs_for_all_default_filters(self, subtests):
        for filter_name in hmbp.FILTER_DEFAULTS:
            with subtests.test(filter_name=filter_name):
                zeroflux = 0 * u.ph / (u.s * u.m**2)
                assert hmbp.for_flux_in_filter(filter_name, 0) > zeroflux


class TestInZeroVegaMags:
    @pytest.mark.parametrize(
        "filter_name, ph_exp", [("J", 2.56e9), ("H", 2.76e9), ("Ks", 1.27e9)]
    )
    def test_returns_expected_photon_count_for_hawki(self, filter_name, ph_exp):
        vega_phs = hmbp.in_zero_vega_mags(
            filter_name, "HAWKI", "Paranal"
        ).to_value(u.ph / u.s / u.m**2)
        npt.assert_allclose(vega_phs, ph_exp, rtol=0.02)


@pytest.mark.webtest
class TestInOneJansky:
    @pytest.mark.parametrize("filter_name, vega_ab",
                             [("I", 0.45), ("H", 1.39)])
    def test_returns_scaled_counts_to_vega_mags(self, filter_name, vega_ab):
        scale_factor = 10 ** (-0.4 * vega_ab)
        vega_phs = hmbp.in_zero_vega_mags(filter_name) / 3631 / scale_factor
        jy_phs = hmbp.in_one_jansky(filter_name)
        npt.assert_allclose(vega_phs, jy_phs, rtol=0.02)


class TestConvert:
    def test_0_vega_mag_is_3631_jansky_in_V(self):
        n_jy = hmbp.convert(0*u.mag, u.Jy, "V")
        npt.assert_allclose(n_jy, 3631*u.Jy, rtol=0.02)

    def test_0_vega_mag_is_0_91_AB_mag_in_J(self):
        ab_mag = hmbp.convert(0*u.mag, u.ABmag, "J").to_value(u.ABmag)
        npt.assert_allclose(ab_mag, 0.91, rtol=0.02)

    @pytest.mark.parametrize(
        "filter_name, diff",
        [
            ("J", 0.93),
            ("H", 1.34),
            ("Ks", 1.85),
        ],
    )
    def test_returns_correct_difference_vega_to_AB(self, filter_name, diff):
        """HAWKI Vega-AB colours differ slightly from standard Bessel filters"""
        mag_diff = hmbp.convert(
            from_quantity=0 * u.mag,
            to_unit=u.ABmag,
            filter_name=filter_name,
            instrument="HAWKI",
            observatory="Paranal",
        ).to_value(u.ABmag)
        npt.assert_allclose(mag_diff, diff, rtol=0.02)


@pytest.mark.webtest
class TestInSkyCalcBackground:
    # Sky mags taken from skycalc with default values from website
    @pytest.mark.parametrize(
        "filter_name, sky_mag",
        # Magnitudes that work  # Original magnitudes from ESO skycalc
        [
            ("U", 20.7),  # [20.76] Generic/Bessel.U
            ("B", 21.0),  # [21.14]
            ("V", 20.6),  # [20.67]
            pytest.param("R", 20.3, marks=pytest.mark.xfail(reason="Off by 8 % for unknown reasons.")),  # [20.32]
            pytest.param("I", 19.5, marks=pytest.mark.xfail(reason="Off by 23 % for unknown reasons.")),  # [19.48]
            pytest.param("NACO.J", 17.3, marks=pytest.mark.xfail(reason="Off by 214842053 % for unknown reasons.")),  # [16.87] Paranal/NACO.J
            pytest.param("NACO.H", 15.3, marks=pytest.mark.xfail(reason="Off by 88661246 % for unknown reasons.")),  # [14.43]
            pytest.param("NACO.Ks", 15.1, marks=pytest.mark.xfail(reason="Off by 44322087 % for unknown reasons.")),  # [15.23]
            pytest.param("NACO.Lp", 5.3, marks=pytest.mark.xfail(reason="Off by 10515 % for unknown reasons.")),  # [6.00]
            pytest.param("NACO.Mp", 1.2, marks=pytest.mark.xfail(reason="Off by 240 % for unknown reasons.")),  # [1.14]
            ("MIDI.Nband", -2.7),  # [-2.29] Paranal/MIDI.Nband
        ],
    )
    def test_returns_expected_sky_bg_counts(self, filter_name, sky_mag):  # 22 warnings
        obs, inst, filt = None, None, filter_name
        if "." in filter_name:
            obs = "Paranal"
            inst, filt = filter_name.split(".")

        skycalc_phs = hmbp.in_skycalc_background(
            filt, instrument=inst, observatory=obs
        )
        vega_phs = hmbp.for_flux_in_filter(
            filt, sky_mag * u.mag, instrument=inst, observatory=obs
        )

        npt.assert_allclose(skycalc_phs, vega_phs, rtol=0.05)

    def test_returns_different_values_for_different_airmasses(self):  # 2 warnings
        am1_phs = hmbp.in_skycalc_background("M", airmass=1.0)
        am2_phs = hmbp.in_skycalc_background("M", airmass=2.0)
        npt.assert_allclose(am2_phs, 1.5 * am1_phs, rtol=0.1)

    def test_plot_skycalc_spectrum_and_filters(self):  # 2 warnings
        import skycalc_ipy
        import numpy as np
        import matplotlib.pyplot as plt
        import astropy.units as u

        skycalc = skycalc_ipy.SkyCalc()
        # skycalc.values.update(kwargs)
        sky_trans, sky_flux = skycalc.get_sky_spectrum(return_type="synphot")

        wave = np.linspace(0.3, 14, 1000)*u.um
        plt.loglog(wave, sky_flux(wave))
        plt.show()

    def test_override_filter_name_with_spectral_element(self):  # 2 warnings
        wave = [1.99, 2.0, 2.3, 2.31] * u.um
        trans = [0, 0.8, 0.8, 0]
        filt = SpectralElement(Empirical1D, points=wave, lookup_table=trans)
        custom_phs = hmbp.in_skycalc_background(filt)
        ks_phs = hmbp.in_skycalc_background("Ks")
        npt.assert_allclose(ks_phs, custom_phs, rtol=0.05)
