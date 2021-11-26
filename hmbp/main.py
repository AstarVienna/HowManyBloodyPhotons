import numpy as np
from astropy import units as u


from scopesim.effects import ter_curves_utils as tcu
from scopesim.source import spectrum_templates as st

FILTER_DEFAULTS = tcu.FILTER_DEFAULTS


def for_flux_in_filter(filter_name, flux, instrument=None, observatory=None):
    """
    Returns the number of photons in a given filter at a given magnitude

    Filter anme, instrument, and observatory are as per the Spanish VO filter
    service name

    Parameters
    ----------
    filter_name : str
    flux : float, Quantity
        [vega mag, AB mag, Jy] If float, Vega mag is assumed
    instrument : str, optional
    observatory : str, optional

    Returns
    -------
    n_ph : expected photon count in the filter
        [ph s-1 m-2 (arcsec-2)]

    Examples
    --------
    Find the number of photons emitted by Vega in V band [ph s-1 m-2]::

        >>> from astropy import units as u
        >>> import hmbp
        >>>
        >>> hmbp.in_zero_vega_mags("V")

    Find the number of photons emitted by a Ks=20 [ABmag] point source through
    the HAWKI Ks filter::

        >>> hmbp.for_flux_in_filter("Ks", 20*u.ABmag,
                                         observatory="Paranal",
                                         instrument="HAWKI")

    See Also
    --------
    http://svo2.cab.inta-csic.es/theory/fps/

    """

    if observatory is None and instrument is None and \
            filter_name in FILTER_DEFAULTS:
        svo_str = FILTER_DEFAULTS[filter_name]
    else:
        svo_str = f"{observatory}/{instrument}.{filter_name}"

    spec_template = st.vega_spectrum
    if isinstance(flux, u.Quantity):
        if flux.unit.physical_type == "spectral flux density":             # ABmag and Jy
            spec_template = st.ab_spectrum
            flux = flux.to(u.ABmag)
        flux = flux.value

    spec = spec_template(flux)
    filt = tcu.download_svo_filter(svo_str)
    wave = filt.waveset
    dwave = 0.5 * (np.r_[[0], np.diff(wave)] + np.r_[np.diff(wave), [0]])

    flux = spec(wave)  # ph/s/cm2/AA
    flux *= filt(wave) * dwave  # ph/s
    n_ph = np.sum(flux.to(u.ph / u.s / u.m ** 2))

    return n_ph


def in_zero_vega_mags(filter_name, instrument=None, observatory=None):
    return for_flux_in_filter(filter_name, 0, instrument, observatory)


def in_zero_AB_mags(filter_name, instrument=None, observatory=None):
    return for_flux_in_filter(filter_name, 0 * u.ABmag, instrument, observatory)


def in_one_jansky(filter_name, instrument=None, observatory=None):
    return for_flux_in_filter(filter_name, 1 * u.Jy, instrument, observatory)


in_zero_vega_mags.__doc__ = for_flux_in_filter.__doc__
in_zero_AB_mags.__doc__ = for_flux_in_filter.__doc__
in_one_jansky.__doc__ = for_flux_in_filter.__doc__


def convert(from_quantity, to_unit, filter_name,
            instrument=None, observatory=None):
    """
    Converts units within a certain instrument filters

    If ``filter_name`` is not a generic name as listed in
    ``hmbp.FILTER_DEFAULTS``, the instrument and observatory must be given as
    per the name definitions in the Spanish VO filter service.

    Parameters
    ----------
    from_quantity : astropy.Quantity
        [u.mag, u.ABmag, u.Jy] The flux quantity to convert
    to_unit : astropy.Unit
        [u.mag, u.ABmag, u.Jy] The Unit to convert flux to
    filter_name: str
    instrument: str, optional
    observatory: str, optional

    Returns
    -------
    new_flux : astropy.Quantity
        FLux is units of ``to_unit``

    Examples
    --------
    Convert from Vega magnitudes to Jansky in the NACO M-prime filter::

        >>> from astropy import units as u
        >>> import hmbp
        >>>
        >>> hmbp.convert(20 * u.mag, u.Jy, filter_name="Mp",
        >>>              observatory="Paranal", instrument="NACO")


    See Also
    --------
    http://svo2.cab.inta-csic.es/theory/fps/

    """

    base_fn = {u.mag: in_zero_vega_mags,
               u.ABmag: in_zero_AB_mags,
               u.Jy: in_one_jansky}[to_unit]
    to_phs = base_fn(filter_name, instrument, observatory)
    from_phs = for_flux_in_filter(filter_name, from_quantity, instrument,
                                  observatory)

    scale_factor = (from_phs / to_phs).value
    if to_unit in [u.mag, u.ABmag]:
        scale_factor = -2.5 * np.log10(scale_factor)

    new_flux = scale_factor * to_unit

    return new_flux
