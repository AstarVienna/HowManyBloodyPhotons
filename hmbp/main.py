import numpy as np
from astropy import units as u


from scopesim.effects import ter_curves_utils as tcu
from scopesim.source import spectrum_templates as st

FILTER_DEFAULTS = tcu.FILTER_DEFAULTS


def for_magnitude_in_filter(filter_name, magnitude,
                            instrument=None, observatory=None):
    """
    Returns the number of photons in a given filter at a given magnitude

    Parameters
    ----------
    filter_name : str
    magnitude : float, Quantity
        [vega mag or AB mag]
    instrument : str, optional
    observatory : str, optional

    Returns
    -------
    n_ph : expected photon count in the filter
        [ph s-1 m-2 (arcsec-2)]

    Examples
    --------
    ::
        >>> from astropy import units as u
        >>> import hmbp

    Find the number of photons emitted by Vega in V band [ph s-1 m-2]::

        >>> hmbp.in_vega_spectrum("V")

    Find the number of photons emitted by a Ks=20 [ABmag] point source through
    the HAWKI Ks filter::

        >>> hmbp.for_magnitude_in_filter("Ks", 20*u.ABmag,
                                         observatory="Paranal",
                                         instrument="HAWKI")

    """

    if observatory is None and instrument is None and \
            filter_name in FILTER_DEFAULTS:
        svo_str = FILTER_DEFAULTS[filter_name]
    else:
        svo_str = f"{observatory}/{instrument}.{filter_name}"

    spec_template = st.vega_spectrum
    if isinstance(magnitude, u.Quantity):
        if magnitude.unit == u.ABmag:
            spec_template = st.ab_spectrum
        magnitude = magnitude.value

    spec = spec_template(magnitude)
    filt = tcu.download_svo_filter(svo_str)
    wave = filt.waveset
    dwave = 0.5 * (np.r_[[0], np.diff(wave)] + np.r_[np.diff(wave), [0]])

    flux = spec(wave)  # ph/s/cm2/AA
    flux *= filt(wave) * dwave  # ph/s
    n_ph = np.sum(flux.to(u.ph / u.s / u.m ** 2).value)

    # print(f"\nVega spectrum over the {filter_name} band "
    #       f"({filter_name}=0mag) has a flux of "
    #       f"{int(n_ph * 1e-6)}e6 ph/s/m2")

    return n_ph


def in_vega_spectrum(filter_name, instrument=None, observatory=None):
    return for_magnitude_in_filter(filter_name, 0, instrument, observatory)


in_vega_spectrum.__doc__ = for_magnitude_in_filter.__doc__


def in_ab_spectrum(filter_name, instrument=None, observatory=None):
    return for_magnitude_in_filter(filter_name, 0*u.ABmag, instrument, observatory)


in_ab_spectrum.__doc__ = for_magnitude_in_filter.__doc__
