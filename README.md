# How Many Photons
Tiny repo to calculate astronomical photon counts for various instruments

Note: the python package must be called using the abbreviation ``hmbp``.

The abbreviation ("how many *bloody* photons") is an expression of the authors frustration with the various
convoluted legacy units still used in modern day Astronomy. 


## Install

    pip install HowManyPhotons 


## Basic usage

HowManyPhotons is mainly used to determine how many photons fall within the
the wavelength range of any instrument filter contained in the spanish VO
filter database.

```python
from astropy import units as u
import hmbp

hmbp.for_flux_in_filter("V", 10*u.ABmag)
hmbp.convert(from_quantity=10*u.mJy, to_unit=u.mag, filter_name="J", 
             instrument="HAWKI", observatory="Paranal")
```


## Input and output units

``hmbp`` accepts and converts between ``u.mag`` (Vega), ``u.ABmag`` and ``u.Jy``.

ALL functions return photons counts in units of ``[ph s-1 m-2]``.

ONLY the filter curve is included in the calclutation. 
The following instrumental transmission profiles are NOT included in the photon 
flux calculation:

- atmospheric transmission
- mirror reflection
- detector quantum efficiency

These spectral profiles may be included in later releases of
HowManyBloodyPhotons, but they are currently NOT considered. 


## Main functions

The are two main functions: ``for_flux_in_filter`` and ``convert``:

### hmbp.for_flux_in_filter

Returns the number of incoming photons through a specific filter. 
If no ``instrument`` and ``observatory`` are provided, ``hmbp`` looks for a 
corresponding ``filter_name`` in the dictionary ``hmbp.FILTER_DEFAULTS``.

The result is an ``astropy.Quantity`` with the units ``[ph s-1 m-2]``.

Function signature: 

``
hmbp.for_flux_in_filter(filter_name, flux, instrument=None, observatory=None)
``

Some short examples:

```python
hmbp.for_flux_in_filter("V", 10*u.ABmag)
hmbp.for_flux_in_filter("Ks", 20*u.mag, instrument="HAWKI", observatory="Paranal")
hmbp.for_flux_in_filter("Si6", 1*u.Jy, instrument="Michelle", observatory="Gemini")
```


### hmbp.convert

Converts one common flux unit into another common flux: 
(``mag``, ``ABmag``, ``Jy``)

Function signature: 

``
hmbp.convert(from_quantity, to_unit, filter_name, instrument=None, observatory=None)
``

Some short examples:

```python
hmbp.convert(10*u.mag, u.Jy, "BrGamma")
hmbp.convert(from_quantity=0*u.mag, to_unit=u.ABmag, filter_name="J",
             instrument="HAWKI", observatory="Paranal")
```


## Convenience functions

We have also provided a few helper functions for several common flux conversions:

- ``hmbp.in_zero_vega_mags``
- ``hmbp.in_zero_AB_mags``
- ``hmbp.in_one_jansky``

The function signatures follow the same pattern as ``hmbp.for_flux_in_filter``,
just without needing to explicitly specify the flux parameter.

Returned units are ``[ph s-1 m-2]``

Some short examples:

```python
hmbp.in_zero_vega_mags("V")
hmbp.in_zero_AB_mags("Ks", "HAWKI", "Paranal")
hmbp.in_one_jansky("NeII", instrument="VISIR", observatory="Paranal")
```
