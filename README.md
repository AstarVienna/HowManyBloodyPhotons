# How Many Bloody Photons
Tiny repo to calculate astronomical photon counts for various instruments

### Install

    pip install git+https://github.com/AstarVienna/HowManyBloodyPhotons.git 


### Basic usage

```python
>>> from astropy import units as u
>>> import hmbp
```

Find the number of photons emitted by Vega in V band [ph s-1 m-2]::
         
```python
>>> hmbp.in_vega_spectrum("V")
```

Find the number of photons emitted by a Ks=20 [ABmag] point source through
the HAWKI Ks filter::
    
```python
>>> hmbp.for_magnitude_in_filter("Ks", 20*u.ABmag, 
                                 observatory="Paranal", 
                                 instrument="HAWKI")
```
