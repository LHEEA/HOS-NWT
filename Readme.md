# HOS-NWT

Numerical Wave Tank based on High-Order Spectral method

[![Travis][buildstatus_image_travis]][travisci]
[![Appveyor][buildstatus_image_appveyor]][appveyorci]

This README file describes the different cases that might be computed with HOS-NWT
and gives instructions to set the numerical parameters `n1` to `p2`

Setting the value of integers `n1`, `n2`, `M`, `p1` and `p2` in `common_vars.f90`

## 2D simulation

For a 2D simulation,
   Compile with `n2=1` `AND` `p2=1` to adjust the memory allocation to minimum

   If partial dealiasing is used, compile with `p1` set to maximal required value
   (total dealiasing is obtained with `p1=M`
   but it can be reduced if `p1` is further set to a value below `M`)

## 3D simulation

For a 3D simulation,
   Compile with `n2\=1` `AND` `p2` set to required value

   If partial dealiasing is used in x-direction,
   compile with `p1` set to maximal required value (total dealiasing is obtained with `p1=M`
   but it can be reduced if `p1` is further set to a value below `M`)

   If partial dealiasing is used in y-direction,
   compile with `p2` set to maximal required value (total dealiasing is obtained with `p2=M`
   but it can be reduced if `p2` is further set to a value below `M`)

* 'input_HOS-NWT.dat' will have to be attached to run this program *

Setting the value of integer `icase` in `input_HOS-NWT.dat`

- `icase = 1` : Sloshing case
    Computation starts with a natural mode of the tank (in x) of a given amplitude

- `icase = 2` : Monochromatic case
    Regular wave is generated in the NWT
    User defines frequency, amplitude, phase (and angle of propagation if 3D simulation)

- `icase = 3` and `31` and `32` and `33` : File case
    Wavemaker movement is deduced from input file named 'filename'
    - `3`  - `filename.cfg` describes the configuration of wavemaker
        `filename.dat` describes the frequency components of wavemaker movement
    - `31` - `filename.txt` is an output of control software used in ECN Wave Basin
    - `32` - `filename.txt` is an output of control software used in ECN Towing Tank
    - `33` - `filename.txt` is an output of control software used in other tanks

- `icase = 4` and `41` : Irregular wave
   Wavemaker movement creates an irregular wave field with a given Hs and Tp
    - `4`  - JONSWAP spectrum
    - `41` - Bretschneider spectrum

Further details about input file, output of the code... may be find at the Wiki page of HOS-NWT project: https://github.com/LHEEA/HOS-NWT/wiki

***

## Solution method

### Numerical tools

    1. fully FFT resolution for spectral and additionnal modes.
        1.1. Fully nonlinear HOS for free surface
        1.2. 1st, 2nd and 3rd order wavemaker modeling
    2. Runge-Kutta 4th-order scheme in time.


### Coordinates system

             ^ Y (the width of wave tank)
             |
             |
             O --------->X (the length of wave tank)


### Discretization of wave field

      Y ^  eta(1,n2) eta(2,n2) --------------- eta(n1,n2)
        |    .         .                          .
        |    .         .                          .
        |  eta(1,3)    .                          .
        |  eta(1,2)  eta(2,2) ----------------- eta(n1,2)
        |  eta(1,1)  eta(2,1) eta(3,1) -------- eta(n1,1)
      O |-------------------------------------------> X


### Nondimensionalizations

    [L] = h
    [T] = 1/sqrt(g/h)


[buildstatus_image_travis]: https://travis-ci.org/LHEEA/HOS-NWT.svg?branch=master
[travisci]: https://travis-ci.org/LHEEA/HOS-NWT

[buildstatus_image_appveyor]: https://ci.appveyor.com/api/projects/status/9pu4d2njd540ffkp?svg=true
[appveyorci]: https://ci.appveyor.com/project/Gjacquenot/hos-nwt
