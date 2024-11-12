`ORG.GLC:`
![GitHub release (latest by date)](https://img.shields.io/github/v/release/Geo-Linux-Calculations/ggmcalc)
![GitHub Release Date](https://img.shields.io/github/release-date/Geo-Linux-Calculations/ggmcalc)
![GitHub repo size](https://img.shields.io/github/repo-size/Geo-Linux-Calculations/ggmcalc)
![GitHub all releases](https://img.shields.io/github/downloads/Geo-Linux-Calculations/ggmcalc/total)
![GitHub](https://img.shields.io/github/license/Geo-Linux-Calculations/ggmcalc)  

# GGMCalc

## Upgrate by Roman Sermiagin (2024-11)

### Upgrate

* Add CLI interface for loading data files
* Add installation to `\opt\bin\` folder for Linux

### ToDO

* Add to CLI interface arguments for the calculation setting 

## About

GGMCalc is a program to compute some components of the gravity field, using
Global Geopotential Models, including

1. Undulation
2. Heigh anomaly
3. Gravity disturbance
4. Classical gravity anomaly
5. Molodensky gravity anomaly

## Reference

Moazezi S. and Zomorrodian H. (2012) GGMCalc a software for calculation of the geoid undulation and the height anomaly using the iteration method, and classical gravity anomaly, Earth Science Informatics, 5(2), pp. 123-136. https://doi.org/10.1007/s12145-012-0102-2

## Building and Installing

`gfortran` and "GNU Make" should be installed on the OS. Using `make` command
one can build `ggmcalc`.

```bash
    make
    sudo make install
```

## Usage

```bash
Usage: ggmcalc [OPTION]... [FILE]...

Reads EGM coefficients and point coordinates from specified FILEs,
then calculates height anomalies, geoid undulations, and other values.

Options:
  --coeffs        Specify the EGM coefficients file in ICGEM format.
  --ellip_conf    Specify the ellipsoid configuration file with the following format:

                  3986004.418d8
                  6378137.0d0
                  6356752.3142d0
                  298.257223563
                  7292115.0d-11
                  62636851.7146d0
                  WGS_84

  --input_files   Specify the input file(s) with coordinates and other required fields.

  --output        Specify the path for the output CSV file with results.

  --help          Display this help message and exit.
```

The global geopotential model (GGM) file by option `--coeffs` is used by the
software. One can download GGM files from ICGEM website at
<http://icgem.gfz-potsdam.de/ICGEM/> and specity path to it.

The ellipsoid parameters file by option `--ellip_conf` is used by the software.

* `dat/ellpsoid.dat.GRS80` contains parameters of GRS80 ellipsoid.
* `dat/ellpsoid.dat.WGS84` contains parameters of WGS84 ellipsoid.

To use a particular ellipsoid, specify path to it.

The `dat/input.dat` text file (option `--input_files`) contains the coordinates of
the calculating points (`longitude latitude height` format).

Descriptions of the functions and the subroutines including corresponding
references can be found in the source codes.

## License

Copyright (C) 2010-2018 Siamak Moazezi

This file is part of GGMCalc.

GGMCalc is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

GGMCalc is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
GGMCalc. If not, see <http://www.gnu.org/licenses/>.

Contact info:

* [http://www.sourceforge.net/projects/xgravity](http://www.sourceforge.net/projects/xgravity)
* Siamak Moazezi <s.moazezi@srbiau.ac.ir>
