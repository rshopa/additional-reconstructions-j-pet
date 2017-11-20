# TomoPy FBP image reconstruction

Implementation of 3D FBP (somewhat simplified, so 2D+1) for the simulated J-PET scanner, using the open-source [Python package for tomographic data processing and image reconstruction](tomopy.readthedocs.io "TomoPy"). 

The scripts from *src/* allow to operate with both the single data file and multiple ones.

## One data file (*get_slices.py*)

The file *get_slices.py* operates with the convenience used for file names (data is in list mode ASCII format) -- the prefix "*PSF_*" and the geometry set at the end such as "*_x10_y0_z0*". Example of the name: *PSF_384strips_370kBq_600s_L050_x10_y0_z1875*. The corresponding ending will be used for naming the output data. The geometry in the example (in *data/* folder) reflects only the scanner length of *L*=50 cm, so that *z*-coordinate of the source is *z*=0 cm or *z*=18.75 cm. 

Assuming that all files are located in the place which the shell command is executed from, the basic usage of the script is as follows:
```
$ python get_slices.py <input_file> [biased] [centered]
```
It is implied that the option *biased* added at the end would mean that *z*=18.75 cm and it is verified in the script. The option *centered* defines the usage of time-of-flight (TOF) information: if present -- *z*-coordinate of each LOR is estimated from its centre, otherwise hit time differences are considered. The order of options *biased* and *centered* in the shell command is not important. 

(!) Note that the word "centered" is written in American English (DO NOT use British spelling "centred").

The usage related to the examples, using TOF data (will lead to more accurate reconstruction):
```
$ python get_slices.py PSF_384strips_370kBq_600s_L050_x10_y0_z0
$ python get_slices.py PSF_384strips_370kBq_600s_L050_x10_y0_z1875 biased
```
Without TOF (using *centered*):
```
$ python get_slices.py PSF_384strips_370kBq_600s_L050_x10_y0_z0 centered
$ python get_slices.py PSF_384strips_370kBq_600s_L050_x10_y0_z1875 centered biased
```
As the result, the script will return three vectors (*x_axis*, *y_axis*, *z_axis*, in ASCII), defining the scale (in centimeters), and the data, stored as ASCII matrices (tab as a separator), for no filtering (*none*) and FBP (*fbp*) used for the reconstruction. The output files for the examples *PSF_384strips_370kBq_600s_L050_x10_y0_z0* above will be:
- x_axis
- y_axis
- z_axis
- xy_fbp_x10_y0_z0
- xz_fbp_x10_y0_z0
- xy_none_x10_y0_z0
- xz_none_x10_y0_z0

If the source was defined correctly and the word *biased* used properly, the script will return such matrices for the planes *XY* and *XZ*, made along the *z*- and *y*-position of the center of the source.  For *PSF_384strips_370kBq_600s_L050_x10_y0_z1875*, only the ending will differ (such as *xz_none_x10_y0_z1875* etc).

## Multiple files (*execute.sh*)

The usage of the script *execute.sh* implies that all the files -- the data, the class *ListModeReconstruction.py* (including compiled *ListModeReconstruction.pyc*), *get_slices.py* -- are located at the same directory which the command is executed from:
```
$ ./execute.sh [centered]
```
Consequently, the script will execute *get_slices.py* over all data files, with the automatic definition, whether the option *biased* is needed, by parsing the name (ending with "_z1875" or "_z0"). The option *centered*, related to TOF and explained above, is applied to all data files. 

The correct execution for this repository requires the files from *src/* to be copied to *data/* and executed from there or vice versa.
