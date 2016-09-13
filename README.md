Spector
------

A package for reducing long-slit IDS data

Usage
-----

*Before running Spector:*
    -   Remove all files that should not be included in reduction
        i.e. bad flats, data from other central wavelengths or gratings etc.
        Place all this data in a folder called 'junk'. Use *.int logfile to
        help weed out unwanted files
    -   Adjust the settings in the **EDIT HERE** section below

*Assumptions:*
   1. This code was written for a particular observing strategy - Acquire, Arc, Science, Arc, repeat. Therefore the script is looking for this structure in the data when matching arcs to spectra. This can be easily changed by modifying the relevant section

*What Spector does:*
   1. Prepares the working environment, making copies of original data
   1. Strips multi extension fits to single extension
   1. Renames files for easier viewing:
      1. i_* = integration (spectra)
      1. a_* = arc
      1. b_* = bias
      1. f_* = flat
   1. Debiases and flat field corrects the data using CCDPROC in ASTROPY
   1. Trims the spectra based on INT's heavily vignetted region
   1. Calcuates mid exposure datetime in BJD\_TDB, HJD and JD formats
   1. Tidies up calibration frames into 'calibs' folder
   1. Interactive spectral analysis with KPNOSLIT in IRAF
   1. Normalises the spectra
   1. Tidies up working directory after reductions are finished
   1. Removes all intermediate data products to save disc space
   1. Logs spectra information to the databse

*What Spector does not do:*
   1. Flux calibration

*What you will end up with:*
   1.  A set of reduced, wavelength calibrated, normalised 1D spectra.

Gathering Spectra for Analysis with iSpec
-----------------------------------------

After reduction and extraction with Spector the files should be combined
into one location for RV analysis using:

```python gatherSpectraFor_iSpec.py```

This gets the unique list of targets from the database and puts all the spectra
together in one folder per object for analysis with iSpec.

Measuring RVs with iSpec
------------------------

Loop over all the objects with the ```iSpec_RVs.py``` script. This will
cross correlate spectra with the first one to determine the RVs.

The barycentric correction is not the most accurate, but this is ok for
initial analysis. The BCV should be calculated accurately before any detailed
analysis.

Determining Outstanding Observations
------------------------------------

With a partially complete dataset we do the following before each run to 
see what observations are required. Check the Q1 and Q2 columns of each 
object. 

If !Q1 and !Q2, then try for both

If !Q1 only, then do Q1 (same for Q2). 

If Q1 and Q2 and dRV > 25 km/s, then flag for phase coverage

If Q1 and Q2 and dRV < 25 km/s, then flag for phase coverage-stabilized

Ignore all SB2, IGNORE, EB, BEB objects

Dealing with blends
-------------------

Need to determine a way to treat blends properly

Schema for database table
-------------------------

```
create table eblm_ids (
image_id                      varchar(40) not null primary key,
swasp_id                      varchar(40),
object_name                   varchar(20) not null,
night                         date,
utmiddle                      datetime,
hjd_mid                       double,
bjd_mid                       double,
jd_mid                        double,
n_traces                      int(2),
sky_pa                        double,
ccf_height                    double,
ccf_fwhm                      double,
atomic_velocity               double,
atomic_velocity_err           double,
telluric_velocity             double,
telluric_velocity_err         double,
barycentric_velocity_iSpec    double,
barycentric_velocity_exo      double,
comment                       varchar(40),
analyse                       tinyint(1)
);
```

Contributors
------------

James McCormac



