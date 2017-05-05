Spectroscopy
-------------

Spectroscopy is carried out using 1 of two codes depending on the instrument.

   1. Spector is a package for reducing long-slit IDS data
   1. For echelle spectra we use the [CERES](https://github.com/rabrahm/ceres) package from Rafael Brahm. *Update these docs to include CAFE/FIES info?*

Spector Usage
-------------

*Before running Spector:*
   1. ```chmod -R ugo+rwx *``` so ```Spector``` has permission to edit the files
   1. Remove all files that should not be included in reduction
      i.e. bad flats, data from other central wavelengths or gratings etc.
      Place all this data in a folder called e.g. 'junk'. Use *.int logfile to
      help weed out unwanted files
   1. Adjust the settings in the **EDIT HERE** section below

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
into one location for RV analysis. The files from multiple nights are
linked by their ```SWASP\_ID```. First populate the database column with
each spectrum's ```SWASP\_ID``` using:

```python database/setSwaspIds.py```

Manually add any missing entries to the database and/or resolve naming
or typos at the telescope. Add the ```--update`` flag to update the database

Next update the Q1 and Q2 keywords so they next observations can be planned:

```python database/updateRvQuadStatus.py```

use the ```--commit``` flag after checking the quads are ok. This will 
commit them to the database

Once all the ```SWASP_ID``` keywords have been updated and quads have been
populated, gather all the spectra into per-object folders for analysis using:

```python gatherSpectraFor_iSpec.py```

This gets the unique list of targets from the database and puts all the 
spectra together in one folder per object ready for analysis with iSpec.
Add the ```--copy``` flag to copy the data to the combined area.

Dealing with blends
-------------------

# WE NEED TO FIX BLENDS HERE! SET THEIR SWASPIDS RIGHT AWAY TO NOT LEAVE THEM BEHIND!

Estimating Spectral Type
------------------------

To use the right CCF mask with iSpec we need to estimate the spectral type
of the target. This can be done by getting the Teff_jh and Teff_vk from Hunter
and then using the Teff-to-Spectral_Type translation to estimate the type
based on the average of the two temperatures. This is not super accurate but
should be enough to distinguish F from K or M stars.

Measuring RVs with iSpec
------------------------

Loop over all the objects with the ```iSpec_RVs.py``` script. This will
cross correlate spectra with the first one to determine the RVs. A second
cross correlation is also done with the first spectrum as a double check.

By default ```iSpec_RVs.py``` will loop over all the objects that have spectra
that have not had their RV measured. You can force it to analyse a given object
by specifiying ```--swasp_id```.

The barycentric correction is not the most accurate, but this is ok for
initial analysis. The BCV should be calculated accurately before any detailed
analysis.

Determining Outstanding Observations
------------------------------------

This is now done by ```updateRvQuadStatus.py``` above. Run this after each
reduction run to always have the correct current state of the world for the
project

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



