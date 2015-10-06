Scripts for analysing spectroscopy, primarily from IDS

Process for reducing IDS spectra: <br/>
1. Run Spector.py on a nightly basis and extract all the 1D spectra <br/>
2. Run INT_TimeCorrection.py to fix the times/airmass in the headers. Values from IRAF are a mess. SetJD and SetAirmass disabled in Spector as of 20151006 <br/>
3. Run LogSpectraToDB.py on each night's data to log the images to the database. This will be added to later when FXCOR is ran <br/>
4. Check for common template RV standards across all nights <br/>
5. Analyse the single stars (n_traces=1) first. Blends need a more carefull analysis afterwards (n_traces>1) <br/>
6. Run FXCOR.py to get RV shifts. FXCOR will log all the values to the database <br/>
7. Repeat 5 and 6 for the BLENDs.  <br/>
<br/>
In step 3 we log the following to the database for easier analysis later <br/>
	1. image_id <br/>
	2. object_id <br/>
	3. hjd_mid <br/>
	4. n_traces <br/>
	5. template_id <br/>
	6. ref_star <br/>
<br/>
Then in step 5 (and 7) we add the following values from FXCOR: <br/>
	1. peak_shift_pix <br/>
	2. correlation_height <br/>   
	3. fwhm_peak_pix <br/>
	4. fwhm_peak_kms  <br/>       
	5. relative_velocity_kms  <br/>
	6. observed_velocity_kms  <br/>
	7. helio_velocity_kms  <br/>  
<br/>
BLENDSs need more attention to extract the right trace ID (e.g. PA | PA+180) <br/>
