Scripts for analysing spectroscopy, primarily from IDS

Process for reducing IDS spectra:
1. Run Spector.py on a nightly basis and extract all the 1D spectra
2. Run INT_TimeCorrection.py to fix the times/airmass in the headers. Values from IRAF are a mess. SetJD and SetAirmass disabled in Spector as of 20151006
3. Run LogSpectraToDB.py on each night's data to log the images to the database. This will be added to later when FXCOR is ran
4. Check for common template RV standards across all nights
5. Analyse the single stars (n_traces=1) first. Blends need a more carefull analysis afterwards (n_traces>1)
6. Run FXCOR.py to get RV shifts. FXCOR will log all the values to the database
7. Repeat 5 and 6 for the BLENDs. 

In step 3 we log the following to the database for easier analysis later
	1. image_id
	2. object_id
	3. hjd_mid
	4. n_traces
	5. template_id
	6. ref_star

Then in step 5 (and 7) we add the following values from FXCOR:
	1. peak_shift_pix       
	2. correlation_height   
	3. fwhm_peak_pix        
	4. fwhm_peak_kms        
	5. relative_velocity_kms
	6. observed_velocity_kms
	7. helio_velocity_kms   


BLENDSs need more attention to extract the right trace ID (e.g. PA | PA+180)
