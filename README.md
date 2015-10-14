Scripts for analysing spectroscopy, primarily from IDS

Process for reducing IDS spectra: <br/>
1. Run Spector.py on a nightly basis and extract all the 1D spectra <br/>
2. Run INT_TimeCorrection.py to fix the times/airmass in the headers. Values from IRAF are a mess. SetJD and SetAirmass disabled in Spector as of 20151006 <br/>
3. Run LogSpectraToDB.py on each night's data to log the images to the database. This will be added to later when FXCOR is ran <br/>
4. Check for common template RV standards across all nights <br/>
5. Analyse the single stars (n_traces=1) first. Blends need a more carefull analysis afterwards (n_traces>1) <br/>
6. Check that FXCOR settings are coorect (osample etc)
7. Run FXCOR.py to get RV shifts. FXCOR will log all the values to the database <br/>
8. Repeat 6 and 7 for the BLENDs.  <br/>
<br/>
In step 3 we log the following to the database for easier analysis later <br/>
	1. image_id <br/>
	2. object_name <br/>
	3. template_image_id <br/>
	4. ref_star_name <br/>
	5. hjd_mid <br/>
	6. utmiddle <br/>
	7. n_traces <br/>	
<br/>
Then in step 7 (and 8) we add the following values from FXCOR: <br/>
	1. peak_shift_pix <br/>
	2. correlation_height <br/>   
	3. fwhm_peak_pix <br/>
	4. fwhm_peak_kms  <br/>       
	5. relative_velocity_kms  <br/>
	6. observed_velocity_kms  <br/>
	7. helio_velocity_kms  <br/> 
	8. commnets <br />
	9. swasp_id - populated later <br/> 
<br/>

We also add the following eblm parameters to the database for easier analysis using the EBLMParamsToDB.py script: <br/>
	1. swasp_id <br/>
	2. period <br/>   
	3. epoch <br/>
	4. Vmag  <br/>       
	5. paramfit_spec_type  <br/>

We then cross match the two tables to add the swasp_id to the spectra <br/>

BLENDSs need more attention to extract the right trace ID (e.g. PA | PA+180) <br/>
