Scripts for analysing spectroscopy, primarily from IDS

Data is stored in the eblm_ids_new table <br/>
<br/>

Process for reducing IDS spectra: <br/>
1. Run Spector.py on a nightly basis and extract all the 1D spectra <br/>
2. Run INT_TimeCorrection.py to fix the times/airmass in the headers. Values from IRAF are a mess. SetJD and SetAirmass disabled in Spector as of 20151006 <br/>
3. Apply the radial velocity correction using RVcorrect.py<br/>
4. Run LogSpectraToDB.py on each night's data to log the images to the database. This will be added to later when FXCOR is ran <br/>
5. Check for common template RV standards across all nights <br/>
6. Analyse the single stars (n_traces=1) first. Blends need a more careful analysis afterwards (n_traces>1) <br/>
7. Check that FXCOR settings are coorect (osample etc)
8. Run FXCOR.py to get RV shifts. FXCOR will log all the values to the database <br/>
9. Repeat 6 and 7 for the BLENDs.  <br/>
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
	8. comments <br />
	9. ccf height <br />
	10. tonry & davis r-value <br />
<br/>

We also add the following eblm parameters to the database for easier analysis using the EBLMParamsToDB.py script: <br/>
	1. swasp_id <br/>
	2. period <br/>   
	3. epoch <br/>
	4. Vmag  <br/>       
	5. paramfit_spec_type  <br/>

We then cross match the two tables to add the swasp_id to the spectra <br/>

BLENDSs need more attention to extract the right trace ID (e.g. PA | PA+180) <br/>
I've now updated FXCOR to account for multiple objects in the slit. The plan is to use the SkyPA to distinguish between the spectra. 

For the blended objects use the correctBlends.py script to fix the trace ids. MeasureRVs.py can then use the updated trace IDs to plot the right spectra with each other.  <br/>

