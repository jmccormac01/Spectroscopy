Scripts for analysing spectroscopy, primarily from IDS

Data is stored in the eblm_ids_new table <br/>
<br/>

Process for reducing IDS spectra: <br/>
   1. Read the night log, remove any bad images etc
   1. Run Spector.py on a nightly basis and extract all the 1D spectra 
   1. Run INT_TimeCorrection.py to fix the times/airmass in the headers. Values from IRAF are a mess. SetJD and SetAirmass disabled in Spector as of 20151006
   1. Apply the radial velocity correction using RVcorrect.py
   1. Run LogSpectraToDB.py on each night's data to log the images to the database. This will be added to later when FXCOR is ran 
   1. Check for common template RV standards across all nights 
   1. Analyse the single stars (n_traces=1) first. Blends need a more careful analysis afterwards (n_traces>1) 
   1. Check that FXCOR settings are coorect (osample etc)
   1. Run FXCOR.py to get RV shifts. FXCOR will log all the values to the database 
   1. Repeat 6 and 7 for the BLENDs. 
<br/>
In step 3 we log the following to the database for easier analysis later: <br/>
   1. image_id 
   1. object_name 
   1. template_image_id 
   1. ref_star_name 
   1. hjd_mid 
   1. utmiddle 
   1. n_traces 	
<br/>
Then in step 7 (and 8) we add the following values from FXCOR: <br/>
   1. peak_shift_pix 
   1. correlation_height    
   1. fwhm_peak_pix 
   1. fwhm_peak_kms        
   1. relative_velocity_kms  
   1. observed_velocity_kms  
   1. helio_velocity_kms  
   1. comments 
<br/>

We also add the following eblm parameters to the database for easier analysis using the EBLMParamsToDB.py script: <br/>
   1. swasp_id 
   1. period 
   1. epoch 
   1. Vmag       
   1. paramfit_spec_type  

We then cross match the two tables to add the swasp_id to the spectra <br/>

BLENDSs need more attention to extract the right trace ID (e.g. PA | PA+180) <br/>
I've now updated FXCOR to account for multiple objects in the slit. The plan is to use the SkyPA to distinguish between the spectra. 

For the blended objects use the correctBlends.py script to fix the trace ids. MeasureRVs.py can then use the updated trace IDs to plot the right spectra with each other.  <br/>

