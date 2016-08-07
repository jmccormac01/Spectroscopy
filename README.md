Scripts for analysing spectroscopy, primarily from IDS

Data is stored in the eblm\_ids\_new table <br/>
<br/>

Process for reducing IDS spectra: <br/>
   1. Read the night log, remove any bad images etc
   1. Run Spector.py on a nightly basis and extract all the 1D spectra 
      1. In indentify use the following settings:
         * thres 700 
         * mark 'm' the two lines around 6506.53 and 6532.88 (to the right of big gap in middle, if cenwave~6600A)
         * use 'l' to find more lines
         * niter 7   
         * 'f' to enter fitting window
         * order 3
         * thres 50 to pick up all the faint lines
   1. Set read/write/execute permissions on all the data. If not IRAF shits itself
      1. chmod -R ugo+rwx *
   1. Run INT\_TimeCorrection.py to fix the times/airmass in the headers. Times from IRAF are a mess. SetJD and SetAirmass disabled in Spector as of 20151006. This has to run after RVcorrect or IRAF ruins the HJD value!!
   1. Apply the radial velocity correction using RVcorrect.py
   1. Run INT\_TimeCorrection.py again to fix the small HJD correction error from IRAF
   1. Run LogSpectraToDB.py on each night's data to log the images to the database. 
   1. Add the following eblm parameters to the database for easier analysis using 
   the EBLMParamsToDB.py script. In cases where the swasp\_id is not found, 
   locate them on hunter or in the wasp db and add their values.  <br/>
       1. swasp\_id 
       1. period 
       1. epoch 
       1. Vmag       
       1. paramfit\_spec\_type  

Once all nights have been reduced and ingested into the database: <br/>
	1. Check for common template RV standards across all nights 
	1. Analyse the single stars (n\_traces=1) first. Blends need a more careful analysis afterwards (n\_traces>1) 
	1. Check that FXCOR settings are coorect (osample etc)
	1. Run FXCOR.py to get RV shifts. FXCOR will log all the values to the database 
	1. Repeat for the BLENDs. 
<br/>
In step 6 we log the following to the database for easier analysis later: <br/>
   1. image\_id 
   1. object\_name 
   1. template\_image\_id 
   1. ref\_star\_name 
   1. hjd\_mid 
   1. utmiddle 
   1. n\_traces 	
<br/>
Then in step 8 (and 9) we add the following values from FXCOR: <br/>
   1. peak\_shift\_pix 
   1. correlation\_height    
   1. fwhm\_peak\_pix 
   1. fwhm\_peak\_kms        
   1. relative\_velocity\_kms  
   1. observed\_velocity\_kms  
   1. helio\_velocity\_kms  
   1. comments 
<br/>



We then cross match the two tables to add the swasp\_id to the spectra <br/>

BLENDSs need more attention to extract the right trace ID (e.g. PA | PA+180) <br/>
I've now updated FXCOR to account for multiple objects in the slit. The plan is to use the SkyPA to distinguish between the spectra. 

For the blended objects use the correctBlends.py script to fix the trace ids. MeasureRVs.py can then use the updated trace IDs to plot the right spectra with each other.  <br/>


UPDATE
------
This process is being redone. 
   1. Check the night log and remove all spectra that have issues
   1. Run Spector.py
      1. In indentify use the following settings:
         * thres 700 
         * mark 'm' the two lines around 6506.53 and 6532.88 (to the right of big gap in middle, if cenwave~6600A)
         * use 'l' to find more lines
         * niter 7   
         * 'f' to enter fitting window
         * order 3
         * thres 50 to pick up all the faint lines
