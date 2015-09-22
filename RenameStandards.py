
# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# RenameStandards.py - Plot scintillation with exposure time 

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 22/07/11	- 	code writen
# 22/07/11	-	code tested 
                  
import os
import commands as cmd

t=cmd.getoutput('ls i_s_r*.fits').split('\n')

for i in range(0,len(t)):
	f1=t[i]
	f2="%s-S.ms.fits" % (t[i].split('.')[0])
	
	os.system('cp %s %s' % (f1, f2))