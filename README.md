# igrowup_update
About WHO IGROWUP

IGROWUP is a stata macro developed by the Department of Nutrition and Food Safety.


More details on the license can be found here: https://www.who.int/about/policies/publishing/copyright


Update of IGROWUP Stata Macro Prevalence Module


**Contact for reporting bugs/ comments: **
Should you encounter any problems with this package, please send an e-mail with a clear description of the identified problem to "anthro2005@who.int", specifying in the subject line that it concerns the igrowup_Stata package , the name of the macro (igrowup_standard or igrowup_restricted) and kindly indicate which version of STATA 
you are using. Thank you.


The package igrowup_stata contains the following items: 
	1. Two macros (igrowup_standard.ado and igrowup_restricted.ado). 
	2. Nine permanent (read-only) Stata data sets containing the WHO Child Growth Standards: weianthro.dta, lenanthro.dta, wflanthro.dta, wfhanthro.dta, bmianthro.dta, 
	hcanthro.dta, acanthro.dta, tsanthro.dta and ssanthro.dta. 
	3. The file Readme.pdf 
	4. An example set, survey.dta. 
	5. Two example do-files, survey_standard.do and survey_restricted.do. 
	6. The example output files: mysurvey_z_st.xls, mysurvey_z_rc.xls, mysurvey_z_st.dta, mysurvey_z_rc.dta, mysurvey_prev_st.xls and mysurvey_prev_rc.xls.

**Pre-requisites **
**Individual level Z-score outputs**
STATA Version 7.0 Stata/SE (Special Edition of Stata) or higher is required to run two macros (igrowup_standard.ado and igrowup_restriced). 
Intercooled Stata has a limit of 2,047 variables and with that the macros will only produce the z-scores output files

**Prevalence Estimates**
STATA Version 13.0 Stata/SE (Special Edition of Stata) or higher is required to generate prevalence estimates

**Precautions: **
1. Avoid any variable names starting with underscore "_" in the input STATA data set; otherwise they may be replaced by the derived ones created by the macro. 
2. Avoid any temporary format names starting with underscore "_"; otherwise they may be replaced by the temporary ones created by the macro. 
3. Avoid any STATA global macro variable names staring with underscore “_”, except those defined by the system. 

**Recommended setup and run **

Step 1. 
Create a sub-directory, for example "D:\WHO igrowup STATA", where you wish to save the package (igrowup_stata.zip). This directory should be reserved only for the references tables (*anthro.dta) and the macros  (igrowup_standard.ado and igrowup_restricted) that are contained in the zip file.

Step 2. 
Create a sub-directory, for example "D:\WHO igrowup workdata", where the example data (survey.dta and pertaining output files) and your STATA input data can be stored and where all the macro output files will be written to.

Step 3. 
It is recommended that you start by loading and running the example code below (also found in survey_standard.do or survey_restricted.do) in the STATA do-file  editor to see how the data should be prepared and to fill in the macros' parameters according to their requirements. Note: The macros run on Stata/SE: The Special Edition version of Stata (type "help SpecialEdition"). For users of Intercooled Stata, the macros will only produce the z scores output files and the user gets the message 
 
 Please wait, programme is calculating prevalences............. 
 .............................................. 
no room to add more variables

