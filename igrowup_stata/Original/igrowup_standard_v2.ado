/*************************************************************************************************************/
/********  	WHO Child Growth Standards                                                         	  ****/
/********  	Department of Nutrition for Health and Development                                 	  ****/
/********  	World Health Organization                                                          	  ****/
/********  	Last modified on 23/05/2007     -     For STATA versions 12.0 and above             	  ****/
/*************************************************************************************************************/

/*************************************************************************************************************/
/********  	A macro / program for calculating the z-scores and prevalences for a nutritional survey  *****/
/*************************************************************************************************************/

*! version 1	  01may2007
program define igrowup_standard
	version 7.0
	discard
	set type float
	args reflib datalib datalab sex age ageunit weight height measure hc armc tris subs oedema sw 
	if ("`age'" == "" | "`ageunit'" == "" | "`weight'" == "" | "`height'" == "" | "`measure'" == ""| "`hc'" == "" | "`armc'" == "" | "`tris'" == "" | "`subs'" == "" | "`oedema'" == "" | "`sw'" == "") {
		di "Error - You must specify 15 arguments: "
		di "	reflib datalib datalab sex age ageunit weight"
		di "	height measure oedema hc armc tris subs oedema sw"
	exit
	}
	/*	Phase I =  generation of z-scores = _zwei, _zlen, _zwfl, _zbmi */
	di as txt _n "Please wait, programme is running............."
	di as txt _n ".............................................."

	/*	Data preparation */
	
	/* quietly begins*/
	
	qui {
	tempvar tsort tsex length2 height2 toedema lorh lenhei2 uselgth wsd3p wsd23p wsd3n wsd23n 
	tempvar bsd3p bsd23p bsd3n bsd23n interp temp1 lenhigh lenlow length9 
	tempvar zwhzll asd3p asd23p asd3n asd23n zwhzlh csd3p csd23p csd3n csd23n
	tempvar height9 interph temp2 hgtlow hgthigh zwhzhl dsd3p dsd23p dsd3n dsd23n zwhzhh esd3p esd23p esd3n esd23n
	tempvar abovel ratiol aboveh ratioh 
	tempvar msd3p msd23p msd3n msd23n
	tempvar tsd3p tsd23p tsd3n tsd23n
	tempvar ssd3p ssd23p ssd3n ssd23n
		 
	gen `tsort'=_n
	inspect `sex'
	if  r(N_unique)==0 {
		gen `tsex'=1 if `sex'=="m" | `sex'=="M" 
		replace `tsex'=2 if `sex'=="f" | `sex'=="F"
		replace `tsex'=. if `sex'== " " 
	}
	if  r(N_unique)~=0 {
		gen `tsex'=1 if `sex'==1
		replace `tsex'=2 if `sex'==2
		replace `tsex'=. if `sex'==.
	}
	qui ta `oedema'
	if  r(r)==0 {
		gen `toedema'=1  if `oedema'==.
	}
	if  r(r)~=0 {
		gen `toedema'= 1 
		replace `toedema'= 2 if `oedema'=="y" | `oedema'=="Y"
	}
	inspect `measure'
	if  r(N_unique)==0 {
		gen `lorh'= 1 if `measure'=="l" | `measure'=="L"
		replace `lorh'= 2 if `measure'=="h" | `measure'=="H"
		replace `lorh'= . if `measure'==" " | `measure'=="."
	}
	
	capture gen _agedays=`age' 
	if _rc==0 {
		replace _agedays=`age'*30.4375 if `ageunit'=="months"
		replace _agedays=round(_agedays,1)
		lab var _agedays "Calculated age in days"	
	}
	if _rc~=0 {
		drop _agedays
		gen _agedays=`age'
		replace _agedays=`age'*30.4375 if `ageunit'=="months"
		replace _agedays=round(_agedays,1)
		lab var _agedays "Calculated age in days for deriving z score"	
	}
	
	tempvar agetemp 
	gen `agetemp'=`age'
	replace `agetemp'=`age'*30.4375 if `ageunit'=="months"	
	
	gen `lenhei2'= `height'
	gen `uselgth'=-99
	replace `uselgth'=-99 if `lenhei2'==.
	replace `lenhei2'= `height'+.7 if (`lorh'==2 & _agedays<731) 
	replace `lenhei2'= `height'-.7 if (`lorh'==1 & _agedays>=731)
	replace `uselgth'=1 if (`lorh'==2 & _agedays<731)
	replace `uselgth'=2 if (`lorh'==1 & _agedays>=731)
	replace `uselgth'=1 if (`lorh'==1 & _agedays<731)
	replace `uselgth'=2 if (`lorh'==2 & _agedays>=731)
	
	* 	if missing the recumbent indicator but have age, we assume they have it right.
	replace `uselgth'=1 if (`lorh'==. &  _agedays<731)
	replace `uselgth'=2 if (`lorh'==. &  _agedays>=731)
	replace `lenhei2'= `height' if (`lorh'==1 & _agedays==.) 
	replace `lenhei2'= `height' if (`lorh'==2 & _agedays==.) 
	replace `uselgth'=1 if (`lorh'==1 & _agedays==.)
	replace `uselgth'=2 if (`lorh'==2 & _agedays==.)
	
	* 	if age missing & indicator missing, use length of child to figure.

	replace `uselgth'=1 if (`lorh'==. & _agedays==. &  `lenhei2'<87)
	replace `uselgth'=2 if (`lorh'==. & _agedays==. &  `lenhei2'>=87)

	* ======= check inputs & create outputs for dbf file ANTHRO check =========
	
	macro def under5 "if _agedays >= 61*30.4375"
	
	capture gen _clenhei= `lenhei2' 
	if _rc==0 {
		lab var _clenhei "Converted length/height for deriving z score (in cms)"	
	}
	if _rc~=0 {
		drop _clenhei
		gen _clenhei= `lenhei2'
		lab var _clenhei "Converted length/height for deriving z score (in cms)"	
	}

	capture gen double _cbmi= `weight'*10000/(`lenhei2'*`lenhei2')
	if _rc==0 {
		lab var _cbmi "Calculated bmi=weight / squared(_clenhei)"
	}
	if _rc~=0 {
		drop _cbmi
		gen double _cbmi= `weight'*10000/(`lenhei2'*`lenhei2')
		lab var _cbmi "Calculated bmi=weight / squared(_clenhei)"
	}
	
	* ======================== ZWEI ============================================
	sort `tsex' _agedays
	local string "xxx\weianthro.dta"
	local i=`reflib'
	global wazfile: subinstr local string "xxx" "`i'"
	
	
	merge `tsex' _agedays using "$wazfile"
	gen double _zwei=(((`weight'/m)^l)-1)/(s*l)
	gen double `wsd3p'=m*((1+l*s*3)^(1/l))
	gen double `wsd23p'=`wsd3p'- m*((1+l*s*2)^(1/l))
	replace _zwei=3+((`weight'-`wsd3p')/`wsd23p') if (_zwei>3 & _zwei~=.) 
	gen double `wsd3n'=m*((1+l*s*(-3))^(1/l))
	gen double `wsd23n'= m*((1+l*s*(-2))^(1/l))-`wsd3n'
	replace _zwei=(-3)-((`wsd3n'-`weight')/`wsd23n') if (_zwei<-3& _zwei~=.) 
	replace _zwei =. $under5
	keep if _merge~=2
	drop l m s _merge

	* ==========================  ZLEN  =======================================

	* 	if not done already, need age in days to compare to norms.
	*	select if (agedays2 < (61*30.4375) or missing(agedays2)).

	sort `tsex' _agedays
	local string "xxx\lenanthro.dta"
	local i=`reflib'
	global hazfile: subinstr local string "xxx" "`i'"
	
	merge `tsex' _agedays using "$hazfile"
	gen double _zlen=(((`lenhei2'/m)^l)-1)/(s*l)
	replace _zlen =. $under5
	keep if _merge~=2
	drop l m s loh _merge 
	
	* =============================== ZBMI ======================================

	sort `tsex' _agedays
	local string "xxx\bmianthro.dta"
	local i=`reflib'
	global bmifile: subinstr local string "xxx" "`i'"
	
	merge `tsex' _agedays using "$bmifile"
	gen double _zbmi=(((_cbmi/m)^l)-1)/(s*l)
	gen double `bsd3p'=m*((1+l*s*3)^(1/l))
	gen double `bsd23p'=`bsd3p'- m*((1+l*s*2)^(1/l))
	replace _zbmi=3+((_cbmi-`bsd3p')/`bsd23p') if (_zbmi>3 & _zbmi~=.)  
	gen double `bsd3n'=m*((1+l*s*(-3))^(1/l))
	gen double `bsd23n'= m*((1+l*s*(-2))^(1/l))-`bsd3n'
	replace _zbmi=(-3)-((`bsd3n'-_cbmi)/`bsd23n') if (_zbmi<-3 & _zbmi~=.) 
	replace _zbmi =. $under5
	keep if _merge~=2
	drop l m s loh _merge
	
	set type float
	* ========================== ZWFL using length ============================

	* ======== interpolate length =======
	* 	determine if length is of 2 decimals significance
	* 	since the LMS charts are only to one, such as 87 point 2
	* 	the `interp' variable =1 if we need to interpolate

	gen `length2'=`lenhei2' if `uselgth'==1
	gen `interp'=0
	gen `temp1'= abs(`length2' - (round(`length2'*10,1))/10)
	replace `interp'=1 if (`temp1' >.001)

	* 	if it does then interpolate on lms table.
	* 	first we will find the lower level on LMS chart.

	gen `lenlow'=int(`length2'*10)/10 if `interp'==1
	replace `lenlow'=round(`length2'*10,1)/10 if `interp'==0

	* 	above makes double sure it will match lms chart.
	* 	next we will find the upper level on LMS chart.

	gen `lenhigh'=round((`length2'+.05)*10,1)/10 if `interp'==1
	replace `lenhigh'=round(`length2'*10,1)/10 if `interp'==0

	* 	we will be doing many calculations with `length2' .
	* 	so we will remember the original values here.

	gen `length9'=`length2'
	
	* 	next we will get the LMS numbers for lenlow & lenhigh

	* ===============lenlow LMS calculations============
	
	replace `length2'=`lenlow'
	sort `tsex' `length2'
	local string "xxx\wflanthro.dta"
	local i=`reflib'
	global wflfile: subinstr local string "xxx" "`i'"
	merge `tsex' `length2' using "$wflfile"

	*	do the ln_hts--> zscore for lenlow where ll=lenlow
	
	gen double `zwhzll'=(((`weight'/m)^l)-1)/(s*l)
	gen double `asd3p'=m*((1+l*s*3)^(1/l))
	gen double `asd23p'=`asd3p'- m*((1+l*s*2)^(1/l))
	replace `zwhzll'=3+((`weight'-`asd3p')/`asd23p') if (`zwhzll'>3 & `zwhzll'~=.)
	gen double `asd3n'=m*((1+l*s*(-3))^(1/l))
	gen double `asd23n'= m*((1+l*s*(-2))^(1/l))-`asd3n'
	replace `zwhzll'=-3-((`asd3n'-`weight')/`asd23n') if (`zwhzll'<-3 & `zwhzll'~=.) 
	keep if _merge~=2
	drop l m s lorh _merge 
	
	
	* ===============lenhigh LMS calculations============
	
	replace `length2'=`lenhigh'
	sort `tsex' `length2'
	local string "xxx\wflanthro.dta"
	local i=`reflib'
	global wflfile: subinstr local string "xxx" "`i'"
	merge `tsex' `length2' using "$wflfile"

	* 	do the ln_hts--> zscore for lenhigh where lh=lenhigh

	gen double `zwhzlh'=(((`weight'/m)^l)-1)/(s*l)
	gen double `csd3p'=m*((1+l*s*3)^(1/l))
	gen double `csd23p'=`csd3p'- m*((1+l*s*2)^(1/l))
	replace `zwhzlh'= 3+((`weight'-`csd3p')/`csd23p') if (`zwhzlh'>3 & `zwhzlh'~=.) 
	gen double `csd3n'=m*((1+l*s*(-3))^(1/l))
	gen double `csd23n'= m*((1+l*s*(-2))^(1/l))-`csd3n'
	replace `zwhzlh'=-3-((`csd3n'-`weight')/`csd23n') if (`zwhzlh'<-3 & `zwhzlh'~=.) 
	keep if _merge~=2
	drop l m s lorh _merge 
	
	* =========== calculations for height used =====.

	* 	use criteria for when height used, height2 is height properly measured

	gen `height2'=`lenhei2' if `uselgth'==2
	
	* 	we will be doing many calculations with height2
	* 	so we will remember the original value here
	
	gen `height9'=`height2'
	
	* 	determine if height is of 2 decimals significance
	* 	since the LMS charts are only to one, such as 87 point 2
	* 	the `interph' variable =1 if we need to interpolate
	
	gen `interph'=0
	gen `temp2'=abs(`height2' - (round(`height2'*10),1)/10)
	replace `interph'=1 if  `temp2' >.001

	* 	if it does then interpolate on lms table
	* 	first we will find the lower level on LMS chart

	gen `hgtlow'=int(`height2'*10)/10 if `interph'==1
	replace `hgtlow'=round(`height2'*10,1)/10 if `interph'==0
	
	
	* 	next we will find the upper level on LMS chart
	
	gen `hgthigh'=round((`height2'+.05)*10,1)/10 if `interph'==1
	replace `hgthigh'=round(`height2'*10,1)/10 if `interph'==0

	* ===================ZWFL using hgtlow===================
	replace `height2'=`hgtlow' 
	sort `tsex' `height2'
	local string "xxx\wfhanthro.dta"
	local i=`reflib'
	global wfhfile: subinstr local string "xxx" "`i'"
	merge `tsex' `height2' using "$wfhfile"
	
	* 	do the heights where hl=hgtlow
	
	gen double `zwhzhl'=(((`weight'/m)^l)-1)/(s*l) if `height2'>0 
	gen double `dsd3p'=m*((1+l*s*3)^(1/l)) if `height2'>0
	gen double `dsd23p'=`dsd3p'- m*((1+l*s*2)^(1/l)) if `height2'>0 
	replace `zwhzhl'= 3+((`weight'-`dsd3p')/`dsd23p') if `height2'>0 & `zwhzhl'>3 & `zwhzhl'~=.
	gen double `dsd3n'=m*((1+l*s*(-3))^(1/l)) if `height2'>0
	gen double `dsd23n'= m*((1+l*s*(-2))^(1/l))-`dsd3n' if `height2'>0
	replace `zwhzhl'=(-3)-((`dsd3n'-`weight')/`dsd23n') if `height2'>0 & `zwhzhl'<-3 & `zwhzhl'~=.
	keep if _merge~=2
	drop l m s lorh _merge
	
	* ===================ZWFL using hgthigh===================
	
	replace `height2'=`hgthigh'
	sort `tsex' `height2'
	local string "xxx\wfhanthro.dta"
	local i=`reflib'
	global wfhfile: subinstr local string "xxx" "`i'"
	merge `tsex' `height2' using "$wfhfile"

	*	do the heights where hh=hgthigh
	
	gen double `zwhzhh'=(((`weight'/m)^l)-1)/(s*l) if `height2'>0 
	gen double `esd3p'=m*((1+l*s*3)^(1/l)) if `height2'>0
	gen double `esd23p'=`esd3p'- m*((1+l*s*2)^(1/l)) if (`height2'>0)
	replace `zwhzhh'=3+((`weight'-`esd3p')/`esd23p') if `height2'>0 & `zwhzhh'>3 & `zwhzhh'~=.
	gen double `esd3n'=m*((1+l*s*(-3))^(1/l)) if `height2'>0
	gen double `esd23n'= m*((1+l*s*(-2))^(1/l))-`esd3n' if `height2'>0
	replace `zwhzhh'=(-3)-((`esd3n'-`weight')/`esd23n') if `height2'>0 & `zwhzhh'<-3 & `zwhzhh'~=.
		
	* =====now do interpolation & choose length or height =====

	* 	length9 is somewhere between lenlow & lenhigh
	* 	find the ratios with #s like 52,20  52,26  & 52,30

	gen `abovel'=`length9'-`lenlow'
	gen `ratiol'=`abovel'/.1

	* 	note that the greater the length, the less the z
	
	gen double _zwfl= `zwhzll'-((`zwhzll'-`zwhzlh')*`ratiol')

	* 	now for height
	* 	height is defined only if uselgth=2 &
	*  	will replace the length calculations if defined
	
	gen `aboveh'=`height9'-`hgtlow'
	gen `ratioh'=`aboveh'/.1
	keep if _merge~=2
	drop l m s lorh _merge 

	* 	note that the greater the height, the less the z

	replace _zwfl=`zwhzhl'-((`zwhzhl'-`zwhzhh')*`ratioh') if `uselgth'==2
		
	
	/**	2nd set additions start here **/
	* ==========================  ZHC  =======================================

	* 	if not done already, need age in days to compare to norms.
	*	select if (agedays2 < (61*30.4375) or missing(agedays2)).

	sort `tsex' _agedays
	local string "xxx\hcanthro.dta"
	local i=`reflib'
	global hczfile: subinstr local string "xxx" "`i'"
	
	merge `tsex' _agedays using "$hczfile"
	gen double _zhc=(((`hc'/m)^l)-1)/(s*l)
	replace _zhc =. $under5
	keep if _merge~=2
	drop l m s _merge 

		* ======================== ZAC ============================================
	sort `tsex' _agedays
	local string "xxx\acanthro.dta"
	local i=`reflib'
	global aczfile: subinstr local string "xxx" "`i'"
	
	
	merge `tsex' _agedays using "$aczfile"
	gen double _zac=(((`armc'/m)^l)-1)/(s*l)
	gen double `msd3p'=m*((1+l*s*3)^(1/l))
	gen double `msd23p'=`msd3p'- m*((1+l*s*2)^(1/l))
	replace _zac=3+((`armc'-`msd3p')/`msd23p') if (_zac>3 & _zac~=.) 
	gen double `msd3n'=m*((1+l*s*(-3))^(1/l))
	gen double `msd23n'= m*((1+l*s*(-2))^(1/l))-`msd3n'
	replace _zac=(-3)-((`msd3n'-`armc')/`msd23n') if (_zac<-3& _zac~=.) 
	replace _zac =. $under5
	keep if _merge~=2
	drop l m s _merge
	
	* ======================== ZTS ============================================
	sort `tsex' _agedays
	local string "xxx\tsanthro.dta"
	local i=`reflib'
	global tszfile: subinstr local string "xxx" "`i'"
	
	
	merge `tsex' _agedays using "$tszfile"
	gen double _zts=(((`tris'/m)^l)-1)/(s*l)
	gen double `tsd3p'=m*((1+l*s*3)^(1/l))
	gen double `tsd23p'=`tsd3p'- m*((1+l*s*2)^(1/l))
	replace _zts=3+((`tris'-`tsd3p')/`tsd23p') if (_zts>3 & _zts~=.) 
	gen double `tsd3n'=m*((1+l*s*(-3))^(1/l))
	gen double `tsd23n'= m*((1+l*s*(-2))^(1/l))-`tsd3n'
	replace _zts=(-3)-((`tsd3n'-`tris')/`tsd23n') if (_zts<-3& _zts~=.) 
	replace _zts =. $under5
	keep if _merge~=2
	drop l m s _merge
	
	* ======================== ZSS ============================================
	sort `tsex' _agedays
	local string "xxx\ssanthro.dta"
	local i=`reflib'
	global sszfile: subinstr local string "xxx" "`i'"
	
	
	merge `tsex' _agedays using "$sszfile"
	gen double _zss=(((`subs'/m)^l)-1)/(s*l)
	gen double `ssd3p'=m*((1+l*s*3)^(1/l))
	gen double `ssd23p'=`ssd3p'- m*((1+l*s*2)^(1/l))
	replace _zss=3+((`subs'-`ssd3p')/`ssd23p') if (_zss>3 & _zss~=.) 
	gen double `ssd3n'=m*((1+l*s*(-3))^(1/l))
	gen double `ssd23n'= m*((1+l*s*(-2))^(1/l))-`ssd3n'
	replace _zss=(-3)-((`ssd3n'-`subs')/`ssd23n') if (_zss<-3& _zss~=.) 
	replace _zss =. $under5
	keep if _merge~=2
	drop l m s _merge

	*===================================================================================
		
	foreach Y in _zwfl _zlen _zwei _zbmi _zhc _zac _zts _zss {
		replace `Y' =round(`Y', 0.01)
	}

	* ====	Set weight-based z-scores to missing for Oedema cases =====
	
	foreach Y in _zwfl _zwei _zbmi {
		replace `Y'=. if `toedema'==2
	}
	
	* ====	Set weight-for-length/height z-scores to missing for children aged 61+ mons =====
	
	foreach Y in _zwfl  {
		replace `Y'=. if _agedays>1856 & _agedays~=.
	}
	
	gen _fwfl=0 if _zwfl~=.
	gen _flen=0 if _zlen~=.
	gen _fwei=0 if _zwei~=.
	gen _fbmi=0 if _zbmi~=.
	gen _fhc=0 if _zhc~=.
	gen _fac=0 if _zac~=.
	gen _fts=0 if _zts~=.
	gen _fss=0 if _zss~=.
	replace _fwfl=1 if (_zwfl< -5 | _zwfl >5) 
	replace _flen=1 if (_zlen< -6 | _zlen >6)
	replace _fwei=1 if (_zwei< -6 | _zwei >5)
	replace _fbmi=1 if (_zbmi< -5 | _zbmi >5)
	replace _fhc=1 if (_zhc< -5 | _zhc >5)
	replace _fac=1 if (_zac< -5 | _zac >5)	
	replace _fts=1 if (_zts< -5 | _zts >5)	
	replace _fss=1 if (_zss< -5 | _zss >5)	
	
	foreach Y in wfl wei bmi len hc ac ts ss {
		replace _f`Y'=. if _z`Y'==.
	}
	
	lab var _zwfl "Weight-for-length/height z-score"
	lab var _zlen "Length/height-for-age z-score"
	lab var _zwei "Weight-for-age z-score"
	lab var _zbmi "BMI-for-age z-score"
	lab var _zhc "Head circumference-for-age z-score"
	lab var _zac "Arm circumference-for-age z-score"
	lab var _zts "Triceps skinfold-for-age z-score"
	lab var _zss "Subscapular skinfold-for-age z-score"

	lab var _fbmi "=1 if (_zbmi < -5 | _zbmi >5)"
	lab var _fwfl "=1 if (_zwfl < -5 | _zwfl >5)"
	lab var _flen "=1 if (_zlen < -6 | _zlen >6)"
	lab var _fwei "=1 if (_zwei < -6 | _zwei >5)"
	lab var _fhc "=1 if (_zhc < -5 | _zhc >5)"
	lab var _fac "=1 if (_zac < -5 | _zac >5)"
	lab var _fts "=1 if (_zts < -5 | _zts >5)"
	lab var _fss "=1 if (_zss < -5 | _zss >5)"

	
	/*	 Clean-up after Phase I = sort data as originally provided	*/
	
	sort `tsort'
	
	} /* quietly ends*/
		
	di "Note: z-scores are flagged according to the following rules:"
	di " "
	di "				_zwfl = . if (_zwfl < -5 or _zwfl >5)"
	di "				_zlen = . if (_zlen < -6 or _zlen >6)"
	di "				_zwei = . if (_zwei < -6 or _zwei >5)"
	di "				_zbmi = . if (_zbmi < -5 or _zbmi >5)"
	di "				_zhc  = . if (_zhc  < -5 or _zhc  >5)"
	di "				_zac  = . if (_zac  < -5 or _zac  >5)"
	di "				_zts  = . if (_zts  < -5 or _zts  >5)"
	di "				_zss  = . if (_zss  < -5 or _zss  >5)"

	di "				"
	
	drop __*
	local string "xxx\yy_z_st.dta"
	local i=`datalib'
	local j=`datalab'
	global outf: subinstr local string "xxx" "`i'" 
	global outfile: subinstr global outf "yy" "`j'"
	
	
	
	
	save "$outfile", replace
	di "Note 1:	Original data plus z-scores are written to"
	di "		$outfile"
	di "				"
	di "				"
	local string "xxx\yy_z_st.xls"
	local i=`datalib'
	local j=`datalab'
	global outs: subinstr local string "xxx" "`i'" 
	global outsheet: subinstr global outs "yy" "`j'"
	outsheet using "$outsheet", replace
	di "Note 2: 	Original data plus z-scores are written to"
	di "		$outsheet"
	
	
	di "				"
	di as txt _n "Please wait, programme is calculating prevalences............."
	di as txt _n ".............................................."
	
	
	/*	Phase II=	Generation of prevalences */
	
	/*	Check the sampling weights before generating prevalences */
	
	qui summ `sw'
	if r(min)<0	{
		di "Error - Negative sampling weights encountered, prevalence tables are not produced!"
	exit
	}
	
	if r(min)>=0 {
			/* quietly begins */
		qui {
	*! version 2	  13jan2019
	
	/// Declaration of relevant files
	
	local string "xxx\yy_z_working_st.dta"
	local i=`datalib'
	local j=`datalab'
	global workf: subinstr local string "xxx" "`i'" 
	global workfile: subinstr global workf "yy" "`j'"
	
	use "$outfile", clear
	save "$workfile", replace
	
	
	/// Generation of age cap (0-59 months)
	gen _age59=0
	replace _age59=1 if _agedays<1826.25
	replace _age59=0 if _agedays==. // 1826.25 is 365.25 (number of days in a year) multiplied by 5
	replace _age59=. if age59==0
	
	gen _age59over=1 if _agedays>=1826/25
	replace _age59over=. if _agedays==.
	/// Generation of age in months
	gen _agemonth=_agedays/365.25*12
	replace _agemonth=. if _age59==0
	
	
	/// Generation of age groups
	gen _agemonth6group=int(_agemonth/6)
	gen _agemonth12group=int(_agemonth/12)
	
	foreach x in _agemonth6group _agemonth12group {
	replace `x'=`x'+1
	}
	
	gen _agemonthgroup=.
	replace _agemonthgroup=1 if _agemonth6group==1 
	replace _agemonthgroup=2 if _agemonth6group==2
	replace _agemonthgroup=3 if _agemonth12group==2
	replace _agemonthgroup=4 if _agemonth12group==3
	replace _agemonthgroup=5 if _agemonth12group==4
	replace _agemonthgroup=6 if _agemonth12group==5
	
	/// Standardized sex	
	
	gen _sex=0
	
	capture confirm string variable sex																// Confirmation loop to confirm sex is string
		if !_rc {	
                       replace _sex=1 if sex=="m" | sex=="M"
					   replace _sex=2 if sex=="f" | sex=="F"
               }
               else {
                       replace _sex=1 if sex==1
					   replace _sex=2 if sex==2
               }
	

	/// Variable Standardization
	
	label variable _age59 "All children under 59 months of age (total sample used for analysis)"
	label variable _agemonth "Age in months"
	label variable _agemonth6group "Age in 6 month groupings"
	label variable _agemonth12group "Age in 12 month groupings"
	label variable _agemonthgroup "Age in 12 month groupings_used for analysis"
	
	/// Value Label Standardization
	
	label define _age59label 1 "Total"
	label define _sexlabel 1 "Male" 2 "Female"
	label define _agemonth6grouplabel 1 "0 to 5 months" 2 "6 to 11 months" 3 "12 to 17 months" 4 "18 to 23 months" 5 "24 to 29 months" 6 "30 to 35 months" 7 "36 to 41 months" 8 "42 to 47 months" 9 "48 to 53 months" 10 "54 to 59 months"
	label define _agemonth12grouplabel 1 "0 to 11 months" 2 "12 to 23 months" 3 "24 to 35 months" 4 "36 to 47 months" 5 "48 to 59 months"	
	label define _agemonthgrouplabel 1 "0 to 5 months" 2 "6 to 11 months" 3 "12 to 23 months" 4 "24 to 35 months" 5 "36 to 47 months" 6 "48 to 59 months"			
	
	/// Apply Value labels to variables
	
	foreach x in _age59 _sex _agemonthgroup {
	label values `x' `x'label
	}
	
	/// Generate Age Sex Disaggregation
	
	foreach x in sex {
		foreach y in _agemonthgroup {
			egen `y'`x'= group(`x' `y'),label
				}
	
	/// Generation of Standardized Anthro Indicators
	
		gen haz=_zlen				// height-for-age
		gen waz=_zwei				// weight-for-age
		gen whz=_zwfl				// weight-for-height	
		gen baz=_zbmi				// bmi-for-age
		
		gen flag_haz=_flen
		gen flag_waz=_fwei
		gen flag_whz=_fwfl
		gen flag_baz=_fbmi
		
		foreach x in haz waz whz baz {
		replace `x'=. if _age59over==1
		}
		
		
	/// Generation of categorical anthropometric variables
	/// A <-3, B <-2, C <-2, D >+1, E >+2, F >+3
		
		foreach x in haz waz whz baz {
		foreach y in _3 _2 _1 1 2 3 {
		gen `x'`y'=0 
		}	
		}
		
		foreach x in haz waz whz baz {
		replace `x'_3=1 if `x'<-3 
		replace `x'_2=1 if `x'<-2 
		replace `x'_1=1 if `x'<-1 
		replace `x'1=1 if `x'>1 
		replace `x'2=1 if `x'>2 
		replace `x'3=1 if `x'>3 
		}	
	
	/// Removal of biologically implausible, children that fall out of the scope
	/// Include children with oedema in the denominator (indicator dependent)
		
		foreach x in haz {
		foreach y in _3 _2 _1 1 2 3 {
		replace `x'`y'=. if `x'==.
		replace `x'`y'=. if `x'>6.00 |`x'<-6.00
		replace `x'=. if `x'>6.00 | `x'<-6.00
		replace `x'`y'=. if _agemonth>60 | _agemonth==60
		}
		}
		
		foreach x in waz  {
		foreach y in _3 _2 _1 1 2 3 {
		replace `x'`y'=. if `x'==.
		replace `x'`y'=. if `x'>5.00 |`x'<-6.00
		replace `x'=. if `x'>5.00 | `x'<-6.00 
		replace `x'`y'=. if _agemonth>60 | _agemonth==60
		replace `x'`y'=0 if oedema=="y" | oedema=="Y"
		}
		}	
	
		foreach x in baz {
		foreach y in _3 _2 _1 1 2 3 {
		replace `x'`y'=. if `x'==.
		replace `x'`y'=. if `x'>5.00 |`x'<-5.00
		replace `x'=. if `x'>5.00 | `x'<-5.00 
		replace `x'`y'=. if _agemonth>60 | _agemonth==60
		replace `x'`y'=0 if oedema=="y" | oedema=="Y"
		}
		}		
		
		foreach x in whz {
		foreach y in _3 _2 _1 1 2 3 {
		replace `x'`y'=. if `x'==.
		replace `x'`y'=. if `x'>5.00 |`x'<-5.00
		replace `x'=. if `x'>5.00 | `x'<-5.00 
		replace `x'`y'=0 if oedema=="y" | oedema=="Y"
		}
		}		
		
		
	/// Include children with oedema in the numerator (indicator dependent)
	
		foreach x in baz waz whz {
			foreach z in _3 _2 _1 {
				replace `x'`z'=1 if oedema=="y" | oedema=="Y"
				}
			}	
	
	
		
	
	
	}
	
end
exit
