/*************************************************************************************************************/
/********  	WHO Child Growth Standards                                                         	  ****/
/********  	Department of Nutrition for Health and Development                                 	  ****/
/********  	World Health Organization                                                          	  ****/
/********  	Last modified on 23/05/2007     -     For STATA versions 7.0 and above             	  ****/
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
	
	tempvar agetemp hcw1 noage agegrp
	
	gen `agetemp'=`age'
	replace `agetemp'=`age'*30.4375 if `ageunit'=="months"
	gen `hcw1'=int(`agetemp'/30.4375)
	gen `noage'=1 if (`hcw1'==.  & `age'==.)
	gen `agegrp'=`hcw1'
	recode `agegrp' 0/5=1 6/11=2 12/23=3 24/35=4 36/47=5 48/60=6 61/max=.
	replace `agegrp'=0 if `noage'==1
	lab def agegrp  0 "no age" 1 "0-5" 2 "6-11" 3 "12-23" 4 "24-35" 5 "36-47" 6 "48-60", modify 
	lab val `agegrp' agegrp
	
	/*	Here I needed to declare the temporary variables for sex, oedema and ovflag again */
	
	tempvar tsex toedema ovflag
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

	qui tab `oedema'
	if  r(r)==0 {
		gen `toedema'=1  if `oedema'==.
	}
	if  r(r)~=0 {
		gen `toedema'= 1 
		replace `toedema'= 2 if `oedema'=="y" | `oedema'=="Y"
	}
	tempvar ovflag _fwfl _flen _fwei _fbmi _fhc _fac _fss _fts
	
	gen `_fwfl'=0 if _zwfl~=.
	gen `_flen'=0 if _zlen~=.
	gen `_fwei'=0 if _zwei~=.
	gen `_fbmi'=0 if _zbmi~=.
	gen `_fhc'=0 if _zhc~=.
	gen `_fac'=0 if _zac~=.
	gen `_fts'=0 if _zts~=.
	gen `_fss'=0 if _zss~=.
		
	replace `_fwfl'=1 if (_zwfl< -5 | _zwfl >5) 
	replace `_flen'=1 if (_zlen< -6 | _zlen >6)
	replace `_fwei'=1 if (_zwei< -6 | _zwei >5)
	replace `_fbmi'=1 if (_zbmi< -5 | _zbmi >5)
	replace `_fhc'=1 if (_zhc< -5 | _zhc >5)
	replace `_fac'=1 if (_zac< -5 | _zac >5)
	replace `_fss'=1 if (_zss< -5 | _zss >5)
	replace `_fts'=1 if (_zts< -5 | _zts >5)
	gen `ovflag'=((`_fwfl'==1) | (`_flen'==1) | (`_fwei'==1) | (`_fbmi'==1))
	
	foreach X in _zwfl _zlen _zwei _zbmi _zhc _zac _zss _zts {
		tempvar `X'p1 `X'p2 `X'p3 `X'm2 `X'm3
	}
		
	/*	For each indicator, declare binary variables to calculate prevalences*/
		*	p1= 1 if zscore above 1SD,  0 otherwise
		*	p2= 1 if zscore above 2SD,  0 otherwise
		*	p3= 1 if zscore above 3SD,  0 otherwise
		*	m2= 1 if zscore below -2SD, 0 otherwise
		*	m3= 1 if zscore below -3SD, 0 otherwise
	
	foreach Y in _zlen _zhc _zac _zts _zss {
		gen ``Y'p1'= `Y' > 1 & `Y' ~=.
		gen ``Y'p2'= `Y' > 2 & `Y' ~=.
		gen ``Y'p3'= `Y' > 3 & `Y' ~=.
		gen ``Y'm2'=  (`Y' < -2 & `Y' ~=.) 
		gen ``Y'm3'=  (`Y' < -3 & `Y' ~=.) 
	}
	
	/*	Special adjustment for the weight-based indicators with oedema for m2 and m3*/
	
	foreach Y in _zwfl _zbmi _zwei  {
		gen ``Y'p1'= `Y' > 1 & `Y' ~=.
		gen ``Y'p2'= `Y' > 2 & `Y' ~=.
		gen ``Y'p3'= `Y' > 3 & `Y' ~=.
		gen ``Y'm2'=  (`Y' < -2 & `Y' ~=.) 
		replace ``Y'm2'= 1 if `toedema'==2 
		gen ``Y'm3'=  (`Y' < -3 & `Y' ~=.) 	
		replace ``Y'm3'= 1 if `toedema'==2 
	}
		
	set type double	
	
		/*	declare flagged values missing by indicator*/

		foreach Y in wfl bmi wei len hc ac ts ss {
			replace _z`Y' =. if _f`Y'==1 
		}
		
		/*	Sexes combined */
		
		/*	Start with the zscores for weight-for-age	*/		
		
		foreach X in _zwei {
		
			/*	Calculates Means and SDs - for total agegroups */
		
			tempvar temp swt y x yx swtsd ysd yxsd
			gen `temp'=1 if `X'~=.				/* temp=1 is a tag for all children with non-missing _zwei */
			gen `swt'=`X'*`sw' if `temp'==1 		/* swt=_zwei * sampling weight for all children if their temp=1 */
			egen `y'=sum(`swt') if `temp'==1		/* y=add-up swt for all children if their temp=1 */
			egen `x'=sum(`sw')  if `temp'==1		/* x=add-up sampling weight of all children if their temp=1 */
			gen `yx'=(`y'/`x')				/* yx=weighted mean of _zwei= y/x	*/
			gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1  	/* swtsd= sum of squared (_zwei - mean) times sampling weight */
			egen `ysd'=sum(`swtsd') if `temp'==1		/* ysd=add-up swtsd of all children if their temp=1 */
			gen `yxsd'=(`ysd'/(`x'-1))^0.5			/* yxsd=standard deviation=ysd divided by x-1*/
			summ `yx'
			if r(N)~=0 {
				matrix mean0= round(r(mean), 0.01)	/* put the mean _zwei in matrix mean0*/
			}
			if r(N)==0 {
				matrix mean0= -69
			}
			summ `yxsd'					/* put the SD _zwei in matrix sd0*/
			if r(N)~=0 {
				matrix sd0= round(r(mean), 0.01)
			}
			if r(N)==0 {
				matrix sd0= -69
			}		
			/*	Calculates Means and SDs - for disaggregated agegroups (same logic as above)*/
			
			foreach Z in 1 2 3 4 5 6 {
				tempvar temp swt y x yx swtsd ysd yxsd
				gen `temp'=1 if `X'~=. & `agegrp'==`Z'
				gen `swt'=`X'*`sw' if `temp'==1 
				egen `y'=sum(`swt') if `temp'==1
				egen `x'=sum(`sw')  if `temp'==1
				gen `yx'=(`y'/`x')
				gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1 
				egen `ysd'=sum(`swtsd') if `temp'==1
				gen `yxsd'=(`ysd'/(`x'-1))^0.5
				summ `yx'
				if r(N)~=0 {
					matrix mean`Z'= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix mean`Z'= -69
				}
				summ `yxsd'
				if r(N)~=0 {
					matrix sd`Z'= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix sd`Z'= -69
				}		
			}
			
			/*	Collate row matrices of total and disaggregate ages
				into matrix mean and matrix sd of weight-for-age zscore	*/
			
			foreach Q in mean sd  {
				matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6
			}
			
			/*	Calculation of prevalences starts here */
			
			foreach Y in m3 m2  {
				
				/*	Calculates prevalences for total agegroups */
				
				tempvar temp swt y x yx swt
				gen `temp'=1 if `X'~=.				 /* temp=1 is a tag for all children with non-missing _zwei */
				recode `temp' .=1 if `toedema'==2 & `agegrp'~=.  /* also temp=1 if children have oedema & younger than 61 months*/
				
				/*	in the "recode" above, is a special adjustment for oedema cases which applies to
					 _zwei _zwfl _zbmi*/
					 
				/*	This special adjustment will be absent when flagsys==1; restricted analysis i.e. listwise deletion*/
				
				/*	Take for example Y=m3, i.e. the prevalence of _zwei below -3SD	*/
				
				gen `swt'=``X'`Y''*`sw' if `temp'==1 		/* swt=_zweim3 times sampling weight for all children if their temp=1*/
				egen `y'=sum(`swt') if `temp'==1		/* y=add-up swt for all children if their temp=1*/
				egen `x'=sum(`sw') if `temp'==1		/* x=add-up sampling weight of all children if their temp=1 */	
				gen `yx'=(`y'/`x')				/* yx=required prevalence of _zwei below -3SD	*/	
				summ `x'					/* x=Weighted N = sum of sampling weights of all children if their temp=1*/
				if r(N)~=0 {
					matrix A0=r(mean)			/* put the weighted N in matrix A0*/
				}
				if r(N)==0 {
					matrix A0=0
				}
				summ `yx'
				if r(N)~=0 {
					matrix B0=r(mean)			/* put the proportion yx in matrix B0*/
				}
				if r(N)==0 {
					matrix B0= -69
				}
				
				/*	Calculates prevalences - for disaggregated agegroups (same logic as above)*/
				
				foreach Z in 1 2 3 4 5 6 {
					tempvar temp swt y x yx swt
					gen `temp'=1 if `X'~=. & `agegrp'==`Z'
					recode `temp' .=1 if `toedema'==2 & `agegrp'==`Z'
					gen `swt'=``X'`Y''*`sw' if `temp'==1 
					egen `y'=sum(`swt') if `temp'==1
					egen `x'=sum(`sw') if `temp'==1
					gen `yx'=(`y'/`x')
					summ `x'
					if r(N)~=0 {
						matrix A`Z'=r(mean)
					}
					if r(N)==0 {
						matrix A`Z'=0
					}
					summ `yx'
					if r(N)~=0 {
						matrix B`Z'=r(mean)
					}
					if r(N)==0 {
						matrix B`Z'= -69
					}
				}
				/* 	For total (Z=0) and disaggregate ages (Z= 1 to 6) 
					matrix A holds the weighted N 
					matrix B holds the proportion or prevalence
					matrix X holds lower 95% CI
					matrix Y holds upper 95% CI	
				*/					
				
				foreach Z in 0 1 2 3 4 5 6 {
					if A`Z'[1,1] ~=0 {
						matrix ga`Z'= inv(A`Z')
						matrix gb`Z'= inv(A`Z'*2)
						matrix D`Z' = 1-B`Z'[1,1]
						matrix F`Z' = 1.96*(B`Z'[1,1]*D`Z'[1,1]*ga`Z'[1,1])^0.5 + gb`Z'[1,1]
						matrix X`Z' = max(B`Z'[1,1] - F`Z'[1,1],0) 
						matrix Y`Z' = B`Z'[1,1] + F`Z'[1,1]
						matrix A`Z' = round(A`Z'[1,1],1)
						matrix B`Z' = round(B`Z'[1,1]*100, 0.1)
						matrix X`Z' = round(X`Z'[1,1]*100, 0.1)
						matrix Y`Z' = round(Y`Z'[1,1]*100, 0.1)
					}
					if A`Z'[1,1] ==0 {
						matrix X`Z' = -69
						matrix Y`Z' = -69
					}
				}
				
				/*	Collates  row matrices of total and disaggregate ages */
				
				foreach Q in A B X Y  {
					matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6
				}
				
				/*	Collates column matrices of prevalence (B) and lower 95% (X) and higher 95% (Y)*/
				matrix `Y'=B, X, Y
			}
			
			/*	Collating together the full matrix for _zwei
				A = Weighted N (which is the same whether the prevalence is m2 or m3)
				m3 = prevalence m3 with 95%C.I.
				m2 = prevalence m2 with 95%C.I.
				mean = mean zscore weight-for-age (_zwei)
				sd =   SD zscore weight-for-age (_zwei)
			*/				
			matrix `X'= A, m3, m2, mean, sd
			matrix colnames `X'=" N" "%<-3SD" "  95%" "  CI " "%<-2SD" "  95%" "  CI "  "Mean" "SD" 
			matrix rownames `X'="Total" "0-5" "6-11" "12-23" "24-35" "36-47" "48-60"
		}
		foreach X in _zlen {
			tempvar temp swt y x yx swtsd ysd yxsd
			gen `temp'=1 if `X'~=. 
			gen `swt'=`X'*`sw' if `temp'==1 
			egen `y'=sum(`swt') if `temp'==1
			egen `x'=sum(`sw')  if `temp'==1
			gen `yx'=(`y'/`x')
			gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1 
			egen `ysd'=sum(`swtsd') if `temp'==1
			gen `yxsd'=(`ysd'/(`x'-1))^0.5
			summ `yx'
			if r(N)~=0 {
				matrix mean0= round(r(mean), 0.01)
			}
			if r(N)==0 {
				matrix mean0= -69
			}
			summ `yxsd'
			if r(N)~=0 {
				matrix sd0= round(r(mean), 0.01)
			}
			if r(N)==0 {
				matrix sd0= -69
			}
			foreach Z in 1 2 3 4 5 6 {
				tempvar temp swt y x yx swtsd ysd yxsd
				gen `temp'=1 if `X'~=. & `agegrp'==`Z'
				gen `swt'=`X'*`sw' if `temp'==1 
				egen `y'=sum(`swt') if `temp'==1
				egen `x'=sum(`sw') if `temp'==1
				gen `yx'=(`y'/`x')
				gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1 
				egen `ysd'=sum(`swtsd') if `temp'==1
				gen `yxsd'=(`ysd'/(`x'-1))^0.5
				summ `yx'
				if r(N)~=0 {
					matrix mean`Z'= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix mean`Z'= -69
				}
				summ `yxsd'
				if r(N)~=0 {
					matrix sd`Z'= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix sd`Z'= -69
				}
			}
		
			foreach Q in mean sd  {
				matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6
			}
		
			foreach Y in m3 m2  {
				tempvar temp swt y x yx swt
				gen `temp'=1 if `X'~=.
				gen `swt'=``X'`Y''*`sw' if `temp'==1 
				egen `y'=sum(`swt') if `temp'==1
				egen `x'=sum(`sw') if `temp'==1
				gen `yx'=(`y'/`x')
				summ `x'
				if r(N)~=0 {
					matrix A0=r(mean)
				}
				if r(N)==0 {
					matrix A0=0
				}
				summ `yx'
				if r(N)~=0 {
					matrix B0=r(mean)
				}
				if r(N)==0 {
					matrix B0= -69
				}
				foreach Z in 1 2 3 4 5 6 {
					tempvar temp swt y x yx swt
					gen `temp'=1 if `X'~=. & `agegrp'==`Z'
					gen `swt'=``X'`Y''*`sw' if `temp'==1 
					egen `y'=sum(`swt') if `temp'==1
					egen `x'=sum(`sw') if `temp'==1
					gen `yx'=(`y'/`x')
					summ `x'
					if r(N)~=0 {
						matrix A`Z'=r(mean)
					}
					if r(N)==0 {
						matrix A`Z'=0
					}
					summ `yx'
					if r(N)~=0 {
						matrix B`Z'=r(mean)
					}
					if r(N)==0 {
						matrix B`Z'= -69
					}
				}
				foreach Z in 0 1 2 3 4 5 6 {
					if A`Z'[1,1] ~=0 {
						matrix ga`Z'= inv(A`Z')
						matrix gb`Z'= inv(A`Z'*2)
						matrix D`Z' = 1-B`Z'[1,1]
						matrix F`Z' = 1.96*(B`Z'[1,1]*D`Z'[1,1]*ga`Z'[1,1])^0.5 + gb`Z'[1,1]
						matrix X`Z' = max(B`Z'[1,1] - F`Z'[1,1],0) 
						matrix Y`Z' = B`Z'[1,1] + F`Z'[1,1]
						matrix A`Z' = round(A`Z'[1,1],1)
						matrix B`Z' = round(B`Z'[1,1]*100, 0.1)
						matrix X`Z' = round(X`Z'[1,1]*100, 0.1)
						matrix Y`Z' = round(Y`Z'[1,1]*100, 0.1)
					}
					if A`Z'[1,1] ==0 {
						matrix X`Z' = -69
						matrix Y`Z' = -69
					}
				}
		
				foreach Q in A B X Y  {
					matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6
				}
				matrix `Y'=B, X, Y
			}
			matrix `X'= A, m3, m2, mean, sd
			matrix colnames `X'=" N" "%<-3SD" "  95%" "  CI " "%<-2SD" "  95%" "  CI "  "Mean" "SD" 
			matrix rownames `X'="Total" "0-5" "6-11" "12-23" "24-35" "36-47" "48-60"
		}
		
		foreach X in _zwfl _zbmi {
			tempvar temp swt y x yx swtsd ysd yxsd
			gen `temp'=1 if `X'~=.
			gen `swt'=`X'*`sw' if `temp'==1 
			egen `y'=sum(`swt') if `temp'==1
			egen `x'=sum(`sw') if `temp'==1
			gen `yx'=(`y'/`x')
			gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1 
			egen `ysd'=sum(`swtsd') if `temp'==1
			gen `yxsd'=(`ysd'/(`x'-1))^0.5
			summ `yx'
			if r(N)~=0 {
				matrix mean0= round(r(mean), 0.01)
			}
			if r(N)==0 {
				matrix mean0= -69
			}
			summ `yxsd'
			if r(N)~=0 {
				matrix sd0= round(r(mean), 0.01)
			}
			if r(N)==0 {
				matrix sd0= -69
			}
			foreach Z in 1 2 3 4 5 6 {
				tempvar temp swt y x yx swtsd ysd yxsd
				gen `temp'=1 if `X'~=. & `agegrp'==`Z'
				gen `swt'=`X'*`sw' if `temp'==1 
				egen `y'=sum(`swt') if `temp'==1
				egen `x'=sum(`sw') if `temp'==1
				gen `yx'=(`y'/`x')
				gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1 
				egen `ysd'=sum(`swtsd') if `temp'==1
				gen `yxsd'=(`ysd'/(`x'-1))^0.5
				summ `yx'
				if r(N)~=0 {
					matrix mean`Z'= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix mean`Z'= -69
				}
				summ `yxsd'
				if r(N)~=0 {
					matrix sd`Z'= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix sd`Z'= -69
				}
			}
		
			foreach Q in mean sd  {
				matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6
			}
			foreach Y in m3 m2 p1 p2 p3 {
				tempvar temp swt y x yx swt
				gen `temp'=1 if `X'~=. 
				recode `temp' .=1 if `toedema'==2 & `agegrp'~=.
				gen `swt'=``X'`Y''*`sw' if `temp'==1 
				egen `y'=sum(`swt') if `temp'==1
				egen `x'=sum(`sw') if `temp'==1
				gen `yx'=(`y'/`x')
				summ `x'
				if r(N)~=0 {
					matrix A0=r(mean)
				}
				if r(N)==0 {
					matrix A0=0
				}
				summ `yx'
				if r(N)~=0 {
					matrix B0=r(mean)
				}
				if r(N)==0 {
					matrix B0= -69
				}
				foreach Z in 1 2 3 4 5 6 {
					tempvar temp swt y x yx swt
					gen `temp'=1 if `X'~=. & `agegrp'==`Z'
					recode `temp' .=1 if `toedema'==2 & `agegrp'==`Z'
					gen `swt'=``X'`Y''*`sw' if `temp'==1 
					egen `y'=sum(`swt') if `temp'==1
					egen `x'=sum(`sw') if `temp'==1
					gen `yx'=(`y'/`x')
					summ `x'
					if r(N)~=0 {
						matrix A`Z'=r(mean)
					}
					if r(N)==0 {
						matrix A`Z'=0
					}
					summ `yx'
					if r(N)~=0 {
						matrix B`Z'=r(mean)
					}
					if r(N)==0 {
						matrix B`Z'= -69
					}
				}
				foreach Z in 0 1 2 3 4 5 6 {
					if A`Z'[1,1] ~=0 {
						matrix ga`Z'= inv(A`Z')
						matrix gb`Z'= inv(A`Z'*2)
						matrix D`Z' = 1-B`Z'[1,1]
						matrix F`Z' = 1.96*(B`Z'[1,1]*D`Z'[1,1]*ga`Z'[1,1])^0.5 + gb`Z'[1,1]
						matrix X`Z' = max(B`Z'[1,1] - F`Z'[1,1],0) 
						matrix Y`Z' = B`Z'[1,1] + F`Z'[1,1]
						matrix A`Z' = round(A`Z'[1,1],1)
						matrix B`Z' = round(B`Z'[1,1]*100, 0.1)
						matrix X`Z' = round(X`Z'[1,1]*100, 0.1)
						matrix Y`Z' = round(Y`Z'[1,1]*100, 0.1)
					}
					if A`Z'[1,1] ==0 {
						matrix X`Z' = -69
						matrix Y`Z' = -69
					}
				}
		
				foreach Q in A B X Y  {
					matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6
				}
				matrix `Y'=B, X, Y
			}
			matrix `X'= A, m3, m2, p1,p2, p3, mean, sd
			matrix colnames `X'=" N" "%<-3SD" "  95%" "  CI " "%<-2SD" "  95%" "  CI " "%>1SD" "  95%" "  CI " "%>2SD" "  95%" "  CI " "%>3SD" "  95%" "  CI " "Mean" "SD"
			matrix rownames `X'="Total" "0-5" "6-11" "12-23" "24-35" "36-47" "48-60"
		}
		foreach X in _zhc _zac _zts _zss {
			tempvar temp swt y x yx swtsd ysd yxsd a
			gen `a'= match("`X'","_zhc")
			gen `temp'=1 if `X'~=. 
			gen `swt'=`X'*`sw' if `temp'==1 
			egen `y'=sum(`swt') if `temp'==1
			egen `x'=sum(`sw')  if `temp'==1
			gen `yx'=(`y'/`x')
			gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1 
			egen `ysd'=sum(`swtsd') if `temp'==1
			gen `yxsd'=(`ysd'/(`x'-1))^0.5
			summ `yx'
			if r(N)~=0 {
				matrix mean0= round(r(mean), 0.01)
			}
			if r(N)==0 {
				matrix mean0= -69
			}
			summ `yxsd'
			if r(N)~=0 {
				matrix sd0= round(r(mean), 0.01)
			}
			if r(N)==0 {
				matrix sd0= -69
			}
			foreach Z in 1 2 3 4 5 6 {
				tempvar temp swt y x yx swtsd ysd yxsd
				gen `temp'=1 if `X'~=. & `agegrp'==`Z'
				gen `swt'=`X'*`sw' if `temp'==1 
				egen `y'=sum(`swt') if `temp'==1
				egen `x'=sum(`sw') if `temp'==1
				gen `yx'=(`y'/`x')
				gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1 
				egen `ysd'=sum(`swtsd') if `temp'==1
				gen `yxsd'=(`ysd'/(`x'-1))^0.5
				summ `yx'
				if r(N)~=0 {
					matrix mean`Z'= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix mean`Z'= -69
				}
				summ `yxsd'
				if r(N)~=0 {
					matrix sd`Z'= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix sd`Z'= -69
				}
			}
		
			foreach Q in mean sd  {
				matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6
			}
		
			foreach Y in m3 m2 p1 p2 p3 {
				tempvar temp swt y x yx swt
				gen `temp'=1 if `X'~=.
				gen `swt'=``X'`Y''*`sw' if `temp'==1 
				egen `y'=sum(`swt') if `temp'==1
				egen `x'=sum(`sw') if `temp'==1
				gen `yx'=(`y'/`x')
				summ `x'
				if r(N)~=0 {
					matrix A0=r(mean)
				}
				if r(N)==0 {
					matrix A0=0
				}
				summ `yx'
				if r(N)~=0 {
					matrix B0=r(mean)
				}
				if r(N)==0 {
					matrix B0= -69
				}
				foreach Z in 1 2 3 4 5 6 {
					tempvar temp swt y x yx swt
					gen `temp'=1 if `X'~=. & `agegrp'==`Z'
					gen `swt'=``X'`Y''*`sw' if `temp'==1 
					egen `y'=sum(`swt') if `temp'==1
					egen `x'=sum(`sw') if `temp'==1
					gen `yx'=(`y'/`x')
					summ `x'
					if r(N)~=0 {
						matrix A`Z'=r(mean)
					}
					if r(N)==0 {
						matrix A`Z'=0
					}
					summ `yx'
					if r(N)~=0 {
						matrix B`Z'=r(mean)
					}
					if r(N)==0 {
						matrix B`Z'= -69
					}
				}
				foreach Z in 0 1 2 3 4 5 6 {
					if A`Z'[1,1] ~=0 {
						matrix ga`Z'= inv(A`Z')
						matrix gb`Z'= inv(A`Z'*2)
						matrix D`Z' = 1-B`Z'[1,1]
						matrix F`Z' = 1.96*(B`Z'[1,1]*D`Z'[1,1]*ga`Z'[1,1])^0.5 + gb`Z'[1,1]
						matrix X`Z' = max(B`Z'[1,1] - F`Z'[1,1],0) 
						matrix Y`Z' = B`Z'[1,1] + F`Z'[1,1]
						matrix A`Z' = round(A`Z'[1,1],1)
						matrix B`Z' = round(B`Z'[1,1]*100, 0.1)
						matrix X`Z' = round(X`Z'[1,1]*100, 0.1)
						matrix Y`Z' = round(Y`Z'[1,1]*100, 0.1)
					}
					if A`Z'[1,1] ==0 {
						matrix X`Z' = -69
						matrix Y`Z' = -69
					}
				}
		
				foreach Q in A B X Y  {
					matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6
				}
				matrix `Y'=B, X, Y
			}
			matrix `X'= A, m3, m2, p1,p2, p3, mean, sd
			matrix colnames `X'=" N" "%<-3SD" "  95%" "  CI " "%<-2SD" "  95%" "  CI " "%>1SD" "  95%" "  CI " "%>2SD" "  95%" "  CI " "%>3SD" "  95%" "  CI " "Mean" "SD"
			matrix rownames `X'="Total" "3-5" "6-11" "12-23" "24-35" "36-47" "48-60"
			if `a'==1 {
				matrix rownames `X'="Total" "0-5" "6-11" "12-23" "24-35" "36-47" "48-60"
			}

		}
		
		/* Sex=1  Males;	 Sex=2 Females */
		foreach S in 1 2  {
			foreach X in _zwei {
				tempvar temp swt y x yx swtsd ysd yxsd
				gen `temp'=1 if `X'~=. & `tsex' ==`S'
				gen `swt'=`X'*`sw' if `temp'==1 
				egen `y'=sum(`swt') if `temp'==1 
				egen `x'=sum(`sw') if `temp'==1 
				gen `yx'=(`y'/`x')
				gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1 
				egen `ysd'=sum(`swtsd') if `temp'==1 
				gen `yxsd'=(`ysd'/(`x'-1))^0.5
				summ `yx'
				if r(N)~=0 {
					matrix mean0= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix mean0= -69
				}
				summ `yxsd'
				if r(N)~=0 {
					matrix sd0= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix sd0= -69
				}
				foreach Z in 1 2 3 4 5 6 {
					tempvar temp swt y x yx swtsd ysd yxsd
					gen `temp'=1 if `X'~=. & `agegrp'==`Z' & `tsex' ==`S'
					gen `swt'=`X'*`sw' if `temp'==1 
					egen `y'=sum(`swt') if `temp'==1 
					egen `x'=sum(`sw') if `temp'==1 
					gen `yx'=(`y'/`x')
					gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1 
			 		egen `ysd'=sum(`swtsd') if `temp'==1 
					gen `yxsd'=(`ysd'/(`x'-1))^0.5
					summ `yx'
					if r(N)~=0 {
						matrix mean`Z'= round(r(mean), 0.01)
					}
					if r(N)==0 {
						matrix mean`Z'= -69
					}
					summ `yxsd'
					if r(N)~=0 {
						matrix sd`Z'= round(r(mean), 0.01)
					}
					if r(N)==0 {
						matrix sd`Z'= -69
					}
				}
		
				foreach Q in mean sd  {
					matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6
				}
				foreach Y in m3 m2  {
					tempvar temp swt y x yx swt
					gen `temp'=1 if `X'~=. & `tsex' ==`S'
					recode `temp' .=1 if `toedema'==2 & `agegrp'~=. & `tsex' ==`S'
					gen `swt'=``X'`Y''*`sw' if `temp'==1 
					egen `y'=sum(`swt') if `temp'==1 
					egen `x'=sum(`sw') if `temp'==1 
					gen `yx'=(`y'/`x')
					summ `x'
					if r(N)~=0 {
						matrix A0=r(mean)
					}
					if r(N)==0 {
						matrix A0=0
					}
					summ `yx'
					if r(N)~=0 {
						matrix B0=r(mean)
					}
					if r(N)==0 {
						matrix B0= -69
					}
					foreach Z in 1 2 3 4 5 6 {
						tempvar temp swt y x yx swt
						gen `temp'=1 if `X'~=. & `agegrp'==`Z' & `tsex' ==`S'
						recode `temp' .=1 if `toedema'==2 & `agegrp'==`Z' & `tsex' ==`S'
						gen `swt'=``X'`Y''*`sw' if `temp'==1 
						egen `y'=sum(`swt') if `temp'==1 
						egen `x'=sum(`sw') if `temp'==1 
						gen `yx'=(`y'/`x')
						summ `x'
						if r(N)~=0 {
							matrix A`Z'=r(mean)
						}
						if r(N)==0 {
							matrix A`Z'=0
						}
						summ `yx'
						if r(N)~=0 {
							matrix B`Z'=r(mean)
						}
						if r(N)==0 {
							matrix B`Z'= -69
						}
					}
					foreach Z in 0 1 2 3 4 5 6 {
						if A`Z'[1,1] ~=0 {
							matrix ga`Z'= inv(A`Z')
							matrix gb`Z'= inv(A`Z'*2)
							matrix D`Z' = 1-B`Z'[1,1]
							matrix F`Z' = 1.96*(B`Z'[1,1]*D`Z'[1,1]*ga`Z'[1,1])^0.5 + gb`Z'[1,1]
							matrix X`Z' = max(B`Z'[1,1] - F`Z'[1,1],0) 
							matrix Y`Z' = B`Z'[1,1] + F`Z'[1,1]
							matrix A`Z' = round(A`Z'[1,1],1)
							matrix B`Z' = round(B`Z'[1,1]*100, 0.1)
							matrix X`Z' = round(X`Z'[1,1]*100, 0.1)
							matrix Y`Z' = round(Y`Z'[1,1]*100, 0.1)
						}
						if A`Z'[1,1] ==0 {
							matrix X`Z' = -69
							matrix Y`Z' = -69
						}				
					}
					foreach Q in A B X Y  {
						matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6
					}
					matrix `Y'=B, X, Y
				}
				matrix `X'`S'= A, m3, m2, mean, sd
				matrix colnames `X'`S'=" N" "%<-3SD" "  95%" "  CI " "%<-2SD" "  95%" "  CI "  "Mean" "SD" 
				matrix rownames `X'`S'="Total" "0-5" "6-11" "12-23" "24-35" "36-47" "48-60"
			}
			foreach X in _zlen {
				tempvar temp swt y x yx swtsd ysd yxsd
				gen `temp'=1 if `X'~=. & `tsex' ==`S'
				gen `swt'=`X'*`sw' if `temp'==1 
				egen `y'=sum(`swt') if `temp'==1
				egen `x'=sum(`sw') if `temp'==1
				gen `yx'=(`y'/`x')
				gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1
				egen `ysd'=sum(`swtsd') if `temp'==1
				gen `yxsd'=(`ysd'/(`x'-1))^0.5
				summ `yx'
				if r(N)~=0 {
					matrix mean0= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix mean0= -69
				}
				summ `yxsd'
				if r(N)~=0 {
					matrix sd0= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix sd0= -69
				}
				foreach Z in 1 2 3 4 5 6 {
					tempvar temp swt y x yx swtsd ysd yxsd
					gen `temp'=1 if `X'~=. & `agegrp'==`Z' & `tsex' ==`S'
					gen `swt'=`X'*`sw' if `temp'==1
					egen `y'=sum(`swt') if `temp'==1
					egen `x'=sum(`sw') if `temp'==1
					gen `yx'=(`y'/`x')
					gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1
					egen `ysd'=sum(`swtsd') if `temp'==1
					gen `yxsd'=(`ysd'/(`x'-1))^0.5
					summ `yx'
					if r(N)~=0 {
						matrix mean`Z'= round(r(mean), 0.01)
					}
					if r(N)==0 {
						matrix mean`Z'= -69
					}
					summ `yxsd'
					if r(N)~=0 {
						matrix sd`Z'= round(r(mean), 0.01)
					}
					if r(N)==0 {
						matrix sd`Z'= -69
					}
				}
			
				foreach Q in mean sd  {
					matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6
				}
			
				foreach Y in m3 m2  {
					tempvar temp swt y x yx swt
					gen `temp'=1 if `X'~=. & `tsex' ==`S'
					gen `swt'=``X'`Y''*`sw' if `temp'==1 
					egen `y'=sum(`swt') if `temp'==1 
					egen `x'=sum(`sw') if `temp'==1 
					gen `yx'=(`y'/`x') 
					summ `x'
					if r(N)~=0 {
						matrix A0=r(mean)
					}
					if r(N)==0 {
						matrix A0=0
					}
					summ `yx'
					if r(N)~=0 {
						matrix B0=r(mean)
					}
					if r(N)==0 {
						matrix B0= -69
					}
					foreach Z in 1 2 3 4 5 6 {
						tempvar temp swt y x yx swt
						gen `temp'=1 if `X'~=. & `agegrp'==`Z' & `tsex' ==`S'
						gen `swt'=``X'`Y''*`sw' if `temp'==1
						egen `y'=sum(`swt') if `temp'==1
						egen `x'=sum(`sw') if `temp'==1
						gen `yx'=(`y'/`x')
						summ `x'
						if r(N)~=0 {
							matrix A`Z'=r(mean)
						}
						if r(N)==0 {
							matrix A`Z'=0
						}
						summ `yx'
						if r(N)~=0 {
							matrix B`Z'=r(mean)
						}
						if r(N)==0 {
							matrix B`Z'= -69
						}
					}
					foreach Z in 0 1 2 3 4 5 6 {
						if A`Z'[1,1] ~=0 {
							matrix ga`Z'= inv(A`Z')
							matrix gb`Z'= inv(A`Z'*2)
							matrix D`Z' = 1-B`Z'[1,1]
							matrix F`Z' = 1.96*(B`Z'[1,1]*D`Z'[1,1]*ga`Z'[1,1])^0.5 + gb`Z'[1,1]
							matrix X`Z' = max(B`Z'[1,1] - F`Z'[1,1],0) 
							matrix Y`Z' = B`Z'[1,1] + F`Z'[1,1]
							matrix A`Z' = round(A`Z'[1,1],1)
							matrix B`Z' = round(B`Z'[1,1]*100, 0.1)
							matrix X`Z' = round(X`Z'[1,1]*100, 0.1)
							matrix Y`Z' = round(Y`Z'[1,1]*100, 0.1)
						}
						if A`Z'[1,1] ==0 {
							matrix X`Z' = -69
							matrix Y`Z' = -69
						}				
					}
					foreach Q in A B X Y  {
						matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6
					}
					matrix `Y'=B, X, Y
				}
				matrix `X'`S'= A, m3, m2, mean, sd
				matrix colnames `X'`S'=" N" "%<-3SD" "  95%" "  CI " "%<-2SD" "  95%" "  CI "  "Mean" "SD" 
				matrix rownames `X'`S'="Total" "0-5" "6-11" "12-23" "24-35" "36-47" "48-60"
			}
			foreach X in _zwfl _zbmi {
				tempvar temp swt y x yx swtsd ysd yxsd
				gen `temp'=1 if `X'~=. & `tsex' ==`S'
				gen `swt'=`X'*`sw' if `temp'==1 
				egen `y'=sum(`swt') if `temp'==1 
				egen `x'=sum(`sw') if `temp'==1
				gen `yx'=(`y'/`x')
				gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1
				egen `ysd'=sum(`swtsd') if `temp'==1
				gen `yxsd'=(`ysd'/(`x'-1))^0.5
				summ `yx'
				if r(N)~=0 {
					matrix mean0= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix mean0= -69
				}
				summ `yxsd'
				if r(N)~=0 {
					matrix sd0= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix sd0= -69
				}
				foreach Z in 1 2 3 4 5 6 {
					tempvar temp swt y x yx swtsd ysd yxsd
					gen `temp'=1 if `X'~=. & `agegrp'==`Z' & `tsex' ==`S'
					gen `swt'=`X'*`sw' if `temp'==1 
					egen `y'=sum(`swt') if `temp'==1 
					egen `x'=sum(`sw') if `temp'==1 
					gen `yx'=(`y'/`x')
					gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1
					egen `ysd'=sum(`swtsd') if `temp'==1
					gen `yxsd'=(`ysd'/(`x'-1))^0.5
					summ `yx'
					if r(N)~=0 {
						matrix mean`Z'= round(r(mean), 0.01)
					}
					if r(N)==0 {
						matrix mean`Z'= -69
					}
					summ `yxsd'
					if r(N)~=0 {
						matrix sd`Z'= round(r(mean), 0.01)
					}
					if r(N)==0 {
						matrix sd`Z'= -69
					}
				}
			
				foreach Q in mean sd  {
					matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6
				}
				foreach Y in m3 m2 p1 p2 p3 {
					tempvar temp swt y x yx swt
					gen `temp'=1 if `X'~=. & `tsex' ==`S'
					recode `temp' .=1 if `toedema'==2 & `agegrp'~=. & `tsex' ==`S'
					gen `swt'=``X'`Y''*`sw' if `temp'==1
					egen `y'=sum(`swt') if `temp'==1 
					egen `x'=sum(`sw') if `temp'==1
					gen `yx'=(`y'/`x')
					summ `x'
					if r(N)~=0 {
						matrix A0=r(mean)
					}
					if r(N)==0 {
						matrix A0=0
					}
					summ `yx'
					if r(N)~=0 {
						matrix B0=r(mean)
					}
					if r(N)==0 {
						matrix B0= -69
					}
	
					foreach Z in 1 2 3 4 5 6 {
						tempvar temp swt y x yx swt
						gen `temp'=1 if `X'~=. & `agegrp'==`Z' & `tsex' ==`S'
						recode `temp' .=1 if `toedema'==2 & `agegrp'==`Z' & `tsex' ==`S'
						gen `swt'=``X'`Y''*`sw' if `temp'==1
						egen `y'=sum(`swt') if `temp'==1
						egen `x'=sum(`sw') if `temp'==1
						gen `yx'=(`y'/`x')
						summ `x'
						if r(N)~=0 {
							matrix A`Z'=r(mean)
						}
						if r(N)==0 {
							matrix A`Z'=0
						}
						summ `yx'
						if r(N)~=0 {
							matrix B`Z'=r(mean)
						}
						if r(N)==0 {
							matrix B`Z'= -69
						}
					}
					foreach Z in 0 1 2 3 4 5 6 {
						if A`Z'[1,1] ~=0 {
							matrix ga`Z'= inv(A`Z')
							matrix gb`Z'= inv(A`Z'*2)
							matrix D`Z' = 1-B`Z'[1,1]
							matrix F`Z' = 1.96*(B`Z'[1,1]*D`Z'[1,1]*ga`Z'[1,1])^0.5 + gb`Z'[1,1]
							matrix X`Z' = max(B`Z'[1,1] - F`Z'[1,1],0) 
							matrix Y`Z' = B`Z'[1,1] + F`Z'[1,1]
							matrix A`Z' = round(A`Z'[1,1],1)
							matrix B`Z' = round(B`Z'[1,1]*100, 0.1)
							matrix X`Z' = round(X`Z'[1,1]*100, 0.1)
							matrix Y`Z' = round(Y`Z'[1,1]*100, 0.1)
						}
						if A`Z'[1,1] ==0 {
							matrix X`Z' = -69
							matrix Y`Z' = -69
						}				
					}
					foreach Q in A B X Y  {
						matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6
					}
					matrix `Y'=B, X, Y
				}
				matrix `X'`S'= A, m3, m2, p1,p2, p3, mean, sd
				matrix colnames `X'`S'=" N" "%<-3SD" "  95%" "  CI " "%<-2SD" "  95%" "  CI " "%>1SD" "  95%" "  CI " "%>2SD" "  95%" "  CI " "%>3SD" "  95%" "  CI " "Mean" "SD"
				matrix rownames `X'`S'="Total" "0-5" "6-11" "12-23" "24-35" "36-47" "48-60"
			}
			foreach X in _zhc _zac _zts _zss {
				tempvar temp swt y x yx swtsd ysd yxsd a
				gen `a'= match("`X'","_zhc")
				gen `temp'=1 if `X'~=. & `tsex' ==`S'
				gen `swt'=`X'*`sw' if `temp'==1 
				egen `y'=sum(`swt') if `temp'==1
				egen `x'=sum(`sw') if `temp'==1
				gen `yx'=(`y'/`x')
				gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1
				egen `ysd'=sum(`swtsd') if `temp'==1
				gen `yxsd'=(`ysd'/(`x'-1))^0.5
				summ `yx'
				if r(N)~=0 {
					matrix mean0= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix mean0= -69
				}
				summ `yxsd'
				if r(N)~=0 {
					matrix sd0= round(r(mean), 0.01)
				}
				if r(N)==0 {
					matrix sd0= -69
				}
				foreach Z in 1 2 3 4 5 6 {
					tempvar temp swt y x yx swtsd ysd yxsd
					gen `temp'=1 if `X'~=. & `agegrp'==`Z' & `tsex' ==`S'
					gen `swt'=`X'*`sw' if `temp'==1
					egen `y'=sum(`swt') if `temp'==1
					egen `x'=sum(`sw') if `temp'==1
					gen `yx'=(`y'/`x')
					gen `swtsd'=(`X'-`yx')^2*`sw' if `temp'==1
					egen `ysd'=sum(`swtsd') if `temp'==1
					gen `yxsd'=(`ysd'/(`x'-1))^0.5
					summ `yx'
					if r(N)~=0 {
						matrix mean`Z'= round(r(mean), 0.01)
					}
					if r(N)==0 {
						matrix mean`Z'= -69
					}
					summ `yxsd'
					if r(N)~=0 {
						matrix sd`Z'= round(r(mean), 0.01)
					}
					if r(N)==0 {
						matrix sd`Z'= -69
					}
				}
			
				foreach Q in mean sd  {
					matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6
				}
			
				foreach Y in m3 m2 p1 p2 p3 {
					tempvar temp swt y x yx swt
					gen `temp'=1 if `X'~=. & `tsex' ==`S'
					gen `swt'=``X'`Y''*`sw' if `temp'==1 
					egen `y'=sum(`swt') if `temp'==1 
					egen `x'=sum(`sw') if `temp'==1 
					gen `yx'=(`y'/`x') 
					summ `x'
					if r(N)~=0 {
						matrix A0=r(mean)
					}
					if r(N)==0 {
						matrix A0=0
					}
					summ `yx'
					if r(N)~=0 {
						matrix B0=r(mean)
					}
					if r(N)==0 {
						matrix B0= -69
					}
					foreach Z in 1 2 3 4 5 6 {
						tempvar temp swt y x yx swt
						gen `temp'=1 if `X'~=. & `agegrp'==`Z' & `tsex' ==`S'
						gen `swt'=``X'`Y''*`sw' if `temp'==1
						egen `y'=sum(`swt') if `temp'==1
						egen `x'=sum(`sw') if `temp'==1
						gen `yx'=(`y'/`x')
						summ `x'
						if r(N)~=0 {
							matrix A`Z'=r(mean)
						}
						if r(N)==0 {
							matrix A`Z'=0
						}
						summ `yx'
						if r(N)~=0 {
							matrix B`Z'=r(mean)
						}
						if r(N)==0 {
							matrix B`Z'= -69
						}
					}
					foreach Z in 0 1 2 3 4 5 6 {
						if A`Z'[1,1] ~=0 {
							matrix ga`Z'= inv(A`Z')
							matrix gb`Z'= inv(A`Z'*2)
							matrix D`Z' = 1-B`Z'[1,1]
							matrix F`Z' = 1.96*(B`Z'[1,1]*D`Z'[1,1]*ga`Z'[1,1])^0.5 + gb`Z'[1,1]
							matrix X`Z' = max(B`Z'[1,1] - F`Z'[1,1],0) 
							matrix Y`Z' = B`Z'[1,1] + F`Z'[1,1]
							matrix A`Z' = round(A`Z'[1,1],1)
							matrix B`Z' = round(B`Z'[1,1]*100, 0.1)
							matrix X`Z' = round(X`Z'[1,1]*100, 0.1)
							matrix Y`Z' = round(Y`Z'[1,1]*100, 0.1)
						}
						if A`Z'[1,1] ==0 {
							matrix X`Z' = -69
							matrix Y`Z' = -69
						}				
					}
					foreach Q in A B X Y  {
						matrix `Q'= `Q'0 \ `Q'1 \ `Q'2 \ `Q'3 \ `Q'4 \ `Q'5 \ `Q'6
					}
					matrix `Y'=B, X, Y
				}
				matrix `X'`S'= A, m3, m2, p1,p2, p3, mean, sd
				matrix colnames `X'`S'=" N" "%<-3SD" "  95%" "  CI " "%<-2SD" "  95%" "  CI " "%>1SD" "  95%" "  CI " "%>2SD" "  95%" "  CI " "%>3SD" "  95%" "  CI " "Mean" "SD"
				matrix rownames `X'`S'="Total" "3-5" "6-11" "12-23" "24-35" "36-47" "48-60"
				if `a'==1 {
					matrix rownames `X'`S'="Total" "0-5" "6-11" "12-23" "24-35" "36-47" "48-60"
				}
			}
		}
		
	} /* quietly ends*/
	
	} 	/* end of loop "if r(min)~=0" */

	/**		On-screen display	**/
			
	di as txt _n "Set 1:		Sexes combined"
	di as txt _n "Age groups 	Weight-for-age"
	matrix list _zwei, noheader 
	di as txt _n "Age groups	Length/height-for-age"
	matrix list _zlen, noheader 
	di as txt _n "Age groups 	Weight-for-length/height"
	matrix list _zwfl, noheader 
	di as txt _n "Age groups	BMI-for-age"
	matrix list _zbmi, noheader 
	
	di as txt _n "Set 2:		Males"
	di as txt _n "Age groups 	Weight-for-age"
	matrix list _zwei1, noheader 
	di as txt _n "Age groups	Length/height-for-age"
	matrix list _zlen1, noheader 
	di as txt _n "Age groups 	Weight-for-length/height"
	matrix list _zwfl1, noheader 
	di as txt _n "Age groups	BMI-for-age"
	matrix list _zbmi1, noheader 
	
	di as txt _n "Set 3:		Females"
	di as txt _n "Age groups 	Weight-for-age"
	matrix list _zwei2, noheader 
	di as txt _n "Age groups	Length/height-for-age"
	matrix list _zlen2, noheader 
	di as txt _n "Age groups 	Weight-for-length/height"
	matrix list _zwfl2, noheader 
	di as txt _n "Age groups	BMI-for-age"
	matrix list _zbmi2, noheader 
	
	di as txt _n "Set 1:		Sexes combined"
	di as txt _n "Age groups 	Head circumference-for-age"
	matrix list _zhc, noheader 
	di as txt _n "Age groups	MUAC-for-age"
	matrix list _zac, noheader 
	di as txt _n "Age groups 	Triceps skinfold-for-age"
	matrix list _zts, noheader 
	di as txt _n "Age groups	Subscapular skinfold-for-age"
	matrix list _zss, noheader 
	
	di as txt _n "Set 2:		Males"
	di as txt _n "Age groups 	Head circumference-for-age"
	matrix list _zhc1, noheader 
	di as txt _n "Age groups	MUAC-for-age"
	matrix list _zac1, noheader 
	di as txt _n "Age groups 	Triceps skinfold-for-age"
	matrix list _zts1, noheader 
	di as txt _n "Age groups	Subscapular skinfold-for-age"
	matrix list _zss1, noheader 
	
	di as txt _n "Set 3:		Females"
	di as txt _n "Age groups 	Head circumference-for-age"
	matrix list _zhc2, noheader 
	di as txt _n "Age groups	MUAC-for-age"
	matrix list _zac2, noheader 
	di as txt _n "Age groups 	Triceps skinfold-for-age"
	matrix list _zts2, noheader 
	di as txt _n "Age groups	Subscapular skinfold-for-age"
	matrix list _zss2, noheader 
	
	/**	Prevalence tables to the excel sheet begin here	**/
	
	/* quietly begins */
	qui {
	#delimit ;
	gen __t1=-1000; gen __t2=-2000; gen __t3=-9000; gen __t4=-10000; gen __t5=-3000; gen __t6=-9000; 
	gen __t7=-10000; gen __t8=-6000; gen __t9=-9000; gen __t10=-10000; gen __t11=-7000; gen __t12=-9000;
	gen __t13=-10000; gen __t14=-8000; gen __t15=-9000; gen __t16=-10000; gen __t17=-4000; gen __t18=-5000;
	#delimit cr
	
	
	forval i = 1/7 {
		gen __u`i'=__t`i'
	}
	gen __u8=-4000
	gen __u9=-5000
	
	mkmat __u1-__u9 in 1 , matrix(U)
	mkmat __t1-__t18 in 1, matrix(T)
	
	matrix emp=J(7,9,-6666)
	matrix uemp=J(1,9,-6666)
	forval i = 1/2 {
		matrix _zwei`i'=_zwei`i', emp
		matrix _zlen`i'=_zlen`i', emp
	}
	matrix _zwei=_zwei, emp
	matrix _zlen=_zlen, emp
	matrix U=U, uemp
	
	*set trace on
	matrix first=J(228,1,0)
	forval i = 1/228 {
		matrix first[`i',1]=`i'
	}
	matrix row1=J(1,18,-6666)
	matrix row1[1,1]=-11	/*	Set1:*/
	matrix row1[1,2]=-22	/*	Sexes*/
	matrix row1[1,3]=-33	/*	combined*/
	
	matrix rowM=J(1,18,-6666)
	matrix rowM[1,1]=-1155	/*	Set2:*/
	matrix rowM[1,2]=-2255	/*	Males*/
	
	matrix rowF=J(1,18,-6666)
	matrix rowF[1,1]=-1166	/*	Set3:*/
	matrix rowF[1,2]=-2266	/*	Females*/
	
	matrix row2=J(1,18,-6666)
	matrix row2[1,1]=-111		/*	Weight*/
	matrix row2[1,2]=-222		/*	-for-*/
	matrix row2[1,3]=-333		/*	age*/
	
	matrix row11=J(1,18,-6666)
	matrix row11[1,1]=-1111		/*	Length*/
	matrix row11[1,2]=-2222		/*	/height*/
	matrix row11[1,3]=-3333		/*	-for-*/
	matrix row11[1,4]=-4444		/*	age*/
	
	matrix row20=J(1,18,-6666)
	matrix row20[1,1]=-11111	/*	Weight*/
	matrix row20[1,2]=-22222	/*	-for-*/
	matrix row20[1,3]=-33333	/*	length*/
	matrix row20[1,4]=-44444	/*	/height*/
	
	matrix row29=J(1,18,-6666)
	matrix row29[1,1]=-111111	/*	BMI*/
	matrix row29[1,2]=-222222	/*	-for-*/
	matrix row29[1,3]=-333333	/*	age*/
	
	matrix row38=J(1,18,-7777)

	matrix srow2=J(1,18,-6666)
	matrix srow2[1,1]=-1111111	/*	Head*/
	matrix srow2[1,2]=-2222222	/*	circumference*/
	matrix srow2[1,3]=-3333333	/*	-for-*/
	matrix srow2[1,4]=-4444444	/*	age*/
	
	matrix srow11=J(1,18,-6666)
	matrix srow11[1,1]=-818181	/*	MUAC*/
	matrix srow11[1,2]=-828282	/*	-for-*/
	matrix srow11[1,3]=-838383	/*	age*/
		
	matrix srow20=J(1,18,-6666)
	matrix srow20[1,1]=-848484	/*	Triceps*/
	matrix srow20[1,2]=-858585	/*	skinfold*/
	matrix srow20[1,3]=-868686	/*	-for-*/
	matrix srow20[1,4]=-878787	/*	age*/
	
	matrix srow29=J(1,18,-6666)
	matrix srow29[1,1]=-888888	/*	Subscapular*/
	matrix srow29[1,2]=-898989	/*	skinfold*/
	matrix srow29[1,3]=-808080	/*	-for-*/
	matrix srow29[1,4]=-919191	/*	age*/
	
	matrix row38=J(1,18,-7777)

	matrix Comb=row1 \ row2 \ U \ _zwei \ row11 \ U \ _zlen \ row20 \ T \ _zwfl \ row29 \ T \ _zbmi \ row38
	matrix Comb1=rowM \ row2 \ U \ _zwei1 \ row11 \ U \ _zlen1 \ row20 \ T \ _zwfl1 \ row29 \ T \ _zbmi1 \ row38
	matrix Comb2=rowF \ row2 \ U \ _zwei2 \ row11 \ U \ _zlen2 \ row20 \ T \ _zwfl2 \ row29 \ T \ _zbmi2 \ row38
	matrix CombA=Comb \ Comb1 \ Comb2 

	matrix Somb=row1 \ srow2 \ T \ _zhc \ srow11 \ T \ _zac \ srow20 \ T \ _zts \ srow29 \ T \ _zss \ row38
	matrix Somb1=rowM \ srow2 \ T \ _zhc1 \ srow11 \ T \ _zac1 \ srow20 \ T \ _zts1 \ srow29 \ T \ _zss1 \ row38
	matrix Somb2=rowF \ srow2 \ T \ _zhc2 \ srow11 \ T \ _zac2 \ srow20 \ T \ _zts2 \ srow29 \ T \ _zss2 \ row38
	matrix SombA=Somb \ Somb1 \ Somb2 
	matrix temp= CombA \ SombA
	matrix __CombA=first, temp
	svmat __CombA
	} 
	/* quietly ends*/
	
	qui gen str3 __xCombA1 = string(__CombA1)
	qui replace __xCombA1 = " " if __xCombA1 == "-6666"
	qui replace __xCombA1 = " " if __xCombA1 == "38" | __xCombA1 == "76" | __xCombA1 == "114" | __xCombA1 == "152" | __xCombA1 == "190" | __xCombA1 == "228"
	
	foreach i in 1 2 11 20 29  39 40 49 58 67  77 78 87 96 105  115 116 125 134 143  153 154 163 172 181  191 192 201 210 219 {
		qui replace __xCombA1=" " if __xCombA1=="`i'" 
	}
	foreach i in 3 12 21 30  41 50 59 68  79 88 97 106  117 126 135 144  155 164 173 182  193 202 211 220 {
		qui replace __xCombA1="Age" if __xCombA1=="`i'" 
	}
	foreach i in 127 136 145   165 174 183  203 212 221 {
		qui replace __xCombA1="(3-60)" if __xCombA1=="`i'" 
	}
	foreach i in 4 13 22 31  42 51 60 69  80 89 98 107  118 156  194 {
		qui replace __xCombA1="(0-60)" if __xCombA1=="`i'" 
	}
	foreach i in 128 137 146  166 175 184 204 213 222 {
		qui replace __xCombA1="(3-5)" if __xCombA1=="`i'" 
	}
	foreach i in 5 14 23 32  43 52 61 70  81 90 99 108  119 157 195 {
		qui replace __xCombA1="(0-5)" if __xCombA1=="`i'" 
	}
	foreach i in 128 137 146  166 175 184 204 213 222 {
		qui replace __xCombA1="(3-5)" if __xCombA1=="`i'" 
	}
	foreach i in 6 15 24 33  44 53 62 71  82 91 100 109  120 129 138 147  158 167 176 185  196 205 214 223 {
		qui replace __xCombA1="(6-11)" if __xCombA1=="`i'" 
	}
	foreach i in 7 16 25 34  45 54 63 72  83 92 101 110  121 130 139 148  159 168 177 186  197 206 215 224 {
		qui replace __xCombA1="(12-23)" if __xCombA1=="`i'" 
	}
	foreach i in 8 17 26 35  46 55 64 73  84 93 102 111  122 131 140 149  160 169 178 187  198 207 216 225 {
		qui replace __xCombA1="(24-35)" if __xCombA1=="`i'" 
	}
	foreach i in 9 18 27 36  47 56 65 74  85 94 103 112  123 132 141 150  161 170 179 188  199 208 217 226 {
		qui replace __xCombA1="(36-47)" if __xCombA1=="`i'" 
	}
	foreach i in 10 19 28 37  48 57 66 75  86 95 104 113 124 133 142 151  162 171 180 189  200 209 218 227 {
		qui replace __xCombA1="(48-60)" if __xCombA1=="`i'" 
	}
	forval i = 2/19 {
		qui gen str10 __xCombA`i' = string(__CombA`i')
		qui replace __xCombA`i' = "N" if __xCombA`i' == "-1000"
		qui replace __xCombA`i' = "%<-3SD" if __xCombA`i' == "-2000"
		qui replace __xCombA`i' = "%<-2SD" if __xCombA`i' == "-3000"
		qui replace __xCombA`i' = "Mean" if __xCombA`i' == "-4000"
		qui replace __xCombA`i' = "SD" if __xCombA`i' == "-5000"
		qui replace __xCombA`i' = "%>+1SD" if __xCombA`i' == "-6000"
		qui replace __xCombA`i' = "%>+2SD" if __xCombA`i' == "-7000"
		qui replace __xCombA`i' = "%>+3SD" if __xCombA`i' == "-8000"
		qui replace __xCombA`i' = "95%" if __xCombA`i' == "-9000"
		qui replace __xCombA`i' = "C.I." if __xCombA`i' == "-10000"
		qui replace __xCombA`i' = "Set 1:" if __xCombA`i' == "-11"
		qui replace __xCombA`i' = "Sexes" if __xCombA`i' == "-22"
		qui replace __xCombA`i' = "combined" if __xCombA`i' == "-33"
		qui replace __xCombA`i' = "Set 2:" if __xCombA`i' == "-1155"
		qui replace __xCombA`i' = "Males" if __xCombA`i' == "-2255"
		qui replace __xCombA`i' = "Set 3:" if __xCombA`i' == "-1166"
		qui replace __xCombA`i' = "Females" if __xCombA`i' == "-2266"
		qui replace __xCombA`i' = "Weight" if __xCombA`i' == "-111"
		qui replace __xCombA`i' = "-for-" if __xCombA`i' == "-222"
		qui replace __xCombA`i' = "age" if __xCombA`i' == "-333"
		qui replace __xCombA`i' = "Length" if __xCombA`i' == "-1111"
		qui replace __xCombA`i' = "/height" if __xCombA`i' == "-2222"
		qui replace __xCombA`i' = "-for-" if __xCombA`i' == "-3333"
		qui replace __xCombA`i' = "age" if __xCombA`i' == "-4444"
		qui replace __xCombA`i' = "Weight" if __xCombA`i' == "-11111"
		qui replace __xCombA`i' = "-for-" if __xCombA`i' == "-22222"
		qui replace __xCombA`i' = "length" if __xCombA`i' == "-33333"
		qui replace __xCombA`i' = "/height" if __xCombA`i' == "-44444"
		qui replace __xCombA`i' = "BMI" if __xCombA`i' == "-111111"
		qui replace __xCombA`i' = "-for-" if __xCombA`i' == "-222222"
		qui replace __xCombA`i' = "age" if __xCombA`i' == "-333333"
		qui replace __xCombA`i' = "Head" if __xCombA`i' == "-1111111"
		qui replace __xCombA`i' = "circumference" if __xCombA`i' == "-2222222"
		qui replace __xCombA`i' = "-for-" if __xCombA`i' == "-3333333"
		qui replace __xCombA`i' = "age" if __xCombA`i' == "-4444444"
		qui replace __xCombA`i' = "MUAC" if __xCombA`i' == "-818181"
		qui replace __xCombA`i' = "-for-" if __xCombA`i' == "-828282"
		qui replace __xCombA`i' = "age" if __xCombA`i' == "-838383"
		qui replace __xCombA`i' = "Triceps" if __xCombA`i' == "-848484"
		qui replace __xCombA`i' = "skinfold" if __xCombA`i' == "-858585"
		qui replace __xCombA`i' = "-for-" if __xCombA`i' == "-868686"
		qui replace __xCombA`i' = "age" if __xCombA`i' == "-878787"
		qui replace __xCombA`i' = "Subscapular" if __xCombA`i' == "-888888"
		qui replace __xCombA`i' = "skinfold" if __xCombA`i' == "-898989"
		qui replace __xCombA`i' = "-for-" if __xCombA`i' == "-808080"
		qui replace __xCombA`i' = "age" if __xCombA`i' == "-919191"
		qui replace __xCombA`i' = " " if __xCombA`i' == "-6666"
		qui replace __xCombA`i' = " " if __xCombA`i' == "-7777"
		qui replace __xCombA`i' = " " if __xCombA`i' == "-69"
	}

	/***	Writing out to the worksheet **/
	local string "xx\yy_prev_st.xls"
	local i=`datalib'
	local j=`datalab'
	global outp: subinstr local string "xx" "`i'"
	global outprev: subinstr global outp "yy" "`j'"
	outsheet __xCombA1-__xCombA19 using "$outprev" in 1/228, nonames nolabel replace
	
	di " "
	di "Note:	Prevalences are written to"
	di "		$outprev"
	
	/* Cleaning up after Phase II*/
	drop _*
	tempvar clear
end
exit
