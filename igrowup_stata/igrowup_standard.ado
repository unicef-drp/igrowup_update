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
	///version 7.0
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
	
	local string "xxx\yy_prev_st.dta"
	local i=`datalib'
	local j=`datalab'
	global finalf: subinstr local string "xxx" "`i'" 
	global finalfile: subinstr global finalf "yy" "`j'"	
	
	local string "xxx\yy_prev_st.csv"
	local i=`datalib'
	local j=`datalab'
	global exportf: subinstr local string "xxx" "`i'" 
	global exportfile: subinstr global exportf "yy" "`j'"	
	
	/// Generation of age cap (0-59 months)
	gen _age59=0
	replace _age59=1 if _agedays<1826.25
	replace _age59=0 if _agedays==. // 1826.25 is 365.25 (number of days in a year) multiplied by 5
	replace _age59=. if _age59==0
	
	gen _age59over=1 if _agedays>=1826.25
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
	
	capture confirm string variable gender																// Confirmation loop to confirm sex is string
		if !_rc {	
                       replace _sex=1 if gender=="m" | gender=="M"
					   replace _sex=2 if gender=="f" | gender=="F"
               }
               else {
                       replace _sex=1 if gender==1
					   replace _sex=2 if gender==2
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
	label define _agemonthgroup_sexlabel 1 "Male 0 to 5 months" 2 "Male 6 to 11 months" 3 "Male 12 to 23 months" 4 "Male 24 to 35 months" 5 "Male 36 to 47 months" 6 "Male 48 to 59 months"	7 "Female 0 to 5 months" 8 "Female 6 to 11 months" 9 "Female 12 to 23 months" 10 "Female 24 to 35 months" 11 "Female 36 to 47 months" 12 "Female 48 to 59 months"
	
	
	/// Generate Age Sex Disaggregation
	
	gen _agemonthgroup_sex=.
	replace _agemonthgroup_sex=1 if _sex==1 & _agemonthgroup==1
	replace _agemonthgroup_sex=2 if _sex==1 & _agemonthgroup==2
	replace _agemonthgroup_sex=3 if _sex==1 & _agemonthgroup==3
	replace _agemonthgroup_sex=4 if _sex==1 & _agemonthgroup==4
	replace _agemonthgroup_sex=5 if _sex==1 & _agemonthgroup==5
	replace _agemonthgroup_sex=6 if _sex==1 & _agemonthgroup==6
	replace _agemonthgroup_sex=7 if _sex==2 & _agemonthgroup==1
	replace _agemonthgroup_sex=8 if _sex==2 & _agemonthgroup==2
	replace _agemonthgroup_sex=9 if _sex==2 & _agemonthgroup==3
	replace _agemonthgroup_sex=10 if _sex==2 & _agemonthgroup==4
	replace _agemonthgroup_sex=11 if _sex==2 & _agemonthgroup==5
	replace _agemonthgroup_sex=12 if _sex==2 & _agemonthgroup==6
	
	/// Apply Value labels to variables
	
	foreach x in _age59 _sex _agemonthgroup _agemonthgroup_sex {
	label values `x' `x'label
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
	

	

	/// Levels 
	gen all=1	
	
	sum all 
	local all_check=r(N)
	
	if `all_check'>0 {
		gen str sta_all="Total"
	}
	else {
		gen sta_all=""
	}
	

	foreach x in _sex _agemonthgroup _agemonthgroup_sex { 
		sum `x'
		local `x'_check=r(N)
		if ``x'_check'>0 {
			decode `x', generate(sta_`x')
						}
		else {
			gen str sta_`x'=" "
			}
	}

	
		save "$workfile", replace	
	
	local string "xxx\yy_z_working_prev_st.dta"
	local i=`datalib'
	local j=`datalab'
	global prevf: subinstr local string "xxx" "`i'" 
	global prevfile: subinstr global prevf "yy" "`j'"
	
	clear
	
	gen str strat_label=" "
	gen str standard_strat=" "

	label variable strat_label "Statification Label"
	label variable standard_strat "Standardized Stratification"
	

	
	foreach x in haz waz whz baz {
	foreach z in _3 _2 _1 1 2 3 {
		foreach y in r se t pvalue ll ul df crit weighted_N unweighted_N weighted_Num unweighted_Num {
				gen `x'`z'_`y'=.
	}
	}
	}
	
	foreach x in haz waz whz baz {
	foreach y in r se t pvalue ll ul df crit weighted_N unweighted_N weighted_Num unweighted_Num {
		gen `x'_`y'=.
	}
}
	
	foreach x in haz waz whz baz {
	foreach y in sd {
		gen `x'_`y'=.
	}
}


/// Label variables 

label variable haz_3_r "height-for-age <-3SD [point estimate]"
label variable haz_3_se "height-for-age <-3SD [standard error]"
label variable haz_3_t "height-for-age <-3SD [t-statistic]"
label variable haz_3_pvalue "height-for-age <-3SD [p value]"
label variable haz_3_ll "height-for-age <-3SD [95% lower confidence limit]"
label variable haz_3_ul "height-for-age <-3SD [95% upper confidence limit]"
label variable haz_3_df "height-for-age <-3SD [degrees of freedom]"
label variable haz_3_crit "height-for-age <-3SD [critical value]"
label variable haz_3_weighted_N "height-for-age <-3SD [weighted denominator]"
label variable haz_3_unweighted_N "height-for-age <-3SD [unweighted denominator]"
label variable haz_3_weighted_Num "height-for-age <-3SD [weighted numerator]"
label variable haz_3_unweighted_Num "height-for-age <-3SD [unweighted numerator]"
label variable haz_2_r "height-for-age <-2SD [point estimate]"
label variable haz_2_se "height-for-age <-2SD [standard error]"
label variable haz_2_t "height-for-age <-2SD [t-statistic]"
label variable haz_2_pvalue "height-for-age <-2SD [p value]"
label variable haz_2_ll "height-for-age <-2SD [95% lower confidence limit]"
label variable haz_2_ul "height-for-age <-2SD [95% upper confidence limit]"
label variable haz_2_df "height-for-age <-2SD [degrees of freedom]"
label variable haz_2_crit "height-for-age <-2SD [critical value]"
label variable haz_2_weighted_N "height-for-age <-2SD [weighted denominator]"
label variable haz_2_unweighted_N "height-for-age <-2SD [unweighted denominator]"
label variable haz_2_weighted_Num "height-for-age <-2SD [weighted numerator]"
label variable haz_2_unweighted_Num "height-for-age <-2SD [unweighted numerator]"
label variable haz_1_r "height-for-age <-1SD [point estimate]"
label variable haz_1_se "height-for-age <-1SD [standard error]"
label variable haz_1_t "height-for-age <-1SD [t-statistic]"
label variable haz_1_pvalue "height-for-age <-1SD [p value]"
label variable haz_1_ll "height-for-age <-1SD [95% lower confidence limit]"
label variable haz_1_ul "height-for-age <-1SD [95% upper confidence limit]"
label variable haz_1_df "height-for-age <-1SD [degrees of freedom]"
label variable haz_1_crit "height-for-age <-1SD [critical value]"
label variable haz_1_weighted_N "height-for-age <-1SD [weighted denominator]"
label variable haz_1_unweighted_N "height-for-age <-1SD [unweighted denominator]"
label variable haz_1_weighted_Num "height-for-age <-1SD [weighted numerator]"
label variable haz_1_unweighted_Num "height-for-age <-1SD [unweighted numerator]"
label variable haz1_r "height-for-age >+1SD [point estimate]"
label variable haz1_se "height-for-age >+1SD [standard error]"
label variable haz1_t "height-for-age >+1SD [t-statistic]"
label variable haz1_pvalue "height-for-age >+1SD [p value]"
label variable haz1_ll "height-for-age >+1SD [95% lower confidence limit]"
label variable haz1_ul "height-for-age >+1SD [95% upper confidence limit]"
label variable haz1_df "height-for-age >+1SD [degrees of freedom]"
label variable haz1_crit "height-for-age >+1SD [critical value]"
label variable haz1_weighted_N "height-for-age >+1SD [weighted denominator]"
label variable haz1_unweighted_N "height-for-age >+1SD [unweighted denominator]"
label variable haz1_weighted_Num "height-for-age >+1SD [weighted numerator]"
label variable haz1_unweighted_Num "height-for-age >+1SD [unweighted numerator]"
label variable haz2_r "height-for-age >+2SD [point estimate]"
label variable haz2_se "height-for-age >+2SD [standard error]"
label variable haz2_t "height-for-age >+2SD [t-statistic]"
label variable haz2_pvalue "height-for-age >+2SD [p value]"
label variable haz2_ll "height-for-age >+2SD [95% lower confidence limit]"
label variable haz2_ul "height-for-age >+2SD [95% upper confidence limit]"
label variable haz2_df "height-for-age >+2SD [degrees of freedom]"
label variable haz2_crit "height-for-age >+2SD [critical value]"
label variable haz2_weighted_N "height-for-age >+2SD [weighted denominator]"
label variable haz2_unweighted_N "height-for-age >+2SD [unweighted denominator]"
label variable haz2_weighted_Num "height-for-age >+2SD [weighted numerator]"
label variable haz2_unweighted_Num "height-for-age >+2SD [unweighted numerator]"
label variable haz3_r "height-for-age >+3SD [point estimate]"
label variable haz3_se "height-for-age >+3SD [standard error]"
label variable haz3_t "height-for-age >+3SD [t-statistic]"
label variable haz3_pvalue "height-for-age >+3SD [p value]"
label variable haz3_ll "height-for-age >+3SD [95% lower confidence limit]"
label variable haz3_ul "height-for-age >+3SD [95% upper confidence limit]"
label variable haz3_df "height-for-age >+3SD [degrees of freedom]"
label variable haz3_crit "height-for-age >+3SD [critical value]"
label variable haz3_weighted_N "height-for-age >+3SD [weighted denominator]"
label variable haz3_unweighted_N "height-for-age >+3SD [unweighted denominator]"
label variable haz3_weighted_Num "height-for-age >+3SD [weighted numerator]"
label variable haz3_unweighted_Num "height-for-age >+3SD [unweighted numerator]"
label variable waz_3_r "weight-for-age <-3SD [point estimate]"
label variable waz_3_se "weight-for-age <-3SD [standard error]"
label variable waz_3_t "weight-for-age <-3SD [t-statistic]"
label variable waz_3_pvalue "weight-for-age <-3SD [p value]"
label variable waz_3_ll "weight-for-age <-3SD [95% lower confidence limit]"
label variable waz_3_ul "weight-for-age <-3SD [95% upper confidence limit]"
label variable waz_3_df "weight-for-age <-3SD [degrees of freedom]"
label variable waz_3_crit "weight-for-age <-3SD [critical value]"
label variable waz_3_weighted_N "weight-for-age <-3SD [weighted denominator]"
label variable waz_3_unweighted_N "weight-for-age <-3SD [unweighted denominator]"
label variable waz_3_weighted_Num "weight-for-age <-3SD [weighted numerator]"
label variable waz_3_unweighted_Num "weight-for-age <-3SD [unweighted numerator]"
label variable waz_2_r "weight-for-age <-2SD [point estimate]"
label variable waz_2_se "weight-for-age <-2SD [standard error]"
label variable waz_2_t "weight-for-age <-2SD [t-statistic]"
label variable waz_2_pvalue "weight-for-age <-2SD [p value]"
label variable waz_2_ll "weight-for-age <-2SD [95% lower confidence limit]"
label variable waz_2_ul "weight-for-age <-2SD [95% upper confidence limit]"
label variable waz_2_df "weight-for-age <-2SD [degrees of freedom]"
label variable waz_2_crit "weight-for-age <-2SD [critical value]"
label variable waz_2_weighted_N "weight-for-age <-2SD [weighted denominator]"
label variable waz_2_unweighted_N "weight-for-age <-2SD [unweighted denominator]"
label variable waz_2_weighted_Num "weight-for-age <-2SD [weighted numerator]"
label variable waz_2_unweighted_Num "weight-for-age <-2SD [unweighted numerator]"
label variable waz_1_r "weight-for-age <-1SD [point estimate]"
label variable waz_1_se "weight-for-age <-1SD [standard error]"
label variable waz_1_t "weight-for-age <-1SD [t-statistic]"
label variable waz_1_pvalue "weight-for-age <-1SD [p value]"
label variable waz_1_ll "weight-for-age <-1SD [95% lower confidence limit]"
label variable waz_1_ul "weight-for-age <-1SD [95% upper confidence limit]"
label variable waz_1_df "weight-for-age <-1SD [degrees of freedom]"
label variable waz_1_crit "weight-for-age <-1SD [critical value]"
label variable waz_1_weighted_N "weight-for-age <-1SD [weighted denominator]"
label variable waz_1_unweighted_N "weight-for-age <-1SD [unweighted denominator]"
label variable waz_1_weighted_Num "weight-for-age <-1SD [weighted numerator]"
label variable waz_1_unweighted_Num "weight-for-age <-1SD [unweighted numerator]"
label variable waz1_r "weight-for-age >+1SD [point estimate]"
label variable waz1_se "weight-for-age >+1SD [standard error]"
label variable waz1_t "weight-for-age >+1SD [t-statistic]"
label variable waz1_pvalue "weight-for-age >+1SD [p value]"
label variable waz1_ll "weight-for-age >+1SD [95% lower confidence limit]"
label variable waz1_ul "weight-for-age >+1SD [95% upper confidence limit]"
label variable waz1_df "weight-for-age >+1SD [degrees of freedom]"
label variable waz1_crit "weight-for-age >+1SD [critical value]"
label variable waz1_weighted_N "weight-for-age >+1SD [weighted denominator]"
label variable waz1_unweighted_N "weight-for-age >+1SD [unweighted denominator]"
label variable waz1_weighted_Num "weight-for-age >+1SD [weighted numerator]"
label variable waz1_unweighted_Num "weight-for-age >+1SD [unweighted numerator]"
label variable waz2_r "weight-for-age >+2SD [point estimate]"
label variable waz2_se "weight-for-age >+2SD [standard error]"
label variable waz2_t "weight-for-age >+2SD [t-statistic]"
label variable waz2_pvalue "weight-for-age >+2SD [p value]"
label variable waz2_ll "weight-for-age >+2SD [95% lower confidence limit]"
label variable waz2_ul "weight-for-age >+2SD [95% upper confidence limit]"
label variable waz2_df "weight-for-age >+2SD [degrees of freedom]"
label variable waz2_crit "weight-for-age >+2SD [critical value]"
label variable waz2_weighted_N "weight-for-age >+2SD [weighted denominator]"
label variable waz2_unweighted_N "weight-for-age >+2SD [unweighted denominator]"
label variable waz2_weighted_Num "weight-for-age >+2SD [weighted numerator]"
label variable waz2_unweighted_Num "weight-for-age >+2SD [unweighted numerator]"
label variable waz3_r "weight-for-age >+3SD [point estimate]"
label variable waz3_se "weight-for-age >+3SD [standard error]"
label variable waz3_t "weight-for-age >+3SD [t-statistic]"
label variable waz3_pvalue "weight-for-age >+3SD [p value]"
label variable waz3_ll "weight-for-age >+3SD [95% lower confidence limit]"
label variable waz3_ul "weight-for-age >+3SD [95% upper confidence limit]"
label variable waz3_df "weight-for-age >+3SD [degrees of freedom]"
label variable waz3_crit "weight-for-age >+3SD [critical value]"
label variable waz3_weighted_N "weight-for-age >+3SD [weighted denominator]"
label variable waz3_unweighted_N "weight-for-age >+3SD [unweighted denominator]"
label variable waz3_weighted_Num "weight-for-age >+3SD [weighted numerator]"
label variable waz3_unweighted_Num "weight-for-age >+3SD [unweighted numerator]"
label variable whz_3_r "weight-for-height <-3SD [point estimate]"
label variable whz_3_se "weight-for-height <-3SD [standard error]"
label variable whz_3_t "weight-for-height <-3SD [t-statistic]"
label variable whz_3_pvalue "weight-for-height <-3SD [p value]"
label variable whz_3_ll "weight-for-height <-3SD [95% lower confidence limit]"
label variable whz_3_ul "weight-for-height <-3SD [95% upper confidence limit]"
label variable whz_3_df "weight-for-height <-3SD [degrees of freedom]"
label variable whz_3_crit "weight-for-height <-3SD [critical value]"
label variable whz_3_weighted_N "weight-for-height <-3SD [weighted denominator]"
label variable whz_3_unweighted_N "weight-for-height <-3SD [unweighted denominator]"
label variable whz_3_weighted_Num "weight-for-height <-3SD [weighted numerator]"
label variable whz_3_unweighted_Num "weight-for-height <-3SD [unweighted numerator]"
label variable whz_2_r "weight-for-height <-2SD [point estimate]"
label variable whz_2_se "weight-for-height <-2SD [standard error]"
label variable whz_2_t "weight-for-height <-2SD [t-statistic]"
label variable whz_2_pvalue "weight-for-height <-2SD [p value]"
label variable whz_2_ll "weight-for-height <-2SD [95% lower confidence limit]"
label variable whz_2_ul "weight-for-height <-2SD [95% upper confidence limit]"
label variable whz_2_df "weight-for-height <-2SD [degrees of freedom]"
label variable whz_2_crit "weight-for-height <-2SD [critical value]"
label variable whz_2_weighted_N "weight-for-height <-2SD [weighted denominator]"
label variable whz_2_unweighted_N "weight-for-height <-2SD [unweighted denominator]"
label variable whz_2_weighted_Num "weight-for-height <-2SD [weighted numerator]"
label variable whz_2_unweighted_Num "weight-for-height <-2SD [unweighted numerator]"
label variable whz_1_r "weight-for-height <-1SD [point estimate]"
label variable whz_1_se "weight-for-height <-1SD [standard error]"
label variable whz_1_t "weight-for-height <-1SD [t-statistic]"
label variable whz_1_pvalue "weight-for-height <-1SD [p value]"
label variable whz_1_ll "weight-for-height <-1SD [95% lower confidence limit]"
label variable whz_1_ul "weight-for-height <-1SD [95% upper confidence limit]"
label variable whz_1_df "weight-for-height <-1SD [degrees of freedom]"
label variable whz_1_crit "weight-for-height <-1SD [critical value]"
label variable whz_1_weighted_N "weight-for-height <-1SD [weighted denominator]"
label variable whz_1_unweighted_N "weight-for-height <-1SD [unweighted denominator]"
label variable whz_1_weighted_Num "weight-for-height <-1SD [weighted numerator]"
label variable whz_1_unweighted_Num "weight-for-height <-1SD [unweighted numerator]"
label variable whz1_r "weight-for-height >+1SD [point estimate]"
label variable whz1_se "weight-for-height >+1SD [standard error]"
label variable whz1_t "weight-for-height >+1SD [t-statistic]"
label variable whz1_pvalue "weight-for-height >+1SD [p value]"
label variable whz1_ll "weight-for-height >+1SD [95% lower confidence limit]"
label variable whz1_ul "weight-for-height >+1SD [95% upper confidence limit]"
label variable whz1_df "weight-for-height >+1SD [degrees of freedom]"
label variable whz1_crit "weight-for-height >+1SD [critical value]"
label variable whz1_weighted_N "weight-for-height >+1SD [weighted denominator]"
label variable whz1_unweighted_N "weight-for-height >+1SD [unweighted denominator]"
label variable whz1_weighted_Num "weight-for-height >+1SD [weighted numerator]"
label variable whz1_unweighted_Num "weight-for-height >+1SD [unweighted numerator]"
label variable whz2_r "weight-for-height >+2SD [point estimate]"
label variable whz2_se "weight-for-height >+2SD [standard error]"
label variable whz2_t "weight-for-height >+2SD [t-statistic]"
label variable whz2_pvalue "weight-for-height >+2SD [p value]"
label variable whz2_ll "weight-for-height >+2SD [95% lower confidence limit]"
label variable whz2_ul "weight-for-height >+2SD [95% upper confidence limit]"
label variable whz2_df "weight-for-height >+2SD [degrees of freedom]"
label variable whz2_crit "weight-for-height >+2SD [critical value]"
label variable whz2_weighted_N "weight-for-height >+2SD [weighted denominator]"
label variable whz2_unweighted_N "weight-for-height >+2SD [unweighted denominator]"
label variable whz2_weighted_Num "weight-for-height >+2SD [weighted numerator]"
label variable whz2_unweighted_Num "weight-for-height >+2SD [unweighted numerator]"
label variable whz3_r "weight-for-height >+3SD [point estimate]"
label variable whz3_se "weight-for-height >+3SD [standard error]"
label variable whz3_t "weight-for-height >+3SD [t-statistic]"
label variable whz3_pvalue "weight-for-height >+3SD [p value]"
label variable whz3_ll "weight-for-height >+3SD [95% lower confidence limit]"
label variable whz3_ul "weight-for-height >+3SD [95% upper confidence limit]"
label variable whz3_df "weight-for-height >+3SD [degrees of freedom]"
label variable whz3_crit "weight-for-height >+3SD [critical value]"
label variable whz3_weighted_N "weight-for-height >+3SD [weighted denominator]"
label variable whz3_unweighted_N "weight-for-height >+3SD [unweighted denominator]"
label variable whz3_weighted_Num "weight-for-height >+3SD [weighted numerator]"
label variable whz3_unweighted_Num "weight-for-height >+3SD [unweighted numerator]"
label variable baz_3_r "BMI-for-age <-3SD [point estimate]"
label variable baz_3_se "BMI-for-age <-3SD [standard error]"
label variable baz_3_t "BMI-for-age <-3SD [t-statistic]"
label variable baz_3_pvalue "BMI-for-age <-3SD [p value]"
label variable baz_3_ll "BMI-for-age <-3SD [95% lower confidence limit]"
label variable baz_3_ul "BMI-for-age <-3SD [95% upper confidence limit]"
label variable baz_3_df "BMI-for-age <-3SD [degrees of freedom]"
label variable baz_3_crit "BMI-for-age <-3SD [critical value]"
label variable baz_3_weighted_N "BMI-for-age <-3SD [weighted denominator]"
label variable baz_3_unweighted_N "BMI-for-age <-3SD [unweighted denominator]"
label variable baz_3_weighted_Num "BMI-for-age <-3SD [weighted numerator]"
label variable baz_3_unweighted_Num "BMI-for-age <-3SD [unweighted numerator]"
label variable baz_2_r "BMI-for-age <-2SD [point estimate]"
label variable baz_2_se "BMI-for-age <-2SD [standard error]"
label variable baz_2_t "BMI-for-age <-2SD [t-statistic]"
label variable baz_2_pvalue "BMI-for-age <-2SD [p value]"
label variable baz_2_ll "BMI-for-age <-2SD [95% lower confidence limit]"
label variable baz_2_ul "BMI-for-age <-2SD [95% upper confidence limit]"
label variable baz_2_df "BMI-for-age <-2SD [degrees of freedom]"
label variable baz_2_crit "BMI-for-age <-2SD [critical value]"
label variable baz_2_weighted_N "BMI-for-age <-2SD [weighted denominator]"
label variable baz_2_unweighted_N "BMI-for-age <-2SD [unweighted denominator]"
label variable baz_2_weighted_Num "BMI-for-age <-2SD [weighted numerator]"
label variable baz_2_unweighted_Num "BMI-for-age <-2SD [unweighted numerator]"
label variable baz_1_r "BMI-for-age <-1SD [point estimate]"
label variable baz_1_se "BMI-for-age <-1SD [standard error]"
label variable baz_1_t "BMI-for-age <-1SD [t-statistic]"
label variable baz_1_pvalue "BMI-for-age <-1SD [p value]"
label variable baz_1_ll "BMI-for-age <-1SD [95% lower confidence limit]"
label variable baz_1_ul "BMI-for-age <-1SD [95% upper confidence limit]"
label variable baz_1_df "BMI-for-age <-1SD [degrees of freedom]"
label variable baz_1_crit "BMI-for-age <-1SD [critical value]"
label variable baz_1_weighted_N "BMI-for-age <-1SD [weighted denominator]"
label variable baz_1_unweighted_N "BMI-for-age <-1SD [unweighted denominator]"
label variable baz_1_weighted_Num "BMI-for-age <-1SD [weighted numerator]"
label variable baz_1_unweighted_Num "BMI-for-age <-1SD [unweighted numerator]"
label variable baz1_r "BMI-for-age >+1SD [point estimate]"
label variable baz1_se "BMI-for-age >+1SD [standard error]"
label variable baz1_t "BMI-for-age >+1SD [t-statistic]"
label variable baz1_pvalue "BMI-for-age >+1SD [p value]"
label variable baz1_ll "BMI-for-age >+1SD [95% lower confidence limit]"
label variable baz1_ul "BMI-for-age >+1SD [95% upper confidence limit]"
label variable baz1_df "BMI-for-age >+1SD [degrees of freedom]"
label variable baz1_crit "BMI-for-age >+1SD [critical value]"
label variable baz1_weighted_N "BMI-for-age >+1SD [weighted denominator]"
label variable baz1_unweighted_N "BMI-for-age >+1SD [unweighted denominator]"
label variable baz1_weighted_Num "BMI-for-age >+1SD [weighted numerator]"
label variable baz1_unweighted_Num "BMI-for-age >+1SD [unweighted numerator]"
label variable baz2_r "BMI-for-age >+2SD [point estimate]"
label variable baz2_se "BMI-for-age >+2SD [standard error]"
label variable baz2_t "BMI-for-age >+2SD [t-statistic]"
label variable baz2_pvalue "BMI-for-age >+2SD [p value]"
label variable baz2_ll "BMI-for-age >+2SD [95% lower confidence limit]"
label variable baz2_ul "BMI-for-age >+2SD [95% upper confidence limit]"
label variable baz2_df "BMI-for-age >+2SD [degrees of freedom]"
label variable baz2_crit "BMI-for-age >+2SD [critical value]"
label variable baz2_weighted_N "BMI-for-age >+2SD [weighted denominator]"
label variable baz2_unweighted_N "BMI-for-age >+2SD [unweighted denominator]"
label variable baz2_weighted_Num "BMI-for-age >+2SD [weighted numerator]"
label variable baz2_unweighted_Num "BMI-for-age >+2SD [unweighted numerator]"
label variable baz3_r "BMI-for-age >+3SD [point estimate]"
label variable baz3_se "BMI-for-age >+3SD [standard error]"
label variable baz3_t "BMI-for-age >+3SD [t-statistic]"
label variable baz3_pvalue "BMI-for-age >+3SD [p value]"
label variable baz3_ll "BMI-for-age >+3SD [95% lower confidence limit]"
label variable baz3_ul "BMI-for-age >+3SD [95% upper confidence limit]"
label variable baz3_df "BMI-for-age >+3SD [degrees of freedom]"
label variable baz3_crit "BMI-for-age >+3SD [critical value]"
label variable baz3_weighted_N "BMI-for-age >+3SD [weighted denominator]"
label variable baz3_unweighted_N "BMI-for-age >+3SD [unweighted denominator]"
label variable baz3_weighted_Num "BMI-for-age >+3SD [weighted numerator]"
label variable baz3_unweighted_Num "BMI-for-age >+3SD [unweighted numerator]"
label variable haz_r "mean height-for-age  [point estimate]"
label variable haz_se "mean height-for-age  [standard error]"
label variable haz_t "mean height-for-age  [t-statistic]"
label variable haz_pvalue "mean height-for-age  [p value]"
label variable haz_ll "mean height-for-age  [95% lower confidence limit]"
label variable haz_ul "mean height-for-age  [95% upper confidence limit]"
label variable haz_df "mean height-for-age  [degrees of freedom]"
label variable haz_crit "mean height-for-age  [critical value]"
label variable haz_weighted_N "mean height-for-age  [weighted denominator]"
label variable haz_unweighted_N "mean height-for-age  [unweighted denominator]"
label variable haz_weighted_Num "mean height-for-age  [weighted numerator]"
label variable haz_unweighted_Num "mean height-for-age  [unweighted numerator]"
label variable waz_r "mean weight-for-age  [point estimate]"
label variable waz_se "mean weight-for-age  [standard error]"
label variable waz_t "mean weight-for-age  [t-statistic]"
label variable waz_pvalue "mean weight-for-age  [p value]"
label variable waz_ll "mean weight-for-age  [95% lower confidence limit]"
label variable waz_ul "mean weight-for-age  [95% upper confidence limit]"
label variable waz_df "mean weight-for-age  [degrees of freedom]"
label variable waz_crit "mean weight-for-age  [critical value]"
label variable waz_weighted_N "mean weight-for-age  [weighted denominator]"
label variable waz_unweighted_N "mean weight-for-age  [unweighted denominator]"
label variable waz_weighted_Num "mean weight-for-age  [weighted numerator]"
label variable waz_unweighted_Num "mean weight-for-age  [unweighted numerator]"
label variable whz_r "mean weight-for-height  [point estimate]"
label variable whz_se "mean weight-for-height  [standard error]"
label variable whz_t "mean weight-for-height  [t-statistic]"
label variable whz_pvalue "mean weight-for-height  [p value]"
label variable whz_ll "mean weight-for-height  [95% lower confidence limit]"
label variable whz_ul "mean weight-for-height  [95% upper confidence limit]"
label variable whz_df "mean weight-for-height  [degrees of freedom]"
label variable whz_crit "mean weight-for-height  [critical value]"
label variable whz_weighted_N "mean weight-for-height  [weighted denominator]"
label variable whz_unweighted_N "mean weight-for-height  [unweighted denominator]"
label variable whz_weighted_Num "mean weight-for-height  [weighted numerator]"
label variable whz_unweighted_Num "mean weight-for-height  [unweighted numerator]"
label variable baz_r "mean BMI-for-age  [point estimate]"
label variable baz_se "mean BMI-for-age  [standard error]"
label variable baz_t "mean BMI-for-age  [t-statistic]"
label variable baz_pvalue "mean BMI-for-age  [p value]"
label variable baz_ll "mean BMI-for-age  [95% lower confidence limit]"
label variable baz_ul "mean BMI-for-age  [95% upper confidence limit]"
label variable baz_df "mean BMI-for-age  [degrees of freedom]"
label variable baz_crit "mean BMI-for-age  [critical value]"
label variable baz_weighted_N "mean BMI-for-age  [weighted denominator]"
label variable baz_unweighted_N "mean BMI-for-age  [unweighted denominator]"
label variable baz_weighted_Num "mean BMI-for-age  [weighted numerator]"
label variable baz_unweighted_Num "mean BMI-for-age  [unweighted numerator]"
label variable haz_sd "mean height-for-age  [standard deviation]"
label variable waz_sd "mean weight-for-age  [standard deviation]"
label variable whz_sd "mean weight-for-height  [standard deviation]"
label variable baz_sd "mean BMI-for-age  [standard deviation]"

/// Category Placement

foreach x of numlist 1/21 {
set obs `x'
}

/// Total
replace standard_strat = "Total" in 1
replace strat_label = "sta_all" in 1
/// Sex
replace strat_label = "sta__sex" in 2
replace standard_strat = "Male" in 2
replace strat_label = "sta__sex" in 3
replace standard_strat = "Female" in 3
/// Age Groups
replace strat_label = "sta__agemonthgroup" in 4
replace standard_strat = "0 to 5 months" in 4
replace strat_label = "sta__agemonthgroup" in 5
replace standard_strat = "6 to 11 months" in 5
replace strat_label = "sta__agemonthgroup" in 6
replace standard_strat = "12 to 23 months" in 6
replace strat_label = "sta__agemonthgroup" in 7
replace standard_strat = "24 to 35 months" in 7
replace strat_label = "sta__agemonthgroup" in 8
replace standard_strat = "36 to 47 months" in 8
replace strat_label = "sta__agemonthgroup" in 9
replace standard_strat = "48 to 59 months" in 9
/// Age Groups x Sex
replace strat_label = "sta__agemonthgroup_sex" in 10
replace standard_strat = "Male 0 to 5 months" in 10
replace strat_label = "sta__agemonthgroup_sex" in 11
replace standard_strat = "Male 6 to 11 months" in 11
replace strat_label = "sta__agemonthgroup_sex" in 12
replace standard_strat = "Male 12 to 23 months" in 12
replace strat_label = "sta__agemonthgroup_sex" in 13
replace standard_strat = "Male 24 to 35 months" in 13
replace strat_label = "sta__agemonthgroup_sex" in 14
replace standard_strat = "Male 36 to 47 months" in 14
replace strat_label = "sta__agemonthgroup_sex" in 15
replace standard_strat = "Male 48 to 59 months" in 15
replace strat_label = "sta__agemonthgroup_sex" in 16
replace standard_strat = "Female 0 to 5 months" in 16
replace strat_label = "sta__agemonthgroup_sex" in 17
replace standard_strat = "Female 6 to 11 months" in 17
replace strat_label = "sta__agemonthgroup_sex" in 18
replace standard_strat = "Female 12 to 23 months" in 18
replace strat_label = "sta__agemonthgroup_sex" in 19
replace standard_strat = "Female 24 to 35 months" in 19
replace strat_label = "sta__agemonthgroup_sex" in 20
replace standard_strat = "Female 36 to 47 months" in 20
replace strat_label = "sta__agemonthgroup_sex" in 21
replace standard_strat = "Female 48 to 59 months" in 21

gen RUN=1
save "$prevfile", replace	

/// Analysis

 use "$prevfile", clear  
  local M=_N                                                  	// Number of lines in the dataset
  count if RUN==1
  local Mtoanalyse=r(N)                                       	// Number of surveys to analyse in the control dataset
  local crv=1  

  forval line = 1 / `M' {
    
    if RUN[`line']!= 1 continue
    

   
	foreach x in standard_strat strat_label {

	local `x'name = `x'[`line']
	}

	 *** Display survey being Prepped
    noisily display _newline(2) as input "Analysis: `crv' out of `Mtoanalyse' - `standard_stratname'" _newline
    capture confirm file "$workfile"
    if _rc > 0 {
      noisily display as error ">>>> Working file not available"
      use "$prevfile" , clear
      local ++crv
	  continue
   }
	
	/// Renaming Script

	foreach x in haz waz whz baz {
		foreach y in _3 _2 _1 1 2 3 {
	
		use "$workfile", clear				
	
		sum strata
  		local strata_check=r(N)
		
		sum cluster
		local cluster_check=r(N)
  
		if `strata_check'>0 & `cluster_check'>0 {
		svyset cluster [pw=sw], strata(strata) singleunit(centered)
		}
		
		if `strata_check'>0 & `cluster_check'==0 {
		svyset [pw=sw], strata(strata) singleunit(centered)
		}
		
		if `strata_check'==0 & `cluster_check'>0 {
		svyset cluster [pw=sw], singleunit(centered)
		}

		if `strata_check'==0 & `cluster_check'==0 {
		svyset [pw=sw], singleunit(centered)
		}
		
				
		sum `x'`y' if `x'`y'==1 & `strat_labelname'=="`standard_stratname'"
		local `x'`y'_Num=r(sum)
		
		sum `x'`y' [iw=sw] if `x'`y'==1 & `strat_labelname'=="`standard_stratname'"
		local `x'`y'_wNum=r(sum)
		
		sum `x'`y' if `strat_labelname'=="`standard_stratname'"
		local `x'`y'_D=r(N)
		
		sum `x'`y' [iw=sw] if `strat_labelname'=="`standard_stratname'"
		local `x'`y'_wD=r(sum_w)
		
		gen pop_int=0
		replace pop_int=1 if `strat_labelname'=="`standard_stratname'"
		replace pop_int=. if `strat_labelname'==""
		
	
		capture svy, subpop(pop_int):prop `x'`y'
		capture mat R = r(table)
		capture mat P = e(_N_subp)
		capture mat N = e(_N)
		
		capture mat RR = R["b".."crit",1...]'
		capture mat PP = P[1,1...]'
		capture mat NN = N[1,1...]'
		
		capture mat def out1 = (RR, PP, NN)
		
		clear
		capture {
		svmat out1
		
		drop in 1
		
		sum out11
		local `x'`y'_r_temp=r(sum)
		
		sum out12
		local `x'`y'_se_temp=r(sum)

		sum out13
		local `x'`y'_t_temp=r(sum)

		sum out14
		local `x'`y'_pvalue_temp=r(sum)
	
		sum out15
		local `x'`y'_ll_temp=r(sum)
	
		sum out16
		local `x'`y'_ul_temp=r(sum)
	
		sum out17
		local `x'`y'_df_temp=r(sum)
	
		sum out18
		local `x'`y'_crit_temp=r(sum)	
	
	
		}
		}
		}
			
		use "$prevfile", clear
		
		foreach x in haz waz whz baz {
				foreach y in _3 _2 _1 1 2 3 {
				foreach z in r se t pvalue ll ul df crit {
					replace `x'`y'_`z'=``x'`y'_`z'_temp' if strat_label=="`strat_labelname'" & standard_strat=="`standard_stratname'"
					replace `x'`y'_weighted_Num=``x'`y'_wNum' if strat_label=="`strat_labelname'" & standard_strat=="`standard_stratname'"
					replace `x'`y'_unweighted_Num=``x'`y'_Num' if strat_label=="`strat_labelname'" & standard_strat=="`standard_stratname'"
					replace `x'`y'_weighted_N=``x'`y'_wD' if strat_label=="`strat_labelname'" & standard_strat=="`standard_stratname'"
					replace `x'`y'_unweighted_N=``x'`y'_D' if strat_label=="`strat_labelname'" & standard_strat=="`standard_stratname'"
										}			
					}
					}
		
	foreach x in haz waz whz baz {
				foreach y in _3 _2 _1 1 2 3 {
				foreach z in r se t pvalue ll ul df crit {
					replace `x'`y'_`z'=0 if `x'`y'_unweighted_Num==0 & `x'`y'_unweighted_N==0
					replace `x'`y'_`z'=. if `x'`y'_unweighted_N==0
									}			
					}
					}	
	
		save "$prevfile", replace	
	
	/// For mean referencing
	
	foreach x in haz waz whz baz {

		use "$workfile", clear				
	
		sum strata
  		local strata_check=r(N)
		
		sum cluster
		local cluster_check=r(N)
  
		if `strata_check'>0 & `cluster_check'>0 {
		svyset cluster [pw=sw], strata(strata) singleunit(centered)
		}
		
		if `strata_check'>0 & `cluster_check'==0 {
		svyset [pw=sw], strata(strata) singleunit(centered)
		}
		
		if `strata_check'==0 & `cluster_check'>0 {
		svyset cluster [pw=sw], singleunit(centered)
		}

		if `strata_check'==0 & `cluster_check'==0 {
		svyset [pw=sw], singleunit(centered)
		}
		
		sum `x'`y' if `x'`y'==1 & `strat_labelname'=="`standard_stratname'"
		local `x'`y'_Num=r(sum)
		
		sum `x'`y' [iw=sw] if `x'`y'==1 & `strat_labelname'=="`standard_stratname'"
		local `x'`y'_wNum=r(sum)
		
		sum `x'`y' if `strat_labelname'=="`standard_stratname'"
		local `x'`y'_D=r(N)
		
		sum `x'`y' [iw=sw] if `strat_labelname'=="`standard_stratname'"
		local `x'`y'_wD=r(sum_w)
		
		gen pop_int=0
		replace pop_int=1 if `strat_labelname'=="`standard_stratname'"
		replace pop_int=. if `strat_labelname'==""
		
	
		capture svy, subpop(pop_int):mean `x'
		capture mat R = r(table)
		capture mat P = e(_N_subp)
		capture mat N = e(_N)
		capture svy, subpop(pop_int):mean `x'		
		capture estat sd, srssubpop
		capture mat O = r(sd)
		capture mat RR = R["b".."crit",1...]'
		capture mat PP = P[1,1...]'
		capture mat NN = N[1,1...]'
		
		
		
		capture mat def out1 = (RR,O)
		
		clear
		capture {
		svmat out1
		
		sum out11
		local `x'_r_temp=r(sum)
		
		sum out12
		local `x'_se_temp=r(sum)

		sum out13
		local `x'_t_temp=r(sum)

		sum out14
		local `x'_pvalue_temp=r(sum)
	
		sum out15
		local `x'_ll_temp=r(sum)
	
		sum out16
		local `x'_ul_temp=r(sum)
	
		sum out17
		local `x'_df_temp=r(sum)
	
		sum out18
		local `x'_crit_temp=r(sum)	
	
		sum out19
		local `x'_sd_temp=r(sum)
		}
		
		}
		
		use "$prevfile", clear
		
		foreach x in haz waz whz baz {
				foreach z in r se t pvalue ll ul df crit sd {
					replace `x'_`z'=``x'_`z'_temp' if strat_label=="`strat_labelname'" & standard_strat=="`standard_stratname'"
					replace `x'_weighted_Num=``x'_wNum' if strat_label=="`strat_labelname'" & standard_strat=="`standard_stratname'"
					replace `x'_unweighted_Num=``x'_Num' if strat_label=="`strat_labelname'" & standard_strat=="`standard_stratname'"
					replace `x'_weighted_N=``x'_wD' if strat_label=="`strat_labelname'" & standard_strat=="`standard_stratname'"
					replace `x'_unweighted_N=``x'_D' if strat_label=="`strat_labelname'" & standard_strat=="`standard_stratname'"
								
					}
					}	
		
				foreach x in haz waz whz baz {
				foreach z in r se t pvalue ll ul df crit sd {
					replace `x'_`z'=0 if `x'_unweighted_Num==0 & `x'_unweighted_N==0	
					replace `x'_`z'=. if `x'_unweighted_N==0	
					}
					}
		
		save "$prevfile", replace	
		

		use "$prevfile" , clear
		local ++crv	
	}
 
  /// For next disaggregation
	
/// Cleanup
		

	
	use "$prevfile" , clear
	drop strat_label RUN

	save "$finalfile", replace
	export delimited using "$exportfile", replace
	
	erase "$prevfile"
	erase "$workfile"
	

	noisily display _newline(2) as input "Note 3:	Prevalence estimates are written to " _newline
	noisily display _newline(2) as input "		$finalfile" _newline
	noisily display _newline(2) as input "Note 2: Prevalence estimates are written to" _newline
	noisily display _newline(2) as input "		$exportfile" _newline
	
		use "$outfile", clear
		clear matrix
		clear mata
		
/// Run programme successfully
}
}
end 


exit
