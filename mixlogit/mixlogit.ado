*! mixlogit 1.0.0 5Dec2006
*! author arh

program mixlogit
	version 9.2
	if replay() {
		if (`"`e(cmd)'"' != "mixlogit") error 301
		Replay `0'
	}
	else	Estimate `0'
end

program Estimate, eclass
	syntax varlist [if] [in], 		///
		GRoup(varname) 			///
		RAND(varlist) [			///
		ID(varname) 			///
		LN(integer 0) 			///
		CORR					///
		NREP(integer 50)			///
		BURN(integer 15)			///
		FRom(string)			/// 
		Level(integer `c(level)')	///
		TRace					///
		GRADient				///
		HESSian				///
		SHOWSTEP				///
		ITERate(passthru)			///
		TOLerance(passthru)		///
		LTOLerance(passthru)		///
		GTOLerance(passthru)		///
		NRTOLerance(passthru)		///
		CONSTraints(passthru)		///
		TECHnique(passthru)		///
		DIFficult				///
	]

	local mlopts `trace' `gradient' `hessian' `showstep' `iterate' `tolerance' ///
	`ltolerance' `gtolerance' `nrtolerance' `constraints' `technique' `difficult'

	if ("`technique'" == "technique(bhhh)") {
		di in red "technique(bhhh) is not allowed."
		exit
	}

	** Mark the estimation sample **
	marksample touse
	markout `touse' `group' `rand' `id'

	** Estimate conditional logit model **
	gettoken lhs fixed : varlist
	local rhs `fixed' `rand'
	qui clogit `lhs' `rhs' if `touse', group(`group')
	local ll = e(ll)
	local k  = e(k)
	qui replace `touse' = e(sample)

	** Drop missing data **
	preserve
	qui keep if `touse'

	** Check that no variables have been specified to have both fixed and random coefficients **
	qui _rmcoll `rhs' 
	if "`r(varlist)'" != "`rhs'" {
		di in red "At least one variable has been specified to have both fixed and random coefficients"
		exit 498
	}

	** Generate individual id **
	if ("`id'" != "") {
		tempvar nchoice pid
		sort `group'
		by `group': gen `nchoice' = cond(_n==_N,1,0)
		sort `id'
		by `id': egen `pid' = sum(`nchoice')		
		qui duplicates report `id'
		mata: _mixl_np = st_numscalar("r(unique_value)")
		mata: _mixl_T = st_data(., ("`pid'"))
	}

	global mixl_panid `pid'

	** Generate choice occacion id **
	tempvar csid
	sort `group'
	by `group': egen `csid' = sum(1)
	qui duplicates report `group'
	mata: _mixl_nobs = st_numscalar("r(unique_value)")

	** Sort data **
	sort `id' `group'

	** Set Mata matrices and scalars to be used in optimisation routine **
	local kfix: word count `fixed'
	local krnd: word count `rand'

	mata: _mixl_X = st_data(., tokens(st_local("rhs")))
	mata: _mixl_Y = st_data(., ("`lhs'"))
	mata: _mixl_CSID = st_data(., ("`csid'"))

	mata: _mixl_nrep = strtoreal(st_local("nrep"))
	mata: _mixl_kfix = strtoreal(st_local("kfix"))
	mata: _mixl_krnd = strtoreal(st_local("krnd"))
	mata: _mixl_krln = strtoreal(st_local("ln"))
	mata: _mixl_burn = strtoreal(st_local("burn"))
	mata: _mixl_totobs = st_nobs()

	** Restore data **
	restore

	** Create macro to define equations for optimisation routine **
	local mean (Mean: `rhs', noconst)
	if ("`corr'" == "") {
		mata: _mixl_corr = 0
		local sd (SD: `rand', noconst)
		local max `mean' `sd' 
	}
	else {
		mata: _mixl_corr = 1
		local cho = `krnd'*(`krnd'+1)/2
		mata: _mixl_ncho = strtoreal(st_local("cho"))
		local max `mean'
		forvalues i = 1(1)`krnd' {
			forvalues j = `i'(1)`krnd' {
				local max `max' /l`j'`i'
			}
		}
	}

	** Create matrix of starting values unless specified **
	if ("`from'" == "") {
		tempname b from
		matrix `b' = e(b)
		if ((`kfix'+`krnd')>`ln') matrix `from' = `b'[1,1..(`kfix'+`krnd'-`ln')]
		forvalues i = 1(1)`ln' {
			if (`b'[1,(`kfix'+`krnd'-`ln'+`i')] <= 0) {
				di in red "Variables specified to have log-normally distributed coefficients should have positive"
				di in red "coefficients in the conditional logit model. Try multiplying the variable by -1."
				exit 498
			}
			if ((`kfix'+`krnd')==`ln' & `i'==1) matrix `from' = ln(`b'[1,1])
			else matrix `from' = `from', ln(`b'[1,(`kfix'+`krnd'-`ln'+`i')])
		} 
		if ("`corr'" == "") matrix `from' = `from', J(1,`krnd',0.1)
		else matrix `from' = `from', J(1,`cho',0.1)
		local copy , copy
	}

	** Run optimisation routine **
	ml model d0 mixlog_d0 							///
		`max' if `touse', search(off) init(`from' `copy') 	///
		`mlopts' maximize	lf0(`k' `ll') missing nopreserve			


	** To be returned as e() **
	ereturn local title "Mixed logit model"
	ereturn local indepvars `rhs'
	ereturn local depvar `lhs'
	ereturn local group `group'
	ereturn scalar kfix = `kfix'
	ereturn scalar krnd = `krnd'
	ereturn scalar krln = `ln'
	ereturn scalar nrep = `nrep'
	ereturn scalar burn = `burn'
	if ("`corr'" != "") {
		ereturn scalar corr = 1
		ereturn scalar k_aux = `cho'
	}
	else ereturn scalar corr = 0
	if ("`id'" != "") ereturn local id `id'
	ereturn local cmd "mixlogit"

	Replay , level(`level')
end

program Replay
	syntax [, Level(integer `c(level)') ]
	ml display , level(`level')
end

exit


