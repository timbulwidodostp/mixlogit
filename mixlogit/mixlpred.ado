*! mixlpred 1.0.1 02Jun2007
*! author arh

program define mixlpred
	version 9.2

	syntax newvarname [if] [in], [NREP(integer 50) BURN(integer 15)]

	if ("`e(cmd)'" != "mixlogit") error 301

	** Mark the prediction sample **
	marksample touse, novarlist
	markout `touse' `e(indepvars)' `e(group)' `e(id)'

	** Generate variables used to sort data **
	tempvar sorder altid
	gen `sorder' = _n
	sort `touse' `e(id)' `e(group)'
	by `touse' `e(id)' `e(group)': gen `altid' = _n 
	
	** Drop data not in prediction sample **
	preserve
	qui keep if `touse'

	** Generate individual id **
	if ("`e(id)'" != "") {
		tempvar nchoice pid
		sort `e(group)'
		by `e(group)': gen `nchoice' = cond(_n==_N,1,0)
		sort `e(id)'
		by `e(id)': egen `pid' = sum(`nchoice')		
		qui duplicates report `e(id)'
		local np = r(unique_value)
		mata: _mixl_T = st_data(., ("`pid'"))
	}

	** Generate choice occacion id **
	tempvar csid
	sort `e(group)'
	by `e(group)': egen `csid' = sum(1)
	qui duplicates report `e(group)'
	local nobs = r(unique_value)

	** Sort data **
	sort `e(id)' `e(group)' `altid'

	** Set Mata matrices to be used in optimisation routine **
	local rhs `e(indepvars)'
	mata: _mixl_X = st_data(., tokens(st_local("rhs")))
	mata: _mixl_CSID = st_data(., ("`csid'"))
	local totobs = _N	

	** Restore data **
	restore
	
	tempname b
	matrix `b' = e(b)
	
	qui gen double `varlist' = .

	if ("`e(id)'" == "") mata: mixlpred_cross("`b'", "`varlist'", "`touse'")	
	else mata: mixlpred_pan("`b'", "`varlist'", "`touse'")

	** Restore sort order **
	sort `sorder'	
end

version 9.2
mata: 
void mixlpred_cross(string scalar B_s, string scalar P_s, string scalar TOUSE_s)
{
	external _mixl_X
	external _mixl_CSID

	nrep = strtoreal(st_local("nrep"))
	nobs = strtoreal(st_local("nobs"))
	totobs = strtoreal(st_local("totobs"))
	kfix = st_numscalar("e(kfix)")
	krnd = st_numscalar("e(krnd)")
	krln = st_numscalar("e(krln)")
	burn = strtoreal(st_local("burn"))
	corr = st_numscalar("e(corr)")
		
	B = st_matrix(B_s)'

	k = kfix + krnd

	if (kfix > 0) {
		MFIX = B[|1,1\kfix,1|]
		MFIX = MFIX :* J(kfix,nrep,1)	
	}

	MRND = B[|(kfix+1),1\k,1|]

	if (corr == 1) {
		ncho = st_numscalar("e(k_aux)")
		SRND = invvech(B[|(k+1),1\(k+ncho),1|]) :* lowertriangle(J(krnd,krnd,1))
	}
	else {
		SRND = diag(B[|(k+1),1\(kfix+2*krnd),1|])
	}

	P = J(totobs,1,0)

	i = 1
	for (n=1; n<=nobs; n++) {
		ERR = invnormal(halton(nrep,krnd,(1+burn+nrep*(n-1)))')
		if (kfix > 0) BETA = MFIX \ (MRND :+ (SRND*ERR))
		else BETA = MRND :+ (SRND*ERR)
		if (krln > 0) {
			if ((k-krln) > 0) { 
				BETA = BETA[|1,1\(k-krln),nrep|]\exp(BETA[|(k-krln+1),1\k,nrep|])
			}
			else {
				BETA = exp(BETA)
			}
		}
		EV = _mixl_X[|i,1\(i+_mixl_CSID[i,1]-1),cols(_mixl_X)|]
		EV = exp(EV*BETA)		
		R = EV :/ colsum(EV)
		P[|i,1\(i+_mixl_CSID[i,1]-1),1|] = mean(R',1)'	
		i = i + _mixl_CSID[i,1]
	}
	st_store(.,P_s,TOUSE_s,P)
}
end

version 9.2
mata: 
void mixlpred_pan(string scalar B_s, string scalar P_s, string scalar TOUSE_s)
{
	external _mixl_X
	external _mixl_T
	external _mixl_CSID

	nrep = strtoreal(st_local("nrep"))
	np = strtoreal(st_local("np"))
	totobs = strtoreal(st_local("totobs"))
	kfix = st_numscalar("e(kfix)")
	krnd = st_numscalar("e(krnd)")
	krln = st_numscalar("e(krln)")
	burn = strtoreal(st_local("burn"))
	corr = st_numscalar("e(corr)")

	B = st_matrix(B_s)'

	k = kfix + krnd

	if (kfix > 0) {
		MFIX = B[|1,1\kfix,1|]
		MFIX = MFIX :* J(kfix,nrep,1)	
	}

	MRND = B[|(kfix+1),1\k,1|]

	if (corr == 1) {
		ncho = st_numscalar("e(k_aux)")
		SRND = invvech(B[|(k+1),1\(k+ncho),1|]) :* lowertriangle(J(krnd,krnd,1))
	}
	else {
		SRND = diag(B[|(k+1),1\(kfix+2*krnd),1|])
	}

	P = J(totobs,1,0)

	i = 1
	for (n=1; n<=np; n++) {
		ERR = invnormal(halton(nrep,krnd,(1+burn+nrep*(n-1)))')
		if (kfix > 0) BETA = MFIX \ (MRND :+ (SRND*ERR))
		else BETA = MRND :+ (SRND*ERR)
		if (krln > 0) {
			if ((k-krln) > 0) { 
				BETA = BETA[|1,1\(k-krln),nrep|]\exp(BETA[|(k-krln+1),1\k,nrep|])
			}
			else {
				BETA = exp(BETA)
			}
		}
		t = 1
		nc = _mixl_T[i,1]
		for (t=1; t<=nc; t++) {
			EV = _mixl_X[|i,1\(i+_mixl_CSID[i,1]-1),cols(_mixl_X)|]
			EV = exp(EV*BETA)
			R = EV :/ colsum(EV)
			P[|i,1\(i+_mixl_CSID[i,1]-1),1|] = mean(R',1)'
			i = i + _mixl_CSID[i,1]
		}
	}
	st_store(.,P_s,TOUSE_s,P)	
}
end	

exit

