*! mixl_d0 1.0.1 28Mar2007
*! author arh

program define mixlog_d0 
	version 9.2
	args todo b lnf

	if "$mixl_panid" != "" {
		mata: mixl_ll_pan("`b'")	
	}
	else {
		mata: mixl_ll_cross("`b'")
	}

	scalar `lnf' = r(ll)
end

version 9.2
mata: 
void mixl_ll_cross(string scalar B_s)
{
	external _mixl_X
	external _mixl_Y
	external _mixl_CSID

	external _mixl_nrep
	external _mixl_nobs
	external _mixl_kfix
	external _mixl_krnd
	external _mixl_krln
	external _mixl_burn
	external _mixl_corr

	B = st_matrix(B_s)'

	k = _mixl_kfix + _mixl_krnd

	if (_mixl_kfix > 0) {
		MFIX = B[|1,1\_mixl_kfix,1|]
		MFIX = MFIX :* J(_mixl_kfix,_mixl_nrep,1)	
	}

	MRND = B[|(_mixl_kfix+1),1\k,1|]

	if (_mixl_corr == 1) {
		external _mixl_ncho
		SRND = invvech(B[|(k+1),1\(k+_mixl_ncho),1|]) :* lowertriangle(J(_mixl_krnd,_mixl_krnd,1))
	}
	else {
		SRND = diag(B[|(k+1),1\(_mixl_kfix+2*_mixl_krnd),1|])
	}

	P = J(_mixl_nobs,1,0)

	i = 1
	for (n=1; n<=_mixl_nobs; n++) {
		ERR = invnormal(halton(_mixl_nrep,_mixl_krnd,(1+_mixl_burn+_mixl_nrep*(n-1)))')
		if (_mixl_kfix > 0) BETA = MFIX \ (MRND :+ (SRND*ERR))
		else BETA = MRND :+ (SRND*ERR)
		if (_mixl_krln > 0) {
			if ((k-_mixl_krln) > 0) { 
				BETA = BETA[|1,1\(k-_mixl_krln),_mixl_nrep|] \ ///
				exp(BETA[|(k-_mixl_krln+1),1\k,_mixl_nrep|])
			}
			else {
				BETA = exp(BETA)
			}
		}
		EV = _mixl_X[|i,1\(i+_mixl_CSID[i,1]-1),cols(_mixl_X)|]
		EV = exp(EV*BETA)
		R = colsum(EV :* _mixl_Y[|i,1\(i+_mixl_CSID[i,1]-1),cols(_mixl_Y)|]) :/ colsum(EV)
		P[n,1] = mean(R',1)	
		i = i + _mixl_CSID[i,1]
	}
	st_numscalar("r(ll)", colsum(ln(P)))	
}
end

version 9.2
mata: 
void mixl_ll_pan(string scalar B_s)
{
	external _mixl_X
	external _mixl_Y
	external _mixl_T
	external _mixl_CSID

	external _mixl_nrep
	external _mixl_np
	external _mixl_kfix
	external _mixl_krnd
	external _mixl_krln
	external _mixl_burn
	external _mixl_corr

	B = st_matrix(B_s)'

	k = _mixl_kfix + _mixl_krnd

	if (_mixl_kfix > 0) {
		MFIX = B[|1,1\_mixl_kfix,1|]
		MFIX = MFIX :* J(_mixl_kfix,_mixl_nrep,1)	
	}

	MRND = B[|(_mixl_kfix+1),1\k,1|]

	if (_mixl_corr == 1) {
		external _mixl_ncho
		SRND = invvech(B[|(k+1),1\(k+_mixl_ncho),1|]) :* lowertriangle(J(_mixl_krnd,_mixl_krnd,1))
	}
	else {
		SRND = diag(B[|(k+1),1\(_mixl_kfix+2*_mixl_krnd),1|])
	}

	P = J(_mixl_np,1,0)

	i = 1
	for (n=1; n<=_mixl_np; n++) {
		ERR = invnormal(halton(_mixl_nrep,_mixl_krnd,(1+_mixl_burn+_mixl_nrep*(n-1)))')
		if (_mixl_kfix > 0) BETA = MFIX \ (MRND :+ (SRND*ERR))
		else BETA = MRND :+ (SRND*ERR)
		if (_mixl_krln > 0) {
			if ((k-_mixl_krln) > 0) {
				BETA = BETA[|1,1\(k-_mixl_krln),_mixl_nrep|] \ ///
				exp(BETA[|(k-_mixl_krln+1),1\k,_mixl_nrep|])
			}
			else {
				BETA = exp(BETA)
			}
		}
		R = J(1,_mixl_nrep,1)
		t = 1
		nc = _mixl_T[i,1]
		for (t=1; t<=nc; t++) {
			EV = _mixl_X[|i,1\(i+_mixl_CSID[i,1]-1),cols(_mixl_X)|]
			EV = exp(EV*BETA)
			EV = colsum(EV :* _mixl_Y[|i,1\(i+_mixl_CSID[i,1]-1),cols(_mixl_Y)|]) :/ colsum(EV)
			R = R :* EV
			i = i + _mixl_CSID[i,1]
		}
		P[n,1] = mean(R',1)
	}
	st_numscalar("r(ll)", colsum(ln(P)))	
}
end	

exit




			


