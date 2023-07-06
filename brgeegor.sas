/* SAS macro: brgeegor.sas                                                   */
/* by Yukio Tada                                                */

/* id	    - 	An ID variable indicating cluster (person)		     */
/* y	    - 	The category response taking on values 1,...,K		     */
/* x	    - 	The covariates which model the categorical response	     */
/* beta0    -   Initial value of regression parameter                        */
/* criteria -	The convergence criterion, which is the maximum absolute     */
/*		value of the difference of parameter estimates between	     */
/*		iterations.						     */
/* maxit    -	The maximum number of iterations allowed.		     */
/* assoc    -   The working covariance matrix: 1=independent, 2=exchangeable */
/* mgee     -   modified GEE based on Firth(1993): 1=No, 2=Yes               */

proc iml worksize=5000;
  reset nolog noprint;

start brgeegor(ID,TP,X,Y,BETA0,CRITERIA,MAXIT,ASSOC,MGEE);

  n = max(id);             /* total number of persons */

  maxy = max(y);                      /* # levels of multinomial */
  maxu = maxy * maxy;
  stop = maxy - 1;
  stop2 = stop * stop;

  trans = j(stop,stop,1);    /* cumulative transform matrix */
  do r = 1 to stop;
    do c = 1 to stop;
      if c > r then do;
      trans[r,c]=0;
      end;
    end;
  end;
  transinv = inv(trans);

/* set initial value of parameter */
  beta = j(stop,1,0)//beta0;
  nbeta = nrow(beta);

  if nrow(loc(beta0^=0)) = 0 then do;
    ysize = nrow(y);
    do yr = 1 to stop;
      cum_p=nrow(y[loc(y<=yr),])/ysize;
      beta[yr,1]=log(cum_p/(1-cum_p));
    end;
  end;

/* To estimate crude odds ratio */
  maxt = max(tp);      /* max of time points */

/**/
if ASSOC=2 then do;
  npair=0;
  twLOR=0;
  mwLOR=0;
  do t1=1 to maxt-1;
    do t2=t1+1 to maxt;
	  wLOR=0;
	  tVarOR=0;

      id_ = id[loc(tp=t1|tp=t2),];
      tp_ = tp[loc(tp=t1|tp=t2),];
      y_ = y[loc(tp=t1|tp=t2),];

	  do L=1 to stop;
        do M=1 to stop;
	      A_LM=0;
          B_LM=0;
          C_LM=0;
          D_LM=0;
	      do i=1 to n;
	        pair=ncol(loc(id_=i));
		    if pair=2 then do;
		      Y1=Y_[loc(ID_=i & tp_=t1),];
		      Y2=Y_[loc(ID_=i & tp_=t2),];
		      A_LM=A_LM+(Y1>L & Y2>M);
		      B_LM=B_LM+(Y1<=L & Y2>M);
		      C_LM=C_LM+(Y1>L & Y2<=M);
		      D_LM=D_LM+(Y1<=L & Y2<=M);
		    end;
	      end;
          LOR_LM=log(A_LM+0.5)+log(D_LM+0.5)-log(B_LM+0.5)-log(C_LM+0.5);
	      VarOR_LM=1/(A_LM+0.5)+1/(B_LM+0.5)+1/(C_LM+0.5)+1/(D_LM+0.5);
	      wLOR=wLOR+LOR_LM/VarOR_LM;
	      tVarOR=tVarOR+1/VarOR_LM;
	    end;
      end;
      wLOR=wLOR/tVarOR;
      npair=npair+1;
      twLOR=twLOR+wLOR;
    end;
  end;
  mwLOR=twLOR/npair;
end;
/**/

crit=1;

Do it=1 to maxit while (crit > criteria);

  U1 = j( nbeta, 1, 0);
  dvd = j( nbeta, nbeta, 0);
  u1sq = j( nbeta, nbeta, 0);
  u1sq_c = j( nbeta, nbeta, 0);
  u1_adj = j( nbeta, 1, 0);
  n_obs = 0;

Do i=1 to n;
  free Y_ie Y_i X_ie X_i mp_ie p_ie p_i A_ie A_i D_ie D_i
    PSI_i Gam_i1 Gam_i2 S_ijk F_ijk PP_ijk Vjk_i V1j_i V1_i V_i;
    times = loc(id=i);
    T_i = ncol(times);
    Yvar_i = Y[times,];
    Xvar_i = X[times,];
	A_i = 0;

  IF T_i > 1 then do;
    Do j=1 to T_i;
      Y_ie = j(maxy,1,0);
      Y_ie[Yvar_i[j,], 1] = 1;
      Y_ie = Y_ie[1:stop,1];
      Y_i = Y_i // Y_ie;

      X_ie = j(stop,1,1) @ Xvar_i[j,];
      X_ie = I(stop) || X_ie;
      X_i = X_i // X_ie;

      mp_ie = logistic(X_ie * beta);

      p_ie = transinv * mp_ie;
      p_i = p_i // p_ie;

/**/
	  A_ie = diag(p_ie) - p_ie*p_ie`;
	  A_i = Block(A_i, A_ie);
/**/

      D_ie = X_ie` * Diag( mp_ie#(1-mp_ie) ) * transinv`;
      D_i = D_i || D_ie;
	END;

    nrA = nrow(A_i);
    A_i = A_i[2:nrA,2:nrA]; 

/**/
	if ASSOC=2 then do;
    Do j=1 to T_i-1;
	  free V1j_i;
	  Do k=j+1 to T_i;
	    Gam_i1 = cusum(p_i[(j-1)*stop+1:j*stop,1]);
		Gam_i2 = cusum(p_i[(k-1)*stop+1:k*stop,1]);
		Gam_i1 = Gam_i1[1:stop,1];
		Gam_i2 = Gam_i2[1:stop,1];

        PSI_i = exp(mwLOR);

        S_ijk = J( stop, stop, 0);
        F_ijk  = J( maxy, maxy, 0);
        PP_ijk = J( stop, stop, 0);
        Vjk_i = J(stop, stop, 0);
        do L=1 to stop;
          do M=1 to stop;
            S_ijk[L,M] = sqrt(((1+(Gam_i1[L,1]+Gam_i2[M,1])*(PSI_i-1))**2)+(4*PSI_i*(1-PSI_i)*Gam_i1[L,1]*Gam_i2[M,1]));
            if PSI_i^=1 then F_ijk[L,M] = (1+(Gam_i1[L,1]+Gam_i2[M,1])*(PSI_i-1)-S_ijk[L,M])/(2*(PSI_i-1));
            else F_ijk[L,M] = Gam_i1[L,1]*Gam_i2[M,1];
                 if L=1 & M=1 then PP_ijk[L,M] = F_ijk[L,M];
            else if L=1 & M>1 then PP_ijk[L,M] = F_ijk[L,M] - F_ijk[L,M-1];
            else if L>1 & M=1 then PP_ijk[L,M] = F_ijk[L,M] - F_ijk[L-1,M];
            else if L>1 & M>1 then PP_ijk[L,M] = F_ijk[L,M] - F_ijk[L,M-1] - F_ijk[L-1,M] + F_ijk[L-1,M-1];
			Vjk_i[L,M] = PP_ijk[L,M] - p_i[(j-1)*stop+L,1] * p_i[(k-1)*stop+M,1];
          end;
        end;
        V1j_i = V1j_i || Vjk_i;
	  END;
      V1j_i = J(stop, j*stop, 0) || V1j_i; 
	  V1_i = V1_i // V1j_i;
    END; 
	V1_i = V1_i // J(stop, T_i*stop, 0);
	V_i = A_i + V1_i + V1_i`;
	end;
	else do;
	V_i = A_i;
	end;

  END;
  
  dvd =  dvd + D_i*ginv(V_i)*D_i`;

  n_obs = n_obs + T_i;

end;               * this ends the loop over each person for B;

Do i=1 to n;
  free Y_ie Y_i X_ie X_i mp_ie mp_i p_ie p_i A_ie A_i D_ie D_i DE_ie DE_i
    PSI_i Gam_i1 Gam_i2 S_ijk F_ijk PP_ijk Vjk_i V1j_i V1_i V_i u1_i;
    times = loc(id=i);
    T_i = ncol(times);
    Yvar_i = Y[times,];
    Xvar_i = X[times,];
	A_i = 0;
	DE_i = 0;

  IF T_i > 1 then do;
    Do j=1 to T_i;
      Y_ie = j(maxy,1,0);
      Y_ie[Yvar_i[j,], 1] = 1;
      Y_ie = Y_ie[1:stop,1];
      Y_i = Y_i // Y_ie;

      X_ie = j(stop,1,1) @ Xvar_i[j,];
      X_ie = I(stop) || X_ie;
      X_i = X_i // X_ie;

      mp_ie = logistic(X_ie * beta);
	  mp_i = mp_i // mp_ie;

      p_ie = transinv * mp_ie;
      p_i = p_i // p_ie;

/**/
	  A_ie = diag(p_ie) - p_ie*p_ie`;
	  A_i = Block(A_i, A_ie);
/**/

      D_ie = X_ie` * Diag( mp_ie#(1-mp_ie) ) * transinv`;
      D_i = D_i || D_ie;

/**/
	  DE_ie = Diag(mp_ie#(1-mp_ie)) * transinv`;
	  DE_i = Block(DE_i,DE_ie);
/**/
	END;

    nrA = nrow(A_i);
    A_i = A_i[2:nrA,2:nrA]; 

	nrDE = nrow(DE_i);
	DE_i = DE_i[2:nrDE,2:nrDE];

/**/
	if ASSOC=2 then do;
    Do j=1 to T_i-1;
	  free V1j_i;
	  Do k=j+1 to T_i;
	    Gam_i1 = cusum(p_i[(j-1)*stop+1:j*stop,1]);
		Gam_i2 = cusum(p_i[(k-1)*stop+1:k*stop,1]);
		Gam_i1 = Gam_i1[1:stop,1];
		Gam_i2 = Gam_i2[1:stop,1];

        PSI_i = exp(mwLOR);

        S_ijk = J( stop, stop, 0);
        F_ijk  = J( maxy, maxy, 0);
        PP_ijk = J( stop, stop, 0);
        Vjk_i = J(stop, stop, 0);
        do L=1 to stop;
          do M=1 to stop;
            S_ijk[L,M] = sqrt(((1+(Gam_i1[L,1]+Gam_i2[M,1])*(PSI_i-1))**2)+(4*PSI_i*(1-PSI_i)*Gam_i1[L,1]*Gam_i2[M,1]));
            if PSI_i^=1 then F_ijk[L,M] = (1+(Gam_i1[L,1]+Gam_i2[M,1])*(PSI_i-1)-S_ijk[L,M])/(2*(PSI_i-1));
            else F_ijk[L,M] = Gam_i1[L,1]*Gam_i2[M,1];
                 if L=1 & M=1 then PP_ijk[L,M] = F_ijk[L,M];
            else if L=1 & M>1 then PP_ijk[L,M] = F_ijk[L,M] - F_ijk[L,M-1];
            else if L>1 & M=1 then PP_ijk[L,M] = F_ijk[L,M] - F_ijk[L-1,M];
            else if L>1 & M>1 then PP_ijk[L,M] = F_ijk[L,M] - F_ijk[L,M-1] - F_ijk[L-1,M] + F_ijk[L-1,M-1];
			Vjk_i[L,M] = PP_ijk[L,M] - p_i[(j-1)*stop+L,1] * p_i[(k-1)*stop+M,1];
          end;
        end;
        V1j_i = V1j_i || Vjk_i;
	  END;
      V1j_i = J(stop, j*stop, 0) || V1j_i; 
	  V1_i = V1_i // V1j_i;
    END; 
	V1_i = V1_i // J(stop, T_i*stop, 0);
	V_i = A_i + V1_i + V1_i`;
	end;
	else do;
	V_i = A_i;
	end;

  END;

/**/
  dd_de_i=J(stop#T_i#stop#T_i, stop#T_i, 0);
  do ss=1 to T_i;
    do tt=1 to stop;
	  do uu=1 to stop;
	    ii=(ss-1)#stop+tt;
		jj=(ss-1)#stop+uu;
		kk=(ii-1)#stop#T_i+jj;
		if ii=jj then dd_de_i[kk,ii]=(1-2#mp_i[ii])#(1-mp_i[ii])#mp_i[ii];
		if ii=jj+1 then dd_de_i[kk,jj]=-(1-2#mp_i[jj])#(1-mp_i[jj])#mp_i[jj];
	  end;
	end;
  end;

  U1_adj_i=J(nbeta,1,0);
  DE_V_i=DE_i*ginv(V_i);
  do r=1 to nbeta;
    do s=1 to stop#T_i;
	  U1_adj_i[r]=U1_adj_i[r]+0.5*trace(X_i*ginv(dvd)*X_i`*(DE_V_i[s,]@I(stop#T_i))*dd_de_i)*X_i[s,r];
	end;
  end;

  u1_i = D_i*ginv(V_i)*( Y_i - p_i );

  U1 = U1 + u1_i;                * the first set of equations;

  u1sq = u1sq + u1_i*u1_i`;      * estimate of sigma 1,1 matrix;

  U1_adj = U1_adj + U1_adj_i;

/**/
/*   h_i = D_i`*ginv(dvd)*D_i*ginv(V_i); */
/*   u1sq_c = u1sq_c + D_i*ginv(V_i)*inv(I((stop) # T_i)-h_i)*( Y_i - p_i )*( Y_i - p_i )`*inv(I((stop) # T_i)-h_i`)*ginv(V_i)*D_i`; */
/**/

end;               * this ends the loop over each person for B;

    if MGEE=2 then U1 = U1 + U1_adj;

    DELTA1 = ginv(dvd)*U1;

    beta = beta + DELTA1;
    crit = max( ABS(DELTA1));

end;                       * this ends the loop over each iteration;

    vb0 = ginv(dvd);

/*robust variance*/
	vb = vb0 * u1sq * vb0;
/*Mancl and DeRouen (2001)*/
    *vb_c = vb0 * u1sq_c * vb0;    *robust variance matrix of beta;

    dvd2 = dvd + U1_adj*U1_adj`;
    u1sq2 = u1sq + U1_adj*U1_adj`;
    if MGEE=1 then do;
      vb02 = vb0;
      vb2 = vb;
    end;
    else if MGEE=2 then do;
      vb02 = ginv(dvd2);
      vb2 = vb02 * u1sq2 * vb02;
    end;


   conv1 = j(nbeta,1,0);
   conv2 = j(nbeta,1,0);
   estimat1 = j(nbeta,1,.);
   se_est1 = j(nbeta,1,.);
   se_est2 = j(nbeta,1,.);
   z1 = j(nbeta,1,.);
   z2 = j(nbeta,1,.);
   p1 = j(nbeta,1,.);
   lower1 = j(nbeta,1,.);
   upper1 = j(nbeta,1,.);
   p2 = j(nbeta,1,.);
   lower2 = j(nbeta,1,.);
   upper2 = j(nbeta,1,.);
   p3 = j(nbeta,1,.);
   lower3 = j(nbeta,1,.);
   upper3 = j(nbeta,1,.);
   p4 = j(nbeta,1,.);
   lower4 = j(nbeta,1,.);
   upper4 = j(nbeta,1,.);
 
   if crit <= criteria then do;
/**/
     if min(vecdiag(vb2)) > -1e-8 then do; 
        conv1 = j(nbeta,1,1); 
        se_est1 = sqrt(choose(vecdiag(vb2) > 0, vecdiag(vb2), 0)); 
 	 end; 
	 else do;
       se_est1 = sqrt(vecdiag(vb02));
	 end;
     if min(vecdiag(vb)) > -1e-8 then do; 
        conv2 = j(nbeta,1,1); 
        se_est2 = sqrt(choose(vecdiag(vb) > 0, vecdiag(vb), 0)); 
 	 end; 
	 else do;
       se_est2 = sqrt(vecdiag(vb0));
	 end;
/**/
	 estimat1 = beta;
	 _nozero1 = loc(se_est1>0);
	 _nozero2 = loc(se_est2>0);
     z1[_nozero1] = beta[_nozero1]/se_est1[_nozero1];    *z-statistics for beta;
     z2[_nozero2] = beta[_nozero2]/se_est2[_nozero2];    *z-statistics for beta;
/**/
	 p1 = choose(z1^=., 2 * (1 - cdf("t",abs(z1),N-nrow(beta))), .);
	 lower1 = beta - tinv(0.975, N-nrow(beta)) * se_est1;
	 upper1 = beta + tinv(0.975, N-nrow(beta)) * se_est1;
/**/
	 p2 = choose(z2^=., 2 * (1 - cdf("t",abs(z2),N-nrow(beta))), .);
	 lower2 = beta - tinv(0.975, N-nrow(beta)) * se_est2;
	 upper2 = beta + tinv(0.975, N-nrow(beta)) * se_est2;
/**/
	 p3 = choose(z1^=., 2 * (1 - cdf("normal",abs(z1))), .);
	 lower3 = beta - probit(0.975) * se_est1;
	 upper3 = beta + probit(0.975) * se_est1;
/**/
	 p4 = choose(z2^=., 2 * (1 - cdf("normal",abs(z2))), .);
	 lower4 = beta - probit(0.975) * se_est2;
	 upper4 = beta + probit(0.975) * se_est2;
   end;

     return (estimat1||se_est1||se_est2||p1||lower1||upper1||p2||lower2||upper2||p3||lower3||upper3||p4||lower4||upper4||conv1||conv2);

finish brgeegor;

store module=_all_;

quit;
