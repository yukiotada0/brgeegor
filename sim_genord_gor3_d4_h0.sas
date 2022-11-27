/*Program Name : sim_genord_gor3_d4.sas*/
/*Programmer : Yukio Tada                                 */

%include "/***/brgeegor.sas";
%include "/***/gen_ordinal_gor3.sas";

%macro sim_genord_gor(N, LGOR, covtype);
	proc iml worksize=5000;
		load module=_all_;
		
		K=3;
		q=K-1;
		beta0={-0.2 0.7 0 -0.9 -0.6 -0.3 0 0 0}`;
		beta_init={0 0 0 0 0 0 0}`;
		design0={0 1 0 0 0 0 0, 0 0 1 0 0 0 0, 0 0 0 1 0 0 0, 0 0 0 0 0 0 0};
		T=nrow(design0);
		design1={1 1 0 0 1 0 0, 1 0 1 0 0 1 0, 1 0 0 1 0 0 1, 1 0 0 0 0 0 0};
		trans=J(K, K, 1);

		/*cumulative transform matrix*/
		do r=1 to K;

			do c=1 to K;

				if c>r then
					do;
						trans[r, c]=0;
					end;
			end;
		end;
		transinv=inv(trans);

		do tt=1 to T;
			Z0_t=I(q)||(J(q, 1, 1)@design0[tt, ]);
			Z1_t=I(q)||(J(q, 1, 1)@design1[tt, ]);
			Gam0_t=logistic(Z0_t*beta0);
			Gam0_t=Gam0_t//J(1, 1, 1);
			Gam1_t=logistic(Z1_t*beta0);
			Gam1_t=Gam1_t//J(1, 1, 1);
			pr0_t=transinv*Gam0_t;
			pr1_t=transinv*Gam1_t;
			pr0=pr0||pr0_t;
			pr1=pr1||pr1_t;
		end;
				
		N=&N.;
		SimNo=5000;
		ALPHA0=&LGOR./2;
        covtype = "&covtype.";
/* 		 */
        CALL JOINT_PI(pr0,ALPHA0,PI0,CODE_PI0,covtype,RC0);
        CALL JOINT_PI(pr1,ALPHA0,PI1,CODE_PI1,covtype,RC1);
        IF RC0=0 & RC1=0 THEN DO;
        CUSUM_PI0=CUSUM(PI0);
        CUSUM_PI1=CUSUM(PI1);
		CALL RANDSEED(12345);
        RanX=J(N,SimNo,.);
        CALL RANDGEN(RanX,"UNIFORM");
        INDATA=J(N#T,SimNo,.);
        DO J=1 TO SimNo;
          DO I=1 TO N;
            TEMP0=LOC(CUSUM_PI0>RanX[I,J]);
            LOCATION0=TEMP0[1];
            TEMP1=LOC(CUSUM_PI1>RanX[I,J]);
            LOCATION1=TEMP1[1];
            if I<=N/2 then INDATA[(I-1)#T+1:I#T,J]=CODE_PI0[LOCATION0,]`;
            else INDATA[(I-1)#T+1:I#T,J]=CODE_PI1[LOCATION1,]`;
          END;
        END;
        END;
/* 		 */
		ID=(1:N)`@J(T, 1, 1);
		TP=J(N, 1, 1)@(1:T)`;
		X=(J(N/2, 1, 1)@design0)//(J(N/2, 1, 1)@design1);

		do assoc0=1 to 2;

			do mgee0=1 to 2;
				free est0 bias0 stderr0 sqerr0 cover_z0 length_z0 test_z0 cover_t0 length_t0 test_t0 conv0;

				do s=1 to SimNo;
					free est00;
					est00=brgeegor(ID, TP, X, INDATA[, s], beta_init, 1e-3, 20, assoc0, mgee0);

					if est00[, 1] > . then
						do;
							est0=est0//est00[, 1]`;
							bias0=bias0//(est00[, 1]-beta0)`;
							stderr0=stderr0//est00[, 2]`;
							sqerr0=sqerr0//((est00[, 1]-beta0)##2)`;
							
							lower_ci_z=est00[, 1]-probit(0.975)*est00[, 2];
							upper_ci_z=est00[, 1]+probit(0.975)*est00[, 2];
							cover_z0=cover_z0//(lower_ci_z<beta0 & beta0<upper_ci_z)`;
							length_z0=length_z0//(upper_ci_z-lower_ci_z)`;
							test_z0=test_z0//(abs(est00[, 3])>probit(0.975))`;
							
							lower_ci_t=est00[,1]-tinv(0.975, N-nrow(beta0))*est00[,2];
							upper_ci_t=est00[,1]+tinv(0.975, N-nrow(beta0))*est00[,2];
							cover_t0=cover_t0//(lower_ci_t<beta0 & beta0<upper_ci_t)`;
							length_t0=length_t0//(upper_ci_t-lower_ci_t)`;
							test_t0=test_t0//(abs(est00[, 3])>tinv(0.975, N-nrow(beta0)))`;

							conv0=conv0//j(nrow(beta0), 1, 1)`;
						end;
					else
						do;
							est0=est0//j(nrow(beta0), 1, .)`;
							bias0=bias0//j(nrow(beta0), 1, .)`;
							stderr0=stderr0//j(nrow(beta0), 1, .)`;
							sqerr0=sqerr0//j(nrow(beta0), 1, .)`;
							
							cover_z0=cover_z0//j(nrow(beta0), 1, .)`;
							length_z0=length_z0//j(nrow(beta0), 1, .)`;
							test_z0=test_z0//j(nrow(beta0), 1, .)`;
							
							cover_t0=cover_t0//j(nrow(beta0), 1, .)`;
							length_t0=length_t0//j(nrow(beta0), 1, .)`;
							test_t0=test_t0//j(nrow(beta0), 1, .)`;
							
							conv0=conv0//j(nrow(beta0), 1, 0)`;
						end;
				end;
				param=param||(1:nrow(beta0));
				size=size||J(1, nrow(beta0), N);
				lgor=lgor||J(1, nrow(beta0), ALPHA0);
				wcov=wcov||J(1, nrow(beta0), assoc0);
				mgee=mgee||J(1, nrow(beta0), mgee0);
				est=est||est0;
				bias=bias||bias0;
				stderr=stderr||stderr0;
				sqerr=sqerr||sqerr0;
				cover_z=cover_z||cover_z0;
				length_z=length_z||length_z0;
				test_z=test_z||test_z0;
				cover_t=cover_t||cover_t0;
				length_t=length_t||length_t0;
				test_t=test_t||test_t0;
				conv=conv||conv0;
			end;
		end;
		_result1_1=(size//lgor//wcov//mgee//param//est);
		_result1_2=(size//lgor//wcov//mgee//param//stderr);

		/*  print _result1;*/
		create result1_1 from _result1_1;
		append from _result1_1;
		close result1_1;
		create result1_2 from _result1_2;
		append from _result1_2;
		close result1_2;
		mean_est=mean(est);
		mean_bias=mean(bias);
		mean_se=mean(stderr);
		mean_sqerr=mean(sqerr);
		simse=std(est);
		se_simse=mean_se/simse;
		mean_cover_z=mean(cover_z);
		mean_length_z=mean(length_z);
		mean_test_z=mean(test_z);
		mean_cover_t=mean(cover_t);
		mean_length_t=mean(length_t);
		mean_test_t=mean(test_t);
		n_conv=conv[+, ];
		_result2_1=(size//lgor//wcov//mgee//param//mean_est//mean_bias//mean_se//mean_sqerr//simse//se_simse//mean_cover_z//mean_length_z//mean_test_z//mean_cover_t//mean_length_t//mean_test_t//n_conv)`;
		varNames2={"N" "LGOR" "WCOV" "mGEE" "param" "est" "bias" "se" "mserr" "simse"
			"se_simse" "coverage_z" "length_z" "power_z" "coverage_t" "length_t" "power_t" "iterations"};
		print _result2_1[colname=varNames2];
		create result2_1 from _result2_1[colname=varNames2];
		append from _result2_1;
		close result2_1;

		/**/
		miss=missing(est);
		rownomiss=loc(miss[, +]=0);
		mean_est=mean(est[rownomiss, ]);
		mean_bias=mean(bias[rownomiss, ]);
		mean_se=mean(stderr[rownomiss, ]);
		mean_sqerr=mean(sqerr[rownomiss, ]);
		simse=std(est[rownomiss, ]);
		se_simse=mean_se/simse;
		mean_cover_z=mean(cover_z[rownomiss, ]);
		mean_length_z=mean(length_z[rownomiss, ]);
		mean_test_z=mean(test_z[rownomiss, ]);
		mean_cover_t=mean(cover_t[rownomiss, ]);
		mean_length_t=mean(length_t[rownomiss, ]);
		mean_test_t=mean(test_t[rownomiss, ]);
		n_conv=conv[rownomiss, ][+, ];
		_result2_2=(size//lgor//wcov//mgee//param//mean_est//mean_bias//mean_se//mean_sqerr//simse//se_simse//mean_cover_z//mean_length_z//mean_test_z//mean_cover_t//mean_length_t//mean_test_t//n_conv)`;
		print _result2_2[colname=varNames2];
		create result2_2 from _result2_2[colname=varNames2];
		append from _result2_2;
		close result2_2;
		
/* 		create indata2 from indata; */
/* 		append from indata; */
/* 		close indata2; */
/* 		 */
/* 		print ranx, indata, pr0, PI0, CUSUM_PI0; */
		quit;
		libname out "P:\desktop6\biostat\SPH\maunscript\output";

	data out.out&N._&LGOR._&covtype._k3d4_est_aSE_h0;
		set result1_1;
	run;

	data out.out&N._&LGOR._&covtype._k3d4_se_aSE_h0;
		set result1_2;
	run;

	data out.out&N._&LGOR._&covtype._k3d4_sum1_aSE_h0;
		set result2_1;
	run;

	data out.out&N._&LGOR._&covtype._k3d4_sum2_aSE_h0;
		set result2_2;
	run;

%mend sim_genord_gor;

/*Exchangeable*/
%sim_genord_gor(20, 2, exch);
%sim_genord_gor(20, 3, exch);
%sim_genord_gor(20, 4, exch);

/**/
%sim_genord_gor(30, 2, exch);
%sim_genord_gor(30, 3, exch);
%sim_genord_gor(30, 4, exch);

/**/
%sim_genord_gor(40, 2, exch);
%sim_genord_gor(40, 3, exch);
%sim_genord_gor(40, 4, exch);

/**/
%sim_genord_gor(50, 2, exch);
%sim_genord_gor(50, 3, exch);
%sim_genord_gor(50, 4, exch);

/*AR-type*/
%sim_genord_gor(20, 2, ar);
%sim_genord_gor(20, 3, ar);
%sim_genord_gor(20, 4, ar);

/**/
%sim_genord_gor(30, 2, ar);
%sim_genord_gor(30, 3, ar);
%sim_genord_gor(30, 4, ar);

/**/
%sim_genord_gor(40, 2, ar);
%sim_genord_gor(40, 3, ar);
%sim_genord_gor(40, 4, ar);

/**/
%sim_genord_gor(50, 2, ar);
%sim_genord_gor(50, 3, ar);
%sim_genord_gor(50, 4, ar);
