/*Program Name : sim_genord_gor3_d4_h0.sas*/
/*Programmer : Yukio Tada                                 */

%include "***\brgeegor.sas";
%include "***\gen_ordinal_gor3.sas";

%macro sim_genord_gor_h0(N, LGOR, covtype);
	proc iml worksize=5000;
		load module=_all_;
		
		K=3;
		q=K-1;
		beta0={1.1 2.2 0 -1.8 -1.2 -0.6 0 0 0}`;
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
				free est0 bias0 stderr10 stderr20 sqerr0 
                     cover_t10 length_t10 test_t10
                     cover_t20 length_t20 test_t20
                     cover_z10 length_z10 test_z10
                     cover_z20 length_z20 test_z20
                     conv10 conv20 conv30 conv40 separate_a0 separate_b0 cseparate0 flag10 flag20 flag30 flag40;

				do s=1 to SimNo;
					free est00;
					est00=brgeegor(ID, TP, X, INDATA[, s], beta_init, 1e-4, 50, assoc0, mgee0);
/**/
					separate_a00 = 0;
					separate_b00 = 0;
					cseparate00 = 0;
					do _t=1 to T;
					  if sum(INDATA[loc(TP=_t & X[,1]=0),s]=1)=N/2 | sum(INDATA[loc(TP=_t & X[,1]=0),s]=K)=N/2 |
                         sum(INDATA[loc(TP=_t & X[,1]=1),s]=1)=N/2 | sum(INDATA[loc(TP=_t & X[,1]=1),s]=K)=N/2 then separate_a00=1; * complete or quasi-complete separation;
					  if sum(INDATA[loc(TP=_t & X[,1]=0),s]=1)=N/2 | sum(INDATA[loc(TP=_t & X[,1]=0),s]=K)=N/2 |
                         sum(INDATA[loc(TP=_t & X[,1]=1),s]=1)=N/2 | sum(INDATA[loc(TP=_t & X[,1]=1),s]=K)=N/2 | 
                         (sum(INDATA[loc(TP=_t & X[,1]=0),s]=1)=0 & sum(INDATA[loc(TP=_t & X[,1]=1),s]=K)=0) | 
                         (sum(INDATA[loc(TP=_t & X[,1]=1),s]=1)=0 & sum(INDATA[loc(TP=_t & X[,1]=0),s]=K)=0) then separate_b00=1; * complete or quasi-complete separation;
					  if (sum(INDATA[loc(TP=_t & X[,1]=0),s]=1)=N/2 & sum(INDATA[loc(TP=_t & X[,1]=1),s]=K)=N/2) |
                         (sum(INDATA[loc(TP=_t & X[,1]=1),s]=1)=N/2 & sum(INDATA[loc(TP=_t & X[,1]=0),s]=K)=N/2) then cseparate00=1; * complete separation;
					end;
/**/
					if est00[, 1] > . then 
						do;
							est0=est0//est00[, 1]`;
							bias0=bias0//(est00[, 1]-beta0)`;
							stderr10=stderr10//est00[, 2]`;
							stderr20=stderr20//est00[, 3]`;
							sqerr0=sqerr0//((est00[, 1]-beta0)##2)`;
							
							cover_t10=cover_t10//(est00[,5]<=beta0 & beta0<=est00[,6])`;
							length_t10=length_t10//(est00[,6]-est00[,5])`;
							test_t10=test_t10//(est00[, 4] < 0.05)`;
							
							cover_t20=cover_t20//(est00[,8]<=beta0 & beta0<=est00[,9])`;
							length_t20=length_t20//(est00[,9]-est00[,8])`;
							test_t20=test_t20//(est00[, 7] < 0.05)`;

							cover_z10=cover_z10//(est00[,11]<=beta0 & beta0<=est00[,12])`;
							length_z10=length_z10//(est00[,12]-est00[,11])`;
							test_z10=test_z10//(est00[, 10] < 0.05)`;

							cover_z20=cover_z20//(est00[,14]<=beta0 & beta0<=est00[,15])`;
							length_z20=length_z20//(est00[,15]-est00[,14])`;
							test_z20=test_z20//(est00[, 13] < 0.05)`;

						end;
					else
						do;
							est0=est0//j(nrow(beta0), 1, .)`;
							bias0=bias0//j(nrow(beta0), 1, .)`;
							stderr10=stderr10//j(nrow(beta0), 1, .)`;
							stderr20=stderr20//j(nrow(beta0), 1, .)`;
							sqerr0=sqerr0//j(nrow(beta0), 1, .)`;
							
							cover_t10=cover_t10//j(nrow(beta0), 1, .)`;
							length_t10=length_t10//j(nrow(beta0), 1, .)`;
							test_t10=test_t10//j(nrow(beta0), 1, .)`;
							
							cover_t20=cover_t20//j(nrow(beta0), 1, .)`;
							length_t20=length_t20//j(nrow(beta0), 1, .)`;
							test_t20=test_t20//j(nrow(beta0), 1, .)`;
							
							cover_z10=cover_z10//j(nrow(beta0), 1, .)`;
							length_z10=length_z10//j(nrow(beta0), 1, .)`;
							test_z10=test_z10//j(nrow(beta0), 1, .)`;
							
							cover_z20=cover_z20//j(nrow(beta0), 1, .)`;
							length_z20=length_z20//j(nrow(beta0), 1, .)`;
							test_z20=test_z20//j(nrow(beta0), 1, .)`;
							
						end;
                      conv10=conv10//((choose(est00[,1]^=., 1, 0)` + (max(abs(est00[,1])) < 10))=2);
                      conv20=conv20//((choose(est00[,1]^=., 1, 0)` + (max(abs(est00[,1])) < 10))=2);
					  conv30=conv30//est00[, 16]`;
					  conv40=conv40//est00[, 17]`;
					  separate_a0=separate_a0//j(1,nrow(beta0),separate_a00);
					  separate_b0=separate_b0//j(1,nrow(beta0),separate_b00);
					  cseparate0=cseparate0//j(1,nrow(beta0),cseparate00);
					  flag10=flag10//((choose(est00[,1]^=., 1, 0)` + (max(abs(est00[,1])) < 10))<2 & separate_a00=0);
					  flag20=flag20//((choose(est00[,1]^=., 1, 0)` + (max(abs(est00[,1])) < 10))<2 & separate_a00=0);
					  flag30=flag30//((choose(est00[,16]=0, 1, 0)` + choose(est00[,1]^=., 1, 0)` + (max(abs(est00[,1])) < 10))=3);
					  flag40=flag40//((choose(est00[,17]=0, 1, 0)` + choose(est00[,1]^=., 1, 0)` + (max(abs(est00[,1])) < 10))=3);

				end;
				param=param||(1:nrow(beta0));
				size=size||J(1, nrow(beta0), N);
				lgor=lgor||J(1, nrow(beta0), ALPHA0);
				wcov=wcov||J(1, nrow(beta0), assoc0);
				mgee=mgee||J(1, nrow(beta0), mgee0);
				est=est||est0;
				bias=bias||bias0;
				stderr1=stderr1||stderr10;
				stderr2=stderr2||stderr20;
				sqerr=sqerr||sqerr0;
				cover_t1=cover_t1||cover_t10;
				length_t1=length_t1||length_t10;
				test_t1=test_t1||test_t10;
				cover_t2=cover_t2||cover_t20;
				length_t2=length_t2||length_t20;
				test_t2=test_t2||test_t20;
				cover_z1=cover_z1||cover_z10;
				length_z1=length_z1||length_z10;
				test_z1=test_z1||test_z10;
				cover_z2=cover_z2||cover_z20;
				length_z2=length_z2||length_z20;
				test_z2=test_z2||test_z20;
				conv1=conv1||conv10;
				conv2=conv2||conv20;
				conv3=conv3||conv30;
				conv4=conv4||conv40;
				separate_a=separate_a||separate_a0;
				separate_b=separate_b||separate_b0;
				cseparate=cseparate||cseparate0;
				flag1=flag1||flag10;
				flag2=flag2||flag20;
				flag3=flag3||flag30;
				flag4=flag4||flag40;

				size_=size_||J(1, 9, N);
				lgor_=lgor_||J(1, 9, ALPHA0);
				wcov_=wcov_||J(1, 9, assoc0);
				mgee_=mgee_||J(1, 9, mgee0);
				param_=param_||{-9 -8 -7 -6 -5 -4 -3 -2 -1};
				info=info||conv10[,1]||stderr10[,<>]||conv20[,1]||stderr20[,<>]||separate_a0[,1]||separate_b0[,1]||cseparate0[,1]||conv30[,1]||conv40[,1];
			end;
		end;
		_result1_1=(size//lgor//wcov//mgee//param//est);
		_result1_2=(size//lgor//wcov//mgee//param//stderr1);
		_result1_3=(size//lgor//wcov//mgee//param//stderr2);
		_result1_4=(size_//lgor_//wcov_//mgee_//param_//info);

		create result1_1 from _result1_1;
		append from _result1_1;
		close result1_1;
		create result1_2 from _result1_2;
		append from _result1_2;
		close result1_2;
		create result1_3 from _result1_3;
		append from _result1_3;
		close result1_3;
		varNames1={"conv1_1" "maxse1_1" "conv2_1" "maxse2_1" "separate_a1" "separate_b1" "cseparate1" "conv3_1" "conv4_1"
                   "conv1_2" "maxse1_2" "conv2_2" "maxse2_2" "separate_a2" "separate_b2" "cseparate2" "conv3_2" "conv4_2"
                   "conv1_3" "maxse1_3" "conv2_3" "maxse2_3" "separate_a3" "separate_b3" "cseparate3" "conv3_3" "conv4_3"
                   "conv1_4" "maxse1_4" "conv2_4" "maxse2_4" "separate_a4" "separate_b4" "cseparate4" "conv3_4" "conv4_4"};
		create result1_4 from _result1_4[colname=varNames1];
		append from _result1_4;
		close result1_4;
		mean_est=mean(est);
		mean_bias=mean(bias);
		mean_sqerr=mean(sqerr);
		median_se1=median(stderr1);
		median_se2=median(stderr2);
		simse=std(est);
		se1_simse=median_se1/simse;
		se2_simse=median_se2/simse;
		mean_cover_t1=mean(cover_t1);
		median_length_t1=median(length_t1);
		mean_test_t1=mean(test_t1);
		mean_cover_t2=mean(cover_t2);
		median_length_t2=median(length_t2);
		mean_test_t2=mean(test_t2);
		mean_cover_z1=median(cover_z1);
		median_length_z1=median(length_z1);
		mean_test_z1=mean(test_z1);
		mean_cover_z2=median(cover_z2);
		median_length_z2=median(length_z2);
		mean_test_z2=mean(test_z2);
		n_conv1=conv1[+, ];
		n_conv2=conv2[+, ];
		n_conv3=conv3[+, ];
		n_conv4=conv4[+, ];
		n_separate_a=separate_a[+, ];
		n_separate_b=separate_b[+, ];
		n_cseparate=cseparate[+, ];
		n_flag1=flag1[+, ];
		n_flag2=flag2[+, ];
		n_flag3=flag3[+, ];
		n_flag4=flag4[+, ];
		_result2_1=(size//lgor//wcov//mgee//param//mean_est//mean_bias//mean_sqerr//median_se1//median_se2//simse//se1_simse//se2_simse//
                    mean_cover_t1//median_length_t1//mean_test_t1//mean_cover_t2//median_length_t2//mean_test_t2//
                    mean_cover_z1//median_length_z1//mean_test_z1//mean_cover_z2//median_length_z2//mean_test_z2//
                    n_conv1//n_conv2//n_conv3//n_conv4//n_separate_a//n_separate_b//n_cseparate//n_flag1//n_flag2//n_flag3//n_flag4)`;
		varNames2={"N" "LGOR" "WCOV" "mGEE" "param" "est" "bias" "mserr" "se1" "se2" "SimSE" "se1/SimSE" "se2/SimSE"
                   "coverage_t1" "length_t1" "power_t1" "coverage_t2" "length_t2" "power_t2"
                   "coverage_z1" "length_z1" "power_z1" "coverage_z2" "length_z2" "power_z2"
                   "convergence1" "convergence2" "convergence3" "convergence4" "separate_a" "separate_b" "cseparate" "flag1" "flag2" "flag3" "flag4"};
		print _result2_1[colname=varNames2];
		create result2_1 from _result2_1[colname=varNames2];
		append from _result2_1;
		close result2_1;

		quit;
		libname out "P:\desktop6\biostat\SPH\maunscript\output\k3d4\5000\&covtype.\h0";

	data out.out&N._&LGOR._&covtype._k3d4_est_h0;
		merge result1_1 result1_4;
	run;

	data out.out&N._&LGOR._&covtype._k3d4_se1_h0;
		merge result1_2 result1_4;
	run;

	data out.out&N._&LGOR._&covtype._k3d4_se2_h0;
		merge result1_3 result1_4;
	run;

	data out.out&N._&LGOR._&covtype._k3d4_summary_h0;
		set result2_1;
	run;

	proc datasets lib = work;
	  delete result:;
	quit;

%mend sim_genord_gor_h0;

/*Exchangeable*/
%sim_genord_gor_h0(20, 2, exch);
%sim_genord_gor_h0(20, 3, exch);
%sim_genord_gor_h0(20, 4, exch);

/**/
%sim_genord_gor_h0(30, 2, exch);
%sim_genord_gor_h0(30, 3, exch);
%sim_genord_gor_h0(30, 4, exch);

/**/
%sim_genord_gor_h0(40, 2, exch);
%sim_genord_gor_h0(40, 3, exch);
%sim_genord_gor_h0(40, 4, exch);

/**/
%sim_genord_gor_h0(50, 2, exch);
%sim_genord_gor_h0(50, 3, exch);
%sim_genord_gor_h0(50, 4, exch);

/*AR-type*/
%sim_genord_gor_h0(20, 2, ar);
%sim_genord_gor_h0(20, 3, ar);
%sim_genord_gor_h0(20, 4, ar);

/**/
%sim_genord_gor_h0(30, 2, ar);
%sim_genord_gor_h0(30, 3, ar);
%sim_genord_gor_h0(30, 4, ar);

/**/
%sim_genord_gor_h0(40, 2, ar);
%sim_genord_gor_h0(40, 3, ar);
%sim_genord_gor_h0(40, 4, ar);

/**/
%sim_genord_gor_h0(50, 2, ar);
%sim_genord_gor_h0(50, 3, ar);
%sim_genord_gor_h0(50, 4, ar);
