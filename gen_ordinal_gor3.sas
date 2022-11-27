/*SAS macro for generating correlated discrete ordinal data
SEE: Ibrahim NA et al.(2011)*/
/* modification by Yukio Tada on 11Feb2021 */

/******************************************/
/*Ver2: To assume AR(1) as covariance structure*/
/******************************************/

/* NOTE: THE linprog MODUL IS TAKEN FROM SAS HELP*/
PROC IML;
OPTIONS NODATE NOCENTER PS=1000;

START PI_GOR(PI1,PI2,ALPHA,RESULT);
PSI=EXP(ALPHA);

P1 = cusum(PI1);
P2 = cusum(PI2);

Nc=NROW(P1);
S_ijk=J(Nc,Nc,0);
F_ijk=J(Nc,Nc,0);
PP_ijk=J(Nc,Nc,0);
DO L=1 to Nc;
DO M=1 to Nc;
S_ijk[L,M]=SQRT((1+(P1[L,1]+P2[M,1])*(PSI-1))**2+4*PSI*(1-PSI)*P1[L,1]*P2[M,1]);
IF PSI^=1 THEN F_ijk[L,M]=(1+(P1[L,1]+P2[M,1])*(PSI-1)-S_ijk[L,M])/(2*(PSI-1));
ELSE F_ijk[L,M]=P1[L,1]*P2[M,1];
     if L=1 & M=1 then PP_ijk[L,M]=F_ijk[L,M];
else if L=1 & M>1 then PP_ijk[L,M]=F_ijk[L,M]-F_ijk[L,M-1];
else if L>1 & M=1 then PP_ijk[L,M]=F_ijk[L,M]-F_ijk[L-1,M];
else if L>1 & M>1 then PP_ijk[L,M]=F_ijk[L,M]-F_ijk[L,M-1]-F_ijk[L-1,M]+F_ijk[L-1,M-1];
END;
END;
RESULT=PP_ijk;
FINISH PI_GOR;

/* Subroutine to solve Linear Programs */
/* names: names of the decision variables */
/* obj: coefficients of the objective function */
/* maxormin: the value �fMAX�f or �fMIN�f, upper or lowercase */
/* coef: coefficients of the constraints */
/* rel: character array of values: �f<=�f or �f>=�f or �f=�f */
/* rhs: right-hand side of constraints */
/* activity: returns the optimal value of decision variables*/
/* SOURCE: SAS HELP */
start linprog(names, obj, maxormin, coef, rel, rhs, activity,rc);
bound=1.0e10;
m=nrow(coef);
n=ncol(coef);
/* Convert to maximization */
if upcase(maxormin)="MIN" then o=-1;
else o=1;
/* Build logical variables */
rev=(rhs<0);
adj=(-1*rev)+^rev;
ge=((rel = ">==") & ^rev) | ((rel = "<=") & rev);
eq=(rel="=");
if max(ge)=1 then do;
sr=I(m);
logicals=-sr[,loc(ge)]||I(m);
artobj=repeat(0,1,ncol(logicals)-m)|(eq+ge)`;
end;
else do;
logicals=I(m);
artobj=eq`;
end;
nl=ncol(logicals);
nv=n+nl+2;
/* Build coef matrix */
a=((o*obj)||repeat(0,1,nl)||{
-1 0
})//
(repeat(0,1,n)||-artobj||{
0 -1
})//
((adj#coef)||logicals||repeat(0,m,2));
/* rhs, lower bounds, and basis */
b={0,0}//(adj#rhs);
L=repeat(0,1,nv-2)||-bound||-bound;
basis=nv-(0:nv-1);
/* Phase 1 - primal feasibility */
call lp(rc,x,y,a,b,nv,,l,basis);
print ({
" ",
"**********Primal infeasible problem************",
" ",
"*********Numerically unstable problem**********",
"*********Singular basis encountered************",
"*******Solution is numerically unstable********",
"***Subroutine could not obtain enough memory***",
"**********Number of iterations exceeded********"
}[rc+1]);
if x[nv] ^=0 then do;
print "**********Primal infeasible problem************";
stop;
end;
if rc>0 then stop;
/* phase 2 - dual feasibility */
u=repeat(.,1,nv-2)||{. 0};
L=repeat(0,1,nv-2)||-bound||0;
call lp(rc,x,y,a,b,nv-1,u,l,basis);
/* Report the solution */
print ({
"*************Solution is optimal***************",
"*********Numerically unstable problem**********",
"**************Unbounded problem****************",
"*******Solution is numerically unstable********",
"*********Singular basis encountered************",
"*******Solution is numerically unstable********",
"***Subroutine could not obtain enough memory***",
"**********Number of iterations exceeded********"
}[rc+1]);
value=o*x [nv-1];
activity= x [1:n];
/* print, "Objective Value" value;
print,"Decision Variables" activity;
lhs=coef*x[1:n];
dual=y[3:m+2];
print,"Constraints" lhs rel rhs dual,
"**********************************************";
*/
finish linprog;

START JOINT_PI(MARGINAL,ALPHA,RESULT,CODE_PI,covtype,rc);
/* MODUL COMPUTING JOINT PROBABILTY GIVEN 
MARGINAL MODEL AND SPECIFIC CORRELATION MATRIX.
CORRELATION IN THE FORM OF GLOBAL ODDS RATIO.
INPUT:
MARGINAL: MATRIX OF MARGINAL DISTRIBUTION FOR ALL LEVEL AND ALL TIME. 
DIMMENSION = L X T; 
ONE COLUMN REPRESENT MARGINAL MODEL FOR RESPECTIVE TIME.
CORRELATION: IS [T(T-1)/2]x[T(T-1)/2] MATRIX, REPRESENT ALL COMBINATION
OF THE ASSOCIATION (GROBAL ODDS RATIO - ALPHA=log(GOR)).
THE RESULT IS A MATRIX RESULT.
covtype="EXCH": Exchangeable, covtype="AR": Autoregressive*/
N=(MARGINAL^=.);
N=N[+,];
T=NCOL(MARGINAL);
/* T = N TIME*/
COUNT=1;
COMBI_T=J(T*(T-1)/2,2,.);
DO I=1 TO T-1;
DO J=I+1 TO T;
COMBI_T[COUNT,1]=I;
COMBI_T[COUNT,2]=J;
COUNT=COUNT+1;
END;
END;
TEMP=1:N[1];
TEMP=TEMP`;
CODE_PI=TEMP;
DO I =2 TO T;
TEMP=1:N[I];
TEMP=TEMP`;
T1=TEMP@J(NROW(CODE_PI),1,1);
CODE_PI=J(N[I],1,1)@CODE_PI;
CODE_PI=T1||CODE_PI;
END;
CODE_PI=NCOL(CODE_PI):1//CODE_PI;
CODE_PI=CODE_PI`;
CALL SORT(CODE_PI,{1},{1});
CODE_PI=CODE_PI`;
CODE_PI=CODE_PI[2:NROW(CODE_PI),];
/* NOTE:
EQUATION AX=B
CODE_PI IS LEVEL CODE FOR X
BELOW IS COMMAND TO CONSTRUCT VECTOR B
*/
CODE_TEMP=1:N[1]-1;
CODE_TEMP=CODE_TEMP`;
CODE_TEMP2=J(N[1]-1,T,0);
DO I=1 TO T;
IF I=1 THEN DO;
CODE_B=CODE_TEMP2;
CODE_B[,1]=CODE_TEMP;
B=MARGINAL[1:N[1]-1,I];
END;
ELSE DO;
TEMP=J(N[I]-1,T,0);
CODE_TEMP=1:N[I]-1;
CODE_TEMP=CODE_TEMP`;
TEMP[,I]=CODE_TEMP;
CODE_B=CODE_B//TEMP;
B=B//MARGINAL[1:N[I]-1,I];
END;
END;
DO I=1 TO NROW(COMBI_T);
T1=COMBI_T[I,1];
T2=COMBI_T[I,2];
P1=MARGINAL[1:N[T1],T1];
P2=MARGINAL[1:N[T2],T2];
/**/
if upcase(covtype)="EXCH" then ALPHA_T=ALPHA;
if upcase(covtype)="AR" then ALPHA_T=ALPHA/(abs(T1-T2));
/**/
CALL PI_GOR(P1,P2,ALPHA_T,MATRIX);
RESULT=MATRIX[1:N[T1]-1, 1:N[T2]-1];
DO J=1 TO N[T2]-1;
TEMP=J(N[T1]-1,T,0);
TEMP[,T1]=(1:N[T1]-1)`;
TEMP[,T2]=J(NROW(CODE_TEMP),1,J);
CODE_B=CODE_B//TEMP;
B=B//RESULT[,J];
END;
END;
/* CONSTRUCTING MATIX A IN AX = B*/
A=J(1,NROW(CODE_PI),1);
DO I=1 TO NROW(B);
TEMP=CODE_B[I,];
T1=LOC(TEMP^=0);
IF NCOL(T1)=1 THEN A=A//(CODE_PI[,T1]=TEMP[,T1])`;
ELSE A=A//((CODE_PI[,T1[1]]=TEMP[,T1[1]])&(CODE_PI[,T1[2]]=TEMP[,T1[2]]))`;
END;
B=1//B;
NAMES=1:NCOL(A);
OBJ=J(1,NCOL(A),1);
RELATION=J(NROW(B),1,"=");
CALL LINPROG(NAMES,OBJ,"MAX",A,RELATION,B,RESULT,RC);
FINISH JOINT_PI;

store module=_all_;

QUIT;
