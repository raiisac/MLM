/* Assignment Concepts of multivariate, longitudinal and mixed models
Group 3*/

/* Create library */
libname MLMM "/home/u49934288/MLMM2";

/* 1. EXPLORATORY DATA ANALYSIS */
/* N, mean, stddev, min, max, median*/
proc means data=MLMM.mmse N mean stddev min max median maxdec=4;
run;

/* Distribution of msse*/
proc univariate data=MLMM.mmse;
   histogram;
run;

/* Differences between groups */
/* Differences for neuro */
proc sort data= mlmm.mmse out=mlmm.mmse_sortneuro;
by neuro;
run;
proc means data=mlmm.mmse_sortneuro min max mean median maxdec=4;
by neuro;
run;

/* Differences for housing */
proc sort data= mlmm.mmse out=mlmm.mmse_sorthousing;
by housing;
run;
proc means data=mlmm.mmse_sorthousing min max mean median maxdec=4;
by housing;
run;

/* Mean and variance structure */
title "Mean Response";
proc sgplot data=mlmm.mmse;
   vline time / response=mmse stat=mean limitstat=stderr datalabel;
   yaxis label='Mean +/- SEM';
run;
/* -> Intercept in random statement */

/* Mean and variance structure for neuro */
title "Mean Response by Neuro";
proc sgplot data=mlmm.mmse;
   vline time / response=mmse group=neuro stat=mean limitstat=stderr datalabel;
   yaxis label='Mean +/- SEM';
run;
/* -> Difference in intercept for neuro -> neuro in model statement */
/* -> No difference in slope for neuro -> NOT neuro*t in model statement */
/* -> Variance differences -> group = neuro in random statement */

/* Mean and variance structure for housing */
title "Mean Response by Housing";
proc sgplot data=mlmm.mmse;
   vline time / response=mmse group=housing stat=mean limitstat=stderr datalabel;
   yaxis label='Mean +/- SEM';
run;
/* -> Difference in intercept for housing -> housing in model statement */

/* Mean and variance structure for age*/
title "Mean Response";
proc sgplot data=mlmm.mmse;
   vline age / response=mmse  stat=mean limitstat=stderr datalabel;
   yaxis label='Mean +/- SEM';
run;
/* -> Difference in intercept for age -> age in model statement */
title "Mean Response";
proc sgplot data=mlmm.mmse;
   vline time / response=mmse group=age stat=mean limitstat=stderr datalabel;
   yaxis label='Mean +/- SEM';
run;
/* -> No difference in slope for age -> NOT age*t in model statement */

/* Plot of individual profiles */
/* Sort data by id and time */
proc sort data=mlmm.mmse out=mlmm.mmse_sortidtime;
	by id time;
run;
/* Transform table from narrow to wide */
data mlmm.mmse_wide;
	set mlmm.mmse_sortidtime;
	by id;
	retain id t_1 t_3 t_5 t_8 t_12 age neuro mmse time housing;
	keep id t_1 t_3 t_5 t_8 t_12 age neuro mmse time housing;
	if time=1 then
		t_1=mmse;
	else if time=3 then
		t_3=mmse;
	else if time=5 then
		t_5=mmse;
	else if time=8 then
		t_8=mmse;
	else if time=12 then
		t_12=mmse;
	if last.id=1 then
		output;
run;
/* Randomly sample 20 patients from the dataset */
proc surveyselect data=mlmm.mmse_wide out=mlmm.mmse_sample method=SRS 
		sampsize=20 seed=1234567;
run;
/* Transform table from wide to narrow */
data mlmm.mmse_samplenarrow;
	set mlmm.mmse_sample;
	by id;
	retain id time age neuro mmse time housing;
	keep id time age neuro mmse time housing;
	time=1;
	mmse=t_1;
	output;
	time=3;
	mmse=t_3;
	output;
	time=5;
	mmse=t_5;
	output;
	time=8;
	mmse=t_8;
	output;
	time=12;
	mmse=t_12;
	output;
run;
/* Plot the individual profiles */
proc sgplot data=mlmm.mmse_samplenarrow noautolegend;
	series x=time y=mmse / group=id lineattrs=(pattern=solid);
	xaxis values=(1, 3, 5, 8, 12);
run;

/* Correlation within the dataset */
proc corr data=mlmm.mmse;
run;
proc corr data=mlmm.mmse_wide;
run;


/* Explore missingness */
proc mi data=mlmm.mmse_wide nimpute=0; 
ods select misspattern; 
run;
proc mi data=mlmm.mmse_wide seed=4000 simple nimpute=0;   
title 'standard EM';
em itprint outem=mlmm.mmse_wide_missingness;
var NEURO mmse id housing age t_1 t_3 t_5 t_8 t_12; 
run;


/* 2. LINEAR MIXED MODEL WITH RANDOM INTERCEPT AND SLOPE */
/* Log transform time + */
data MLMM.mmse1;
	set MLMM.mmse;
	t=log(time);
	age=age-65;
run;

/* Parametrization 1: estimates for both groups*/
proc mixed data=MLMM.mmse1;
	class id neuro(ref="0");
	model mmse = neuro neuro*t / noint solution;
	random intercept t /type=un subject=id group=neuro g gcorr solution;
	ods exclude solutionr;
	ods output solutionr=MLMM.mmse1_lin;
	contrast 'neuro contrast 1 -1' neuro 1 -1;
	contrast 'neuro*t contrast 1 -1' neuro*t 1 -1;
run;

/* Parametrization 2: one group as baseline and the difference between the groups is given */
proc mixed data=MLMM.mmse1;
	class id neuro(ref="0");
	model mmse = t neuro neuro*t / solution;
	random intercept t /type=un subject=id group=neuro g gcorr solution;
	ods exclude solutionr;
	ods output solutionr=MLMM.mmse1_lin1;
run;
/* Estimated G matrix is not positive definite -> delete t or neuro from random statement */

/* (2.1) Random intercept + slope */
proc mixed data=MLMM.mmse1;
	class id neuro(ref="0");
	model mmse = t neuro neuro*t/ solution;
	random intercept t /type=un subject=id g gcorr solution;
	ods exclude solutionr;
	ods output solutionr=MLMM.mmse1_lin2;
run;
/* Reduced model */
proc mixed data=MLMM.mmse1;
	class id neuro(ref="0");
	model mmse = t neuro / solution;
	random intercept t /type=un subject=id g gcorr solution;
	ods exclude solutionr;
	ods output solutionr=MLMM.mmse1_lin2R;
run;

/* (2.2) Random intercept + group effect of neuro */
proc mixed data=MLMM.mmse1;
	class id neuro(ref="0");
	model mmse = t neuro neuro*t/ solution;
	random intercept /type=un subject=id group=neuro g gcorr solution;
	ods exclude solutionr;
	ods output solutionr=MLMM.mmse1_lin3;
run;
/* Reduced model */
proc mixed data=MLMM.mmse1;
	class id neuro(ref="0");
	model mmse = t neuro / solution;
	random intercept /type=un subject=id group=neuro g gcorr solution;
run;

/* t + neuro = significant */

/* 3. SCATTERPLOT */
/* Full model random intercept + slope */
proc transpose data=MLMM.mmse1_lin2 out=MLMM.mmse1_lin_trans (drop=_NAME_ 
		rename=(COL1=Intercept COL2=Slope));
	var Estimate;
	by id;
run;
proc sgscatter DATA=MLMM.mmse1_lin_trans;
	title "Random Effects Estimates";
	plot Intercept*Slope ;
run;
/* Negative correlation between the intercept and the slope */


/* 4. ADD BASELINE CHARACTERISTICS */
/* Parametrization 1: estimates for both groups*/
proc mixed data=MLMM.mmse1;
	class id neuro(ref="0") housing(ref="3");
	model mmse = neuro neuro*t age age*t housing housing*t/ noint solution;
	random intercept t/type=un subject=id g group=neuro gcorr solution;
	ods exclude solutionr;
run;

/* Parametrization 2: one group as baseline and the difference between the groups is given */
proc mixed data=MLMM.mmse1;
	class id neuro(ref="0") housing(ref="3");
	model mmse = t neuro neuro*t age age*t housing housing*t / solution;
	random intercept t/type=un subject=id group=neuro g gcorr solution;
	ods exclude solutionr;
run;

/* (4.1) Random intercept + slope */
proc mixed data=MLMM.mmse1 covtest;
	class id neuro(ref="0") housing(ref="3");
	model mmse = t neuro neuro*t age age*t housing housing*t / solution;
	random intercept t/type=un subject=id g gcorr solution;
	ods exclude solutionr;
	ods output solutionr=MLMM.mmse1_cov1;
run;

/* Reduced model */
proc mixed data=MLMM.mmse1;
	class id neuro(ref="0") housing(ref="3");
	model mmse = t neuro age housing / solution;
	random intercept t/type=un subject=id g gcorr solution;
	ods exclude solutionr;
run;

/* (4.2) Random intercept + group effect of neuro */
proc mixed data=MLMM.mmse1;
	class id neuro(ref="0") housing(ref="3");
	model mmse = t neuro neuro*t age age*t housing housing*t / solution;
	random intercept /type=un subject=id group=neuro g gcorr solution;
	ods exclude solutionr;
run;

/* Reduced model */
proc mixed data=MLMM.mmse1;
	class id neuro(ref="0") housing(ref="3");
	model mmse = t neuro age housing / solution;
	random intercept /type=un subject=id group=neuro g gcorr solution;
	ods exclude solutionr;
run;

/* t + neuro + age + housing significant */


/* 5. DICHOTOMIZATION */
data MLMM.mmse2;
	set MLMM.mmse1;
	if mmse>23 then
		abnormal=0;
	else
		abnormal=1;
run;

/* 6. LOGISTIC RANDOM INTERCEPTS MODEL */
/* Only neuro*/
/* Parametrization 1: estimates for both groups*/
proc glimmix data=MLMM.mmse2;
	class id neuro(ref="0") housing(ref="3");
	model abnormal (event='1')= neuro neuro*t / noint dist=binary 
		solution;
	random intercept /type=un subject=id group=neuro;
run;
/* Parametrization 2: one group as baseline and the difference between the groups is given */
proc glimmix data=MLMM.mmse2;
	class id neuro(ref="0") housing(ref="3");
	model abnormal (event='1')= t neuro neuro*t / dist=binary 
		solution;
	random intercept t /type=un subject=id group=neuro;
run;

/* (6.1) Random intercept + slope */
proc glimmix data=MLMM.mmse2;
	class id neuro(ref="0") housing(ref="3");
	model abnormal (event='1')= t neuro neuro*t / dist=binary 
		solution;
	random intercept t /type=un subject=id;
run;
/* "t neuro neuro*t" -> no convergence */
/* "t neuro" -> Estimated G matrix is not positive definite. -> remove the random slope */

/* (6.2) Random intercept + group effect of neuro */
proc glimmix data=MLMM.mmse2;
	class id neuro(ref="0") housing(ref="3");
	model abnormal (event='1')= t neuro neuro*t / dist=binary 
		solution;
	random intercept /type=un subject=id group=neuro;
run;

/* Reduced model */
proc glimmix data=MLMM.mmse2;
	class id neuro(ref="0") housing(ref="3");
	model abnormal (event='1')= neuro / dist=binary 
		solution;
	random intercept /type=un subject=id group=neuro;
run;
/* Neuro significant */


/* Additional covariates*/
/* Parametrization 1: estimates for both groups*/
proc glimmix data=MLMM.mmse2;
	class id neuro(ref="0") housing(ref="3");
	model abnormal (event='1')= neuro neuro*t age age*t housing housing*t / noint dist=binary 
		solution;
	random intercept /type=un subject=id group=neuro;
run;

/* Parametrization 2: one group as baseline and the difference between the groups is given */
proc glimmix data=MLMM.mmse2;
	class id neuro(ref="0") housing(ref="3");
	model abnormal (event='1')= t neuro neuro*t age age*t housing housing*t / dist=binary 
		solution;
	random intercept t /type=un subject=id group=neuro;
run;

/* (6.1) Random intercept + slope */
proc glimmix data=MLMM.mmse2;
	class id neuro(ref="0") housing(ref="3");
	model abnormal (event='1')= t neuro neuro*t age housing / dist=binary 
		solution;
	random intercept t /type=un subject=id;
run;
/* "t neuro neuro*t age age*t housing housing*t" -> no convergence */
/* "t neuro age housing" -> Estimated G matrix is not positive definite. -> remove the random slope */

/* (6.2) Random intercept + group effect of neuro */
proc glimmix data=MLMM.mmse2;
	class id neuro(ref="0") housing(ref="3");
	model abnormal (event='1')= t neuro neuro*t age age*t housing housing*t / dist=binary 
		solution;
	random intercept /type=un subject=id group=neuro;
run;

/* Reduced model */
proc glimmix data=MLMM.mmse2;
	class id neuro(ref="0") housing(ref="3");
	model abnormal (event='1')= neuro age housing / dist=binary 
		solution;
	random intercept /type=un subject=id group=neuro;
run;
/* Neuro + age + housing significant */


