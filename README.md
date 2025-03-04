This file describes the structure of the files in the zip folder. This contains the code to generate the figures from the paper and to sustain the conjectures.

-- The code is written in Matlab based on the Performance Estimation Toolbox (PESTO).
-- Required packages:
	1. Performance Estimation Toolbox (PESTO), downloaded from: https://github.com/PerformanceEstimation/Performance-Estimation-Toolbox
		-- on this link are available the installation instructions, with the required packages:
	2. Yalmip
	3. An SDP solver (we use Mosek)
-- The PESTO toolbox is updated as follows:
	-> File Hypoconvex.m in folder Functions_clsses
		-- this includes interpolation inequalities for weakly convex (or hypoconvex) functions
	-> File pep.m includes the class of Hypoconvex functions
	The easiest way to use the code: add the folder AIStats_Code to the path and then run the scripts.			
	
Details on the scripts:
--> main folder:
	--> PEP_DCA_no_opt.m: a function to solve the (PEP-DCA) problem to find numerically the worst-case value for a given setup. It doesn't incorporate conditions including a stationary point F star.
	
	--> PEP_DCA_with_opt.m: similar to PEP_DCA_no_opt.m, but incorporating conditions to include a stationary point F star.
	
	--> PEP_simulation_N_steps.m: a script to solve the (PEP-DCA) problem for every iteration k up to some maximum number N. In the default arguments it includes different setups covering all possible cases (all 6 regimes). These setups were used to generate Figures 3,4,5,6 from Section D to check numerically the tightness for each of the 6 regimes.
	
	--> compute_DCA_rate_6_regimes.m: a function which determines the expression of the denominator pi after one iteration, covering all possible 6 regimes. 


--> Describe_regimes:
	--> All_regimes_no_contours.m: it generates a similar plot with Figure 1, but without the contour lines. It can be adapted to other curvature choices as well. 
	
	--> All_regimes_with_contours.m: this is the script to generate Figure 1.
	
	--> Illustrate_wc_func_p1.m: this script builds a worst-case function example corresponding to regime p1 when both functions are smooth. This is described in the Appendix in Section D, following Proposition D.1. It generates Figure 7, which also provides some intuition on how the DCA iterations look like.

	--> Plot_regimes_f1_nonsmooth.m: a sketch of the possible regimes when L1=Inf. For scaling reasons, we use L1=10 x L2.
	--> Plot_regimes_f2_nonsmooth.m: same, but having L2=Inf and L1<Inf.

--> Best_splitting_problem: it includes functions to compute the optimal curvature subtraction in both functions to improve the convergence rate (larger denominator pi)
	--> get_p_lambda.m: it computes the corresponding denominator after subtracting both curvatures by lambda
	--> get_opt_lambda.m: it computes the lambda corresponding to the largest improvement
	--> plot_denominators_as_func_of_shifting_param_lambda.m: this is a script to play with the curvature parameters and plots the denominators as function of lambda. In particular, this script was used to generate Table 3. Figures associated with the different setups are in the folder Figures/Best_splitting_lambdas.
	
	--> best_splittings_for_smooth_nonconvex_objective_F.m: it is used to generate the plot in the Appendix in Section F, extending Section 4. This shows all possible regimes when fixing curvature muF and LF, and illustrates how to shift the curvature for several setups in terms of L2 and mu2. 
	
	--> compare_GD_vs_smooth_DCA.m: This script was used to generate Example 4.1, which shows that the best worst-case for smooth DCA is better than the best worst-case when applying Gradient Descent (GD) with optimal stepsize, for smooth nonconvex functions.

--> SPCA_numerical_experiments: this includes the implementation of the numerical example described in Section 7
	--> experiment_SPCA_Case_1.m, experiment_SPCA_Case_2.m: the script to run the experiments.
	--> rng_1.mat: the seed used to generate the random initializations x0.
	--> Sigma_vals.mat: the randomly generated covariance matrix Sigma
	
--> Figures: This folder includes the plots generated with the provided code.
	--> Best_splitting_lambdas: see description of best_splittings_for_smooth_nonconvex_objective_F.m
	--> PEP_worst_case_values_simulations: see description of PEP_simulation_N_steps.m
	--> Regimes_nonsmooth_case: plots from Appendix D.
	--> Worst_case_functions_p1_p2: plots from Appendices B.1 and B.2 on worst-case examples.
	--> SPCA: figures with the SPCA simulations.