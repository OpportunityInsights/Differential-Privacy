********************************************************************************
*
*	 A Practical Method to Reduce Privacy Loss
*		when Disclosing Statistics Based on Small Samples 
*
*    Example of how to implement the noise-infusion algorithm 
* 		outlined in Appendix A of Chetty and Friedman (2019).
*
*	 Updated: 1/10/2019
*
********************************************************************************
* 
* Objective: apply the noise-infusion algorithm in order to publicly release 
* the predicted value of child income rank at the 25th percentile of the parental
* income distribution in each of the cells of a simulated dataset of 10,000 individuals.

use "private_data_by_cells.dta", clear

*Choose distribution of omega (see Step 4 of algorithm) from laplace or normal
local omega_dist normal

* ------------------------------------------------------------------------------
* 1. Calculate Local Sensitivity
* ------------------------------------------------------------------------------
	
	*Generate variables to fill in with key parameters in each cell g:
	gen theta_g = .     //estimate of interest
	gen SE_theta_g = .  //SE of estimate of interest
	gen theta_g_d = . 	//estimate obtained when adding one observation
	gen LS_g = .		//local sensitivity
	
	levelsof cell, local(cells)  //save list of cells in local

	*Compute true statistic, SE, and local sensitivity for each cell
	foreach g of local cells { 
		
		*Compute true statistic: 
		qui reg kid_rank parent_rank  if cell == `g'

		*Save statistic and SE:
		replace theta_g = _b[parent_rank]*0.25 + _b[_cons] if cell == `g'
		lincom _b[parent_rank]*0.25 + _b[_cons]
		replace SE_theta_g = r(se) if cell == `g'

		*Compute Local Sensitivity (LS) for each cell:
			
			*Add one additional observation in the cell:
			count
			local additional_obs = r(N) + 1
			set obs `additional_obs'
			replace cell = `g' if  _n == `additional_obs'
			
			*loop over 4 corners of the rank-rank space (0,0), (0,1), (1,0), (1,1)
			replace LS_g = 0 if cell == `g' 
			forvalues i = 0/3 {
			
				replace parent_rank = floor(`i'/2) if _n==`additional_obs'
				replace kid_rank    = mod(`i',2)   if _n==`additional_obs'

				qui reg kid_rank parent_rank if cell==`g'
				replace theta_g_d = _b[parent_rank]*0.25 + _b[_cons] if cell == `g'
				
				*compute LS as the max absolute difference between theta_g_d and theta_g:
				replace LS_g = abs(theta_g_d - theta_g) if abs(theta_g_d - theta_g) > LS_g & cell == `g'
				
				}
		
		drop if  _n==`additional_obs'
	}

* ------------------------------------------------------------------------------
* 2. Compute Maximum Observed Sensitivity (chi)
* ------------------------------------------------------------------------------
	
	*compute number of observations per cell
	bys cell: gen N_g = _N
	
	*compute N_g * LS_g
	gen N_g_LS_g = N_g * LS_g
	
	*Find it's max:
	egen chi = max(N_g_LS_g)
	
* ------------------------------------------------------------------------------
* 3. Determine Privacy Parameter (epsilon)
* ------------------------------------------------------------------------------
	
	*Collapse at the cell level: 
	collapse theta_g SE_theta_g N_g chi, by(cell)
	
	*Compute MSE for multiple epsilons: 
	set seed 419
	forval epsilon = 1(1)10 {
		
		*Generate several draws of random noise omega from Normal or Laplace distribution:
		local draws = 500
		forval d=1(1)`draws' {
		
		if "`omega_dist'"=="normal"		gen omega_`d' = rnormal(0, 1)       		// Normal
		if "`omega_dist'"=="laplace"	gen omega_`d' = rlaplace(0, 1/sqrt(2)) 		// Laplace
		
		*Compute noise infused estimate (as in step 4):
		gen noise_infused_theta_g_`d' = theta_g + sqrt(2)*(chi / (`epsilon' * N_g)) * omega_`d'
		gen diff_true_noise_`d' = (noise_infused_theta_g _`d' - theta_g) ^ 2
		}
	
	*Compute Mean Squared Error
	egen MSE_eps_`epsilon' = rowmean(diff_true_noise_*)
	
	drop omega_*  noise_infused_theta_g_* diff_true_noise_*
	}
	
	*Plot Mean Squared Error as a function of epsilon: 
	preserve
	collapse MSE_eps_*
	gen id=1
	reshape long MSE_eps_, i(id) j(epsilon)
	twoway (connected MSE_eps_ epsilon), ///
	xtitle("Epsilon") ytitle( "Mean Squared Error (Average)" " ")
	restore
	
	drop MSE_eps_*
	
	*Say that, after discussing with the data provider, you choose an epsilon of 4: 
	local epsilon = 4
	
* ------------------------------------------------------------------------------
* 4. Add Random Noise to Estimate
* ------------------------------------------------------------------------------
	
	*A. Add random noise: 
	set seed 5711
	if "`omega_dist'"=="normal"		gen omega = rnormal(0, 1)       		// Normal
	if "`omega_dist'"=="laplace"	gen omega = rlaplace(0, 1/sqrt(2)) 		// Laplace
	
	gen noise_infused_theta_g = theta_g + sqrt(2)*(chi / (`epsilon' * N_g))*omega
	
	*B.Compute standard error of noised-infused estimate: 
	gen SE_noise_infused_theta_g = sqrt(SE_theta_g^2 + 2*((chi / (`epsilon' * N_g))^2))

	* ! Repeat steps 1 - 3 treating SE_noise_infused_theta_g as the relevant parameter (i.e. theta)
		*(For the sake of brevity that is omitted in this example do-file)
	
	*C. Construct noise infused estimates of the counts in each cell: 
	gen noise_infused_N_g = N_g + sqrt(2)*(omega / `epsilon')
	
	*D. Compute summary statistics of amount of noise added:
	
	*Standard deviation of noise in cell of interest: 
	gen SD_noise_g = sqrt(2) * (chi / (`epsilon' * N_g))
	sum SD_noise_g

	*Share of variance due to noise: 
	*total variance
	sum noise_infused_theta_g
	local total_variance = `r(Var)'

	*noise variance: 
	gen Var_noise_g = SD_noise_g ^ 2 
	sum Var_noise_g
	local noise_var = `r(mean)'

	local share_noise_variance = `noise_var' / `total_variance'
	dis %4.3f `share_noise_variance'

* ------------------------------------------------------------------------------
* 5. Release Noise-Infused Estimates
* ------------------------------------------------------------------------------
	
	*Release only noise-infused estimates:
	drop theta_g SE_theta_g N_g omega chi Var_noise_g
	
	*Export dataset:
	export excel using "example_cell_public_estimates_p25.xlsx", replace firstr(var)
