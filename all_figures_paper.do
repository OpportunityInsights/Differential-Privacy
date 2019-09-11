********************************************************************************
* Figures of: A Practical Method to Reduce Privacy Loss
* when Disclosing Statistics Based on Small Samples
********************************************************************************

clear all
set maxvar 10000
global title "title("")"

/* Need to set these paths to the files in the replication package (see README) */
global figs "C:\Users\nit518\Opportunity Insights Dropbox\Research files\outside\Differential_Privacy\paper\aea_short_paper\package\"
local atlas_data_filepath "C:\Users\nit518\Opportunity Insights Dropbox\Research files\outside\Differential_Privacy\paper\aea_short_paper\package\tract_outcomes_early_dta_dp.dta"
local covariates_filepath "C:\Users\nit518\Opportunity Insights Dropbox\Research files\outside\Differential_Privacy\paper\aea_short_paper\package\tract_covariates.dta"
****************************************************************************
*Figure 1: Local Sensitivity
****************************************************************************

*Define the tract to have 20 children
set obs 20
set seed 45

*Generate parents' income ranks up to p99
g par = runiform()*0.99

*Generate their kid's income rank
g kid = max(min(0.01 + 0.85 * par + 0.03* rnormal(),1),0)

*Add a 21st child in the tract
set obs 21

*Set this child to be an outlier with low parent rank but high child rank
replace par = 0 in 21
replace kid = 1 in 21


twoway scatter kid par || lfit kid par
qui reg kid par

*Generate the p25 prediction with this outlier added
g pred_25 = _b[par]*0.25 + _b[_cons]
di _b[par]*0.25 + _b[_cons]

* Calculate sensitivity
g sens = 0
g new_pred_25 = .
g n_rep = .
forvalues i = 1/21 {
	preserve
	replace par = . if _n==`i'
	replace kid = . if _n==`i'
	reg kid par
	restore
	local new_est = _b[par]*0.25 + _b[_cons]
	if abs(`new_est' - pred_25)>sens{
	replace sens = abs(`new_est' - pred_25)
	replace new_pred_25 = `new_est'
	replace n_rep = `i'
	}
	qui sum sens
	global sens = r(mean)
	}
di sens
di n_rep

*This is the prediction with the addtion of the outlier
sum pred_25
local pred_25= round(r(mean), .001)

*This is the prediction that occurs when the outlier is removed to return to original data
sum new_pred_25
local new_pred_25= round(r(mean), .001)

local diff : dis %4.3fc round(abs(`pred_25' - `new_pred_25'), .001)

twoway (scatter kid par if _n!=21, msize(small) text( .1 .8 "t = `new_pred_25'", size(vsmall)) text( .8 .30 "Local Sensitivity = `diff'", size(vsmall))) ///
(lfit kid par, range(0 1) lcolor(maroon) lpattern(dash) ylabel(0 "0.0" .2 "0.2" .4 "0.4" .6 "0.6" .8 "0.8" 1 "1.0") xlabel(0 "0.0" .2 "0.2" .4 "0.4" .6 "0.6" .8 "0.8" 1 "1.0")) ///
(lfit kid par if _n!=21, range(0 1) lcolor(navy) xtitle("Parents' Income Rank") ytitle("Child's Income Rank")) ///
(scatteri 1 0, mcolor(navy) msymbol(circle_hollow) msize(small)), ///
$title legend(order(3 "OLS Estimate in Original Data" 2 "OLS Estimate with Addition of Outlier") size(small))

graph export "${figs}/fig_1.wmf", replace

****************************************************************************
*Figure 2: MOSE
****************************************************************************

clear all

*Set parameters for simulation:
set seed 42
global N_per = 3
global N_low = 20
global N_high = 300
global par_range = 0.75
global slope_low = 0
global slope_high = 0.5
global int_low = 0
global int_high = 0.5
global noise = 0.18
global obs = $N_per * ($N_high - $N_low)

*Create dataset to store local sensibility for each cell:
set obs $obs
g tract = _n
g N = .
g sens = .

*For each cell, generate data and compute local sensibility:
forvalues n = 1/$obs {
	preserve
	clear

	*set number of observatons in cell:
	global ob = floor((`n' - 1)/$N_per) + $N_low
	set obs $ob

	*generate data:
	qui g par_range = runiform()*$par_range
	qui sum par_range
	global pr = r(mean)
	qui g start = runiform()*(1-$pr)
	qui sum start
	global st = r(mean)
	qui g par = runiform()*$pr + $st
	drop par_range start

	qui g intercept = runiform()*($int_high - $int_low ) + $int_low
	qui g slope = runiform()*($slope_high - $slope_low) + $slope_low
	qui g kid = max(min(intercept + slope * par + $noise * rnormal(),1),0)

	*Compute p25 prediction:
	qui reg kid par
	global p25 = _b[par]*0.25 + _b[_cons]

	*Compute local sensibility:
	global ob1 = $ob + 1
	set obs $ob1
	qui g sens = 0
	qui g new_pred_25 = .
		forvalues i = 0/3 {
		qui replace par = floor(`i'/2) if _n==($ob + 1)
		qui replace kid = mod(`i',2) if _n==($ob + 1)
		qui reg kid par
		*di _b[par]*0.25 + _b[_cons]
		qui replace sens = abs(_b[par]*0.25 + _b[_cons] - $p25) if abs(_b[par]*0.25 + _b[_cons] - $p25)>sens
		qui replace new_pred_25 = _b[par]*0.25 + _b[_cons] if abs(_b[par]*0.25 + _b[_cons] - $p25)==sens
		qui sum sens
		global sens = r(mean)
		}
	restore

	*Store local sensibility:
	qui replace sens = $sens if _n== `n'
	qui replace N = $ob if _n == `n'
	}

*Compute MOS:
g cap = N*sens
qui sum cap
g ub = r(max)/N
local mos = round(r(max), .01)

gen N_rounded= round(N,5)

sum sens if N_rounded ==50
local max_at_50 = round(r(max), .01)
di `max_at_50'

twoway (scatter sens N_rounded, yscale(log range(.015 .6)) xscale(log) msize(tiny) xtitle("Number of Individuals in Tract (Log Scale)") ytitle("Sensitivity of p25 Predictions (Log Scale)") xlabel(20 50 100 200 300) ylabel(0.02 0.05 0.1 0.2 0.5, format(%5.2f)) legend(off)) ///
(scatter ub N, m(i) c(l)), text( .4 200 "{&delta} = `mos'", size(medium)) $title

graph export "${figs}/fig_2.wmf", replace

****************************************************************************
*Figure 3a: Teenage Birth and Two-Parent Share (Noise-Infused)
****************************************************************************

	* switches
	local var			teenbrth
	local g 			female

	* load data
	use "`atlas_data_filepath'", clear

	//renvars state county tract, suff(10)

	merge 1:1 state county tract using ///
		"`covariates_filepath'", ///
		nogen keep(3) keepusing(pop* poor_share2000 singleparent_share2000)

	* winsorize to deal with outliers
	foreach race in white black {
		qui: su `var'_`race'_`g'_p25 [w=`var'_`race'_`g'_n], d
		qui: replace `var'_`race'_`g'_p25 = `r(p99)' if ///
			`var'_`race'_`g'_p25 > `r(p99)' & ///
			`var'_`race'_`g'_p25  !=.
		qui: replace `var'_`race'_`g'_p25 = `r(p1)' if ///
			`var'_`race'_`g'_p25 < `r(p1)' & ///
			`var'_`race'_`g'_p25  !=.
	}

	tempfile data
	save `data'

	* suppress data and save tempfiles
	* loop over different intervals for numerator suppression
	foreach i of numlist 0(1)10 {
		use `data', clear
		foreach race in black white {
			qui: gen num_`var'_`race'_`g'_`i' = `var'_`race'_`g'_p25 * ///
				`var'_`race'_`g'_n * frac_below_median_`race'_`g'
			qui: gen om_`var'_`race'_`g'_`i' = 0
			local z = `i'+0.5
			if `i' > 0 {
				qui: replace om_`var'_`race'_`g'_`i' = 1 if ///
					inrange(num_`var'_`race'_`g'_`i', 0.51, `z')
				qui: replace num_`var'_`race'_`g'_`i' = . if ///
					inrange(num_`var'_`race'_`g'_`i', 0.51, `z')
				}
			gen frac_`var'_`race'_`g'_`i' = num_`var'_`race'_`g'_`i' / ///
				(`var'_`race'_`g'_n * frac_below_median_`race'_`g')
			tempfile collapsed_`i'
			save `collapsed_`i''
			}
		}

	keep state county tract frac_* om_*

	* merge tempfiles of different suppressions
	use `collapsed_0', clear
	foreach i of numlist 1(1)10 {
		merge 1:1 state county tract using `collapsed_`i'', nogen
		}

	* merge initial data
	merge 1:1 state county tract using `data', ///
		keepusing(*_n) nogen assert(3)

	* rescale rate with 100 to get %
	replace  frac_teenbrth_black_female_0 = frac_teenbrth_black_female_0 * 100
	replace  frac_teenbrth_black_female_4 = frac_teenbrth_black_female_4 * 100

	replace  singleparent_share2000 = singleparent_share2000 * 100


	*Regressions
	preserve
	drop if mi(frac_teenbrth_black_female_0) | mi(singleparent_share2000) | mi(kid_black_female_blw_p50_n)
	keep if poor_share2000 < 0.07
	reg frac_teenbrth_black_female_0 singleparent_share2000 [w=kid_black_female_blw_p50_n]
	local beta_0 : di %4.3f _b[singleparent_share2000]
	local se_0 : di %4.3f _se[singleparent_share2000]
	restore


	preserve
	drop if mi(frac_teenbrth_black_female_4) | mi(singleparent_share2000) | mi(kid_black_female_blw_p50_n)
	keep if poor_share2000 < 0.07
	reg frac_teenbrth_black_female_4 singleparent_share2000 [w=kid_black_female_blw_p50_n]
	local beta_4 : di %4.3f _b[singleparent_share2000]
	local se_4 : di %4.3f _se[singleparent_share2000]
	restore


	*BINSCATTERS:

		*1. teenage birth rate and singleparent share (black female)
		binscatter frac_teenbrth_black_female_0 singleparent_share2000 if poor_share2000 < 0.07 [w=kid_black_female_blw_p50_n], ///
		xtitle("Single Parent Share in Tract (2000)") ytitle("Teen Birth Rates for Black Women" "Raised in Low-Income Families") ///
		$title ylabel(32 "32%" 34 "34%" 36 "36%" 38 "38%" 40 "40%") ///
		xlabel(10 "10%" 20 "20%"  30 "30%" 40 "40%" 50 "50%") xscale(range(8 53)) yscale(range(35.3 41)) ///
		text(35 43 "Slope = `beta_0'") text(34.2 45.5 "(`se_0')")
		graph export "${figs}/fig_3a_v_final.wmf", replace

****************************************************************************
*Figure 3b: Teenage Birth and Two-Parent Share (Count-Suppressed Data)
****************************************************************************

		binscatter frac_teenbrth_black_female_4 singleparent_share2000 if poor_share2000 < 0.07 [w=kid_black_female_blw_p50_n], ///
		xtitle("Single Parent Share in Tract (2000)") ytitle("Teen Birth Rates for Black Women" "Raised in Low-Income Families") ///
		text(35 43 "Slope = `beta_4'") text(34.2 45.5 "(`se_4')") ///
		ylabel(32 "32%" 34 "34%" 36 "36%" 38 "38%" 40 "40%") xlabel(10 "10%" 20 "20%"  30 "30%" 40 "40%" 50 "50%") ///
		$title xscale(range(8 53)) yscale(range(35.3 41))
		graph export "${figs}/fig_3b_v_final.wmf", replace
