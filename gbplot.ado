*
*		PROGRAM GBPLOT
*
*		version 2.1
*		Date: 10/20/26
*		Allyson Benton, Soren Jordan and Andrew Philips
*
*
* Changes from v2.0:
*	[v2.1 N1] Added three simulation options (no delta method allowed):
*	          fcast_epu         expected path uncertainty: each bootstrap b recurses on 
*	                            its OWN sigma2_{t-1,b}, not E[sigma2_{t-1,b}]
*	          fcast_arch        epu + ARCH innovations ~ N(0,1)
*	          fcast_arch_resid  epu + ARCH innovations ~ whitened residuals saved by
*	                            gbmake, fcast_arch (now gb_innov in using file)
*
* Changes from v1.4:
*	[v2 N2]  Default plot now spikeplot---areaplot now option
*	[v2 S1]  speed: results collected in matrix -> dataset once. preivously redundancy. 
*				same with nlcom blocks.
*	[v2 S2]  speed: using path keeps 2 working variables instead of 1 per period.
*	[v2 F1]  mean(Z) uses e(sample), respects if/in, not full dataset
*	[v2 F2]  shockvariable() checked against the het() variables/using dataset at start
*	[v2 F3]  nlcom now quiet in every period
*	[v2 F4]  Matrices/scalars now tempnames
*	[v2 F5]  burnin(0) and psv_usematrix(1) previously broke program
*	[v2 F6]  psv_usematrix without usematrix now ignored, and shocktime <= range 
*	           only enforced when shock is actually applied.
*	[v2 F7]  moved 'replace' override to top to avoid grumpy user. allow spaces
*			 in file paths too.
*	[v2 F8]  related to S2, using dataset previously read twice, now once 
*	[v2 F9]  duplicate plotting code now new helper function _gbplot_draw
*
* Changes from v1.3:
*	[Fix D1] Add level(`level') to all nlcom step2 calls
*	[Fix D2] Adjust preshocktime for burnin reindexing in using path PSV block
*	[Fix D3] Same preshocktime adjustment for delta path PSV block.
*	[Fix D4] noplot informational note changed to `di in g`
*	[Fix D5] If we use ,post in nlcom should replay saved garch results
*
* Changes from v1.2:
* 	[Fix C1] .dta was appended always, even if user specified .dta
*	[Fix C2] nlcom fix attempt when it fails
*	[Fix C3] add in a perm option, else shock is temporary
*	[Fix C4] replaced r(table) with r(b)/r(V) for delta method
*
* Changes from v1.1:
*	[FIX B2] Ensured we're using locals, not scalars for nlcom that seemed to
*				cause error under certain OS's. Also upped iterations = 1000
*	[FIX B3] removed garch*0 when t=1 b/c that was reading as garch*. for
*				certain OS's.
*
* Changes from v1.0:
*   [FIX B4] Added user_burnin flag to avoid spurious "burnin ignored" warning
*            when usematrix specified but burnin not set
*   [FIX B6] Added guard against preshocktime < 1 in PSV calculation
*   [FIX B7] Replaced data variable reference inside nlcom with explicit scalar
*            (su + local) for the GARCH lag term in the delta method
*
* -------------------------------------------------------------------------
* -------------------------------------------------------------------------
* -------------------------------------------------------------------------
cap prog drop gbplot
prog define gbplot, rclass
syntax [if] [in], [SHOCKVARiable(string) SHOCKVALue(numlist max=1)]  ///
    [delta using(string) burnin(numlist integer >=0 max=1)] ///
    [usematrix(string) shocktime(numlist integer >1 max=1)] ///
    [range(numlist integer >1 max=1) level(numlist max=1)] ///
    [customvalues(numlist) psv SAVing(string)] ///
	[noplot psv_usematrix(numlist integer >=1 max=1)] ///
	[areaplot spikeplot replace PERManent] ///
	[fcast_epu fcast_arch fcast_arch_resid]

* delta: uses arch call and uses delta method, not using
* using: dataset with simulated data to plot impulse responses from
* burnin: number of burn in observations to drop from simulated data
* usematrix: name of matrix with shocks to apply to simulated data
* shockvariable: variable to shock
* shockvalue: size of shock to apply
* shocktime: time period to apply shock (default is 5)
* range: number of periods to plot (default is 15)
* level: significance level for confidence intervals (default is 0.95).
* customvalues: list of values corresponding exactly to variables in the het equation
*	(for delta method) or all variables in the simulated dataset (excluding __cons, arch1, garch1)
*	to set when simulating shocks
* psv: express plot as % of pre-shock value
* saving: if exporting simulated results for additional plotting is desired
* noplot: suppress showing plot (advised to save then since nothing gets stored to memory)
* psv_usematrix: when using usematrix along with psv, specify the time period to use as
* 	pre-shock value (default is 1)
* areaplot: show areaplot instead of the default spikeplot
* spikeplot: (default; still accepted for backward compatibility, but ignored otherwise)
* replace: overwrite using saving if file exists
* permanent: default is impulse, this option allows permanent shock
* fcast_epu: [v2.1] (using only) expected path uncertainty. For EACH bootstrap draw b:
*		s2_{1,b} = exp(Z_1 g_b)
*		s2_{t,b} = exp(Z_t g_b) + garch_b * s2_{t-1,b}
*	i.e. the lag is draw b's own previous variance, not the mean across draws (default).
*	Plotted: median and percentile CIs of s2_{t,b} across b.
* fcast_arch: [v2.1] (using only) epu plus ARCH innovations. Draw u_{t-1,b} ~ N(0,1), set
*	e_{t-1,b} = sigma_{t-1,b} * u_{t-1,b}, so
*		s2_{t,b} = exp(Z_t g_b) + arch_b * s2_{t-1,b} * u^2 + garch_b * s2_{t-1,b}
*	Requires arch1 in using file.
* fcast_arch_resid: [v2.1] (resid bs only) same as fcast_arch, but u resampled from
*	whitened residuals saved in variable gb_innov by gbmake (residual bootstrap), fcast_arch.

/*
general approach:
	1. First decide whether to use delta method or simulated dataset
		- delta method:
			- ensure arch previously run, ensure GARCH(1,0) or GARCH(1,1) only, ensure het() specified
			- get parameters and standard errors from e()
			- simulate shocks via delta method
		- simulated dataset:
			- load in simulated dataset
			- get means/modes of all variables in the dataset (excluding __cons, arch1, garch1)
			- simulate shocks by changing specified variable at specified time period
	2. Then plot impulse responses based on the simulated shocks
		- default: plot all variables over default range of periods (15)
		- custom: plot specified variables over specified range of periods
		- confidence intervals based on specified significance level (default 0.95)

	## Additional options:
		- usematrix allows for custom matrix to be specified for the shocks
		- customvalues allow for setting specific values for variables when simulating the shocks
		- plot types: default spikeplot, areaplot optional
		- fcast_epu / fcast_arch / fcast_arch_resid: bootstrap path-specific variance recursion,
	 		with ARCH innovations (using only)
*/

* ------------ Checks -----------------------
* assert one of either using or delta
* catch if nothing is specified:
if "`using'" == "" & "`delta'" == "" {
	di in r _n "Error: must specify either 'using' or 'delta'"
	exit 198
}

* If using, assert delta not specified:
if "`using'" != "" {
	cap confirm file "`using'"
	if _rc {
    	di in r _n "Error: file '`using'' not found."
    	exit 601
	}
	if "`delta'" != "" {
		di in r _n "Error: 'using' and 'delta' must be specified exclusively; you provided both"
		exit 198
	}
}

* [v2.1 N1] fcast_* checks: at most one, using only
local fmode ""
local nf = ("`fcast_epu'" != "") + ("`fcast_arch'" != "") + ("`fcast_arch_resid'" != "")
if `nf' > 1 {
	di in r _n "Error: specify at most one of 'fcast_epu', 'fcast_arch', 'fcast_arch_resid'."
	exit 198
}
if `nf' == 1 {
	if "`delta'" != "" { // assert delta not used for fcast
		di in r _n "Error: 'fcast_*' options are only available with 'using' (bootstrap draws), not 'delta'"
		exit 198
	}
	if "`fcast_epu'" != "" local fmode "epu"
	if "`fcast_arch'" != "" local fmode "arch"
	if "`fcast_arch_resid'" != "" local fmode "resid"
}

* make sure shockvariable and shockvalue specified (unless usematrix specified)
if "`usematrix'" == "" {
	if "`shockvariable'" == "" {
		di in r _n "Error: must specify variable to shock via 'shockvariable' option (unless using 'usematrix')."
		exit 198
	}
	if "`shockvalue'" == "" {
		di in r _n "Error: must specify shock value via 'shockvalue' option (unless using 'usematrix')."
		exit 198
	}
}

* assert only 1 variable specified for shockvariable:
if "`shockvariable'" != "" {
	local shockvariable_count: word count `shockvariable'
	if `shockvariable_count' != 1 {
		di in r _n "Error: Only one variable may be specified using 'shockvariable'. See 'usematrix' for multiple shocks."
		exit 198
	}
}

* is there a shock applied by gbplot itself (as opposed to one in a usematrix)?
local hasshock = ("`shockvariable'" != "" & "`usematrix'" == "")

* burnin: [FIX B4] track whether user specified it; defaults are 0 for usematrix, 20 otherwise
local user_burnin = ("`burnin'" != "")
if !`user_burnin' {
	if "`usematrix'" != "" {
		loc burnin 0
	}
	else {
		loc burnin 20
	}
}

* level: number is logical, and tag for below
local user_level = 0
if "`level'" != "" {
	if `level' <= 0 | `level' >= 1 {
		di in r _n "Error: 'level' must be between 0 and 1 (e.g., 0.95 = 95% CIs)."
		exit 198
	}
	local user_level = 1
	local level = round(`level'*100, .01)
}
else {
	local level = 95 // default
}

* warn if usematrix specified along with incompatible options
if "`usematrix'" != "" {
	if "`shockvariable'" != "" | "`shockvalue'" != "" {
		di in r _n "Warning: 'usematrix' option specified, 'shockvariable' and 'shockvalue' options ignored."
	}
	if "`shocktime'" != "" {
		di in r _n "Warning: 'usematrix' option specified, 'shocktime' option ignored."
	}
	if "`range'" != "" {
		di in r _n "Warning: 'usematrix' option specified, 'range' option ignored."
	}
	if "`customvalues'" != "" {
		di in r _n "Warning: 'usematrix' option specified, 'customvalues' option ignored."
	}
}
* set defaults for range and shocktime
if "`range'" == "" {
	local range = 15
}
if "`shocktime'" == "" {
	local shocktime = 5
}
* assert that shocktime is within range (only matters if gbplot applies a shock)
if `hasshock' & `shocktime' > `range' {
	di in r _n "Error: 'shocktime' must be less than or equal to 'range'."
	exit 198
}

* set totalsimtime and startshock based on burnin
local totalsimtime = `range' + `burnin'
local startshock = `shocktime' + `burnin' // when to start shock

* throw warnings/errors for psv_usematrix
if "`psv_usematrix'" != "" {
	if "`usematrix'" == "" {
		di in g _n "Warning: 'psv_usematrix' option specified without providing 'usematrix'; 'psv_usematrix' ignored."
		local psv_usematrix ""
	}
	else {
		* if psv_usematrix specified, set psv accordingly, regardless of whether psv was called
		local psv "psv"
		* error if it's earlier than burnin stops
		if `psv_usematrix' <= `burnin' {
			di in r _n "Error: 'psv_usematrix' is specifed during the burnin period but must come after t = `burnin'."
			exit 198
		}
	}
}

* psv needs a shock or usematrix
if "`psv'" != "" & "`usematrix'" == "" & !`hasshock' {
	di in r _n "Error: 'psv' requires a shock (shockvariable/shockvalue) or usematrix."
	exit 198
}

* [v2 N2] plot type
if "`areaplot'" != "" & "`spikeplot'" != "" {
	di in r _n "Error: specify only one of 'areaplot' or 'spikeplot'."
	exit 198
}
local plottype = cond("`areaplot'" != "", "area", "spike")

* strip .dta from saving if user added it:
if "`saving'" != "" {
	if lower(substr("`saving'", -4, .)) == ".dta" {
		local saving = substr("`saving'", 1, length("`saving'") - 4)
	}
	* [v2 F7] check the output file NOW
	cap confirm file "`saving'.dta"
	if _rc == 0 & "`replace'" == "" {
		di in r _n "Error: File `saving'.dta already exists (use 'replace' to overwrite)"
		exit 602
	}
}

* observations for means of Z: user's if/in, now intersected e(sample) [v2 F1]
marksample touse, novarlist
if "`e(cmd)'" == "arch" {
	qui replace `touse' = 0 if !e(sample)
}

tempname M R Q bm Vm // matrix/scalar tempnames
preserve


* ------------ Delta Method -----------------
if "`delta'" != "" {

	* assert arch previously run (same checks as in gbmake)
	if "`e(cmd)'" != "arch"	{
		di in r _n "Error: you must first run arch before using the 'delta' option."
		exit 198
	}

	* ensure GARCH(1,0) or GARCH(1,1) only, and record for below
	if "`e(garch)'" == "1" {
		if "`e(arch)'" == "1" {
			local is_arch = 1
		}
		else if "`e(arch)'" == "0" {
			local is_arch = 0
		}
		else {
			di in r _n "Error: only ARCH(1) or ARCH(0) process currently supported"
			exit 198
		}
	}
	else {
		di in r _n "Error: only GARCH(1,1) or GARCH(1,0) process currently supported"
		exit 198
	}
	* ensure no weird GARCH flavors:
	if "`e(saarch)'" != "" | "`e(tarch)'" != "" | "`e(aarch)'" != "" | "`e(narch)'" != "" | ///
 	"`e(narchk)'" != "" | "`e(abarch)'" != "" | "`e(atarch)'" != "" | "`e(sdgarch)'" != "" | ///
 	"`e(earch)'" != "" | "`e(parch)'" != "" | "`e(tparch)'" != "" | "`e(aparch)'" != "" | ///
  	"`e(nparch)'" != "" | "`e(nparchk)'" != "" | "`e(pgarch)'" != "" {
		di in r _n "Error: only GARCH(1,1) or GARCH(1,0) process currently supported"
		exit 198
	}

	* ensure no ARIMA terms
	if "`e(ar)'" != "" | "`e(ma)'" != "" {
		di in r _n "Error: no ARIMA(p,d,q) terms currently supported"
		exit 198
	}

	* ensure a HET term is actually specified:
	if "`e(mhet)'" == "" {
		di in r _n "Error: no het() equation specified; nothing to simulate"
		exit 198
	}
	* ensure no ARCH-in-Mean, or just 1 contemporaneous
	if "`e(archm)'" != "" {
		if "`e(archm)'" != "0" {
			di in r _n "Error: only contemporaneous ARCH-in-Mean currently supported"
			exit 198
		}
	}

	* if 'nlcom , post' triggers, must restore original garch estimates
	quietly estimates store _gbplot_original_arch

	* now do the namepeeler (similar to gbmake but simpler).
	loc eqnames `e(eqnames)'
	matrix `bm' = e(b)
	local bnames : colnames `bm'
	local old_dv "`e(depvar)'" // for assertion of no DV in het eq
	local index 1 // this will let us loop through bnames (list of variable names)
	local Bvar_Z "" // list of variables in het eq (excluding _cons)
	* next, go through each word `eq_name' in eqnames:
	foreach eq_name of local eqnames {

		* first, obtain the variable name:
		loc vname: word `index' of `bnames'
		if "`eq_name'" == "HET" { // if eqnames == HET, we're in the variance eq
			if regexm("`vname'", "^[^.]*\.`old_dv'\$") { // we don't want any form of the DV in this equation
				di in r _n "Error: the variance eq. cannot contain the dependent variable"
				exit 198
			}
			else if "`vname'" != "_cons" {
				loc Bvar_Z "`Bvar_Z' `vname'" // needed for below
			}
		}
		local index = `index' + 1
	}
	local Bvar_Z = trim("`Bvar_Z'")
	local het_count: word count `Bvar_Z'

	* [v2 F2] shockvariable must exist in het() variables
	if `hasshock' & !`: list shockvariable in Bvar_Z' {
		di in r _n "Error: shockvariable '`shockvariable'' is not a variable in the het() equation (`Bvar_Z')."
		exit 198
	}

	* next use delta method to loop through from 1 to range (or matrix rows), implementing any shocks

	* If usematrix is specified, derive totalsimtime from its rows and ensure column count matches het variables
	if "`usematrix'" != "" {
		cap matrix `M' = `usematrix'
		if _rc {
			di in r _n "Error: could not load matrix specified in usematrix(); provide a valid Stata matrix or literal."
			exit 198
		}
		_gbplot_chkmat `M' `het_count' "het() variable"
		local totalsimtime = rowsof(`M')
		local startshock = 1
		noi di "usematrix provided: rows=`=rowsof(`M')', cols=`=colsof(`M')'. Ignoring shock/range/customvalues options."
		noi di "Het values will be read from usematrix() each time period; Bvar_Z list is: `Bvar_Z'"
	}
	else {
		* Build het equation expression string (for use in nlcom).
		* If customvalues are provided, use them in the order of the het() variables; otherwise use means
		local het_expr "_b[HET:_cons]"
		if "`customvalues'" != "" {
			local custom_count: word count `customvalues'
			if `custom_count' != `het_count' {
				di in r _n "Error: number of customvalues (`custom_count' specified) must equal number of variables in het(), which is `het_count'."
				exit 198
			}
			local idx = 1
			foreach var of local Bvar_Z {
				local thisval: word `idx' of `customvalues'
				local het_expr "`het_expr' + _b[HET:`var']*(`thisval')"
				local idx = `idx' + 1
			}
		}
		else {
			foreach var of local Bvar_Z {
				quietly summarize `var' if `touse', meanonly
				* [v2.1] missing mean (e.g., data reloaded after arch, so e(sample) is empty)
				if missing(r(mean)) {
					di in r _n "Error: mean of '`var'' is missing (no usable observations). If you reloaded the data after running arch, re-run arch first, or pass customvalues()."
					exit 198
				}
				local het_expr "`het_expr' + _b[HET:`var']*`r(mean)'"
			}
		}
		noi di "het expression is: `het_expr'"
		noi di "Bvar_Z list is: `Bvar_Z'"
	}

	* GARCH coefficient expression (we've already asserted it's GARCH(1) only)
	local garch_expr "_b[ARCH:L.garch]"
	local z = invnormal(0.5 + `level'/200) // critical value for the delta-method CIs

	* Main simulation loop: loop through time (range or usematrix rows)
	* results go in a matrix: point estimate, lower, upper
	matrix `R' = J(`totalsimtime', 3, .)
	local pastsigma2 = 0
	forv t = 1/`totalsimtime' {

		* period-specific het expression
		if "`usematrix'" != "" {
			local het_t "_b[HET:_cons]"
			local colidx = 1
			foreach var of local Bvar_Z {
				local curval = `M'[`t', `colidx']
				local het_t "`het_t' + _b[HET:`var']*(`curval')"
				local colidx = `colidx' + 1
			}
		}
		else {
			local het_t "`het_expr'"
			* shock at shocktime; continues afterwards only if permanent
			if `hasshock' & (`t' == `startshock' | (`t' > `startshock' & "`permanent'" != "")) {
				local het_t "`het_t' + _b[HET:`shockvariable']*(`shockvalue')"
			}
		}

		* E[sigma2_t] = exp(Z g) + garch*E[sigma2_(t-1)]. At t=1 there is no lag.
		* [FIX B7] the lag is the previous point estimate, passed in as a number
		local expr "exp(`het_t')"
		if `t' > 1 {
			local expr "`expr' + `garch_expr'*`pastsigma2'"
		}

		* Fix [C2]: catch for nlcom fail
		capture qui nlcom `expr', level(`level') iterate(1000)
		if _rc {
			qui nlcom (step1: 1e+3 + `expr'), level(`level') iterate(1000) post
			qui nlcom (step2: (_b[step1] - 1e+3)), level(`level')
			quietly estimates restore _gbplot_original_arch // restore OG garch estimates if ', post' triggered
		}

		* Fix C4: r(b)/r(V)
		matrix `bm' = r(b)
		matrix `Vm' = r(V)
		matrix `R'[`t', 1] = `bm'[1,1]
		matrix `R'[`t', 2] = `bm'[1,1] - `z'*sqrt(`Vm'[1,1])
		matrix `R'[`t', 3] = `bm'[1,1] + `z'*sqrt(`Vm'[1,1])
		local pastsigma2 = `bm'[1,1]
	}

	* matrix -> dataset, dropping the burnin periods
	matrix colnames `R' = sigma2_pred_ sigma2_ll_ sigma2_ul_
	qui clear
	qui svmat double `R', names(col)
	if `burnin' > 0 {
		qui drop in 1/`burnin'
	}
	qui gen time = _n

	local cvar "sigma2_pred_"
	local bands "sigma2_ll_ sigma2_ul_"
	local ytitle "E[Conditional Error Variance]"
	local xtitle "time"
	local note "`level'% confidence intervals shown"

	quietly estimates drop _gbplot_original_arch // clean up OG garch est from memory
}
* ----------------------------------




* ------------ Using Simulated Data, `using` -----------------
else {
	* Get variable names from using dataset, excluding __cons, arch1, garch1
	* [v2 F8] describe using: read var names without loading
	quietly describe using "`using'", varlist
	local allvars `r(varlist)'
	local usingvars ""
	* check to make sure const and garch1 are in using (don't worry about arch)
	local __chk_cons 0
	local __chk_garch 0
	local has_arch1 0
	local has_innov 0
	foreach var of local allvars {
		if "`var'" == "__cons" {
			local __chk_cons 1
		}
		else if "`var'" == "garch1" {
			local __chk_garch 1
		}
		else if "`var'" == "arch1" {
			local has_arch1 1
		}
		else if "`var'" == "gb_innov" { // [v2.1] whitened residuals saved by gbmake, fcast_arch; not a Z variable
			local has_innov 1
		}
		else {
			local usingvars `usingvars' `var'
		}
	}

	* assert checks:
	if `__chk_cons' != 1 {
		di in r _n "Error: 'using' does not contain a vector of bootstrapped constants, named __cons"
		exit 198
	}
	if `__chk_garch' != 1 {
		di in r _n "Error: 'using' does not contain a vector of bootstrapped GARCH(1) parameters, named garch1"
		exit 198
	}
	local using_count: word count `usingvars'

	* [v2.1 N1] requirements for the ARCH-innovation options
	if inlist("`fmode'", "arch", "resid") & !`has_arch1' {
		di in r _n "Error: 'fcast_arch' / 'fcast_arch_resid' need ARCH(1) parameters (arch1) in 'using'; this model has none."
		exit 198
	}
	if "`fmode'" == "resid" & !`has_innov' {
		di in r _n "Error: 'fcast_arch_resid' needs the whitened innovations variable gb_innov in 'using'. Re-run gbmake with btype(residual) and fcast_arch."
		exit 198
	}

	* [v2 F2] shockvariable must be one of the simulated-data variables
	if `hasshock' & !`: list shockvariable in usingvars' {
		di in r _n "Error: shockvariable '`shockvariable'' is not a variable in the simulated dataset (`usingvars')."
		exit 198
	}

	* If usematrix is provided, validate dimensions and derive simulation length
	if "`usematrix'" != "" {
		cap matrix `M' = `usematrix'
		if _rc {
			di in r _n "Error: could not load matrix specified in usematrix(); provide valid Stata matrix."
			exit 198
		}
		_gbplot_chkmat `M' `using_count' "simulated-data variable"
		local totalsimtime = rowsof(`M')
		local startshock = 1
		noi di "usematrix provided: rows=`=rowsof(`M')', cols=`=colsof(`M')'. Ignoring shock/range/customvalues options."
		loc messageidx = 1
		foreach var of local usingvars {
			noi di "	values for `var' will come from matrix column `messageidx' each period."
			loc messageidx = `messageidx' + 1
		}
	}
	else {
		* Allow overriding means with customvalues (only when not using matrix)
		* Note: means computed from the original estimation data (before switching to bootstrap data below).
		if "`customvalues'" != "" {
			local custom_count: word count `customvalues'
			if `custom_count' != `using_count' {
				di in r _n "Error: number of customvalues (`custom_count') must equal number of variables in simulated data (`using_count')."
				exit 198
			}
		}

		* now create local of constant Z-Gamma values that's going to be used in every sim
		local zgamma ""
		local idx = 1
		foreach var of local usingvars {
			if "`customvalues'" != "" {
				local thisval: word `idx' of `customvalues'
				local mean_`var' = `thisval'
				local idx = `idx' + 1
			}
			else {
				cap confirm numeric variable `var'
				if _rc {
					di in r _n "Error: variable '`var'' from the simulated dataset is not in the data in memory (needed for its mean)."
					exit 111
				}
				quietly summarize `var' if `touse', meanonly
				* [v2.1] missing mean (e.g., data reloaded after arch, so e(sample) is empty)
				if missing(r(mean)) {
					di in r _n "Error: mean of '`var'' is missing (no usable observations). If you reloaded the data after running arch, re-run arch first, or pass customvalues()."
					exit 198
				}
				local mean_`var' = r(mean)
			}
			noi di "mean of `var' is: `mean_`var''"
			local zgamma "`zgamma' + `var'*`mean_`var''"
		}
		noi di "z-gamma list is: `zgamma'"
	}

	* Switch to the using dataset
	use "`using'", clear

	* [v2.1 N1] innovation pool for fcast_arch_resid: held in Mata. gbmake pads the file with
	* missing parameter rows when there are more residuals than reps; drop those.
	if "`fmode'" == "resid" {
		mata: _gbplot_pool = st_data(., "gb_innov"); _gbplot_pool = select(_gbplot_pool, _gbplot_pool :< .)
		mata: st_local("npool", strofreal(rows(_gbplot_pool)))
		if `npool' == 0 {
			mata: mata drop _gbplot_pool
			di in r _n "Error: gb_innov in 'using' has no non-missing values."
			exit 198
		}
	}
	qui drop if missing(__cons)

	* [v2.1] warn about draws that are negative or explosive (they are used as-is)
	qui count if garch1 < 0
	if r(N) > 0 di in g _n "Warning: `r(N)' bootstrap draw(s) have a negative garch1; variance paths for these draws can go negative."
	if `has_arch1' {
		qui count if arch1 < 0
		if r(N) > 0 di in g _n "Warning: `r(N)' bootstrap draw(s) have a negative arch1; variance paths for these draws can go negative."
		qui count if arch1 + garch1 >= 1
		if r(N) > 0 di in g _n "Warning: `r(N)' bootstrap draw(s) have arch1 + garch1 >= 1 (explosive variance paths)."
	}
	else {
		qui count if garch1 >= 1
		if r(N) > 0 di in g _n "Warning: `r(N)' bootstrap draw(s) have garch1 >= 1 (explosive variance paths)."
	}

	* percentiles to compute, and names for the columns of the results matrix
	if `user_level' == 1 {
		local user_ll = (100 - `level')/2
		local user_ul = 100 - `user_ll'
		local plist "`user_ll', 50, `user_ul'"
		local cnames "esigma2_ll_ esigma2_med_ esigma2_ul_"
		local cvar "esigma2_med_"
		local bands "esigma2_ll_ esigma2_ul_"
		local note "`level'% CIs shown"
	}
	else {
		local plist "2.5, 5, 12.5, 50, 87.5, 95, 97.5"
		local cnames "esigma2_ll95_ esigma2_ll90_ esigma2_ll75_ esigma2_med_ esigma2_ul75_ esigma2_ul90_ esigma2_ul95_"
		local cvar "esigma2_med_"
		local bands "esigma2_ll95_ esigma2_ul95_ esigma2_ll90_ esigma2_ul90_ esigma2_ll75_ esigma2_ul75_"
		local note "75, 90 and 95% CIs shown"
	}
	local ncol: word count `cnames'
	if `totalsimtime' <= `burnin' {
		if "`fmode'" == "resid" mata: mata drop _gbplot_pool // [v2.1] don't leak pool on this exit
		di in r _n "Error: usematrix() has `totalsimtime' rows, which are all consumed by burnin(`burnin')."
		exit 198
	}
	matrix `Q' = J(`totalsimtime' - `burnin', `ncol', .)

	* loop through simulations and apply shocks. Two working variables only:
	* scur = this period's sigma2 (one value per bootstrap draw); slag = last period's
	tempvar scur slag uu
	qui gen double `scur' = .
	qui gen double `slag' = .
	qui gen double `uu' = . // [v2.1] ARCH innovation u for this period (fcast_arch / fcast_arch_resid)
	local pastsigma2 = 0
	forv t = 1/`totalsimtime' {

		* Z-gamma for period t
		if "`usematrix'" != "" {
			local zg ""
			local colidx = 1
			foreach var of local usingvars {
				local curval = `M'[`t', `colidx']
				local zg "`zg' + `var'*(`curval')"
				local colidx = `colidx' + 1
			}
		}
		else {
			local zg "`zgamma'"
			* only incorporate shock starting at the shocktime; continues only if permanent specified
			if `hasshock' & (`t' == `startshock' | (`t' > `startshock' & "`permanent'" != "")) {
				local zg "`zg' + `shockvariable'*(`shockvalue')"
			}
		}

		* lag term
		if "`fmode'" != "" {
			* [v2.1 N1] each bootstrap b carries its OWN previous variance (slag), s2_0 = 0 so no lag at t=1.
			*	epu:         garch_b*s2_{t-1,b}
			*	arch/resid:  (garch_b + arch_b*u^2)*s2_{t-1,b}, since e_{t-1,b}^2 = s2_{t-1,b}*u^2
			if `t' == 1 local lag "0"
			else if "`fmode'" == "epu" local lag "garch1*`slag'"
			else {
				if "`fmode'" == "arch" qui replace `uu' = rnormal()
				else mata: st_store(., "`uu'", _gbplot_pool[ceil(runiform(st_nobs(), 1) :* rows(_gbplot_pool))])
				local lag "(garch1 + arch1*`uu'^2)*`slag'"
			}
		}
		else {
			* original approach: arch effects ignored, garch applied to the mean of last period's sigma2
			* (for t = 1, we don't know prior garch value so it's 0)
			local lag "garch1*`pastsigma2'"
		}
		qui replace `scur' = exp(__cons `zg') + `lag'

		* percentile CIs, once we're through burnin
		if `t' > `burnin' {
			_pctile `scur', p(`plist')
			forv j = 1/`ncol' {
				matrix `Q'[`t' - `burnin', `j'] = r(r`j')
			}
		}

		qui replace `slag' = `scur'
		qui su `scur', meanonly
		local pastsigma2 = r(mean)
	}

	* [v2.1] missing percentiles mean the paths overflowed (explosive draws)
	mata: st_local("nmissQ", strofreal(sum(st_matrix("`Q'") :>= .)))
	if `nmissQ' > 0 {
		di in g _n "Warning: `nmissQ' band value(s) are missing because simulated variances were explosive draws; those periods will not be plotted."
	}

	* matrix -> dataset
	matrix colnames `Q' = `cnames'
	qui clear
	qui svmat double `Q', names(col)
	qui gen time = _n

	if "`fmode'" == "resid" mata: mata drop _gbplot_pool
	local ytitle = cond(inlist("`fmode'", "arch", "resid"), "Conditional Error Variance", "E[Conditional Error Variance]")
	local xtitle "time"
	if "`fmode'" == "epu" local note "`note' (parameter + expected path uncertainty)"
	if "`fmode'" == "arch" local note "`note' (parameter + path + ARCH innovation uncertainty, N(0,1) draws)"
	if "`fmode'" == "resid" local note "`note' (parameter + path + ARCH innovation uncertainty, resampled residuals)"
}

* ------------ Common: PSV, plot, save -----------------
if "`psv'" != "" { // is pre-shock value specified? If so rescale as % of preshock value
	if "`usematrix'" != "" { // if usematrix, preshock time is either specified or defaults to 1
		if "`psv_usematrix'" != "" {
			* [Fix D2/D3] after dropping burnin, time=1..(totalsimtime-burnin), so shift the raw matrix-row index
			local preshocktime = `psv_usematrix' - `burnin'
			noi di "Note: 'psv_usematrix' specified; using time period `preshocktime' as pre-shock value (adjusted for burnin=`burnin')."
		}
		else {
			noi di "Note: 'psv_usematrix' option not specified; using time period 1 as pre-shock value."
			local preshocktime = 1
		}
	}
	else {
		local preshocktime = `shocktime' - 1
	}
	* [FIX B6] Guard against preshocktime < 1
	if `preshocktime' < 1 {
		di in r _n "Error: 'shocktime' must be > 1 when using 'psv' (pre-shock reference period would be t=`preshocktime')."
		exit 198
	}
	su `cvar' if time == `preshocktime', meanonly
	* [OPT] Guard against zero/missing pre-shock value, since would divide by 0
	if missing(r(mean)) | r(mean) == 0 {
		di in r _n "Error: pre-shock expected variance is zero or missing at t=`preshocktime'; cannot compute PSV."
		exit 198
	}
	local preshock_val = r(mean)
	foreach var in `cvar' `bands' {
		qui replace `var' = (`var'*100)/`preshock_val'
	}
	local psvtxt "(% of pre-shock value)"
}

* Plot unless requested not to
if "`noplot'" != "" {
	di in g _n "Note: 'noplot' option specified; estimates calculated but graph suppressed"
}
else {
	_gbplot_draw, cvar(`cvar') bands(`bands') type(`plottype') ytitle(`ytitle') xtitle(`xtitle') note(`note') psvtxt(`psvtxt')
}

* save if requested (file check already done above)
if "`saving'" != "" {
	noi save "`saving'.dta", replace
}

restore
end
* ----------------------------------


* ---- helper: validate the usematrix matrix (no missing values, right number of columns) ----
cap prog drop _gbplot_chkmat
prog define _gbplot_chkmat
args mat ncols what
mata: st_local("hasmiss", strofreal(any(st_matrix("`mat'") :>= .)))
if `hasmiss' {
	di in r _n "Error: usematrix() contains missing values."
	exit 198
}
if colsof(`mat') != `ncols' {
	di in r _n "Error: usematrix() must have `ncols' columns (one per `what'); found `=colsof(`mat')'."
	exit 198
}
end

* ---- helper: draw spike (default) or area plot; bands = lo1 hi1 [lo2 hi2 lo3 hi3], widest first ----
cap prog drop _gbplot_draw
prog define _gbplot_draw
syntax, cvar(string) bands(string) type(string) ytitle(string) xtitle(string) note(string) [psvtxt(string)]
if "`psvtxt'" != "" local ytitle `""`ytitle'" "`psvtxt'""' // two-line title so it doesn't get clipped
else local ytitle `""`ytitle'""'
local nb: word count `bands'
local colors = cond(`nb' == 2, "eltblue", "eltblue ebblue edkblue")
local plots ""
forv j = 1/`=`nb'/2' {
	local lo: word `=2*`j'-1' of `bands'
	local hi: word `=2*`j'' of `bands'
	local col: word `j' of `colors'
	if "`type'" == "spike" local plots `plots' (rspike `lo' `hi' time, lwidth(medthick) color(`col'))
	else local plots `plots' (rarea `lo' `hi' time, color(`col'))
}
if "`type'" == "spike" local plots `plots' (scatter `cvar' time, color(black))
else local plots `plots' (line `cvar' time, lwidth(medthick) lcolor(black))
twoway `plots', legend(off) ytitle(`ytitle') xtitle("`xtitle'") title("") note("`note'") xscale(range(1 .))
end
