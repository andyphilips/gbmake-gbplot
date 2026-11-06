*
*		PROGRAM GBMAKE
*
*		version 2.1
*		date: 10/20/26
*		authors: Allyson Benton, Soren Jordan and Andrew Philips
*
*
*
* Changes from v2.0:
*	[v2.1 N1] added fcast_arch option (for btype(residual) only). Saves the whitened residuals used by
*			  the bootstrap as variable gb_innov in the results file, for use with gbplot's 
*	          ARCH innovations. Rows 1..N_resid are filled; if N_resid > number of saved reps, file
*		 	  padded with missing parameter rows.
*
* Changes from v1.4:
*	[v2 S1]  speed: residual bootstrap defaults to warm-start each iteration using original estimated
*	          coefficients (from(b, copy)). should 1) make bs's closer 2) lead
*	 		  to faster convergence. `coldstart` new option to turn this off. (also suggest now that
*			  users limit iterate(#) in og arch call to increase speed, since gbmake uses those options).
*	[v2 S2]  speed: bootstrapped epsilons now drawn from a Mata call from a pre-cleaned
*	          pool
*
*	[v2 F1]  Parametric path names columns by position, not variable name. v1 crashed
*	          (corr2data dimension mismatch) with L2.y or a ts-operator such as L.x in
*	          the mean equation (allowed in parametric only), so might've misnamed cols if
*			  L.y not typed first.
*	[v2 F2]  moved replace data already exists to before bs loop to avoid grumpy user. allow spaces
*		     in file paths too.
*	[v2 F3]  Vars created inside residual bootstrap now tempvars. v1 used varnames (so potential conflict
*			 if they already existed)  
*	[v2 F4]  Similarly, matrices now tempnames (B and VCV)
*	[v2 F5]  Models with noconstant in mean eq will no longer produce a syntax error.  
*	[v2 F6]  Ctrl-Break in residual loop now closes postfile, restores og arch estimates, and exits
*			 with rc=1 (v1 didn't do first 2, and returned rc=602).
*	[v2 F7]  seems unlikely, but resid bs w/ zero converged runs exits with an error instead of
*	          saving empty dataset; and simulated series now stored in double precision.
*	[v2 F8]  if/in removed from syntax (previously was ignored anyway).
*
* Changes from v1.3:
*	[fix D1] while-loop redraw: changed & to | for y_t0/sigma2_t0 (was only
*	          redrawing when BOTH missing; should redraw if EITHER is missing)
*	[Fix D2] Convert e(tmin)/e(tmax) time-variable values to row numbers before
*	          using `in` qualifier; e(tmin)/e(tmax) are time values (e.g. 1990),
*	          not row indices — `in 1990` silently did nothing for most real data
*	[Fix D3] Guard garch1 bounds check in parametric path: skip if Bvar_garch=""
*	          (ARCH-only models have no garch1 variable; was crashing)
*	[Fix D4] Move residual-path garch1 bounds check to after `use coefs_temp`:
*	          the check was running on the original preserved dataset (no garch1),
*	          crashing every residual run; also added Bvar_garch guard
*
* Changes from v1.2:
*	[Fix C1] .dta was appended always, even if user specified .dta
*	[Fix C2] quietly examine GARCH terms and throw warning if any outside logical bounds
*
* Changes from v1.1:
*	[Fix B1] exit if resid and gaps exist
*
* Changes from v1.0:
*   [FIX A1] Added 'replace' to syntax line (was missing; option was non-functional)
*   [FIX A2] Added total_attempts counter in residual bootstrap while-loop
*             to cap iterations at minreps_try total attempts (not successes)
*   [FIX A3] Changed reps(string) to reps(numlist integer max=1 >0) for cleaner validation
*   [FIX A4] Added manual loop kill to resid bootstrap
*   [FIX A5] fixed the initial t=1 and t=2 issue whereby we weren't using the same (though
*			 lagged) error to generate both sigma2 and y in the resid bs
*   [OPT]    Clarified ARCH/GARCH t=1 formula comment and notation in error messages
*
* -----------------------------------------------------------------------------
* -----------------------------------------------------------------------------
* -----------------------------------------------------------------------------


cap prog drop gbmake
prog define gbmake, rclass
syntax, [Btype(string) Reps(numlist integer >0 max=1) Saving(string) replace coldstart fcast_arch]
* gbmake: Generate bootstrap samples from a previously estimated ARCH/GARCH model
* btype: parametric (default) or residual
* reps: number of bootstrap replications (default 1000 for parametric, 100 for residual)
* saving: filename to save bootstrap data (default gbmake_results.dta)
* replace: overwrite 'saving' if file already exists
* coldstart: residual bootstrap only; let arch choose its own starting values each fit
*	(default is to start at the original estimates, which is much faster)
* fcast_arch: residual bootstrap only; saves whitened residuals as gb_innov (for gbplot, fcast_arch_resid)

* ------------ Checks -----------------
* assert arch previously run
if "`e(cmd)'" != "arch"	{
	di in r _n "Error: you must first run arch before using gbmake"
	exit 198
}

* validate btype option
if "`btype'" == "" {
	loc btype "parametric"
}
else if !inlist("`btype'", "parametric", "residual") {
	di in r _n "Error: btype() must be either 'parametric' or 'residual'"
	exit 198
}

* ensure GARCH(1,0) or GARCH(1,1) only, and record for below
if "`e(garch)'" == "1" {
	local is_garch = 1
	if "`e(arch)'" == "1" {
		local is_arch = 1 // pass
	}
	else if "`e(arch)'" == "0" {
		local is_arch = 0 // pass
	}
	else {
		di in r _n "Error: only arch(1) or arch(0) currently supported"
		exit 198
	}
}
else {
	di in r _n "Error: only arch(1) garch(1) or arch(0) garch(1) currently supported"
	exit 198
}

* ensure no weird GARCH flavors:
if "`e(saarch)'" != "" | "`e(tarch)'" != "" | "`e(aarch)'" != "" | "`e(narch)'" != "" | ///
 "`e(narchk)'" != "" | "`e(abarch)'" != "" | "`e(atarch)'" != "" | "`e(sdgarch)'" != "" | ///
 "`e(earch)'" != "" | "`e(parch)'" != "" | "`e(tparch)'" != "" | "`e(aparch)'" != "" | ///
  "`e(nparch)'" != "" | "`e(nparchk)'" != "" | "`e(pgarch)'" != "" {
	  di in r _n "Error: only arch(1) garch(1) or arch(0) garch(1) currently supported"
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

* [v2.1 N1] fcast_arch only makes sense with the residual bootstrap
if "`fcast_arch'" != "" & "`btype'" != "residual" {
	di in r _n "Error: 'fcast_arch' requires btype(residual) (parametric bootstrap has no residuals to save)."
	exit 198
}

* strip .dta from saving if user added it, then set output file name:
if "`saving'" != "" {
	if lower(substr("`saving'", -4, .)) == ".dta" {
		local saving = substr("`saving'", 1, length("`saving'") - 4)
	}
}
else local saving "gbmake_results"

* [v2 F2] check output file doesn't already exist
cap confirm file "`saving'.dta"
if _rc == 0 & "`replace'" == "" {
	di in r _n "Error: File `saving'.dta already exists (use 'replace' to overwrite)"
	exit 602
}
* ----------------------------------


* ------------ Obtain important values -----------------
* set number of boots (parametric, default = 1000) / minimum successful boots (residual, default = 100)
* (reps() is already validated as a single positive integer by syntax)
if "`reps'" != "" {
	loc bs `reps'
}
else {
	if "`btype'" == "residual" {
		loc bs 100
	}
	else {
		loc bs 1000
	}
}
if "`btype'" == "residual" {
	loc minreps_try = `bs' * 3 // set target to 3x desired successful reps (= max total attempts)
}



* ---------- namepeeler: vectorize all relevant variables for use below --------
* go through column `index' of e(b); the equation name tells us where it lives:
	* if eqnames == DV name: [then it's in the mean EQ]
		* if bname == L.DV, then it's an LDV and should be stored accordingly
		* else, if bname != _cons, then its an X and should be stored accordingly
		* else, it's the constant
	* if eqnames == ARCHM, it's an archm term to store
	* if eqnames == HET:
		* store as het term. also assert there's no DV lurking here...too much to deal with
	* if eqnames == ARCH:
		* if bname == L.arch, store as arch term
		* if bname == L.garch, store as garch term

tempname B VCV b0 // [v2 F4] added tempnames
mat `B' = e(b)
loc bnames: colnames `B'
loc eqnames `e(eqnames)'
local old_dv "`e(depvar)'" // store for later use in residual bootstrap and to avoid calls to new depvar
local t_min = `e(tmin)' // store for later use in residual bootstrap
local t_max = `e(tmax)'

* assert no gaps, which cause an error in the resid bootstrap:
if `e(N_gaps)' != 0 & "`btype'" == "residual" {
	di in r _n "Error: Residual bootstrap cannot be used when gaps exist in the ARCH model"
	exit 198
}

* [v2 F5] linear predictors start at 0 so a missing constant can't leave a dangling "+"
loc meqlist_B "0"
loc hetlist_Z "0"
* [v2 F1] keepidx/keepnam: specify which columns of e(b) we keep for parametric bootstrap, and what to call them
loc keepidx ""
loc keepnam ""

local index 1 // loop through bnames (list of variable names)
foreach eq_name of local eqnames {

	* first, obtain the variable name:
	loc vname: word `index' of `bnames'

	* if eq_name == DV name: [then it's in the mean EQ]
	if "`eq_name'" == "`old_dv'" {
		* assert no higher-order LDVs if doing residual bootstrap
		if regexm("`vname'", "^L([2-9]|[1-9][0-9]+)\.`old_dv'$") {
			if "`btype'" == "residual" {
				di in r _n "Error: lagged dependent variables greater than t-1 only supported using parametric bootstrap"
		    	exit 198
			}
		}
		* if bname == L.DV, then it's an LDV and should be stored accordingly
		else if "`vname'" == "L.`old_dv'" | "`vname'" == "l.`old_dv'" | "`vname'" == "L1.`e(depvar)'" | "`vname'" == "l1.`old_dv'" {
			loc ldv_B "_b[`old_dv':`vname']" // store parameter for resid bs
			loc ldv_vname "`vname'"	// store variable name for coefficient extraction
		}
		else if "`vname'" != "_cons" { // all X's get stored accordingly
			loc meqlist_B "`meqlist_B' + _b[`old_dv':`vname']*`vname'" // store for resid bs
		}
		else { // its a constant
			loc meqlist_B "`meqlist_B' + _b[`old_dv':_cons]" // store for resid bs
		}
	}
	if "`eq_name'" == "ARCHM" { // if eqnames == ARCHM, it's an archm term to store.
		* We've already asserted that we don't have anything silly like higher order ARCHM
		loc archm_B "_b[ARCHM:sigma2]" // store for resid bs
	}

	if "`eq_name'" == "HET" { // if eqnames == HET, we're in variance eq
		if regexm("`vname'", "^[^.]*\.`old_dv'\$") { // we don't want any form of the DV in this equation
			di in r _n "Error: the variance eq. cannot contain the dependent variable"
			exit 198
		}
		else if "`vname'" != "_cons" {
			* [v2 F1] since these names become variable names in the saved dataset, no ts/factor operators allowed
			if strpos("`vname'", ".") {
				di in r _n "Error: het() variable `vname' uses a time-series or factor-variable operator; not currently supported"
				exit 198
			}
			loc Bvar_Z "`Bvar_Z' `vname'" // list for param bs
			loc hetlist_Z "`hetlist_Z' + _b[HET:`vname']*`vname'" // list for resid bs
			loc keepnam "`keepnam' `vname'"
		}
		else {
			loc Bvar_Z "`Bvar_Z' __cons" // stata won't allow word "_cons" in parametric bs?!
			loc hetlist_Z "`hetlist_Z' + _b[HET:_cons]" // list for resid bs
			loc keepnam "`keepnam' __cons"
		}
		loc keepidx "`keepidx' `index'"
	}

	if "`eq_name'" == "ARCH" {
		if "`vname'" == "L.arch" { // obtain any ARCH terms
			loc Bvar_arch "arch1"
			loc keepidx "`keepidx' `index'"
			loc keepnam "`keepnam' arch1"
		}
		if "`vname'" == "L.garch" { // obtain any GARCH terms
			loc Bvar_garch "garch1"
			loc keepidx "`keepidx' `index'"
			loc keepnam "`keepnam' garch1"
		}
	}

	local index = `index' + 1
}

* ----------------------------------



* ----------------------------------
* DO PARAMETRIC [also default]
if "`btype'" == "parametric" {
	preserve
	mat `VCV' = e(V)
	clear

	* [v2 F1] draw ALL parameters with generic positional names (__b1, __b2, ...), then
	* rename the ones we keep. Works for any mean equation (L2.y, L.x, ...).
	* when using parametric
	loc k = colsof(`B')
	forv i = 1/`k' {
		loc allb "`allb' __b`i'"
	}
	qui corr2data `allb', n(`bs') means(`B') cov(`VCV')

	* keep only relevant params; just the het, arch and garch
	loc n_keep: word count `keepidx'
	forv j = 1/`n_keep' {
		loc i: word `j' of `keepidx'
		loc nm: word `j' of `keepnam'
		rename __b`i' `nm'
	}
	keep `keepnam'
	qui compress

	* Fix C2: warn if GARCH outside logical bounds
	* Fix D3: added Bvar_garch wrapper just in case someone specifies ARCH-only and no GARCH
	if "`Bvar_garch'" != "" {
		qui count if garch1 <= 0 | garch1 >= 1
		if r(N) != 0 {
			noi di in g "WARNING: `r(N)' GARCH parameters outside (0,1) logical bounds; consider dropping illogical bootstraps or respecifying model"
		}
	}

	* print and save (file check already done above)
	noisily di in g "`bs' parametric bootstraps completed, saved as:"
	noi save "`saving'.dta", replace
	restore
}
* ----------------------------------

* ----------------------------------
* DO RESIDUAL
if "`btype'" == "residual" {
	* Store the original arch estimation results before bootstrap (estimates will be overwritten)
	quietly estimates store _gbmake_original_arch
	mat `b0' = e(b) // [v2 S1] starting values for every re-estimation

	preserve
	cap tsset // need to reset tsset for some strange reason
	if _rc {
		di in r _n "Error: data must be tsset for residual bootstrap"
		restore
		exit 198
	}

	* [Fix D2] Capture timevar name from r(timevar) immediately after tsset.
	* e(tmin)/e(tmax) store actual time-variable VALUES (e.g. calendar year 1990),
	* not row numbers.
	local tvar "`r(timevar)'"
	qui count if `tvar' < `t_min'
	local row_min = r(N) + 1
	qui count if `tvar' <= `t_max'
	local row_max = r(N)
	local row_plus1 = `row_min' + 1

	* [v2 F3] tempvars now to avoid existing names
	tempvar nu sigma2h epsilon sig2 ysynth epsb xb0 het0

	* predict vars that are going to be replaced:
	qui predict double `nu', resid // nu_t = \sigma_t\epsilon_t
	qui predict double `sigma2h', variance // sigma^2_t
	qui gen double `epsilon' = `nu'/sqrt(`sigma2h') // rescale, nu_t/\sigma_t = \epsilon_t
	su `epsilon', meanonly // mean-center
	qui replace `epsilon' = `epsilon'-r(mean)

	* generate quantities that get replaced every bootstrap
	qui gen double `sig2' = . // for storing sigma2.
	qui gen double `ysynth' = . // for storing synthetic y.
	qui gen double `epsb' = . // for storing bootstrapped epsilons.

	* need these values from the OG arch call
	loc cmd `"`e(cmdline)'"'

	* --- Capture original predicted Xb and heteroskedastic term using original estimates ---
	* evaluate `meqlist_B' and `hetlist_Z' now (before any re-estimation) and store as variables
	qui gen double `xb0' = `meqlist_B'
	qui gen double `het0' = exp(`hetlist_Z')

	* store numeric copies of coefficients that will be overwritten by re-estimation
	if "`Bvar_arch'" != "" {
		local archcoef = _b[ARCH:L.arch]
	}
	else local archcoef = 0
	local garchcoef = _b[ARCH:L.garch]
	if "`archm_B'" != "" {
		local archmcoef = _b[ARCHM:sigma2]
	}
	else local archmcoef = 0
	if "`ldv_B'" != "" {
		local ldvcoef = _b[`old_dv':`ldv_vname']
	}
	else local ldvcoef = 0

	* [v2 S2] build re-estimation call once: obtain arch call, replacing all instances
	* of DV name with synthetic DV (captures any LDVs)
	loc new_arch = ustrregexra(`"`cmd'"', "\b`old_dv'\b", "`ysynth'")
	* [v2 S1] warm start: add from(original b, copy) to options. Need to know whether the
	* call already has an options comma (ignore commas nested inside parentheses).
	if "`coldstart'" == "" & !strpos(lower(`"`cmd'"'), "from(") {
		loc tmp `"`new_arch'"'
		while ustrregexm(`"`tmp'"', "\([^()]*\)") {
			loc tmp = ustrregexra(`"`tmp'"', "\([^()]*\)", "")
		}
		if strpos(`"`tmp'"', ",") loc new_arch `"`new_arch' from(`b0', copy)"'
		else loc new_arch `"`new_arch', from(`b0', copy)"'
	}

	* [v2 S2] now using pool of non-missing std residuals, held in Mata for speed
	mata: _gbmake_eps = st_data(., "`epsilon'"); _gbmake_eps = select(_gbmake_eps, _gbmake_eps :< .)

	* set up postfile to store coefficients from each successful convergence
	cap postclose coeffs
	qui tempfile coefs_temp
	qui postfile coeffs `Bvar_Z' `Bvar_arch' `Bvar_garch' using `coefs_temp', replace

	* [v2 S2] expression to post is the same every time; build once
	foreach var of local Bvar_Z {
		if "`var'" == "__cons" local post_expr "`post_expr' (_b[HET:_cons])"
		else local post_expr "`post_expr' (_b[HET:`var'])"
	}
	if "`Bvar_arch'" != "" local post_expr "`post_expr' (_b[ARCH:L.arch])"
	if "`Bvar_garch'" != "" local post_expr "`post_expr' (_b[ARCH:L.garch])"

	* loop across bootstraps:
	local successful_sims = 0
	local total_attempts = 0
	noisily _dots 0

	while `successful_sims' < `bs' & `total_attempts' < `minreps_try' {
		local total_attempts = `total_attempts' + 1

		* Initialize locals for this iteration
		local resid_t0 = .
		local resid_t1 = .
		local y_t0 = .
		local sigma2_t0 = .
		loc fullmeanlist = ""
		loc add_arch_t1 = ""
		loc add_arch = ""

		* Initialize: draw t=0 and t=1 resids, for sigma2 and y, respectively
		while missing(`resid_t0') { // redraw if NA. Note these while calls only work for scalars.
			local resid_t0 = `epsilon'[runiformint(1, _N)] // t=0 error needed for sigma2
		}
		while missing(`resid_t1') {
			local resid_t1 = `epsilon'[runiformint(1, _N)] // t=1 error for y
		}

		* Initialize: draw y and sigma2 jointly.
		* These are t=0 values for y (if LDV = TRUE) and sigma2 (for GARCH), respectively
		* [Fix D1] Changed & to | so we redraw if EITHER initial value is missing
		while missing(`y_t0') | missing(`sigma2_t0') {
			local rand_index = runiformint(1, _N) // index random row
			local y_t0 = `old_dv'[`rand_index']
			local sigma2_t0 = `sigma2h'[`rand_index']
		}

		* [v2 S2] bootstrapped epsilons for every row, drawn from clean pool in one call
		* (so no missing values and no redraw loop)
		mata: st_store(., "`epsb'", _gbmake_eps[ceil(runiform(st_nobs(), 1) :* rows(_gbmake_eps))])

		* Do we have ARCH effects? Needed for DGP below (possibly in future add other arch flavors here).
		* (catches above rule out no GARCH effects)
		if "`Bvar_arch'" != "" {
			loc add_arch_t1 = "`archcoef'*`resid_t0'^2*`sigma2_t0' + " // arch effect for t=1 (requires arbitrary t=0 resid value)
		}

		* -- create t=1 values of sigma2 and y --
		* [Fix D2] Use row_min (row number) instead of t_min (time value) for `in`
		qui replace `sig2' = `add_arch_t1' `het0' + `garchcoef'*`sigma2_t0' in `row_min'
		if "`archm_B'" != "" {
			loc fullmeanlist "`fullmeanlist' + `archmcoef'*`sig2'"
		}
		if "`ldv_B'" != "" {
			loc fullmeanlist "`fullmeanlist' + `ldvcoef'*`y_t0'"
		}
		* use stored Xb from original estimates; full mean = original Xb + any LDV or archm terms
		qui replace `ysynth' = `xb0' `fullmeanlist' + `resid_t1'*sqrt(`sig2') in `row_min'

		* [FIX A5] Overwrite epsb at row_min with resid_t1.
		qui replace `epsb' = `resid_t1' in `row_min'

		* generate the rest of synthetic series from row_min+1 to row_max
		if "`Bvar_arch'" != "" {
			loc add_arch = "`archcoef'*l.`epsb'^2*l.`sig2' + " // ARCH for t=2/T
		}
		qui replace `sig2' = `add_arch' `het0' + `garchcoef'*l.`sig2' in `row_plus1'/`row_max'

		loc fullmeanlist ""
		if "`archm_B'" != "" {
			loc fullmeanlist "`fullmeanlist' + `archmcoef'*`sig2'"
		}
		if "`ldv_B'" != "" {
			loc fullmeanlist "`fullmeanlist' + `ldvcoef'*l.`ysynth'"
		}
		qui replace `ysynth' = `xb0' `fullmeanlist' + `epsb'*sqrt(`sig2') in `row_plus1'/`row_max'

		* Estimate the ARCH model, using the synthetic y as the new DV:
		cap qui `new_arch'
		local rc = _rc

		* [FIX A4] / [v2 F6] manual loop kill: clean up before leaving
		if `rc' == 1 {
			postclose coeffs
			mata: mata drop _gbmake_eps
			restore
			qui estimates restore _gbmake_original_arch
			qui estimates drop _gbmake_original_arch
			di in r _n "Error: User manually stopped bootstrap loop"
			exit 1
		}

		if !`rc' & "`e(converged)'" == "1" {
			local successful_sims = `successful_sims' + 1
			post coeffs `post_expr'
			noisily _dots 1 0
		}
		else {
			noisily _dots 1 1  // display 'x' for failure
		}

	}

	* close postfile, clear Mata pool, load coefficient results
	postclose coeffs
	noi di ""
	noi di in g "Achieved `successful_sims' successful simulations out of desired `bs'."

	* [v2 F7] nothing converged: stop here
	if `successful_sims' == 0 {
		mata: mata drop _gbmake_eps
		restore
		qui estimates restore _gbmake_original_arch
		qui estimates drop _gbmake_original_arch
		di in r _n "Error: no bootstrap replication converged; nothing saved. Try coldstart, or respecify model."
		exit 430
	}

	* [Fix D4] Load coefs_temp BEFORE running bounds check
	use `coefs_temp', clear

	* [v2.1 N1] attach whitened innovations (the pool the bootstrap resampled from), then free Mata
	if "`fcast_arch'" != "" {
		mata: st_local("npool", strofreal(rows(_gbmake_eps)))
		if `npool' > _N qui set obs `npool'
		qui gen double gb_innov = .
		mata: st_store((1::rows(_gbmake_eps)), "gb_innov", _gbmake_eps)
		noi di in g "Saved `npool' whitened innovations as gb_innov (use gbplot, fcast_arch_resid)."
	}
	mata: mata drop _gbmake_eps
	
	* warn if bad GARCH params (save 1 < ARCH < 0 or ARCH + GARCH >= 1 for gbplot)
	if "`Bvar_garch'" != "" {
		qui count if !missing(garch1) & (garch1 <= 0 | garch1 >= 1)
		if r(N) != 0 {
			noi di in g "WARNING: `r(N)' GARCH parameters outside (0,1) logical bounds; consider dropping illogical bootstraps or respecifying model"
		}
	}

	noisily di in g "`successful_sims' residual bootstrap simulations completed, saved as:"
	noi save "`saving'.dta", replace // (file check already done above so overwrite always)

	restore

	* Restore the original arch estimation results
	quietly estimates restore _gbmake_original_arch
	quietly estimates drop _gbmake_original_arch

	* Final summary
	noisily di in g "Residual bootstrap completed: `successful_sims' successful out of `total_attempts' total attempted simulations (max allowed attempts: `minreps_try')."
	if `successful_sims' < `bs' {
		noisily di in r "Warning: desired number of `bs' successful simulations not achieved."
	}
}
end
* ----------------------------------
