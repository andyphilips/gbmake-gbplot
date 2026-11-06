{smcl}
{* *! version 2.1  22oct2026}{...}
{findalias asfradohelp}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "[R] help" "help help"}{...}
{viewerjumpto "Syntax" "gbplot##syntax"}{...}
{viewerjumpto "Description" "gbplot##description"}{...}
{viewerjumpto "Options" "gbplot##options"}{...}
{viewerjumpto "Authors" "gbplot##authors"}{...}
{viewerjumpto "References" "gbplot##references"}
{viewerjumpto "Citations" "gbplot##citations"}
{viewerjumpto "Examples" "gbplot##examples"}{...}
{viewerjumpto "Version" "gbplot##version"}{...}

{title:Title}

{phang}
{bf:gbplot} {hline 2} Create dynamic simulations of ARCH-GARCH models

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:gbplot} 
{cmd:[if] [in],} [{cmdab:shockvar:iable(string)} {cmdab:shockval:ue(#)} {cmdab:delta} {cmdab:using(string)} {cmdab:burnin(#)} {cmdab:usematrix(string)} {cmdab:shocktime(#)} {cmdab:range(#)} {cmdab:level(#)} {cmdab:customvalues(numlist)} {cmdab:psv} {cmdab:sav:ing(string)} {cmdab:noplot} {cmdab:psv_usematrix(#)} {cmdab:areaplot} {cmdab:replace} {cmdab:perm:anent} {cmdab:fcast_epu} {cmdab:fcast_arch} {cmdab:fcast_arch_resid}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opth shockvar:iable(string)}}Z variable to impulse/shock{p_end}
{synopt:{opth shockval:ue(numlist)}}size of shock/impulse{p_end}
{synopt:{opt delta}}use delta method instead of bootstraps{p_end}
{synopt:{opth using(string)}}dataset containing bootstraps from {cmd:gbmake}{p_end}
{synopt:{opth burnin(#)}}initial burn-in period{p_end}
{synopt:{opth usematrix(string)}}Stata matrix of shocks/impulses for all Z{p_end}
{synopt:{opth shocktime(#)}}simulation time for impulse/shock{p_end}
{synopt:{opth range(#)}}total simulation range{p_end}
{synopt:{opth level(#)}}significance level{p_end}
{synopt:{opth customvalues(numlist)}}set custom values other than means for Z{p_end}
{synopt:{opt psv}}express predictions as % of pre-shock value{p_end}
{synopt:{opth sav:ing(string)}}save predictions as Stata dataset{p_end}
{synopt:{opt noplot}}suppress plot creation{p_end}
{synopt:{opth psv_usematrix(#)}}time to use for {cmd:psv} calculation if {cmd:usematrix} specified{p_end}
{synopt:{opt areaplot}}create areaplot rather than spikeplot{p_end}
{synopt:{opt replace}}overwrite Stata dataset when saving{p_end}
{synopt:{opt permanent}}specify permanent rather than default impulse shock{p_end}
{synopt:{opt fcast_epu}}forecast ({cmd:using()} only), expected path uncertainty{p_end}
{synopt:{opt fcast_arch}}forecast ({cmd:using()} only), expected path uncertainty with N(0,1) ARCH innovations{p_end}
{synopt:{opt fcast_arch_resid}}forecast ({cmd:using()} from residual bootstrap only), expected path uncertainty with sample-specific ARCH innovations{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}

{marker description}{...}
{title:Description}

{pstd}
{cmd:gbplot} uses saved bootstraps from {cmd:gbmake} to create dynamic impulse-response style plots from Z regressors in the heteroskedastic variance equation from ARCH-GARCH models. The full approach is outlined in Benton, Jordan and Philips (2026).

{pstd}
Users can find and download the most up-to-date version of {cmd:gbplot} at:
 http://andyphilips.github.io/FOO/.

{pstd} 
{cmd:gbplot} can be used in two main ways (Approach 1 and Approach 2): 

{pstd}
In Approach 1, the user specifies the saved bootstrap dataset using {cmd:using()} (which was created using {cmd:gbmake} and contains a (b x k) dataset of b bootstraps for k regressors, Z, from the heteroskedastic variance equation). When using {cmd:using()}, the program warns (but does not remove) bootstrap simulations where ARCH/GARCH terms fall outside logical bounds (e.g., < 0). Alternatively, users can first estimate a model using {cmd:arch} and then use option {cmd:delta} to use the delta method approximation for confidence intervals in lieu of bootstraps. 

{pstd}
Approach 1 requires specifying the following: {cmd:shockvariable()}, the variable z observe an impulse/shock of size {cmd:shockvalue()} at time {cmd:shocktime()} out of total range time {cmd:range()}. By default, the program sets all variables Z to their means throughout the period (only {cmd:shockvariable()} experiences a change at {cmd:shocktime()}), although {cmd:customvalue()} allows users to specify other values instead of means. Optionally, {cmd:burnin()} adds additional time at the start of the simulation before {cmd:range()} that might create more stable predictions. Default shocks are impulses, although permanent shocks can be specified using {cmd:permanent}.

{pstd}
Approach 2 is much more customizable as it allows users to set Z to any values throughout the period of simulation, and allows for multiple, complex shocks to occur. Either a bootstrapped dataset or the delta method must be provided as in Approach 1. In addition, a Stata matrix of size (m x k) must be provided, where each column k shows the desired value for regressor z_k at time (row) m. Options {cmd:shockvariable()}, {cmd:shockvalue()}, {cmd:shocktime()}, {cmd:range()}, {cmd:customvalue()} are ignored under Approach 2. Burn-in periods can be incorporated into the matrix by the user, therefore using the {cmd:burnin()} option for Approach 2 will not ignore the 1st through m'th row specified in {cmd:burnin()} when plotting.

{pstd}
{cmd:gbplot} supports three forms of forecasting. {cmd:fcast_epu} incorporates 'expected path uncertainty' (EPU) among the bootstrap replicants by using each bootstrap's own t-1  conditional variance when calculating t's conditional variance, instead of the default of averaging across bootstraps for form the expected conditional error variance at t-1. {cmd:fcast_arch} incorporates EPU as well as a N(0,1) innovation which, when squared and multiplied by the lagged conditional error variance, enters as the ARCH effect in each period. {cmd:fcast_arch_resid} does the same but uses whitened sample-specific innovations from the residual-based bootstrap. For all three approaches, median expected values are shown along with percentile confidence intervals. {cmd:fcast_epu} and {cmd:fcast_arch} only work if {cmd:using()} is specified, and {cmd:fcast_arch_resid} only if {cmd:fcast_arch} was specified for the residual bootstrap when using {cmd:gbmake} (see that help file for details).

{pstd}
Additional options for both approaches include {cmd:saving()}, which saves the simulations as a Stata dataset (e.g., for customizable plots), and {cmd:psv}, which expresses predictions not in expected variance (default) but as a percentage of the pre-shock expected variance. If using {cmd:psv} for Approach 2, it is advisable to also specify {cmd:psv_usematrix()} which is used as the baseline time point to calculate all future percentage changes from. {cmd:noplot} suppresses plot creation and {cmd:areaplot} creates areaplots rather than the default spikeplot.

{marker options}{...}
{title:Options}
 
{phang}
{opth shockvar:iable(string)} is the z variable in Z to impulse/shock using Approach 1. Only one z may be specified.

{phang}
{opth shockval:ue(#)} is the size of the shock/impulse to {cmd:shockvariable}; thus the level of z will change to (mean + {cmd:shockvalue}) at time {cmd:shocktime}. This option is required when using Approach 1. If not specified this defaults to 1.

{phang}
{opth shocktime(#)} is the simulation time for the impulse/shock to occur. If not specified this defaults to t=5. This option is for Approach 1.

{phang}
{opth range(#)} specifies the total simulation range to create and plot. If not specified this defaults to t=15. This option is for Approach 1.

{phang}
{opth customvalues(numlist)} allows for custom values other than the default means for Z. One value must be provided for each Z, respectively (i.e., in the variable order specified in {cmd:arch}). This option is for Approach 1.

{phang}
{opt delta} uses the delta method for approximating confidence intervals rather than bootstraps. This option is required unless using {cmd:using} for both Approach 1 and 2, and requires first estimating a ARCH-GARCH model using {cmd:arch}.

{phang}
{opth using(string)} dataset containing bootstraps from {cmd:gbmake}. This option is required unless using {cmd:delta} for both Approach 1 and 2.

{phang}
{opth burnin(#)} optional burn-in period before starting the simulation time observed in the plot. Specifying burn-in is optional but is highly suggested as it will likely produce more stable simulations. Under Approach 2 this will discard the first row up to the row specified when plotting.

{phang}
{opth usematrix(string)} Stata matrix of shocks/impulses for all Z. This option is required when using Approach 2. Users might especially favor this option if the het equation uses many mutually exclusive dummy variables (i.e. the means are uninformative OR cannot be simultaneously met).

{phang}
{opth level(numlist)} optionally set significance level (in form of 1-alpha).

{phang}
{opt psv} express predictions as a percentage relative to the the expected pre-shock value of variance rather than the default, which shows the expected value of variance. For Approach 1 {cmd:psv} will use the time period just before {cmd:shocktime()} as the reference period. For Approach 2 the default is the first time period of {cmd:usematrix()}, so it is highly suggested to change this to the user's desired reference period using {cmd:psv_usematrix()}.

{phang}
{opth sav:ing(string)} Stata dataset name to save predictions as for future use, possibly for plotting in other software or custom plots with custom colors.

{phang}
{opt replace} Overwrite existing dataset when using {cmd:saving()}, if it exists.

{phang}
{opt noplot} suppresses plot creation (e.g., if preferring to save predictions and make a future custom plot).

{phang}
{opth psv_usematrix(#)} time to use for {cmd:psv} calculation if {cmd:usematrix} specified. This option is for Approach 2.

{phang}
{opt areaplot} plots an areaplot rather than the default spikeplot.

{phang}
{opt permanent} requests a permanent rather than temporary (the default) shock when using Approach 1.

{phang}
{opt fcast_epu} incorporates 'expected path uncertainty' (EPU) among the bootstrap replicants by using each bootstrap's own t-1  conditional variance when calculating t's conditional variance, instead of the default of averaging across bootstraps for form the expected conditional error variance at t-1. {cmd:using()} must be specified.

{phang}
{opt fcast_arch} incorporates {cmd:fcast_epu} as well as a N(0,1) innovation which, when squared and multiplied by the lagged conditional error variance, enters as the ARCH effect in each period. {cmd:using()} must be specified.

{phang}
{opt fcast_arch_resid} is the same as {cmd:fcast_arch} but uses whitened sample-specific innovations from the residual-based bootstrap. {cmd:using()} must be specified, and must have been created by  {cmd:gbmake} with option {cmd:fcast_arch} specified (see that help file for details).

{marker authors}{...}
{title:Authors}

{pstd}
Allyson L. Benton {break}
Department of Government {break}
University of Essex {break}
allyson.benton@essex.ac.uk {p_end}

{pstd}
Soren Jordan {break}
Department of Political Science {break}
Texas A&M University {break}
sorenjordanpols@gmail.com {p_end}

{pstd}
Andrew Q. Philips {break}
Department of Political Science {break}
University of Colorado Boulder {break}
andrew.philips@colorado.edu {p_end}

{marker references}{...}
{title:References}

If you use {cmd:gbplot}, please cite us:

{phang}Benton, Allyson L, Soren Jordan and Andrew Q Philips. 2026. “Taking Variance Seriously: Visualizing the Statistical and Substantive Significance of ARCH-GARCH Models.” The Journal of Politics 88(1):18–35.{p_end}

and:

{phang}Benton, Allyson L, Soren Jordan and Andrew Q Philips. 2026. “Dynamic inferences from models for conditional autoregressive heteroskedasticity.” Working Paper.{p_end}

{marker citations}{...}
{title:Citations}

{phang}Benton, Allyson L, Soren Jordan and Andrew Q Philips. 2026. “Dynamic inferences from models for conditional autoregressive heteroskedasticity.” Working Paper.{p_end}

{phang}Benton, Allyson L, Soren Jordan and Andrew Q Philips. 2026. “Taking Variance Seriously: Visualizing the Statistical and Substantive Significance of ARCH-GARCH Models.” The Journal of Politics 88(1):18–35.{p_end}

{marker examples}{...}
{title:Examples}

{phang}Open fictional Acme/Anvil stock price data and estimate GARCH(1,1) with heteroskedastic variance regressors Z:{p_end}
{phang}{cmd:use https://www.stata-press.com/data/r18/acmeh}{p_end}
{phang}{cmd:arch acme l.acme afrelated afinputs apex, arch(1) garch(1) het(afrelated apex)}{p_end}

{phang}Create b=2000 bootstraps using the parametric bootstrap (see gbmake for details){p_end}
{phang}{cmd:gbmake, btype(parametric) reps(2000)}{p_end}

{phang}Create an impulse of +1 increase to afrelated:{p_end}
{phang}{cmd:gbplot, shockvar(afrelated) shockval(1) using("gbmake_results.dta")}{p_end}

{phang}Same impulse but increasing the range of the simulation, the burn-in and holding afrelated at -1 pre-shock and apex at 2{p_end}
{phang}{cmd:gbplot, shockvar(afrelated) shockval(1) using("gbmake_results.dta") range(25) burnin(50) customvalues(-1 2)}{p_end}

{phang}Change to the delta method and expressing expected values as % of pre-shock value:{p_end}
{phang}{cmd:gbplot, shockvar(afrelated) shockval(1) delta range(25) burnin(50) customvalues(-1 2) psv}{p_end}

{phang}Add expected path uncertainty and ARCH effects (using N(0,1) innovations), and request an areaplot:{p_end}
{phang}{cmd:gbplot, shockvar(afrelated) shockval(1) using("gbmake_results.dta") fcast_arch areaplot}{p_end}

{marker version}{...}
{title:Version}

version 2.1, October 22, 2026.