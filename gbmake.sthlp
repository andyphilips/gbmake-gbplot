{smcl}
{* *! version 2.1  22oct2026}{...}
{findalias asfradohelp}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "[R] help" "help help"}{...}
{viewerjumpto "Syntax" "gbmake##syntax"}{...}
{viewerjumpto "Description" "gbmake##description"}{...}
{viewerjumpto "Options" "gbmake##options"}{...}
{viewerjumpto "Authors" "gbmake##authors"}{...}
{viewerjumpto "References" "gbmake##references"}
{viewerjumpto "Citations" "gbmake##citations"}
{viewerjumpto "Examples" "gbmake##examples"}{...}
{viewerjumpto "Version" "gbmake##version"}{...}

{title:Title}

{phang}
{bf:gbmake} {hline 2} Generate bootstrap samples of the conditional variance equation from a previously estimated GARCH model

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:gbmake} 
{cmd:,} [{cmdab:b:type(string)} {cmdab:r:eps(#)} {cmdab:s:aving(string)} {cmd:replace} {cmd:coldstart} {cmd:fcast_arch}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opth b:type(string)}}bootstrap type{p_end}
{synopt:{opth r:eps(#)}}number of bootstraps{p_end}
{synopt:{opth s:aving(string)}}custom filename for bootstrapped parameter dataset{p_end}
{synopt:{opt replace}}overwrite Stata dataset when saving{p_end}
{synopt:{opt coldstart}}let arch choose its own starting values each residual bootstrap replication{p_end}
{synopt:{opt fcast_arch}}save whitened innovations from residual bootstrap for forecasting{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}

{marker description}{...}
{title:Description}

{pstd}
{cmd:gbmake} implements two bootstrapping techniques---parametric and residual-based---and saves bootstrapped values for future use in {cmd:gbplot}. The full approach is outlined in Benton, Jordan and Philips (2026).

{pstd}
Users can find and download the most up-to-date version of {cmd:gbmake} at:
 http://andyphilips.github.io/FOO/.

{pstd} 
{cmd:gbmake} first requires estimating a GARCH(1,0) or GARCH(1,1) model using {cmd:arch}, with a list of Z regressors in the heteroskedastic variance equation (i.e., {cmd:het()} is specified). By default, the program will then simulate b bootstrap values of each parameter in the heteroskedastic variance equation using a parametric bootstrap. Values are saved in a Stata dataset for future use in {cmd:gbplot}. Instead of the default parametric bootstrap, users can also perform a residual-based bootstrap (see Benton, Jordan and Philips 2026), which is noticeably slower but may increase precision of predictions. For both, the program warns (but does not remove) bootstrap simulations where ARCH/GARCH terms fall outside logical bounds (e.g., < 0).

{pstd}
If the residual-based bootstrap is desired, {cmd:gbmake} output will display information on the number of simulations attempted, which will always be greater than the number of simulations requested (due to challenges of GARCH models converging). (These failures of convergence are denoted by a red X.) This does not have an impact on the observed results. {cmd:gbmake} can also save whitened innovations from the residual-based bootstrap for future use in forecasting using {cmd:gbplot}. 

{marker options}{...}
{title:Options}
 
{phang}
{opth b:type(string)} is either {cmd:parametric} (also the default) for the parametric bootstrap, or {cmd:residual} for the residual-based bootstrap.

{phang}
{opth r:eps(#)} is the number of b bootstraps taken. By default, b=100 for the residual bootstrap and b=1000 for the parametric bootstrap. Since not all residual bootstraps may converge, the program will attempt to create b successful bootstraps (it will try up to reps x 3 iterations). Higher values should increase precision but increase computational time.

{phang}
{opth s:aving(string)} specifies the filename to save the bootstrapped parameter values for future use in {cmd:gbplot}. By default this dataset is named "gbmake_results.dta".

{phang}
{opt replace} Overwrite existing dataset when using {cmd:saving()}, if it exists.

{phang}
{opt coldstart} To speed up residual-based bootstrap simulations, {cmd:gbmake} uses original estimates from {cmd:arch} as starting values during estimation; specifying this option reverts to using whatever starting values {cmd:arch} uses. This option cannot be used for the {cmd:parametric} bootstrap.

{phang}
{opt fcast_arch} saves the whitened innovations used in the residual-based bootstrap for future use in forecasting using {cmd:gbplot}. These are called "gb_innov" in the saved dataset. This option cannot be used for the {cmd:parametric} bootstrap.

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

If you use {cmd:gbmake}, please cite us:

{phang}Benton, Allyson L, Soren Jordan and Andrew Q Philips. 2026. “Taking Variance Seriously: Visualizing the Statistical and Substantive Significance of ARCH-GARCH Models.” The Journal of Politics 88(1):18–35.{p_end}

and

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

{phang}Create b=2000 bootstraps using the parametric bootstrap{p_end}
{phang}{cmd:gbmake, btype(parametric) reps(2000)}{p_end}

{phang}Create b=20 bootstraps using the residual-based bootstrap, changing the name of the saved dataset, and saving the innovations for future forecasting{p_end}
{phang}{cmd:gbmake, btype(residual) reps(20) saving("gbmake_resid.dta") replace fcast_arch}{p_end}

{marker version}{...}
{title:Version}

version 2.1, Oct 22, 2026.