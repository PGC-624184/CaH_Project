Here's a collection of model atmospheres for use in the different projects.
They are available for
The Sun:
	ap00t5750g45k0odfnew.dat
	t5750_g+4.5_z+0.00.mod
	lte058-4.5-0.0a+0.0.BT-Settl.6
Canopus - Mike has used examples with 6500 K but you may also want to consider the 7500 K solutions:
	ap00t6500g15k0odfnew.dat
	ap00t7500g15k0odfnew.dat
	lte065-1.5-0.0a+0.0.BT-Settl.6
	lte076-1.5-0.0a+0.0.BT-Settl.6
	t6500_g+1.5_z+0.00.mod
	t7500_g+1.5_z-0.25.mod
Proxima Cen:
	lte027-5.5-0.0a+0.0.BT-Settl.6
	t2700_g+5.5_z+0.00.mod
A0 stars:
	ap00t10000g40k0odfnew.dat
	lte100-4.0-0.0a+0.0.BT-Settl.6
B0 stars:
	ap00t30000g40k0odfnew.dat
	lte300-4.0-0.0a+0.0.BT-Settl.6


Three different sources were used. All files are plain text files that you can read in any text editor. Your OS may try to open them in Quicktime or similar. Formats are described briefly below.

Castelli & Kurucz https://wwwuser.oats.inaf.it/fiorella.castelli/grids.html - https://ui.adsabs.harvard.edu/abs/2003IAUS..210P.A20C/abstract
	ap00t10000g40k0odfnew.dat
	ap00t30000g40k0odfnew.dat
	ap00t5750g45k0odfnew.dat
	ap00t6500g15k0odfnew.dat
	ap00t7500g15k0odfnew.dat
Filenames give stellar parameters such that e.g. the last file has [M/H]=0.0, Teff=7500 K, log(g)=1.5.
The input parameters are listed without any comments, but you don't need them. 
Abundances are listed for the first 99 elemements: N_X/N_tot for H and He, and then log(N_X/N_tot) for all other elements.
The thermodynamic variables use cgs units: RHOX is the column mass i.e. the integral of rho*dx from top of atmosphere to current point [g / cm^2], P is the gas pressure [dyne / cm^2], XNE the number density of free electrons [cm^-3], ABROSS is the Rosseland opacity per mass [cm^2 / g], ACCRAD is the radiative acceleration [cm / s^2], and the others are typically just zeros. The final column, VELSND is the local sound speed [cm / s].



PHOENIX / BT-Settl: https://phoenix.ens-lyon.fr/Grids/BT-Settl/AGSS2009/STRUCTURES/ - https://ui.adsabs.harvard.edu/abs/2014IAUS..299..271A
	lte027-5.5-0.0a+0.0.BT-Settl.6
	lte058-4.5-0.0a+0.0.BT-Settl.6
	lte065-1.5-0.0a+0.0.BT-Settl.6
	lte076-1.5-0.0a+0.0.BT-Settl.6
	lte100-4.0-0.0a+0.0.BT-Settl.6
	lte300-4.0-0.0a+0.0.BT-Settl.6
Filenames give stellar parameters such that the final file has Teff=7600 K, log(g)=1.5, [Fe/H]=0.0 and [alpha/Fe]=0.0.
Files are EXTREMELY verbose, and are actually the log output from the program run. They contain the last few iterations towards converging the atmosphere, so atmospheric quantities are printed several times. Make sure you grab the LAST occurrence. Units are cgs.
Basic thermodynamic properties: search for "Pressure and temperature structure". tstd or tau is the reference optical depth. pgas (or pg) and pe are gas and free electron pressure [dyn / cm^2]. density (or rho) [g / cm^3]. mu is the mean molecular weight. radius and extension in cm. kapq is the Rosseland mean opaciy [cm^2 / g]


MARCS: https://marcs.astro.uu.se/ - https://ui.adsabs.harvard.edu/abs/2008A%26A...486..951G/abstract
	t2700_g+5.5_z+0.00.mod
	t5750_g+4.5_z+0.00.mod
	t6500_g+1.5_z+0.00.mod
	t7500_g+1.5_z-0.25.mod
Filenames give stellar parameters such that the final file has Teff=7500 K, log(g)=1.5, [Fe/H]=-0.25.
Files contain several tables stacked one after another. Here's a Fortran program that can read them, for inspiration: https://marcs.astro.uu.se/documents/auxiliary/readmarcs.f
Rather than code a complicated reader, consider cut-and-pasting the files to keep just the parts that you need.
First come the basic properties of the model, with physical units for the input parameters. 
Abundances are listed for the first 92 elements, on a log(N_X/N_H)+12 scale.
Thermodynamic quantities use cgs units: lgTauR and lgTau5 are log(Rosseland optical depth) and log(tau_500nm). Note that Depth [cm] increases inward. Pe, Pg, Prad and Pturb are the electron, gas, radiation and turbulent pressures; all use [dyne / cm^2].
KappaRoss is the Rosseland opacity [cm^2 / g]. Density [g / cm^3]. Mu is the mean molecular weight [amu]. Vconv [cm/s] and Fconv/F are the convective velocity and flux. RHOX is the column mass [g / cm^2].
Partial pressures are again in cgs units [dyn / cm^2]. 




