! ----------------------------------------------------------------
! Project: Complete One loop Predictions for B->mu+mu- in the MSSM
! Please cite: Dedes, Rosiek, Tanedo (2007)
! ----------------------------------------------------------------
! Main file for B->mu+mu- calculation
! filename: README.TXT
! Current version: 1.00 (SM)
! By Flip Tanedo (philip.tanedo@gmail.com)
! Date: 29 November 2007
! Latest version available at: n/a
! ----------------------------------------------------------------
! This is the README FILE
! ----------------------------------------------------------------
! RECORD OF REVISIONS
!	Date	 Programmer	Description
!	====	 ==========	==================================
!	18/11/07 Flip Tanedo	First version (running)
!	20/11/07 Flip Tanedo	Error messages, Troubleshooting
!	28/11/07 Flip Tanedo	Matching
!	29/11/07 Flip Tanedo	1.00 revision
! ----------------------------------------------------------------
! CONTENTS
! 0. Production-notes
! 1. About
! 2. Files
! 3. How to compile
! 4. Hacking the code
! 5. THINGS TO DO STILL
! 6. List of error messages
! 7. Troubleshooting
! 8. Important notes from debugging
! 9. Version notes
! ----------------------------------------------------------------

! ----------------------------------------------------------------
! 0. PRODUCTION NOTES (NOT FOR GENERAL PUBLIC)
! ----------------------------------------------------------------

1.	Current code is Standard-Model ONLY
2. 	ckm_init must be run in driver.f to initialise CKM matrix
3. 	I've turned off the switches in the diagrams since they're
	redundant in the SM-only calculations. These are 
	straightforward to include when I extend to SUSY.
	Write a module with the switches (ih,ig, etc... in Janusz's
	code).

! ----------------------------------------------------------------
! 1. ABOUT
! ----------------------------------------------------------------

This is Flip's version of the B->mu+mu- code.

The purpose is to provide an independent check of the programming
in the larger code.

If used externally, please cite Dedes, Rosiek, Tanedo (2007/8).
If you don't, I'll never cite any of your papers. Ever. o_o

! ----------------------------------------------------------------
! 2. FILES
! ----------------------------------------------------------------

Files listed below the single dashed line are non-essential.

FILENAME		DESCRIPTION
======================	==========================================
driver.f95		Main driver file (name subject to change)
program_parameters.f95	Contains non-physics parameters necessary
			for the program to run
physics_data.f95	Contains particle data
tuning_data.f95		Non-physical parameters (e.g. 't Hooft
			mass, divergence, `small difference')
integrals.f95		Passarino-Veltman Integrals
sorter.f95		Array sorting (for integrals.95)
errorlog.f95		Module for error-tracking
diagrams_box.f95	Box diagram contributions
diagrams_zpenguin.f95	Z-penguin contributions
diagrams_
	dselfenergy.f95 down quark self-energy contributions
amplitudes.f95		Calculates amplitude from diagrams_%%%.f95 files
branchingratio.f95	Calculates the branching ratio
----------------------  ------------------------------------------
Makefile		Makefile
README.TXT		This file.

! ----------------------------------------------------------------
! 3. HOW TO COMPILE AND RUN
! ----------------------------------------------------------------

Type in: "Make" (no quotations).
Type in: ./Bmumu

My Makefile is very rudimentary and needs updating.

! ----------------------------------------------------------------
! 4. HACKING THE CODE
! ----------------------------------------------------------------

BASIC STRUCTURE: 

	PROGRAM: bmumu
	--------------
	call	ckm_init (physics_data.f95)
		[troubleshooting write statements]
	func	br (branchingratio.f95)
	
	FUNCTION: br
	------------
	func	Fv,Fa,Fs,Fp
	func	amp2


PARAMETERS:
-----------

PROGRAM_PARAMETERS.F95 contains non-physics information relevant
to the running of the code. For example, it contains the definition
of the REAL KIND (KIND=DBL) for processor-independence. It also
contains the definition of pi and zero as a REAL(KIND=DBL).

PHYSICS_DATA.F95 contains all relevant data from the Particle Data Group,
as well as mixing and SUSY parameters and their initialization subroutines.
The current values are set to match the parameters in Janusz's code
(after RG evolution).

TUNING_DATA.F95 contains details about the 2/(d-4) divergence and 't Hooft
mass scale from dimensional regularisation. It also contains a definition
of `small difference' for the integral to distinguish between identical
arguments. The final result should be INDEPENDENT of these values.
NOTE: if you want to scan these parameters, turn off the `parameter' 
(ironic, isn't it?) tag in this file, and ALSO in the definitions in
integrals.f95. The definitions in this latter file are so that the
integral functions can use these values without having them listed
explicitly in their arguments. (This is one of the benefits of modules.) 

The non-obvious parameters here are 
	mu_tH, the 't Hooft scale (which should cancel)
	del, the pole from divergent integrals (should also cancel)
	eps, definition of a 'small' mass difference for integrals

STANDARD MODEL VS MSSM. Version 1 of this code only contains SM values.

CALCULATING THE AMPLITUDE. The amplitude is calculated from the
diagrams_%%%.f95 files/modules. Some of these modules contain complete
diagrams (contributions to the Wilson coefficients), such as 
the box diagrams. However, others only contain pieces of a diagram,
such as the penguins or self energies. These have to be pieced together
into complete diagrams within the amplitude module. I don't really like
this tangled-up structure, perhaps I'll make this more compartamentalised
in the future so that the amplitude module only sums complete diagrams.

BRANCHING RATIO. branchingratio.f95 contains the Br calculation.
If you're looking at this file for the first time, start at the last 
function and work your way up as each function calls a function
that was defined just above it.


GOING FROM FORTRAN 77 to FORTRAN 95
-----------------------------------
* Instead of complex*16 or double precision, use a processor-independent
  definition for data type. In program_parameters.f95 I define the
  integer DBL as the real number 'kind' that produces a specified
  precision for whatever processor is running the code. (i.e. the
  definition asks the processor to identify the correct integer that
  refers to a real number with the specified precision.) When 
  declaring variables, use e.g. COMPLEX(KIND=DBL).
* Instead of common blocks, use modules. Each module is given its
  own file with the same name. For example, all of my physics data
  are in physics_data.f95 in a module called physics_data. This is
  nice because it prevents name ambiguities. 
* Instead of $ at column 5 to indicate a spill-over line, use
  & at the end of the previous line. I also use a & at column 1 of the
  spill-over line, but this is purely optional. F95 has no line 
  length maximum, so you don't even have to spill over.
* Instead of c, use ! to indicate comments
* No WRITE statements in functions! (Recursive call error.)
* Instead of "implicit double precision (a-h,o-z)", use "IMPLICIT NONE"
  since this helps catch errors.


AN IMPORTANT NOTE ABOUT PORTING JANUSZ CODE: The main task is to 
change the F77 formatting, i.e. get rid of commons, implicits,
and include appropriate types (KIND=DBL). HOWEVER, I've made one
further change: Janusz's CKM matrix actually CKM^\dagger. My
CKM matrix is the true CKM matrix, so I've made the following
replacements when porting Jansuz's code:

	Janusz				This Code
	---------------			-----------------
	CKM(i,j)		--->	DCONJG(CKM(j,i))
	DCONJG(CKM(i,j))	--->	CKM(j,i)

! ----------------------------------------------------------------
! 5. THINGS TO DO STILL
! ----------------------------------------------------------------

1. Check values of physical data
	a. Check the values from PDG

2. I've put in Janusz's b-functions, but I don't know what all of this
	dbstat stuff is.

	Instead of going crazy over this (i think most of the options
	are irrelevant), I'll just put an error message if one of the
	unactivated integrals are called. (See section 6 below.)

3. Running masses
	Current inputs are at top mass scale. I want to fix this so that
	one inputs pole masses and the values are run to the top mass.
	Also make sure other parameters that are initialised 
	appropriately, e.g. Yukawas which are defined in terms of
	the quark masses.

4. Include SUSY terms... lots of baggage. This is the next big step.


5. This is easy and important: put in error checks for
	flavor integers greater than 3 and less than 1
	function: checkindices(i,j,k,l, `filename')

! ----------------------------------------------------------------
! 6. List of error messages
! ----------------------------------------------------------------

In addition to the usual FORTRAN 95 error calls, I have inserted my
own 'error' messages in the code. These messages do not halt the
code (though I suspect they lead to fatal errors), but just print out
tagged error notices. The dictionary of error tags is provided below.
UPDATE: No need for a dictionary anymore, the error log has the errors
listed explicitly.

The error handling module is a standard linked list set up.
It is located in errorlog.f95.



The following errors are old and should no longer be used.

integrals.f95
-------------
ERROR Int.dp0.1: two or more masses = 0 in dp0
ERROR Int.dp1.1: 3 or more masses=0 in dp1
ERROR Int.cp0.1: 2 or more masses = 0 in cp0
ERROR Int.cp1.1: 3 x 0 masses in cp1

ERROR Int.b0.1:	 Trying to access a sector of the b-functions
		 that should not be relevant for this calculation.
		 (See troubleshooting note #2 below)

ERROR Int.b1.1:	 Trying to access a sector of the b-functions
		 that should not be relevant for this calculation.
		 (See troubleshooting note #2 below)
ERROR Int.b1.2:  b1 called for s = m1 = m2 = 0

branchingratio.f95
------------------
ERROR Br.Fs.1: 	Invalid array argument: d=1, s=2
		The argument 'a' must either be 1 or 2
		corresponding to the B_d or B_s meson.
ERROR Br.Fp.1: 	See above.
ERROR Br.Fv.1: 	See above.
ERROR Br.Fa.1: 	See above.

! ----------------------------------------------------------------
! 7. Troubleshooting
! ----------------------------------------------------------------

If the above error messages don't solve your problem, here are some
thoughts I've had while writing up the code.

1.  I've used REAL(KIND=DBL) instead of REAL. This is for double precision.
	DBL is defined in program_parameters.f95. This should be a 
	processor-independent definition. A unfortunate consequence, however,
	is that REAL constants (e.g. 1.D0) are NOT automatically typecasted into
	REAL(KIND=DBL) constants when they are fed into the argument of a function.
	(It is still ok to define a REAL(KIND=DBL) with a REAL and the conversion
	is automatic.)

	This should not really be a problem since none of my functions should be
	called with non-names constants. However, if one is getting type errors
	with some of the functions, try putting in an explicit type conversion:
		REAL(c,DBL) turns the real part of c into a REAL(KIND=DBL)

	One should be careful with entering zero, however. Be sure to enter
	0.D0, not just 0. If this becomes tedious, I should define a 
	parameter REAL(KIND=DBL)::ZERO=0.D0.

2.  The b-functions (self energies) are a very rough hack. (`kludge') They're not
	general and I've cut out several options from Janusz Rosiek's more complete
	set of options. Hopefully this hasn't caused any errors (some of which are
	tagged above). There are a few things that I've removed that I don't really
	understand. I don't think they're relevant for the code. These are the weird
	dbstat and infstat logical switches. I don't know what they do, but I think
	the code will survive without them. Just like the human appendix.

3.  If there are any I/O errors, I list the I/O usage below:
	subroutine createlogfile uses IOUNIT=25
	   This should only be called at the end of a program, after
	   all other IO streams have been closed. There shouldn't
	   be any errors. 

4.  The b-functions in Janusz's code are complex, while I've left them real.
	I don't think this will cause problems, but who knows.

5.  NOTE: Janusz's CKM matrix is my (and the rest of the world's) CKM^\dagger
	So Janusz: CKM(i,j) = Flip: DCONJG(CKM(j,i))

! ----------------------------------------------------------------
! 8. Troubleshooting
! ----------------------------------------------------------------

1. Janusz's code has several versions of CKM floating around. The
	CKM in common/ckm is the true CKM. The CKM in common/km_mat/
	is the dagger of the true CKM. Note that the common/km_mat/
	is what is used in the calculations, even though the
	common/ckm is what is listed in bll.f (Why?)

2. Janusz's fermion masses are defined at some scale (call it S1) in
	bll.f. However, the masses are re-initialised when qmass_nlo
	(in qcd_fun.f) is called by b_ll (in phen_2q.f). The mass that
	is outputted after running b_ll is different from the masses
	inputted at the top of bll.f. However, the Yukawa couplings:
		yu(3)	=  um*sqrt(2.D0)/v2
	are defined in terms of the pre-re-initialised (!) masses
	BEFORE qmass_nlo is run. (This is probably right,b ut shouldn't
	these values run also?)  By the way, where are the scales 
	amuu, amud defined?

3. Known discrepancies with the code. The current version of the code 
	has physics_data taken from Janusz's values. (After appropriate
	running of quark masses.) However, there are still some 
	discrepancies with my code. A few comments on this are
	in order.

	(a) The discrepancy is at the 5th (or higher) decimal
	(b) The main source is a discrepancy in yh_eff_l/r
	    I believe this is from some hidden higher order
	    redefinition in Janusz's code that I didn't find.
	(c) Janusz's branching ratio changes slightly with 
	    del (the divergence from dimensional regularisation
	    which cancels in a physical process). I believe
	    this is just due to rounding error since there is
	    no dependence on amiu ('t Hooft mass). Recall
	    that del and amiu always appear together in the
	    combination del + log(amiu). By the way, my
	    del and amiu = mu_tH are defined in tuning_data.f95,
	    and the code has no del or mu_tH dependence.
	(d) Janusz's code seems to give slightly different values
	    if you run it consecutively. I'm not sure why this is,
	    but I suspect it partially has to do with the running of
	    eta_2q_evol only during the first call of b_ll. Here's
	    an example of two consecutive runs:

		Br(B_s -> ll) =   4.507764590846329E-009
		Br(B_s -> ll) =   4.507764640148487E-009


! ----------------------------------------------------------------
! 9. Version Notes
! ----------------------------------------------------------------


Bmumu_29Nov07.tar
-----------------
This is a `full' matching with Janusz's code. See the corresponding
JcodeWorking_29Nov07Special.tar file. If you run both programs, the 
output should be identical.

The results show minor  disagreement in the yh_eff_l/r functions, 
which I presume come from some residual definition of yh_eff_l/r 
that I couldn't find. Supporting this hypothesis is the fact that
the diagrams without yh_eff_l/r dependence appear to have exact
agreement.

Also, there are occasionally some differences between the values
of divergent diagrams. However, these differences seem to cancel
consistently and the amplitude matches to the accuracy that one
would expect given the differences in yh_eff_l/r.

See the files MatchingResults_28Nov07J/F.TXT for comparisons
of the outputs.
