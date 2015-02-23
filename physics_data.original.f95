! ----------------------------------------------------------------
! Project: Complete One loop Predictions for B->mu+mu- in the MSSM
! Please cite: Dedes, Rosiek, Tanedo (2007)
! ----------------------------------------------------------------
! Main file for B->mu+mu- calculation
! filename: physics_data.f95
! Current version: 0.01
! By Flip Tanedo (philip.tanedo@gmail.com)
! Date: 16 November 2007
! Latest version available at: http://www.dur.ac.uk/philip.tanedo
! ----------------------------------------------------------------
! This file contains physics data. Should also include subroutines
! for the running couplings and masses. 
! All data using Planck Units. (h-bar = c = G = k =1)
! Data is from the PDG unless otherwise noted.
!
! But I still need to work on running values. (alpha, masses, ...)
! ----------------------------------------------------------------
! RECORD OF REVISIONS
!	Date	 Programmer	Description
!	====	 ==========	==================================
!	18/11/07 Flip Tanedo	Pre-Original Code 
!	24/11/07 Flip Tanedo	Agreement with Janusz's Code
! ----------------------------------------------------------------
! DATA INCLUDED
!
!	NAME		DESCRIPTION
!	==============	==========================================
!	hbar		hbar, for conversion of from seconds to GeV
!	um(3)		up-type quark masses (@ appropriate scale)
!	dm(3)		down-type quark masses (@ appropriate scale)
!	em(3		charged lepton masses (@ appropriate scale)
!	wm		W boson mass
!	zm		Z boson mass
!	st		sin(\theta_W)
!	ct		cos(\theta_W)	(Comes from wm and zm)
!	e		unit of electric charge (from alpha_fine)
!	ckm(3,3)	CKM matrix (careful: janusz defines as CKM^dag)
!	 phi		Phase for CKM
!	 s12		Angle for CKM
!	 s23		Angle for CKM
!	 s13		Angle for CKM
!	v		SM Higgs mass
!	v1		SUSY Higgs 1 mass (= v in SM)
!	v2		SUSY Higgs 2 mass (= v in SM)
!	cm(2)		Charged higgs masses (Currently SM value only)
!	zh(2,2)		Charged higgs mixing from CP to R/L basis
!	fB(2)		B(d,s) decay constants, $f_{B_{d,s}}$
!	MB(2)		Mass of the $B^0_{d,s}$ meson.
!	tauB(2)		B(d,s) lifetime
! ----------------------------------------------------------------
! DATA INCLUDED, DERIVED QUANTITIES
! (For porting amplitudes from Janusz's code)
!
!	NAME		DESCRIPTION
!	==============	==========================================
!	st2		sin(\theta_W)^2
!	ct2		cos(\theta_W)^2
!	sct		sin(\theta_W)cos(\theta_W)
!	sct2		sin(\theta_W)^2 cos(\theta_W)^2
!	e2		electric charge ^2
!	alpha		Fine structure constant
!	wm2		W mass squared
!	zm2		Z mass squared
!	yu,yd,yl(3)	Up-, down-, lepton- yukawa couplings
!	
! ----------------------------------------------------------------
! SUBROUTINES
!
!	NAME		DESCRIPTION
!	==============	==========================================
!	ckm_init	Initialises CKM matrix, driver MUST run this
!
! ----------------------------------------------------------------
! FUNCTIONS (Though these might as well be redefined as arrays)
!
!	NAME		DESCRIPTION
!	==============	==========================================
!	yh_eff_l,r	`Effective' left-right charged Higgs Yukawa,
!			contains CKM matrix. It's kind of an odd
!			quantity, but used for Janusz-compatibility.
!
! ----------------------------------------------------------------
! DATA NOT INCLUDED, see TUNING_DATA.F95
!	mu_tH	'	t Hooft mass scale, should cancel 
!	del		2/(d-4), should cancel (for integrals)
!	eps		"small mass difference" for integrals
! ----------------------------------------------------------------
! DATA NOT INCLUDED, see PROGRAM_PARAMETERS.F95
!	PI	'	3.14... (for Integrals, but not really)
!	sq2		SQRT(REAL(2,DBL))
! ----------------------------------------------------------------
! PRODUCTION NOTES
! 1. Note: fB(2) and MB(2) are arrays.
!	element 1 corresponds to the B_d system
!	element 2 corresponds to the B_s system
! 2. CURRENT VERSION: 
! ----------------------------------------------------------------


MODULE physics_data

USE program_parameters		! for DBL and PI
! USE errorlog
IMPLICIT NONE
SAVE

! The SAVE statement is important, it allows the driver file
! to call CKM_INIT and for the values to be maintained 
! when the other modules access PHYSICS_DATA


!	Planck units used. All mass dimensions are given in GeV
	REAL(KIND=DBL), PARAMETER:: hbar	= 6.58211915D-25 ! GeV sec

!	Electric coupling (dimensionless)
! 	----------------------------------------------------------------
	REAL(KIND=DBL), PARAMETER:: alpha 	= 1.D0/127.861285D0      ! @ Z scale
	REAL(KIND=DBL), PARAMETER:: e2		= alpha*4.D0*PI
	REAL(KIND=DBL), PARAMETER:: e		= sqrt(e2)

!	Up and down type masses, check this for 1S vs MS-bar mass
! 	----------------------------------------------------------------
!	Important note: see README.TXT section 8
!	TOP AND BOTTOM MASSES (To match Janusz's numbers)
	REAL(KIND=DBL), PARAMETER:: um(3)	= (/0.004D0, 1.3D0, 161.221993663823D0/)
	REAL(KIND=DBL), PARAMETER:: dm(3)	= (/0.007D0, 0.11D0, 2.96479470613322D0/)
	REAL(KIND=DBL), PARAMETER:: em(3)	= (/0.000511D0, 0.1056590D0, 1.7769D0/)

!	Up and down type masses FOR YUKAWAS ONLY
! 	----------------------------------------------------------------
!	Important note: see README.TXT section 8
!	REAL(KIND=DBL), PARAMETER:: umy(3)	= (/0.004D0,1.3D0,165.0D0/)
!	REAL(KIND=DBL), PARAMETER:: dmy(3)	= (/0.007D0, 0.11D0,2.72D0/)


!	Gauge Boson Masses
! 	----------------------------------------------------------------
	REAL(KIND=DBL), PARAMETER:: wm		= 80.4D0
	REAL(KIND=DBL), PARAMETER:: zm		= 91.1863D0
	REAL(KIND=DBL), PARAMETER:: wm2	= wm**2
	REAL(KIND=DBL), PARAMETER:: zm2	= zm**2

!	Higgs Boson Masses
! 	----------------------------------------------------------------
	REAL(KIND=DBL) 	   	:: cm(2)	= (/50000.06464D0,wm/)	! Charged Higgs
!	  This is a TEMPORARY, SM-ONLY assignment!
!	  We'll need to include an initialisation subroutine for SUSY.

!	Weak mixing angle
! 	----------------------------------------------------------------
	REAL(KIND=DBL), PARAMETER:: ct		= wm/zm
	REAL(KIND=DBL), PARAMETER:: ct2		= ct*ct
	REAL(KIND=DBL), PARAMETER:: st2		= 1.D0 - ct2
	REAL(KIND=DBL), PARAMETER:: st		= sqrt(st2)
	REAL(KIND=DBL), PARAMETER:: sct		= st*ct
	REAL(KIND=DBL), PARAMETER:: sct2	= sct*sct


!	CKM Matrix Declaration, Format: CKM (row, col)
! 	----------------------------------------------------------------
!	Fill this in using ckm_init subroutine below
	COMPLEX(KIND=DBL):: CKM(3,3) 		= &
	  RESHAPE((/0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0/),(/3,3/))
!	NOTE: This is not a PARAMETER. It is initialised by
!	  the subroutine ckm_init, which you must run in the driver!

	REAL(KIND=DBL), PARAMETER :: phi	= 7.D0/18*PI ! from Janusz	
	REAL(KIND=DBL), PARAMETER :: s12	= 0.2272D0   ! from Janusz	
	REAL(KIND=DBL), PARAMETER :: s23	= 0.04221D0  ! from Janusz	
	REAL(KIND=DBL), PARAMETER :: s13	= 0.00396D0  ! from Janusz


!	STANDARD MODEL (only) VEVS -- REPLACE THIS WITH SUSY VEVS
! 	----------------------------------------------------------------
	REAL(KIND=DBL), PARAMETER :: v		= wm*2*st/e*SIN(PI/4)
!	The factor of SIN(PI/4) = sin(beta) = cos(beta) in SM
!	in mh_diag.f of Janusz's code,
!		v1 = 2*wm*st/e*cb
!		v2 = 2*wm*st/e*sb
!		sb = zh(1,1)
!		cb = zh(2,1)

!	SUSY HIGGS SECTOR
! 	----------------------------------------------------------------
	REAL(KIND=DBL), PARAMETER:: beta	= PI/4

!	SUPERSYMMETRIC VEVS -- SM VALUES FOR NOW
! 	----------------------------------------------------------------
	REAL(KIND=DBL), PARAMETER :: v1		= v	     ! SM value
	REAL(KIND=DBL), PARAMETER :: v2		= v	     ! SM value

!	CHARGED HIGGS MIXING MATRIX (CP to RL basis)
! 	----------------------------------------------------------------
	REAL(KIND=DBL), PARAMETER :: zh(2,2)	= &
&	  RESHAPE((/v2,v1,-v1,v2/), (/2,2/))/sqrt(v1*v1 + v2*v2)


!	YUKAWA COUPLINGS (cf vf_def.f in Janusz's code)
!	Important note: see README.TXT section 8
!	umy, dmy, emy are temporary values to put in the right
!	values for quark mass running, which I haven't built in yet.
! 	----------------------------------------------------------------
!	REAL(KIND=DBL), PARAMETER :: yu(3)	=  umy*sqrt(2.D0)/v2
!	REAL(KIND=DBL), PARAMETER :: yd(3)	= -dmy*sqrt(2.D0)/v1
	REAL(KIND=DBL), PARAMETER :: yu(3)	= (/0.330591253D-4, 0.107442157D-1, 0.133246452D+1/)
	REAL(KIND=DBL), PARAMETER :: yd(3)	= (/-0.578534761D-4, -0.909126053D-3, -0.245033828D-1/)
	REAL(KIND=DBL), PARAMETER :: yl(3)	= -em*sqrt(2.D0)/v1


!	B-meson properties (nonperturbative)
! 	----------------------------------------------------------------
	REAL(KIND=DBL), PARAMETER :: fB(2)	= (/0.202D0, 0.232D0/)
	REAL(KIND=DBL), PARAMETER :: MB(2)	= (/5.2794D0, 5.368D0/)
!		this is amb in phen_4q.f
	REAL(KIND=DBL), PARAMETER :: tauB(2)	= (/1.530D-12/hbar,2227246220542.82D0/)
!	! tauB(2) SET TO MATCH JANUSZ's parameters

! ----------------------------------------------------------------
CONTAINS
! ----------------------------------------------------------------

!	Initialize CKM
!	Based on Janusz's CKM_INIT function in mh_init.f
	SUBROUTINE ckm_init

!		COMPLEX(KIND=DBL) :: ckm_herm(3,3)
!		I'll make this a global module variable

		REAL(KIND=DBL) :: sith = s12
		REAL(KIND=DBL) :: coth = sqrt(1 - s12*s12)
		REAL(KIND=DBL) :: sibe = s13
		REAL(KIND=DBL) :: cobe = sqrt(1 - s13*s13)
		REAL(KIND=DBL) :: siga = s23
		REAL(KIND=DBL) :: coga = sqrt(1 - s23*s23)
		REAL(KIND=DBL) :: vkm1
		REAL(KIND=DBL) :: vkm2
		REAL(KIND=DBL) :: vkm3
		COMPLEX(KIND=DBL) :: cphi = exp(cmplx(0.d0,phi))
		INTEGER ::  i,j

		vkm1 = sith*cobe
		vkm2 = coth*cobe
		vkm3 = siga*cobe

      		CKM(1,1) = cobe*coth
      		CKM(1,2) = cobe*sith
      		CKM(1,3) = sibe/cphi
      		CKM(2,1) = - siga*coth*sibe*cphi - sith*coga
      		CKM(2,2) = coga*coth - siga*sibe*sith*cphi
      		CKM(2,3) = siga*cobe
      		CKM(3,1) = - sibe*coga*coth*cphi + siga*sith
      		CKM(3,2) = - coga*sibe*sith*cphi - siga*coth
      		CKM(3,3) = coga*cobe

!	UPDATE: I now have discovered a very subtle point in
!	Janusz's code!! 
!
!	 tCKM: the true CKM		
!	 jCKM: from common/ckm		= tCKM	
!	 kCKM: from common/km_mat	= tCKM^\dagger
!
!	 kCKM is what's ACTUALLY used in the calculation
!	 CKM in this subroutine (above) is tCKM
		RETURN

	END SUBROUTINE ckm_init

! ----------------------------------------------------------------

!	'EFFECTIVE' CHARGED HIGGS YUKAWA COUPLING
!	This translates the charged Higgs Yukawa couplings from the
!	  CP-even/odd basis to the Left/Right helicity basis.
!	  (i.e. from a basis where the fermion couplings are
!	  proportional to 1 and \gamma_5 to a basis where those
!	  couplings are proportional to P_R and P_L.)
!	Arguments: 
!	  i and j are the CKM matrix indices
!	  k is the index for the CP even and CP odd Higgs
!	NOTE: recall that Janusz's CKM is my CKM^\dag
!	NOTE: alternately, I could have defined this as an array
	COMPLEX(KIND=DBL) FUNCTION yh_eff_l(i,j,k)
		! Flavour indices i,j,k
		INTEGER, INTENT(IN)::i,j,k

		yh_eff_l = yu(j)*zh(2,k)*CKM(j,i) 
	END FUNCTION yh_eff_l
	
! ----------------------------------------------------------------

!	'EFFECTIVE' CHARGED HIGGS YUKAWA COUPLING
!	See notes above for yh_eff_r
	COMPLEX(KIND=DBL) FUNCTION yh_eff_r(i,j,k)
		! Flavour indices i,j,k
		INTEGER, INTENT(IN)::i,j,k

		yh_eff_r = -yd(i)*zh(1,k)*CKM(j,i)
	END FUNCTION yh_eff_r

! ----------------------------------------------------------------




END MODULE physics_data
