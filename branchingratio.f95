! ----------------------------------------------------------------
! Project: Complete One loop Predictions for B->mu+mu- in the MSSM
! Please cite: Dedes, Rosiek, Tanedo (2007)
! ----------------------------------------------------------------
! Main file for B->mu+mu- calculation
! filename: amplitude.f95
! Current version: 0.01
! By Flip Tanedo (philip.tanedo@gmail.com)
! Date: 22 November 2007
! Latest version available at: http://www.dur.ac.uk/philip.tanedo
! ----------------------------------------------------------------
! This module calculates the B->mu+mu- branching ratio
!	SM only for now! 
!
!	Arguments: 
!	  i is the antiquark flavour (s or d)
! 	  j is the quark flavour (b)
!	  k and l are the lepton flavors
!
!	NOTE: error if i is not 1 or 2
!	      error if j is not 3
!	(Calculating B decay only.)
!
! See our paper for descriptions of the functions Fs, Fp, Fv, Fa.
! Do NOT confuse with the decay constants fB(2).
!
! ----------------------------------------------------------------
! RECORD OF REVISIONS
!	Date	 Programmer	Description
!	====	 ==========	==================================
!	22/11/07 Flip Tanedo	Original code
! ----------------------------------------------------------------
! USAGE
! USE branchingratio
! ----------------------------------------------------------------
! PRODUCTION NOTES
! 1. In the calculation of the F's, the arguments 'i,j,k,l' are
!	fermion flavour indices as usual.
!	The integer a labels the B -d or -s (1 or 2) meson
! 2. There may be a better way to do this,
!	if so, I should do it that way instead.
! 3. I can get rid of a = i (assuming i is the antiquark
!	and the B meson is b\bar{q} q = d,s
! ----------------------------------------------------------------

MODULE branchingratio
USE program_parameters		! for DBL, Pi, i = img
USE physics_data		! for coefficients
! USE integrals			! for integrals, of course
! USE tuning_data		! for mu_tH, 't Hooft mass, and del
! USE sorter			! for array sorting
! USE errorlog			! for error logging
! USE diagrams_box		! for box diagrams
! USE diagrams_zpenguin		! for Z-penguin piece
! USE diagrams_dselfenergy	! for d-self energy piece
USE amplitude
IMPLICIT NONE

CONTAINS



! ================================================================
! 	Calculate the F coefficients
! ================================================================

	COMPLEX(KIND=DBL) FUNCTION Fs(i,j,k,l)
!	Scalar coefficient
	INTEGER, INTENT(IN) :: i,j,k,l

	Fs=(0.D0,0.D0)

	Fs = dl_sll(i,j,k,l) + dl_slr(i,j,k,l) &
&		- dl_srr(i,j,k,l) - dl_srl(i,j,k,l)
	Fs = Fs * img * Mb(j)**2 * fB(j)/4/(dm(j)+dm(i))

!	Temp: troubleshooting
!	Fs = (-0.365258930D-13, 0.745703829D-15)

	END FUNCTION Fs

! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION Fp(i,j,k,l)
!	Pseudoscalar coefficient
	INTEGER, INTENT(IN) :: i,j,k,l

	Fp=(0.D0,0.D0)

	Fp = -dl_sll(i,j,k,l) + dl_slr(i,j,k,l) &
&		- dl_srr(i,j,k,l) + dl_srl(i,j,k,l)
	Fp = Fp * img * Mb(j)**2 * fB(j)/4/(dm(j)+dm(i))

!	TEMP
!	Fp = (0.393407001D-13, -0.803170253D-15)


	END FUNCTION Fp

! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION Fv(i,j,k,l)
!	Vector coefficient
	INTEGER, INTENT(IN) :: i,j,k,l
	COMPLEX(KIND=DBL) dls_vll, dls_vrr, dls_vlr, dls_vrl

	Fv=(0.D0,0.D0)

	Fv = dl_vll(i,j,k,l) + dl_vlr(i,j,k,l) &
&		- dl_vrr(i,j,k,l) - dl_vrl(i,j,k,l)
	Fv = -Fv*img*fB(j)/4

!	TEMP
!	Fv = (0.123719347D-09, -0.252593635D-11)

	END FUNCTION Fv


! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION Fa(i,j,k,l)
!	Axial coefficient
	INTEGER, INTENT(IN) :: i,j,k,l

	Fa=(0.D0,0.D0)

	Fa = - dl_vll(i,j,k,l) + dl_vlr(i,j,k,l) &
&		- dl_vrr(i,j,k,l) + dl_vrl(i,j,k,l)
	Fa = -img * Fa *fB(j)/4

!	TEMP
!	Fa = (-0.460921442D-09, 0.941108211D-11)


	END FUNCTION Fa


! ================================================================
! 	Calculate the |M|^2
! ================================================================

	COMPLEX(KIND=DBL) FUNCTION amp2(i,j,k,l)
	INTEGER, INTENT(IN) :: i,j,k,l
	COMPLEX(KIND=DBL)::amp2t

	amp2t = (0.D0,0.D0)
	amp2t = amp2t + 2*abs(Fs(i,j,k,l))**2*(MB(j)**2-(em(l)+em(k))**2)
	amp2t = amp2t + 2*abs(Fp(i,j,k,l))**2*(MB(j)**2-(em(l)-em(k))**2)	
	amp2t = amp2t + 2*abs(Fv(i,j,k,l))**2*(MB(j)**2*(em(k)-em(l))**2 - (em(k)**2-em(l)**2)**2)
	amp2t = amp2t + 2*abs(Fa(i,j,k,l))**2*(MB(j)**2*(em(k)+em(l))**2 - (em(k)**2-em(l)**2)**2)
	amp2t = amp2t + 4*REAL(Fs(i,j,k,l)*DCONJG(Fv(i,j,k,l)), KIND=DBL) &
&			*( MB(j)**2*(em(l)-em(k)) + (em(k)**2 - em(l)**2)*(em(l)+em(k)) )
	amp2t = amp2t + 4*REAL(Fp(i,j,k,l)*DCONJG(Fa(i,j,k,l)), KIND=DBL) &
&			*( MB(j)**2*(em(l)+em(k)) + (em(k)**2 - em(l)**2)*(em(l)-em(k)) )

	amp2 = amp2t

	return
	END FUNCTION amp2

! ================================================================
! 	Calculate the Branching Ratio
! ================================================================

	REAL(KIND=DBL) FUNCTION br(i,j,k,l)
	INTEGER, INTENT(IN) :: i,j,k,l
!	TEMP:
!	REAL(KIND=DBL)::tauBee

	COMPLEX(KIND=DBL) :: brc ! `complex' branching ratio for kind-fitting
	REAL(KIND=DBL) :: temp1 
	REAL(KIND=DBL) :: temp2

	temp1 = 1.D0 - ((em(k)+em(l))/MB(j))**2	
	temp2 = 1.D0 - ((em(k)-em(l))/MB(j))**2

!	Matching (TEMP): to use the same lifetime as Janusz
!	tauBee =  2227246220542.82D0
!	I've now put this value into physics_data by hand.

!	Matching: to use the same lifetime as Janusz
	brc = tauB(j)/16/PI*amp2(i,j,k,l)/MB(j)*sqrt(temp1)*sqrt(temp2)
!	brc = tauBee/16/PI*amp2(i,j,k,l)/MB(j)*sqrt(temp1)*sqrt(temp2)

	br = abs(brc)
	return
	END FUNCTION br	

END MODULE branchingratio
