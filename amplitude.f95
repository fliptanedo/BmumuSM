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
! This module calculates the B->mu+mu- amplitude
!	SM only for now! 
!
!	Arguments: 
!	  i and j are the quark flavours
!	  k and l are the lepton flavors
!
!	This code is mainly from dd_ll.f (see bottom of the file)
!		Janusz Rosiek, 02-08-2007
!
! ----------------------------------------------------------------
! RECORD OF REVISIONS
!	Date	 Programmer	Description
!	====	 ==========	==================================
!	22/11/07 Flip Tanedo	Original code
! ----------------------------------------------------------------
! USAGE
! USE amplitude
! ----------------------------------------------------------------
! PRODUCTION NOTES
! 1. Diagrams are located in the diagrams_%%%.f95 modules
! 2. Some modules have complete amplitudes (e.g. diagrams_box.f95)
! 3. Some modules only have pieces of amplitudes (e.g. Z-penguin)
!	For these we stick on the other factors in this file.
! ----------------------------------------------------------------
! DIRECTORY OF FUNCTIONS
!
! Wilson coefficient format: dl_%**
!  where % = vector, scalar, tensor
!  and  ** = {L,R} x {L,R}
!
!	dl_vll	
!	dl_vlr	
!	dl_vrl	
!	dl_vrr	
!	dl_sll	
!	dl_slr	
!	dl_srl	
!	dl_srr	
!	dl_tl	
!	dl_tr	
!
! ----------------------------------------------------------------

MODULE amplitude
USE program_parameters		! for DBL, Pi
USE physics_data		! for coefficients
! USE integrals			! for integrals, of course
! USE tuning_data		! for mu_tH, 't Hooft mass, and del
! USE sorter			! for array sorting
!USE errorlog			! for error logging
USE diagrams_box		! for box diagrams
USE diagrams_zpenguin		! for Z-penguin piece
USE diagrams_dselfenergy	! for d-self energy piece
IMPLICIT NONE

CONTAINS

! ================================================================
! 	Sum of the box, Z- and Higgs penguin contributions 
! ================================================================

!	Huh??? Where are the Higgs penguin contributions?

! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_vll(i,j,k,l)
!     	Full A^V_LL formfactor
	INTEGER, INTENT(IN) :: i,j,k,l
      	dl_vll = dl_vll_box(i,j,k,l)
      	if (k.eq.l) then
      	   dl_vll = dl_vll - e*(1.D0 - 2.D0*st2)/2.D0/sct/zm2*(zdd_vl(i,j) &
&     	        - e/2.D0/sct*(1.D0 - 2.d0/3.D0*st2)*(dv_sig0(i,j) - da_sig0(i,j)))
      	end if

!	Troubleshooting: (JcodeWorking value)
!	dl_vll = (-0.503579849D-08, 0.102909249D-09)

      	return
      	END FUNCTION dl_vll

! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_vlr(i,j,k,l)
!     	Full A^V_LR formfactor
	INTEGER, INTENT(IN) :: i,j,k,l
      	dl_vlr = dl_vlr_box(i,j,k,l)
      	if (k.eq.l) then
      	   dl_vlr = dl_vlr + e*st/ct/zm2*(zdd_vl(i,j) &
&     	        - e/2.D0/sct*(1.D0 - 2.d0/3.D0*st2)*(dv_sig0(i,j) - da_sig0(i,j)))
      	end if

!	Troubleshooting
!	dl_vlr = (0.290544226D-08, -0.593574207D-10)

      	return
      	END FUNCTION dl_vlr

! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_vrl(i,j,k,l)
!     	Full A^V_RL formfactor
	INTEGER, INTENT(IN) :: i,j,k,l
      	dl_vrl = dl_vrl_box(i,j,k,l)
      	if (k.eq.l) then 
      	   dl_vrl = dl_vrl - e*(1.D0 - 2.D0*st2)/2.D0/sct/zm2*(zdd_vr(i,j) &
&     	        + e*st/3.D0/ct*(dv_sig0(i,j) + da_sig0(i,j)))
      	end if

!	Troubleshooting
!	dl_vrl = (-0.153947670D-13, 0.314828671D-15)

      	return
      	END FUNCTION dl_vrl

! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_vrr(i,j,k,l)
!     	Full A^V_RR formfactor
	INTEGER, INTENT(IN) :: i,j,k,l
      	dl_vrr = dl_vrr_box(i,j,k,l)
      	if (k.eq.l) then
      	   dl_vrr = dl_vrr + e*st/ct/zm2*(zdd_vr(i,j) &
&     	        + e*st/3.D0/ct*(dv_sig0(i,j) + da_sig0(i,j)))
      	end if

!	Troubleshooting
!	dl_vrr = (0.123520522D-13, -0.252604030D-15)

      	return
      	END FUNCTION dl_vrr

! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_sll(i,j,k,l)
!     	Full A^S_LL formfactor
	INTEGER, INTENT(IN) :: i,j,k,l
      	dl_sll = dl_sll_box(i,j,k,l)
      	return
      	END FUNCTION dl_sll

! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_slr(i,j,k,l)
!     	Full A^S_LR formfactor
	INTEGER, INTENT(IN) :: i,j,k,l
      	dl_slr = dl_slr_box(i,j,k,l)

!	Troubleshooting
!	dl_slr = (0.254354361D-14, -0.520164478D-16)

      	return
      	END FUNCTION dl_slr

! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_srl(i,j,k,l)
!     	Full A^S_RL formfactor
	INTEGER, INTENT(IN) :: i,j,k,l
      	dl_srl = dl_srl_box(i,j,k,l)

!	Troubleshooting
!	dl_srl = (0.685553149D-13, -0.140198263D-14)

      	return
      	END FUNCTION dl_srl

! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_srr(i,j,k,l)
!     	Full A^S_RR formfactor
	INTEGER, INTENT(IN) :: i,j,k,l
      	dl_srr = dl_srr_box(i,j,k,l)
      	return
      	END FUNCTION dl_srr

! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_tl(i,j,k,l)
!     	Full A^T_L formfactor
	INTEGER, INTENT(IN) :: i,j,k,l
      	dl_tl = dl_tl_box(i,j,k,l)
      	return
      	END FUNCTION dl_tl

! ----------------------------------------------------------------

	DOUBLE COMPLEX FUNCTION dl_tr(i,j,k,l)
!     	Full A^T_R formfactor
	INTEGER, INTENT(IN) :: i,j,k,l
      	dl_tr = dl_tr_box(i,j,k,l)
      	return
      	END FUNCTION dl_tr


! ----------------------------------------------------------------

	
END MODULE amplitude
