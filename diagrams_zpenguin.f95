! ----------------------------------------------------------------
! Project: Complete One loop Predictions for B->mu+mu- in the MSSM
! Please cite: Dedes, Rosiek, Tanedo (2007)
! ----------------------------------------------------------------
! Main file for B->mu+mu- calculation
! filename: diagrams_zpenguin.f95
! Current version: 1.00
! By Flip Tanedo (philip.tanedo@gmail.com)
! Date: 22 November 2007
! Latest version available at: http://www.dur.ac.uk/philip.tanedo
! ----------------------------------------------------------------
! This module contains the coefficients for the 
! Zdd penguin diagram. It is translated from dd_ll.f in 
! Janusz Rosiek's original code. 
!
! I would like to rewrite this calculation in the future, but
! for now I'll just port it directly since it's been checked
! several times already. (It will be easy to check new functions
! once I have this file set up.)
!
!	Arguments: 
!	  i and j are the quark flavours
!
!    From Janusz's header:
!
!      General form of the vertex is:                                  
!      i\bar d^J (A^V_L \gamma_mu P_L + A^V_R \gamma_mu P_R) d^I Z^mu  
!                         
!    By the way, I helped rewrite the original zdd_vert0.f in
!    Janusz's code. Notation for diagrams: 3 letter codes refer
!    to the particles running in the loop. For exampple, Huu
!    refers to a higgs and two up quarks in the triangle.
!	Original code: zdd_vert0.f, Rosiek et al.  22-08-2007
!
! The functions here jsut account for the actual penguin. The 
! full diagram for B->mu+mu- (i.e. adding Z propagator and
! lepton lines) is calculated in amplitude.f95
!
! ----------------------------------------------------------------
! RECORD OF REVISIONS
!	Date	 Programmer	Description
!	====	 ==========	==================================
!	22/11/07 Flip Tanedo	Original code
!	22/11/07 Flip Tanedo	Completed
!	29/11/07 Flip Tanedo	Cosmetic changes
! ----------------------------------------------------------------
! USAGE
! USE diagrams_penguin
! zdd_vl(i,j) or zdd_vr(i,j) gives the form factor contributions
! (amplitudes.f95 calls these and stick on the leptonic part)
! ----------------------------------------------------------------
! PRODUCTION NOTES
! 1. These are lifted from Janusz's code with minimal reformatting.
! 2. Physical constants are as defined in physics_data.f95
! 3. Functions are complex-valued (CKM matrix is complex)
! 4. Recall Janusz's CKM is the dagger of the 'true' CKM
! ----------------------------------------------------------------
! DIRECTORY OF FUNCTIONS
!	zdd_vl		Left-handed vertex
! 	zdd_vr 		Right-handed vertex
!	zdd_vl,r_*	Diagrams: * = g,h,hg
! ----------------------------------------------------------------

MODULE diagrams_zpenguin
USE program_parameters		! for DBL, Pi
USE physics_data		! For coefficients
USE integrals			! for integrals, of course
USE tuning_data			! for mu_tH, 't Hooft mass, and del
!USE sorter			! for array sorting
!USE errorlog			! for error logging
IMPLICIT NONE

CONTAINS

! ================================================================
! 	Formfactor A^V_L
! ================================================================

      	COMPLEX(KIND=DBL) FUNCTION zdd_vl_g(i,j)
!     	Gauge contribution
	INTEGER, INTENT(IN)::i,j
	INTEGER::m 

      	zdd_vl_g = (0.d0,0.d0)

      	do m=1,3
      	   zdd_vl_g = zdd_vl_g + CKM(m,i)*DCONJG(CKM(m,j)) &
!     	Wuu contribution
&     	        * (((1 - 4.d0/3*st2)*(cp1(wm,um(m),um(m)) - 0.5d0) &
!&     	        * (((1 - 4.d0/3*st2)*(1.D0 - 0.5d0) &
&     	        + 8.d0/3*st2*um(m)*um(m)*cp0(wm,um(m),um(m))) &
!     	WWu contribution
&     	        + ct2*(6*cp1(wm,wm,um(m)) - 1))
!&     	        + ct2*(6*1.D0 - 1))
      	end do

      	zdd_vl_g = e*e2/4/sct/st2*zdd_vl_g

!	TEMP BELOW
!	zdd_vl_g = CKM(3,i)*DCONJG(CKM(3,j))

      	return
      	END FUNCTION zdd_vl_g

! ----------------------------------------------------------------

      	COMPLEX(KIND=DBL) FUNCTION zdd_vl_hg(i,j)
!     	Higgs-gauge contribution
	INTEGER, INTENT(IN)::i,j
	INTEGER::m 

      	zdd_vl_hg = (0.d0,0.d0)

      	do m=1,3
      	   zdd_vl_hg = zdd_vl_hg + (CKM(m,i)*DCONJG(yh_eff_l(j,m,2)) &
&     	        + DCONJG(CKM(m,j))*yh_eff_l(i,m,2))*um(m)*cp0(wm,wm,um(m))
      	end do
!	Note: changed overall sign as per monthlong sign error
!	debacle from 9 Oct 2007.

      	zdd_vl_hg = e2*wm/sq2/ct*zdd_vl_hg
      	return
      	END FUNCTION zdd_vl_hg

! ----------------------------------------------------------------

      	COMPLEX(KIND=DBL) FUNCTION zdd_vl_h(i,j)
!     	Higgs contribution
	INTEGER, INTENT(IN)::i,j
	INTEGER::m,l 

      	zdd_vl_h = (0.d0,0.d0)

      	do m=1,3
!	   Loop truncated for SM only, for SUSY l=1,2 (see README)
      	   do l=1,2
      	     zdd_vl_h = zdd_vl_h - yh_eff_l(i,m,l)*DCONJG(yh_eff_l(j,m,l)) &
!     	Huu contribution
&     	           * ((2*st2/3.d0*(cp1(cm(l),um(m),um(m)) - 0.5d0) &
!&     	           * ((2*st2/3.d0*(1.D0 - 0.5d0) &
&     	           + (1 - 4*st2/3.d0)*um(m)*um(m)*cp0(cm(l),um(m),um(m))) &
!     	HHu contribution
&     	           + (ct2 - st2)/2*(cp1(cm(l),cm(l),um(m)) + 0.5d0))
!&     	           + (ct2 - st2)/2*(1.D0 + 0.5d0))
      	  end do
      	end do

      	zdd_vl_h = e/2/sct*zdd_vl_h
!	FLIP's next line
!	zdd_vl_h = yh_eff_l(j,3,1)

      	return
      	END FUNCTION zdd_vl_h

! ----------------------------------------------------------------
!	I'm skipping the SUSY contributions for now
!	* Chargino
!	* Neutralino
!	* Gluino
! ----------------------------------------------------------------

! ----------------------------------------------------------------
!	Full Z penguin A^V_L formfactor
! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION zdd_vl(i,j)
!     	Full Z penguin A^V_L formfactor
	INTEGER, INTENT(IN)::i,j

! 	I haven't put in the switches for now
      	zdd_vl = (zdd_vl_g(i,j) + zdd_vl_hg(i,j) + zdd_vl_h(i,j) &
&			)/16/pi/pi
      	return
      	END FUNCTION zdd_vl


! ================================================================
! 	Formfactor A^V_R
! ================================================================


	COMPLEX(KIND=DBL) FUNCTION zdd_vr_g(i,j)
!	Gaugino contribution
	INTEGER, INTENT(IN) :: i,j

	zdd_vr_g = (0.D0, 0.D0)
	return
	END FUNCTION zdd_vr_g

! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION zdd_vr_hg(i,j)
!	Higgs-gauge contribution
	INTEGER, INTENT(IN) :: i,j

	zdd_vr_hg = (0.D0, 0.D0)
	return
	END FUNCTION zdd_vr_hg

! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION zdd_vr_h(i,j)
!	Higgs contribution
	INTEGER, INTENT(IN) :: i,j
	INTEGER :: m,l

      	zdd_vr_h = (0.D0,0.D0)

      	do m=1,3
!	   Loop truncated for SM only, for SUSY l=1,2 (see README)
      	   do l=1,2
      	     zdd_vr_h = zdd_vr_h + yh_eff_r(i,m,l)*dconjg(yh_eff_r(j,m,l)) &
!     	Huu contribution
&     	           * (((1 - 4*st2/3.d0)*(cp1(cm(l),um(m),um(m)) - 0.5d0) &
&     	           + 2*st2/3.d0*um(m)*um(m)*cp0(cm(l),um(m),um(m))) &
!     	HHu contribution
&     	           - (ct2 - st2)*(cp1(cm(l),cm(l),um(m)) + 0.5d0))
      	  end do
      	end do

      	zdd_vr_h = e/4/sct*zdd_vr_h

      	return
      	END FUNCTION zdd_vr_h

! ----------------------------------------------------------------
!	I'm skipping the SUSY contributions for now
!	* Chargino
!	* Neutralino
!	* Gluino
! ----------------------------------------------------------------

! ----------------------------------------------------------------
!	Full Z penguin A^V_L formfactor
! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION zdd_vr(i,j)
!	Full Z penguin A^V_L formfactor
	INTEGER, INTENT(IN) :: i,j

! 	I haven't put in the switches for now
      	zdd_vr = (zdd_vr_g(i,j) + zdd_vr_hg(i,j) + zdd_vr_h(i,j) &
&     			)/16/pi/pi
      	return
      	END FUNCTION zdd_vr


END MODULE diagrams_zpenguin
