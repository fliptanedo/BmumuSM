! ----------------------------------------------------------------
! Project: Complete One loop Predictions for B->mu+mu- in the MSSM
! Please cite: Dedes, Rosiek, Tanedo (2007)
! ----------------------------------------------------------------
! Main file for B->mu+mu- calculation
! filename: diagrams_box.f95
! Current version: 0.9
! By Flip Tanedo (philip.tanedo@gmail.com)
! Latest version available at: http://www.dur.ac.uk/philip.tanedo
! ----------------------------------------------------------------
! This module contains the Wilson coefficients for the B->mu+mu-
! box diagrams. It is translated from dd_ll.f in Janusz Rosiek's
! original code. 
!
!	Arguments: 
!	  i and j are the quark flavours
!	  k and l are the lepton flavours
!
!     General form of the Hamiltonian is:                         
!     H = A^V_LL H^V_LL + A^V_RR H^V_RR                           
!       + A^V_LR H^V_LR + A^V_LR H^V_RL                           
!       + A^S_LL H^S_LL + A^S_RR H^S_RR                           
!       + A^S_LR H^S_LR + A^S_LR H^S_RL                           
!       + A^T_L  H^T_L  + A^T_R  H^T_R    
!
!     The original code: 
!	dd_ll.f Rosiek, 2-08-2007
!
! ----------------------------------------------------------------
! RECORD OF REVISIONS
!	Date	 Programmer	Description
!	====	 ==========	==================================
!	20/11/07 Flip Tanedo	Original code
!	29/03/08 Flip Tanedo	LFV corrections
! ----------------------------------------------------------------
! USAGE
! USE diagrams_box
! ----------------------------------------------------------------
! PRODUCTION NOTES
! 1. These are lifted from Janusz's code with minimal reformatting.
! 2. Physical constants are as defined in physics_data.f95
! 3. Functions are complex-valued (CKM matrix is complex)
! 4. Recall Janusz's CKM is the dagger of the 'true' CKM
! ----------------------------------------------------------------
! DIRECTORY OF FUNCTIONS
!	dl_v**_box =	A^V_** form factor, * = L,R
!	dl_s**_box = 	A^S_** form factor, * = L,R
!	dl_v,s**_% = 	Contribution from % = g, h, hg
! ----------------------------------------------------------------

MODULE diagrams_box
USE program_parameters		! for DBL, Pi
USE physics_data		! For coefficients
USE integrals			! for integrals, of course
USE tuning_data			! for mu_tH, 't Hooft mass, and del
IMPLICIT NONE

CONTAINS


! ================================================================
! 	Formfactor A^V_LL
! ================================================================

      	COMPLEX(KIND=DBL) FUNCTION dl_vll_g(i,j,k,l)
!     	Gauge contribution
	INTEGER, INTENT(IN)::i,j,k,l
	INTEGER::m 

     	dl_vll_g = (0.d0,0.d0)

      	if (k.ne.l) return
      	DO m=1,3
      	   dl_vll_g = dl_vll_g + CKM(m,i)*DCONJG(CKM(m,j)) &
&     	        * dp1(wm,wm,um(m),0.d0)
      	END DO
      	dl_vll_g = e2*e2/4/st2/st2*dl_vll_g
      	return
      	END FUNCTION dl_vll_g

! ----------------------------------------------------------------

      	COMPLEX(KIND=DBL) FUNCTION dl_vll_h(i,j,k,l)
!     	Double Higgs contribution
	INTEGER, INTENT(IN)::i,j,k,l

      	dl_vll_h = (0.d0,0.d0)
      	return
      	END FUNCTION dl_vll_h

! ----------------------------------------------------------------

      	COMPLEX(KIND=DBL) FUNCTION dl_vll_hg(i,j,k,l)
!     	W-Higgs contribution
	INTEGER, INTENT(IN)::i,j,k,l

      	dl_vll_hg = (0.d0,0.d0)

      	return
      	END FUNCTION dl_vll_hg

! ----------------------------------------------------------------
!	I'm skipping the SUSY contributions for now
!	* Chargino
!	* Neutralino
! ----------------------------------------------------------------

! ----------------------------------------------------------------
!	Full box A^V_LL formfactor
! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_vll_box(i,j,k,l)
!     	Full box A^V_LL formfactor
	INTEGER, INTENT(IN)::i,j,k,l

! 	I haven't put in the switches for now
	dl_vll_box = (dl_vll_g(i,j,k,l) + dl_vll_h(i,j,k,l) &
&          + dl_vll_hg(i,j,k,l))/16/pi/pi
	return
	END FUNCTION dl_vll_box


! ================================================================
! 	Formfactor A^V_RR
! ================================================================

	COMPLEX(KIND=DBL) FUNCTION dl_vrr_g(i,j,k,l)
!	Gauge contribution
	INTEGER, INTENT(IN)::i,j,k,l
	dl_vrr_g = (0.d0,0.d0)
	return
	END FUNCTION dl_vrr_g

! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_vrr_h(i,j,k,l)
!	Double Higgs contribution
	INTEGER, INTENT(IN)::i,j,k,l
	INTEGER::m,kk,ll
      	dl_vrr_h = (0.d0,0.d0)
      	if (k.ne.l) return
      	do m=1,3
!	   Loop truncated for SM only, for SUSY kk=1,2 (see README)
      	   do kk=1,2
!	      Loop truncated for SM only, for SUSY ll=1,2
      	      do ll=1,2
      	         dl_vrr_h = dl_vrr_h + zh(1,kk)*zh(1,ll)*yh_eff_r(i,m,kk) &
&     	              * DCONJG(yh_eff_r(j,m,ll)) &
&     	              * dp0(cm(kk),cm(ll),um(m),0.D0)
      	      end do
      	   end do
      	end do
      	dl_vrr_h = yl(k)*yl(l)/4*dl_vrr_h
	return
	END FUNCTION dl_vrr_h

! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_vrr_hg(i,j,k,l)
!     	W-Higgs contribution
	INTEGER, INTENT(IN)::i,j,k,l
      	dl_vrr_hg = (0.d0,0.d0)
      	return
      	END FUNCTION dl_vrr_hg

! ----------------------------------------------------------------
!	I'm skipping the SUSY contributions for now
!	* Chargino
!	* Neutralino
! ----------------------------------------------------------------

! ----------------------------------------------------------------
!	Full box A^V_RR formfactor
! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_vrr_box(i,j,k,l)
!	Full box A^V_RR formfactor
	INTEGER, INTENT(IN)::i,j,k,l

! 	I haven't put in the switches for now
	dl_vrr_box = (dl_vrr_g(i,j,k,l) + dl_vrr_h(i,j,k,l) &
&          + dl_vrr_hg(i,j,k,l))/16/pi/pi
      	return
      	END FUNCTION dl_vrr_box

! ================================================================
! 	Formfactor A^V_LR
! ================================================================

	COMPLEX(KIND=DBL) FUNCTION dl_vlr_g(i,j,k,l)
!	Gauge contribution
	INTEGER, INTENT(IN)::i,j,k,l
      	dl_vlr_g = (0.d0,0.d0)
      	return
      	END FUNCTION dl_vlr_g

! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_vlr_h(i,j,k,l)
!	Double Higgs contribution
	INTEGER, INTENT(IN)::i,j,k,l
	INTEGER::m,kk,ll
	dl_vlr_h = (0.d0,0.d0)
      	if (k.ne.l) return
      	do m=1,3
!	   Loop truncated for SM only, for SUSY kk=1,2 (see README)
      	   do kk=1,2
!	      Loop truncated for SM only, for SUSY ll=1,2 (see README)
      	      do ll=1,2
      	         dl_vlr_h = dl_vlr_h + zh(1,kk)*zh(1,ll)*yh_eff_l(i,m,kk) &
&     	              * dconjg(yh_eff_l(j,m,ll)) &
&     	              * dp0(cm(kk),cm(ll),um(m),0.d0)
      	      end do
      	   end do
      	end do
      	dl_vlr_h = yl(k)*yl(l)/4*dl_vlr_h
      	return
      	END FUNCTION dl_vlr_h

! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_vlr_hg(i,j,k,l)
!     	W-Higgs contribution
	INTEGER, INTENT(IN)::i,j,k,l

	dl_vlr_hg = (0.d0,0.d0)
      	return
      	END FUNCTION dl_vlr_hg

! ----------------------------------------------------------------
!	I'm skipping the SUSY contributions for now
!	* Chargino
!	* Neutralino
! ----------------------------------------------------------------

! ----------------------------------------------------------------
!	Full box A^V_LR formfactor
! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_vlr_box(i,j,k,l)
!     	Full box A^V_LR formfactor
	INTEGER, INTENT(IN)::i,j,k,l

! 	I haven't put in the switches for now
      	dl_vlr_box = (dl_vlr_g(i,j,k,l) + dl_vlr_h(i,j,k,l) &
&     	     + dl_vlr_hg(i,j,k,l))/16/pi/pi
      	return
      	END FUNCTION dl_vlr_box

! ================================================================
! 	Formfactor A^V_RL
! ================================================================

	COMPLEX(KIND=DBL) FUNCTION dl_vrl_g(i,j,k,l)
!     	Gauge contribution
	INTEGER, INTENT(IN)::i,j,k,l
      	dl_vrl_g = (0.d0,0.d0)
      	return
      	END FUNCTION dl_vrl_g

! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_vrl_h(i,j,k,l)      
!     	Double Higgs contribution
	INTEGER, INTENT(IN)::i,j,k,l
      	dl_vrl_h = (0.d0,0.d0)
      	return
      	END FUNCTION dl_vrl_h

! ----------------------------------------------------------------
      
	COMPLEX(KIND=DBL) FUNCTION dl_vrl_hg(i,j,k,l)
!     	W-Higgs contribution
	INTEGER, INTENT(IN)::i,j,k,l
      	dl_vrl_hg = (0.d0,0.d0)
      	return
      	END FUNCTION dl_vrl_hg

! ----------------------------------------------------------------
!	I'm skipping the SUSY contributions for now
!	* Chargino
!	* Neutralino
! ----------------------------------------------------------------

! ----------------------------------------------------------------
!	Full box A^V_LR formfactor
! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_vrl_box(i,j,k,l)
!     	Full box A^V_RL formfactor
	INTEGER, INTENT(IN)::i,j,k,l

! 	I haven't put in the switches for now
      	dl_vrl_box = (dl_vrl_g(i,j,k,l) + dl_vrl_h(i,j,k,l) &
&     	     + dl_vrl_hg(i,j,k,l))/16/pi/pi
      	return
      	END FUNCTION dl_vrl_box

! ================================================================
! 	Formfactor A^S_LL
! ================================================================

	COMPLEX(KIND=DBL) FUNCTION dl_sll_g(i,j,k,l)
!     	Gauge contribution
	INTEGER, INTENT(IN)::i,j,k,l
      	dl_sll_g = (0.d0,0.d0)
      	return
      	END FUNCTION dl_sll_g

! ----------------------------------------------------------------
      
	COMPLEX(KIND=DBL) FUNCTION dl_sll_h(i,j,k,l)
!     	Double Higgs contribution
	INTEGER, INTENT(IN)::i,j,k,l
      	dl_sll_h = (0.d0,0.d0)
      	return
      	END FUNCTION dl_sll_h

! ----------------------------------------------------------------
      
	COMPLEX(KIND=DBL) FUNCTION dl_sll_hg(i,j,k,l)
!     	W-Higgs contribution
	INTEGER, INTENT(IN)::i,j,k,l
      	dl_sll_hg = (0.d0,0.d0)
      	return
      	END FUNCTION dl_sll_hg

! ----------------------------------------------------------------
!	I'm skipping the SUSY contributions for now
!	* Chargino
!	* Neutralino
! ----------------------------------------------------------------

! ----------------------------------------------------------------
!	Full box A^S_LL formfactor
! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_sll_box(i,j,k,l)
!     	Full box A^S_LL formfactor
	INTEGER, INTENT(IN)::i,j,k,l

! 	I haven't put in the switches for now
      	dl_sll_box = (dl_sll_g(i,j,k,l) + dl_sll_h(i,j,k,l) &
&     	     + dl_sll_hg(i,j,k,l))/16/pi/pi
      	return
      	END FUNCTION dl_sll_box

! ================================================================
! 	Formfactor A^S_RR
! ================================================================

	COMPLEX(KIND=DBL) FUNCTION dl_srr_g(i,j,k,l)
!     	Gauge contribution
	INTEGER, INTENT(IN)::i,j,k,l
      	dl_srr_g = (0.d0,0.d0)
      	return
      	END FUNCTION dl_srr_g

! ----------------------------------------------------------------
      
	COMPLEX(KIND=DBL) FUNCTION dl_srr_h(i,j,k,l)
!     	Double Higgs contribution
	INTEGER, INTENT(IN)::i,j,k,l
      	dl_srr_h = (0.d0,0.d0)
      	return
      	END FUNCTION dl_srr_h

! ----------------------------------------------------------------
      
	COMPLEX(KIND=DBL) FUNCTION dl_srr_hg(i,j,k,l)
!     	W-Higgs contribution
	INTEGER, INTENT(IN)::i,j,k,l
      	dl_srr_hg = (0.d0,0.d0)
      	return
      	END FUNCTION dl_srr_hg

! ----------------------------------------------------------------
!	I'm skipping the SUSY contributions for now
!	* Chargino
!	* Neutralino
! ----------------------------------------------------------------

! ----------------------------------------------------------------
!	Full box A^S_RR formfactor
! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_srr_box(i,j,k,l)
!     	Full box A^S_RR formfactor
	INTEGER, INTENT(IN)::i,j,k,l

! 	I haven't put in the switches for now
      	dl_srr_box = (dl_sll_g(i,j,k,l) + dl_srr_h(i,j,k,l) &
&     	     + dl_srr_hg(i,j,k,l))/16/pi/pi
      	return
      	END FUNCTION dl_srr_box

! ================================================================
! 	Formfactor A^S_LR
! ================================================================

	COMPLEX(KIND=DBL) FUNCTION dl_slr_g(i,j,k,l)
!     	Gauge contribution
	INTEGER, INTENT(IN)::i,j,k,l
      	dl_slr_g = (0.d0,0.d0)
      	return
      	END FUNCTION dl_slr_g

! ----------------------------------------------------------------
      
	COMPLEX(KIND=DBL) FUNCTION dl_slr_h(i,j,k,l)
!     	Double Higgs contribution
	INTEGER, INTENT(IN)::i,j,k,l
      	dl_slr_h = (0.d0,0.d0)
      	return
      	END FUNCTION dl_slr_h

! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_slr_hg(i,j,k,l)
!     	W-Higgs contribution
	INTEGER, INTENT(IN)::i,j,k,l
	INTEGER::m,n
      	dl_slr_hg = (0.d0,0.d0)
      	do m=1,3
!	   Loop truncated for SM only, for SUSY n=1,2 (see README)
      	   do n=1,2
!	   Careful with the CKM and such.
      	      dl_slr_hg = dl_slr_hg + CKM(m,i)*DCONJG(yh_eff_r(j,m,n)) &
&     	           * zh(1,n)*dp1(um(m),cm(n),wm,0.d0)
      	   end do
      	end do
      	dl_slr_hg = e2/st2/2*yl(k)*dl_slr_hg
      	return
      	END FUNCTION dl_slr_hg

! ----------------------------------------------------------------
!	I'm skipping the SUSY contributions for now
!	* Chargino
!	* Neutralino
! ----------------------------------------------------------------

! ----------------------------------------------------------------
!	Full box A^S_LR formfactor
! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_slr_box(i,j,k,l)
!     	Full box A^S_LR formfactor
	INTEGER, INTENT(IN)::i,j,k,l

! 	I haven't put in the switches for now
      	dl_slr_box = (dl_slr_g(i,j,k,l) + dl_slr_h(i,j,k,l) &
&     	     + dl_slr_hg(i,j,k,l))/16/pi/pi
      	return
      	END FUNCTION dl_slr_box

! ================================================================
! 	Formfactor A^S_RL
! ================================================================

	COMPLEX(KIND=DBL) FUNCTION dl_srl_g(i,j,k,l)
!     	Gauge contribution
	INTEGER, INTENT(IN)::i,j,k,l
      	dl_srl_g = (0.d0,0.d0)
      	return
      	END FUNCTION dl_srl_g

! ----------------------------------------------------------------
      
	COMPLEX(KIND=DBL) FUNCTION dl_srl_h(i,j,k,l)
!     	Double Higgs contribution
	INTEGER, INTENT(IN)::i,j,k,l
      	dl_srl_h = (0.d0,0.d0)
      	return
      	END FUNCTION dl_srl_h

! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_srl_hg(i,j,k,l)
!     	W-Higgs contribution
	INTEGER, INTENT(IN)::i,j,k,l
	INTEGER::m,n
      	dl_srl_hg = (0.d0,0.d0)
      	do m=1,3
!	   Loop truncated for SM only, for SUSY n=1,2 (see README)
      	   do n=1,2
      	      dl_srl_hg = dl_srl_hg + DCONJG(CKM(m,j))*yh_eff_r(i,m,n)*zh(1,n) &
&     	           * dp1(um(m),cm(n),wm,0.d0)
      	   end do
      	end do
      	dl_srl_hg = e2/st2/2*yl(l)*dl_srl_hg
      	return
      	END FUNCTION dl_srl_hg

! ----------------------------------------------------------------
!	I'm skipping the SUSY contributions for now
!	* Chargino
!	* Neutralino
! ----------------------------------------------------------------

! ----------------------------------------------------------------
!	Full box A^V_RL formfactor
! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_srl_box(i,j,k,l)
!     	Full box A^S_RL formfactor
	INTEGER, INTENT(IN)::i,j,k,l

! 	I haven't put in the switches for now
      	dl_srl_box = (dl_srl_g(i,j,k,l) + dl_srl_h(i,j,k,l) &
&     	     + dl_srl_hg(i,j,k,l))/16/pi/pi
      	return
      	END FUNCTION dl_srl_box

! ================================================================
! 	Formfactor A^T_L
! ================================================================

	COMPLEX(KIND=DBL) FUNCTION dl_tl_g(i,j,k,l)
!     	Gauge contribution
	INTEGER, INTENT(IN)::i,j,k,l

      	dl_tl_g = (0.d0,0.d0)
      	return
      	END FUNCTION dl_tl_g

! ----------------------------------------------------------------
      
	COMPLEX(KIND=DBL) FUNCTION dl_tl_h(i,j,k,l)
!     	Double Higgs contribution
	INTEGER, INTENT(IN)::i,j,k,l

      	dl_tl_h = (0.d0,0.d0)
      	return
      	END FUNCTION dl_tl_h

! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_tl_hg(i,j,k,l)
!     	W-Higgs contribution
	INTEGER, INTENT(IN)::i,j,k,l

      	dl_tl_hg = (0.d0,0.d0)
      	return
      	END FUNCTION dl_tl_hg

! ----------------------------------------------------------------
!	I'm skipping the SUSY contributions for now
!	* Chargino
!	* Neutralino
! ----------------------------------------------------------------

! ----------------------------------------------------------------
!	Full box A^T_L formfactor
! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_tl_box(i,j,k,l)
!     	Full box A^T_L formfactor
	INTEGER, INTENT(IN)::i,j,k,l

! 	I haven't put in the switches for now
      	dl_tl_box = (dl_tl_g(i,j,k,l) + dl_tl_h(i,j,k,l) &
&          + dl_tl_hg(i,j,k,l))/16/pi/pi
      	return
      	END FUNCTION dl_tl_box

! ----------------------------------------------------------------

! ================================================================
! 	Formfactor A^T_R
! ================================================================

	COMPLEX(KIND=DBL) FUNCTION dl_tr_g(i,j,k,l)
!     	Gauge contribution
	INTEGER, INTENT(IN)::i,j,k,l

      	dl_tr_g = (0.d0,0.d0)
      	return
      	END FUNCTION dl_tr_g

! ----------------------------------------------------------------
      
	COMPLEX(KIND=DBL) FUNCTION dl_tr_h(i,j,k,l)
!     	Double Higgs contribution
	INTEGER, INTENT(IN)::i,j,k,l

      	dl_tr_h = (0.d0,0.d0)
      	return
      	END FUNCTION dl_tr_h

! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_tr_hg(i,j,k,l)
!     	W-Higgs contribution
	INTEGER, INTENT(IN)::i,j,k,l

      	dl_tr_hg = (0.d0,0.d0)
      	return
      	END FUNCTION dl_tr_hg

! ----------------------------------------------------------------
!	I'm skipping the SUSY contributions for now
!	* Chargino
!	* Neutralino
! ----------------------------------------------------------------

! ----------------------------------------------------------------
!	Full box A^T_R formfactor
! ----------------------------------------------------------------

	COMPLEX(KIND=DBL) FUNCTION dl_tr_box(i,j,k,l)
!     	Full box A^T_R formfactor
	INTEGER, INTENT(IN)::i,j,k,l

! 	I haven't put in the switches for now
      	dl_tr_box = (dl_tr_g(i,j,k,l) + dl_tr_h(i,j,k,l) &
&          + dl_tr_hg(i,j,k,l))/16/pi/pi
      	return
      	END FUNCTION dl_tr_box

! ----------------------------------------------------------------


END MODULE diagrams_box
