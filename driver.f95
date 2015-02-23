! ----------------------------------------------------------------
! Project: Complete One loop Predictions for B->mu+mu- in the MSSM
! Please cite: Dedes, Rosiek, Tanedo (2007)
! ----------------------------------------------------------------
! Main file for B->mu+mu- calculation
! Current version: 0.03
! By Flip Tanedo (philip.tanedo@gmail.com)
! Date: 29 November 2007
! Latest version available at: n/a
! ----------------------------------------------------------------
! This is the main driver file for the calculation of the 
! branching ratio of $\B^0_{s,d} \rightarrow \ell^+\ell'^-$
! ----------------------------------------------------------------
! RECORD OF REVISIONS
!	Date	 Programmer	Description
!	====	 ==========	==================================
!	17/11/07 Flip Tanedo	Pre-Original Code
!	29/11/07 Flip Tanedo	Original Code, SM only
! ----------------------------------------------------------------
! PRODUCTION NOTES
! 1. I use REAL(KIND=DBL) for double precision. 
!	DBL is defined in program_parameters.f95
! ----------------------------------------------------------------


PROGRAM bmumu
USE program_parameters		! defines DBL for DOUBLE precision
USE physics_data		! Masses, particle properties, etc.
USE tuning_data			! divergences, `small' vales, etc.
USE integrals
!USE sorter			! For integrals, not needed here
USE errorlog			! For error tracking
USE diagrams_box
USE diagrams_zpenguin
USE diagrams_dselfenergy
USE amplitude
USE branchingratio
IMPLICIT NONE  	

REAL(KIND=DBL)::temp1, temp2, temp1u, temp1d, temp2u, temp2d		! Temporary

! DATA DICTIONARY: DECLARATION OF CONSTANTS
! DATA DICTIONARY: LIST OF LOCAL VARIABLES
! PROGRAM CODE GOES BELOW
! MANDATORY SUBROUTINES 
	CALL ckm_init		! Initialise CKM matrix, see physics_data.f95

! I should have some error flag that comes up if ckm_init isn't run.

100	FORMAT ('', 1A20, (3E20.9)) 

	WRITE(*,*)


!	24/11/07 Debugging

!	WRITE(*,*) 'DEL    = ', DEL
!	WRITE(*,*) 'ST     = ', ST
!	WRITE(*,*) 'CT     = ', CT
!	WRITE(*,*) 'ST2    = ', ST2
!	WRITE(*,*) 'CT2    = ', CT2
!	WRITE(*,*) 'SCT    = ', SCT
!	WRITE(*,*) 'SCT2   = ', SCT2
!	WRITE(*,*) 'E      = ', E
!	WRITE(*,*) 'E2     = ', E2
!	WRITE(*,*) 'ALPHA  = ', ALPHA
!	WRITE(*,*) 'WM     = ', WM
!	WRITE(*,*) 'WM2    = ', WM2
!	WRITE(*,*) 'ZM     = ', ZM
!	WRITE(*,*) 'ZM2    = ', ZM2
!	WRITE(*,*) 'PI     = ', PI
!	WRITE(*,*) 'SQ2    = ', SQ2	

!	WRITE(*,*)

!	WRITE(*,100) 'u(m) = ', um(1), um(2), um(3)
!	WRITE(*,100) 'd(m) = ', dm(1), dm(2), dm(3)
!	WRITE(*,100) 'e(m) = ', em(1), em(2), em(3)
!	WRITE(*,100) 'yu(m) = ', yu(1), yu(2), yu(3)
!	WRITE(*,100) 'yd(m) = ', yd(1), yd(2), yd(3)
!	WRITE(*,100) 'cm(m) = ', cm(1), cm(2)


!	WRITE(*,*)
!	WRITE(*,100) 'zh(1,*) = ', zh(1,1), zh(1,2)
!	WRITE(*,100) 'zh(2,*) = ', zh(2,1), zh(2,2)

!	WRITE(*,*)

!	WRITE(*,100) 'cp1(1.1,1.1, 10) = ', cp1(1.1D1,1.1D1, 10.0D0)
!	WRITE(*,100) 'cp1(.1,.1, 100) = ', cp1(0.1D0,0.1D0, 100.0D0)


!	WRITE(*,*)

!	WRITE(*,100)  'yh_eff_l(3,1,1) = ', yh_eff_l(3,1,1)
!	WRITE(*,100)  'yh_eff_l(3,2,1) = ', yh_eff_l(3,2,1)
!	WRITE(*,100)  'yh_eff_l(3,3,1) = ', yh_eff_l(3,3,1)
!	WRITE(*,100)  'yh_eff_l(3,1,2) = ', yh_eff_l(3,1,2)
!	WRITE(*,100)  'yh_eff_l(3,2,2) = ', yh_eff_l(3,2,2)
!	WRITE(*,100)  'yh_eff_l(3,3,2) = ', yh_eff_l(3,3,2)
!	WRITE(*,100)  'yh_eff_r(3,1,1) = ', yh_eff_r(3,1,1)
!	WRITE(*,100)  'yh_eff_r(3,2,1) = ', yh_eff_r(3,2,1)
!	WRITE(*,100)  'yh_eff_r(3,3,1) = ', yh_eff_r(3,3,1)
!	WRITE(*,100)  'yh_eff_r(3,1,2) = ', yh_eff_r(3,1,2)
!	WRITE(*,100)  'yh_eff_r(3,2,2) = ', yh_eff_r(3,2,2)
!	WRITE(*,100)  'yh_eff_r(3,3,2) = ', yh_eff_r(3,3,2)

!	WRITE(*,*)

!	WRITE(*,100)  'yh_eff_l(2,1,1) = ', yh_eff_l(2,1,1)
!	WRITE(*,100)  'yh_eff_l(2,2,1) = ', yh_eff_l(2,2,1)
!	WRITE(*,100)  'yh_eff_l(2,3,1) = ', yh_eff_l(2,3,1)
!	WRITE(*,100)  'yh_eff_l(2,1,2) = ', yh_eff_l(2,1,2)
!	WRITE(*,100)  'yh_eff_l(2,2,2) = ', yh_eff_l(2,2,2)
!	WRITE(*,100)  'yh_eff_l(2,3,2) = ', yh_eff_l(2,3,2)
!	WRITE(*,100)  'yh_eff_r(2,1,1) = ', yh_eff_r(2,1,1)
!	WRITE(*,100)  'yh_eff_r(2,2,1) = ', yh_eff_r(2,2,1)
!	WRITE(*,100)  'yh_eff_r(2,3,1) = ', yh_eff_r(2,3,1)
!	WRITE(*,100)  'yh_eff_r(2,1,2) = ', yh_eff_r(2,1,2)
!	WRITE(*,100)  'yh_eff_r(2,2,2) = ', yh_eff_r(2,2,2)
!	WRITE(*,100)  'yh_eff_r(2,3,2) = ', yh_eff_r(2,3,2)

!	WRITE(*,*)

!	WRITE(*,100) 'CKM(1,1)^d = ', DCONJG(CKM(1,1))
!	WRITE(*,100) 'CKM(1,2)^d = ', DCONJG(CKM(2,1))
!	WRITE(*,100) 'CKM(1,3)^d = ', DCONJG(CKM(3,1))
!	WRITE(*,100) 'CKM(2,1)^d = ', DCONJG(CKM(1,2))
!	WRITE(*,100) 'CKM(2,2)^d = ', DCONJG(CKM(2,2))
!	WRITE(*,100) 'CKM(2,3)^d = ', DCONJG(CKM(3,2))
!	WRITE(*,100) 'CKM(3,1)^d = ', DCONJG(CKM(1,3))
!	WRITE(*,100) 'CKM(3,2)^d = ', DCONJG(CKM(2,3))
!	WRITE(*,100) 'CKM(3,3)^d = ', DCONJG(CKM(3,3))

!	WRITE(*,*)

!	WRITE(*,*)

!	WRITE(*,100) 'dl_vll_g(3,2,2,2) = ', dl_vll_g(3,2,2,2)
!	WRITE(*,100) 'dl_vll_box(3,2,2,2) = ', dl_vll_box(3,2,2,2)
!	WRITE(*,*)
!	WRITE(*,100) 'dl_vrr_h(3,2,2,2) = ', dl_vrr_h(3,2,2,2)
!	WRITE(*,100) 'dl_vrr_box(3,2,2,2) = ', dl_vrr_box(3,2,2,2)
!	WRITE(*,*)
!	WRITE(*,100) 'dl_vlr_h(3,2,2,2) = ', dl_vlr_h(3,2,2,2)
!	WRITE(*,100) 'dl_vlr_box(3,2,2,2) = ', dl_vll_box(3,2,2,2)
!	WRITE(*,100) 'dl_vrl_box(3,2,2,2) = ', dl_vrl_box(3,2,2,2)
!	WRITE(*,100) 'dl_sll_box(3,2,2,2) = ', dl_sll_box(3,2,2,2)
!	WRITE(*,100) 'dl_srr_box(3,2,2,2) = ', dl_srr_box(3,2,2,2)
!	WRITE(*,*)
!	WRITE(*,100) 'dl_slr_hg(3,2,2,2) = ', dl_slr_hg(3,2,2,2)
!	WRITE(*,100) 'dl_slr_box(3,2,2,2) = ', dl_slr_box(3,2,2,2)
!	WRITE(*,*)
!	WRITE(*,100) 'dl_srl_hg(3,2,2,2) = ', dl_srl_hg(3,2,2,2)
!	WRITE(*,100) 'dl_srl_box(3,2,2,2) = ', dl_srl_box(3,2,2,2)
!	WRITE(*,*)
!	WRITE(*,100) 'dl_tl_box(3,2,2,2) = ', dl_tl_box(3,2,2,2)
!	WRITE(*,100) 'dl_tr_box(3,2,2,2) = ', dl_tr_box(3,2,2,2)

!	WRITE(*,*)

!	WRITE(*,100) 'zdd_vl_g(3,2) = ', zdd_vl_g(3,2)
!	WRITE(*,100) 'zdd_vl_hg(3,2) = ', zdd_vl_hg(3,2)
!	WRITE(*,100) 'zdd_vl_h(3,2) = ', zdd_vl_h(3,2)
!	WRITE(*,100) 'zdd_vl(3,2) = ', zdd_vl(3,2)
!	WRITE(*,*)
!	WRITE(*,100) 'zdd_vr_h(3,2) = ', zdd_vr_h(3,2)
!	WRITE(*,100) 'zdd_vr(3,2) = ', zdd_vr(3,2)
!	WRITE(*,*)
!	WRITE(*,100) 'dv_sig0_1(3,2) = ', dv_sig0_1(3,2)
!	WRITE(*,100) 'dv_sig0_4(3,2) = ', dv_sig0_4(3,2)
!	WRITE(*,100) 'dv_sig0(3,2) = ', dv_sig0(3,2)
!	WRITE(*,*)
!	WRITE(*,100) 'dl_vll(3,2,2,2) = ', dl_vll(3,2,2,2)
!	WRITE(*,100) 'dl_vlr(3,2,2,2) = ', dl_vlr(3,2,2,2)
!	WRITE(*,100) 'dl_vrl(3,2,2,2) = ', dl_vrl(3,2,2,2)
!	WRITE(*,100) 'dl_vrr(3,2,2,2) = ', dl_vrr(3,2,2,2)
!	WRITE(*,100) 'dl_sll(3,2,2,2) = ', dl_sll(3,2,2,2)
!	WRITE(*,100) 'dl_slr(3,2,2,2) = ', dl_slr(3,2,2,2)
!	WRITE(*,100) 'dl_srl(3,2,2,2) = ', dl_srl(3,2,2,2)
!	WRITE(*,100) 'dl_srr(3,2,2,2) = ', dl_srr(3,2,2,2)
!	WRITE(*,100) 'dl_tl(3,2,2,2) = ', dl_tl(3,2,2,2)
!	WRITE(*,100) 'dl_tr(3,2,2,2) = ', dl_tr(3,2,2,2)

!	WRITE(*,*)
!	WRITE(*,100) 'fb(2) = ', fb(2)
!	WRITE(*,100) 'MB(2) = ', MB(2)
!	WRITE(*,*)
!	WRITE(*,100) 'fv = ', Fv(3,2,2,2)
!	WRITE(*,100) 'fa = ', Fa(3,2,2,2)
!	WRITE(*,100) 'fs = ', Fs(3,2,2,2)
!	WRITE(*,100) 'fp = ', Fp(3,2,2,2)
!	WRITE(*,*)
!	WRITE(*,100) '|M|^2 = ', amp2(3,2,2,2)
!	WRITE(*,100) 'hbar = ', hbar
!	WRITE(*,100) 'tauB(2) = ', tauB(2)
!	WRITE(*,100) 'gam_b(2) = ', tauB(2)*hbar

	fB(1)	= 0.21565D0
	fB(2)	= 0.23958D0

	WRITE(*,100) 'fBd = ', fB(1)
	WRITE(*,100) 'fBs = ', fB(2) 

	WRITE(*,100) 'BR(s) = ', Br(3,2,2,2)
	WRITE(*,100) 'BR(d) = ', Br(3,1,2,2)

	WRITE(*,*)

	fB(1)	= 0.21635D0
	fB(2)	= 0.24042D0

	WRITE(*,100) 'fBd = ', fB(1)
	WRITE(*,100) 'fBs = ', fB(2) 

	WRITE(*,100) 'BR(s) = ', Br(3,2,2,2)
	WRITE(*,100) 'BR(d) = ', Br(3,1,2,2)

	WRITE(*,*)

	fB(1)	= 0.21565D0
	fB(2)	= 0.23958D0

	temp1d = Br(3,1,2,2)
	temp2d = Br(3,2,2,2)

	fB(1)	= 0.21635D0
	fB(2)	= 0.24042D0

	temp1u = Br(3,1,2,2)
	temp2u = Br(3,2,2,2)

	temp1 = (temp1u + temp1d)/2.D0
	temp2 = (temp2u + temp2d)/2.D0

	WRITE(*,100) 'BR(d) = ', temp1
	WRITE(*,100) '+ or - ', temp1u - temp1
	WRITE(*,100) 'BR(s) = ', temp2
	WRITE(*,100) '+ or - ', temp2u - temp2


	CALL printerrors
	CALL createlogfile('errors.log')
	WRITE(*,*) 'Any errors are listed in errors.log'

END PROGRAM

