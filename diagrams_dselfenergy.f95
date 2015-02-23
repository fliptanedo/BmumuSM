! ----------------------------------------------------------------
! Project: Complete One loop Predictions for B->mu+mu- in the MSSM
! Please cite: Dedes, Rosiek, Tanedo (2007)
! ----------------------------------------------------------------
! Main file for B->mu+mu- calculation
! filename: diagrams_zpenguin.f95
! Current version: 0.01
! By Flip Tanedo (philip.tanedo@gmail.com)
! Date: 22 November 2007
! Latest version available at: http://www.dur.ac.uk/philip.tanedo
! ----------------------------------------------------------------
! This module contains the self energies for the d quark
!
! I would like to rewrite this calculation in the future, but
! for now I'll just port it directly since it's been checked
! several times already. (It will be easy to check new functions
! once I have this file set up.)
!
!	Arguments: 
!	  i and j are the quark flavours
!
!
!    The original code is from d_self0.f (Chankowski, Rosiek 28-1-2000)
!
! The functions here juts account for the actual penguin. The 
! full diagram for B->mu+mu- (i.e. adding Z propagator and
! lepton lines) is calculated in amplitude.f95
!
! ----------------------------------------------------------------
! RECORD OF REVISIONS
!	Date	 Programmer	Description
!	====	 ==========	==================================
!	22/11/07 Flip Tanedo	Original code
! ----------------------------------------------------------------
! USAGE
! USE diagrams_dselfenergy
! dv_sig0(i,j) gives the vectror self energy
! da_sig0(i,j) gives the axial self energy
! amplitudes.f95 calls these and sticks on vertices and lepton parts
! ----------------------------------------------------------------
! PRODUCTION NOTES
! 1. These are lifted from Janusz's code with minimal reformatting.
! 2. Physical constants are as defined in physics_data.f95
! 3. Functions are complex-valued (CKM matrix is complex)
! 4. Recall Janusz's CKM is the dagger of the 'true' CKM
! 5. The naming convention for the contributions and the pieces
!	is weird, but I'm keeping them to stick with Janusz's notation.
! ----------------------------------------------------------------
! DIRECTORY OF FUNCTIONS
! 	dv_sig0(i,j) gives the vectror self energy
! 	da_sig0(i,j) gives the axial self energy
!	ds_sig0(i,j) gives the scalar self energy
!	dp_sig0(i,j) gives the pseudoscalar self energy
!
!	d*_sig0_1    gives the up-W vector contribution
!	d*_sig_07    gives the neutralino contribution (SUSY)
!	d*_sig_08    gives the chargino contribution (SUSY)
!	d*_sig_09    gives the gluino contribution (SUSY)
!
!	... why these strange numbers? I dunno.
! ----------------------------------------------------------------

MODULE diagrams_dselfenergy
USE program_parameters		! for DBL, Pi
USE physics_data		! For coefficients
USE integrals			! for integrals, of course
USE tuning_data			! for mu_tH, 't Hooft mass, and del
!USE sorter			! for array sorting
!USE errorlog			! for error logging
IMPLICIT NONE

CONTAINS

! ================================================================
! 	Vector self-energy (proportional to k(mu)gamma(mu) = G(k)))                                
! ================================================================

      	COMPLEX(KIND=DBL) FUNCTION dv_sig0_1(i,j)
!     	up quark + W in loop
	INTEGER, INTENT(IN)::i,j
	INTEGER::l

!      	dv_sig0_1 = 0
	dv_sig0_1 = (0.D0, 0.D0)
      	do l=1,3
      	   dv_sig0_1 = dv_sig0_1 + e2/2.D0/st2*CKM(l,i)*DCONJG(CKM(l,j)) &
&     	        * b1(0.d0,um(l),wm)
      	end do
!	teMP:
!	dv_sig0_1 = e2/2.D0/st2*CKM(1,i)*DCONJG(CKM(1,j))*bc1(0.d0,um(1),wm)
!	dv_sig0_1 = e2/2.D0/st2*CKM(1,i)*DCONJG(CKM(1,j))*b1(0.d0,um(1),wm)
      	return
      	END FUNCTION dv_sig0_1

! ----------------------------------------------------------------

      	COMPLEX(KIND=DBL) FUNCTION dv_sig0_4(i,j)
!     	up quark + charged Higgs (Goldstone) in loop
	INTEGER, INTENT(IN)::i,j
	INTEGER::l,k

      	dv_sig0_4 = 0
      	do l=1,3
!	   Loop truncated for SM only, for SUSY k=1,2 (see README)
      	   do k=1,2
      	      dv_sig0_4 = dv_sig0_4 + ((zh(2,k)*yu(l))**2 &
&     	           + yd(i)*yd(j)*zh(1,k)**2)/2 &
&     	           * CKM(l,i)*DCONJG(CKM(l,j))*b1(0.d0,um(l),cm(k))
      	   end do
      	end do
      	return
      	END FUNCTION dv_sig0_4

! ----------------------------------------------------------------
!	I'm skipping the SUSY contributions for now
!	* Neutralino + down squark
!	* Chargino + up squark
!	* Gluino + down squark
! ----------------------------------------------------------------

! ----------------------------------------------------------------
!	Full bare down quark self-denergy, vector part
! ----------------------------------------------------------------

      	COMPLEX(KIND=DBL) FUNCTION dv_sig0(i,j)
!     	Full bare down quark self-denergy, vector part
	INTEGER, INTENT(IN)::i,j
	
! 	I haven't put in the switches for now
      	dv_sig0 = (0.d0,0.d0)
      	dv_sig0 = dv_sig0 + dv_sig0_1(i,j) + dv_sig0_4(i,j)
      	dv_sig0 = dv_sig0/16/pi/pi
      	return
      	END FUNCTION dv_sig0


! ================================================================
! 	Axial self-energy (proportional to k(mu)gamma(mu)gamma(5)
! ================================================================

      	COMPLEX(KIND=DBL) FUNCTION da_sig0_1(i,j)
!     	up quark + W in loop
	INTEGER, INTENT(IN)::i,j
	INTEGER::l

      	da_sig0_1 = 0
      	do l=1,3
      	   da_sig0_1 = da_sig0_1 - e2/2/st2*CKM(l,i)*DCONJG(ckm(l,j)) &
&     	        * b1(0.d0,um(l),wm)
      	end do
      	return
      	END FUNCTION da_sig0_1

! ----------------------------------------------------------------

      	COMPLEX(KIND=DBL) FUNCTION da_sig0_4(i,j)
!     	up quark + charged Higgs in loop
	INTEGER, INTENT(IN)::i,j
	INTEGER::l,k

      	da_sig0_4 = 0
      	do l=1,3
!	   Loop truncated for SM only, for SUSY k=1,2 (see README)
      	   do k=1,2
      	      da_sig0_4 = da_sig0_4 - ((zh(2,k)*yu(l))**2 &
&     	           - yd(i)*yd(j)*zh(1,k)**2)/2 &
&     	           * CKM(l,i)*DCONJG(CKM(l,j))*b1(0.d0,um(l),cm(k))
      	   end do
      	end do
      	return
      	END FUNCTION da_sig0_4

! ----------------------------------------------------------------
!	I'm skipping the SUSY contributions for now
!	* Neutralino + down squark
!	* Chargino + up squark
!	* Gluino + down squark
! ----------------------------------------------------------------

! ----------------------------------------------------------------
!	Full bare down quark self-denergy, axial part
! ----------------------------------------------------------------

      	COMPLEX(KIND=DBL) FUNCTION da_sig0(i,j)
!     	Full bare down quark self-denergy, vector part
	INTEGER, INTENT(IN)::i,j

      	da_sig0 = (0.d0,0.d0)

! 	I haven't put in the switches for now
      	da_sig0 = da_sig0 + da_sig0_1(i,j) + da_sig0_4(i,j) 
      	da_sig0 = da_sig0/16/pi/pi
      	return
      	END FUNCTION da_sig0

! ================================================================
! 	Scalar self-energy 
! ================================================================

      	COMPLEX(KIND=DBL) FUNCTION ds_sig0_1(i,j)
!     	up quark + W in loop
	INTEGER, INTENT(IN)::i,j

	ds_sig0_1 = (0.D0, 0.D0)
	return
	END FUNCTION ds_sig0_1

! ----------------------------------------------------------------


      	COMPLEX(KIND=DBL) FUNCTION ds_sig0_4(i,j)
!     	up quark + charged Higgs in loop
	INTEGER, INTENT(IN)::i,j
	INTEGER::l,k

      	ds_sig0_4 = 0
      	do l=1,3
!	   Loop truncated for SM only, for SUSY k=1,2 (see README)
      	   do k=1,2
      	      ds_sig0_4 = ds_sig0_4 + zh(2,k)*zh(1,k)*um(l)*yu(l)/2 &
&     	           * (yd(i) + yd(j)) &
&     	           * CKM(l,i)*DCONJG(CKM(l,j))*b0(0.d0,um(l),cm(k))
      	   end do
      	end do
      	return
      	END FUNCTION ds_sig0_4


! ----------------------------------------------------------------
!	I'm skipping the SUSY contributions for now
!	* Neutralino + down squark
!	* Chargino + up squark
!	* Gluino + down squark
! ----------------------------------------------------------------

! ----------------------------------------------------------------
!	Full bare down quark self-denergy, scalar part
! ----------------------------------------------------------------

      	COMPLEX(KIND=DBL) FUNCTION ds_sig0(i,j)
!     	Full bare down quark self-denergy, scalar part
	INTEGER, INTENT(IN)::i,j

      	ds_sig0 = (0.d0,0.d0)

! 	I haven't put in the switches for now
      	ds_sig0 = ds_sig0 + ds_sig0_1(i,j) + ds_sig0_4(i,j) 
      	ds_sig0 = ds_sig0/16/pi/pi
      	return
      	END FUNCTION ds_sig0

! ================================================================
! 	Pseudoscalar self-energy (proportional to G(5))
! ================================================================

      	COMPLEX(KIND=DBL) FUNCTION dp_sig0_1(i,j)
!     	up quark + W in loop
	INTEGER, INTENT(IN)::i,j

	dp_sig0_1 = (0.D0, 0.D0)
	return
	END FUNCTION dp_sig0_1

! ----------------------------------------------------------------

      	COMPLEX(KIND=DBL) FUNCTION dp_sig0_4(i,j)
!     	up quark + charged Higgs in loop
	INTEGER, INTENT(IN)::i,j
	INTEGER::l,k

      	dp_sig0_4 = 0
      	do l=1,3
!	   Loop truncated for SM only, for SUSY k=1,2 (see README)
      	   do k=1,2
              dp_sig0_4 = dp_sig0_4 + zh(2,k)*zh(1,k)*um(l)*yu(l) &
&                  * (yd(i) - yd(j))/2 &
&                  * CKM(l,i)*DCONJG(CKM(l,j))*b0(0.d0,um(l),cm(k))
      	   end do
      	end do
      	return
      	END FUNCTION dp_sig0_4


! ----------------------------------------------------------------
!	I'm skipping the SUSY contributions for now
!	* Neutralino + down squark
!	* Chargino + up squark
!	* Gluino + down squark
! ----------------------------------------------------------------

! ----------------------------------------------------------------
!	Full bare down quark self-denergy, pseudoscalar part
! ----------------------------------------------------------------

      	COMPLEX(KIND=DBL) FUNCTION dp_sig0(i,j)
!     	Full bare down quark self-denergy, scalar part
	INTEGER, INTENT(IN)::i,j

      	dp_sig0 = (0.d0,0.d0)

! 	I haven't put in the switches for now
      	dp_sig0 = dp_sig0 + dp_sig0_1(i,j) + dp_sig0_4(i,j) 
      	dp_sig0 = dp_sig0/16/pi/pi
      	return
      	END FUNCTION dp_sig0

! ----------------------------------------------------------------

	
END MODULE diagrams_dselfenergy
