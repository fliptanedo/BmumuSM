! ----------------------------------------------------------------
! Project: Complete One loop Predictions for B->mu+mu- in the MSSM
! Please cite: Dedes, Rosiek, Tanedo (2007)
! ----------------------------------------------------------------
! Main file for B->mu+mu- calculation
! filename: tuning_data.f95
! Current version: 0.01
! By Flip Tanedo (philip.tanedo@gmail.com)
! Date: 19 November 2007
! Latest version available at: http://www.dur.ac.uk/philip.tanedo
! ----------------------------------------------------------------
! This file contains non-physical data. Final results should not
! depend on the particular values in this file. I.e. final results
! should be insensitive to `[fine]-tuning.' 
! ----------------------------------------------------------------
! RECORD OF REVISIONS
!	Date	 Programmer	Description
!	====	 ==========	==================================
!	19/11/07 Flip Tanedo	Pre-Original Code 
!	29/11/07 Flip Tanedo	SM release
! ----------------------------------------------------------------
! DATA INCLUDED
!	NAME		DESCRIPTION
!	==============	==========================================
!	mu_tHooft'	t Hooft mass scale, should cancel 
!	del_diverge	2/(d-4), should cancel (for integrals)
!	eps		"small mass difference" for integrals
! ----------------------------------------------------------------
! PRODUCTION NOTES
! 1. This is the file where you change the numerical values of
! these parameters. However, please note that in the file
! integrals.f95, which is the only file to call these values, 
! mu_tHooft is renamed mu_tH and del_divergence is renamed del.
! This is only so that the funcctions in integrals.f95 can use
! this values without having them listed explicitly in the
! function arguments.
! ----------------------------------------------------------------


MODULE tuning_data
USE program_parameters		! for DBL

!	't Hooft Mass Scale, big; should cancel
	REAL(KIND=DBL), PARAMETER:: mu_tH	= 1.D6
!	Used by divergent integrals in INTEGRALS.F95 and calling subroutines

!	Integral pole, big; should cancel
	REAL(KIND=DBL), PARAMETER:: del		= 1.D6
!	Used by divergent integrals in INTEGRALS.F95 and calling subroutines

!	"Small mass difference" ratio (dimensionless)
	REAL(KIND=DBL), PARAMETER:: eps		= 1.D-20
! 	Used by SORTER.F95


END MODULE tuning_data
