! ----------------------------------------------------------------
! Project: Complete One loop Predictions for B->mu+mu- in the MSSM
! Please cite: Dedes, Rosiek, Tanedo (2007)
! ----------------------------------------------------------------
! Main file for B->mu+mu- calculation
! filename: program_parameters.f95
! Current version: 0.02
! By Flip Tanedo (philip.tanedo@gmail.com)
! Date: 16 November 2007
! Latest version available at: http://www.dur.ac.uk/philip.tanedo
! ----------------------------------------------------------------
! This file contains non-physics data for running the program in
! a processor-independent way.
! ----------------------------------------------------------------
! RECORD OF REVISIONS
!	Date	 Programmer	Description
!	====	 ==========	==================================
!	18/11/07 Flip Tanedo	Pre-Original Code 
!	20/11/07 Flip Tanedo	Added PI (used in integrals.f95)
!	20/11/07 Flip Tanedo	Added sq2 (used in diagrams)
! ----------------------------------------------------------------
! PRODUCTION NOTES
! 1. Accuracy is Pi is to match Janusz's code
! ----------------------------------------------------------------

MODULE program_parameters
IMPLICIT NONE

INTEGER, PARAMETER :: DBL = KIND(1.0D0)
! Defines the the kind DBL, used for double precision
! Example usage: 
!	REAL(kind=DBL) :: a = 1.

REAL(KIND=DBL), PARAMETER :: PI = 3.1415926536D0
REAL(KIND=DBL), PARAMETER :: sq2 = SQRT(REAL(2.D0,DBL))
COMPLEX(KIND=DBL), PARAMETER :: img = (0.D0, 1.D0)

END MODULE program_parameters
