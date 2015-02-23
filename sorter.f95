! ----------------------------------------------------------------
! Project: Complete One loop Predictions for B->mu+mu- in the MSSM
! Please cite: Dedes, Rosiek, Tanedo (2007)
! ----------------------------------------------------------------
! Main file for B->mu+mu- calculation
! filename: sorter.f95
! Current version: 0.01
! By Flip Tanedo (philip.tanedo@gmail.com)
! Date: 19 November 2007
! Latest version available at: http://www.dur.ac.uk/philip.tanedo
! ----------------------------------------------------------------
! This module contains a basic selection sort subroutine for
! sorting a REAL(KIND=DBL) array x of dimension n.
! The final sort ordering is such that x(1) <= x(2) <= ... <= x(n) 
! Additional feature: Following Janusz Rosiek's code, I'm also
! building into soft an automatic 'smoother' for masses that are
! close to one another. I.e. Integrals with arguments (m, m+eps, ...)
! end up having an awkward dependence on this small mass difference.
!
! A subroutine averager(avgs,food,kill) is also included to make this second
! feature work properly.
! 
! This is used by integrals.f95.
!
! ----------------------------------------------------------------
! RECORD OF REVISIONS
!	Date	 Programmer	Description
!	====	 ==========	==================================
!	19/11/07 Flip Tanedo	Original Code
! ----------------------------------------------------------------
! USAGE
! CALL sort_arr(x,n)
! ----------------------------------------------------------------

MODULE sorter
USE program_parameters		! for DBL
USE tuning_data			! for eps (small mass difference ratio)
IMPLICIT NONE

CONTAINS

!	Sort an n-dim REAL(KIND=DBL) array from low to high
!	Also smooth out small mass differences of order eps
	SUBROUTINE sort_arr(x,n)

!		PSEUDOCODE:
!		Loop(i=1,n)
!		   Loop(j=i,n)
!			IF(x(i)>x(j)), swap their values
!		   END loop
!		END loop

		IMPLICIT NONE
		INTEGER, INTENT(IN)::n
		REAL(kind=DBL), INTENT(INOUT)::x(n)
		INTEGER::i,j,k
		REAL(kind=DBL)::temp
		REAL(kind=DBL)::temp2=1.0

		DO i=1,n
		   DO j=i,n
			IF(x(i)>x(j)) THEN
			   temp=x(i)
			   x(i)=x(j)
			   x(j)=temp
			END IF
		   END DO
		END DO	

!		CONTINUE TO MASS DIFFERENCE CHECKER!

!		PSEUDOCODE:
!		i=1
!		Kill any old average (reset)
!		DO WHILE(i<n)
!		   Define counter j=i+1
!		   Call running average
!		   DO WHILE([|average-x(j)|/max(|average|,|x(j)|)< eps].AND.[j<=n])
!			Feed average, keep going
!			j = j+1
!		   END DO
!		   Loop(k=i,j-1)
!			x(k) = average
!		   END LOOP
!		   i = j
!		END DO

		i=1				! reset counter	
		temp2=0.			! food = 0
		DO WHILE(i<n)			! nothing to average when i=n
		   CALL averager(temp,temp2,1)	! reset the averager
		   j=i+1
		   CALL averager(temp,x(i),0)	! Take first value
			DO WHILE((ABS(temp-x(j))/MAX(ABS(temp),ABS(x(j)))<eps).AND.(j<=n))
			CALL averager(temp,x(j),0)
			j=j+1
			END DO

			DO k=i,j-1
			x(k) = temp
			END DO
		   i=j
		END DO

	END SUBROUTINE sort_arr

! ----------------------------------------------------------------

!	AVERAGER subroutine
	SUBROUTINE averager(avg,food,kill)

	REAL(KIND=DBL), INTENT(OUT)::avg	! average
	REAL(KIND=DBL), INTENT(IN)::food	! new value
	INTEGER, INTENT(IN)::kill		! 1 means reset, 0 means go on
	INTEGER, SAVE::n			! number of data points
	REAL(KIND=DBL),SAVE::thesum		! running sum

	IF(kill==1) THEN
		thesum	= 0
		n 	= 0
!		WRITE(*,*) 'averager killed'

	ELSE
		n = n+1
		thesum = thesum + food
		avg = thesum/n
!		WRITE(*,*) 'n = ', n, ', avg = ', avg

	END IF

	END SUBROUTINE averager	

! ----------------------------------------------------------------

END MODULE sorter


! ----------------------------------------------------------------

!	JANUSZ'S CODE

!      subroutine sort_arr(x,n)
!      implicit double precision (a-h,o-z)
!      dimension x(n)
!      common/cp_acc/eps
!      external init_cd_fun
!      if (n.lt.2) return
!c     Sort array. Final ordering x(1) <= x(2) <= ... <= x(n)
!      do i=1,n-1
!        do j=1,n-i
!          if (x(j).gt.x(j+1)) then
!             tmp     =  x(j+1)
!             x(j+1)  =  x(j)
!             x(j)    =  tmp
!           end if
!        end do
!      end do
!c     Replace very close masses by the equal ones, i.e.
!c     x(i) = x(i+1) = ... = x(i+n) = (x(1) + ... x(n))/n
!      i = 1
! 10   if (abs(x(i+1) - x(i))/max(abs(x(i+1)),abs(x(i))).le.eps) then
!         avg = x(i)
!         j = 1
! 20      avg = (j*avg + x(i+j))/(j+1)
!         j = j + 1 
!         if ((i+j).gt.n) goto 30 
!         if (abs(x(i+j) - avg)/max(abs(x(i+j)),abs(avg)).le.eps) goto 20
! 30      j = j - 1
!         do k=0,j
!            x(i+k) = avg
!         end do
!         i = i + j
!      else
!         i = i + 1
!      end if      
!      if (i.ge.n) return
!      goto 10 
!      end
