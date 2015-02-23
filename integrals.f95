! ----------------------------------------------------------------
! Project: Complete One loop Predictions for B->mu+mu- in the MSSM
! Please cite: Dedes, Rosiek, Tanedo (2007)
! ----------------------------------------------------------------
! Main file for B->mu+mu- calculation
! filename: integrals.f95
! Current version: 1.00
! By Flip Tanedo (philip.tanedo@gmail.com)
! Date: 16 November 2007
! Latest version available at: http://www.dur.ac.uk/philip.tanedo
! ----------------------------------------------------------------
! This module contains the Passarino-Veltman integrals relevant
! to B->mu+mu-. Expressions were taken from Misiak and Rosiek.
! See, for example, Axelrod Nucl.Phys.B209(1982)p.349. 
! NOTE: Following Rosiek et al., we use a different sign for the
! b functions. (Sorry, I know this is kind of annoying.)
! 	Comparison to JR's FORTRAN 77 code:
!		Box, Z integrals: 	cd_fun.f
!		Self-energy:		b_fun.f
! ----------------------------------------------------------------
! RECORD OF REVISIONS
!	Date	 Programmer	Description
!	====	 ==========	==================================
!	18/11/07 Flip Tanedo	Pre-Original Code 
!	19/11/07 Flip Tanedo	Box and Z-penguin integrals
!	20/11/07 Flip Tanedo	Self-energy integrals
!	22/11/07 Flip Tanedo	No explicit divergence arguments
! ----------------------------------------------------------------
! USAGE
! For example:
!	dp1(mass1, mass2, mass3, mass4)
! Return value is a REAL(kind=DBL); see production note #2
!
! USE integrals
! Call functions as above.
! ----------------------------------------------------------------
! PRODUCTION NOTES
! 1. These are lifted from Janusz's code with minimal reformatting.
!	They have been checked by Flip on 20/11/07. (not extensively)
!	Not all of Janusz's functions are included here, just the
!	ones that we need for this program.
! 2. Janusz's b functions were originally defined to be complex,
!	but I believe the sector that we care about is manifestly
!	real, so I've left the b-functions as REAL(KIND=DBL).
! 3. Originally I had the divergent integrals and convergent integrals
!	behave take different arguments (del and mu_tH). However,
!	I decided to change this and define del and mu_tH within
!	the main part of the module. Thus all integrals take the
!	same arguments as they did in Janusz's code.
! ----------------------------------------------------------------
! DIRECTORY OF FUNCTIONS
! 1. Box diagrams: 	dp0, dp1
! 2. Z-penguins:	cp0, cp1
! 3. Self energies:	b0,  b1 
! ----------------------------------------------------------------


MODULE integrals
USE program_parameters		! for DBL, Pi
USE tuning_data			! for mu_tH, 't Hooft mass, and del
USE sorter			! for array sorting
USE errorlog			! for error logging
IMPLICIT NONE

! Do not change these. Instead, change the values in tuning_data.f95
! The reason that I've put this (trivial) statement here is so that
! the functions in this module can call them without having them
! explicitly listed in their arguments. 
! (This might be an extraneous step, but it compartamentalises things)
! REAL(KIND=DBL), PARAMETER :: del = del_diverge
! REAL(KIND=DBL), PARAMETER :: mu_tH = mu_tHooft
! 22/11/07: I removed this and just left the parameters in tuning_data.f95

CONTAINS

! ----------------------------------------------------------------
!	Four point integral
!		/   d^4k	        1
!	dp0 = 	|  ------- -------------------------
!		/  (2pi)^4 (k^2-am^2) .... (k^2-dm^2)
! ----------------------------------------------------------------

      	REAL(KIND=DBL) FUNCTION dp0(am,bm,cm,dm)
!     	Four point scalar function (1 in the numerator)
	REAL(KIND=DBL), INTENT(IN)::am,bm,cm,dm
	REAL(KIND=DBL)::x(4)

	REAL(KIND=DBL)::a,b,c,d

      	x(1) = am*am
      	x(2) = bm*bm
      	x(3) = cm*cm
      	x(4) = dm*dm

      	call sort_arr(x,4)

      	a = x(1)
      	b = x(2)
      	c = x(3)
      	d = x(4)

!	---------------------
      	if (a.eq.0) then
!     	Mass combination 00cd
      	  if (b.eq.0) then
		call logerror('ERROR Int.dp0.1: 2 or more masses=0 in dp0')
		dp0=0
!		The error message below is redundant
!		stop 'ERROR Int.dp0.1: 2 or more masses=0 in dp0'
		return
	  end if

      	  if (b.eq.d) then
!     	Mass combination 0bbb
      	    dp0 = 1.d0/2/b/b
      	  else if (b.eq.c) then
!     	Mass combination 0bbd
      	    dp0 = (1/b - log(d/b)/(d - b))/(d - b)
      	  else if (c.eq.d) then
!     	Mass combination 0bcc
      	    dp0 = (1/c - log(b/c)/(b - c))/(b - c)
      	  else
!     	Mass combination 0bcd
      	    dp0 = (log(c/b)/(c - b) - log(d/b)/(d - b))/(d - c)
      	  end if
      	return
      	end if

!	---------------------
      	if (a.eq.d) then
!     	Mass combination aaaa
      	  dp0 = 1.d0/6/a/a
      	else if (a.eq.c) then
!     	Mass combination aaad
      	  dp0 = ((1 + d/a)/2 - log(d/a)/(1 - a/d))/(d - a)/(d - a)
      	else if (a.eq.b) then
      	  if (c.eq.d) then
!     	Mass combination aacc
      	    dp0 = - (2 + (a + c)/(a - c)*log(c/a))/(c - a)/(c - a)
      	  else
!     	Mass combination aacd
      	    dp0 = - 1/(d - a)/(c - a) - log(c/a)/(1 - d/c)/(c - a)/(c - a) &
&     	                              - log(d/a)/(1 - c/d)/(d - a)/(d - a)
      	  end if
      	else
      	  if (b.eq.d) then
!     	Mass combination abbb
      	    dp0 = ((1 + a/b)/2 - log(a/b)/(1 - b/a))/(a - b)/(a - b)
      	  else if (b.eq.c) then
!     	Mass combination abbd
      	    dp0 = - 1/(d - b)/(a - b) - log(a/b)/(1 - d/a)/(a - b)/(a - b) &
&     	                              - log(d/b)/(1 - a/d)/(d - b)/(d - b)
      	  else if (c.eq.d) then
!     	Mass combination abcc
      	    dp0 = - 1/(b - c)/(a - c) - log(a/c)/(1 - b/a)/(a - c)/(a - c) &
&     	                              - log(b/c)/(1 - a/b)/(b - c)/(b - c)
      	  else 
!     	Mass combination abcd
      	    dp0 = - log(b/a)/(1 - a/b)/(b - c)/(b - d) &
&     	          - log(c/a)/(1 - a/c)/(c - b)/(c - d) &
&     	          - log(d/a)/(1 - a/d)/(d - b)/(d - c)
      	  end if
      	end if
      	return
      	END FUNCTION dp0

! ----------------------------------------------------------------
!	Four Point Integral
!		/   d^4k	      k^2
!	dp1 = 	|  ------- -------------------------
!		/  (2pi)^4 (k^2-am^2) .... (k^2-dm^2)
! ----------------------------------------------------------------

	REAL(KIND=DBL) FUNCTION dp1(am,bm,cm,dm)
!     	Four point scalar function (k^2 in the numerator)
	REAL(KIND=DBL), INTENT(IN)::am,bm,cm,dm
	REAL(KIND=DBL)::x(4)
	REAL(KIND=DBL)::a,b,c,d

      	x(1) = am*am
      	x(2) = bm*bm
      	x(3) = cm*cm
      	x(4) = dm*dm
	
      	call sort_arr(x,4)
      
	a = x(1)
      	b = x(2)
      	c = x(3)
      	d = x(4)

!	--------------------------------

      	if (a.eq.0) then
!     	Mass combination 000d
      	  if (c.eq.0) then
		call logerror('ERROR Int.dp1.1: 3 or more masses = 0 in dp1')
		dp1=0
!		stop 'ERROR Int.dp1.1: 3 or more masses = 0 in dp1'
		return
	  end if
      	  if (b.eq.0) then
      	    if (c.eq.d) then
!     	Mass combination 00cc
      	      dp1 = - 1/c
      	    else
!     	Mass combination 00cd
      	      dp1 = - log(d/c)/(d - c)
      	    end if
      	  else if (b.eq.d) then
!     	Mass combination 0bbb
      	    dp1 = - 1.d0/2/b
      	  else if (b.eq.c) then
!     	Mass combination 0bbd
      	    dp1 = (1 - log(d/b)/(1 - b/d))/(d - b)
      	  else if (c.eq.d) then
!     	Mass combination 0bcc
      	    dp1 = (1 - log(b/c)/(1 - c/b))/(b - c)
      	  else
!     	Mass combination 0bcd
      	    dp1 = (log(b/d)/(1 - d/b) - log(c/d)/(1 - d/c))/(c - b)
      	  end if
      	return
      	end if

!	--------------------------------

      	if (a.eq.d) then
!     	Mass combination aaaa
      	  dp1 = - 1.d0/3/a
      	else if (a.eq.c) then
!     	Mass combination aaad
      	  dp1 = ((3*d - a)/2 - d*log(d/a)/(1 - a/d))/(d - a)/(d - a)
      	else if (a.eq.b) then
      	  if (c.eq.d) then
!     	Mass combination aacc
      	    dp1 = - (a + c + 2*c*log(c/a)/(1 - c/a))/(c - a)/(c - a)
      	  else
!     	Mass combination aacd
      	    dp1 = - a/(d - a)/(c - a) &
&    	          - c*log(c/a)/(1 - d/c)/(c - a)/(c - a) &
&     	          - d*log(d/a)/(1 - c/d)/(d - a)/(d - a)
      	  end if
      	else
      	  if (b.eq.d) then
!     	Mass combination abbb
      	    dp1 = ((3*a - b)/2 - a*log(a/b)/(1 - b/a))/(a - b)/(a - b)
      	  else if (b.eq.c) then
!     	Mass combination abbd
      	    dp1 = - b/(d - b)/(a - b) &
&     	          - a*log(a/b)/(1 - d/a)/(a - b)/(a - b) &
&     	          - d*log(d/b)/(1 - a/d)/(d - b)/(d - b)
      	  else if (c.eq.d) then
!     	Mass combination abcc
      	    dp1 = - c/(b - c)/(a - c) &
&     	          - a*log(a/c)/(1 - b/a)/(a - c)/(a - c) &
&     	          - b*log(b/c)/(1 - a/b)/(b - c)/(b - c)
      	  else
!     	Mass combination abcd
      	    dp1 = - b*log(b/a)/(1 - a/b)/(b - c)/(b - d) &
&     	          - c*log(c/a)/(1 - a/c)/(c - b)/(c - d) &
&     	          - d*log(d/a)/(1 - a/d)/(d - b)/(d - c)
      	  end if
      	end if
      	return
      	END FUNCTION dp1


! ----------------------------------------------------------------
!	Three Point Functions
!		/   d^4k	        1
!	cp0 = 	|  ------- -------------------------
!		/  (2pi)^4 (k^2-am^2) .... (k^2-cm^2)
! ----------------------------------------------------------------

      	REAL(KIND=DBL) FUNCTION cp0(am,bm,cm)
!     	Three point scalar function (1 in the numerator)

	REAL(KIND=DBL),INTENT(IN)::am,bm,cm
	REAL(KIND=DBL)::x(3)
	REAL(KIND=DBL)::a,b,c
	
      	x(1) = am*am
      	x(2) = bm*bm
      	x(3) = cm*cm

      	call sort_arr(x,3)

      	a = x(1)
      	b = x(2)
      	c = x(3)

!	-------------------------------------------

      	if (a.eq.0) then
      	  if (b.eq.0) then
		call logerror('ERROR Int.cp0.1: 2 or more masses = 0 in cp0')
		cp0=0
!		stop 'ERROR Int.cp0.1: 2 or more masses = 0 in cp0'
		return
	  end if
      	  if (b.eq.c) then
      	    cp0 = - 1/b
      	  else 
      	    cp0 = log(b/c)/(c - b)
      	  end if
      	  return
      	end if
      	if (a.eq.c) then
      	  cp0 = - 1.d0/2/a
      	else if (a.eq.b) then
      	  cp0 =  (1  - log(c/a)/(1 - a/c))/(c - a)
      	else if (b.eq.c) then
      	  cp0 =  (1  - log(a/c)/(1 - c/a))/(a - c)
      	else
      	  cp0 = - (log(b/a)/(1 - a/b) - log(c/a)/(1 - a/c))/(b - c)
      	end if
      	return
      	END FUNCTION cp0


! ----------------------------------------------------------------
!	Three Point Integral (NOTE ARGUMENTS!)
!	This integral diverges, include del and 't Hooft mass
!
!		/   d^4k	      k^2
!	cp1 = 	|  ------- -------------------------
!		/  (2pi)^4 (k^2-am^2) .... (k^2-cm^2)
! ----------------------------------------------------------------

      	REAL(KIND=DBL) FUNCTION cp1(am,bm,cm)

!      	REAL(KIND=DBL) FUNCTION cp1(am,bm,cm, del,mu_tH)
!     	Three point scalar function (k^2 in the numerator)
!      	common/renorm/del,amiu2,infstat
!      	external init_spence
	REAL(KIND=DBL), INTENT(IN)::am,bm,cm !,del,mu_tH
	REAL(KIND=DBL)::x(3)
	REAL(KIND=DBL)::a,b,c

      	x(1) = am*am
      	x(2) = bm*bm
      	x(3) = cm*cm

      	call sort_arr(x,3)

      	a = x(1)
      	b = x(2)
      	c = x(3)

! 	--------------------------------------------

      	if (a.eq.0) then
      	  if (c.eq.0) then
		call logerror('ERROR: Int.cp1.1: 3 x 0 masses in cp1')
		cp1 = 0
!		stop 'ERROR: Int.cp1.1: 3 x 0 masses in cp1'
		return
	  end if
      	  cp1 = del + log(mu_tH/c)
      	  if (b.eq.0) then
      	    cp1 = cp1 + 1
      	    return
      	  end if
      	  if (b.ne.c) cp1 = cp1 + 1 + log(c/b)/(1 - c/b)
      	  return
      	end if
      	if (a.eq.c) then
      	  cp1 = del + log(mu_tH/c) - 0.5d0
      	else if (a.eq.b) then
      	  cp1 = del + log(mu_tH/a) + (1 - log(c/a)/(1 - a/c))/(1 - a/c)
      	else if (b.eq.c) then
      	  cp1 = del + log(mu_tH/c) + (1 - log(a/c)/(1 - c/a))/(1 - c/a)
      	else
      	  cp1 = del + log(mu_tH/a) + 1 - log(b/a)/(a/b - 1)/(c/b - 1) &
&     	      - log(c/a)/(a/c - 1)/(b/c - 1)
      	end if
      	return
      	END FUNCTION cp1


! ----------------------------------------------------------------
!	SELF ENERGY INTEGRAL
!	This integral diverges, include del and 't Hooft mass
!
!		/   d^4k	        1
!	b0 = 	|  ------- -------------------------
!		/  (2pi)^4 ((k-p)^2-a1^2) (k^2-a2^2)
!
!	Depends on a bunch of switches and auxilliary functions.
!	For now I'll just put up an error message any time one
!	of these functions is called. I don't think they're necessary
!	but I'll add them in if the error messages come up
!	during these calculations.
!
!	In Janusz's code, b0 and b1 are complex. Why?
! ----------------------------------------------------------------

      	REAL(KIND=DBL) FUNCTION b0(s,a1,a2)
!      	REAL(KIND=DBL) FUNCTION b0(s,a1,a2,del,mu_tH)
!     	B0 Veltman function
	REAL(KIND=DBL), INTENT(IN)::s,a1,a2 !,del,mu_tH
	REAL(KIND=DBL)::summ,dif
!	COMMENT: sum -> summ from Janusz's code, sum is a keyword

!	COMMENTED OUT: what the hell is dbstat?
!      	if (dbstat) then
!      	  b0 = db0(s,a1,a2)
!      	  return
!      	end if

!	COMMENTED OUT: what the hell is infstat?
!      	if (infstat) then
!      	  b0 = del
!      	  return
!      	end if


!	This is the case that we care about
      	if (s.eq.0.d0) then

      	  if (abs(a1) + abs(a2).eq.0.d0) then
      	    b0 = del
      	    return
      	  end if

      	  if (a1*a2.eq.0.d0) then
      	    b0 = del - log((a1 + a2)*(a1 + a2)/mu_tH) + 1
      	  else
      	    b0 = del - log(a1*a2/mu_tH)
      	    if (a1.ne.a2) then
      	      summ = a1*a1 + a2*a2
      	      dif = a1*a1 - a2*a2
     	      b0 = b0 + 1 - summ/dif*log(a1/a2)
      	    end if
      	  end if

      	  return

!	ERROR HANDLING if first argument isn't 0
	else
		CALL logerror('ERROR Int.b0.1: bad parameters, see README')
		b0 = 0
		return

      	end if

!	COMMENT: The code below was commented out because it shouldn't be
!	necessary and was producing weird recursive errors. I replaced it
!	with the 'else' statement above.

!!	COMMENT: Error flag to say that we're accessing an unexpected part
!!	of this function. We shouldn't need this for Bmumu.
!!	WRITE(*,*) 'ERROR Int.b0.1: see README.txt'

!!	If (s!=0), not a case we care about in Bmumu
!      	if ((a1*a2).eq.0.d0) then
!      	  if (a1.eq.a2) then
!      	    b0 = del - log(abs(s)/mu_tH) + 2
!!      	    if (s.gt.0.d0) b0 = b0 + (0,1)*pi
!!	    COMMENT: the above line adds a complex part!
!!	    if (s.gt.0.d0) WRITE(*,*) 'ERROR Int.b0.2: Imaginary part of integral.'
!      	  else
!!      	    b0 = del - log((a1 + a2)*(a1 + a2)/mu_tH) + 1 + f(s,a1,a2)
!!	    COMMENT: I don't think we need the above line. Error otherwise.
!!	    WRITE(*,*) 'ERROR Int.b0.3: Trying to access f functions'
!      	  end if
!      	  return
!      	end if

!!      	b0 = del - log(a1*a2/mu_tH) + f(s,a1,a2)
!!	COMMENT: I don't think we need the above line. Error otherwise.
!!	WRITE(*,*) 'ERROR Int.b0.3: Trying to access f functions'
	

!!	COMMENT: I'm commenting out the block belowk see error above
!!      	if (a1.ne.a2) then
!!      	  sum = a1*a1 + a2*a2
!!      	  dif = a1*a1 - a2*a2
!!      	  b0 = b0 + 1 - sum/dif*log(a1/a2)
!!      	end if

      	return
      	END FUNCTION b0
	

! ----------------------------------------------------------------
!	SELF ENERGY INTEGRAL
!	This integral diverges, include del and 't Hooft mass
!
!		/   d^4k	       p*k
!	b1 = 	|  ------- -------------------------
!		/  (2pi)^4 ((k-p)^2-a1^2) (k^2-a2^2)
!
!	Depends on function f, see below
! ----------------------------------------------------------------


      	REAL(KIND=DBL) FUNCTION b1(s,a1,a2)
!      	REAL(KIND=DBL) FUNCTION b1(s,a1,a2,del,mu_tH)
!     	B1 Veltman function
	REAL(KIND=DBL), INTENT(IN)::s,a1,a2 !,del,mu_tH
	REAL(KIND=DBL)::dif

!	COMMENT: I have no idea what dbstat is for
!      	if (dbstat) then
!      	  b1 = db1(s,a1,a2)
!      	  return
!      	end if

!	COMMENT: I have no idea what infstat is for
!      	if (infstat) then
!      	  b1 = - del/2
!      	  return
!      	end if

!	COMMENT: This is the case we care about
      	if (s.eq.0.d0) then

      	  if (abs(a1) + abs(a2).eq.0.d0) then
		call logerror('ERROR Int.b1.2: b1 called for s=m1=m2=0')
		b1=0
!		stop 'ERROR Int.b1.2: b1 called for s=m1=m2=0'
		return
	  end if

      	  if (a1.eq.0.d0) then
      	    b1 = - del/2 + log(a2*a2/mu_tH)/2 - 0.25d0
      	  else if (a2.eq.0.d0) then
      	    b1 = - del/2 + log(a1*a1/mu_tH)/2 - 0.75d0
      	  else
      	    b1 = - del/2 + log(a1*a2/mu_tH)/2
      	    if (a1.ne.a2) then
      	      dif = a1*a1 - a2*a2
      	      b1 = b1 - 0.75d0 - a2*a2/2/dif &
&     	         + (a1**4/dif/dif - 0.5d0)*log(a1/a2)
      	    end if
 
     	  end if
      	  return

!	ERROR HANDLING if first argument isn't 0
	else
		call logerror('ERROR Int.b1.1: bad parameters, see README')
		b1 = 0
		return

      	end if

!	COMMENT: The code below was commented out because it shouldn't be
!	necessary and was producing weird recursive errors. I replaced it
!	with the 'else' statement above.


!!	COMMENT: Error flag to say that we're accessing an unexpected part
!!	of this function. We shouldn't need this for Bmumu.
!!	WRITE(*,*) 'ERROR Int.b1.1: see README.txt'
	

!!	COMMENT: These cases are irrelevant  
!    	if ((a1*a2).eq.0.d0) then
!      	  if (a1.eq.a2) then
!      	    b1 = - (del - log(abs(s)/mu_tH))/2 - 1
!!      	    if (s.gt.0.d0) b1 = b1 - (0,1)*pi/2
!!	    COMMENT: The above adds an imaginary part
!!	    if (s.gt.0.d0) WRITE(*,*) 'ERROR Int.b1.2: Imaginary Part!'
!      	  else
!!      	    if (a1.eq.0.d0) then
!!      	      b1 = - (del - log(a2*a2/mu_tH) + 1 &
!!&     	         + (1 - a2*a2/s)*f(s,a1,a2))/2
!!      	    else
!!      	      b1 = - (del - log(a1*a1/mu_tH) + 1 &
!!&     	         + (1 + a1*a1/s)*f(s,a1,a2))/2
!!      	    end if
!!	    WRITE(*,*) 'ERROR Int.b1.2: Imaginary Part!'
!      	  end if
!      	  return
!      	end if

!!      	b1 = - (del - log(a1*a2/mu_tH) &
!!&     	   + (s  + a1*a1 - a2*a2)/s*f(s,a1,a2))/2
!!	COMMENT: I don't think we need the above, error otherwise
!!	WRITE(*,*) 'ERROR Int.b1.3: b1 trying to access f functions.'

!!	COMMENT: I'm commenting out the stuff below too
!!      	if (a1.ne.a2) then
!!      	  sum = a1*a1 + a2*a2
!!      	  dif = a1*a1 - a2*a2
!!      	  b1 = b1 - (1 - sum/dif*log(a1/a2))/2
!!      	end if
!!      	return

      	END FUNCTION b1
	

! ----------------------------------------------------------------
!	AUXILLIARY FUNCTIONS
!	I don't really care what these correspond to!
!	The b functions seem to depend on them
! ----------------------------------------------------------------

!	Define auxilliary functions here, e.g. f functions


! ----------------------------------------------------------------
!	SELF ENERGY INTEGRAL
!	This integral diverges, include del and 't Hooft mass
!
!		/   d^4k	       p*k
!	b1 = 	|  ------- -------------------------
!		/  (2pi)^4 ((k-p)^2-a1^2) (k^2-a2^2)
!
!	Depends on function f, see below
! ----------------------------------------------------------------


      	COMPLEX(KIND=DBL) FUNCTION bc1(s,a1,a2)
	REAL(KIND=DBL), INTENT(IN)::s,a1,a2 !,del,mu_tH
	REAL(KIND=DBL)::dif

!	COMMENT: I have no idea what dbstat is for
!      	if (dbstat) then
!      	  b1 = db1(s,a1,a2)
!      	  return
!      	end if

!	COMMENT: I have no idea what infstat is for
!      	if (infstat) then
!      	  b1 = - del/2
!      	  return
!      	end if

!	COMMENT: This is the case we care about
      	if (s.eq.0.d0) then

      	  if (abs(a1) + abs(a2).eq.0.d0) then
		call logerror('ERROR Int.b1.2: b1 called for s=m1=m2=0')
		bc1=0
!		stop 'ERROR Int.b1.2: b1 called for s=m1=m2=0'
		return
	  end if

      	  if (a1.eq.0.d0) then
      	    bc1 = - del/2 + log(a2*a2/mu_tH)/2 - 0.25d0
      	  else if (a2.eq.0.d0) then
      	    bc1 = - del/2 + log(a1*a1/mu_tH)/2 - 0.75d0
      	  else
      	    bc1 = - del/2 + log(a1*a2/mu_tH)/2
      	    if (a1.ne.a2) then
      	      dif = a1*a1 - a2*a2
      	      bc1 = bc1 - 0.75d0 - a2*a2/2/dif &
&     	         + (a1**4/dif/dif - 0.5d0)*log(a1/a2)
      	    end if
 
     	  end if
      	  return

!	ERROR HANDLING if first argument isn't 0
	else
		call logerror('ERROR Int.b1.1: bad parameters, see README')
		bc1 = 0
		return

      	end if


      	END FUNCTION bc1

	

END MODULE integrals



