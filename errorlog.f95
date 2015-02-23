! ----------------------------------------------------------------
! Project: Complete One loop Predictions for B->mu+mu- in the MSSM
! Please cite: Dedes, Rosiek, Tanedo (2007)
! ----------------------------------------------------------------
! Main file for B->mu+mu- calculation
! filename: errorlog.f95
! Current version: 0.01
! By Flip Tanedo (philip.tanedo@gmail.com)
! Date: 20 November 2007
! Latest version available at: http://www.dur.ac.uk/philip.tanedo
! ----------------------------------------------------------------
! This module keeps track of any error flags encountered at
! runtime. Points in the code that may generate errors automatically
! call the errorlog function to signal when an improper part of the
! code is accessed. For example, some of the Passarino-Veltman
! functions are only defined for certain values of parameters.
! If an invalid parameter is fed into the program, then the
! integral is set to zero (hopefully to make the error obvious) 
! and an error flag is waved. This error is passed through this
! module and the resulting error messages are stored in a linked
! list which are then outputted at the end of runtime or to a
! log file.
!
! ----------------------------------------------------------------
! RECORD OF REVISIONS
!	Date	 Programmer	Description
!	====	 ==========	==================================
!	20/11/07 Flip Tanedo	Original Code
! ----------------------------------------------------------------
! USAGE
! USE errorlog
! CALL logerror(...)		! Flags an error
! CALL printerrors		! Outputs all flagged errors
!
! Note: do not confuse module name (errorlog) with subroutine (logerror)
! ----------------------------------------------------------------
! LINKED LIST STRUCTURE
!	Each errorflag data type contains an error message 
!	(description) and a pointer to the next error.
!	The end of the list is signified by a nullified
!	pointer. 
! ----------------------------------------------------------------


MODULE errorlog
USE program_parameters		! for DBL
IMPLICIT NONE

!	Declare errorflag linked list data type
	TYPE :: errorflag
		CHARACTER(len=70) :: description
		TYPE(errorflag), POINTER :: next
	END TYPE errorflag


!	Data dictionary
	TYPE(errorflag), POINTER::errorhead	! Points to first error
	TYPE(errorflag), POINTER::errortail	! Points to the last error
	INTEGER::istat				! Status: 0 for success
						! for pointer allocation
	CHARACTER(len=70)::terrormessage	! Temp error message

! ----------------------------------------------------------------

CONTAINS

! ----------------------------------------------------------------

	SUBROUTINE logerror(message)
	IMPLICIT NONE
	CHARACTER(len=*), INTENT(IN)::message
	
	IF (.NOT. ASSOCIATED(errorhead)) THEN	! If the linked list is empty
	   ALLOCATE (errorhead,STAT=istat)	! Allocate new value (node)
	   errortail => errorhead		! Tail points to new value, too
	   NULLIFY(errortail%next)		! Nullify link in new value
	   errortail%description = message 	! Store error message in tail

	ELSE					! Values already in list
	   ALLOCATE(errortail%next,STAT=istat)	! Allocate a new value
	   errortail => errortail%next		! Tail points to latest node
	   errortail%description = message	! Store error message in tail

	END IF

	END SUBROUTINE logerror

! ----------------------------------------------------------------

	SUBROUTINE printerrors			! To print error list
	IMPLICIT NONE
	TYPE(errorflag), POINTER::ptr		! Temporary pointer value
	
	ptr => errorhead			! pointer points to first error

	IF (.NOT. ASSOCIATED(ptr)) THEN		! If the linked list is empty
	   WRITE(*,*) 'Error checker: No errors flagged.'
	   WRITE(*,*) '... Flip is a good programmer!'

	ELSE
	   DO
		IF(.NOT. ASSOCIATED(ptr)) EXIT	! Exit loop if no more data
		WRITE(*,*) ptr%description	! Write error
		ptr => ptr%next			! Go to next pointer
	   END DO

	END IF
	END SUBROUTINE printerrors

! ----------------------------------------------------------------

	SUBROUTINE createlogfile(filename)	! To create error logfile
	IMPLICIT NONE
	CHARACTER(len=*), INTENT(IN)::filename	! filename for logfile
	TYPE(errorflag), POINTER::ptr		! Temporary pointer value
	INTEGER, PARAMETER::IOUNIT=25		! Arbitrary I/O UNIT
	INTEGER::ioerror			! For file I/O

	OPEN(UNIT=IOUNIT, FILE=filename, &	! Open file
&	   STATUS='REPLACE', ACTION ='WRITE', &
&	   IOSTAT=ioerror)		
	
	ptr => errorhead			! pointer points to first error

	IF (.NOT. ASSOCIATED(ptr)) THEN		! If the linked list is empty
	   WRITE(IOUNIT,*) 'Error checker: No errors flagged.'
	   WRITE(IOUNIT,*) '... Flip is a good programmer!'

	ELSE
	   DO
		IF(.NOT. ASSOCIATED(ptr)) EXIT	! Exit loop if no more data
		WRITE(IOUNIT,*) ptr%description	! Write error
		ptr => ptr%next			! Go to next pointer
	   END DO

	END IF

	CLOSE(UNIT=IOUNIT)			! Release I/O UNIT

	END SUBROUTINE createlogfile

! ----------------------------------------------------------------



END MODULE errorlog
