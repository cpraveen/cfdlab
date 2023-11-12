      SUBROUTINE GetArg (n,arg)
      CHARACTER*64 arg
! Get first command line argument, arg, using the Lahey non-standard 
! intrinsic GetCl. Works with the GetCl from Lahey F77L, LF90 and LF95. 
! n is a placeholder variable with no meaning
! Ajit J. Thakkar, Chemistry Dept., U. of New Brunswick
      CALL GetCl(arg)
      CALL Crop (arg) 
      RETURN
      END
      SUBROUTINE Crop (Word)
!  Crop tidies up the character string Word by left-justifying
!   and blanking out everything after the first trailing blank.
! Ajit J. Thakkar, Chemistry Dept., U. of New Brunswick
      IMPLICIT NONE
      CHARACTER*(*) Word
      INTEGER I,J,K
      CHARACTER*1 Blank
      PARAMETER (Blank=' ')
      DO I=1,LEN(Word)
        IF(Word(I:I).NE.Blank)THEN
          J=I
          GOTO 20
        END IF
      END DO
      RETURN
   20 K=0
      DO I=J,LEN(Word)
        IF(Word(I:I).NE.Blank)THEN
          K=K+1
          Word(K:K)=Word(I:I)
        ELSE
          EXIT
        END IF
      END DO
      DO I=K+1,LEN(Word)
       Word(I:I)=Blank
      END DO
      RETURN
      END
