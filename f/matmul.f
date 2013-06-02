      SUBROUTINE MULTALLNN(U1, U2, HELP, MAXSIZE, SIZE,
     &           NEWSITE1, NEWSITE2)

**********************************************************************
*This subroutine caclulates the product of two normal links U1*U2      
*for all lattice sites
**********************************************************************

      IMPLICIT NONE
      
      INTEGER MAXSIZE, SIZE, I, NEWSITE1(SIZE), NEWSITE2(SIZE),
     &        SITE1,  SITE2
      COMPLEX U1(2,MAXSIZE), U2(2,MAXSIZE), HELP(2,MAXSIZE)

***********************************************************************

      DO  I = 1, SIZE

        SITE1 = NEWSITE1(I)
        SITE2 = NEWSITE2(I)
        
        HELP(1,I)=U1(1,SITE1)*U2(1,SITE2)-U1(2,SITE1)*CONJG(U2(2,SITE2))
        HELP(2,I)=U1(1,SITE1)*U2(2,SITE2)+U1(2,SITE1)*CONJG(U2(1,SITE2))

      ENDDO

      RETURN
 
      END

      SUBROUTINE MULTALLDN(U1, U2, HELP, MAXSIZE, SIZE,
     &           NEWSITE1, NEWSITE2)

**********************************************************************
*This subroutine caclulates the product of U1^+*U2
*for all lattice sites
**********************************************************************

      IMPLICIT NONE

      INTEGER MAXSIZE, SIZE, I, NEWSITE1(SIZE), NEWSITE2(SIZE),
     &        SITE1,  SITE2
      COMPLEX U1(2,MAXSIZE), U2(2,MAXSIZE), HELP(2,MAXSIZE)

***********************************************************************

      DO  I = 1, SIZE

        SITE1 = NEWSITE1(I)
        SITE2 = NEWSITE2(I)

        HELP(1,I)=CONJG(U1(1,SITE1))*U2(1,SITE2)+
     &  U1(2,SITE1)*CONJG(U2(2,SITE2))

        HELP(2,I)=CONJG(U1(1,SITE1))*U2(2,SITE2)-
     &  U1(2,SITE1)*CONJG(U2(1,SITE2))

      ENDDO
    
      RETURN
      END


      SUBROUTINE MULTALLND(U1, U2, HELP, MAXSIZE, SIZE,
     &           NEWSITE1, NEWSITE2)

**********************************************************************
*This subroutine caclulates the product of two links U1*U2^+
*for all lattice sites
**********************************************************************

      IMPLICIT NONE

      INTEGER MAXSIZE, SIZE, I, NEWSITE1(SIZE), NEWSITE2(SIZE),
     &        SITE1,  SITE2
      COMPLEX U1(2,MAXSIZE), U2(2,MAXSIZE), HELP(2,MAXSIZE)

***********************************************************************

      DO  I = 1, SIZE

        SITE1 = NEWSITE1(I)
        SITE2 = NEWSITE2(I)

        HELP(1,I)=U1(1,SITE1)*CONJG(U2(1,SITE2))+
     &  U1(2,SITE1)*CONJG(U2(2,SITE2))

        HELP(2,I)=-U1(1,SITE1)*U2(2,SITE2)+
     &  U1(2,SITE1)*U2(1,SITE2)
     
      ENDDO

      RETURN
    
      END

      SUBROUTINE MULTALLDD(U1, U2, HELP, MAXSIZE, SIZE,
     &           NEWSITE1, NEWSITE2)

**********************************************************************
*This subroutine caclulates the product of two  links U1^+*U2^+
*for all lattice sites
**********************************************************************

      IMPLICIT NONE

      INTEGER MAXSIZE, SIZE, I, NEWSITE1(SIZE), NEWSITE2(SIZE),
     &        SITE1,  SITE2
      COMPLEX U1(2,MAXSIZE), U2(2,MAXSIZE), HELP(2,MAXSIZE)

***********************************************************************

      DO  I = 1, SIZE

        SITE1 = NEWSITE1(I)
        SITE2 = NEWSITE2(I)

        HELP(1,I)=CONJG(U1(1,SITE1))*CONJG(U2(1,SITE2))-
     &  U1(2,SITE1)*CONJG(U2(2,SITE2))

        HELP(2,I)=-CONJG(U1(1,SITE1))*U2(2,SITE2)-
     &  U1(2,SITE1)*U2(1,SITE2)
     
      ENDDO

      RETURN
      END

      SUBROUTINE MULTALLAD(A, U, HELP, MAXSIZE, SIZE,
     &           NEWSITE1, NEWSITE2)

**********************************************************************
*This subroutine caclulates the product of A field and daggered link U^+
*for all lattice sites
**********************************************************************

      IMPLICIT NONE

      INTEGER MAXSIZE, SIZE, I, NEWSITE1(SIZE), NEWSITE2(SIZE),
     &        SITE1,  SITE2
      COMPLEX U(2,MAXSIZE), HELP(2,MAXSIZE)
      REAL A(3,MAXSIZE)
***********************************************************************

      DO  I = 1, SIZE

        SITE1 = NEWSITE1(I)
        SITE2 = NEWSITE2(I)

        HELP(1,I)=CMPLX(0.0,A(3,SITE1))*CONJG(U(1,SITE2))+
     &  CMPLX(A(2,SITE1),A(1,SITE1))*CONJG(U(2,SITE2))

        HELP(2,I)=-CMPLX(0.0,A(3,SITE1))*U(2,SITE2)+
     &  CMPLX(A(2,SITE1),A(1,SITE1))*U(1,SITE2)  

      ENDDO

      RETURN
      
      END

      SUBROUTINE MULTALLAN(A, U, HELP, MAXSIZE, SIZE,
     &           NEWSITE1, NEWSITE2)

**********************************************************************
*This subroutine caclulates the product of A-field and an SU(2) matrix
*for all lattice sites
**********************************************************************

      IMPLICIT NONE

      INTEGER MAXSIZE, SIZE, I, NEWSITE1(SIZE), NEWSITE2(SIZE),
     &        SITE1,  SITE2
      COMPLEX U(2,MAXSIZE),  HELP(2,MAXSIZE)
      REAL A(3,MAXSIZE)
***********************************************************************

      DO  I = 1, SIZE

        SITE1 = NEWSITE1(I)
        SITE2 = NEWSITE2(I)
      
        HELP(1,I)=CMPLX(0.0,A(3,SITE1))*U(1,SITE2)-
     &  CMPLX(A(2,SITE1),A(1,SITE1))*CONJG(U(2,SITE2))

        HELP(2,I)=CMPLX(0.0,A(3,SITE1))*U(2,SITE2)+
     &  CMPLX(A(2,SITE1),A(1,SITE1))*CONJG(U(1,SITE2))


       ENDDO
     

      RETURN
      END
