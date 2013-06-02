c **********************************************************************
c     Kennedy Pendleton Heat Bath
c **********************************************************************
c
c     input varibales : b, a0_old,1/2 Tr(U stap)
c     ouput variables : v, nacch
c **********************************************************************

      subroutine heat_kpi( b, v, a0_old, nacch )

      implicit none
*
*     include global definitions
* 
      INCLUDE 'parallel_parameter.inc'

      INTRINSIC COUNT
*     
*    definition of the local variables
*
      complex v(2,ownhalfvol(3))
      real a0(ownhalfvol(3)),a0_old(ownhalfvol(3)),a1(ownhalfvol(3)),
     &     a2(ownhalfvol(3)),a3(ownhalfvol(3)),b(ownhalfvol(3))
      real tpi,rx
      integer ier,site,nacch
      logical lacc(ownhalfvol(3))


      tpi=2.0*acos(-1.0)
c
c create trials for 1-a0
c
      DO site=1,ownhalfvol(3)
               call random_number(rx)
               a1(site)=rx
               call random_number(rx)
               a2(site)=rx
               call random_number(rx)
               a3(site)=rx
      ENDDO

      where ( a1.eq.0.0 )
      a1=0.0000033
      endwhere
      where ( a2.eq.0.0 )
      a2=0.0000066
      endwhere

      a0=-(log(a1)+log(a2)*cos(tpi*a3)**2)*b
c
c decide which ones to accept
c
      DO site=1,ownhalfvol(3)
          call random_number(rx)
               a1(site)=rx
      ENDDO

      lacc = a1**2 .le. (1.0 - 0.5*a0)
      nacch = count( lacc )
      where ( lacc )
      a0 = 1.0 - a0
      elsewhere
      a0 = a0_old
      endwhere
c
c now create the angles randomly and form SU(2) matrix
c
      a1 = 1.0 - a0**2
      where ( a1 .ge. 0.0 )
      a1 = sqrt( a1 )
      elsewhere
      a1 = 0.0
      endwhere

      DO site=1,ownhalfvol(3)
               call random_number(rx)
               a2(site)=rx
      ENDDO

      a2 = 2.0 * ( a2 - 0.5 )
      a3 = a1 * a2
      v(1,:) = cmplx( a0, a3)
      a0 = 1.0 - a2**2
      where ( a0 .ge. 0.0 )
      a0 = a1 * sqrt( a0 )
      elsewhere
      a0 = 0.0
      endwhere

            DO site=1,ownhalfvol(3)
               call random_number(rx)
               a3(site)=rx
            ENDDO

      a3 = tpi * a3
      a1 = a0 * cos ( a3 )
      a2 = a0 * sin ( a3 )
      v(2,:) = cmplx( a2, a1)

      return
      end













