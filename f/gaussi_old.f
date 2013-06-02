**********************************************************************
*                    Gaussian Random numbers                         *
**********************************************************************

      subroutine gauss ( grndm , width )

      implicit none

      INCLUDE 'parameter.inc'
      INCLUDE 'parallel_parameter.inc'     

      real grndm(3,ownhalfvol(3)),help1(ownhalfvol(3)),
     &     help2(ownhalfvol(3)),width(ownhalfvol(3)),
     &     rand1(ownhalfvol(3)),rand2(ownhalfvol(3))
      real tpi,rx
      integer ier,ig,is


      tpi=2.0*acos(-1.0)

      do 100 ig=1,3

        do is=1,ownhalfvol(3)
               call random_number(rx)
               rand1(is)=rx
               call random_number(rx)
               rand2(is)=rx
        enddo

      where(rand1.eq.0.0)
      rand1=0.0000055
      endwhere
      help2= 1.0 / width
      help1=-1.0*help2*log(rand1)
      grndm(ig,:)=sqrt(help1)*cos(tpi*rand2)
 100  continue

      return
      end
  










