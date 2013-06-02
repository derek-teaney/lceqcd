  
**********************************************************************
*                   normalize the gauge fields                       *
**********************************************************************

      subroutine renorm

      implicit none

*
*including all global MESSAGE PASSING features and variables
*
      INCLUDE 'parallel_parameter.inc'

*
*including global variables, parameters and common blocks
*
      INCLUDE 'parameter.inc'
      INCLUDE 'gaugefield.inc'

      integer eo,mu,is
      real tmp(ownhalfvol(3))
      
      do is=1,ownhalfvol(3)
       do mu=1,3
        do eo=0,1
         tmp(is)=real( u(1,is,eo,mu)*conjg(u(1,is,eo,mu))+
     &   u(2,is,eo,mu)*conjg(u(2,is,eo,mu)))
         tmp(is)=1.0/sqrt(tmp(is))
         u(1,is,eo,mu)=tmp(is)*u(1,is,eo,mu)
         u(2,is,eo,mu)=tmp(is)*u(2,is,eo,mu)
        enddo
       enddo
      enddo

      return
      end

