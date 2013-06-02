****************************************************************
*   subroutine for normalization of gauge fixed links          * 
*   to unitary SU(2) matricies                                 * 
****************************************************************
     
      subroutine u_f_renorm

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
      integer mu,eo,is
      real tmp(ownhalfvol(3))

      do mu=1,3
       do eo=0,1
        do is=1,ownhalfvol(3)
         tmp(is)=u_f(1,is,eo,mu)*conjg(u_f(1,is,eo,mu))+
     &   u_f(2,is,eo,mu)*conjg(u_f(2,is,eo,mu))
         tmp(is)=1.0/sqrt(tmp(is))
         u_f(1,is,eo,mu)=tmp(is)*u_f(1,is,eo,mu)        
         u_f(2,is,eo,mu)=tmp(is)*u_f(2,is,eo,mu)
        enddo
       enddo
      enddo

      return
      end

