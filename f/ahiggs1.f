*******************************************************************
*   subroutine for higgs field updating with Metropolis algorithm *
*   input  : delta, measflag                                      *
*   output : accept                                               *
*******************************************************************
     
      subroutine ahiggs1(delta,measflag,accept)
      
      implicit none

*******************************************************************

* include message passing variables
 
      INCLUDE 'mpif.h'
      INCLUDE 'parallel_parameter.inc'
      INCLUDE 'parallel_mpi_types.inc'

* include global variables and common blocks
 
      INCLUDE 'parameter.inc'
      INCLUDE 'input_parameter.inc'
      INCLUDE 'pointer.inc'
      INCLUDE 'gaugefield.inc'
      INCLUDE 'measurement.inc'

* definitions of local variables 

      integer eo,site,i,measflag,acc
      real accept,rx,delta,
     & a_n(3,ownhalfvol(3)),sf_o(ownhalfvol(3)),sf_n(ownhalfvol(3)),
     & rho_o(ownhalfvol(3)),rho_n(ownhalfvol(3)),boltz(ownhalfvol(3))

*********************************************************************

      a0kin=0.0
      a02=0.0
      a04=0.0
      acc=0.0

c     start loop over colour

      do eo=0,1
      
*   calculate the force

       call a_force(eo) 
       
*   propose a new configuration for the Higgs field       

       do site=1,ownhalfvol(3)
        do i=1,3
        call random_number(rx)
        a_n(i,site)=a(i,site,eo)+(rx-0.5)*delta
        enddo
       enddo

*   calculate the Higgs field length and cinetic term
*   for new and old configurations

       sf_o=0.0
       sf_n=0.0
       rho_o=0.0
       rho_n=0.0

       do site=1,ownhalfvol(3)
         do i=1,3
          sf_o(site)=sf_o(site)+a(i,site,eo)*force(i,site)
          sf_n(site)=sf_n(site)+a_n(i,site)*force(i,site)
          rho_o(site)=rho_o(site)+a(i,site,eo)**2
          rho_n(site)=rho_n(site)+a_n(i,site)**2
         enddo
       enddo
  
*   calculate the ratio of boltzmann factors

       boltz=exp(-beta_a*(sf_o-sf_n)+beta_2*(rho_o-rho_n)+
     &           beta_4*(rho_o**2-rho_n**2))

*   accept reject procedure

       do site=1,ownhalfvol(3)
          call random_number(rx)
          if (boltz(site).ge.rx) then
             a(1,site,eo)=a_n(1,site)  
             a(2,site,eo)=a_n(2,site)
             a(3,site,eo)=a_n(3,site)
             rho_o(site)=rho_n(site)
             sf_o(site)=sf_n(site)
             acc=acc+1
          endif
       enddo
       
*      measurements
       
       if (measflag.eq.1)then
       a02=a02+sum(rho_o)
       a04=a04+sum(rho_o**2)       
       a0kin=a0kin+sum(sf_o)
       endif
 
       enddo
c  end loop over colour 

*  normalize measurements
      
      if (measflag.eq.1)then
        a02=a02/real(ownvol(3))
        a04=a04/real(ownvol(3))
        a0kin=a0kin/real(6*ownvol(3))
      endif

       accept=real(acc)/real(ownvol(3))

*  adjust delta to keep acceptance rate close to 0.5
   
      if (abs(accept-0.5).gt.0.1)then
      if (accept.gt.0.5)then
      delta=1.1131241*delta
      else
      delta=0.8989669*delta
      endif
      endif

      return
      
      end









































