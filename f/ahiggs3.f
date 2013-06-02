********************************************************************
**   subroutine for higgs field updating with  overrelaxation      *
**   input  : measflag                                             *
**   output : accept                                               *
********************************************************************
     
      subroutine ahiggs3(measflag,accept)
      
      implicit none

*******************************************************************

** include message passing variables

      INCLUDE 'mpif.h'
      INCLUDE 'parallel_parameter.inc'
      INCLUDE 'parallel_mpi_types.inc' 

** include global variables and common blocks
 
      INCLUDE 'parameter.inc'
      INCLUDE 'input_parameter.inc'
      INCLUDE 'pointer.inc'
      INCLUDE 'gaugefield.inc'
      INCLUDE 'measurement.inc'

** definitions of local variables 

      integer eo,site,i,measflag,acc
      real accept,rx,
     & a_n(3,ownhalfvol(3)),g1(ownhalfvol(3)),
     & rho_o(ownhalfvol(3)),rho_n(ownhalfvol(3)),boltz(ownhalfvol(3)),
     & conf2(ownhalfvol(3)),alph(ownhalfvol(3)),width(ownhalfvol(3)),
     & cxval(ownhalfvol(3)),g2(ownhalfvol(3)),g3(ownhalfvol(3))

*********************************************************************
**     start loop over colour

      a0kin=0.0
      a02=0.0 
      a04=0.0
      acc=0.0

      do eo=0,1
      
**   calculate the force

       call a_force(eo) 

**------- plane version of the heatbath ---------**
    
      if(imode.eq.0)then
      alph=0.0
      width=beta_2
      cxval=1/beta_2/2
      endif
**-----------------------------------------------**

**---- sophisticated version of the heatbath -----**

      if(imode.eq.1)then

      conf2=0.0
      do i=1,3
       do site=1,ownhalfvol(3)
        conf2(site)=conf2(site)+beta_a*force(i,site)**2
       enddo
      enddo
    
      do site=1,ownhalfvol(3)
        alph(site)=beta_4*conf2(site)
        alph(site)=alph(site)/(2*(beta_2+beta_4*conf2(site))**2) 
        width(site)=beta_2+alph(site)
        cxval(site)=1/2.0/width(site)
      enddo 

      endif
**--------------------------------------------------**     
       
**   propose a new configuration for the Higgs field       


       do site=1,ownhalfvol(3)
        a_n(1,site)=-a(1,site,eo)+2.0*beta_a*force(1,site)*cxval(site)
        a_n(2,site)=-a(2,site,eo)+2.0*beta_a*force(2,site)*cxval(site)
        a_n(3,site)=-a(3,site,eo)+2.0*beta_a*force(3,site)*cxval(site)
       enddo   

**   calculate the Higgs field length and cinetic term
**   for new and old configurations

       rho_o=0.0
       rho_n=0.0

         do i=1,3
          do site=1,ownhalfvol(3)
           rho_o(site)=rho_o(site)+a(i,site,eo)**2
           rho_n(site)=rho_n(site)+a_n(i,site)**2
          enddo 
         enddo
  
**   calculate the ratio of boltzmann factors

       boltz=exp(-alph*(rho_o-rho_n)+
     &           beta_4*(rho_o**2-rho_n**2))

**   accept reject procedure

       do site=1,ownhalfvol(3)
        call random_number(rx)
        if (boltz(site).ge.rx)then
        a(1,site,eo)=a_n(1,site)  
        a(2,site,eo)=a_n(2,site)
        a(3,site,eo)=a_n(3,site)
        rho_o(site)=rho_n(site)
        acc=acc+1
        endif
       enddo
       
**      measurements

       if(measflag.eq.1)then
       a02=a02+sum(rho_o)
       a04=a04+sum(rho_o**2)       
       
       do i=1,3
        do site=1,ownhalfvol(3)
          a0kin=a0kin+a(i,site,eo)*force(i,site)
        enddo 
       enddo
       endif 

       enddo
**  end loop over colour 

**  normalize measurements

      if (measflag.eq.1)then
      a02=a02/real(ownvol(3))
      a04=a04/real(ownvol(3))
      a0kin=a0kin/real(6*ownvol(3))
      accept=real(acc)/real(ownvol(3))
      endif


      return
      
      end







