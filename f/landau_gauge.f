*********************************************************
*                 Landau gauge fixing                   *
* output: iter_or, da                                   *
*********************************************************

      subroutine landau_gauge(da1,iter_or1)

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
      INCLUDE 'measurement.inc'
      INCLUDE 'input_parameter.inc'
*
*local VARIABLES
*
      integer eo,mu,is,iter_or1,ig
      real da1,rx, max

      da1=10000.0
      iter_or1=0

* copy unfixed gauge configuration to the buffer
      
      do is=1,ownhalfvol(3)
       do eo=0,1
        do mu=1,3
         u_f(1,is,eo,mu)=u(1,is,eo,mu)
         u_f(2,is,eo,mu)=u(2,is,eo,mu)
        enddo
       enddo
      enddo


      do is=1,ownhalfvol(3)
       do eo=0,1
         a_f(1,is,eo)=a(1,is,eo)
         a_f(2,is,eo)=a(2,is,eo)
         a_f(3,is,eo)=a(3,is,eo)
       enddo
      enddo

* perform overrelaxation till condition |d_mu A_mu|^2 is satisfied with
* the precision limit_or

130   if(da1.gt.limit_or) then
       iter_or1=iter_or1+1
        call or_gauge_fixing(da1)
        call u_f_renorm
        goto 130
       endif

      return 
       
      end






