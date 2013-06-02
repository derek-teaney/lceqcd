*********************************************************
*              Testing   Landau gauge fixing            *
* output: max, min                                      *
*********************************************************

      subroutine testmax(max,min)

      implicit none

      include 'mpif.h'

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
      include 'pointer.inc'
       

      real max, max1, min 
      integer eo, oe, is, mu 
      complex atmp(2,ownhalfvol(3),0:1,3), tmp(2,ownhalfvol(3),0:1)
      max1=0.0

      do eo=0,1
       do is=1,ownhalfvol(3)
        max1=max1+real(u_f(1,is,eo,3))
        do mu=1,2 
         max1=max1+real(u_f(1,is,eo,mu))
        enddo
       enddo
      enddo

      max=max1
c     CALL MPI_REDUCE(max1,max,1,MPI_REAL,MPI_SUM,root,
c     &                comm,ierror)

* check the condition d_mu A_mu=0 if nproc=1

      min=0.0

      if (nproc .gt. 1) goto 11


* define the vector fields

      do eo=0,1
       do is=1,ownhalfvol(3)
        do mu=1,3
         atmp(1,is,eo,mu)=0.5*( u_f(1,is,eo,mu)-conjg(u_f(1,is,eo,mu)) )
         atmp(2,is,eo,mu)=u_f(2,is,eo,mu)
        enddo
       enddo
      enddo


* calculate d_mu A_mu      

      tmp=0.0

      do eo=0,1
       oe=1-eo
       do is=1,ownhalfvol(3)

        do mu=1,2
         tmp(1,is,eo)=tmp(1,is,eo)+(atmp(1,is,eo,mu) -
     &                    atmp(1,sitedn(is,eo,mu),oe,mu) )
         tmp(2,is,eo)=tmp(2,is,eo)+(atmp(2,is,eo,mu) -
     &                    atmp(2,sitedn(is,eo,mu),oe,mu) )
        enddo

         tmp(1,is,eo)=tmp(1,is,eo) + alpha*(atmp(1,is,eo,3) -
     &                    atmp(1,sitedn(is,eo,3),oe,3) )
         tmp(2,is,eo)=tmp(2,is,eo) + alpha*(atmp(2,is,eo,3) -
     &                    atmp(2,sitedn(is,eo,3),oe,3) )
       enddo
      enddo

* calculate sum | d_mu A_mu |^2

      min=0.0

      do eo=0,1
       do is=1,ownhalfvol(3)
       min=min+tmp(1,is,eo)*conjg(tmp(1,is,eo))+
     &         tmp(2,is,eo)*conjg(tmp(2,is,eo))
       enddo
      enddo

      min = min/ownvol(3)

11    end
