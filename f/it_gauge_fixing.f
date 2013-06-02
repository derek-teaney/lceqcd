***************************************************************
*              Iterative gauge fixing                         * 
* ioutput: da1                                                *
***************************************************************
      subroutine it_gauge_fixing(da1)

      implicit none

*
*including all global MESSAGE PASSING features and variables
*
      INCLUDE 'mpif.h'
      INCLUDE 'parallel_parameter.inc'
      INCLUDE 'parallel_mpi_types.inc'

*
*including global variables, parameters and common blocks
*
      INCLUDE 'parameter.inc'
      INCLUDE 'input_parameter.inc'
      INCLUDE 'pointer.inc'
      INCLUDE 'gaugefield.inc'
     
      real    da1,da11
      real    t1(ownhalfvol(3)),t2(ownhalfvol(3))
      integer eo,oe,mu,ig,is
      complex Omega(2,maxarray,0:1), g1(2,ownhalfvol(3)), g2(2,ownhalfvol(3)),
     &       g(2,ownhalfvol(3)),u_tmp1(2,ownhalfvol(3)),u_tmp2(2,ownhalfvol(3))

*           
* Calculate i d_mu A_mu, do it in 5 steps
*                  
* 1) Send and recive U (x-mu) from oe
*                     mu       
*                           +     
* 2) Calculate 1/2 ( U (x)-U (x) )=g1
*                         mu    mu    
*
* 3) Wait to complete the transfer
*
*                              +
* 4) Calculate 1/2 ( U (x-mu)-U (x-mu) )=g2
*                     mu       mu
*
* 5) Calculate the matrix g= i d_mu A_mu= g1-g2
*


      da11=0.0

      do eo=0,1 
*     start lopp over even-odd

      oe=mod(eo+1,2)

      g1=cmplx(0.0,0.0)
      g2=cmplx(0.0,0.0)

      do mu=1,3
* 1) 
      if (dimproc(mu).gt.1) then
      call MPI_ISSEND(u_f(1,1,oe,mu),1,u_VECTOR_upTYPE(mu,oe),
     &                neigh(mu),1000+oe+2*mu,comm,send_request,ierror)
      call MPI_IRECV(u_f(1,nsiteshalf+1,oe,mu),1,u_CONT_TYPE,
     &               neigh(-mu),1000+oe+2*mu,comm,rec_request,ierror)
      endif

* 2)    
      do is=1,ownhalfvol(3)
        g1(1,is)=g1(1,is)+cmplx(0.0, aimag(u_f(1,is,eo,mu)) )
        g1(2,is)=g1(2,is)+u_f(2,is,eo,mu)
      enddo
* 3)
      if (dimproc(mu).gt.1) then
      call MPI_WAIT(send_request,send_status,ierror)
      call MPI_WAIT(rec_request,rec_status,ierror)
      endif

* 4) 
      do is=1,ownhalfvol(3)
       g2(1,is)=g2(1,is)+cmplx(0.0,aimag(u_f(1,sitedn(is,eo,mu),oe,mu)) )
       g2(2,is)=g2(2,is)+u_f(2,sitedn(is,eo,mu),oe,mu)
      enddo

      enddo

* 5)  
      g(1,:)=g1(1,:)-g2(1,:)
      g(2,:)=g1(1,:)-g2(2,:)

 
*
* calculate the average norm of the g matrix sum |d_mu A_mu|**2= sum |g|
*

      if (eo.eq.1) then
        da11=sum( aimag(g(1,:))**2+g(2,:)*conjg(g(2,:)))
        da11=da11/real(ownhalfvol(3))
      endif
*
* calculate  alpha*d_mu A_mu
*

       g(1,:)=alpha*g(1,:)
       g(2,:)=alpha*g(2,:)

*
* calculate the gauge transformation Omega(x)=exp(i alpha d_mu A_mu)
* 
* Omega(x)=1*cos(|g|)+ g*sin(|g|)/|g|
*  

      
      t1=sqrt( aimag(g(1,:))**2+real(g(2,:))**2+aimag(g(2,:))**2 )
      
        do is=1,ownhalfvol(3)
          t2(is)=sin(t1(is))/t1(is)
          Omega(1,is,eo)=g(1,is)*t2(is)+cos(t1(is))
          Omega(2,is,eo)=g(2,is)*t2(is)
        enddo

* end loop over even odd
*      enddo

*
*  Perform the gauge transformation
*                                                                +
*  U (x)  ->Omega(x) U  (x) Omega(x+mu) 
*   mu                mu   
*                           +
*  A->Omega(x) A(x) Omega(x)  
*
*  Do this in 5 steps:
*  1) send Omega(x+mu) in -mu direction
*  2) Omega(x) U(x)=u_tmp1
*  3) wait to complete the transfer 
*  4) u_tmp1 Omega(x+mu)
*  5) Omega(x) A(x) Omega(x)^dagg

*      do eo=0,1
* start loop over even odd
*      oe=1-eo

      do mu=1,3
* 1)       
       if (dimproc(mu).gt.1) then
       call MPI_ISSEND(Omega(1,1,oe),1,u_VECTOR_dnTYPE(mu,oe),
     &                neigh(-mu),1000+oe+2*mu,comm,send_request,ierror)
       call MPI_IRECV(Omega(1,nsiteshalf+1,oe),1,u_CONT_TYPE,
     &               neigh(mu),1000+oe+2*mu,comm,rec_request,ierror)
       endif

* 2)
       call multallnn(Omega(1,1,eo),u_f(1,1,eo,mu),u_tmp1(1,1),
     &                maxarray,ownhalfvol(3),normalsites(1),normalsites(1))
        
         do is=1,ownhalfvol(3)
           u_f(1,is,eo,mu)=u_tmp1(1,is)
           u_f(2,is,eo,mu)=u_tmp1(2,is)
         enddo 
     
* 3)
       if (dimproc(mu).gt.1) then
       call MPI_WAIT(send_request,send_status,ierror)
       call MPI_WAIT(rec_request,rec_status,ierror)
       endif

* 4)       
       call multallnd(u_f(1,1,oe,mu),Omega(1,1,eo),u_tmp2(1,1),
     &                maxarray,ownhalfvol(3),normalsites(1),siteup(1,oe,mu))
      
         do is=1,ownhalfvol(3)
           u_f(1,is,oe,mu)=u_tmp2(1,is)
           u_f(2,is,oe,mu)=u_tmp2(2,is)
         enddo

       enddo

* 5)
       call multallad(a_f(1,1,eo),Omega(1,1,eo),u_tmp1(1,1),
     &                maxarray,ownhalfvol(3),normalsites(1),normalsites(1))
       call multallnn(Omega(1,1,eo),u_tmp1(1,1),u_tmp2(1,1),
     &                maxarray,ownhalfvol(3),normalsites(1),normalsites(1))

       do is=1,ownhalfvol(3)
        a_f(1,is,eo)=aimag(u_tmp2(2,is))
        a_f(2,is,eo)=real(u_tmp2(2,is))
        a_f(3,is,eo)=aimag(u_tmp2(1,is))
       enddo

       enddo 
* end loop over even odd

      if (nproc.gt.1) then
       call MPI_ALLREDUCE(da11,da1,1,MPI_REAL,MPI_SUM,
     &                    comm,ierror)
       da1=da1/nproc
      else
       da1=da11
      endif

       

       return
      
       end








