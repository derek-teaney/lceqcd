
******************************************************************
* subroutine for the calculation of the force                    *
*                                                                *
******************************************************************

      subroutine a_force(eo)
      
      implicit none

*  include MPI definition and parameters
       
      INCLUDE 'mpif.h'
      INCLUDE 'parallel_parameter.inc'
      INCLUDE 'parallel_mpi_types.inc'


*   include global variables and parameters

      INCLUDE 'parameter.inc'
      INCLUDE 'input_parameter.inc'
      INCLUDE 'pointer.inc'
      INCLUDE 'gaugefield.inc'

*   local variables 
       
      integer mu,eo,oe,is
      complex u_tmp1(2,ownhalfvol(3)),u_tmp2(2,maxarray),
     &         u_tmp3(2,ownhalfvol(3))

*****************************************************************      
*          +                 +  
*         U(x-mu)           U(x)
*    /-------<-------\/-------<-------\               /\ nu 
*   oe------->-------eo------->-------oe               | 
* A(x-mu)  U(x-mu)   A(x)    U(x)     A(x+mu)          |
*                                                      ------> 
*                                                            mu
* 1) transfer A(x+mu) in -mu direction
* 2) calculate A(x-mu)*U(x-mu)=u_tmp1(x-mu)
* 3) calculate U^dagger(x-mu)*u_tmp1(x-mu)=u_tmp2(x-mu)
* 4) wait for completing transfer 
* 5) transfer u_tmp2(x-mu) in mu direction
* 6) calculate A(x+mu)*U^dagger(x)=u_tmp3(x)
* 7) calculate U(x)*u_tmp3(x)=u_tmp1(x)
* 8) wait for completing transfer
* 9) calculate the force
*****************************************************************

      oe=1-eo

      force=0.0

      do mu=1,3

* 1)
         if(dimproc(mu).gt.1)then
            call MPI_ISSEND(a(1,1,oe),1,a_VECTOR_dnTYPE(mu,oe),neigh(-mu),
     &                      1000+oe+2*mu,comm,send_request,ierror)
            call MPI_IRECV(a(1,nsiteshalf+1,oe),1,a_CONT_TYPE,neigh(mu),
     &                     1000+oe+2*mu,comm,rec_request,ierror)
         endif
* 2)      
         call multallan(a(1,1,oe),u(1,1,oe,mu),u_tmp1(1,1),maxarray,
     &                  ownhalfvol(3),normalsites(1),normalsites(1))
* 3)
         call multalldn(u(1,1,oe,mu),u_tmp1(1,1),u_tmp2(1,1),maxarray,
     &                  ownhalfvol(3),normalsites(1),normalsites(1))
* 4)      
         if(dimproc(mu).gt.1)then
            call MPI_WAIT(send_request,send_status,ierror)
            call MPI_WAIT(rec_request,rec_status,ierror)
         endif

*         call barrier(comm,ierror)
* 5)      
         if(dimproc(mu).gt.1)then
            call MPI_ISSEND(u_tmp2(1,1),1,u_VECTOR_upTYPE(mu,oe),neigh(mu),
     &                      1000+oe+2*mu,comm,send_request,ierror)
            call MPI_IRECV(u_tmp2(1,nsiteshalf+1),1,u_CONT_TYPE,neigh(-mu),
     &                     1000+oe+2*mu,comm,rec_request,ierror)
         endif

* 6)
         call multallad(a(1,1,oe),u(1,1,eo,mu),u_tmp3(1,1),maxarray,
     &                  ownhalfvol(3),siteup(1,eo,mu),normalsites(1))
* 7)
         call multallnn(u(1,1,eo,mu),u_tmp3(1,1),u_tmp1(1,1),maxarray,
     &                  ownhalfvol(3),normalsites(1),normalsites(1))

* 8)       
         if (dimproc(mu).gt.1)then
            call MPI_WAIT(send_request,send_status,ierror)
            call MPI_WAIT(rec_request,rec_status,ierror)
         endif
      
*         call barrier(comm,ierror)
* 9)
         do is=1,ownhalfvol(3)
            force(3,is)=force(3,is) +
     &            aimag(u_tmp1(1,is) + u_tmp2(1,sitedn(is,eo,mu)))

            force(2,is)=force(2,is) +  
     &            real(u_tmp1(2,is) + u_tmp2(2,sitedn(is,eo,mu)))

            force(1,is)=force(1,is) + 
     &           aimag(u_tmp1(2,is) + u_tmp2(2,sitedn(is,eo,mu)))
         enddo
      
      enddo


      return

      end


