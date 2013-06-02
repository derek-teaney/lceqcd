****************************************************************
*                                                              *
*     construction of the staples for a defined mu direction   *
*     and half of the sites                                    *
*                                                              *
****************************************************************


      subroutine staples(eo,oe,mu)


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

*
*local VARIABLES
*
      integer nu,nnu,nnnu(2),i,j,site,eo,oe,mu,newsite,
     &        l,m
      complex helpbt1(2,ownhalfvol(3)), helpbt2(2,maxarray),
     &        helpup1(2,ownhalfvol(3)), helpup2(2,ownhalfvol(3))

***********************************************************************


      i=0

*              
*construction of the nu direction with nu not equal mu
*so alltoghether 4 staples can be added together
*
      do 5 nnu=1,3
         if (nnu.eq.mu) go to 5
         i=i+1
         nnnu(i)=nnu
 5    continue


*
*initialisation of the staple
*
      do 15 site=1,ownhalfvol(3)
         v(1,site,eo,mu)=(0.,0.)
         v(2,site,eo,mu)=(0.,0.)
 15   continue



      do 100 nnu=1,2
         nu=nnnu(nnu)

*
*     o--3---e
*     |      |   or (even => odd) and vice versa
*     2      4
*     |      |
*     e--1---o   staple to be computed 
*     |      |   therefore compute up and bottom staple separately
*     5      7   on the even and the odd lattice
*     |      |   and then pass the bottom staple one up 
*     o--6---e   helpbt2 has therefore a tail to store the sent data
*
* do it in five parts:
*


*
*  1a) transfer 3 in -nu direction
*   b) compute 6\dad * 5
*   c) wait for completing the transfer
*
c -- a -----------------------------------------------------------
         if (dimproc(nu).gt.1) then
            CALL MPI_ISSEND(u(1,1,oe,mu),1,u_VECTOR_dnTYPE(nu,oe),
     1           neigh(-nu),1000+oe+2*mu+8*nu,comm,
     2           send_request,ierror)      
            CALL MPI_IRECV(u(1,nsiteshalf+1,oe,mu),1,u_CONT_TYPE,
     1           neigh(nu),1000+oe+2*mu+8*nu,comm,
     2           rec_request,ierror)      
         endif
c -- b --------------------------------------------------------------

         call multalldn(u(1,1,oe,mu),u(1,1,oe,nu),helpbt1(1,1),
     &                  maxarray,ownhalfvol(3),
     &                  normalsites(1),normalsites(1))

c -- c -------------------------------------------------------------         
         if (dimproc(nu).gt.1) then
            CALL MPI_WAIT(send_request,send_status,ierror)      
            CALL MPI_WAIT(rec_request,rec_status,ierror)       
         endif
 
         call mpi_barrier(comm,ierror)
        
*---------------------------------------------------------------------
*
*  2a) transfer 7 in -mu direction
*   b) compute 3\dag * 2\dag
*   c) wait for completing the transfer
*
      
c -- a ----------------------------------------------------------------

         if (dimproc(mu).gt.1) then
            CALL MPI_ISSEND(u(1,1,eo,nu),1,u_VECTOR_dnTYPE(mu,eo),
     1           neigh(-mu),1000+oe+2*mu+8*nu,comm,
     2           send_request,ierror)       
            CALL MPI_IRECV(u(1,nsiteshalf+1,eo,nu),1,u_CONT_TYPE,
     1           neigh(mu),1000+oe+2*mu+8*nu,comm,
     2           rec_request,ierror)       
         endif
c -- b ----------------------------------------------------------------

         call multalldd(u(1,1,oe,mu),u(1,1,eo,nu),helpup1(1,1),
     &                  maxarray,ownhalfvol(3),
     &                  siteup(1,eo,nu),normalsites(1))

c -- c ----------------------------------------------------------------
         
         if (dimproc(mu).gt.1) then
            CALL MPI_WAIT(send_request,send_status,ierror)       
            CALL MPI_WAIT(rec_request,rec_status,ierror)       
         endif

         call mpi_barrier(comm,ierror)

*-----------------------------------------------------------------------
*
*  3a) transfer 4 in -mu direction
*   b) compute 7\dag * (6\dag * 5)
*   c) wait for completing the transfer
*
      
c -- a ----------------------------------------------------------------
         if (dimproc(mu).gt.1) then
            CALL MPI_ISSEND(u(1,1,oe,nu),1,u_VECTOR_dnTYPE(mu,oe),
     1           neigh(-mu),1000+oe+2*mu+8*nu,comm,
     2           send_request,ierror)       
            CALL MPI_IRECV(u(1,nsiteshalf+1,oe,nu),1,u_CONT_TYPE,
     1           neigh(mu),1000+oe+2*mu+8*nu,comm,
     2           rec_request,ierror)      
         endif
c -- b ---------------------------------------------------------------


         call multalldn(u(1,1,eo,nu),helpbt1(1,1),helpbt2(1,1),
     &                  maxarray,ownhalfvol(3),
     &                  siteup(1,oe,mu),normalsites(1))
         
c -- c --------------------------------------------------------------- 
         if (dimproc(mu).gt.1) then
            CALL MPI_WAIT(send_request,send_status,ierror)       
            CALL MPI_WAIT(rec_request,rec_status,ierror)       
         endif

         call mpi_barrier(comm,ierror)

*---------------------------------------------------------------------- 
*
*  4a) transfer the bottom staple (7\dag * 6\dag * 5) in nu direction
*   b) compute 4 * (3\dag * 2\dag)
*   c) wait for completing the transfer
*

c -- a -----------------------------------------------------------------
         if (dimproc(nu).gt.1) then
            CALL MPI_ISSEND(helpbt2(1,1),1,
     1           u_VECTOR_upTYPE(nu,oe),neigh(nu),
     2           2000+oe+2*mu+8*nu,comm,
     3           send_request,ierror)       
            CALL MPI_IRECV(helpbt2(1,nsiteshalf+1),1,
     1           u_CONT_TYPE,neigh(-nu),
     2           2000+oe+2*mu+8*nu,comm,
     3           rec_request,ierror)       
         endif

c -- b -----------------------------------------------------------------
         call multallnn(u(1,1,oe,nu),helpup1(1,1),helpup2(1,1),
     &                  maxarray,ownhalfvol(3),
     &                  siteup(1,eo,mu),normalsites(1))

c -- c -----------------------------------------------------------------
         if (dimproc(nu).gt.1) then
            CALL MPI_WAIT(send_request,send_status,ierror)       
            CALL MPI_WAIT(rec_request,rec_status,ierror)       
         endif

         call mpi_barrier(comm,ierror)

c ------------------------------------------------------------------------
*
*  5) add both staples together
*
         do 50 site=1,ownhalfvol(3)
            newsite=sitedn(site,eo,nu)
            v(1,site,eo,mu)=v(1,site,eo,mu)+helpup2(1,site)+
     &                        helpbt2(1,newsite)
            v(2,site,eo,mu)=v(2,site,eo,mu)+helpup2(2,site)+
     &                        helpbt2(2,newsite)
50       continue




 100  continue


      return

      end          

