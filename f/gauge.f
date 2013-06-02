****************************************************************
*                                                              *
*     local or global gauge of all the link variables          *
*     u`(x,mu) = g(x) * u(x,mu) * g(x+mu)\dag  (g\dag=g**-1)   *
*     and femionfields                                         *
*     phi'(x) = g(x) * phi(x)                                  *
*                                                              *
****************************************************************


      subroutine gauge(subgauge)


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
      INCLUDE 'gaugefield.inc'
      INCLUDE 'gauge.inc'
      INCLUDE 'pointer.inc'
      INCLUDE 'momenta.inc'
      

*
*local VARIABLES
*
      integer eo,oe,mu,subgauge,i,j,site,newsite,count
      complex help2(3,3),help(3,3),helpphi(nc)


***********************************************************************

      if (subgauge.eq.1) then
*
*global gauge
*
*only for ROOT
*first create a random matrix
*
         if (myrank.eq.root) then
            help(1,1)=cmplx(ranf(),ranf())
            help(2,1)=cmplx(ranf(),ranf())
            help(3,1)=cmplx(ranf(),ranf())
            help(1,2)=cmplx(ranf(),ranf())
            help(2,2)=cmplx(ranf(),ranf())
            help(3,2)=cmplx(ranf(),ranf())
            help(1,3)=cmplx(ranf(),ranf())
            help(2,3)=cmplx(ranf(),ranf())
            help(3,3)=cmplx(ranf(),ranf())
*
*unitarise it and broadcast it to all the other processes
*
            call unit(help(1,1))
         endif
         CALL MPI_BCAST(help(1,1),1,matrixtype,root,comm,ierror)

*
*then all processes gauge the links
*
         do 20 mu=1,nmu
            do 20 eo=0,1
               do 20 site=1,ownhalfvol(4)
                  call multnnd(help(1,1),u(1,1,site,eo,mu),
     1                 help(1,1),help2(1,1))
                  u(1,1,site,eo,mu)=help2(1,1)
                  u(2,1,site,eo,mu)=help2(2,1)
                  u(3,1,site,eo,mu)=help2(3,1)
                  u(1,2,site,eo,mu)=help2(1,2)
                  u(2,2,site,eo,mu)=help2(2,2)
                  u(3,2,site,eo,mu)=help2(3,2)
                  u(1,3,site,eo,mu)=help2(1,3)
                  u(2,3,site,eo,mu)=help2(2,3)
                  u(3,3,site,eo,mu)=help2(3,3)
 20      continue
*
* then all processes gauge the fermionfields
*
         do eo = even,odd
            do site = 1, ownhalfvol(4)
               do mu = 1, nmu
                  call gauge_phi(help(1,1),phi(1,mu,site,eo),helpphi(1))
                  phi(:,mu,site,eo)=helpphi(:)
               end do
            end do
         end do


      else
*
*local gauge
*


*       
*compute random matrices      
*       
         do 30 eo=0,1
            do 30 site=1,ownhalfvol(4)
               g(1,1,site,eo)=cmplx(ranf(),ranf())
               g(2,1,site,eo)=cmplx(ranf(),ranf())
               g(1,2,site,eo)=cmplx(ranf(),ranf())
               g(2,2,site,eo)=cmplx(ranf(),ranf())
               g(1,3,site,eo)=cmplx(ranf(),ranf())
               g(2,3,site,eo)=cmplx(ranf(),ranf())
*
*and unitarize them
*         
               call unit(g(1,1,site,eo))

 30      continue


*
*and transfer the halo to all the other processes which need them
*the necessary points are always the nearest neighbours to the
*gauged link
*so one transfer in each direction (if number of processes in that 
*direction is greater than one)
*         
*transfer first the x/odd direction (it will be needed for x/even 
*computation)
*then do loop : 1) send/receive the next direction : x/even, y/odd,
*                  y/even,...
*               2) gauge the links of the sent before : x/even,
*                  x/odd, y/even,...
*               3) end the transfer
*
*mat_VECTOR_TYPE(direction,even/odd)
*request/status(direction to send, even/odd to send, dummi)
*


         
*x-direction / odd
         if (dimproc(1).gt.1) then
            CALL MPI_ISSEND(g(1,1,1,1),1,mat_VECTOR_dnTYPE(1,1),
     1           neigh(-1),101,comm,send_request,ierror)       
            CALL MPI_IRECV(g(1,1,nsiteshalf+1,1),1,mat_CONT_TYPE,
     1           neigh(1),101,comm,rec_request,ierror)       
            CALL MPI_WAIT(send_request,send_status,ierror)       
            CALL MPI_WAIT(rec_request,rec_status,ierror)       
         endif

         call mpi_barrier(comm,ierror)

         do 50 mu=1,nmu
            do 50 eo=0,1

               oe=mod((eo+1),2)
               if ((mu+eo).lt.5) then
                  if (dimproc(mu+eo).gt.1) then
                     CALL MPI_ISSEND(g(1,1,1,eo),1,
     1                    mat_VECTOR_dnTYPE(mu+eo,eo),neigh(-(mu+eo)),
     2                    (100+2*mu+eo),comm,send_request,ierror)
                     CALL MPI_IRECV(g(1,1,nsiteshalf+1,eo),1,
     1                    mat_CONT_TYPE,neigh(mu+eo),
     2                    (100+2*mu+eo),comm,rec_request,ierror)
                  endif
               endif

               do 55 site=1,ownhalfvol(4)
                  newsite=siteup(site,eo,mu)
                  call multnnd(g(1,1,site,eo),u(1,1,site,eo,mu),
     1                 g(1,1,newsite,oe),help(1,1))
                  u(1,1,site,eo,mu)=help(1,1)
                  u(2,1,site,eo,mu)=help(2,1)
                  u(3,1,site,eo,mu)=help(3,1)
                  u(1,2,site,eo,mu)=help(1,2)
                  u(2,2,site,eo,mu)=help(2,2)
                  u(3,2,site,eo,mu)=help(3,2)
                  u(1,3,site,eo,mu)=help(1,3)
                  u(2,3,site,eo,mu)=help(2,3)
                  u(3,3,site,eo,mu)=help(3,3)
*
* then all processes gauge the fermionfields
*
                  call gauge_phi(g(1,1,site,eo),phi(1,mu,site,eo),helpphi(1))
                  phi(:,mu,site,eo)=helpphi(:)
 55            continue

               if ((mu+eo).lt.5) then
                  if (dimproc(mu+eo).gt.1) then
                     CALL MPI_WAIT(send_request,send_status,ierror)       
                     CALL MPI_WAIT(rec_request,rec_status,ierror)       
                  endif
               endif

         call mpi_barrier(comm,ierror)

 50      continue
                  
      endif



      end subroutine gauge






