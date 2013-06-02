***************************************************************
*              Gauge fixing with overrelaxation               * 
*  output: da1                                                *
***************************************************************
      subroutine or_gauge_fixing(da1)

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


      real  da1,da11 
      real  k(ownhalfvol(3)),theta(ownhalfvol(3)),theta_o(ownhalfvol(3)),
     &       tmp1(ownhalfvol(3)),tmp2(ownhalfvol(3))
      integer eo,oe,mu,ig,is
      complex Omega(2,maxarray), u_tmp1(2,ownhalfvol(3)),
     &        u_tmp2(2,ownhalfvol(3)),u_tmp3(2,ownhalfvol(3)),
     &        u_tmp(2,maxarray), Omega1(2,maxarray),
     &        Omega2(2,maxarray), Omega3(2,maxarray)

*           __         +
* calculate > ( U (x)+U (x-mu) ) = k(x) Omega(x)
*           --   mu    mu
* do it in 6 steps :
*
* 1) send U (x-mu) in +mu direction for mu=1,2,3
*          mu
* 2) calculate sum U  (x) over direction
*                   mu            
* 3) wait for end of the transfer for mu=1 and add  U_1(x-1)
*                  
* 4) wait for end of the transfer for mu=2 and add  U_2(x-2)
*
* 5) wait for end of the transfer for mu=3 and add  U_3(x-3)
*
* 6) sum( U_mu(x)+U_mu(x-mu)^dagg )
*

      da11=0.0

      do eo=0,1 
*     start loop over even-odd

      oe=mod(eo+1,2)

      u_tmp1=cmplx(0.0,0.0)
      u_tmp2=cmplx(0.0,0.0)



* 1) 
      if (dimproc(1).gt.1) then
       call MPI_ISSEND(u_f(1,1,oe,1),1,u_VECTOR_upTYPE(1,oe),
     &                 neigh(1),1000+oe+2*1,comm,send_request,ierror)
      endif

      if (dimproc(2).gt.1) then
       call MPI_ISSEND(u_f(1,1,oe,2),1,u_VECTOR_upTYPE(2,oe),
     &                 neigh(2),1000+oe+2*2,comm,send_request1,ierror1)
      endif
      
      if (dimproc(3).gt.1) then
       call MPI_ISSEND(u_f(1,1,oe,3),1,u_VECTOR_upTYPE(3,oe),
     &                 neigh(3),1000+oe+2*3,comm,send_request2,ierror2)
      endif
      
      if (dimproc(1).gt.1) then
       call MPI_IRECV(u_f(1,nsiteshalf+1,oe,1),1,u_CONT_TYPE,
     &               neigh(-1),1000+oe+2*1,comm,rec_request,ierror)
      endif

      if (dimproc(2).gt.1) then
       call MPI_IRECV(u_f(1,nsiteshalf+1,oe,2),1,u_CONT_TYPE,
     &                neigh(-2),1000+oe+2*2,comm,rec_request1,ierror1)
      endif


      if (dimproc(3).gt.1) then
       call MPI_IRECV(u_f(1,nsiteshalf+1,oe,3),1,u_CONT_TYPE,
     &               neigh(-3),1000+oe+2*3,comm,rec_request2,ierror2)
      endif

* 2)  
      do mu=1,2  
       do is=1,ownhalfvol(3)
        u_tmp1(1,is)=u_tmp1(1,is)+u_f(1,is,eo,mu)
        u_tmp1(2,is)=u_tmp1(2,is)+u_f(2,is,eo,mu)
        enddo
      enddo
      
      do is=1,ownhalfvol(3)
       u_tmp1(1,is)=u_tmp1(1,is) + alpha*u_f(1,is,eo,3)
       u_tmp1(2,is)=u_tmp1(2,is) + alpha*u_f(2,is,eo,3)
      enddo

* 3)
      if (dimproc(1).gt.1) then
       call MPI_WAIT(send_request,send_status,ierror)
       call MPI_WAIT(rec_request,rec_status,ierror)
      endif

* 4) 
      do is=1,ownhalfvol(3)
       u_tmp2(1,is)=u_tmp2(1,is)+u_f(1,sitedn(is,eo,1),oe,1)
       u_tmp2(2,is)=u_tmp2(2,is)+u_f(2,sitedn(is,eo,1),oe,1)
      enddo

* 5)
      if (dimproc(2).gt.1) then
       call MPI_WAIT(send_request1,send_status1,ierror1)
       call MPI_WAIT(rec_request1,rec_status1,ierror1)
      endif

      do is=1,ownhalfvol(3)
       u_tmp2(1,is)=u_tmp2(1,is)+u_f(1,sitedn(is,eo,2),oe,2)
       u_tmp2(2,is)=u_tmp2(2,is)+u_f(2,sitedn(is,eo,2),oe,2)
      enddo


* 6)
      if (dimproc(3).gt.1) then
       call MPI_WAIT(send_request2,send_status2,ierror2)
       call MPI_WAIT(rec_request2,rec_status2,ierror2)
      endif

      do is=1,ownhalfvol(3)
       u_tmp2(1,is)=u_tmp2(1,is) + alpha*u_f(1,sitedn(is,eo,3),oe,3)
       u_tmp2(2,is)=u_tmp2(2,is) + alpha*u_f(2,sitedn(is,eo,3),oe,3)
      enddo
* 7)  
      
      u_tmp3(1,:)=u_tmp1(1,:)+conjg(u_tmp2(1,:))
      u_tmp3(2,:)=u_tmp1(2,:)-u_tmp2(2,:)

 
*
* calculate |d_mu A_mu|**2
*

      if (eo.eq.1) then
         tmp1=real( (aimag(u_tmp1(1,:)-u_tmp2(1,:)))**2+
     &              u_tmp3(2,:)*conjg(u_tmp3(2,:)) )
        da11=sum(tmp1)
        da11=da11/real(ownhalfvol(3))
      endif
*
* calculate the determinant k(x)
*
      k(:)=real( u_tmp3(1,:)*conjg(u_tmp3(1,:))+
     &            u_tmp3(2,:)*conjg(u_tmp3(2,:)) )     

*
*  calculate theta , Omega(x)=1*cos(theta)+i sin(theta) a/|a| sigma
*  and perform the overralaxation step  theta->theta_o
*
      do is=1,ownhalfvol(3)
       Omega(1,is)=u_tmp3(1,is)/sqrt(k(is))
       Omega(2,is)=u_tmp3(2,is)/sqrt(k(is))
  

       theta(is)=acos(real(Omega(1,is)))

       if(theta(is).ne.0) then
        theta_o(is)=lambda*theta(is)
        tmp1(is)=cos(theta_o(is)) 
        tmp2(is)=sin(theta_o(is))/sin(theta(is))*aimag(Omega(1,is))
        Omega(1,is)=cmplx(tmp1(is),-tmp2(is))
        Omega(2,is)=-sin(theta_o(is))/sin(theta(is))*Omega(2,is)
       endif

      enddo
*                                                                +
*  U (x)  ->Omega(x) U  (x), for eo and U (x-mu)=U (x-mu) Omega(x), oe
*   mu                mu                 mu       mu
*                           +
*  A->Omega(x) A(x) Omega(x)  for eo and A-> A for oe
*
*  Do this in 6 steps:
*  0) copy Omega to Omega1, Omega2, Omega3
*  1) send Omega(x+mu) in -mu direction, for mu=1,2,3
*  2) perform a gauge transformation on A for even (eo) points
*  3) Omega1(x) U_1(x) for eo and U_1(x)*Omega(x+1) for oe
*  4) Omega2(x) U_2(x) for eo and U_2(x)*Omega(x+2) for oe 
*  5) Omega3(x) U_3(x) for eo and U_3(x)*Omega(x+3) for oe
*

* 0)
       Omega1=Omega
       Omega2=Omega
       Omega3=Omega

* 1)       
       if (dimproc(1).gt.1) then
       call MPI_ISSEND(Omega1(1,1),1,u_VECTOR_dnTYPE(1,eo),
     &                 neigh(-1),1000+eo+2*1,comm,send_request,ierror)
       endif

       if (dimproc(2).gt.1) then
       call MPI_ISSEND(Omega2(1,1),1,u_VECTOR_dnTYPE(2,eo),
     &                 neigh(-2),1000+eo+2*2,comm,send_request1,ierror1)
       endif


       if (dimproc(3).gt.1) then
       call MPI_ISSEND(Omega3(1,1),1,u_VECTOR_dnTYPE(3,eo),
     &                 neigh(-3),1000+eo+2*3,comm,send_request2,ierror2)
       endif


       if (dimproc(1).gt.1) then
       call MPI_IRECV(Omega1(1,nsiteshalf+1),1,u_CONT_TYPE,
     &                neigh(1),1000+eo+2*1,comm,rec_request,ierror)
       endif

       
       if (dimproc(2).gt.1) then
       call MPI_IRECV(Omega2(1,nsiteshalf+1),1,u_CONT_TYPE,
     &                neigh(2),1000+eo+2*2,comm,rec_request1,ierror1)
       endif
      

       if (dimproc(3).gt.1) then
       call MPI_IRECV(Omega3(1,nsiteshalf+1),1,u_CONT_TYPE,
     &                neigh(3),1000+eo+2*3,comm,rec_request2,ierror2)
       endif


* 2)

       call multallad(a_f(1,1,eo),Omega(1,1),u_tmp1(1,1),
     &                maxarray,ownhalfvol(3),normalsites(1),normalsites(1))
       call multallnn(Omega(1,1),u_tmp1(1,1),u_tmp2(1,1),
     &                maxarray,ownhalfvol(3),normalsites(1),normalsites(1))

       do is=1,ownhalfvol(3)
        a_f(1,is,eo)=aimag(u_tmp2(2,is))
        a_f(2,is,eo)=real(u_tmp2(2,is))
        a_f(3,is,eo)=aimag(u_tmp2(1,is))
       enddo


*3a)
       call multallnn(Omega1(1,1),u_f(1,1,eo,1),u_tmp1(1,1),
     &                maxarray,ownhalfvol(3),normalsites(1),normalsites(1))

       do is=1,ownhalfvol(3)
        u_f(1,is,eo,1)=u_tmp1(1,is)
        u_f(2,is,eo,1)=u_tmp1(2,is)
       enddo
*3b)       
       if (dimproc(1).gt.1) then
        call MPI_WAIT(send_request,send_status,ierror)
        call MPI_WAIT(rec_request,rec_status,ierror)
       endif

*3c)
       call multallnd(u_f(1,1,oe,1),Omega1(1,1),u_tmp2(1,1),
     &                maxarray,ownhalfvol(3),normalsites(1),siteup(1,oe,1))

       do is=1,ownhalfvol(3)
        u_f(1,is,oe,1)=u_tmp2(1,is)
        u_f(2,is,oe,1)=u_tmp2(2,is)
       enddo
       
*4a)
       call multallnn(Omega2(1,1),u_f(1,1,eo,2),u_tmp1(1,1),
     &                maxarray,ownhalfvol(3),normalsites(1),normalsites(1))

       do is=1,ownhalfvol(3)
        u_f(1,is,eo,2)=u_tmp1(1,is)
        u_f(2,is,eo,2)=u_tmp1(2,is)
       enddo
*4b)       
       if (dimproc(2).gt.1) then
        call MPI_WAIT(send_request1,send_status1,ierror1)
        call MPI_WAIT(rec_request1,rec_status1,ierror1)
       endif

*4c)
       call multallnd(u_f(1,1,oe,2),Omega2(1,1),u_tmp2(1,1),
     &                maxarray,ownhalfvol(3),normalsites(1),siteup(1,oe,2))

       do is=1,ownhalfvol(3)
        u_f(1,is,oe,2)=u_tmp2(1,is)
        u_f(2,is,oe,2)=u_tmp2(2,is)
       enddo



*5a)
       call multallnn(Omega3(1,1),u_f(1,1,eo,3),u_tmp1(1,1),
     &                maxarray,ownhalfvol(3),normalsites(1),normalsites(1))

       do is=1,ownhalfvol(3)
        u_f(1,is,eo,3)=u_tmp1(1,is)
        u_f(2,is,eo,3)=u_tmp1(2,is)
       enddo
*5b)       
       if (dimproc(3).gt.1) then
        call MPI_WAIT(send_request2,send_status2,ierror2)
        call MPI_WAIT(rec_request2,rec_status2,ierror2)
       endif

*5c)
       call multallnd(u_f(1,1,oe,3),Omega3(1,1),u_tmp2(1,1),
     &                maxarray,ownhalfvol(3),normalsites(1),siteup(1,oe,3))

       do is=1,ownhalfvol(3)
        u_f(1,is,oe,3)=u_tmp2(1,is)
        u_f(2,is,oe,3)=u_tmp2(2,is)
       enddo



       enddo 
* end loop over even-odd
       
       if (nproc.gt.1) then
        call MPI_ALLREDUCE(da11,da1,1,MPI_REAL,MPI_SUM,
     &                    comm,ierror)
        da1=da1/nproc
       else
        da1=da11
       endif

       return
      
       end













