******************************************************************
** subroutine for updating gauge fields.                        **
** input  : measflag                                            **
** output : accept                                              **
******************************************************************

      subroutine agauge(measflag,accept)

      implicit none

******************************************************************* 

**
**  including all global MESSAGE PASSING features and variables 
**
      INCLUDE 'mpif.h'
      INCLUDE 'parallel_parameter.inc'
      INCLUDE 'parallel_mpi_types.inc'

**
**  including global variables, parameters and common blocks
**
      INCLUDE 'parameter.inc'
      INCLUDE 'input_parameter.inc'
      INCLUDE 'pointer.inc'
      INCLUDE 'gaugefield.inc'
      INCLUDE 'measurement.inc'


**
**   local VARIABLES
**
      integer eo,oe,mu,i,site,measflag,acc,ihit,nhit
      real    accept,tmp1(ownhalfvol(3)),tmp2(ownhalfvol(3)),
     &        sa_old(ownhalfvol(3)),sa_new(ownhalfvol(3)),rx,
     &        etmp,acc1,norm(ownhalfvol(3))
      complex u_tmp1(2,ownhalfvol(3)),u_tmp2(2,ownhalfvol(3)),
     &        u_tmp3(2,ownhalfvol(3)),u_tmp4(2,ownhalfvol(3)),
     &        h(2,ownhalfvol(3),0:1,3),nstaple(2,ownhalfvol(3))
      parameter(nhit=3)
***********************************************************************


**
**  the whole sweep is done in six parts 
**  one part for even-odd and every mu-direction
**  because after one iteration of the heat bath the links for the next
**  staples of another mu direction have changed
**  part I : pure gauge update   
**  part II : propose new configuration and accept/reject step for
**  gauge-higgs interaction,  these
**  will be done 3 times to get higher acceptance rates
**  part III : calculate the average plaq. and acceptance rates for
**  higgs-gauge interaction 
**  part IV  : normalize measurements
**
      accept=0.0
      acc1=0.0
      plaq0=0.0

      do 10 eo=0,1

        oe=mod(eo+1,2)

         do 20 mu=1,3

******************************************************* 
**                Pure gauge update                  **
*******************************************************

**
**   calculate the staple
**
            call staples(eo,oe,mu)
**
**   normalize the staple v to su(2) matrix nstaple, calculate 
**   calculate 1/(beta tmp1)=norm
**   and do this  for all sites i
**
           do site=1,ownhalfvol(3)

            tmp1(site)=real(v(1,site,eo,mu)*conjg(v(1,site,eo,mu))+
     &           v(2,site,eo,mu)*conjg(v(2,site,eo,mu)))
            tmp1(site)=1.0/sqrt(tmp1(site))

            nstaple(1,site)=tmp1(site)*v(1,site,eo,mu)
            nstaple(2,site)=tmp1(site)*v(2,site,eo,mu)

            norm(site)=tmp1(site)/beta

            enddo

*****************************************************************
**  Generate new configuration according the pure gauge action **
**  and accept/reject step for the gauge-Higgs interaction     **
**  this is reapeted ihit=3 times                              **
***************************************************************** 

**  start loop over hits

        do ihit=1,nhit
 
**  calculate U*nstaple

           do site=1,ownhalfvol(3)
            tmp2(site)=real(u(1,site,eo,mu)*nstaple(1,site)-
     &          u(2,site,eo,mu)*conjg(nstaple(2,site)))
           enddo

**   create SU(2) matrices with probability distribution
**   exp{ 1/(2*a) Tr V },    with V = u_tmp2 and a = norm.
**   with Kennedy-Pendleton heatbath.
**   input : tmp1, tmp2
**   output : u_tmp2, acc
**

          call heat_kpi(norm,u_tmp2,tmp2,acc)
          accept=accept+real(acc)
**  
**   new U's are:  U = u_tmp2 * nstaple^dagger
**
      call multallnd(u_tmp2(1,1),nstaple(1,1),u_tmp3(1,1),maxarray,
     &           ownhalfvol(3),normalsites(1),normalsites(1))

      do site=1,ownhalfvol(3)
       h(1,site,eo,mu)=u_tmp3(1,site)
       h(2,site,eo,mu)=u_tmp3(2,site)
      enddo

******     goto 9

**
**  transfer a(x+mu) from the neighbouring processor
**

         if (dimproc(mu).gt.1) then
            CALL MPI_ISSEND(a(1,1,oe),1,a_VECTOR_dnTYPE(mu,oe),
     1           neigh(-mu),1000+oe+2*mu,comm,
     2           send_request,ierror)
            CALL MPI_IRECV(a(1,nsiteshalf+1,oe),1,a_CONT_TYPE,
     1           neigh(mu),1000+oe+2*mu,comm,
     2           rec_request,ierror)
         endif

**  calculate a(x) U_mu(x) for new and old U configuration

         call multallan(a(1,1,eo),u(1,1,eo,mu),u_tmp1(1,1),maxarray,
     &                   ownhalfvol(3),normalsites(1),normalsites(1))
         call multallan(a(1,1,eo),h(1,1,eo,mu),u_tmp2(1,1),maxarray,
     &                   ownhalfvol(3),normalsites(1),normalsites(1))

**  wait until all transfer is completed

         if (dimproc(mu).gt.1) then
            CALL MPI_WAIT(send_request,send_status,ierror)
            CALL MPI_WAIT(rec_request,rec_status,ierror)
         endif

         call mpi_barrier(comm,ierror)

**  calculate the old adjoint Higgs action

         call multallad(a(1,1,oe),u(1,1,eo,mu),u_tmp3(1,1),maxarray,        
     &                   ownhalfvol(3),siteup(1,eo,mu),normalsites(1))

         call multallnn(u_tmp1(1,1),u_tmp3(1,1),u_tmp4(1,1),maxarray,
     &                   ownhalfvol(3),normalsites(1),normalsites(1))

         sa_old=beta_a*real(u_tmp4(1,:))

**  calculate the new adjoint Higgs action

         call multallad(a(1,1,oe),h(1,1,eo,mu),u_tmp3(1,1),maxarray,
     &                   ownhalfvol(3),siteup(1,eo,mu),normalsites(1))

         call multallnn(u_tmp2(1,1),u_tmp3(1,1),u_tmp4(1,1),maxarray,
     &                   ownhalfvol(3),normalsites(1),normalsites(1))


         sa_new=beta_a*real(u_tmp4(1,:))

**  generate random number for accept/reject

         

**  accept reject step
         
         do site=1,ownhalfvol(3)
           call random_number(rx)
           if (exp(sa_old(site)-sa_new(site)).ge.rx) then
           u(1,site,eo,mu)=h(1,site,eo,mu)
           u(2,site,eo,mu)=h(2,site,eo,mu)
           acc1=acc1+1.0
           endif
         enddo

**   end loop over hits
      
      enddo

*************************************************
**      Calculate the average plaquette        **
*************************************************

9     if(measflag.eq.1)then
      call multallnn(u(1,1,eo,mu),v(1,1,eo,mu),u_tmp1(1,1),maxarray,
     &           ownhalfvol(3),normalsites(1),normalsites(1))
      tmp1=real(u_tmp1(1,:))
      etmp=sum(tmp1)
      plaq0=plaq0+etmp
      endif

**  calling a barrier to make sure that no staple will be computed
**  untill not all links have been updated
   
      call MPI_BARRIER(comm,ierror)
       
**   end loop over even-odd and mu

20    continue
10    continue

*******************************************************  
**             Normalize measurements                **
******************************************************* 

      if (measflag.eq.1) then
      plaq0=plaq0/(12.0*ownvol(3))
      accept=accept/(3*ownvol(3))
      acc1=acc1/(12.0*ownvol(3))/nhit
      accept=accept*acc1
      endif


      return

      end
















