********************************************************************
********************************************************************
*                                                                  *
* Written by Peter Schmidt for an SU(3) pure gauge code and a      *
* Wilson code to invert the fermion matrix                         *
*                                                                  *
* Adapted for an SU(2) adjoint Higgs model by Manfred Oevers       *
* November 1997                                                    *
*                                                                  *
********************************************************************
********************************************************************
********************************************************************
*                                                                  *
*       routine for the generation of all the linkvariables        *
*                                                                  *
********************************************************************

      subroutine confget

      implicit none
*
*include MESSAGE PASSING features and variables
*
      INCLUDE 'mpif.h'
      INCLUDE 'parallel_parameter.inc'

*
*include global variables, parameters and common blocks
*
      INCLUDE 'input_parameter.inc'
      INCLUDE 'filenames.inc'
      INCLUDE 'times.inc'

*
*local variables
*



***********************************************************************
     

*
* 0) randomly chosen gaugefields and higgsfields 
* 1) same random gauge and higgsfield on all links/sites
* 2) identity for all gauge and higgsfields 
* 3) read the old configuration
*
      firsttime(8)=MPI_WTIME()

      if (init.eq.0) call randomconf

      if (init.eq.1) call onerandconf

      if (init.eq.2) call unitconf

      if (init.eq.3) call oldconf

      if (init.eq.4) call testconf

      if ((init.gt.4).OR.(init.lt.0)) then
         if (myrank.eq.root) then
             write(*,*)'break : wrong gauge initialisation'
             call flush(6)
             call abort()
         endif
      endif

      lasttime(8)=MPI_WTIME()
      needtime(8)=lasttime(8)-firsttime(8)

      if (myrank.eq.0) then 
         write(*,*)'Time needed in confget: ', needtime(8)
      endif


      end subroutine confget

****************************************************************
*                                                              *
*     generates random SU(3) matrices for  all links and       *
*     a random higgs field of length one for  each site        *
*                                                              *
****************************************************************

      subroutine randomconf

      implicit none
*
*including all global MESSAGE PASSING features and variables 
*
      INCLUDE 'mpif.h'
      INCLUDE 'parallel_parameter.inc'
*
*including global variables, parameters and common blocks
*
      INCLUDE 'parameter.inc'
*      INCLUDE 'input_parameter.inc'
      INCLUDE 'gaugefield.inc'
      INCLUDE 'const.inc'
*
*local VARIABLES
*
      complex    uhelp(2)
      real       ahelp(3)
      real       xx(4)
      real       s1,s2,s3,s4,s5,r
      integer site,eo,mu


*
* first for all sites, directions, even/odd 
* a random su(2) matrix is created 
*
      do mu=1,nmu
         do eo=0,1
            do site=1,ownhalfvol(3)
               call random_number(r)
               xx(1) = 2.0 * r - 1.0
               s1=sqrt(1.0-xx(1)*xx(1))
               call random_number(r)
               s2 = 2.0 * r - 1.0
               xx(2) = s1*s2
               s3 = sqrt(1.0 - s2*s2)
               s4 = s1*s3
               call random_number(r)
               s5=twopi * r
               xx(4)=s4*cos(s5)
               xx(3)=s4*sin(s5)

               u(1,site,eo,mu) = cmplx(xx(1),xx(2))
               u(2,site,eo,mu) = cmplx(xx(3),xx(4))
            end do 
         end do
      end do
*
* then a random higgs field is created
*
      do site =1, ownhalfvol(3)
         do eo = even, odd
            call random_number(r)
            a(3,site,eo) = 2.0 * r - 1.0
            s1=sqrt(1.0 - a(3,site,eo)*a(3,site,eo))
            call random_number(r)
            s2=twopi * r
            a(1,site,eo)=s1*cos(s2)
            a(2,site,eo)=s1*sin(s3)
         end do
      end do

      end subroutine randomconf






****************************************************************
*                                                              *
*     generates one random SU(2) matrix and copies it to       *
*     all links and creates one higgs field of length one      *
*     pointing in a random direction                           *
*                                                              *
**************************************************************** 

      subroutine onerandconf

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
      INCLUDE 'const.inc'
*
*local VARIABLES
*
      complex    uhelp(2)
      real       ahelp(3)
      real       xx(4)
      real       s1,s2,s3,s4,s5,r
      integer site,eo,mu,i

* create random su(2) matrix via 
* U= cos(phi)*1 + i*sin(phi)*\sum_i n_i * \sigma_i
* with n a unitvector in R^3 on root

      if (myrank.eq.root) then

         call random_number(r)
         xx(1) = 2.0 * r - 1.0
         s1=sqrt(1.0-xx(1)*xx(1))
         call random_number(r)
         s2 = 2.0 * r - 1.0
         xx(2) = s1*s2
         s3 = sqrt(1.0 - s2*s2)
         s4 = s1*s3
         call random_number(r)
         s5=twopi * r
         xx(4)=s4*cos(s5)
         xx(3)=s4*sin(s5)

         uhelp(1)= cmplx(xx(1),xx(2))
         uhelp(2)= cmplx(xx(3),xx(4))

         s1=xx(1)**2 + xx(2)**2 + xx(3)**2 + xx(4)**2
         
         write(*,*) 'Length of S1=',s1
      endif

*
* then broadcast the matrix to all other processes
*
      CALL MPI_BCAST(uhelp(1),1,matrixtype,root,comm,ierror)

*
*then copy it to all the links
*
      do site=1,maxarray
      do eo=even,odd
      do mu=1,nmu
         u(:,site,eo,mu) = uhelp
      enddo
      enddo
      enddo


*
* create a higgsfield with random direction and length 1 on root
*
      if (myrank.eq.root) then

         call random_number(r)
         ahelp(3) = 2.0 * r - 1.0
         s1=sqrt(1.0 - ahelp(3)*ahelp(3))
         call random_number(r)
         s2=twopi * r
         ahelp(1)=s1*cos(s2)
         ahelp(2)=s1*sin(s3)

      endif

*
* then broadcast the matrix to all other processes
*
      CALL MPI_BCAST(ahelp(1),1,higgstype,root,comm,ierror)

*
*then copy it to all the links
*
      do site=1,maxarray
      do eo=even,odd
         a(:,site,eo) = ahelp
      enddo
      enddo

      end subroutine onerandconf

****************************************************************
*                                                              *
*     initializes u to the unitmatrix and a to sigma_3         *
*                                                              *
**************************************************************** 


      subroutine unitconf


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

*
*local variables
*
      integer site, eo, mu,i

      do site=1,maxarray
      do eo=even,odd
      do mu=1,nmu
         u(:,site,eo,mu) = (/(1.0, 0.0), (0.0, 0.0)/)
      enddo
      enddo
      enddo

      do site=1,maxarray
      do eo=even,odd
         a(:,site,eo) = (/0.0, 0.0, 1.0/)
      enddo
      enddo

      end subroutine unitconf

****************************************************************
*                                                              *
*     read the old configurations                              *
*                                                              *
**************************************************************** 


      subroutine oldconf


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
      INCLUDE 'times.inc'
      INCLUDE 'filenames.inc'
*     
*local VARIABLES
*
      integer  u_rl, a_rl, recnumber,coords(npdim),dest,
     &         mu, eo, x_proc, y_proc, z_proc, x_coord,y_coord,z_coord,
     &         slice, globalsite, localsite, offset,  sendsite, receivesite

      complex  u_helparray(nc,n_lattice(1)*n_lattice(2)/2)
      real     a_helparray(npdim,n_lattice(1)*n_lattice(2)/2)


***********************************************************************

      if(myrank==root)then

* record length for direct access, one global xy slice for fixed eo and mu 
* is read read in and distributed among the processors

         u_rl=16 * nc    * n_lattice(1) * n_lattice(2) / 2
         a_rl= 8 * npdim * n_lattice(1) * n_lattice(2) / 2
*
* open the files
*
         open(Unit=u_unit,file=u_conf_file,status='unknown',
     &                 access='direct',
     &                 form='unformatted',recl=u_rl)

         open(Unit=a_unit,file=a_conf_file,status='unknown',
     &                 access='direct',
     &                 form='unformatted',recl=a_rl)
      end if
      
* read in the gauge field

      recnumber=1
      do mu = 1, nmu
         do eo = even, odd
            do z_proc = 0, dimproc(3)-1
               do slice = 0, extension(3)-1
                  if (myrank==root) then

* read in xy global slice

                     read(Unit=u_unit,rec=recnumber,err=150) u_helparray
                     goto 200
 150                 call abort()
 200                 continue

*     reorder links in the force field for transfer
                     do y_proc = 0, dimproc(2) - 1
                        do x_proc = 0, dimproc(1) - 1
                           do y_coord = 0, extension(2) - 1
                              do x_coord = 0, ownhalfvol(1) -1

                                globalsite = 1 + x_proc * ownhalfvol(1) + x_coord +
     &                          (y_proc * extension(2) + y_coord)*ownhalfvol(1)*dimproc(1)

                                localsite  = 1 + x_coord + y_coord * ownhalfvol(1)

                                offset     = ( x_proc + y_proc * dimproc(1) ) *
     &                                        ownhalfvol(2)
                              
                                u_f(:,offset + localsite,eo,mu) =
     &                                                u_helparray(:,globalsite)
                             end do
                          end do
                       end do
                    end do
                 end if

* distribute xy slice among processors

                 do x_proc = 0, dimproc(1) - 1
                    do y_proc = 0, dimproc(2) - 1
                       coords(1)=x_proc
                       coords(2)=y_proc
                       coords(3)=z_proc
                       call MPI_CART_RANK(comm,coords,dest,ierror)
                       
                       if (myrank==root) then
                          
                          sendsite = ( x_proc + y_proc * dimproc(1) ) *
     &                                 ownhalfvol(2) + 1
                          call MPI_SEND(u_f(1,sendsite,eo,mu),1,
     &                                  load_u_type,dest,
     &                                  1000*slice+100*x_proc+10*y_proc+z_proc,
     &                                  comm,ierror)
                       end if
                       if (myrank==dest) then

                          receivesite = slice * ownhalfvol(2) + 1

                          call MPI_RECV(u(1,receivesite,eo,mu),1,
     &                                  load_u_type,root,
     &                                  1000*slice+100*x_proc+10*y_proc+z_proc,
     &                                  comm,save_rec_status,ierror)
                       end if
                       call MPI_BARRIER(comm,ierror1)
                    end do
                 end do
                 recnumber=recnumber+1
              end do
           end do
        end do
      end do

* read in the higgs field
      recnumber=1
      do eo=even, odd
         do z_proc = 0, dimproc(3)-1
            do slice = 0, extension(3)-1
               
               if (myrank==root) then

* read in global xy slice

                  read(Unit=a_unit,rec=recnumber,err=250) a_helparray
                  goto 300
 250              call abort()
 300              continue

* reorder higgsfield in the force field for transfer

                  do y_proc = 0, dimproc(2) - 1
                     do x_proc = 0, dimproc(1) - 1
                        do y_coord = 0, extension(2) - 1
                           do x_coord = 0, ownhalfvol(1) - 1

                              globalsite = 1 + x_proc * ownhalfvol(1) + x_coord +
     &                          (y_proc * extension(2) + y_coord)*ownhalfvol(1)*dimproc(1)

                              localsite  = 1 + x_coord + y_coord * ownhalfvol(1)
                           
                              offset     = ( x_proc + y_proc * dimproc(1) ) *
     &                                        ownhalfvol(2)
                              
                              a_f(:,offset + localsite,eo) =
     &                                                a_helparray(:,globalsite)
                           end do
                        end do
                     end do
                  end do
               end if

* distribute xy slice among processors

               do x_proc = 0, dimproc(1) - 1
                  do y_proc = 0, dimproc(2) - 1
                     coords(1)=x_proc
                     coords(2)=y_proc
                     coords(3)=z_proc
                     call MPI_CART_RANK(comm,coords,dest,ierror)

                     if (myrank==root)then

                        sendsite = ( x_proc + y_proc * dimproc(1) ) *
     &                               ownhalfvol(2) + 1

                        call MPI_SEND(a_f(1,sendsite,eo),1,
     &                                 load_a_type,dest,
     &                                 1000*slice+100*x_proc+10*y_proc+z_proc,
     &                                 comm,ierror)
                     end if
                     if (myrank==dest)then

                        receivesite = slice * ownhalfvol(2) + 1

                        call MPI_RECV(a(1,receivesite,eo),1,
     &                                load_a_type,root,
     &                                1000*slice+100*x_proc+10*y_proc+z_proc,
     &                                 comm,save_rec_status,ierror)
                     end if
                     call MPI_BARRIER(comm,ierror1)
                  end do
               end do

               recnumber=recnumber+1
            end do
         end do
      end do
*
* close the files for read in
*
      if (myrank==root)then
         close(20)
         close(21)
      end if

      return

      end subroutine oldconf

      subroutine testconf
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
      INCLUDE 'const.inc'
*
*local variables
*
      integer site, eo, mu, slice, start

      do mu=1,nmu
         do eo = even, odd
            start=0
            do slice = 0, extension(3)
               do site= 1, ownhalfvol(2)
                 u(1,start+site,eo,mu)=zi*(my_z_column_rank*extension(3)+slice)
                 u(2,start+site,eo,mu)=zi*(my_z_column_rank*extension(3)+slice)
               end do
               start=start+ownhalfvol(2)
            end do
         end do
         u(:,:,:,mu)=real(mu)*u(:,:,:,mu)
      end do


      do eo = even, odd
         start=0
         do slice = 0, extension(3)-1
            do site= 1, ownhalfvol(2)
               a(1,start+site,eo)=1./sqrt(3.)
               a(2,start+site,eo)=1./sqrt(3.)
               a(3,start+site,eo)=1./sqrt(3.)
            end do
            start=start+ownhalfvol(2)
         end do
      end do

      end subroutine testconf
