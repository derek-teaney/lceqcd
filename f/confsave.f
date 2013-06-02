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
* routine to save the configuration of the link and higgs field    *
*                                                                  *
********************************************************************

      subroutine confsave

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
      integer  u_rl, a_rl, recnumber,coords(npdim),source,
     &         mu, eo, x_proc, y_proc, z_proc, x_coord,y_coord,z_coord,
     &         slice, globalsite, localsite, offset,  sendsite, receivesite,
     &         betaint, xint, yint

      complex  u_helparray(nc,n_lattice(1)*n_lattice(2)/2)
      real     a_helparray(npdim,n_lattice(1)*n_lattice(2)/2)


***********************************************************************


      firsttime(7)=MPI_WTIME()

      if(myrank==root)then

* record length for direct access, one global xy slice for fixed eo and mu 
* is read read in and distributed among the processors

         u_rl=16 * nc    * n_lattice(1) * n_lattice(2) / 2
         a_rl= 8 * npdim * n_lattice(1) * n_lattice(2) / 2
*
* open the files
*
         betaint    = nint(1000 * beta)
         xint  = nint(10000 * x)
         yint  = nint(10000 * y)

         write (u_conf_file, fmt = 1000) n_lattice(3),  
     &        betaint, xint, yint, (confstart+nsweep)
         write (a_conf_file, fmt = 2000) n_lattice(3),  
     &        betaint, xint, yint, (confstart+nsweep)

         open(Unit=u_unit,file=u_conf_file,status='unknown',access='direct',
     &                 form='unformatted',recl=u_rl)

         open(Unit=a_unit,file=a_conf_file,status='unknown',access='direct',
     &                 form='unformatted',recl=a_rl)
      end if
      
* write the gauge field

      recnumber=1
      do mu = 1, nmu
         do eo = even, odd
            do z_proc = 0, dimproc(3)-1
               do slice = 0, extension(3)-1

* get xy slice from processors

                  do x_proc = 0, dimproc(1) - 1
                     do y_proc = 0, dimproc(2) - 1
                        coords(1)=x_proc
                        coords(2)=y_proc
                        coords(3)=z_proc
                        call MPI_CART_RANK(comm,coords,source,ierror)

                        if (myrank==source) then

                           sendsite = slice * ownhalfvol(2) + 1

                           call MPI_SEND(u(1,sendsite,eo,mu),1,
     &                                  save_u_type,root,
     &                                  1000*slice+100*x_proc+10*y_proc+z_proc,
     &                                  comm,ierror)
                        end if

                        if (myrank==root) then
                           receivesite = ( x_proc + y_proc * dimproc(1) ) *
     &                          ownhalfvol(2) + 1

                           call MPI_RECV(u_f(1,receivesite,eo,mu),1,
     &                                  save_u_type,source,
     &                                  1000*slice+100*x_proc+10*y_proc+z_proc,
     &                                  comm,save_rec_status,ierror)
                        end if 
                        call MPI_BARRIER(comm,ierror1)
                     end do
                  end do

* reorder links in the force field for transfer

                  if (myrank==root) then

                     do y_proc = 0, dimproc(2) - 1
                        do x_proc = 0, dimproc(1) - 1
                           do y_coord = 0, extension(2)-1
                              do x_coord = 0, ownhalfvol(1)-1

                                globalsite=1 + x_proc * ownhalfvol(1) + x_coord
     &                        + (y_proc * extension(2) + y_coord)*ownhalfvol(1)*dimproc(1)

                                localsite=1 + x_coord + y_coord * ownhalfvol(1)

                                offset   = ( x_proc + y_proc * dimproc(1) ) *
     &                                        ownhalfvol(2)

                                u_helparray(:,globalsite)=
     &                                  u_f(:,offset + localsite,eo,mu)

                             end do
                          end do
                       end do
                    end do



* write global xy slice

                    write(Unit=u_unit,rec=recnumber,err=150) u_helparray

                    goto 200
 150                call abort()
 200                continue

                    recnumber=recnumber+1
                 end if
              end do
           end do
        end do
      end do

* write the higgs field

      recnumber=1
      do eo=even, odd
         do z_proc = 0, dimproc(3)-1
            do slice = 0, extension(3)-1
            
* distribute xy slice among processors

               do x_proc = 0, dimproc(1) - 1
                  do y_proc = 0, dimproc(2) - 1
                     coords(1)=x_proc
                     coords(2)=y_proc
                     coords(3)=z_proc
                     call MPI_CART_RANK(comm,coords,source,ierror)
                     
                     if(myrank==source) then
                     
                        sendsite = slice * ownhalfvol(2) + 1
                        call MPI_SEND(a(1,sendsite,eo),1,
     &                                 save_a_type,root,
     &                                 1000*slice+100*x_proc+10*y_proc+z_proc,
     &                                 comm,ierror)
                     end if
                     if (myrank==root) then
                        receivesite = ( x_proc + y_proc * dimproc(1) ) *
     &                                  ownhalfvol(2) + 1
                        call MPI_RECV(a_f(1,receivesite,eo),1,
     &                       save_a_type,source,
     &                       1000*slice+100*x_proc+10*y_proc+z_proc,
     &                       comm,save_rec_status,ierror)
                     end if
                     call MPI_BARRIER(comm,ierror1)
                  end do
               end do

               if (myrank==root) then

* reorder higgsfield in the force field for transfer

                  do y_proc = 0, dimproc(2) - 1
                     do x_proc = 0, dimproc(1) - 1
                        do y_coord = 0, extension(2) -1
                           do x_coord = 0, ownhalfvol(1) -1

                              globalsite = 1 + x_proc * ownhalfvol(1) + x_coord +
     &                          (y_proc * extension(2) + y_coord)*ownhalfvol(1)*dimproc(1)

                              localsite  = 1 + x_coord + y_coord * ownhalfvol(1)
                           
                              offset     = ( x_proc + y_proc * dimproc(1) ) *
     &                                        ownhalfvol(2)
                              


                             a_helparray(:,globalsite) =
     &                               a_f(:,offset + localsite,eo)
                           end do
                        end do
                     end do
                  end do

* write global xy slice

                  write(Unit=a_unit,rec=recnumber,err=250) a_helparray
                  goto 300
 250              call abort()
 300              continue
               end if

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

      lasttime(7)=MPI_WTIME()
      needtime(7)=lasttime(7)-firsttime(7)

      if (myrank.eq.0) then 
         write(*,*)'Time needed in confsave: ' ,needtime(7)
      endif

1000  format('u_conf_s', i2.2,  '_BE', i5.5, '_x', i5.5, '_y', i5.5,
     &                                        '_U',  i6.6)

2000  format('a_conf_s', i2.2,  '_BE', i5.5, '_x', i5.5, '_y', i5.5,
     &                                        '_U',  i6.6)


      return

      end subroutine confsave


















