*********************************************************************
*********************************************************************
**  Subroutine for calculation of A_0 and A-i correlators          **
**  Written for an SU(2) adjoint Higgs model by Manfred Oevers     **
**  December 1997                                                  **
**  Modified by P.P on 09.04.98                                    **
**  input : meas, mstart                                           **
*********************************************************************
*********************************************************************

      subroutine correlator(meas,mstart)

      implicit none
**
** including all global MESSAGE PASSING features and variables 
**
      INCLUDE 'mpif.h'
      INCLUDE 'parallel_parameter.inc'
      INCLUDE 'parallel_mpi_types.inc'
**
** including global variables, parameters and common blocks
**
      INCLUDE 'parameter.inc'
      INCLUDE 'input_parameter.inc'
      INCLUDE 'gaugefield.inc'
      INCLUDE 'times.inc'
      INCLUDE 'filenames.inc'
      INCLUDE 'const.inc'
**
** local variables
**
      integer sites, eo, slice,start,a_count,u_count,x0,z,meas,mstart
      real    a_local_sum(nmu,0:extension(3)-1),  
     &        a_global_sum(nmu,0:extension(3)-1), 
     &        a_slice(nmu,0:n_lattice(3)-1),      
     &        corr_a(0:n_lattice(3)-1), corr_u(0:n_lattice(3)-1),
     &        corr_3(0:n_lattice(3)-1)            
      complex u_local_sum(nc,0:extension(3)-1,nmu),
     &        u_global_sum(nc,0:extension(3)-1,nmu),
     &        u_slice(nc,0:n_lattice(3)-1,nmu)

**
** sum the higgsfield on the local lattice on each timeslice
**
      a_local_sum = 0.0
      u_local_sum = cmplx(0.0,0.0)
**
** A_mu=1/2i * (U_mu - U_mu^\dag) in our notation we then have:
** U_mu=(a,b)   ===>    U_mu^\dag=(a*,-b)
** ====>        A_mu=(Im(a),-ib) !
**

      do eo = even, odd
         start=0
         do slice = 0,extension(3)-1
            do sites= 1, ownhalfvol(2)
** higgs field               
               a_local_sum(:,slice)   = a_local_sum(:,slice)   + 
     &                                          a_f(:,start + sites,eo)
** gauge field 1-direction
               u_local_sum(1,slice,1) = u_local_sum(1,slice,1) + 
     &                                    aimag(u_f(1,start + sites,eo,1))
               u_local_sum(2,slice,1) = u_local_sum(2,slice,1)  
     &                                    -zi * u_f(2,start + sites,eo,1)
** gauge field 2-direction
               u_local_sum(1,slice,2) = u_local_sum(1,slice,2) + 
     &                                    aimag(u_f(1,start + sites,eo,2))
               u_local_sum(2,slice,2) = u_local_sum(2,slice,2)  
     &                                    -zi * u_f(2,start + sites,eo,2)
**  gauge field 3-direction
               u_local_sum(1,slice,3) = u_local_sum(1,slice,3) + 
     &                                    aimag(u_f(1,start + sites,eo,3))
               u_local_sum(2,slice,3) = u_local_sum(2,slice,3)  
     &                                    -zi * u_f(2,start + sites,eo,3)
                                                
            end do
            start = start + ownhalfvol(2)
         end do
      end do
**
** define the number of data that have to be send
** for the gauge field its for one direction !

      a_count=nmu * extension(3)
      u_count=nc  * extension(3)
**
** sum the results on the roots of the timeslices
**
      call MPI_REDUCE(a_local_sum,a_global_sum,a_count,MPI_REAL,MPI_SUM,
     &                 xy_slice_root,xy_slice_comm)

      call MPI_REDUCE(u_local_sum,u_global_sum,u_count*3,
     &                        MPI_COMPLEX,MPI_SUM,xy_slice_root,xy_slice_comm)
**
** normalize the result on xy_slice_root
**
      if(my_xy_slice_rank==xy_slice_root)then
         a_global_sum = a_global_sum / (n_lattice(1)*n_lattice(2))
         u_global_sum = u_global_sum / (n_lattice(1)*n_lattice(2))
      end if
**
** gather the results from the xy_slice_roots on  z_column_root
**

      if(my_xy_slice_rank==xy_slice_root)then

         call MPI_GATHER(a_global_sum,a_count,MPI_REAL,
     &                   a_slice     ,a_count,MPI_REAL,
     &                   z_column_root,z_column_comm,ierror)
**
** for the gaugefield the directions have to be send separately,
** otherwise one would have the wrong order on z_column_root
**
         call MPI_GATHER(u_global_sum(1,0,1),u_count,MPI_COMPLEX,
     &                   u_slice(1,0,1)     ,u_count,MPI_COMPLEX,
     &                   z_column_root,z_column_comm,ierror)
         call MPI_GATHER(u_global_sum(1,0,2),u_count,MPI_COMPLEX,
     &                   u_slice(1,0,2)     ,u_count,MPI_COMPLEX,
     &                   z_column_root,z_column_comm,ierror)
         call MPI_GATHER(u_global_sum(1,0,3),u_count,MPI_COMPLEX,
     &                   u_slice(1,0,3)     ,u_count,MPI_COMPLEX,
     &                   z_column_root,z_column_comm,ierror)

         
      end if
**
** calculate correlators according to
**
** corr_a(z)=
** 1/n_lattice(3)*sum_{x=0}^n_lattice(3) a(x) \dot a(x+z)
** corr_u(z)=
** 1/n_lattice(3)*sum_{x=0}^n_lattice(3)sum_{mu=1,2}1/2*TR(U_mu(x)*U_mu(x+z))
**
** the trace of the product of two SU(2) matrices U1=(a,b) and U2=(c,d) is
** 1/2*TR(U1*U2)= Re(a*c+b*d^\dag)
** 

      if (myrank==root) then
         do z = 0, n_lattice(3)-1
            corr_a(z)= 0.0
            corr_u(z)= 0.0
            corr_3(z)= 0.0 
            do x0 = 0, n_lattice(3)-1
** higgs field
               corr_a(z)= corr_a(z) + 
     &           dot_product(a_slice(:,x0),a_slice(:,mod(x0+z,n_lattice(3))))
** gaugefield
               corr_u(z) = corr_u(z) + real(
     &           u_slice(1,x0,1) *       u_slice(1,mod(x0+z,n_lattice(3)),1)  +
     &           u_slice(2,x0,1) * conjg(u_slice(2,mod(x0+z,n_lattice(3)),1)) +
     &           u_slice(1,x0,2) *       u_slice(1,mod(x0+z,n_lattice(3)),2)  +
     &           u_slice(2,x0,2) * conjg(u_slice(2,mod(x0+z,n_lattice(3)),2)) )
               corr_3(z) = corr_3(z) + real(
     &           u_slice(1,x0,3) * u_slice(1,mod(x0+z,n_lattice(3)),3)  +
     &           u_slice(2,x0,3) * conjg(u_slice(2,mod(x0+z,n_lattice(3)),3)) ) 

            end do
         end do

** normalize correlation function

         corr_a = corr_a / n_lattice(3)
         corr_u = corr_u / (n_lattice(3) * 2)
         corr_3 = corr_3 / n_lattice(3)

** output results to the file

         open(unit=corr_unit, file=corr_file, status='unknown',position='append' )
         write(corr_unit, '(i3.3)')mstart+meas
         do z=0,n_lattice(3)-1
          write(corr_unit, '(2(e23.16,1x))')corr_a(z),corr_u(z)
         enddo
         close(corr_unit)


         open(unit=checkp_unit, file=checkp_file, status='unknown',position='append' )
         write(checkp_unit, '(i3.3)')mstart+meas
         do z=0,n_lattice(3)-1
          write(checkp_unit, '(e23.16)')corr_3(z)
         enddo
         close(checkp_unit)

      end if


      end subroutine correlator






