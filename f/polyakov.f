********************************************************************
********************************************************************
*                                                                  *
* Polyakov loop correlators by Peter Petreczky                     *
* May 1998                                                         *
*                                                                  *
********************************************************************
********************************************************************

      subroutine polyakov(meas,mstart)

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
      INCLUDE 'const.inc'
*
* local variables
*
      integer sites, eo, slice,start,a_count,x0,z,meas,mstart
      real    a_local_sum(0:extension(3)-1),  
     &        a_global_sum(0:extension(3)-1), 
     &        a_slice(0:n_lattice(3)-1),      
     &        polyak(0:n_lattice(3)-1)

*
* sum the polyakov line  on the local lattice on each z-slice
*
      a_local_sum = 0.0
     
*
* The Polyakov loop in the effective theory is defined by 
* P(x,y,z)=Tr exp(i A_0^a t^a)
* For SU(2) P^dagger=P and P=cos(A_0), A_0 = || A_0^a ||

      do eo = even, odd
         start=0
         do slice = 0,extension(3)-1
            do sites= 1, ownhalfvol(2)
               
               a_local_sum(slice)   = a_local_sum(slice)   + 
     &         cos(4*sqrt( a(1,sites+start,eo)**2 + 
     &                    a(2,sites+start,eo)**2 +
     &                    a(3,sites+start,eo)**2 ) )     
                             
            end do
            start = start + ownhalfvol(2)
         end do
      end do
*
* define the number of data that have to be send
* for the gauge field its for one direction !

      a_count=extension(3)
      
*
* sum the results on the roots of the timeslices
*
      call MPI_REDUCE(a_local_sum,a_global_sum,a_count,MPI_REAL,MPI_SUM,
     &                 xy_slice_root,xy_slice_comm)


*
* normalize the result on xy_slice_root
*
      if(my_xy_slice_rank==xy_slice_root)then
         a_global_sum = a_global_sum / (n_lattice(1)*n_lattice(2))
      end if
*
* gather the results from the xy_slice_roots on  z_column_root
*

      if(my_xy_slice_rank==xy_slice_root)then

         call MPI_GATHER(a_global_sum,a_count,MPI_REAL,
     &                   a_slice     ,a_count,MPI_REAL,
     &                   z_column_root,z_column_comm,ierror)
      endif

*
* calculate correlators according to
*
* polyak(z)=
* 1/n_lattice(3)*sum_{x=0}^n_lattice(3) polyak(x)*polyak(x+z)
* 

      if (myrank==root) then

         do z = 0, n_lattice(3)-1
            polyak(z)= 0.0
            do x0 = 0, n_lattice(3)-1
               polyak(z) = polyak(z) + 
     &         a_slice(x0)*a_slice( mod(x0+z,n_lattice(3)) )
            end do
         end do

         polyak = polyak / n_lattice(3)
   
         open(unit=polyak_unit, file=polyak_file, status='unknown',
     &                                        position='append' )
         write(polyak_unit, '(i5.5)')mstart+meas
         write(polyak_unit, '(e23.16)')(sum(a_slice)/n_lattice(3))**2
         do z=0,n_lattice(3)-1
          write(polyak_unit, '(e23.16)')polyak(z)
         enddo
         close(polyak_unit)


      end if


      end subroutine polyakov






