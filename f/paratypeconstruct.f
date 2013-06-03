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
*         initialisation of the types for MESSAGE PASSING          *
*                                                                  *
********************************************************************

      subroutine paratypeconstruct

*      USE mpi
      implicit none
      INCLUDE 'mpif.h'
*
*including all global MESSAGE PASSING variables 
*
      INCLUDE 'parallel_parameter.inc'
      INCLUDE 'parallel_mpi_types.inc'
*
*including global variables, parameters and common blocks
*
      INCLUDE 'parameter.inc'
      INCLUDE 'input_parameter.inc'
*
*definition of the VARIABLES used for MPI
*
      integer a_of_b_paramet(2),a_of_t_paramet(2)

      integer (kind=MPI_ADDRESS_KIND) :: a_of_d_paramet(2),
     &                                   length_block, ad2, ad1, extent

***********************************************************************
*
* first create three new TYPES to transfer the common block 
* parameter config and latconf for above data and then commit the 
* types to all processes
*
* paramet: 
* 7 X real for couplings and parameters 6 X integer for inorm,ntraj_flag,
* ntime,nsweep,imeas,imode 
* defined in include_parameter.inc

* config : 
* 4 integer for init,rand,save,confstart
* defined in include_parameter.inc

* latconf: 
* 9 integer for n_lattice,dimproc,extension 
* defined in parallel_parameter.inc
*
* a_of_b = array of blocklength
* a_of_d = array of displacements (in bytes)
* a_of_t = array of types of each block
*

* paramet

      length_block=1
      call MPI_TYPE_GET_EXTENT(MPI_REAL, length_block, extent, ierror)
      call mpi_barrier(comm,ierror)
      CALL MPI_GET_ADDRESS(beta,ad1,ierror)
      call mpi_barrier(comm,ierror)
      CALL MPI_GET_ADDRESS(inorm,ad2,ierror)
      CALL mpi_barrier(comm,ierror)

      a_of_b_paramet(1)=7
      a_of_b_paramet(2)=6

      a_of_d_paramet(1)=0
      a_of_d_paramet(2)=ad2-ad1

      a_of_t_paramet(1)=MPI_REAL
      a_of_t_paramet(2)=MPI_INTEGER

      call mpi_barrier(comm,ierror)
      CALL MPI_TYPE_CREATE_STRUCT(2,a_of_b_paramet,a_of_d_paramet,
     &                     a_of_t_paramet,paramettype,ierror)
      call mpi_barrier(comm,ierror)
      CALL MPI_TYPE_COMMIT(paramettype,ierror)

* config

      call mpi_barrier(comm,ierror)
      CALL MPI_TYPE_VECTOR(1,4,4,MPI_INTEGER,configtype,ierror)
      call mpi_barrier(comm,ierror)
      CALL MPI_TYPE_COMMIT(configtype,ierror)

* latconf

      call mpi_barrier(comm,ierror)
      CALL MPI_TYPE_VECTOR(1,9,9,MPI_INTEGER,latconftype,ierror)
      call mpi_barrier(comm,ierror)
      CALL MPI_TYPE_COMMIT(latconftype,ierror)
      call mpi_barrier(comm,ierror)


      return
     
      end












































































