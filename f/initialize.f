      subroutine initialize

***********************************************************
* initializes linkvariables, indexvectors, topology, etc. *
***********************************************************

      implicit none
*
* including MESSAGE PASSING variables
*
      INCLUDE 'mpif.h'
      INCLUDE 'parallel_parameter.inc'
      INCLUDE 'parallel_mpi_types.inc'

*
* including global variables parameters and common blocks
*
      INCLUDE 'input_parameter.inc'
      INCLUDE 'times.inc'
*
* call the MPI Program to start message passing features
*
      call MPI_INIT(ierror)
*
* time for starting the program
*
      firsttime(1) = MPI_WTIME()

* time for starting initialisation

      firsttime(2) = MPI_WTIME()
*
* initialisation of parameters including lattice size etc.
*
      call paraminit
*
* Then build a virtual TOPOLOGY for the 4-D lattice with the given 
* parameters and build the appropriate types for Message passing
*
      call topology
      call typeconstruct
*
* initialisaation of pointers
*
      call adresses
*
* initialisation of the measured quantities and outfiles
*
      call measinit
*
* initialisation of the links
*
      call confget
*
* measure time
*
      lasttime(2)=MPI_WTIME()
      needtime(2)=lasttime(2)-firsttime(2)

      if (myrank.eq.0) then 
         write(*,*)'Time needed for initialisation ', needtime(2)
      endif


      end subroutine initialize





