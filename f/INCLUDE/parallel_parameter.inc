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
****************************************************************
*                                                              *
*     Include file for Message Passing variables               *
*                                                              *
*     the Parameter wit local information                      *
*                                                              *
****************************************************************


* the error handling
      INTEGER ierror,ierror1,ierror2

* the communicator
      INTEGER comm, xy_slice_comm, z_column_comm

* number and dimension of processes
      INTEGER nproc,npdim     

* rank of the root processor
      INTEGER root, xy_slice_root, z_column_root
 
* some Parameters
      Parameter (npdim=3, root=0, xy_slice_root=0, z_column_root=0)

* processes in each dimension
      INTEGER dimproc(npdim) 

* rank of the process 
      INTEGER myrank, my_xy_slice_rank, my_z_column_rank

* coordinates of the process
      INTEGER mycoords(npdim), my_xy_slice_coords(npdim-1), 
     &                         my_z_column_coords(npdim-2)

* neighbours of a process (  1,2,3=x,y,z ; -1,-2,-3=-x,-y,-z ; 0=root  )
      INTEGER neigh(-npdim:npdim) 

* left and right boundary
      INTEGER lowbound(npdim),upbound(npdim)

* extension in mu-direction
      INTEGER extension(npdim)

* volume of the local lattice
      INTEGER ownvol(npdim),ownhalfvol(npdim)

* global lattice size
      INTEGER n_lattice(npdim)

*
*      n_lattice(1:3)   physical dimensions of lattice, see paraminit
*
*      dimproc(1:3)     processors per dimension, see paraminit
*
*      extension(1:3)   extent of each direction for sublattice, so
*                       extension(1)*dimproc(1) = n_lattice(1).  
*                       Each extentsion should be greater than two
*                       and should be even.   see topology
*
*      ownvol(1:3)      Volume of individual lattice sites:
*                       ownvol(1) volume (i.e. length) of x direction 
*                       ownvol(2) volume (i.e. area)  of xy direction
*                       ownvol(3) volume of xyz direction
*                        
*      ownvolhalf(1:3)  Same as ownvol but 1/2,1/4,1/8 as much for each
*                       direction.
*
*      lowbound(1:3) and upbound(1:3)
*
*                        Array limits for a given processor site. One
*                        loops over the lattice via
*
*                        do i=lowbound(1), upbound(1) ...
*
*      comm             Basic cartesian communicator.  Created
*                       with MPI_CART_CREATE. To get the rank of a 
*                       process within this communicator group
* 
*                       CALL MPI_COMM(comm, rank, ierror)
* 
*                       This information is already stored in myrank.
*
*                       Knowledge about how the processes are 
*                       connected can be found by querying this 
*                       communicator group, e.g.
*
*                       CALL MPI_CART_SHIFT(comm, three_minus_one=2, 1,         *                               neighbor_down, neighbor_up, ierror)
*
*                       Gives the ranks of the neighbors up and down
*                       in the third dimension.
*                       
*                       Most of the useful info is already stored
*                       in neigh(-3:3) and mycoords(1:3). See below
*    
*      
*      myrank           The rank of the processor in the mpi-cartesian
*                       communicator group comm.
*
*      mycoords(1:3)    Your coordinates in the grid-computer grid, zero based.
*                       See example
*
*      neigh(-3:3)      Gives the "myrank" of  neighboring process.
*                       neigh(-3) gives the rank of the nearest z-neighbor
*                       below you, assuming periodic boundary conditions.
*                       neigh(0) gives the root 
*      
* Example:
*
* layout a lattice 16x16x16=n_lattice(1:3) and a computer layout 
* of four processors (1,2,4) = ndimproc(1:3). Then the coordinates
* of the different processors, which is recorded by mycoords(1:3) run over
*
* (0,0,0) = mycoords for myrank = 0
* (0,0,1) = mycoords for myrank = 1
* (0,0,2) = mycoords for myrank = 2
* (0,0,3) = mycoords for myrank = 3
* (0,1,0) = mycoords for myrank = 4
* (0,1,1) = mycoords for myrank = 5
* (0,1,2) = mycoords for myrank = 6
* (0,1,3) = mycoords for myrank = 7
*
* For these eight processes we have the following layout
*-------------------------------
* 0 crds 0      0      0
* 0 negh 3      4      0      0      0      4      1
* 0 bndr 1      1      1      1      8      4
*-------------------------------
* 1 crds 0      0      1
* 1 negh 0      5      1      0      1      5      2
* 1 bndr 1      1      5      1       8      8
*-------------------------------
* 2 crds 0      0      2
* 2 negh 1      6      2      0      2      6      3
* 2 bndr 1      1      9      1       8     12
*-------------------------------
* 3 crds 0      0      3
* 3 negh 2      7      3      0      3      7      0
* 3 bndr 1      1      1      1       8     16
*-------------------------------
* 4 crds 0      1      0
* 4 negh 7      0      4      0      4      0      5
* 4 bndr 1      9      1      1      1       4
*-------------------------------
* 5 crds 0      1      1
* 5 negh 4      1      5      0      5      1      6
* 5 bndr 1      9      5      1      1       8
*-------------------------------
* 6 crds 0      1      2
* 6 negh 5      2      6      0      6      2      7
* 6 bndr 1      9      9      1      1      12
*-------------------------------
* 7 crds 0      1      3
* 7 negh 6      3      7      0      7      3      4
* 7 bndr 1      9      1      1      1      16
* 

***********************************************************************
*
*                     THE COMMON BLOCK
*
***********************************************************************

      common /mpivar   / ierror,ierror1,ierror2,
     &                   comm,xy_slice_comm,z_column_comm,nproc
      common /mpiownvar/ myrank,my_xy_slice_rank,my_z_column_rank,
     &                   mycoords,my_xy_slice_coords,my_z_column_coords,
     &                   neigh,lowbound,upbound,
     &                   ownvol,ownhalfvol
      common /lattice  / n_lattice,dimproc,extension







