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







