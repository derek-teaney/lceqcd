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
*     build a virtual TOPOLOGY for the 4-D lattice                 *
*     store the information in the common block mpivar             *
*     included in parallel.inc                                     *
*     make sure that the minimum number of points in the local     *
*     lattice is equal or greater than two and does not            *
*     contain more than 8 points in each direction                 *
*                                                                  *
********************************************************************

      subroutine topology

      IMPLICIT none
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
*
*local VARIABLES
*
      INTEGER mu,nproc2,i
      integer dim,dimproc1(npdim-1),dimproc2(npdim-2)
      LOGICAL periods(npdim),reorder,remain_dim_xy_slice(npdim),
     1        remain_dim_z_column(npdim)

**********************************************************************
*
*first find out how many PROCESSES are used
*
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc2,ierror)

      if ((nproc.ne.nproc2).AND.(myrank.eq.root)) then
         write(*,*)'break : number of processes does not match'
         call abort()
      endif
*
* set the TOPOLOGY of the PROCESSES with PERIODIC boundaries
* they are given by hand with the input data
* but it has to be checked out wether they are fitting or not
* the biggest local lattice is allowed to contain not morer than
* nx*ny*nz lattice sites
*
      extension(1)=n_lattice(1)/dimproc(1)
      extension(2)=n_lattice(2)/dimproc(2)
      extension(3)=n_lattice(3)/dimproc(3)
     
      ownvol(1)=extension(1)
      ownhalfvol(1)=extension(1)/2
      ownvol(2)=extension(2)*ownvol(1)
      ownhalfvol(2)=ownvol(2)/2
      ownvol(3)=extension(3)*ownvol(2)
      ownhalfvol(3)=ownvol(3)/2
*
*break if something doesn't fit
*      
      if (extension(1).gt.nx) then
         if (myrank.eq.root) write(*,*)'break: x extension > nx'
         call abort()
      endif
      if (extension(2).gt.ny) then
         if (myrank.eq.root) write(*,*)'break: y extension > ny'
         call abort()
      endif
      if (extension(3).gt.nz) then
         if (myrank.eq.root) write(*,*)'break: z extension > nz'
         call abort()
      endif
*
* it can happen, that ther is not enough allocated storage to do the
* save of the configurations:
*
      if ( dimproc(1)*dimproc(2) > extension(3)) then
         if (myrank==root) then
            write(*,*)'Not enough space to read/write the configuration!!!!'
            write(*,*)'dimproc(1)*dimproc(2) must be <= extension(3)'
         end if
         call abort()
      end if
      do mu=1,npdim

         if (extension(mu)*dimproc(mu).NE.n_lattice(mu)) then
            if (myrank.eq.root) write(*,*)'break: n_lattice/processes
     1                                     not an integer'
            call abort()
         endif

         if (extension(mu).lt.2) then
            If (myrank.eq.root) write(*,*)'break : extension < 2'
            call abort()
         endif

      end do


*
*define the periodicity of the lattice
*
      do  mu=1,npdim
         periods(mu)=.TRUE.
      end do
      reorder=.TRUE.
*
* then define the lattice COMMUNICATOR with a topology given by the
* input data
*
      CALL MPI_CART_CREATE(MPI_COMM_WORLD,npdim,dimproc,periods,
     1                     reorder,comm,ierror)

*
* tell the process who he is and who the neighbours are
* neigh(mu) is the rank of the process
* the numbers of the dimensions and the processes in each dimension
* begin with 0
*
      CALL MPI_COMM_RANK(comm,myrank,ierror)
      do 30 mu=1,npdim
         CALL MPI_CART_SHIFT(comm,mu-1,1,neigh(-mu),neigh(mu),
     1                       ierror)
 30   continue
      neigh(0)=root

*
*compute the coordinates of the process begining  with (0,0,0,0)
*     
      CALL MPI_CART_COORDS(comm,myrank,npdim,mycoords,ierror)
ccc      write(*,*)myrank,'mycoords',(mycoords(i),i=1,4)
ccc      write(*,*)myrank,'neigh',(neigh(i),i=-4,4)

*
*compute the bounddaries of the local lattice
*
      do  mu=1,npdim
         lowbound(mu)=mycoords(mu)*extension(mu)+1
         upbound(mu)=lowbound(mu)+extension(mu)-1
      end do
ccc      write(*,*)myrank,'boundaries',(lowbound(i),i=1,4),
ccc     1                              (upbound(i),i=1,4)

*
* remove z-coordinate from the topology and build xy-slices to compute
* xy_slice averages of gauge and higgs fields
*
      remain_dim_xy_slice(1)=.true.
      remain_dim_xy_slice(2)=.true.
      remain_dim_xy_slice(3)=.false.

      CALL MPI_CART_SUB(comm,remain_dim_xy_slice,xy_slice_comm,ierror)

      CALL MPI_CARTDIM_GET(xy_slice_comm,dim,ierror)

      CALL MPI_CART_GET(xy_slice_comm,dim,dimproc1,periods,my_xy_slice_coords,
     &                                                                ierror)
      CALL MPI_CART_RANK(xy_slice_comm,my_xy_slice_coords,my_xy_slice_rank,
     &                                                                ierror)
*
* build a topology to gather xy_slice-averaged gauge and higgs fields
*
      remain_dim_z_column(1)=.false.
      remain_dim_z_column(2)=.false.
      remain_dim_z_column(3)=.true.

      CALL MPI_CART_SUB(comm,remain_dim_z_column,z_column_comm,ierror)

      CALL MPI_CARTDIM_GET(z_column_comm,dim,ierror)

      CALL MPI_CART_GET(z_column_comm,dim,dimproc2,periods,my_z_column_coords,
     &                                                                ierror)
      CALL MPI_CART_RANK(z_column_comm,my_z_column_coords,my_z_column_rank,
     &                                                                ierror)

*      do i = 0, nproc-1
*         if(myrank==i)then
*            write(*,*)             'myrank             =',myrank
*            write(*,'(A,3(1x,I2))')'mycoords           =',mycoords
*            write(*,*)'------------------------------------------------------'
*            write(*,'(A,1(1x,I2))')'my_xy_slice_rank   =',my_xy_slice_rank
*            write(*,'(A,2(1x,I2))')'my_xy_slice_coords =',my_xy_slice_coords
*            write(*,*)'------------------------------------------------------'
*            write(*,*)             'my_z_column_rank   =',my_z_column_rank
*            write(*,'(A,1(1x,I2))')'my_z_column_coords =',my_z_column_coords
*            write(*,*)'******************************************************'
*         end if
*         call MPI_BARRIER(comm)
*      end do
*
      return

      end






