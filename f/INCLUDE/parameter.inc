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
*                                                                  *
*     INCLUDE FILE                                                 *
*                                                                  *
*     the PARAMETERS                                               *
*     for the adjoint higgs code                                   *
*                                                                  *
*                                                                  *
********************************************************************
********************************************************************
*
*                   DEFINITION OF THE PARAMETERS
*these parameters are used to define the maximum size of allocated
*memory for a the maximum amount of grid points on one local lattice
***********************************************************************

      integer nc,nx,nxhalf,ny,nz,nsites,nsiteshalf,nlink,
     1        Volxy,Volxyz,Volxyhalf,maxarray,
     2        nmu,even,odd

      PARAMETER (nc   = 2)  ! number of colours
      PARAMETER (nmu  = 3)  ! dimension, and also used as nc^2-1 in spots
      PARAMETER (even = 0)  ! loops over even, and odd
      PARAMETER (odd  = 1)  !

      PARAMETER (nx = 24)    ! maximal number of sites in
      PARAMETER (ny = 24)    ! x,y,z direction this should be an 
      PARAMETER (nz = 24)    ! even number

      PARAMETER (nxhalf     = nx/2)    ! number of even slots in x
      PARAMETER (Volxy      = nx*ny)   ! number of xy slots
      PARAMETER (Volxyhalf  = Volxy/2) ! number of even xy slots
      PARAMETER (nsites     = nx*ny*nz) ! Volume
      PARAMETER (nsiteshalf = nsites/2) ! number of even sits
      PARAMETER (nlink      = nsites*nmu)    ! number of links

* maximal spatial size of the linkfield with even/odd decomposition
      PARAMETER (maxarray   = nsiteshalf+Volxyhalf+4)

C 
