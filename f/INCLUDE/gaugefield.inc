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
*     INCLUDE FILE                                                 *
*                                                                  *
*     the gauge field and adjoint Higgs field                      *
*                                                                  *
*                                                                  *
********************************************************************

      COMPLEX    u(nc, maxarray,even:odd,nmu)  ! linkfield
      REAL       a(nmu,maxarray,even:odd)      ! higgsfield
      COMPLEX    u_f(nc,maxarray,even:odd,nmu) ! gauge fixed linkfiled
      REAL       a_f(nmu,maxarray,even:odd)    ! guage fixed higgsfield
      COMPLEX    v(nc,maxarray, even:odd,nmu)  ! staple
      REAL       force(3,maxarray)             ! force
********************************************************************
*
*                     THE COMMON BLOCK
*
***********************************************************************

cdir$ cache_align /linkvar/
cdir$ cache_align /helpfields/
      common /linkvar/ u,a,u_f,a_f
      common /helpfields/v,force 
***********************************************************************
***********************************************************************
***********************************************************************














