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


c  u(nc, maxarray,even:odd,nmu)   Array of links 
c
c      Each link is an SU2 matrix which is represented by two complex numbers.
c      There are even and odd links.
c
c      u(1,...) = a
c      u(2,...) = b   
c
c      SU(2) matrix = {{a, b},{-b*, a}} = {{m_11,m12},{m21,m22}}
c
c      So in terms of paul matrices c0 = Real(a), c1 = Imag(b), c2=Real(b), c3=Imag(a)
c
c      Then U = c0 I + i c_i \sigma_i.
c
c  a(nmu,maxarray,even:odd)   Array containing the gauge fields
c
c      A = a(1) sigma(1) + a(2) * sigma(2) + a(3) * sigma(3)
c
c      Here sigma(2) = {{0, -I}, {I, 0} } //MatrixForm
c
c u_f(nc,maxarray,even:odd,nmu) ! gauge fixed version of u
c
c a_f(nmu,maxarray,even:odd)    ! guage fixed version of a
c











