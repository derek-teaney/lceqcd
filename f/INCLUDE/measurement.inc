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
* INCLUDE FILE                                                     *
*                                                                  *
* contains measured quantities                                     *
*                                                                  *
********************************************************************
********************************************************************
      real plaq0,a0kin,a02,a04,a0f(3),plaq,akin,a2,a4,af(3)
      real acco,acch,accg,accepto,accepth,acceptg
      real da
      integer iter_or,iter_or_g
********************************************************************
*                                                                  *
*                                                                  *
*                     THE COMMON BLOCK                             *
*                                                                  *
********************************************************************
      common/aver_loc/plaq0,a0kin,a02,a04,a0f
      common/aver_glob/plaq,akin,a2,a4,af
      common/acc_loc/acco,acch,accg
      common/acc_glob/accepto,accepth,acceptg
      common/or_loc/iter_or
      common/or_glob/da,iter_or_g
********************************************************************
********************************************************************

