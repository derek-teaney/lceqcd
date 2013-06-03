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
********************************************************************
*                                                                  *
*     INCLUDE FILE                                                 *
*                                                                  *
*     the input PARAMETER for the adjoint higgs model code         *
*                                                                  *
*                                                                  *
********************************************************************
********************************************************************
      ! gfortran needs twelve seeds to start rd nb generator
      INTEGER, parameter :: nseed=16       
      ! Dimensions of array storing 
      INTEGER, parameter :: nseedblock=4095 

c !!! WARNING !!!!
c
c Do not reoorder common blocks without great care to modify 
c paratypeconstruct.f
c

      REAL     beta               ! gauge coupling          
      REAL     x                  ! continuum parameter 1
      REAL     y                  ! continuum parameter 2
      REAL     lambda             ! overrelaxation parameter
      REAL     alpha              ! parameter for iter. gauge fixing
      REAL     limit_or           ! limit on (D_mu A_mu)^2 = 0
      REAL     limit_it           ! limit on (D_mu A_mu)^2 = 0

      INTEGER  inorm              ! number of steps before normalizing u_f
      INTEGER  nsweep_flag        ! flag to toggle between nsweep and ntime
      INTEGER  ntime              ! maximum of time for the program
      INTEGER  nsweep             ! number of nsweep                
      INTEGER  imeas              ! number of steps without measurements
      INTEGER  imode              ! heatbath mode

      INTEGER  iseed(nseed)       ! local seed for the rd nb generator
      INTEGER  iseedblock(0:nseedblock) ! all seeds
      REAL     beta_a,beta_2,beta_4 ! aux. parameters

      INTEGER  init               ! controls gauge initialisation
      INTEGER  rand               ! old(=0) or new(=1) random numbers
      INTEGER  save               ! save on disc yes or no (1,0)
      INTEGER  confstart          ! starting configuration


***********************************************************************
*
*                     THE COMMON BLOCK
*
***********************************************************************


      common /paramet/ beta,x,y,lambda,alpha,limit_or,limit_it,
     &                 inorm,nsweep_flag,ntime,
     &                 nsweep,imeas,imode,iseed,iseedblock,
     &                 beta_a,beta_2,beta_4
      common /config / init,rand,save,confstart


***********************************************************************
***********************************************************************
***********************************************************************
