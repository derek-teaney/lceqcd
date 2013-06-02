****************************************************************
*                                                              *
*     read the data from the disc using a array with the       *
*     length for one global xy-plane                           *
*     savepara defines how many times vectors of length        *
*     ownhalfvol(1) are read in                                *
*                                                              *
**************************************************************** 


      subroutine readdata(helparray,ownunit,rl,readpara)


      implicit none


*
*including global variables, parameters and common blocks
*
      INCLUDE 'parameter.inc'


*
*local VARIABLES
*
      integer ownunit,rl,readpara
      complex helparray(readpara)


***********************************************************************


      read(Unit=ownunit,rec=rl,err=150)helparray

      goto 200
 150  call abort()
 200  continue


      return

      end







