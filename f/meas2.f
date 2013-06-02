**********************************************************************
*        subroutine for gathering the measured data on root
*        processor and saving to disk           
*  input : meas
**********************************************************************

      subroutine meas2(meas)

      implicit none
      
* include message passing variables

      INCLUDE 'mpif.h'
      INCLUDE 'parallel_parameter.inc'
      INCLUDE 'parallel_mpi_types.inc'

* include global variables and common blocks

      INCLUDE 'parameter.inc'
      INCLUDE 'input_parameter.inc'
      INCLUDE 'measurement.inc'
      INCLUDE 'filenames.inc'

* local variables     

      integer meas,nmeas,mstart

* measure correlation function
  
      mstart=int(confstart/imeas)

      call correlator(meas,mstart)

* calculate sum of gauge fixing iteration for all processors

      if(nproc.gt.1)then
       CALL MPI_REDUCE(iter_or,iter_or_g,1,MPI_INTEGER,MPI_SUM,root,
     &                comm,ierror)
      else
       iter_or_g=iter_or
      endif
  
*  and compute the average quantities

       if (myrank.eq.root) then
          nmeas=int(nsweep/imeas)
          if(meas.eq.nmeas)then
           iter_or_g=iter_or_g/real(nproc)/nmeas
           da=da/nmeas

           open (joblog_unit, file = joblog_file, status = 'unknown',
     &                                  position ='append')
           write(joblog_unit,'(a,1x,i20)')'number of iteration for or. gauge fixing :', 
     &                                     iter_or_g
           write(joblog_unit,'(a,1x,e20.10)')'precision of the condition |d_mu A_mu|^2=0 :',
     &                                        da 
           close(joblog_unit)

          endif

        
       endif

       return
       
       end











