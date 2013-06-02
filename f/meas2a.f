**********************************************************************
**        subroutine for performing measurements with gauge fixed   **
**        quantities, gathering the measured data on root           ** 
**        processor and saving to disk                              **
**  input : meas                                                    ** 
**  P.P  1998.02.26                                                 **
**********************************************************************

      subroutine meas2(meas)

      implicit none
      
** include message passing variables

      INCLUDE 'mpif.h'
      INCLUDE 'parallel_parameter.inc'
      INCLUDE 'parallel_mpi_types.inc'

** include global variables and common blocks

      INCLUDE 'parameter.inc'
      INCLUDE 'input_parameter.inc'
      INCLUDE 'measurement.inc'
      INCLUDE 'gaugefield.inc'
      INCLUDE 'filenames.inc'

** local variables     

      integer meas,nmeas,mstart,eo,is,i

** measure correlation function
  
      mstart=int(confstart/imeas)

      call correlator(meas,mstart)

** measure average of a_f fields on each processor

      a0f=0.0

      do eo=0,1
       do is=1,ownhalfvol(3)
        a0f(1)=a0f(1)+a_f(1,is,eo)
        a0f(2)=a0f(2)+a_f(2,is,eo)
        a0f(3)=a0f(3)+a_f(3,is,eo)
       enddo
      enddo
      
      a0f=a0f/real(ownvol(3))

** calculate sum of gauge fixing iteration for all processors

      if(nproc.gt.1)then
       CALL MPI_REDUCE(iter_or,iter_or_g,1,MPI_INTEGER,MPI_SUM,root,
     &                comm,ierror)
      else
       iter_or_g=iter_or
      endif

** caluculate the global average

      if(nproc.gt.1)then
       do i=1,3
         CALL MPI_REDUCE(a0f(i),af(i),1,MPI_REAL,MPI_SUM,root,
     &                   comm,ierror)
       enddo
      else
        do i=1,3
         af(i)=a0f(i)
        enddo
      endif  

      af=af/real(nproc)
  
**  and compute the average quantities

       if (myrank.eq.root) then
          nmeas=int(nsweep/imeas)
          if(meas.eq.nmeas)then
           iter_or_g=iter_or_g/real(nproc)/(nmeas)
           da=da/nmeas

           open (joblog_unit, file = joblog_file, status = 'unknown',
     &                                  position ='append')
           write(joblog_unit,'(a,1x,i20)')'number of iteration for or. gauge fixing :', 
     &                                     iter_or_g
           write(joblog_unit,'(a,1x,e20.10)')'precision of the condition |d_mu A_mu|^2=0 :',
     &                                        da 
           close(joblog_unit)
          endif
       

        open (afield_unit, file = afield_file, status='unknown',
     &                                         position='append')
            write(afield_unit,'(i5.5,1x,3(e23.16,1x))') meas+mstart, af(1), af(2), af(3)
        
        close(afield_unit)
        
       endif

       return
       
       end











