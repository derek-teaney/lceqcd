**********************************************************************
*        subroutine for gathering the measured data on root          *
*        processor and saving to disk                                *
*  input : isweep                                                    *
**********************************************************************
      
      subroutine meas1(isweep)

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
      
      integer isweep

*   gather the acceptance rates, plaq,akin, a2, a4 on root

      if(isweep.eq.nsweep)then
       if (nproc.gt.1) then
        CALL MPI_REDUCE(acco,accepto,1,MPI_REAL,MPI_SUM,root,
     &                comm,ierror)
        CALL MPI_REDUCE(acch,accepth,1,MPI_REAL,MPI_SUM,root,
     &                comm,ierror)
        CALL MPI_REDUCE(accg,acceptg,1,MPI_REAL,MPI_SUM,root,
     &                comm,ierror)
       else
        accepto=acco
        accepth=acch
        acceptg=accg
       endif
      endif

      if (nproc.gt.1) then
       CALL MPI_REDUCE(plaq0,plaq,1,MPI_REAL,MPI_SUM,root,
     &                comm,ierror)
       CALL MPI_REDUCE(a0kin,akin,1,MPI_REAL,MPI_SUM,root,
     &                comm,ierror)
       CALL MPI_REDUCE(a02,a2,1,MPI_REAL,MPI_SUM,root,
     &                comm,ierror)
       CALL MPI_REDUCE(a04,a4,1,MPI_REAL,MPI_SUM,root,
     &                comm,ierror)
      else
       plaq=plaq0
       akin=a0kin
       a2=a02
       a4=a04
      endif

*  and compute the average quantities on root and save them on disk

       if (myrank.eq.root) then

*  acceptance rates

          if(isweep.eq.nsweep)then
           accepto=accepto/real(nproc)/nsweep
           accepth=accepth/real(nproc)/nsweep
           acceptg=acceptg/real(nproc)/nsweep
           open (joblog_unit, file = joblog_file, status = 'unknown',
     &                                  position ='append')
           write(joblog_unit,'(a,1x,e10.5)')'Overrelaxation acceptance :',
     &                                      accepto
           write(joblog_unit,'(a,1x,e10.5)')'heat-bath acceptance : ', accepth
           write(joblog_unit,'(a,1x,e10.5)')'Gauge  acceptance    : ', acceptg
           close(joblog_unit)
          endif

*  averages of local observables

          plaq=plaq/real(nproc)
          akin=akin/real(nproc)
          a2=a2/real(nproc)
          a4=a4/real(nproc)
        
        if(mod(isweep,100).eq.1)then
         open (action_unit, file=action_file, status='unknown',
     &                              position='append')
        endif

        write(action_unit,'(i5,1x,4(e23.16,1x))')isweep+confstart,plaq,akin,a2,a4
        
        if(mod(isweep,100).eq.0)then
        close(action_unit)
        endif
         
         

       endif

 
       return
       end


























