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
* finishes up by saving configuration and measured quantities      *
*                                                                  *
********************************************************************

      subroutine finalize

      implicit none
*
*including all global MESSAGE PASSING features and variables 
*
      INCLUDE 'mpif.h'
      INCLUDE 'parallel_parameter.inc'
      INCLUDE 'parallel_mpi_types.inc'
*
* include some global variables
*
      INCLUDE 'parameter.inc'
      INCLUDE 'input_parameter.inc'
      INCLUDE 'times.inc'
      INCLUDE 'filenames.inc'

      integer          betaint,xint, yint, i

***************************************************************

      firsttime(2)=MPI_WTIME()
*
* save configuration to disk
*
      call confsave
*
* each process gets its random number seed
* then root gets the seeds from the other processes
*
      call random_seed(GET=iseed)

      CALL MPI_GATHER(iseed,1,MPI_INTEGER,iseedblock(0),1,MPI_INTEGER,
     &                root,comm,ierror)
*
* dump the random numbers on disc
* new random_number_file

      if (myrank.eq.root) then

         betaint    = nint(1000 * beta)
         xint  = nint(10000 * x)
         yint  = nint(10000 * y)

         write (random_number_file, fmt = 1000) n_lattice(3),  
     &        betaint, xint, yint, (confstart+nsweep)

         open(random_number_unit,file=random_number_file,
     &                                status='unknown',err=160)
         write(random_number_unit,*)(iseedblock(i),i=0,255)
         close(random_number_unit)
      endif

      if (myrank.eq.root) then

         write (lastconfig_file, fmt = 1200) n_lattice(3),  
     &        betaint, xint, yint

         open(666,file=lastconfig_file,
     &                                status='unknown',err=160)
         write(666,'(i6.6)')nsweep+confstart
         close(666)
      endif



*
*save the new input with new starting configuration 
*
      if (myrank.eq.root) then

         if (save.eq.1) then
            init=3
            rand=1
         endif


       open(Unit=input_unit,file=inputfile,status='unknown')

       write(input_unit,'(3(I2.2,1x),3x,A)',err=120) n_lattice,
     &        '!nx/ny/nz global lattice size'

       write(input_unit,'(3(I2.2,1x),3x,A)',err=120) dimproc,
     &        '!dimproc     numbers of processes in every dimension mu = 1-3'

       write(input_unit,'(I3.3,9x,A)',err=120) nproc,
     &        '!nproc       total numbers of processes'

       write(input_unit,'(I1.1,11x,A)',err=120) init,
     &        '!init        init for gaugefields 0=rand,1=onerand,2=unit,3=old'

       write(input_unit,'(I1.1,11x,A)',err=120) rand,
     &         '!rand        old(=1) or new(=0) random numbers'

       write(input_unit,'(I1.1,11x,A)',err=120) save,
     &          '!save        save on disc yes or no (1,0)'

       write(input_unit,'(I6.6,6x,A)',err=120) confstart + nsweep,
     &          '!confstart   starting configuration'

       write(input_unit,'(F6.3,6x,A)',err=120) beta,    
     &          '!beta'

       write(input_unit,'(F6.4,6x,A)',err=120) x,    
     &          '!x'

       write(input_unit,'(F6.4,6x,A)',err=120) y,    
     &          '!y'      

       write(input_unit,'(F6.4,6x,A)',err=120) lambda,    
     &          '!lambda overrelaxation parameter'

       write(input_unit,'(F8.6,4x,A)',err=120) alpha,
     &          '!gauge transformation parameter'

       write(input_unit,'(e9.3,3x,A)',err=120) limit_or,    
     &          '!limit_or precision of the condition d_mu A_mu=0'
              
       write(input_unit,'(e9.3,3x,A)',err=120) limit_it,   
     &          '!limit_it precision of the condition d_mu A_mu=0'

       write(input_unit,'(I3.3,9x,A)',err=120) inorm,   
     &          '!inorm   number of iteration  followed by norm. '     

       write(input_unit,'(I1,11x,A)',err=120) nsweep_flag,   
     &          '!nsweep_flag  (0=number of sweeps 1=max time in min) '

       write(input_unit,'(I3.3,9x,A)',err=120) ntime,   
     &          '!ntime       maximal time for the program in minuntes '

       write(input_unit,'(I5.5,7x,A)',err=120) nsweep,   
     &          '!nsweeps       number of sweeps '

      write(input_unit,'(I5.5,7x,A)',err=120) imeas,   
     &          '!imeas  number of sweeps followed by measurement'

      write(input_unit,'(I3.3,9x,A)',err=120) imode,   
     &          '!imode  mode for the heatbath algorithm'





       close(input_unit)


         go to 130
 120     write(*,*)'break: data input'
         call abort()
 130     continue
         
      end if

      lasttime(2)=MPI_WTIME()
      needtime(2)=lasttime(2)-firsttime(2)

      if (myrank.eq.0) then 
         write(*,*)'Time needed for finishing up ', needtime(2)
      endif
*
*error message 
*
      go to 155

 160  if (myrank.eq.root) then
         write(*,*)'break : writing the random numbers'
         call abort()
      endif

 155  continue
*
* timing for the whole propgram
*
      needtime(1) = MPI_WTIME() - firsttime(1)
*
* write timings to stdout
*
      if (myrank == root) then
            write(*,*)'Program             :  ',needtime(1) 
         if (nsweep/=0) then
            write(*,*)'time per trajectory :  ',
     &           (needtime(2) / real(nsweep))
         end if
      endif
*
* formats for the output files
*
 1000 format('random_s', i2.2,  '_BE', i5.5, '_x', i5.5, '_y', i5.5,
     &                                        '_U',  i6.6)
 1100 FORMAT(I5.5)

 1200 format('lastconfig_s', i2.2, '_BE', i5.5, '_x', i5.5, '_y', i5.5)
            
      end subroutine finalize








