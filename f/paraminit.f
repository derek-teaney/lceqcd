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
*     initialisation of the necessary parameters                   *
*                                                                  *
********************************************************************


      subroutine paraminit

      USE mpi
      implicit none
*
*including all global MESSAGE PASSING variables 
*
      INCLUDE 'parallel_parameter.inc'
      INCLUDE 'parallel_mpi_types.inc'
*
*including global variables, parameters and common blocks
*
      INCLUDE 'parameter.inc'
      INCLUDE 'input_parameter.inc'
      INCLUDE 'const.inc'
      INCLUDE 'filenames.inc'
*
*local VARIABLES
*
      integer i
      real    sigma
***********************************************************************
*
* initialize the constants
*
      zi = (0.0,1.0)
      pi = 2.0 * acos(0.0)
      twopi = 2.0 * pi
      sigma=3.1759114

*
*     only for ROOT
*     read the input from inputfile
*
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierror)
      if (myrank.eq.root) then

******************************************************************************
*
*      nx/ny/nz_lattice    ! is the global lattice size
*      dimproc             ! numbers of processes in every dimension mu = 1-4  
*      nproc               ! total numbers of processes                        
*      init                ! controls initialisation of gauge fields           
*      rand                ! old(=1) or new(=0) random numbers
*      save                ! save on disc yes or no (1,0)
*      confstart           ! starting configuration                            
*      beta                ! gauge coupling                
*      x                   ! continuum parameter1
*      y                   ! continuum parameter2
*      lambda              ! overrelaxation parameter
*      limit_or            ! precision of (D_mu A_mu)^2 = 0
*      inorm               ! number of sweeps after which to reunitarize links
*      nsweep_flag         ! (0=number of sweeps 1=time in min)
*      ntime               ! maximum of time for the program
*      nsweep              ! number of trajectories
*      imeas               ! number of sweeps after which to measure things
*      imode               ! heatbath mode    
*******************************************************************************
         input_unit=1
         inputfile='inputfile'
         
         open(Unit=input_unit,file=inputfile,status='unknown')

         read(input_unit,*,err=20) n_lattice(1), n_lattice(2), n_lattice(3)
         read(input_unit,*,err=20) dimproc(1),   dimproc(2),   dimproc(3)   
         read(input_unit,*,err=20) nproc
         read(input_unit,*,err=20) init
         read(input_unit,*,err=20) rand
         read(input_unit,*,err=20) save
         read(input_unit,*,err=20) confstart
         read(input_unit,*,err=20) beta    
         read(input_unit,*,err=20) x    
         read(input_unit,*,err=20) y    
         read(input_unit,*,err=20) lambda
         read(input_unit,*,err=20) alpha
         read(input_unit,*,err=20) limit_or
         read(input_unit,*,err=20) limit_it
         read(input_unit,*,err=20) inorm
         read(input_unit,*,err=20) nsweep_flag   
         read(input_unit,*,err=20) ntime   
         read(input_unit,*,err=20) nsweep   
         read(input_unit,*,err=20) imeas
         read(input_unit,*,err=20) imode   

         close(input_unit)
         
*     
*     error message 
*     
         go to 30
 20      write(*,*)'break: data input',n_lattice,dimproc,nproc,init,rand,
     &        save,confstart,beta,x,y,lambda,limit_or,inorm,nsweep_flag,ntime,
     &        nsweep,imeas,imode
         call abort()
 30      continue
         
         
*     
*     the parameters are given for the user
*     
         write(6,1000)n_lattice(1), n_lattice(2), n_lattice(3),
     &        dimproc(1), dimproc(2), dimproc(3), nproc
         write(*,*)'init 0=rand,1=onerand,2=unit,3=old   :', init        
         write(*,*)'rand 1=old, 0=new random numbers     :', rand        
         write(*,*)'save on disc yes or no (1,0)         :', save        
         write(*,*)'start configuration                  :', confstart   
         write(*,*)'beta                                 :', beta        
         write(*,*)'x                                    :', x
         write(*,*)'y                                    :', y
         write(*,*)'overrelaxation parameter             :', lambda
         write(*,*)'gauge transformation parameter       :', alpha
         write(*,*)'limit_or                             :', limit_or
         write(*,*)'limit_it parameters                  :', limit_it
         write(*,*)'inorm                                :', inorm
         write(*,*)'flag (0=nb of sweeps, 1=time in min) :', nsweep_flag  
         write(*,*)'max time for the program in min      :', ntime   
         write(*,*)'number of sweeps                     :', nsweep       
         write(*,*)'imeas                                :', imeas        
         write(*,*)'heatbath mode                        :', imode        
      endif
*
*     create new types for the broadcast
*     
      call paratypeconstruct
      
      CALL MPI_BCAST(beta,1,paramettype,root,MPI_COMM_WORLD,ierror)
      CALL MPI_BCAST(init,1,configtype,root,MPI_COMM_WORLD,ierror)
      CALL MPI_BCAST(n_lattice(1),1,latconftype,root,MPI_COMM_WORLD,
     &     ierror)
      CALL MPI_BCAST(nproc,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierror)

      beta_a=beta
      beta_2=beta_a*( 3+8*y/beta**2-sigma*(4+5*x)/(2*pi*beta)-
     &       ( (20*x-10*x**2)*(log(1.5*beta)+0.09)+8.7+11.6*x)/(2*pi**2*beta**2) )
      beta_4=x*beta
          
      do i = 0, nproc-1
         if(myrank==i)then
            write(*,*)'********************************************'
            write(*,*)'myrank = ', myrank
            write(*,*)'********************************************'
            write(*,1000) n_lattice(1), n_lattice(2), n_lattice(3),
     &                  dimproc(1),   dimproc(2),   dimproc(3)   
            write(*,*)'init 0=rand,1=onerand,2=unit,3=old   :', init        
            write(*,*)'rand 1=old, 0=new random numbers     :', rand        
            write(*,*)'save on disc yes or no (1,0)         :', save        
            write(*,*)'start configuration                  :', confstart   

            write(*,*)'beta                                 :', beta        
            write(*,*)'x                                    :', x
            write(*,*)'y                                    :', y
            write(*,*)'overrelaxation parameter             :', lambda
            write(*,*)'gauge transformation parameter       :', alpha
            write(*,*)'limit_or                             :', limit_or
            write(*,*)'limit_it parameters                  :', limit_it
            write(*,*)'inorm                                :', inorm
            write(*,*)'flag (0=nb of sweeps, 1=time in min) :', nsweep_flag  
            write(*,*)'max time for the program in min      :', ntime   
            write(*,*)'number of sweeps                     :', nsweep       
            write(*,*)'imeas                                :', imeas        
            write(*,*)'heatbath mode                        :', imode        

         end if
         call MPI_BARRIER(comm, ierror)
      end do
      
      return

1000  format('Adjoint Higgs Modell on a '/
     &        6x,i2,' X 'i2,' X ',i2,' lattice'/
     &     'on a  ',i2,' X ',i2,' X ',i2,' processor lattice'/
     &        'number of processors :',i3/)

      end
      




























