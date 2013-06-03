********************************************************************
********************************************************************
*                                                                  *
* Written by Peter Schmidt for an SU(3) pure gauge code and a      *
* Wilson code to invert the fermion matrix                         *
*                                                                  *
* Adapted for an SU(2) adjoint higgs model by Manfred Oevers       *
* November 1997                                                    *
*                                                                  *
********************************************************************
********************************************************************
********************************************************************
*                                                                  *
*         initialisation of the measured quantities                *
*                                                                  *
********************************************************************

      subroutine measinit

      implicit none
*
*include message passing features
*
      include 'mpif.h'
      include 'parallel_parameter.inc'
      include 'parallel_mpi_types.inc'
*
*including global variables, parameters and common blocks
*
      include 'parameter.inc'
      include 'measurement.inc'
      include 'input_parameter.inc'
      include 'pointer.inc'
      include 'filenames.inc'
*
*local variables
*
      integer i, betaint,xint,yint, nproc2, skip(1), n,
     &        mu, site
      real rvec(nseed)

*************************************************************************

*
* initialisation of pointer to own site
*
      do  i = 1, nsiteshalf
         normalsites (i) = i
      end do

      betaint    = nint(1000 * beta)
      xint  = nint(10000 * x)
      yint  = nint(10000 * y)
*
* open the files for measurement
*

* action_file

      if (myrank.eq.root) then

         write (action_file, fmt = 1000) n_lattice(3),
     &                betaint, xint, yint, confstart

         action_unit = 3

         open (action_unit, file=action_file, status='unknown',
     &                              position='append')

         write(action_unit,1001) n_lattice
         write(action_unit,1002) beta,x,y
         write(action_unit,1003) nsweep

         close(action_unit)

* joblog_file

         write (joblog_file, fmt = 1100) n_lattice(3),
     &                betaint, xint, yint, confstart

         joblog_unit = 4

         open (joblog_unit, file = joblog_file, status = 'unknown',
     &                                  position ='append')

         write(joblog_unit,1001) n_lattice
         write(joblog_unit,1002) beta, x, y
         write(joblog_unit,1003) nsweep

* u_conf_file

         write (u_conf_file, fmt = 1200) n_lattice(3),
     &                betaint, xint, yint, confstart

         u_unit = 20

* a_conf_file

         write (a_conf_file, fmt = 1300) n_lattice(3),
     &                betaint,xint,yint,confstart

         a_unit = 21

* corr_file

         write (corr_file, fmt = 1500) n_lattice(3),
     &                betaint,xint,yint,confstart

         corr_unit = 31
         open(corr_unit, file=corr_file, status='unknown',
     &                                   position='append' )

           write(corr_unit,1001) n_lattice
           write(corr_unit,1002) beta, x, y
           write(corr_unit,1003) nsweep
         close(corr_unit)

* afield_file

         write (afield_file, fmt = 1600) n_lattice(3),
     &                betaint,xint,yint,confstart

         afield_unit = 32
         open(afield_unit, file=afield_file, status='unknown',
     &                                   position='append' )

           write(afield_unit,1001) n_lattice
           write(afield_unit,1002) beta, x, y
           write(afield_unit,1003) nsweep
         close(afield_unit)

* ckeckp_file

         write (checkp_file, fmt = 1700) n_lattice(3),
     &                betaint,xint,yint,confstart

         checkp_unit = 33
         open(checkp_unit, file=checkp_file, status='unknown',
     &                                   position='append' )

           write(checkp_unit,1001) n_lattice
           write(checkp_unit,1002) beta, x, y
           write(checkp_unit,1003) nsweep
         close(checkp_unit)

* polyak_file

         write (polyak_file, fmt = 1800) n_lattice(3),
     &                betaint,xint,yint,confstart

         polyak_unit = 34
         open(polyak_unit, file=polyak_file, status='unknown',
     &                                   position='append' )

           write(polyak_unit,1001) n_lattice
           write(polyak_unit,1002) beta, x, y
           write(polyak_unit,1003) nsweep
         close(polyak_unit)


* corra2_file

         write (corra2_file, fmt = 1900) n_lattice(3),
     &                betaint,xint,yint,confstart

         corra2_unit = 35
         open(corra2_unit, file=corra2_file, status='unknown',
     &                                   position='append' )

           write(corra2_unit,1001) n_lattice
           write(corra2_unit,1002) beta, x, y
           write(corra2_unit,1003) nsweep
         close(corra2_unit)
      end if

*
*read the random numbers from the disc if rand = 1
*
      if (rand.eq.1) then
        if (myrank.eq.root) then
          write (random_number_file, fmt = 1400) n_lattice(3),
     &                betaint, xint, yint,  confstart

          random_number_unit = 2

          open (random_number_unit,file=random_number_file,status ='old')
          read (random_number_unit, *, err = 150) 
     &                       (iseedblock (i), i = 0,nseedblock)
          close (random_number_unit)
        endif

*
*then scatter them to the other processes
*
        call mpi_scatter (iseedblock(0), nseed, mpi_integer,
     &      iseed(1),nseed, mpi_integer, root, comm, ierror)

*
*take new random numbers
*
      else
        iseed = 0
      endif

*
*first find out how many processes are used
*
      call mpi_comm_size (comm, nproc2, ierror)
      call mpi_comm_rank (comm, myrank, ierror)

*
*then initialise the random number generator on each processor
*
      call rluxgo(3,314159+myrank,0,0)
      call ranlux(rvec,nseed) 
      iseed=int(10**6*rvec)
      call random_seed(put=iseed)
      rand=1
*
*error message
*
      goto 155

  150 if (myrank.eq.root) print *,
     &     'break :while reading the random numbers'
           stop

  155 continue

*
*format
*

 1001 format('# ',3(i2.2,1x))

 1002 format('# ',3(f7.4,1x))

 1003 format('# ',i5.5)

 1000 format('action_s', i2.2,  '_be', i5.5, '_x', i5.5, '_y', i5.5,
     &                                        '_u',  i6.6)

 1100 format('joblog_s', i2.2,  '_be', i5.5, '_x', i5.5, '_y', i5.5,
     &                                        '_u',  i6.6)

 1200 format('u_conf_s', i2.2,  '_be', i5.5, '_x', i5.5, '_y', i5.5,
     &                                        '_u',  i6.6)

 1300 format('a_conf_s', i2.2,  '_be', i5.5, '_x', i5.5, '_y', i5.5,
     &                                        '_u',  i6.6)

 1400 format('random_s', i2.2,  '_be', i5.5, '_x', i5.5, '_y', i5.5,
     &                                        '_u',  i6.6)

 1500 format('correl_s', i2.2,  '_be', i5.5, '_x', i5.5, '_y', i5.5,
     &                                        '_u',  i6.6)

 1600 format('afield_s', i2.2,  '_be', i5.5, '_x', i5.5, '_y', i5.5,
     &                                        '_u',  i6.6)

 1700 format('checkp_s', i2.2,  '_be', i5.5, '_x', i5.5, '_y', i5.5,
     &                                        '_u',  i6.6)

 1800 format('polyak_s', i2.2,  '_be', i5.5, '_x', i5.5, '_y', i5.5,
     &                                        '_u',  i6.6)

 1900 format('corra2_s', i2.2,  '_be', i5.5, '_x', i5.5, '_y', i5.5,
     &                                        '_u',  i6.6)

      return
      end









