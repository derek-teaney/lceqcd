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


      subroutine typeconstruct

c Purpose:
c                                                                  
c     To initialisation of the types for MESSAGE PASSING 
c     filling up the common blocks described in parallel_mpi_types.inc         
c
c Inputs
c
c     nc=2      Number of colors from parameter.inc
c
c     nx,ny     Maximum size of the xy planes. See parameter.inc
c
c     nmu=3     Number of dimensions. See parameter.inc
c 
c     ownhalfvol  See parallel_parameter.inc. Describes the local
c                 lattice
c
c     npdim     See parallel_parameter.inc Used as a proxy for nc^2 -1. CHANGEME
c
c     myrank/root
c
c                See parallel_parameter.inc Used for debugging only. When
c                printing s commented out (as in production runs) these
c                parameters
c
c     dimproc(1:3)  
c
c                See parallel_parameter.inc. Dimesnions of computer grid
c
c     extension  Size of domains 
c
c Outputs:
c 
c     All of the parameters of the common blocks defined in
c     parallel_mpi_types.inc. See this file for a full description
c
c     ierror    Generic MPI error code. See paralllel_parameter.inc

      implicit none


*
*including all global MESSAGE PASSING variables 
*
      INCLUDE 'mpif.h'
      INCLUDE 'parallel_parameter.inc'
      INCLUDE 'parallel_mpi_types.inc'

*
*including global variables, parameters and common blocks
*
      INCLUDE 'parameter.inc'
      INCLUDE 'gaugefield.inc'

*
*local PARAMETERS
*
      integer    count
      parameter (count=nx*ny/2)

*
*definition of the VARIABLES used for MPI
*
      integer a_of_b_TYPE(count,0:1,nmu,0:1),
     1        a_of_d_TYPE(count,0:1,nmu,0:1),
     2        a_of_t_TYPE(count,0:1,nmu,0:1),
     6        a_of_b(8),
     7        a_of_d(8),
     8        a_of_t(8)

      integer eo,ad2,ad1,extent2,extent


*
*local VARIABLES
*
      integer iter,iter0,iter1,jump,x,y,z,t,mu


*
*create new types
*1) to transfer a matrix for diverse subroutines
*

      CALL MPI_TYPE_CONTIGUOUS(nc,MPI_COMPLEX,matrixtype,ierror)
      CALL MPI_TYPE_COMMIT(matrixtype,ierror)
      CALL MPI_TYPE_EXTENT(matrixtype,extent2,ierror)

      CALL MPI_ADDRESS(u(1,1,0,1),ad1,ierror)
      CALL MPI_ADDRESS(u(1,2,0,1),ad2,ierror)

      extent=ad2-ad1

*      if(myrank==root)then
*         write (*,*) 'matrixtype'
*         write (*,*) 'MPI_EXTENT=',extent2
*         write (*,*) 'EXTENT=',extent
*      end if


****************************************************************************
*
*2) to transfer xy-slize of gauge field for read and write from/to disk
*
****************************************************************************

      CALL MPI_TYPE_CONTIGUOUS( ownhalfvol(2),
     &                          matrixtype,load_u_type,ierror)
      CALL MPI_TYPE_COMMIT(load_u_type,ierror)

      CALL MPI_TYPE_CONTIGUOUS( ownhalfvol(2),
     &                          matrixtype,save_u_type,ierror)
      CALL MPI_TYPE_COMMIT(save_u_type,ierror)

****************************************************************************

*
*3) to send a selected part of the links or the gauge field
*   ((usefull),.....,(usefull),.....,(usefull),.....)
*   and to receive it in a contiguous buffer
*   ((usefull),(usefull),(usefull))
*

* type to receive linkfields in a contiguous buffer

      CALL MPI_TYPE_CONTIGUOUS(count,matrixtype,u_CONT_TYPE,ierror)
      CALL MPI_TYPE_COMMIT(u_CONT_TYPE,ierror)


* now define the types for sending the linkfields
*
* consider first the movement down and then up
* for this two different types of transfer patterns are necessary
*
* I think Peter has assumed in what follows, that the local lattice is EVEN 
* in every direction. A runtime check of this is given in topology.
*
* indices : a_of_b/d/t_TYPE(i,even/odd,direction,up/down)
*          i=1 to count; eo=0,1; dir=1 to 4; down=0,up=1
*
* x-direction for down:
*    even:
*         1) A=(M, , , , , , , ) and B=( , , , ,M, , , ) x
*         2) C=(A,A,A,A)         and D=(B,B,B,B) xy
*         3) E=(C,D,C,D,C,D,C,D) and F=(D,C,D,C,D,C,D,C) xyz
*         4) Vector=(E,F,E,F,E,F,E,F) xyzt
*    odd :   A=>B and B=>A
*
*            for up : A=( , , , , , , ,M) and B=( , , ,M, , , , )
*

      if (dimproc(1).gt.1) then

         jump=0
         iter0=0
         iter1=0
         do  z=1,extension(3)
            do  y=1,extension(2)
               if (mod((y+z),2).eq.0) then 
                  iter0=iter0+1
                  a_of_d_TYPE(iter0,0,1,0)=jump*extent
                  a_of_b_TYPE(iter0,0,1,0)=1
                  a_of_t_TYPE(iter0,0,1,0)=matrixtype

                  a_of_d_TYPE(iter0,1,1,1)=(jump+ownhalfvol(1)-1)*extent
                  a_of_b_TYPE(iter0,1,1,1)=1
                  a_of_t_TYPE(iter0,1,1,1)=matrixtype
               else
                  iter1=iter1+1
                  a_of_d_TYPE(iter1,1,1,0)=jump*extent
                  a_of_b_TYPE(iter1,1,1,0)=1
                  a_of_t_TYPE(iter1,1,1,0)=matrixtype

                  a_of_d_TYPE(iter1,0,1,1)=(jump+ownhalfvol(1)-1)*extent
                  a_of_b_TYPE(iter1,0,1,1)=1
                  a_of_t_TYPE(iter1,0,1,1)=matrixtype
               endif
               jump=jump+ownhalfvol(1)
            end do
         end do

* to send the even points down in x-direction

         CALL MPI_TYPE_STRUCT(iter0,a_of_b_TYPE(1,0,1,0),
     &                              a_of_d_TYPE(1,0,1,0),
     &                              a_of_t_TYPE(1,0,1,0),
     &        u_VECTOR_dnTYPE(1,0),ierror)
         CALL MPI_TYPE_COMMIT(u_VECTOR_dnTYPE(1,0),ierror)
            
* to send the odd points down in x-direction
         
         CALL MPI_TYPE_STRUCT(iter0,a_of_b_TYPE(1,1,1,0),
     &                              a_of_d_TYPE(1,1,1,0),
     &                              a_of_t_TYPE(1,1,1,0),
     &        u_VECTOR_dnTYPE(1,1),ierror)
         CALL MPI_TYPE_COMMIT(u_VECTOR_dnTYPE(1,1),ierror)
         
* to send the even points up in x-direction
         
         CALL MPI_TYPE_STRUCT(iter0,a_of_b_TYPE(1,0,1,1),
     &                              a_of_d_TYPE(1,0,1,1),
     &                              a_of_t_TYPE(1,0,1,1),
     &        u_VECTOR_upTYPE(1,0),ierror)
         CALL MPI_TYPE_COMMIT(u_VECTOR_upTYPE(1,0),ierror)
            
* to send the odd points up in x-direction
         
         CALL MPI_TYPE_STRUCT(iter0,a_of_b_TYPE(1,1,1,1),
     &                              a_of_d_TYPE(1,1,1,1),
     &                              a_of_t_TYPE(1,1,1,1),
     &        u_VECTOR_upTYPE(1,1),ierror)
         CALL MPI_TYPE_COMMIT(u_VECTOR_upTYPE(1,1),ierror)


      endif
****************************************************************************
*
* y-direction (for down)
*
*    even and odd :
*         1) A=(M,M,M,M)
*         2) B=(A,...,A,...,A,...,A,...,A,...,A,...,A,...,A,...) xy
*         3) C=(B,B,B,B,B,B,B,B) xyz
*         4) Vector=(C,C,C,C,C,C,C,C) xyzt
*

      if(dimproc(2).gt.1) then

         jump=0
         iter=0
         do z=1,extension(3)
            do  x=1,ownhalfvol(1)
               iter=iter+1
               a_of_d_TYPE(iter,0,2,0)=jump*extent
               a_of_b_TYPE(iter,0,2,0)=1
               a_of_t_TYPE(iter,0,2,0)=matrixtype

               a_of_d_TYPE(iter,0,2,1)=
     &                    (jump + ownhalfvol(2) - ownhalfvol(1))*extent
               a_of_b_TYPE(iter,0,2,1)=1
               a_of_t_TYPE(iter,0,2,1)=matrixtype
               jump=jump+1
            end do
            jump=jump-ownhalfvol(1)+ownhalfvol(2)
         end do

* to send the even points in y-direction down

         CALL MPI_TYPE_STRUCT(iter,a_of_b_TYPE(1,0,2,0),
     &                             a_of_d_TYPE(1,0,2,0),
     &                             a_of_t_TYPE(1,0,2,0),
     &        u_VECTOR_dnTYPE(2,0),ierror)
         CALL MPI_TYPE_COMMIT(u_VECTOR_dnTYPE(2,0),ierror)

* to send the odd points in y-direction down

         u_VECTOR_dnTYPE(2,1)=u_VECTOR_dnTYPE(2,0)

* to send the even points in y-direction up

         CALL MPI_TYPE_STRUCT(iter,a_of_b_TYPE(1,0,2,1),
     &                             a_of_d_TYPE(1,0,2,1),
     &                             a_of_t_TYPE(1,0,2,1),
     &        u_VECTOR_upTYPE(2,0),ierror)
         CALL MPI_TYPE_COMMIT(u_VECTOR_upTYPE(2,0),ierror)

* to send the odd points in y-direction up

         u_VECTOR_upTYPE(2,1)=u_VECTOR_upTYPE(2,0)

      endif
****************************************************************************
*
* z-direction
*
*    even and odd :
*         1) A=(M,M,M,M)
*         2) B=(A,A,A,A,A,A,A,A) xy
*         3) C=(B, , , , , , , ) xyz
*         4) Vector=(C,C,C,C,C,C,C,C) xyzt
*

      if (dimproc(3).gt.1) then

         jump=0
         iter=0
         do y=1,extension(2)
            do x=1,ownhalfvol(1)
               iter=iter+1
               a_of_d_TYPE(iter,0,3,0)=jump*extent
               a_of_b_TYPE(iter,0,3,0)=1
               a_of_t_TYPE(iter,0,3,0)=matrixtype
                
               a_of_d_TYPE(iter,0,3,1)=
     &                    (jump - ownhalfvol(2) + ownhalfvol(3))*extent
               a_of_b_TYPE(iter,0,3,1)=1
               a_of_t_TYPE(iter,0,3,1)=matrixtype
               jump=jump+1
            end do
         end do



* to send the even points in z-direction down

         CALL MPI_TYPE_STRUCT(iter,a_of_b_TYPE(1,0,3,0),
     &                             a_of_d_TYPE(1,0,3,0),
     &                             a_of_t_TYPE(1,0,3,0),
     &        u_VECTOR_dnTYPE(3,0),ierror)
         CALL MPI_TYPE_COMMIT(u_VECTOR_dnTYPE(3,0),ierror)

* to send the odd points in z-direction down

         u_VECTOR_dnTYPE(3,1)=u_VECTOR_dnTYPE(3,0)

* to send the even points in z-direction up

         CALL MPI_TYPE_STRUCT(iter,a_of_b_TYPE(1,0,3,1),
     &                             a_of_d_TYPE(1,0,3,1),
     &                             a_of_t_TYPE(1,0,3,1),
     &        u_VECTOR_upTYPE(3,0),ierror)
         CALL MPI_TYPE_COMMIT(u_VECTOR_upTYPE(3,0),ierror)

* to send the odd points in z-direction up

         u_VECTOR_upTYPE(3,1)=u_VECTOR_upTYPE(3,0)

      endif


***********************************************************************
*                                                                     *
*            ******   the adjoint higgs    *****                      *
*                                                                     *
***********************************************************************
*
*4) to transfer a higgs field
*

      ! CHANGEME this is really an abuse this npdim=nc^2 - 1
      CALL MPI_TYPE_CONTIGUOUS(npdim,MPI_REAL,higgstype,ierror)
      CALL MPI_TYPE_COMMIT(higgstype,ierror)
      CALL MPI_TYPE_EXTENT(higgstype,extent2,ierror)

      CALL MPI_ADDRESS(a(1,1,0),ad1,ierror)
      CALL MPI_ADDRESS(a(1,2,0),ad2,ierror)
      extent=ad2-ad1

*      if(myrank==root)then
*         write (*,*) 'higgstype'
*         write (*,*) 'MPI_EXTENT=',extent2
*         write (*,*) 'EXTENT=',extent
*      end if

***************************************************************************
*
*5) to transfer xy-slize of higgs field to read and write from/to disk
*
***************************************************************************

      CALL MPI_TYPE_CONTIGUOUS( ownhalfvol(2),
     &                          higgstype,load_a_type,ierror)
      CALL MPI_TYPE_COMMIT(load_a_type,ierror)

      CALL MPI_TYPE_CONTIGUOUS( ownhalfvol(2),
     &                          higgstype,save_a_type,ierror)
      CALL MPI_TYPE_COMMIT(save_a_type,ierror)

***************************************************************************
*
*6) to send a selected part of the higgs field 
*   ((usefull),.....,(usefull),.....,(usefull),.....)
*   and to receive it in a contiguous buffer
*   ((usefull),(usefull),(usefull))

      ! why have count here and and not ownhafvol(2)?
      CALL MPI_TYPE_CONTIGUOUS(count,higgstype,a_CONT_TYPE,ierror)
      CALL MPI_TYPE_COMMIT(a_CONT_TYPE,ierror)

*
*consider first the movement down and then up
*for this two different types of transfer patterns are necessary
*
*indices : a_of_b/d/t_TYPE(i,even/odd,direction,up/down)
*          i=1 to count; eo=0,1; dir=1 to 4; down=0,up=1
*

*
*x-direction for down:
*    even:
*         1) A=(psi, , , , , , , ) and B=( , , , ,psi, , , ) x
*         2) C=(A,A,A,A)         and D=(B,B,B,B) xy
*         3) E=(C,D,C,D,C,D,C,D) and F=(D,C,D,C,D,C,D,C) xyz
*         4) Vector=(E,F,E,F,E,F,E,F) xyzt
*    odd :   A=>B and B=>A
*
*            for up : A=( , , , , , , ,psi) and B=( , , ,psi, , , , )
*

      if (dimproc(1).gt.1) then

         jump=0
         iter0=0
         iter1=0
         do z=1,extension(3)
            do y=1,extension(2)
               if (mod((y+z),2).eq.0) then 
                  iter0=iter0+1
                  a_of_d_TYPE(iter0,0,1,0)=jump*extent
                  a_of_b_TYPE(iter0,0,1,0)=1
                  a_of_t_TYPE(iter0,0,1,0)=higgstype

                  a_of_d_TYPE(iter0,1,1,1)=(jump + ownhalfvol(1) - 1)*extent
                  a_of_b_TYPE(iter0,1,1,1)=1
                  a_of_t_TYPE(iter0,1,1,1)=higgstype
               else
                  iter1=iter1+1
                  a_of_d_TYPE(iter1,1,1,0)=jump*extent
                  a_of_b_TYPE(iter1,1,1,0)=1
                  a_of_t_TYPE(iter1,1,1,0)=higgstype

                  a_of_d_TYPE(iter1,0,1,1)=(jump + ownhalfvol(1) - 1)*extent
                  a_of_b_TYPE(iter1,0,1,1)=1
                  a_of_t_TYPE(iter1,0,1,1)=higgstype
               endif
               jump=jump+ownhalfvol(1)
            end do
         end do
            
* to send the even points in x-direction down   

         CALL MPI_TYPE_STRUCT(iter0,a_of_b_TYPE(1,0,1,0),
     &                              a_of_d_TYPE(1,0,1,0),
     &                              a_of_t_TYPE(1,0,1,0),
     &        a_VECTOR_dnTYPE(1,0),ierror)
         CALL MPI_TYPE_COMMIT(a_VECTOR_dnTYPE(1,0),ierror)

* to send the odd points in x-direction down   
         
         CALL MPI_TYPE_STRUCT(iter0,a_of_b_TYPE(1,1,1,0),
     &                              a_of_d_TYPE(1,1,1,0),
     &                              a_of_t_TYPE(1,1,1,0),
     &        a_VECTOR_dnTYPE(1,1),ierror)
         CALL MPI_TYPE_COMMIT(a_VECTOR_dnTYPE(1,1),ierror)

* to send the even points in x-direction up   

         CALL MPI_TYPE_STRUCT(iter0,a_of_b_TYPE(1,0,1,1),
     &                              a_of_d_TYPE(1,0,1,1),
     &                              a_of_t_TYPE(1,0,1,1),
     &        a_VECTOR_upTYPE(1,0),ierror)
         CALL MPI_TYPE_COMMIT(a_VECTOR_upTYPE(1,0),ierror)
         
* to send the odd points in x-direction up   

         CALL MPI_TYPE_STRUCT(iter0,a_of_b_TYPE(1,1,1,1),
     &                              a_of_d_TYPE(1,1,1,1),
     &                              a_of_t_TYPE(1,1,1,1),
     &        a_VECTOR_upTYPE(1,1),ierror)
         CALL MPI_TYPE_COMMIT(a_VECTOR_upTYPE(1,1),ierror)


      endif
*************************************************************************
*
* y-direction (for down)
*
*    even and odd :
*         1) A=(psi,psi,psi,psi)
*         2) B=(A,...,A,...,A,...,A,...,A,...,A,...,A,...,A,...) xy
*         3) C=(B,B,B,B,B,B,B,B) xyz
*         4) Vector=(C,C,C,C,C,C,C,C) xyzt
*

      if (dimproc(2).gt.1) then

         jump=0
         iter=0
         do z=1,extension(3)
            do x=1,ownhalfvol(1)
               iter=iter+1
               a_of_d_TYPE(iter,0,2,0)=jump*extent
               a_of_b_TYPE(iter,0,2,0)=1
               a_of_t_TYPE(iter,0,2,0)=higgstype

               a_of_d_TYPE(iter,0,2,1)=
     &                    (jump - ownhalfvol(1) + ownhalfvol(2))*extent
               a_of_b_TYPE(iter,0,2,1)=1
               a_of_t_TYPE(iter,0,2,1)=higgstype
               jump=jump+1
            end do
            jump=jump-ownhalfvol(1)+ownhalfvol(2)
         end do
         
* to send the even points in y-direction down   

         CALL MPI_TYPE_STRUCT(iter,a_of_b_TYPE(1,0,2,0),
     &                             a_of_d_TYPE(1,0,2,0),
     &                             a_of_t_TYPE(1,0,2,0),
     &        a_VECTOR_dnTYPE(2,0),ierror)
         CALL MPI_TYPE_COMMIT(a_VECTOR_dnTYPE(2,0),ierror)

* to send the odd points in y-direction down 

         a_VECTOR_dnTYPE(2,1)=a_VECTOR_dnTYPE(2,0)

* to send the even points in y-direction up   

         CALL MPI_TYPE_STRUCT(iter,a_of_b_TYPE(1,0,2,1),
     &                             a_of_d_TYPE(1,0,2,1),
     &                             a_of_t_TYPE(1,0,2,1),
     &        a_VECTOR_upTYPE(2,0),ierror)
         CALL MPI_TYPE_COMMIT(a_VECTOR_upTYPE(2,0),ierror)

* to send the odd points in y-direction up   

         a_VECTOR_upTYPE(2,1)=a_VECTOR_upTYPE(2,0)

      endif


*
*z-direction
*
*    even and odd :
*         1) A=(psi,psi,psi,psi)
*         2) B=(A,A,A,A,A,A,A,A) xy
*         3) C=(B, , , , , , , ) xyz
*         4) Vector=(C,C,C,C,C,C,C,C) xyzt
*

      if (dimproc(3).gt.1) then

         jump=0
         iter=0
         do y=1,extension(2)
            do x=1,ownhalfvol(1)
               iter=iter+1
               a_of_d_TYPE(iter,0,3,0)=jump*extent
               a_of_b_TYPE(iter,0,3,0)=1
               a_of_t_TYPE(iter,0,3,0)=higgstype

               a_of_d_TYPE(iter,0,3,1)=
     &                    (jump - ownhalfvol(2) + ownhalfvol(3))*extent
               a_of_b_TYPE(iter,0,3,1)=1
               a_of_t_TYPE(iter,0,3,1)=higgstype
               jump=jump+1
            end do
         end do

* to send the even points in z-direction down

         CALL MPI_TYPE_STRUCT(iter,a_of_b_TYPE(1,0,3,0),
     &                             a_of_d_TYPE(1,0,3,0),
     &                             a_of_t_TYPE(1,0,3,0),
     &        a_VECTOR_dnTYPE(3,0),ierror)
         CALL MPI_TYPE_COMMIT(a_VECTOR_dnTYPE(3,0),ierror)

* to send the odd points in z-direction down

         a_VECTOR_dnTYPE(3,1)=a_VECTOR_dnTYPE(3,0)

* to send the even points in z-direction up

         CALL MPI_TYPE_STRUCT(iter,a_of_b_TYPE(1,0,3,1),
     &                             a_of_d_TYPE(1,0,3,1),
     &                             a_of_t_TYPE(1,0,3,1),
     &        a_VECTOR_upTYPE(3,0),ierror)
         CALL MPI_TYPE_COMMIT(a_VECTOR_upTYPE(3,0),ierror)

* to send the odd points in z-direction up

         a_VECTOR_upTYPE(3,1)=a_VECTOR_upTYPE(3,0)

      endif
  
      return

      end
