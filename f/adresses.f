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
*         adresses for the neighbouring sites                      *
*                                                                  *
********************************************************************

      subroutine adresses

      implicit none
*
*including all global MESSAGE PASSING features and variables 
*
      INCLUDE 'parallel_parameter.inc'
*
*including global variables, parameters and common blocks
*
      INCLUDE 'parameter.inc'
      INCLUDE 'pointer.inc'
*
*definition of the local VARIABLES
*
      integer i0,i1,eo,x,y,z,mu,itail(0:1,-nmu:nmu)

***********************************************************************
* 
* the sites are ordered seriell for even and for odd numbers 
* (x+y+z+1):
* so even goes from 1 to ownhalfvol and odd too
* 
*      now construct the neighbouring site in the mu direction
* 
* 
* for y,z and t it is easier than for the x direction :
* 
* if you change one step in y,z or t you have to add extension(1)/2, 
* the xy volume or the xyz volume and change even to odd and vice versa 
* because you move to the dual lattice 
* but be careful at the surface of the lattice => wrap around for
* periodic boundaries in case only one process is existing in that
* direction or you have to link to the buffer where the sent data
* will be stored
* 
* 
* now the x direction
* 
* two neighbouring sites can differ in the lattice coordinate i 
* about one unit or they can be the same because of the duality 
* of the lattice 
* that's why it is a bit tricky with x
* they allways differ in the even/odd coordinate, so this 
* coordinate must be changed for a single step
* 
* for i you have to look whether x+y+z+t is even or odd and if y+z+t 
* is even or odd and then it depends on the direction (up or down) 
* of the step wether i is changed by one or not 
* 
* for example if x+y+z+t is even and y+z+t too (take 2+2+1+1) 
* and you want to go up in x direction then even changes 
* to odd and i to i+1 
* (if i not equal ownhalfvol(1), then you have to go to the buffer or in
* case of periodic boundaries because of only one process being in 
* that dimension you have to wrap around => aint and mod)
* 
* siteup/dn ( local site, even-odd, direction to go)
* 
* 

      i0=0
      i1=0

*
*starting coordinate of the buffer for sent data is always nsiteshalf+1
*
      do  eo=0,1
         do  mu=-nmu,nmu
            itail(eo,mu)=nsiteshalf
         end do
      end do

      do  z=lowbound(3),upbound(3)
         do  y=lowbound(2),upbound(2)
            do  x=lowbound(1),upbound(1)

               if ((mod((x+y+z),2)).eq.1) then

                  i0=i0+1

* x-direction even

* only one process in x-direction => periodic boundaries

                  if (dimproc(1).eq.1) then
                        
                     if ((mod((y+z),2)).eq.1) then
                        siteup(i0,0,1)=mod(i0,ownhalfvol(1))+1+
     &                       aint((real(i0-1))/(real
     &                       (ownhalfvol(1))))*ownhalfvol(1)
                        sitedn(i0,0,1)=i0
                     else
                        siteup(i0,0,1)=i0
                        sitedn(i0,0,1)=mod((i0+ownhalfvol(1)-2),
     &                       ownhalfvol(1))+1+
     &                       aint((real(i0-1))/(real
     &                       (ownhalfvol(1))))*ownhalfvol(1)
                     endif

* more than one process 
                     
                  else

                     if ((mod((y+z),2)).eq.1) then
                        siteup(i0,0,1)=i0+1
                        sitedn(i0,0,1)=i0
                     else
                        siteup(i0,0,1)=i0
                        sitedn(i0,0,1)=i0-1
                     endif
                     
                     if (x.eq.lowbound(1)) then
                        itail(0,-1)=itail(0,-1)+1
                        sitedn(i0,0,1)=itail(0,-1)
                     endif
                     
                     if (x.eq.upbound(1)) then
                        itail(0,1)=itail(0,1)+1
                        siteup(i0,0,1)=itail(0,1)
                     endif

                  endif


* y-direction even

* only one process in y-direction => periodic boundaries

                     
                  if (dimproc(2).eq.1) then

                     siteup(i0,0,2)=mod((i0-1+ownhalfvol(1)),
     &                    ownhalfvol(2))+1+aint((real(i0-1))/
     &                    (real(ownhalfvol(2))))*ownhalfvol(2)
                     
                     sitedn(i0,0,2)=mod((i0-1+ownhalfvol(1)*
     &                    (extension(2)-1)),ownhalfvol(2))+1+
     &                    aint((real(i0-1))/(real(ownhalfvol(2))))*
     &                    ownhalfvol(2)

* more than one process

                  else

                     siteup(i0,0,2)=i0+ownhalfvol(1)
                     sitedn(i0,0,2)=i0-ownhalfvol(1)

                     if (y.eq.lowbound(2)) then
                        itail(0,-2)=itail(0,-2)+1
                        sitedn(i0,0,2)=itail(0,-2)
                     endif

                     if (y.eq.upbound(2)) then
                        itail(0,2)=itail(0,2)+1
                        siteup(i0,0,2)=itail(0,2)
                     endif
                     
                  endif


* z-direction even

* only one process in z-direction => periodic boundaries

                  if (dimproc(3).eq.1) then

                     siteup(i0,0,3)=mod((i0-1+ownhalfvol(2)),
     &                    ownhalfvol(3))+1+
     &                    aint((real(i0-1))/(real(ownhalfvol(3))))*
     &                    ownhalfvol(3)
                     
                     sitedn(i0,0,3)=mod((i0-1+ownhalfvol(2)*
     &                    (extension(3)-1)),ownhalfvol(3))+1+
     &                    aint((real(i0-1))/(real(ownhalfvol(3))))*
     &                    ownhalfvol(3)

*more than one process 

                  else

                     siteup(i0,0,3)=i0+ownhalfvol(2)
                     sitedn(i0,0,3)=i0-ownhalfvol(2)
                     
                     if (z.eq.lowbound(3)) then
                        itail(0,-3)=itail(0,-3)+1
                        sitedn(i0,0,3)=itail(0,-3)
                     endif

                     if (z.eq.upbound(3)) then
                        itail(0,3)=itail(0,3)+1
                        siteup(i0,0,3)=itail(0,3)
                     endif
                        
                  endif

               else
                  i1=i1+1
                  
* x-direction odd 

* only one process in x-direction => periodic boundaries

                  if (dimproc(1).eq.1) then

                     if ((mod((y+z),2)).eq.0) then
                        siteup(i1,1,1)=mod(i1,ownhalfvol(1))+1+
     &                       aint((real(i1-1))/(real
     &                       (ownhalfvol(1))))*ownhalfvol(1)
                        sitedn(i1,1,1)=i1
                     else
                        siteup(i1,1,1)=i1
                        sitedn(i1,1,1)=mod((i1+ownhalfvol(1)-2),
     &                       ownhalfvol(1))+1+
     &                       aint((real(i1-1))/(real
     &                       (ownhalfvol(1))))*ownhalfvol(1)
                     endif
 
* more than one process

                  else
  
                     if ((mod((y+z),2)).eq.0) then
                        siteup(i1,1,1)=i1+1
                        sitedn(i1,1,1)=i1
                     else
                        siteup(i1,1,1)=i1
                        sitedn(i1,1,1)=i1-1
                     endif
                        
                     if (x.eq.lowbound(1)) then
                        itail(1,-1)=itail(1,-1)+1
                        sitedn(i1,1,1)=itail(1,-1)
                     endif
                        
                     if (x.eq.upbound(1)) then
                        itail(1,1)=itail(1,1)+1
                        siteup(i1,1,1)=itail(1,1)
                     endif
                     
                  endif


* y-direction odd

* only one process in y-direction => periodic boundaries

                  if (dimproc(2).eq.1) then

                     siteup(i1,1,2)=mod((i1-1+ownhalfvol(1)),
     &                    ownhalfvol(2))+1+aint((real(i1-1))/
     &                    (real(ownhalfvol(2))))*ownhalfvol(2)
                     
                     sitedn(i1,1,2)=mod((i1-1+ownhalfvol(1)*
     &                    (extension(2)-1)),ownhalfvol(2))+1+
     &                    aint((real(i1-1))/(real(ownhalfvol(2))))*
     &                    ownhalfvol(2)

* more than one process 

                  else
                     
                     siteup(i1,1,2)=i1+ownhalfvol(1)
                     sitedn(i1,1,2)=i1-ownhalfvol(1)
                        
                     if (y.eq.lowbound(2)) then
                        itail(1,-2)=itail(1,-2)+1
                        sitedn(i1,1,2)=itail(1,-2)
                     endif 
                        
                     if (y.eq.upbound(2)) then
                        itail(1,2)=itail(1,2)+1
                        siteup(i1,1,2)=itail(1,2)
                     endif 
 
                  endif


* z-direction odd                       

* only one process in z-direction => periodic boundaries

                  if (dimproc(3).eq.1) then

                     siteup(i1,1,3)=mod((i1-1+ownhalfvol(2)),
     &                    ownhalfvol(3))+1+
     &                    aint((real(i1-1))/(real(ownhalfvol(3))))*
     &                    ownhalfvol(3)

                     sitedn(i1,1,3)=mod((i1-1+ownhalfvol(2)*
     &                    (extension(3)-1)),ownhalfvol(3))+1+
     &                    aint((real(i1-1))/(real(ownhalfvol(3))))*
     &                    ownhalfvol(3)

* more than one process 

                  else

                     siteup(i1,1,3)=i1+ownhalfvol(2)
                     sitedn(i1,1,3)=i1-ownhalfvol(2)
                        
                     if (z.eq.lowbound(3)) then
                        itail(1,-3)=itail(1,-3)+1
                        sitedn(i1,1,3)=itail(1,-3)
                     endif
                        
                     if (z.eq.upbound(3)) then
                        itail(1,3)=itail(1,3)+1
                        siteup(i1,1,3)=itail(1,3)
                     endif 
                     
                  endif


*even/odd
               endif 

            end do
         end do
      end do


      return
      end


	



