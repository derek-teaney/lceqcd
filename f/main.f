      program main


      implicit none

      INCLUDE 'input_parameter.inc'
      INCLUDE 'measurement.inc'
      
      integer isweep, meas, iter_or1, corr_flag, mstart
      real    acc1, acc2, acc3, da1

      call initialize

* set the flag for correlation functions measurements

      corr_flag=1

* loop to  thermalize the system

      do isweep=1,100
         call ahiggs2(0,acc2)
         call ahiggs3(0,acc3)
         call agauge(0,acc1)
         call renorm
      enddo
        
       
      acco=0.0
      acch=0.0
      accg=0.0
      meas=0.0
      iter_or=0
      da=0.0

*  do loop over MC iteration and perform measurements

      do isweep=1,nsweep
         call ahiggs2(1,acc2)
         call ahiggs3(1,acc3)
         call agauge(1,acc1)
         call renorm
         acco=acco+acc3
         acch=acch+acc2
         accg=accg+acc1

*  perform measurements of local quantities
*  colect measurements on root processor and save on disk after  100 sweep

         call meas1(isweep) 

**                                    *
**  correlation function measurements * 
**                                    *
            
         mstart=int(confstart/imeas)      
         
                 
         if(mod(isweep,imeas).eq.0)then
            meas=meas+1
            
            if(corr_flag.eq.1) then
             call landau_gauge(da1,iter_or1)
             da=da+da1
             iter_or=iter_or+iter_or1
             call meas2(meas)
            endif
            
            if(corr_flag.eq.2) then
             call polyakov(meas,mstart)
             call scalarcorr(meas,mstart)
            endif

            if(corr_flag.eq.3) then
             call landau_gauge(da1,iter_or1)
             da=da+da1
             iter_or=iter_or+iter_or1
             call meas2(meas)
             call polyakov(meas,mstart)
            endif

         endif
 
      enddo 

      call finalize

      end program main
































