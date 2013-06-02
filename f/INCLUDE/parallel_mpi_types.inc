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
*     Include file for Message Passing variables                   *
*                                                                  *
*     the MESSAGE PASSING TYPES                                    *
*                                                                  *
********************************************************************


* internal identifier for transfered MPI-Data
      INTEGER matrixtype,morematrixtype,
     &        higgstype,morehiggstype,
     &	      load_u_type,save_u_type,
     &	      load_a_type,save_a_type,
     &        paramettype,configtype,latconftype,
     &        u_CONT_TYPE,a_CONT_TYPE,
     &        u_VECTOR_dnTYPE(npdim,0:1),u_VECTOR_upTYPE(npdim,0:1),
     &        a_VECTOR_dnTYPE(npdim,0:1),a_VECTOR_upTYPE(npdim,0:1)

* communication handles
      INTEGER send_request,send_request1,send_request2,
     1        send_status(MPI_STATUS_SIZE),send_status1(MPI_STATUS_SIZE),
     1        send_status2(MPI_STATUS_SIZE)

      INTEGER rec_request,rec_request1,rec_request2,
     1        rec_status(MPI_STATUS_SIZE),rec_status1(MPI_STATUS_SIZE),
     1        rec_status2(MPI_STATUS_SIZE)

      INTEGER save_rec_status(MPI_STATUS_SIZE)

      INTEGER load_rec_status(MPI_STATUS_SIZE)




***********************************************************************
*
*                     THE COMMON BLOCK
*
***********************************************************************

      common /vartypes / matrixtype,morematrixtype,
     &                   higgstype,morehiggstype,
     &                   load_u_type,save_u_type,
     &                   load_a_type,save_a_type,
     &                   paramettype,configtype,latconftype,
     &                   u_CONT_TYPE,a_CONT_TYPE,
     &                   u_VECTOR_dnTYPE,u_VECTOR_upTYPE,
     &                   a_VECTOR_dnTYPE,a_VECTOR_upTYPE




