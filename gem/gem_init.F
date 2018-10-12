      SUBROUTINE GEM_INIT(irun)
c
c_begin_doc
c  RCS ID string
c  $Id: gem_init.F,v 1.3.4.4 2013/05/06 15:27:38 claseg6 Exp $
c
c  Documentation for subroutine GEM_INIT
c
c  Purpose: GEM supplied routined called at the BEGINNING of EVERY run
c  --------
c
c  Calling Sequence:
c  ----------------
c
c  Input Parameters:  NONE
c  ----------------
c
c  Output Parameters:  NONE
c  -----------------
c
c  Called from:
c  ------------
c
c  Other routines:
c  ---------------
c
c  Notes:
c  ------
c
c  Author: gabriel niculescu, 2005
c  -------
c
c  Major revisions:  created - 2005 GN.
c  ----------------
c_end_doc
c
      IMPLICIT NONE
      SAVE
c
c_begin_inc
c  include files :
c  ---------------------
c BOS common block  uncomment the next line for BOS include file
#include "bcs.inc"
c_end_inc
c
c_begin_var
c  input/output variables:
c  -----------------------
c      
c
c  Local pre-defined variables:
c  ---------------------------
c  RCS information: 
      CHARACTER*(*)  CFILE, CREVIS, CSTATE, CDATE, CAUTHO 
      PARAMETER (CFILE=  '$RCSfile: gem_init.F,v $')
      PARAMETER (CREVIS= '$Revision: 1.3.4.4 $')
      PARAMETER (CSTATE= '$State: Exp $')
      PARAMETER (CDATE=  '$Date: 2013/05/06 15:27:38 $')
      PARAMETER (CAUTHO= '$Author: claseg6 $')
c  Module information:
      CHARACTER*(*)  CRNAME, CRAUTH
      CHARACTER*132  CRMESS
      PARAMETER (CRNAME='GEM_INIT')
      PARAMETER (CRAUTH='Gabriel Niculescu')
c
c  Local gem defined variables:
c  ----------------------------
      logical ok
      integer ind,mlink,evtype
      integer irun,ievnt,irec
      integer nwrite

c Common group for output
c Note: this must match ana/tpccnt.inc!
      integer tpc_hits_nb
      integer chan_nb(1000)
      integer adc(1000)
      real    cadc(1000)
      integer tdc(1000)
      integer track_nb(1000)
      common /tpc_hits/ tpc_hits_nb,chan_nb,adc,cadc,tdc,track_nb

      tpc_hits_nb = 0

c_end_var
c
c  executable code for routine GEM_BRUN:
c  ------------------------------------
c
      IF (NWRITE .LT. 1) THEN
        NWRITE = NWRITE + 1
        CRMESS='Reading GEM calibration constants'
        CALL RECMES(CRNAME,'I',CRMESS)
      ENDIF

      write(*,*) 'gem_init: Reading GEM calibration constants'
      call gem_init_c()

      RETURN
      END
