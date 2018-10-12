      SUBROUTINE GEM_BRUN(irun)
c
c_begin_doc
c  RCS ID string
c  $Id: gem_brun.F,v 1.2.4.4 2013/05/06 15:27:37 claseg6 Exp $
c
c  Documentation for subroutine GEM_BRUN
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
      logical ok
      integer ind,mlink,evtype
      integer irun,ievnt,irec
      integer nwrite

c  ---------------------------
c  RCS information: 
      CHARACTER*(*)  CFILE, CREVIS, CSTATE, CDATE, CAUTHO 
      PARAMETER (CFILE=  '$RCSfile: gem_brun.F,v $')
      PARAMETER (CREVIS= '$Revision: 1.2.4.4 $')
      PARAMETER (CSTATE= '$State: Exp $')
      PARAMETER (CDATE=  '$Date: 2013/05/06 15:27:37 $')
      PARAMETER (CAUTHO= '$Author: claseg6 $')
c  Module information:
      CHARACTER*(*)  CRNAME, CRAUTH
      CHARACTER*132  CRMESS
      PARAMETER (CRNAME='GEM_BRUN')
      PARAMETER (CRAUTH='Gabriel Niculescu')
*      external gem_init_c_
c
c  Local gem defined variables:
c  ----------------------------

c_end_var
c
c  executable code for routine GEM_BRUN:
c  ------------------------------------
c
      IF (NWRITE .LT. 1) THEN
        NWRITE = NWRITE + 1
        CRMESS='gem_brun.F:  Reading GEM calibration constants'
        CALL RECMES(CRNAME,'I',CRMESS)
      ENDIF  
c      write(*,*)'\e[0;32m gem_brun: Initializing GEM data structures \e[0;37m'
      write(*,*) 'gem_brun.F: Initializing GEM data structure'

      call gem_tcl_init()
c     call gem_init_c()
      call gem_read_database(irun)
c     initialize geometry
*   read stuff from map, etc
c
      RETURN
      END
c
c------------------------------------------------------------------------------




