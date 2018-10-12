      SUBROUTINE GEM_TCL_INIT
c
c_begin_doc
c  RCS ID string
c  $Id: gem_tcl_init.F,v 1.2.2.2 2010/04/23 14:43:48 claseg6 Exp $
c
c  Documentation for subroutine GEM_TCL_INIT
c
c  Purpose:
c  --------
c
c  Define variables that have to be accessible to the user via the
c  TCL interface.
c
c  Input Parameters:  (Name - Type - Meaning)
c  ----------------
c
c   None
c
c  Output Parameters:  (Name - Type - Meaning)
c  -----------------
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

c     
c
c_end_doc
c
      IMPLICIT NONE
      SAVE
c
c_begin_inc
c  include files :
c  ---------------------
#include "gem_tcl.inc"
c_end_inc
c
c_begin_var
c  input/output variables:
c  -----------------------
c
c  Local pre-defined variables:
c  ---------------------------
c  RCS information: 
c
c  Local User defined variables:
c  -----------------------------
c_end_var
c
c  executable code for routine DC_TCL_INIT:
c----6----------------------------------------------------------------72
c
c set tcl defaults
c      call gem_set_default()
c      write(37,*)'now in gem_tcl_init.F'
*** tcl_linkvar to your heart's content here
***      call tcl_linkvar('gem_efficiency',gem_efficiency,itcl_double)
      RETURN
      END
c
c------------------------------------------------------------------------------




