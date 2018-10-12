c      include file GEM_tcl
c
c_begin_doc
c  RCS ID string
c  $Id: gem_tcl.inc,v 1.1.1.1.4.2 2010/04/23 14:43:48 claseg6 Exp $
c
c  Documentation for include file GEM_TCL
c
c  Purpose:
c  --------
c  Contains tcl variables for gem.
c
c  Used in routines:
c  ----------------
c  gem_tcl_init
c  gem_getcal_geom
c  gem_xvst_get_dist
c  gem_xvst_get_tdrift
c  Notes:
c  ------
c
c  Author: gabriel niculescu, 2005
c  -------
c
c  Major revisions:  created - 2005 GN.
c  ----------------
c
c_end_doc
c
c_begin_var
c  variables for include file GEM_TCL:

c  Variables contained in GEM_TCL common block
c  ------------------------------------------

      real*8 my_gem_tcl_var
c_end_var      
      COMMON/GEM_TCL/ my_gem_tcl_var
      save /GEM_tcl/
c----6----------------------------------------------------------------72



















