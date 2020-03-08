#######################################################################
#  This is the makefile to create a library of the test matrix
#  generators used in LAPACK.  The files are organized as follows:
#
#     SCATGEN  -- Auxiliary routines called from both REAL and COMPLEX
#     DZATGEN  -- Auxiliary routines called from both DOUBLE PRECISION
#                 and COMPLEX*16
#     SMATGEN  -- Single precision real matrix generation routines
#     CMATGEN  -- Single precision complex matrix generation routines
#     DMATGEN  -- Double precision real matrix generation routines
#     ZMATGEN  -- Double precision complex matrix generation routines
#
#  The library can be set up to include routines for any combination
#  of the four precisions.  To create or add to the library, enter make
#  followed by one or more of the precisions desired.  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments creates a library of all four precisions.
#  The library is called
#       tmglib.a
#  and is created at the LAPACK directory level.
#
#  To remove the object files after the library is created, enter
#       make clean
#  On some systems, you can force the source files to be recompiled by
#  entering (for example)
#       make single FRC=FRC
#
#######################################################################
 
set(SCATGEN  slatm1.c slaran.c slarnd.c)

set(SMATGEN  slatms.c slatme.c slatmr.c slatmt.c 
   slagge.c slagsy.c slakf2.c slarge.c slaror.c slarot.c slatm2.c 
   slatm3.c slatm5.c slatm6.c slatm7.c slahilb.c)

set(CMATGEN  clatms.c clatme.c clatmr.c clatmt.c 
   clagge.c claghe.c clagsy.c clakf2.c clarge.c claror.c clarot.c 
   clatm1.c clarnd.c clatm2.c clatm3.c clatm5.c clatm6.c clahilb.c)

set(DZATGEN  dlatm1.c dlaran.c dlarnd.c)

set(DMATGEN  dlatms.c dlatme.c dlatmr.c dlatmt.c 
   dlagge.c dlagsy.c dlakf2.c dlarge.c dlaror.c dlarot.c dlatm2.c 
   dlatm3.c dlatm5.c dlatm6.c dlatm7.c dlahilb.c)

set(ZMATGEN  zlatms.c zlatme.c zlatmr.c zlatmt.c 
  zlagge.c zlaghe.c zlagsy.c zlakf2.c zlarge.c zlaror.c zlarot.c 
  zlatm1.c zlarnd.c zlatm2.c zlatm3.c zlatm5.c zlatm6.c zlahilb.c)

set(ALLOBJ ${SMATGEN} ${CMATGEN} ${SCATGEN} ${DMATGEN} ${ZMATGEN}
  ${DZATGEN})
if(BUILD_SINGLE)
  set(ALLOBJ $(SMATGEN) $(SCATGEN))
endif()
if(BUILD_DOUBLE)
  set(ALLOBJ $(DMATGEN) $(DZATGEN))
endif()
if(BUILD_COMPLEX)
  set(ALLOBJ  $(CMATGEN) $(SCATGEN))
endif()
if(BUILD_COMPLEX16)
  set(ALLOBJ $(ZMATGEN) $(DZATGEN))
endif()
add_library(tmglib ${ALLOBJ} )

