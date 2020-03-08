########################################################################
#  This is the makefile for the eigenvalue test program from LAPACK.
#  The test files are organized as follows:
#
#     AEIGTST -- Auxiliary test routines used in all precisions
#     SCIGTST -- Auxiliary test routines used in REAL and COMPLEX
#     DZIGTST -- Auxiliary test routines used in DOUBLE PRECISION and
#                COMPLEX*16
#     SEIGTST -- Single precision real test routines
#     CEIGTST -- Single precision complex test routines
#     DEIGTST -- Double precision real test routines
#     ZEIGTST -- Double precision complex test routines
#
#  Test programs can be generated for all or some of the four different
#  precisions.  Enter make followed by one or more of the data types
#  desired.  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments creates all four test programs.
#  The executable files are called
#       xeigtsts, xeigtstd, xeigtstc, and xeigtstz
#  and are created in the next higher directory level.
#
#  To remove the object files after the executable files have been
#  created, enter
#       make clean
#  On some systems, you can force the source files to be recompiled by
#  entering (for example)
#       make single FRC=FRC
#
########################################################################

set(AEIGTST  
   alahdg.c 
   alasum.c 
   alasvm.c 
   alareq.c 
   ilaenv.c 
   xerbla.c 
   xlaenv.c 
   chkxer.c)

set(SCIGTST  slafts.c slahd2.c slasum.c slatb9.c sstech.c sstect.c 
   ssvdch.c ssvdct.c ssxt1.c)

set(SEIGTST  schkee.c 
   sbdt01.c sbdt02.c sbdt03.c 
   schkbb.c schkbd.c schkbk.c schkbl.c schkec.c 
   schkgg.c schkgk.c schkgl.c schkhs.c schksb.c schkst.c 
   sckglm.c sckgqr.c sckgsv.c scklse.c 
   sdrges.c sdrgev.c sdrgsx.c sdrgvx.c 
   sdrvbd.c sdrves.c sdrvev.c sdrvgg.c sdrvsg.c 
   sdrvst.c sdrvsx.c sdrvvx.c 
   serrbd.c serrec.c serred.c serrgg.c serrhs.c serrst.c 
   sget02.c sget10.c sget22.c sget23.c sget24.c sget31.c 
   sget32.c sget33.c sget34.c sget35.c sget36.c 
   sget37.c sget38.c sget39.c sget51.c sget52.c sget53.c 
   sget54.c sglmts.c sgqrts.c sgrqts.c sgsvts.c 
   shst01.c slarfy.c slarhs.c slatm4.c slctes.c slctsx.c slsets.c sort01.c 
   sort03.c ssbt21.c ssgt01.c sslect.c sspt21.c sstt21.c 
   sstt22.c ssyt21.c ssyt22.c)

set(CEIGTST  cchkee.c 
   cbdt01.c cbdt02.c cbdt03.c 
   cchkbb.c cchkbd.c cchkbk.c cchkbl.c cchkec.c 
   cchkgg.c cchkgk.c cchkgl.c cchkhb.c cchkhs.c cchkst.c 
   cckglm.c cckgqr.c cckgsv.c ccklse.c 
   cdrges.c cdrgev.c cdrgsx.c cdrgvx.c 
   cdrvbd.c cdrves.c cdrvev.c cdrvgg.c cdrvsg.c 
   cdrvst.c cdrvsx.c cdrvvx.c 
   cerrbd.c cerrec.c cerred.c cerrgg.c cerrhs.c cerrst.c 
   cget02.c cget10.c cget22.c cget23.c cget24.c 
   cget35.c cget36.c cget37.c cget38.c cget51.c cget52.c 
   cget54.c cglmts.c cgqrts.c cgrqts.c cgsvts.c 
   chbt21.c chet21.c chet22.c chpt21.c chst01.c 
   clarfy.c clarhs.c clatm4.c clctes.c clctsx.c clsets.c csbmv.c 
   csgt01.c cslect.c 
   cstt21.c cstt22.c cunt01.c cunt03.c)

set(DZIGTST  dlafts.c dlahd2.c dlasum.c dlatb9.c dstech.c dstect.c 
   dsvdch.c dsvdct.c dsxt1.c)

set(DEIGTST  dchkee.c 
   dbdt01.c dbdt02.c dbdt03.c 
   dchkbb.c dchkbd.c dchkbk.c dchkbl.c dchkec.c 
   dchkgg.c dchkgk.c dchkgl.c dchkhs.c dchksb.c dchkst.c 
   dckglm.c dckgqr.c dckgsv.c dcklse.c 
   ddrges.c ddrgev.c ddrgsx.c ddrgvx.c 
   ddrvbd.c ddrves.c ddrvev.c ddrvgg.c ddrvsg.c 
   ddrvst.c ddrvsx.c ddrvvx.c 
   derrbd.c derrec.c derred.c derrgg.c derrhs.c derrst.c 
   dget02.c dget10.c dget22.c dget23.c dget24.c dget31.c 
   dget32.c dget33.c dget34.c dget35.c dget36.c 
   dget37.c dget38.c dget39.c dget51.c dget52.c dget53.c 
   dget54.c dglmts.c dgqrts.c dgrqts.c dgsvts.c 
   dhst01.c dlarfy.c dlarhs.c dlatm4.c dlctes.c dlctsx.c dlsets.c dort01.c 
   dort03.c dsbt21.c dsgt01.c dslect.c dspt21.c dstt21.c 
   dstt22.c dsyt21.c dsyt22.c)

set(ZEIGTST  zchkee.c 
   zbdt01.c zbdt02.c zbdt03.c 
   zchkbb.c zchkbd.c zchkbk.c zchkbl.c zchkec.c 
   zchkgg.c zchkgk.c zchkgl.c zchkhb.c zchkhs.c zchkst.c 
   zckglm.c zckgqr.c zckgsv.c zcklse.c 
   zdrges.c zdrgev.c zdrgsx.c zdrgvx.c 
   zdrvbd.c zdrves.c zdrvev.c zdrvgg.c zdrvsg.c 
   zdrvst.c zdrvsx.c zdrvvx.c 
   zerrbd.c zerrec.c zerred.c zerrgg.c zerrhs.c zerrst.c 
   zget02.c zget10.c zget22.c zget23.c zget24.c 
   zget35.c zget36.c zget37.c zget38.c zget51.c zget52.c 
   zget54.c zglmts.c zgqrts.c zgrqts.c zgsvts.c 
   zhbt21.c zhet21.c zhet22.c zhpt21.c zhst01.c 
   zlarfy.c zlarhs.c zlatm4.c zlctes.c zlctsx.c zlsets.c zsbmv.c 
   zsgt01.c zslect.c 
   zstt21.c zstt22.c zunt01.c zunt03.c)

macro(add_eig_executable name )
  add_executable(${name} ${ARGN})
  target_link_libraries(${name} tmglib lapack )
endmacro(add_eig_executable)

add_eig_executable(xeigtsts ${SEIGTST} ${SCIGTST} ${AEIGTST} 
  ${SECOND_SRC} )

add_eig_executable(xeigtstc ${CEIGTST} ${SCIGTST} ${AEIGTST}
  ${SECOND_SRC} )

add_eig_executable(xeigtstd ${DEIGTST} ${DZIGTST} ${AEIGTST} 
  ${DSECOND_SRC} )

add_eig_executable(xeigtstz ${ZEIGTST} ${DZIGTST} ${AEIGTST}
  ${DSECOND_SRC} )
