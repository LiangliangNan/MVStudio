set(ALINTST  
   aladhd.c alaerh.c alaesm.c alahd.c alareq.c 
   alasum.c alasvm.c chkxer.c icopy.c ilaenv.c xlaenv.c xerbla.c)

set(SCLNTST slaord.c)

set(DZLNTST dlaord.c )

set(SLINTST  schkaa.c 
   schkeq.c schkgb.c schkge.c schkgt.c 
   schklq.c schkpb.c schkpo.c schkps.c schkpp.c 
   schkpt.c schkq3.c schkql.c schkqp.c schkqr.c schkrq.c 
   schksp.c schksy.c schktb.c schktp.c schktr.c 
   schktz.c 
   sdrvgt.c sdrvls.c sdrvpb.c 
   sdrvpp.c sdrvpt.c sdrvsp.c sdrvsy.c 
   serrgt.c serrlq.c serrls.c 
   serrpo.c serrps.c serrql.c serrqp.c serrqr.c 
   serrrq.c serrsy.c serrtr.c serrtz.c serrvx.c 
   sgbt01.c sgbt02.c sgbt05.c sgelqs.c sgeqls.c sgeqrs.c 
   sgerqs.c sget01.c sget02.c 
   sget03.c sget04.c sget06.c sget07.c sgtt01.c sgtt02.c 
   sgtt05.c slaptm.c slarhs.c slatb4.c slatb5.c slattb.c slattp.c 
   slattr.c slavsp.c slavsy.c slqt01.c slqt02.c 
   slqt03.c spbt01.c spbt02.c spbt05.c spot01.c 
   spot02.c spot03.c spot05.c spst01.c sppt01.c 
   sppt02.c sppt03.c sppt05.c sptt01.c sptt02.c 
   sptt05.c sqlt01.c sqlt02.c sqlt03.c sqpt01.c 
   sqrt01.c sqrt02.c sqrt03.c sqrt11.c sqrt12.c 
   sqrt13.c sqrt14.c sqrt15.c sqrt16.c sqrt17.c 
   srqt01.c srqt02.c srqt03.c srzt01.c srzt02.c 
   sspt01.c ssyt01.c 
   stbt02.c stbt03.c stbt05.c stbt06.c stpt01.c 
   stpt02.c stpt03.c stpt05.c stpt06.c strt01.c 
   strt02.c strt03.c strt05.c strt06.c 
   stzt01.c stzt02.c sgennd.c)

if(USEXBLAS)
  list(APPEND SLINTST sdrvgex.c serrgex.c sdrvgbx.c sdrvpox.c sebchvxx.c)
else()
  list(APPEND SLINTST sdrvge.c serrge.c sdrvgb.c sdrvpo.c)
endif()

set(CLINTST  cchkaa.c 
   cchkeq.c cchkgb.c cchkge.c cchkgt.c 
   cchkhe.c cchkhp.c cchklq.c cchkpb.c 
   cchkpo.c cchkps.c cchkpp.c cchkpt.c cchkq3.c cchkql.c cchkqp.c 
   cchkqr.c cchkrq.c cchksp.c cchksy.c cchktb.c 
   cchktp.c cchktr.c cchktz.c 
   cdrvgt.c cdrvhe.c cdrvhp.c 
   cdrvls.c cdrvpb.c cdrvpp.c cdrvpt.c 
   cdrvsp.c cdrvsy.c 
   cerrgt.c cerrhe.c cerrlq.c 
   cerrls.c cerrps.c cerrql.c cerrqp.c 
   cerrqr.c cerrrq.c cerrsy.c cerrtr.c cerrtz.c 
   cerrvx.c 
   cgbt01.c cgbt02.c cgbt05.c cgelqs.c cgeqls.c cgeqrs.c 
   cgerqs.c cget01.c cget02.c 
   cget03.c cget04.c cget07.c cgtt01.c cgtt02.c 
   cgtt05.c chet01.c chpt01.c claipd.c claptm.c clarhs.c clatb4.c clatb5.c 
   clatsp.c clatsy.c clattb.c clattp.c clattr.c 
   clavhe.c clavhp.c clavsp.c clavsy.c clqt01.c 
   clqt02.c clqt03.c cpbt01.c cpbt02.c cpbt05.c 
   cpot01.c cpot02.c cpot03.c cpot05.c cpst01.c 
   cppt01.c cppt02.c cppt03.c cppt05.c cptt01.c 
   cptt02.c cptt05.c cqlt01.c cqlt02.c cqlt03.c 
   cqpt01.c cqrt01.c cqrt02.c cqrt03.c cqrt11.c 
   cqrt12.c cqrt13.c cqrt14.c cqrt15.c cqrt16.c 
   cqrt17.c crqt01.c crqt02.c crqt03.c crzt01.c crzt02.c 
   csbmv.c  cspt01.c 
   cspt02.c cspt03.c csyt01.c csyt02.c csyt03.c 
   ctbt02.c ctbt03.c ctbt05.c ctbt06.c ctpt01.c 
   ctpt02.c ctpt03.c ctpt05.c ctpt06.c ctrt01.c 
   ctrt02.c ctrt03.c ctrt05.c ctrt06.c 
   ctzt01.c ctzt02.c sget06.c cgennd.c)

if(USEXBLAS)
  list(APPEND 
    CLINTST cdrvgex.c cdrvgbx.c cerrgex.c cdrvpox.c cerrpox.c cebchvxx.c)
else()
  list(APPEND CLINTST cdrvge.c cdrvgb.c cerrge.c cdrvpo.c cerrpo.c)
endif()

set(DLINTST  dchkaa.c 
   dchkeq.c dchkgb.c dchkge.c dchkgt.c 
   dchklq.c dchkpb.c dchkpo.c dchkps.c dchkpp.c 
   dchkpt.c dchkq3.c dchkql.c dchkqp.c dchkqr.c dchkrq.c 
   dchksp.c dchksy.c dchktb.c dchktp.c dchktr.c 
   dchktz.c 
   ddrvgt.c ddrvls.c ddrvpb.c 
   ddrvpp.c ddrvpt.c ddrvsp.c ddrvsy.c 
   derrgt.c derrlq.c derrls.c 
   derrps.c derrql.c derrqp.c derrqr.c 
   derrrq.c derrsy.c derrtr.c derrtz.c derrvx.c 
   dgbt01.c dgbt02.c dgbt05.c dgelqs.c dgeqls.c dgeqrs.c 
   dgerqs.c dget01.c dget02.c 
   dget03.c dget04.c dget06.c dget07.c dgtt01.c dgtt02.c 
   dgtt05.c dlaptm.c dlarhs.c dlatb4.c dlatb5.c dlattb.c dlattp.c 
   dlattr.c dlavsp.c dlavsy.c dlqt01.c dlqt02.c 
   dlqt03.c dpbt01.c dpbt02.c dpbt05.c dpot01.c 
   dpot02.c dpot03.c dpot05.c dpst01.c dppt01.c 
   dppt02.c dppt03.c dppt05.c dptt01.c dptt02.c 
   dptt05.c dqlt01.c dqlt02.c dqlt03.c dqpt01.c 
   dqrt01.c dqrt02.c dqrt03.c dqrt11.c dqrt12.c 
   dqrt13.c dqrt14.c dqrt15.c dqrt16.c dqrt17.c 
   drqt01.c drqt02.c drqt03.c drzt01.c drzt02.c 
   dspt01.c dsyt01.c 
   dtbt02.c dtbt03.c dtbt05.c dtbt06.c dtpt01.c 
   dtpt02.c dtpt03.c dtpt05.c dtpt06.c dtrt01.c 
   dtrt02.c dtrt03.c dtrt05.c dtrt06.c 
   dtzt01.c dtzt02.c dgennd.c)

if(USEXBLAS)
  list(APPEND 
    DLINTST  ddrvgex.c ddrvgbx.c derrgex.c ddrvpox.c derrpox.c debchvxx.c)
else()
  list(APPEND
    DLINTST  ddrvge.c ddrvgb.c derrge.c ddrvpo.c derrpo.c)
endif()

set(ZLINTST  zchkaa.c 
   zchkeq.c zchkgb.c zchkge.c zchkgt.c 
   zchkhe.c zchkhp.c zchklq.c zchkpb.c 
   zchkpo.c zchkps.c zchkpp.c zchkpt.c zchkq3.c zchkql.c zchkqp.c 
   zchkqr.c zchkrq.c zchksp.c zchksy.c zchktb.c 
   zchktp.c zchktr.c zchktz.c 
   zdrvgt.c zdrvhe.c zdrvhp.c 
   zdrvls.c zdrvpb.c zdrvpp.c zdrvpt.c 
   zdrvsp.c zdrvsy.c 
   zerrgt.c zerrhe.c zerrlq.c 
   zerrls.c zerrps.c zerrql.c zerrqp.c 
   zerrqr.c zerrrq.c zerrsy.c zerrtr.c zerrtz.c 
   zerrvx.c 
   zgbt01.c zgbt02.c zgbt05.c zgelqs.c zgeqls.c zgeqrs.c 
   zgerqs.c zget01.c zget02.c 
   zget03.c zget04.c zget07.c zgtt01.c zgtt02.c 
   zgtt05.c zhet01.c zhpt01.c zlaipd.c zlaptm.c zlarhs.c zlatb4.c zlatb5.c 
   zlatsp.c zlatsy.c zlattb.c zlattp.c zlattr.c 
   zlavhe.c zlavhp.c zlavsp.c zlavsy.c zlqt01.c 
   zlqt02.c zlqt03.c zpbt01.c zpbt02.c zpbt05.c 
   zpot01.c zpot02.c zpot03.c zpot05.c zpst01.c 
   zppt01.c zppt02.c zppt03.c zppt05.c zptt01.c 
   zptt02.c zptt05.c zqlt01.c zqlt02.c zqlt03.c 
   zqpt01.c zqrt01.c zqrt02.c zqrt03.c zqrt11.c 
   zqrt12.c zqrt13.c zqrt14.c zqrt15.c zqrt16.c 
   zqrt17.c zrqt01.c zrqt02.c zrqt03.c zrzt01.c zrzt02.c 
   zsbmv.c  zspt01.c 
   zspt02.c zspt03.c zsyt01.c zsyt02.c zsyt03.c 
   ztbt02.c ztbt03.c ztbt05.c ztbt06.c ztpt01.c 
   ztpt02.c ztpt03.c ztpt05.c ztpt06.c ztrt01.c 
   ztrt02.c ztrt03.c ztrt05.c ztrt06.c 
   ztzt01.c ztzt02.c dget06.c zgennd.c)

if(USEXBLAS)
  list(APPEND
    ZLINTST  zdrvgex.c zdrvgbx.c zerrgex.c zdrvpox.c zerrpox.c zebchvxx.c)
else()
  list(APPEND
    ZLINTST  zdrvge.c zdrvgb.c zerrge.c zdrvpo.c zerrpo.c)
endif()

set(DSLINTST  dchkab.c 
	ddrvab.c ddrvac.c derrab.c derrac.c dget08.c          
	alaerh.c alahd.c  aladhd.c alareq.c 
	chkxer.c dlarhs.c dlatb4.c xerbla.c 
	dget02.c dpot06.c)

set(ZCLINTST  zchkab.c 
	zdrvab.c zdrvac.c zerrab.c zerrac.c zget08.c          
	alaerh.c alahd.c  aladhd.c alareq.c 
	chkxer.c zget02.c zlarhs.c zlatb4.c 
	zsbmv.c  xerbla.c zpot06.c zlaipd.c)

set(SLINTSTRFP  schkrfp.c sdrvrfp.c sdrvrf1.c sdrvrf2.c sdrvrf3.c sdrvrf4.c serrrfp.c 
	slatb4.c slarhs.c sget04.c spot01.c spot03.c spot02.c 
	chkxer.c xerbla.c alaerh.c aladhd.c alahd.c alasvm.c )

set(DLINTSTRFP  dchkrfp.c ddrvrfp.c ddrvrf1.c ddrvrf2.c ddrvrf3.c ddrvrf4.c derrrfp.c 
	dlatb4.c dlarhs.c dget04.c dpot01.c dpot03.c dpot02.c 
	chkxer.c xerbla.c alaerh.c aladhd.c alahd.c alasvm.c )

set(CLINTSTRFP  cchkrfp.c cdrvrfp.c cdrvrf1.c cdrvrf2.c cdrvrf3.c cdrvrf4.c cerrrfp.c 
  claipd.c clatb4.c clarhs.c csbmv.c cget04.c cpot01.c cpot03.c cpot02.c 
  chkxer.c xerbla.c alaerh.c aladhd.c alahd.c alasvm.c )

set(ZLINTSTRFP  zchkrfp.c zdrvrfp.c zdrvrf1.c zdrvrf2.c zdrvrf3.c zdrvrf4.c zerrrfp.c 
  zlatb4.c zlaipd.c zlarhs.c zsbmv.c zget04.c zpot01.c zpot03.c zpot02.c 
  chkxer.c xerbla.c alaerh.c aladhd.c alahd.c alasvm.c )

macro(add_lin_executable name )
  add_executable(${name} ${ARGN})
  target_link_libraries(${name} tmglib lapack)
endmacro(add_lin_executable)

add_lin_executable(xlintsts ${ALINTST} ${SCLNTST} ${SLINTST} 
  ${SECOND_SRC} )

add_lin_executable(xlintstc  ${ALINTST} ${CLINTST} ${SCLNTST}
  ${SECOND_SRC} )
 
add_lin_executable(xlintstd  ${ALINTST} ${DLINTST} ${DZLNTST}
  ${DSECOND_SRC})
add_lin_executable(xlintstz  ${ALINTST} ${ZLINTST} ${DZLNTST}
  ${DSECOND_SRC})
 
add_lin_executable(xlintstds ${DSLINTST}
		${SECOND_SRC} 
		${DSECOND_SRC} )
add_lin_executable(xlintstzc ${ZCLINTST}
		${SECOND_SRC} 
		${DSECOND_SRC} )
 
add_lin_executable(xlintstrfs ${SLINTSTRFP}
		${SECOND_SRC})
 
add_lin_executable(xlintstrfd ${DLINTSTRFP}
  ${DSECOND_SRC})
 
add_lin_executable(xlintstrfc  ${CLINTSTRFP}
  ${SECOND_SRC})
add_lin_executable(xlintstrfz ${ZLINTSTRFP}
  ${DSECOND_SRC})

