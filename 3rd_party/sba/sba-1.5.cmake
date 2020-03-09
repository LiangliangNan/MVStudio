#
# This file was genared by ./proj2cmake and will be overwritten on it's next run!
# Please put all configurations in the cmake_conf/*.cmake files.
#

SET(3rd_sba_SRC
    "sba_chkjac.c"
    "sba_crsm.c"
    "sba_lapack.c"
    "sba_levmar.c"
    "sba_levmar_wrap.c"
   )

SET(3rd_sba_DEPS
        blas
        f2c
        clapack
   )
