project(CLAPACK C)

# The f2c-generated sources use old-style function pointer declarations such as
# `typedef logical (*L_fp)();`. C23 interprets `()` as a zero-argument prototype,
# so recent GCC/Clang releases reject calls through those pointers. Keep this
# bundled third-party code on C11, where `()` still means unspecified arguments.
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS ON)

#enable_testing()
#include(CTest)

if(WIN32 AND NOT CYGWIN)
  set(SECOND_SRC  ${CLAPACK_SOURCE_DIR}/INSTALL/winsecond.c)
  set(DSECOND_SRC  ${CLAPACK_SOURCE_DIR}/INSTALL/windsecnd.c)
  add_definitions(-DNO_ISATTY -DMSDOS -DUSE_CLOCK)
else()
  set(SECOND_SRC  ${CLAPACK_SOURCE_DIR}/INSTALL/second.c)
  set(DSECOND_SRC  ${CLAPACK_SOURCE_DIR}/INSTALL/dsecnd.c)
endif()
enable_testing()
option(USE_BLAS_WRAP "pre-pend f2c_ to each function in blas" OFF)
if(NOT USE_BLAS_WRAP)
# _zrotg_ seems to be missing in the wrap header
  add_definitions(-DNO_BLAS_WRAP)
endif()
include_directories(${CLAPACK_SOURCE_DIR}/INCLUDE)
add_subdirectory(F2CLIBS)
add_subdirectory(BLAS)
add_subdirectory(SRC)
#add_subdirectory(TESTING)
set(CLAPACK_VERSION 3.2.1)
set(CPACK_PACKAGE_VERSION_MAJOR 3)
set(CPACK_PACKAGE_VERSION_MINOR 2)
set(CPACK_PACKAGE_VERSION_PATCH 1)
include(CPack)
export(TARGETS 3rd_f2c 3rd_blas 3rd_clapack FILE clapack-targets.cmake)
configure_file(${CLAPACK_SOURCE_DIR}/clapack-config-version.cmake.in
  ${CLAPACK_BINARY_DIR}/clapack-config-version.cmake @ONLY)
configure_file(${CLAPACK_SOURCE_DIR}/clapack-config.cmake.in
  ${CLAPACK_BINARY_DIR}/clapack-config.cmake @ONLY)
