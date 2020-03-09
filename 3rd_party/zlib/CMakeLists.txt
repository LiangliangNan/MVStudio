#
# This file was genared by ./proj2cmake and will be overwritten on it's next run!
# Please put all configurations in the cmake_conf/*.cmake files.
#

cmake_minimum_required(VERSION 2.8)


INCLUDE("3rd_zlib.cmake")

ADD_LIBRARY(zlib STATIC
            ${3rd_zlib_SRC})

TARGET_LINK_LIBRARIES(zlib
            ${3rd_zlib_DEPS}
            ${3rd_zlib_ADDITIONAL_DEPS}
            ${SOLUTION_STATIC_LIB_DEPS}
            ${SOLUTION_GENERAL_DEPS})

