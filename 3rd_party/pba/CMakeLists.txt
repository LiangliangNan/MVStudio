#
# This file was genared by ./proj2cmake and will be overwritten on it's next run!
# Please put all configurations in the cmake_conf/*.cmake files.
#

cmake_minimum_required(VERSION 2.8)


INCLUDE("pba.cmake")

ADD_LIBRARY(3rd_pba STATIC
            ${3rd_pba_SRC})

TARGET_LINK_LIBRARIES(3rd_pba
            ${3rd_pba_DEPS}
            ${3rd_pba_ADDITIONAL_DEPS}
            ${SOLUTION_STATIC_LIB_DEPS}
            ${SOLUTION_GENERAL_DEPS})

