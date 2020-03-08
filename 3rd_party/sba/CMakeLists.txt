#
# This file was genared by ./proj2cmake and will be overwritten on it's next run!
# Please put all configurations in the cmake_conf/*.cmake files.
#

cmake_minimum_required(VERSION 2.8)


INCLUDE("sba-1.5.cmake")

ADD_LIBRARY(sba STATIC
            ${3rd_sba_SRC})

TARGET_LINK_LIBRARIES(sba
            ${3rd_sba_DEPS}
            ${3rd_sba_ADDITIONAL_DEPS}
            ${SOLUTION_STATIC_LIB_DEPS}
            ${SOLUTION_GENERAL_DEPS})

