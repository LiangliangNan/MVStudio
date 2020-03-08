#
# This file was genared by ./proj2cmake and will be overwritten on it's next run!
# Please put all configurations in the cmake_conf/*.cmake files.
#

cmake_minimum_required(VERSION 2.8)


INCLUDE("sfm.cmake")

ADD_LIBRARY(sfm SHARED
            ${sfm_SRC})

TARGET_LINK_LIBRARIES(sfm
            ${sfm_DEPS}
            ${sfm_ADDITIONAL_DEPS}
            ${SOLUTION_SHARED_LIB_DEPS}
            ${SOLUTION_GENERAL_DEPS})

