#
# This file was genared by ./proj2cmake and will be overwritten on it's next run!
# Please put all configurations in the cmake_conf/*.cmake files.
#

cmake_minimum_required(VERSION 2.8)

INCLUDE("mvglib.cmake")

ADD_LIBRARY(mvglib STATIC
            ${mvglib_SRC})

TARGET_LINK_LIBRARIES(mvglib
            ${mvglib_DEPS}
            ${mvglib_ADDITIONAL_DEPS}
            ${SOLUTION_STATIC_LIB_DEPS}
            ${SOLUTION_GENERAL_DEPS})


set(CIMG_INCLUDE_DIR ${MVSTUDIO_THIRD_PARTY}/cimg)
target_include_directories(mvglib PRIVATE ${CIMG_INCLUDE_DIR} ${MVSTUDIO_INCLUDE_DIR})
