#
# This file was genared by ./proj2cmake and will be overwritten on it's next run!
# Please put all configurations in the cmake_conf/*.cmake files.
#

cmake_minimum_required(VERSION 2.10)


INCLUDE("pmvs.cmake")

ADD_LIBRARY(cmvs-pmvs STATIC
            ${3rd_cmvs-pmvs_SRC})

TARGET_LINK_LIBRARIES(cmvs-pmvs
            ${3rd_cmvs-pmvs_DEPS}
            ${3rd_cmvs-pmvs_ADDITIONAL_DEPS}
            ${SOLUTION_STATIC_LIB_DEPS}
            ${SOLUTION_GENERAL_DEPS})



set(CIMG_INCLUDE_DIR ${MVSTUDIO_THIRD_PARTY}/cimg)
set(METIS_INCLUDE_DIR ${MVSTUDIO_THIRD_PARTY}/graclus/metisLib)
set(JPEG_INCLUDE_DIR ${MVSTUDIO_THIRD_PARTY}/jpeg ${MVSTUDIO_BUILD_DIR}/3rd_party/jpeg)
set(MINIBOOST_INCLUDE_DIR ${MVSTUDIO_THIRD_PARTY}/miniBoost)
set(TINYCTHREAD_INCLUDE_DIR ${MVSTUDIO_THIRD_PARTY}/tinycthread)
set(NLOPT_INCLUDE_DIR ${MVSTUDIO_THIRD_PARTY}/nlopt/src/api ${MVSTUDIO_BUILD_DIR}/3rd_party/nlopt/src/api)

target_include_directories(cmvs-pmvs PRIVATE
        ${METIS_INCLUDE_DIR}
        ${CIMG_INCLUDE_DIR}
        ${JPEG_INCLUDE_DIR}
        ${MINIBOOST_INCLUDE_DIR}
        ${TINYCTHREAD_INCLUDE_DIR}
        ${NLOPT_INCLUDE_DIR}
        )

set_target_properties(cmvs-pmvs PROPERTIES FOLDER "3rd_party")