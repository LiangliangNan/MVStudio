#
# This file was genared by ./proj2cmake and will be overwritten on it's next run!
# Please put all configurations in the cmake_conf/*.cmake files.
#

cmake_minimum_required(VERSION 2.8)


INCLUDE("algo.cmake")

ADD_LIBRARY(algo SHARED
            ${algo_SRC})

TARGET_LINK_LIBRARIES(algo
            ${algo_DEPS}
            ${algo_ADDITIONAL_DEPS}
            ${SOLUTION_SHARED_LIB_DEPS}
            ${SOLUTION_GENERAL_DEPS})

set(MINIBOOST_INCLUDE_DIR ${MVSTUDIO_THIRD_PARTY}/miniBoost)
set(TINYCTHREAD_INCLUDE_DIR ${MVSTUDIO_THIRD_PARTY}/tinycthread)
target_include_directories(algo PRIVATE
        ${MINIBOOST_INCLUDE_DIR}
        ${TINYCTHREAD_INCLUDE_DIR}
        )