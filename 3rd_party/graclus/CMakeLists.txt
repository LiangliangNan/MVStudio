#
# This file was genared by ./proj2cmake and will be overwritten on it's next run!
# Please put all configurations in the cmake_conf/*.cmake files.
#

cmake_minimum_required(VERSION 2.10)


INCLUDE("graclus_1.2.cmake")

ADD_LIBRARY(graclus STATIC
            ${3rd_graclus_SRC})

TARGET_LINK_LIBRARIES(graclus
            ${3rd_graclus_DEPS}
            ${3rd_graclus_ADDITIONAL_DEPS}
            ${SOLUTION_STATIC_LIB_DEPS}
            ${SOLUTION_GENERAL_DEPS})

set(METIS_INCLUDE_DIR ${MVSTUDIO_THIRD_PARTY}/graclus/metisLib)
target_include_directories(graclus PRIVATE ${METIS_INCLUDE_DIR})

set_target_properties(graclus PROPERTIES FOLDER "3rd_party")