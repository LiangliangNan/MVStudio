cmake_minimum_required(VERSION 3.1)

get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})


set(${PROJECT_NAME}_HEADERS
        dense_reconstruction.h
        image_matching.h
        project.h
        sparse_reconstruction.h
        )

set(${PROJECT_NAME}_SOURCES
        dense_reconstruction.cpp
        image_matching.cpp
        project.cpp
        sparse_reconstruction.cpp
        )


add_library(${PROJECT_NAME} STATIC ${${PROJECT_NAME}_SOURCES} ${${PROJECT_NAME}_HEADERS} )

target_include_directories(${PROJECT_NAME} PRIVATE ${MVSTUDIO_INCLUDE_DIR})

target_link_libraries(${PROJECT_NAME}
        clapack
        nlopt
        sift_gpu
        blas
        cmvs-pmvs
        f2c
        tinycthread
        image
        sfm
        mvglib
        )

set(MINIBOOST_INCLUDE_DIR ${MVSTUDIO_THIRD_PARTY}/miniBoost)
set(TINYCTHREAD_INCLUDE_DIR ${MVSTUDIO_THIRD_PARTY}/tinycthread)
target_include_directories(${PROJECT_NAME} PRIVATE
        ${MINIBOOST_INCLUDE_DIR}
        ${TINYCTHREAD_INCLUDE_DIR}
        )

set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "libs")