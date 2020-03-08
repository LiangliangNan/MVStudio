cmake_minimum_required(VERSION 3.1)

get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})


set(${PROJECT_NAME}_HEADERS
        5point.h
        fmatrix.h
        homography.h
        horn.h
        poly1.h
        poly3.h
        qsort.h
        triangulate.h
        )

set(${PROJECT_NAME}_SOURCES
        5point.cpp
        fmatrix.cpp
        homography.cpp
        horn.cpp
        poly1.cpp
        poly3.cpp
        qsort.cpp
        triangulate.cpp
        )


add_library(${PROJECT_NAME} STATIC ${${PROJECT_NAME}_SOURCES} ${${PROJECT_NAME}_HEADERS} )


set(CIMG_INCLUDE_DIR ${MVSTUDIO_THIRD_PARTY}/cimg)
target_include_directories(mvglib PRIVATE ${CIMG_INCLUDE_DIR})

