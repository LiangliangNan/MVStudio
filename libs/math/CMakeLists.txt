cmake_minimum_required(VERSION 3.1)

get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})


set(${PROJECT_NAME}_HEADERS
        complex.h
        factor.h
        line.h
        math_types.h
        matrix.h
        matrix_driver.h
        oriented_line.h
        plane.h
        polynomial.h
        quaternion.h
        segment.h
        svd.h
        vecg.h
        vecgd.h
        vectord.h
        )

set(${PROJECT_NAME}_SOURCES
        factor.cpp
        math_types.cpp
        matrix_driver.c
        quaternion.cpp
        svd.cpp
        vectord.cpp
        )


add_library(${PROJECT_NAME} STATIC ${${PROJECT_NAME}_SOURCES} ${${PROJECT_NAME}_HEADERS} )

target_link_libraries(${PROJECT_NAME} basic cminpack f2c blas clapack)

set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "libs")