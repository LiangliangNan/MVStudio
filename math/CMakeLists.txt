cmake_minimum_required(VERSION 3.1)

get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})


set(${PROJECT_NAME}_HEADERS
        factor.h
        matrix.h
        matrix_driver.h
        svd.h
        )

set(${PROJECT_NAME}_SOURCES
        factor.cpp
        matrix_driver.c
        svd.cpp
        )



add_library(${PROJECT_NAME} STATIC ${${PROJECT_NAME}_SOURCES} ${${PROJECT_NAME}_HEADERS} )

target_link_libraries(${PROJECT_NAME} cminpack)

target_include_directories(${PROJECT_NAME} PRIVATE ${MVSTUDIO_INCLUDE_DIR})
