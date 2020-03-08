cmake_minimum_required(VERSION 3.1)

get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})


set(${PROJECT_NAME}_HEADERS
        image_io.h
        ply_reader_writer.h
        point_cloud_io.h
        resources.h
        )

set(${PROJECT_NAME}_SOURCES
        image_io.cpp
        ply_reader_writer.cpp
        point_cloud_io.cpp
        point_cloud_io_ply.cpp
        resources.cpp
        )


add_library(${PROJECT_NAME} STATIC ${${PROJECT_NAME}_SOURCES} ${${PROJECT_NAME}_HEADERS})

set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "easy3d")

target_include_directories(${PROJECT_NAME} PRIVATE ${MVSTUDIO_INCLUDE_DIR})

set(LASTOOLS_INCLUDE_DIR ${EASY3D_THIRD_PARTY}/LAStools/LASzip/src ${EASY3D_THIRD_PARTY}/LAStools/LASlib/inc)
target_include_directories(${PROJECT_NAME} PRIVATE ${LASTOOLS_INCLUDE_DIR})

target_link_libraries(${PROJECT_NAME} core util rply)

if (MSVC)
    target_compile_definitions(${PROJECT_NAME} PRIVATE _CRT_SECURE_NO_DEPRECATE)
endif ()
