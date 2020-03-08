cmake_minimum_required(VERSION 3.1)

get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})


set(${PROJECT_NAME}_HEADERS
        point_as_sphere.h
        point_set.h
        point_set_io.h
        point_set_render.h
        point_set_serializer_ply.h
        rply.h
        )

set(${PROJECT_NAME}_SOURCES
        point_as_sphere.cpp
        point_set.cpp
        point_set_io.cpp
        point_set_render.cpp
        point_set_serializer_ply.cpp
        rply.c
        )


add_library(${PROJECT_NAME} STATIC ${${PROJECT_NAME}_SOURCES} ${${PROJECT_NAME}_HEADERS} )

target_include_directories(${PROJECT_NAME} PRIVATE ${MVSTUDIO_INCLUDE_DIR})

target_link_libraries(${PROJECT_NAME}
        opengl
        basic
        math
        )

#set(MINIBOOST_INCLUDE_DIR ${MVSTUDIO_THIRD_PARTY}/miniBoost)
#set(TINYCTHREAD_INCLUDE_DIR ${MVSTUDIO_THIRD_PARTY}/tinycthread)
#target_include_directories(${PROJECT_NAME} PRIVATE
#        ${MINIBOOST_INCLUDE_DIR}
#        ${TINYCTHREAD_INCLUDE_DIR}
#        )
