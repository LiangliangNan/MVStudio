cmake_minimum_required(VERSION 3.1)

get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})


set(${PROJECT_NAME}_HEADERS
        opengl_info.h
        quick_text.h
        shader.h
        shader_factory.h
        shader_program.h
        )

set(${PROJECT_NAME}_SOURCES
        opengl_info.cpp
        quick_text.cpp
        shader.cpp
        shader_factory.cpp
        shader_program.cpp
        )


add_library(${PROJECT_NAME} STATIC ${${PROJECT_NAME}_SOURCES} ${${PROJECT_NAME}_HEADERS} )

set(GLEW_INCLUDE_DIR ${MVSTUDIO_THIRD_PARTY}/glew/include)
target_include_directories(${PROJECT_NAME} PRIVATE ${MVSTUDIO_INCLUDE_DIR} ${GLEW_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} glew)

target_compile_definitions(${PROJECT_NAME} PRIVATE GLEW_STATIC)

if (MSVC)
    target_compile_definitions(${PROJECT_NAME} PRIVATE _CRT_SECURE_NO_DEPRECATE)
endif ()