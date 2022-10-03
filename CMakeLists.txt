cmake_minimum_required(VERSION 3.1)

if (APPLE)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "10.12" CACHE STRING "Minimum OS X deployment version" FORCE)
    message(STATUS "Minimum OS X deployment version: ${CMAKE_OSX_DEPLOYMENT_TARGET}")
endif ()

get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)


### Configuration
set(MVSTUDIO_ROOT         ${CMAKE_CURRENT_LIST_DIR})
set(MVSTUDIO_BUILD_DIR    ${CMAKE_CURRENT_BINARY_DIR})
set(MVSTUDIO_THIRD_PARTY  ${MVSTUDIO_ROOT}/3rd_party)

ADD_SUBDIRECTORY("3rd_party/clapack")
ADD_SUBDIRECTORY("3rd_party/cminpack")
ADD_SUBDIRECTORY("3rd_party/cmvs-pmvs")
ADD_SUBDIRECTORY("3rd_party/graclus")
ADD_SUBDIRECTORY("3rd_party/jpeg")
ADD_SUBDIRECTORY("3rd_party/nlopt")
ADD_SUBDIRECTORY("3rd_party/rply")
ADD_SUBDIRECTORY("3rd_party/png")
ADD_SUBDIRECTORY("3rd_party/QGLViewer")
ADD_SUBDIRECTORY("3rd_party/sba")
ADD_SUBDIRECTORY("3rd_party/tinycthread")
ADD_SUBDIRECTORY("3rd_party/zlib")
ADD_SUBDIRECTORY("3rd_party/glew")
ADD_SUBDIRECTORY("3rd_party/SiftGPU")

ADD_SUBDIRECTORY("libs/algo")
ADD_SUBDIRECTORY("libs/basic")
ADD_SUBDIRECTORY("libs/image")
ADD_SUBDIRECTORY("libs/math")
ADD_SUBDIRECTORY("libs/mvglib")
ADD_SUBDIRECTORY("libs/opengl")
ADD_SUBDIRECTORY("libs/pointset")
ADD_SUBDIRECTORY("libs/sfm")

ADD_SUBDIRECTORY("MVStudio")