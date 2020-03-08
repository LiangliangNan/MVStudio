cmake_minimum_required(VERSION 3.1)

get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})


set(${PROJECT_NAME}_HEADERS
        arrays.h
        assertions.h
        basic_types.h
        binary_stream.h
        canvas.h
        counted.h
        file_utils.h
        line_stream.h
        logger.h
        object.h
        progress.h
        real_timer.h
        smart_pointer.h
        stop_watch.h
        text_utils.h
        timer.h
        )

set(${PROJECT_NAME}_SOURCES
        assertions.cpp
        basic_types.cpp
        binary_stream.cpp
        counted.cpp
        file_utils.cpp
        logger.cpp
        object.cpp
        progress.cpp
        real_timer.cpp
        stop_watch.cpp
        text_utils.cpp
        timer.cpp
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