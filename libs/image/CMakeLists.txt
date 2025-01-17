get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})


set(${PROJECT_NAME}_HEADERS
        color.h
        image.h
        image_io.h
        image_serializer.h
        image_serializer_bmp.h
        image_serializer_jpeg.h
        image_serializer_png.h
        image_store.h
        )

set(${PROJECT_NAME}_SOURCES
        image.cpp
        image_io.cpp
        image_serializer.cpp
        image_serializer_bmp.cpp
        image_serializer_jpeg.cpp
        image_serializer_png.cpp
        image_store.cpp
        )


add_library(${PROJECT_NAME} STATIC ${${PROJECT_NAME}_SOURCES} ${${PROJECT_NAME}_HEADERS} )

target_include_directories(${PROJECT_NAME} PRIVATE ${MVSTUDIO_INCLUDE_DIR})

target_link_libraries(${PROJECT_NAME}
        3rd_jpeg
        3rd_png
        3rd_zlib
        basic
        )

set(PNG_INCLUDE_DIR ${MVSTUDIO_THIRD_PARTY}/png)
set(JPEG_INCLUDE_DIR ${MVSTUDIO_THIRD_PARTY}/jpeg ${MVSTUDIO_BUILD_DIR}/3rd_party/jpeg)
target_include_directories(image PRIVATE ${PNG_INCLUDE_DIR} ${JPEG_INCLUDE_DIR})

set_target_properties(image PROPERTIES FOLDER "libs")