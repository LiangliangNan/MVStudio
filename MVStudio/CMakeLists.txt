get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})

include(../cmake/UseQt5.cmake)
if (NOT Qt5_FOUND)
    message(WARNING "${PROJECT_NAME} requires Qt (Qt6 or Qt5), but Qt was not found. You can set 'Qt*_DIR' to the "
            "directory containing 'Qt*Config.cmake' or 'qt*-config.cmake' (* is 5 or 6) to specify the location of Qt.")
    return()
endif ()


set(${PROJECT_NAME}_HEADERS
        main_window.h
        paint_canvas.h
        widget_image_list.h
        image_item.h
        )

set(${PROJECT_NAME}_SOURCES
        main_window.cpp
        main.cpp
        paint_canvas.cpp
        widget_image_list.cpp
        )

set(${PROJECT_NAME}_FORMS
        main_window.ui
        )

set(${PROJECT_NAME}_RESOURCES
        Resources/MVStudio.qrc
        Resources/MVStudio.rc
        )

if (APPLE)
    # This tells CMake where to find and install the ICON file
    set(${PROJECT_NAME}_ICON ${CMAKE_CURRENT_SOURCE_DIR}/Resources/MVStudio.icns)

    #Append the new files to the resources file list
    list(APPEND ${PROJECT_NAME}_RESOURCES ${${PROJECT_NAME}_ICON})

    add_executable(${PROJECT_NAME} MACOSX_BUNDLE ${${PROJECT_NAME}_ICON} ${${PROJECT_NAME}_HEADERS} ${${PROJECT_NAME}_SOURCES} ${${PROJECT_NAME}_FORMS} ${${PROJECT_NAME}_RESOURCES})

    set_source_files_properties(${${PROJECT_NAME}_ICON} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")

    ## Note: On macOS, high-DPI support is enabled by settings in the Info.plist file. Make sure they are present.
    ## For more information, see: http://doc.qt.io/qt-5/highdpi.html
    ## Enable high DPI support by using my Info.plist file in which I added NSHighResolutionCapable to be true.
    set_target_properties(${PROJECT_NAME} PROPERTIES
            INSTALL_RPATH "@executable_path/../Frameworks"
            MACOSX_BUNDLE_BUNDLE_NAME ${PROJECT_NAME}
            MACOSX_BUNDLE TRUE
            MACOSX_BUNDLE_INFO_PLIST "${CMAKE_CURRENT_SOURCE_DIR}/Mac/MVStudio.plist"
            MACOSX_BUNDLE_ICON_FILE "MVStudio.icns"
            MACOSX_BUNDLE_GUI_IDENTIFIER "MVStudio"
            )
else ()
    add_executable(${PROJECT_NAME} ${${PROJECT_NAME}_HEADERS} ${${PROJECT_NAME}_SOURCES} ${${PROJECT_NAME}_FORMS} ${${PROJECT_NAME}_RESOURCES})
endif ()


set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "MVStudio")


# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
target_compile_definitions(${PROJECT_NAME} PRIVATE QT_DEPRECATED_WARNINGS)


target_link_libraries(${PROJECT_NAME} ${QtLibs})

target_link_libraries(${PROJECT_NAME} basic math pointset algo 3rd_QGLViewer 3rd_sift_gpu)

set(GLEW_INCLUDE_DIR ${MVSTUDIO_THIRD_PARTY}/glew/include)
target_include_directories(${PROJECT_NAME} PRIVATE ${GLEW_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} 3rd_glew)
target_compile_definitions(${PROJECT_NAME} PRIVATE GLEW_STATIC)