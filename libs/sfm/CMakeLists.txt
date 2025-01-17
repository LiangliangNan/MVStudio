get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})


set(${PROJECT_NAME}_HEADERS
        bundle_adjustment.h
        camera.h
        epipolar.h
        geometry.h
        image_data.h
        keys.h
        match_table.h
        register.h
        sfm.h
        sfm_option.h
        sfm_util.h
        )

set(${PROJECT_NAME}_SOURCES
        bundle_adjustment.cpp
        camera.cpp
        epipolar.cpp
        geometry.cpp
        image_data.cpp
        keys.cpp
        match_table.cpp
        register.cpp
        sfm.cpp
        sfm_bundle_adjustment.cpp
        sfm_file_management.cpp
        sfm_geometric_constraints.cpp
        sfm_matches.cpp
        sfm_option.cpp
        sfm_tracks.cpp
        sfm_util.cpp
        )


add_library(${PROJECT_NAME} STATIC ${${PROJECT_NAME}_SOURCES} ${${PROJECT_NAME}_HEADERS} )

target_include_directories(${PROJECT_NAME} PRIVATE ${MVSTUDIO_INCLUDE_DIR})

target_link_libraries(${PROJECT_NAME}
        3rd_clapack
        3rd_blas
        3rd_sba
        3rd_f2c
        basic
        math
        image
        mvglib
        pointset
)

set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "libs")