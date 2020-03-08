#
# This file was genared by ./proj2cmake and will be overwritten on it's next run!
# Please put all configurations in the cmake_conf/*.cmake files.
#

SET(sfm_SRC
        "bundle_adjustment.cpp"
        "camera.cpp"
        "epipolar.cpp"
        "geometry.cpp"
        "image_data.cpp"
        "keys.cpp"
        "match_table.cpp"
        "register.cpp"
        "sfm.cpp"
        "sfm_bundle_adjustment.cpp"
        "sfm_file_management.cpp"
        "sfm_geometric_constraints.cpp"
        "sfm_matches.cpp"
        "sfm_option.cpp"
        "sfm_tracks.cpp"
        "sfm_util.cpp"
        )

SET(sfm_DEPS
        pointset
        clapack
        image
        mvglib
        blas
        basic
        f2c
        math
        3rd_sba
        )
