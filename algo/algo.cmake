#
# This file was genared by ./proj2cmake and will be overwritten on it's next run!
# Please put all configurations in the cmake_conf/*.cmake files.
#

SET(algo_SRC
        "dense_reconstruction.cpp"
        "image_matching.cpp"
        "project.cpp"
        "sparse_reconstruction.cpp"
        )

SET(algo_DEPS
        pointset
        clapack
        opengl
        png
        jpeg
        image
        nlopt
        sift_gpu
        blas
        cmvs-pmvs
        basic
        f2c
        sfm
        tinycthread
        )
