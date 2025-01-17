if(NOT IS_MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif()

# prefer GLVND
set(OpenGL_GL_PREFERENCE GLVND)
find_package(OpenGL REQUIRED)
if (OPENGL_opengl_LIBRARY)
    set(FOUND_GL_LIBRARIES ${OPENGL_opengl_LIBRARY})
else ()
    set(FOUND_GL_LIBRARIES ${OPENGL_gl_LIBRARY})
endif ()
if (OPENGL_glx_LIBRARY)
    list(APPEND FOUND_GL_LIBRARIES ${OPENGL_glx_LIBRARY})
endif ()
message(STATUS "OpenGL libraries: ${FOUND_GL_LIBRARIES}")

add_definitions("-DSIFTGPU_NO_DEVIL")

set(SIFT_GPU_SOURCE_FILES
    FrameBufferObject.cpp
    FrameBufferObject.h
    GlobalUtil.cpp
    GlobalUtil.h
    GLTexImage.cpp
    GLTexImage.h
    ProgramGLSL.cpp
    ProgramGLSL.h
    ProgramGPU.h
    PyramidGL.cpp
    PyramidGL.h
    ShaderMan.cpp
    ShaderMan.h
    SiftGPU.cpp
    SiftGPU.h
    SiftMatch.cpp
    SiftMatch.h
    SiftPyramid.cpp
    SiftPyramid.h
)

if(CUDA_ENABLED)
    add_definitions("-DCUDA_SIFTGPU_ENABLED")

    set(SIFT_GPU_SOURCE_FILES
        ${SIFT_GPU_SOURCE_FILES}
        CuTexImage.cpp
        CuTexImage.h
        ProgramCU.cu
        ProgramCU.h
        PyramidCU.cpp
        PyramidCU.h
        SiftMatchCU.cpp
        SiftMatchCU.h
    )

    COLMAP_ADD_STATIC_CUDA_LIBRARY(3rd_sift_gpu ${SIFT_GPU_SOURCE_FILES})
else()
    add_library(3rd_sift_gpu STATIC ${SIFT_GPU_SOURCE_FILES})
endif()

target_link_libraries(3rd_sift_gpu
    ${SIFT_GPU_LIBRARIES}
    ${FOUND_GL_LIBRARIES}
)

set(GLEW_INCLUDE_DIR ${MVSTUDIO_THIRD_PARTY}/glew/include)
target_include_directories(3rd_sift_gpu PRIVATE ${GLEW_INCLUDE_DIR} ${OPENGL_INCLUDE_DIR})
target_link_libraries(3rd_sift_gpu 3rd_glew)
target_compile_definitions(3rd_sift_gpu PRIVATE GLEW_STATIC)

if (WIN32)
    target_link_libraries(3rd_sift_gpu Winmm.lib)
endif (WIN32)

set_target_properties(3rd_sift_gpu PROPERTIES FOLDER "3rd_party")