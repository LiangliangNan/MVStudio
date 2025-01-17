# Define the source files for the graclus library
file(GLOB 3RD_GRACLUS_SRC
        "metisLib/*.c"
        "multilevelLib/*.c"
)

# Add the graclus library
add_library(3rd_graclus STATIC ${3RD_GRACLUS_SRC})

# Set the include directory for the graclus library
set(METIS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/metisLib)
target_include_directories(3rd_graclus PRIVATE ${METIS_INCLUDE_DIR})

# Organize the target in the "3rd_party" folder (for IDEs like Visual Studio)
set_target_properties(3rd_graclus PROPERTIES FOLDER "3rd_party")