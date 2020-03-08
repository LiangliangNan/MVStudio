#
# This file was genared by ./proj2cmake and will be overwritten on it's next run!
# Please put all configurations in the cmake_conf/*.cmake files.
#

SET(image_SRC
    "image.cpp"
    "image_io.cpp"
    "image_serializer.cpp"
    "image_serializer_bmp.cpp"
    "image_serializer_jpeg.cpp"
    "image_serializer_png.cpp"
    "image_store.cpp"
   )

SET(image_DEPS
    png
    jpeg
    basic
    zlib
   )
