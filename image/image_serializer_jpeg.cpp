
#include "image_serializer_jpeg.h"
#include "image.h"
#include <easy3d/util/logging.h>

#include <stdio.h>
extern "C" {
#include "jpeglib.h"
}

#include <setjmp.h>


using namespace easy3d;

struct my_error_mgr {
	struct jpeg_error_mgr pub;	/* "public" fields */
	jmp_buf setjmp_buffer;	/* for return to caller */
};

typedef struct my_error_mgr * my_error_ptr;

/*
* Here's the routine that will replace the standard error_exit method:
*/
METHODDEF(void) my_error_exit(j_common_ptr cinfo) {
	static char buffer[1024];
	(*cinfo->err->format_message)(cinfo, buffer);
	LOG(ERROR) << buffer << std::endl;
	/* cinfo->err really points to a my_error_mgr struct, so coerce pointer */
	my_error_ptr myerr = (my_error_ptr)cinfo->err;
	/* Always display the message. */
	/* We could postpone this until after returning, if we chose. */
	(*cinfo->err->output_message) (cinfo);
	/* Return control to the setjmp point */
	longjmp(myerr->setjmp_buffer, 1);
}


Image* ImageSerializer_jpeg::serialize_read(const std::string& filename_in) {
	//Logger::out("ImageSerializer_jpeg") << "Loading " << filename_in << std::endl ;
	Image* result = nullptr;
	const char* filename = filename_in.c_str();
	/* We use our private extension JPEG error handler.
	* Note that this struct must live as long as the main JPEG parameter
	* struct, to avoid dangling-pointer problems.
	*/
	struct my_error_mgr jerr;
	/* This struct contains the JPEG decompression parameters and pointers to
	* working space (which is allocated as needed by the JPEG library).
	*/
	struct jpeg_decompress_struct cinfo;

	FILE * infile = nullptr;		/* source file */
	JSAMPARRAY buffer;		/* Output row buffer */
	int row_stride;		/* physical row width in output buffer */

	if ((infile = fopen(filename, "rb")) == NULL) {
		LOG(ERROR)
			<< "could not open file:" << filename_in << std::endl;
		return nullptr;
	}

	/* Step 1: allocate and initialize JPEG decompression object */

	/* We set up the normal JPEG error routines, then override error_exit. */
	cinfo.err = jpeg_std_error(&jerr.pub);
	jerr.pub.error_exit = my_error_exit;
	/* Establish the setjmp return context for my_error_exit to use. */
	if (setjmp(jerr.setjmp_buffer)) {
		/* If we get here, the JPEG code has signaled an error.
		* We need to clean up the JPEG object, close the input file, and return.
		*/
		jpeg_destroy_decompress(&cinfo);
		fclose(infile);
		delete result;
		return nullptr;
	}

	/* Now we can initialize the JPEG decompression object. */
	jpeg_create_decompress(&cinfo);

	/* Step 2: specify data source (eg, a file) */
	jpeg_stdio_src(&cinfo, infile);

	/* Step 3: read file parameters with jpeg_read_header() */
	(void)jpeg_read_header(&cinfo, TRUE);
	/* We can ignore the return value from jpeg_read_header since
	*   (a) suspension is not possible with the stdio data source, and
	*   (b) we passed TRUE to reject a tables-only JPEG file as an error.
	* See libjpeg.doc for more info.
	*/

	/* Step 4: set parameters for decompression */

	/* In this example, we don't need to change any of the defaults set by
	* jpeg_read_header(), so we do nothing here.
	*/

	/* Step 5: Start decompressor */

	(void)jpeg_start_decompress(&cinfo);
	/* We can ignore the return value since suspension is not possible
	* with the stdio data source.
	*/

	/* We may need to do some setup of our own at this point before reading
	* the data.  After jpeg_start_decompress() we have the correct scaled
	* output image dimensions available, as well as the output colormap
	* if we asked for color quantization.
	* In this example, we need to make an output work buffer of the right size.
	*/
	/* JSAMPLEs per row in output buffer */
	row_stride = cinfo.output_width * cinfo.output_components;

	// 	Logger::out("ImageSerializer_jpeg") << "JPEG image: " << cinfo.output_width << " x " << cinfo.output_height
	// 		<< " - " << cinfo.output_components << std::endl ;

	switch (cinfo.output_components) {
	case 3:
		result = new Image(Image::RGB, cinfo.output_width, cinfo.output_height);
		break;
	default:
		return nullptr;
		break;
	}

	/* Make a one-row-high sample array that will go away when done with image */
	buffer = (*cinfo.mem->alloc_sarray)
		((j_common_ptr)&cinfo, JPOOL_IMAGE, row_stride, 1);

	/* Step 6: while (scan lines remain to be read) */
	/*           jpeg_read_scanlines(...); */

	/* Here we use the library's state variable cinfo.output_scanline as the
	* loop counter, so that we don't have to keep track ourselves.
	*/

	while (cinfo.output_scanline < cinfo.output_height) {
		/* jpeg_read_scanlines expects an array of pointers to scanlines.
		* Here the array is only one element long, but you could ask for
		* more than one scanline at a time if that's more convenient.
		*/
		(void)jpeg_read_scanlines(&cinfo, buffer, 1);
		// Note: current line is cinfo.output_scanline-1 since jpeg_read_scanlines
		// just incremented it !!
		Memory::copy(result->pixel_base(0, cinfo.output_height - cinfo.output_scanline), buffer[0], row_stride);
	}

	/* Step 7: Finish decompression */

	(void)jpeg_finish_decompress(&cinfo);
	/* We can ignore the return value since suspension is not possible
	* with the stdio data source.
	*/

	/* Step 8: Release JPEG decompression object */

	/* This is an important step since it will release a good deal of memory. */
	jpeg_destroy_decompress(&cinfo);

	/* After finish_decompress, we can close the input file.
	* Here we postpone it until after no more JPEG errors are possible,
	* so as to simplify the setjmp error logic above.  (Actually, I don't
	* think that jpeg_destroy can do an error exit, but why assume anything...)
	*/
	fclose(infile);

	/* At this point you may want to check to see whether any corrupt-data
	* warnings occurred (test whether jerr.pub.num_warnings is nonzero).
	*/
	
	return result;
}


bool ImageSerializer_jpeg::serialize_write(
	const std::string& file_name, const Image* image
	)
{
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;

	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);


	FILE * outfile = fopen(file_name.c_str(), "wb");
	if (outfile  == NULL) {
		LOG(ERROR) << "could not open file:" << file_name << std::endl;
		return false;
	}

	jpeg_stdio_dest(&cinfo, outfile);

	int w = image->width();
	int h = image->height();
	cinfo.image_width = w;     /* image width and height, in pixels */
	cinfo.image_height = h;
	cinfo.input_components = 3;     /* # of color components per pixel */
	cinfo.in_color_space = JCS_RGB; /* colorspace of input image */
	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo, 98, TRUE);

	jpeg_start_compress(&cinfo, TRUE);

	JSAMPROW row = new JSAMPLE[3 * w];
	for (int y = 0; y < h; y++) {
		// JSAMPROW row_pointer[1];        /* pointer to a single row */
		int row_stride;                 /* physical row width in buffer */
		row_stride = w * 3;        /* JSAMPLEs per row in image_buffer */

		for (int x = 0; x < w; x++) {
			Memory::pointer c = const_cast<Image*>(image)->pixel_base(x, h - y - 1); 
			row[3 * x + 0] = c[0];
			row[3 * x + 1] = c[1];
			row[3 * x + 2] = c[2];
		}

		jpeg_write_scanlines(&cinfo, &row, 1);
	}

	jpeg_finish_compress(&cinfo);
	jpeg_destroy_compress(&cinfo);

	delete[] row;

	fclose(outfile);

	return true;
}


bool ImageSerializer_jpeg::read_supported() const {
	return true;
}

bool ImageSerializer_jpeg::streams_supported() const {
	return false;
}


bool ImageSerializer_jpeg::query_image_size(const std::string& file_name, int& width, int& height) {
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;

	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&cinfo);

	FILE *f = fopen(file_name.c_str(), "rb");

	if (f == NULL) {
		LOG(ERROR)
			<< "could not open file:" << file_name << std::endl;
		return false;
	}

	jpeg_stdio_src(&cinfo, f);
	if (jpeg_read_header(&cinfo, TRUE) == JPEG_HEADER_OK) {
		width = cinfo.image_width;
		height = cinfo.image_height;
		jpeg_destroy_decompress(&cinfo);
		fclose(f);
		return true;
	}

	fclose(f);
	return false;
}