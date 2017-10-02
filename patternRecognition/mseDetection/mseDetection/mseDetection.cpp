// imageSizeReducer.cpp : Defines the entry point for the console application.
//

#include "io_bmp.h"
//#include "image_comps.h"

#include "imageMath.h"

#pragma warning(disable:4996)


/*****************************************************************************/
/*                  my_image_comp::perform_boundary_extension                */
/*****************************************************************************/

void my_image_comp::perform_boundary_extension()
{
	int r, c;

	// First extend upwards
	float *first_line = buf;
	for (r = 1; r <= border; r++)
		for (c = 0; c < width; c++)
			first_line[-r*stride + c] = first_line[c];

	// Now extend downwards
	float *last_line = buf + (height - 1)*stride;
	for (r = 1; r <= border; r++)
		for (c = 0; c < width; c++)
			last_line[r*stride + c] = last_line[c];

	// Now extend all rows to the left and to the right
	float *left_edge = buf - border*stride;
	float *right_edge = left_edge + width - 1;
	for (r = height + 2 * border; r > 0; r--, left_edge += stride, right_edge += stride)
		for (c = 1; c <= border; c++)
		{
			left_edge[-c] = left_edge[0];
			right_edge[c] = right_edge[0];
		}
}


/* ========================================================================= */
/*                              Global Functions                             */
/* ========================================================================= */

/*****************************************************************************/
/*                                apply_filter                               */
/*****************************************************************************/

void apply_filter(my_image_comp *in, my_image_comp *out)
{
#define FILTER_EXTENT 4
#define FILTER_DIM (2*FILTER_EXTENT+1)
#define FILTER_TAPS (FILTER_DIM*FILTER_DIM)

	// Create the filter kernel as a local array on the stack, which can accept
	// row and column indices in the range -FILTER_EXTENT to +FILTER_EXTENT.
	float filter_buf[FILTER_TAPS];
	float *mirror_psf = filter_buf + (FILTER_DIM*FILTER_EXTENT) + FILTER_EXTENT;
	// `mirror_psf' points to the central tap in the filter
	int r, c;
	for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++)
		for (c = -FILTER_EXTENT; c <= FILTER_EXTENT; c++)
			mirror_psf[r*FILTER_DIM + c] = 1.0F / FILTER_TAPS;

	// Check for consistent dimensions
	assert(in->border >= FILTER_EXTENT);
	assert((out->height <= in->height) && (out->width <= in->width));

	// Perform the convolution
	for (r = 0; r < out->height; r++)
		for (c = 0; c < out->width; c++)
		{
			float *ip = in->buf + r*in->stride + c;
			float *op = out->buf + r*out->stride + c;
			float sum = 0.0F;
			for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++)
				for (int x = -FILTER_EXTENT; x <= FILTER_EXTENT; x++)
					sum += ip[y*in->stride + x] * mirror_psf[y*FILTER_DIM + x];
			*op = sum;
		}
}
/*****************************************************************************/
/*                                zeroIntensity                               */
/*****************************************************************************/
void zeroIntensity(my_image_comp *out) {
	int r, c;
	float sum = 0.0F;
	float average = 0.0F;

		for (r = 0; r < out->height; r++) {
			for (c = 0; c < out->width; c++) {
				float *op = out->buf + r*out->stride + c;
				sum += op[0];
			}
		}
	average = sum / (out->height * out->width);

	for (r = 0; r < out->height; r++) {
		for (c = 0; c < out->width; c++) {
			float *op = out->buf + r*out->stride + c;
			*op = *op -average;
		}
	}
}


int
main(int argc, char *argv[])
{
	if (argc != 4)
	{
		fprintf(stderr, "Usage: %s <in bmp file> <small sample file> <mode>\n", argv[0]);
		return -1;
	}

	int err_code = 0;
	try {
		// Read the large image
		bmp_in in;
		if ((err_code = bmp_in__open(&in, argv[1])) != 0)
			throw err_code;

		int width = in.cols, height = in.rows;

		int n, num_comps = in.num_components;
		my_image_comp *input_comps = new my_image_comp[num_comps];
		for (n = 0; n < num_comps; n++)
			input_comps[n].init(height, width, 0); // Leave a border of 0

		int r, c; // Declare row index
		io_byte *line = new io_byte[width*num_comps];
		for (r = height - 1; r >= 0; r--)
		{ // "r" holds the true row index we are reading, since the image is
		  // stored upside down in the BMP file.
			if ((err_code = bmp_in__get_line(&in, line)) != 0)
				throw err_code;
			for (n = 0; n < num_comps; n++)
			{
				io_byte *src = line + n; // Points to first sample of component n
				float *dst = input_comps[n].buf + r * input_comps[n].stride;
				for (int c = 0; c < width; c++, src += num_comps) {
					dst[c] = (float)*src; // The cast to type "float" is not
										  // strictly required here, since bytes can always be
										  // converted to floats without any loss of information.

				}
			}
		}
		bmp_in__close(&in);



		printf("Opening rectangle \n");
		//Read the sample Comps
		bmp_in rectangle;
		if ((err_code = bmp_in__open(&rectangle, argv[2])) != 0)
			throw err_code;
		printf("I'm here \n");
		width = rectangle.cols, height = rectangle.rows;

		n, num_comps = rectangle.num_components;
		my_image_comp *rectangle_comps = new my_image_comp[num_comps];
		for (n = 0; n < num_comps; n++)
			rectangle_comps[n].init(height, width, 0); // Leave a border of 0

		r, c; // Declare row index
		io_byte *line2 = new io_byte[width*num_comps];
		for (r = height - 1; r >= 0; r--)
		{ // "r" holds the true row index we are reading, since the image is
	  // stored upside down in the BMP file.
			if ((err_code = bmp_in__get_line(&rectangle, line2)) != 0)
				throw err_code;
			for (n = 0; n < num_comps; n++)
			{
				io_byte *src = line + n; // Points to first sample of component n
				float *dst = rectangle_comps[n].buf + r * rectangle_comps[n].stride;
				for (int c = 0; c < width; c++, src += num_comps) {
					dst[c] = (float)*src; // The cast to type "float" is not
									  // strictly required here, since bytes can always be
									  // converted to floats without any loss of information.

				}
			}
		}
		bmp_in__close(&in);

		int idealCoord[2];
		int mode = atoi(argv[3]);
		if (mode == 1) {
			//Fill out MSE here
			float smallMSE;
			smallMSE = smallestMSE(input_comps, rectangle_comps, idealCoord);
			printf("The smallest MSE is %f", smallMSE);
		}
		else if(mode == 2) {
			float largeCorr;
			largeCorr = maxCorrelation(input_comps, rectangle_comps, idealCoord);
			printf("The largest correlation is %f", largeCorr);
		}

		printf("The ideal coordinates are at: X:%d Y:%d", idealCoord[0], idealCoord[1]);

		delete[] line;
		delete[] line2;
		delete[] input_comps;
		delete[] rectangle_comps;
	}
	catch (int exc) {
		if (exc == IO_ERR_NO_FILE)
			fprintf(stderr, "Cannot open supplied input or output file.\n");
		else if (exc == IO_ERR_FILE_HEADER)
			fprintf(stderr, "Error encountered while parsing BMP file header.\n");
		else if (exc == IO_ERR_UNSUPPORTED)
			fprintf(stderr, "Input uses an unsupported BMP file format.\n  Current "
				"simple example supports only 8-bit and 24-bit data.\n");
		else if (exc == IO_ERR_FILE_TRUNC)
			fprintf(stderr, "Input or output file truncated unexpectedly.\n");
		else if (exc == IO_ERR_FILE_NOT_OPEN)
			fprintf(stderr, "Trying to access a file which is not open!(?)\n");
		return -1;
	}
	return 0;
}
