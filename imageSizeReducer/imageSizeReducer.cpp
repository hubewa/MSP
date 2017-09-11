// imageSizeReducer.cpp : Defines the entry point for the console application.
//

#include "io_bmp.h"
#include "image_comps.h"

#include "imageMath.h"

//Boundary Extension Function
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

//Apply Filter
void apply_filter(my_image_comp *in, my_image_comp *out, int filterExtent)
{
	//#define FILTER_EXTENT 4
	//#define FILTER_DIM (2*FILTER_EXTENT+1)
	//#define filterTaps (FILTER_DIM*FILTER_DIM)

	int filterDim = 2 * filterExtent + 1;
	int filterTaps = filterDim*filterDim;

	// Create the filter kernel as a local array on the stack, which can accept
	// row and column indices in the range -filterExtent to +filterExtent.
	float *filter_buf;

	filter_buf = new float[filterDim];
	float *mirror_psf = filter_buf + (filterDim*filterExtent) + filterExtent;
	// `mirror_psf' points to the central tap in the filter

	float *temp_buf;
	temp_buf = new float[out->height * in->stride];

	int r, c; //Row and Column

	
	int oldHeight = in->height;
	int newWidth = out->width;

	int m, n; //Row and column of PSF
	int row2;

	int level;

	// Perform the convolution - rows
	for (r = 0; r < out->height; r++) {
		//Build PSF for rows
		if (r % 2 == 1) {
			row2 = (r - 1) / 2;
			for (n = -filterExtent; n <= filterExtent; n++) {
				mirror_psf[n + filterExtent] = hanning(n, filterExtent) * q(n - 5 * row2 - 2.5F);
			}
		}
		else {
			row2 = r / 2;
			for (n = -filterExtent; n <= filterExtent; n++) {
				level = n + filterExtent;
				mirror_psf[n + filterExtent] = hanning(n, filterExtent) * q(n - 5 * row2);
			}
		}
		for (c = 0; c < out->width; c++)
		{
			float *ip = in->buf + r*in->stride + c;
			float *temp = temp_buf + r*in->stride + c;
			float sum = 0.0F;
			//Inner Product :D
			if (r*in->stride + c < 0 || r*in->stride + c >= out->height * out->stride) {
				sum += 0;
			}
			else {
				sum += ip[r*in->stride + c] * mirror_psf[c + filterExtent];
			}
			*temp = sum;
		}
	}
	//Convolution for Columns
	for (c = 0; c < out->width; c++) {
		if (c % 2 == 1) {
			row2 = (c - 1) / 2;
			for (n = -filterExtent; n <= filterExtent; n++) {
				mirror_psf[n + filterExtent] = q(n - 5 * row2 - 2.5F);
			}
		}
		else {
			row2 = c / 2;
			for (n = -filterExtent; n <= filterExtent; n++) {
				mirror_psf[n + filterExtent] = q(n - 5 * row2);
			}
		}
		for (r = 0; r < out->height-1; r++){
			float *ip = temp_buf + r*out->stride + c;
			float *op = out->buf + r*out->stride + c;
			float sum = 0.0F;
			//Inner Product :D
			for (int y = -filterExtent; y <= filterExtent; y++)
				// That is, we're looking for even terms 
				if(r*in->stride + y < 0 || r*in->stride + y >= out->height * in->stride){
					sum += 0;
				} else {
					sum += ip[y] * mirror_psf[y + filterExtent];
				}
			*op = sum;
		}
	}
	delete [] filter_buf;
	delete [] temp_buf;
}

int
main(int argc, char *argv[])
{
	//Here we look for Four Arguments
	if (argc != 4)
	{
		fprintf(stderr, "Usage: %s <in bmp file> <out bmp file> <H>\n", argv[0]);
		return -1;
	}

	int err_code = 0;
	try {
		// Read the input image
		bmp_in in;
		int support = atoi(argv[3]);
		if ((err_code = bmp_in__open(&in, argv[1])) != 0)
			throw err_code;

		int width = in.cols, height = in.rows;

		int outWidth = ceil(0.4F * width);
		int outHeight = ceil(0.4F * height);
		int copy;

		int n, num_comps = in.num_components;
		my_image_comp *input_comps = new my_image_comp[num_comps];
		for (n = 0; n < num_comps; n++)
			input_comps[n].init(height, width, 4); // Leave a border of 4

		int r,c; // Declare row index
		io_byte *line = new io_byte[width*num_comps];
		io_byte *outLine = new io_byte[outWidth*num_comps];
		for (r = outHeight - 1; r >= 0; r--)
		{ // "r" holds the true row index we are reading, since the image is
		  // stored upside down in the BMP file.
			if ((err_code = bmp_in__get_line(&in, line)) != 0)
				throw err_code;
			for (n = 0; n < num_comps; n++)
			{
				io_byte *src = line + n; // Points to first sample of component n
				float *dst = input_comps[n].buf + r * input_comps[n].stride;
				for (int c = 0; c < outWidth; c++, src += num_comps) {
					dst[c] = (float)*src; // The cast to type "float" is not
										  // strictly required here, since bytes can always be
										  // converted to floats without any loss of information.

				}
			}
		}
		bmp_in__close(&in);

		// Allocate storage for the filtered output
		my_image_comp *output_comps = new my_image_comp[num_comps];
		for (n = 0; n < num_comps; n++)
			output_comps[n].init(outHeight, outWidth, 0); // Don't need a border for output

													// Process the image, all in floating point (easy)
		for (n = 0; n < num_comps; n++)
			input_comps[n].perform_boundary_extension();
		for (n = 0; n < num_comps; n++)
			apply_filter(input_comps + n, output_comps + n, support);

		// Write the image back out again
		bmp_out out;
		if ((err_code = bmp_out__open(&out, argv[2], outWidth, outHeight, num_comps)) != 0)
			throw err_code;
		for (r = outHeight - 1; r >= 0; r--)
		{ // "r" holds the true row index we are writing, since the image is
		  // written upside down in BMP files.
			for (n = 0; n < num_comps; n++)
			{
				io_byte *dst = line + n; // Points to first sample of component n
				float *src = output_comps[n].buf + r * output_comps[n].stride;
				for (int c = 0; c < outWidth; c++, dst += num_comps) {
					*dst = (io_byte)src[c]; // The cast to type "io_byte" is
											// required here, since floats cannot generally be
											// converted to bytes without loss of information.  The
											// compiler will warn you of this if you remove the cast.
											// There is in fact not the best way to do the
											// conversion.  You should fix it up in the lab.
					copy = r*outWidth + c;
				}
			}
			bmp_out__put_line(&out, line);
		}
		bmp_out__close(&out);
		delete[] line;
		delete[] input_comps;
		delete[] output_comps;
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
