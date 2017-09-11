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

float fract(float x) {
	return x - trunc(x);
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
	float *mirror_psf = filter_buf + filterExtent;
//	float *mirror_psf = filter_buf + (filterDim*filterExtent) + filterExtent;
	// `mirror_psf' points to the central tap in the filter

	float *temp_buf;
	temp_buf = new float[out->height * in->stride];

	int r, c; //Row and Column

	
	int oldHeight = in->height;
	int newWidth = out->width;

	int m, n; //Row and column of PSF
	int row2;

	int level;
	float *qq = new float[filterDim];

	for (r = 0; r < in->height; r++) {
		printf("r=%d, in->buf[r,0]=%f\n", r, in->buf[r*in->stride]);
	}


	// Perform the convolution - rows
	for (c = 0; c < in->width; c++) {
		for (r = 0; r < out->height; r++) {
		//Build PSF for rows
		int is_even = r % 2 == 0;
		float ip_r = 2.5F * r;
		float *psf = mirror_psf;
		float sum = 0.0F;

		if (!is_even) {
			row2 = (r - 1) / 2;
			for (n = -filterExtent; n <= filterExtent; n++, psf++) {
//				mirror_psf[n] = hanning(n, filterDim) * q(n - 5 * row2 - 2.5F);
				if (n == 0) {
//					*psf = 0.;
					mirror_psf[n] = 0.;
//					printf("n=%d, psf=%f\n", n, mirror_psf[n]);
				}
				else if (n < 0) {
//					*psf = hanning(n, filterDim) * q(n + fract(ip_r));
					mirror_psf[n] = hanning(n, filterDim) * q(n + fract(ip_r));
//					printf("n=%d, hann=%f, q=%f, psf=%f\n", n, hanning(n, filterDim), q(n+fract(ip_r)), mirror_psf[n]);
//					qq[n + filterExtent] = q(n + fract(ip_r));
//					qq[n + filterExtent] = q(n - 5 * row2 - 2.5F);
				}
				else {
					mirror_psf[n] = hanning(n, filterDim) * q(n - fract(ip_r));
//					*psf = hanning(n, filterDim) * q(n - fract(ip_r));
//					printf("n=%d, hann=%f, q=%f, psf=%f\n", n, hanning(n, filterDim), q(n - fract(ip_r)), mirror_psf[n]);
					//					qq[n + filterExtent] = q(n + fract(ip_r));
//					qq[n + filterExtent] = q(n - 5 * row2 - 2.5F);
				}
			}
		}
		else {
			row2 = r / 2;
			for (n = -filterExtent; n <= filterExtent; n++) {
				level = n + filterExtent;
//				mirror_psf[n] = hanning(n, filterExtent) * q(n - 5 * row2);
//				*psf = hanning(n, filterDim) * q(n + fract(ip_r));
				mirror_psf[n] = hanning(n, filterDim) * q(n + fract(ip_r));
//				printf("n=%d, hann=%f, q=%f, psf=%f\n", n, hanning(n, filterDim), q(n + fract(ip_r)), mirror_psf[n]);
				//				qq[n + filterExtent] = q(n + fract(ip_r));
//				qq[n + filterExtent] = q(n - 5 * row2);
			}
		}
//		printf("ip_r=%f, qq=9*%f, other = 9*%f \n", fract(ip_r), q(n - 5 * row2), q(n + fract(ip_r)));
//		for (c = 0; c < in->width; c++){
//			for (m = -filterExtent; m <= filterExtent; m++) {
//				float *ip = in->buf + r*in->stride + c;
		{
				float *ip = in->buf + c;
				float *temp = temp_buf + r*in->stride + c;
				float *psf = mirror_psf;
				int k;

				//Inner Product :D
//				if ((r+m)*in->stride + c < 0 || (r + m)*in->stride + c >= out->height * out->stride) {
//					sum += 0;
//				}
//				else 
				{
					if (is_even)
						k = round(ip_r);
					else
						k = round(ip_r + 0.5F);
					//					printf("ip_size=%d, out->height=%d, in->height=%d, r=%d, ip_r=%d, ip_idx=%d, psf_idx=%d\n", in->stride*(in->height+2*in->border), out->height, in->height, r, (int) ip_r, (k + m)*in->stride, m + filterExtent);
				}
				for (m = -filterExtent; m <= filterExtent; m++) {

//					printf("ip=%f, in_stride=%d,  m_psf=%f, r=%d, k=%d, m=%d, sum=%f\n", ip[(k + m)*in->stride], in->stride, mirror_psf[m], r, k, m, sum);
					sum += ip[(k + m)*in->stride] * mirror_psf[m];
//					sum += ip[(k + m)*in->stride] * *(psf+m +filterExtent);
				}
				*temp = sum;
			}
		}
	}
	//Convolution for Columns
	for (r = 0; r < out->height; r++) {
		for (c = 0; c < out->width; c++) {
		int is_even = c % 2 == 0;
		float ip_c = 2.5F * c;
		float sum = 0.0F;

		if (!is_even) {
			row2 = (r - 1) / 2;
			for (n = -filterExtent; n <= filterExtent; n++) {
//				mirror_psf[n] = hanning(n, filterExtent) *q(n - 5 * row2 - 2.5F);
				if (n == 0)
					mirror_psf[n] = 0.;
				else if (n < 0)
					mirror_psf[n] = hanning(n, filterDim) * q(n + fract(ip_c));
				else
					mirror_psf[n] = hanning(n, filterDim) * q(n - fract(ip_c));
			}
		}
		else {
			row2 = r / 2;
			for (n = -filterExtent; n <= filterExtent; n++) {
				level = n + filterExtent;
//				mirror_psf[n] = hanning(n, filterExtent) *q(n - 5 * row2);
				mirror_psf[n] = hanning(n, filterDim) * q(n + fract(ip_c));
			}
		}
//		for (r = 0; r < out->height; r++){
//		float *ip = temp_buf + r*out->stride + c;
		float *ip = temp_buf;
		float *op = out->buf + r*out->stride + c;
//			float sum = 0.0F;
			//Inner Product :D
			for (int y = -filterExtent; y <= filterExtent; y++)
				// That is, we're looking for even terms 
//				if(r*in->stride + y < 0 || r*in->stride + y >= out->height * in->stride){
//					sum += 0;
//				}
//				else 
				{
				    int k;
					if (is_even)
						k = round(ip_c);
					else
						k = round(ip_c + 0.5F);

//					printf("ip=%f, out_stride=%d,  m_psf=%f, r=%d, c=%d, k=%d, y=%d, sum=%f\n", ip[r *out->stride + (k + y)], out->stride, mirror_psf[y], r, c, k, y, sum);
					sum += ip[r *in->stride + (k + y)] * mirror_psf[y];
				}
			*op = sum;
		}
	}
//	delete [] filter_buf;
//	delete [] temp_buf;
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
//		for (r = outHeight - 1; r >= 0; r--)
		for (r = height - 1; r >= 0; r--)
			{ // "r" holds the true row index we are reading, since the image is
		  // stored upside down in the BMP file.
			if ((err_code = bmp_in__get_line(&in, line)) != 0)
				throw err_code;
			for (n = 0; n < num_comps; n++)
			{
				io_byte *src = line + n; // Points to first sample of component n
				float *dst = input_comps[n].buf + r * input_comps[n].stride;
//				for (int c = 0; c < outWidth; c++, src += num_comps) {
				for (int c = 0; c < width; c++, src += num_comps) {
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
