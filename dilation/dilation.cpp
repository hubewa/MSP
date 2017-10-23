/*****************************************************************************/
// File: filtering_main.cpp
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include "io_bmp.h"
#include "image_comps.h"

#define _SECURE_SCL_DEPRECATE 0
#include <algorithm>
#include <math.h>

/* ========================================================================= */
/*                 Implementation of `my_image_comp' functions               */
/* ========================================================================= */

/*****************************************************************************/
/*                  my_image_comp::perform_boundary_extension                */
/*****************************************************************************/

void my_image_comp::perform_boundary_extension()
{
  int r, c;

  // First extend upwards
  float *first_line = buf;
  for (r=1; r <= border; r++)
    for (c=0; c < width; c++)
      first_line[-r*stride+c] = first_line[r*stride+c];

  // Now extend downwards
  float *last_line = buf+(height-1)*stride;
  for (r=1; r <= border; r++)
    for (c=0; c < width; c++)
      last_line[r*stride+c] = last_line[-r*stride+c];

  // Now extend all rows to the left and to the right
  float *left_edge = buf-border*stride;
  float *right_edge = left_edge + width - 1;
  for (r=height+2*border; r > 0; r--, left_edge+=stride, right_edge+=stride)
    for (c=1; c <= border; c++)
      {
        left_edge[-c] = left_edge[c];
        right_edge[c] = right_edge[-c];
      }
}


/* ========================================================================= */
/*                              Global Functions                             */
/* ========================================================================= */

/*****************************************************************************/
/*                                apply_filter                               */
/*****************************************************************************/
//Try to find the largest value inside a to set the correct boundary.
int arrayLargest(int * a, int numElements) {
	int maxValue = a[0];
	int minValue = a[0];
	int i;
	int largestValue;

	for (i = 1; i < numElements; i++) {
		if (a[i] > maxValue) {
			maxValue = a[i];
		}
		else if (a[i] < minValue) {
			minValue = a[i];
		}
	}

	if (abs(maxValue) > abs(minValue)) {
		largestValue = abs(maxValue);
	}
	else {
		largestValue = abs(minValue);
	}
	return largestValue;
}


void apply_filter(my_image_comp *in, my_image_comp *out, int * a, int numElements)
{
  // Create the filter kernel as a local array on the stack, which can accept
  // row and column indices in the range -FILTER_EXTENT to +FILTER_EXTENT.
  //float filter_buf[FILTER_TAPS];
  //float *mirror_psf = filter_buf+(FILTER_DIM*FILTER_EXTENT)+FILTER_EXTENT;
          // `mirror_psf' points to the central tap in the filter
  int r, c;
 /* for (r=-FILTER_EXTENT; r <= FILTER_EXTENT; r++)
    for (c=-FILTER_EXTENT; c <= FILTER_EXTENT; c++)
      mirror_psf[r*FILTER_DIM+c] = 1.0F / FILTER_TAPS;*/

  // Check for consistent dimensions
 // assert(in->border >= FILTER_EXTENT);
  assert((out->height <= in->height) && (out->width <= in->width));

  //Do Erosion
  for (r=0; r < out->height; r++)
    for (c=0; c < out->width; c++)
      {
		float *ip;
		float val = 0;
        float *op = out->buf + r*out->stride + c;
        float sum = 0.0F;

		int arrayIndex;
		int currA;
		int vectorNum = numElements / 2; //Find the number of vectors involved

		for (arrayIndex = 0; arrayIndex < vectorNum; arrayIndex++) {
			currA = a[arrayIndex * 2] * in->stride + a[arrayIndex * 2 + 1]; //Find the value of the vector
			//printf("Array Index: (%d,%d) \n", a[arrayIndex * 2], a[arrayIndex * 2 + 1]);
			//printf("Current coordinate check: %d \n", r*in->stride + c + currA);
			if (r*in->stride + c + currA < 0) { // Check if it exceeds top
				val = 255.0F;
				//printf("Out of bounds \n");
			}
			else if (r*in->stride + c + currA > (in->height*in->stride)){ //Check it doesn't exceed end of the array
				val = 255.0F;
				//printf("Out of bounds \n");
			}
			else {
				ip = in->buf + r*in->stride + c + currA;
				if (*ip == 255) {
					val = 255.0F;
				}
				//printf("Zero Entry \n");
			}
			//printf("End Check \n");
			//printf("Value: %d, val");
		}

		*op = val;
      }
}

/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int
  main(int argc, char *argv[])
{
  if (argc < 5)
    {
      fprintf(stderr,"Usage: %s <in bmp file> <out bmp file> mode vector\n",argv[0]);
      return -1;
    }

  int err_code=0;
  try {
      //Get threshold value
	  
	  // Read the input image
      bmp_in in;
      if ((err_code = bmp_in__open(&in,argv[1])) != 0)
        throw err_code;

      int width = in.cols, height = in.rows;
      int n, num_comps = in.num_components;
      my_image_comp *input_comps = new my_image_comp[num_comps];

	  //sets up the vector a
	  int i;

	  int mode = atoi(argv[3]);
	  int *a;

	  int numElements = 0;

	  if (mode == 0) { //Vector mode 
		  a = new int[argc - 4];
		  for (i = 0; i < argc - 4; i++) {
			  a[i] = atoi(argv[i + 4]);
		  }
		  numElements = argc - 4;
	  }
	  else { //Radius mode
		  int radius = atoi(argv[4]);
		  printf("Radius: %d \n", radius);
		  a = new int[(2*radius+1) * (2 * radius + 1) * 2]; //Times 2 is for the x and y coordinates
		  int r;
		  int c;
		  int counter = 0;
		  for (r = -radius; r <= radius; r++) {
			  for (c = -radius; c <= radius; c++) {
				  if (pow(r, 2) + pow(c, 2) <= pow(radius,2)) {
					  a[counter * 2] = r;
					  a[counter * 2 + 1] = c;
					  printf("(r,c) = (%d,%d) \n", a[counter * 2], a[counter * 2+1]);
					  counter++;
				  }
			  }
		  }
		  numElements = counter*2;
	  }
	  printf("Number of Elements: %d \n", numElements);
	  int maxValue = arrayLargest(a, numElements);
	  printf("max Value = %d", maxValue);

      for (n=0; n < num_comps; n++)
        input_comps[n].init(height,width,4); // Leave a border of 4
      
      int r; // Declare row index
      io_byte *line = new io_byte[width*num_comps];
      for (r=height-1; r >= 0; r--)
        { // "r" holds the true row index we are reading, since the image is
          // stored upside down in the BMP file.
          if ((err_code = bmp_in__get_line(&in,line)) != 0)
            throw err_code;
          for (n=0; n < num_comps; n++)
            {
              io_byte *src = line+n; // Points to first sample of component n
              float *dst = input_comps[n].buf + r * input_comps[n].stride;
              for (int c=0; c < width; c++, src+=num_comps)
                dst[c] = (float) *src; // The cast to type "float" is not
                      // strictly required here, since bytes can always be
                      // converted to floats without any loss of information.
            }
        }
      bmp_in__close(&in);

      // Allocate storage for the filtered output
      my_image_comp *output_comps = new my_image_comp[num_comps];
      for (n=0; n < num_comps; n++)
        output_comps[n].init(height,width,maxValue); // Don't need a border for output

      // Process the image, all in floating point (easy)
      for (n=0; n < num_comps; n++)
        input_comps[n].perform_boundary_extension();
      for (n=0; n < num_comps; n++)
        apply_filter(input_comps+n,output_comps+n, a, numElements);

      // Write the image back out again
      bmp_out out;
      if ((err_code = bmp_out__open(&out,argv[2],width,height,num_comps)) != 0)
        throw err_code;
      for (r=height-1; r >= 0; r--)
        { // "r" holds the true row index we are writing, since the image is
          // written upside down in BMP files.
          for (n=0; n < num_comps; n++)
            {
              io_byte *dst = line+n; // Points to first sample of component n
              float *src = output_comps[n].buf + r * output_comps[n].stride;
              for (int c=0; c < width; c++, dst+=num_comps)
                *dst = (io_byte) src[c]; // The cast to type "io_byte" is
                      // required here, since floats cannot generally be
                      // converted to bytes without loss of information.  The
                      // compiler will warn you of this if you remove the cast.
                      // There is in fact not the best way to do the
                      // conversion.  You should fix it up in the lab.
            }
          bmp_out__put_line(&out,line);
        }
      bmp_out__close(&out);
      delete[] line;
      delete[] input_comps;
      delete[] output_comps;
    }
  catch (int exc) {
      if (exc == IO_ERR_NO_FILE)
        fprintf(stderr,"Cannot open supplied input or output file.\n");
      else if (exc == IO_ERR_FILE_HEADER)
        fprintf(stderr,"Error encountered while parsing BMP file header.\n");
      else if (exc == IO_ERR_UNSUPPORTED)
        fprintf(stderr,"Input uses an unsupported BMP file format.\n  Current "
                "simple example supports only 8-bit and 24-bit data.\n");
      else if (exc == IO_ERR_FILE_TRUNC)
        fprintf(stderr,"Input or output file truncated unexpectedly.\n");
      else if (exc == IO_ERR_FILE_NOT_OPEN)
        fprintf(stderr,"Trying to access a file which is not open!(?)\n");
      return -1;
    }
  return 0;
}
