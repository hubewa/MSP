#include <math.h>
#include "io_bmp.h"
#include "image_comps.h"

float mse(my_image_comp *big, my_image_comp *small, int xAxis, int yAxis);

float smallestMSE(my_image_comp *big, my_image_comp * small, int * smallestCoord);

float correlation(my_image_comp * big, my_image_comp * small, int xAxis, int yAxis);

float maxCorrelation(my_image_comp * big, my_image_comp * small, int * maxCorr);

float correlationFloat(my_image_comp * big, my_image_comp * small, int xAxis, int yAxis, float * bigCorr, float * smallCorr);

float maxFloatCorrelation(my_image_comp * big, my_image_comp * small, int * maxCorr, float * bigCorr, float * smallCorr);
