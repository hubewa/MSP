#include "imageMath.h"
#include <stdio.h>

# define M_PI           3.14159265358979323846 

float sinc(float x) {
	return (x == 0.0) ? 1. : sin(x*M_PI) / (x*M_PI);
}

//Interpolation kernel, in this case, m and n are two co-ordinates for the input vector
float q(float m) {
	float ans = (2.5F) * sinc(m * 2.5F);
//	printf("%f, m: %f \n", ans, m);
	return ans;
}

float hanning(float s, float tau) {
//	return (1.0F - cos((M_PI * s)/tau)) / 2.0F;
	float hann = (1.0F + cos((M_PI * s) / tau)) / 2.0F;
//	printf("%f, s: %f, t: %f \n", hann, s, tau);
	return hann;
}