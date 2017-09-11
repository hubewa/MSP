#include "imageMath.h"

# define M_PI           3.14159265358979323846 

float sinc(float x) {
	float y = sin(x)/x;
	return y;
}

//Interpolation kernel, in this case, m and n are two co-ordinates for the input vector
float q(float m) {
	float ans = (0.4F) * m;
	return ans;
}

float hanning(float s, float tau) {
	return (1.0F + cos((M_PI * s)/tau) / 2.0F);
}