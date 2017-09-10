#include "imageMath.h"

float sinc(float x) {
	float y = sin(x)/x;
	return y;
}

//Interpolation kernel, in this case, m and n are two co-ordinates for the input vector
float q(float m, float n) {
	float sOne = (0.4F) * m;
	float sTwo = (0.4F) * n;
	float ans = (0.16F) * sinc(m) * sinc(n);
	return ans;
}