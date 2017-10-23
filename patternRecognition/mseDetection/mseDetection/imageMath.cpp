#include "imageMath.h"

//Calculates the 
float mse(my_image_comp * big, my_image_comp * small, int xAxis, int yAxis) {
	int xCounter = 0;
	int yCounter = 0;

	int width = small->width;
	int height = small->height;

	int start = yAxis * big->width + xAxis;
	//printf("Start: %d", start);
	float sum = 0;
	float diff;

	float picOne;
	float picTwo;
	//printf("Big Width: %d \n", big->width);
	//printf("Big height: %d \n", big->height);
	
	for (yCounter = 0; yCounter < height; yCounter++) {
		//printf("yCounter: %d \n", yCounter);
		for (xCounter = 0; xCounter < width; xCounter++) {
			picOne = big->buf[(yAxis + yCounter) * big->width + xCounter + xAxis];
			picTwo = small->buf[yCounter * width + xCounter];
			//printf("Big location: %d, small location %d \n", (yAxis + yCounter) * big->width + xCounter + xAxis, yCounter * width + xCounter);
			//printf("PicOne: %f, PicTwo: %f", picOne, picTwo);
			diff = picOne - picTwo;
			//printf(" Entry %d: %f \n", yCounter * width + xCounter, diff);
			sum = sum + (float) pow((double) diff, 2);
		}
	}
	//printf("%f \n", sum);
	//printf("Big Buffer: %f, Small Buffer: %f \n", big->buf[0], small->buf[0]);
	return sum;
}

//Finds all the MSEs and then
float smallestMSE(my_image_comp * big, my_image_comp * small, int * smallestCoord) {
	int xCounter;
	int yCounter;

	float minMSE = 100000000;
	float MSE;

	for (yCounter = 0; yCounter < (big->height - small->height); yCounter++) {
		for (xCounter = 0; xCounter < (big->width - small->width); xCounter++) {
			//printf("X: %d, Y: %d \n", xCounter, yCounter);
			//printf("Big Buffer: %f, Small Buffer: %f \n", big->buf[0], small->buf[0]);
			MSE = mse(big, small, xCounter, yCounter);
			//printf("Big Buffer 3: %f, Small Buffer 3: %f \n", big.buf[0], small.buf[0]);
			if (MSE < minMSE) {
				smallestCoord[0] = xCounter;
				smallestCoord[1] = yCounter;
				minMSE = MSE;
			}
			//printf("Big Buffer 2: %f, Small Buffer 2: %f \n", big.buf[0], small.buf[0]);
		}
	}
	return minMSE;
};

float correlation(my_image_comp * big, my_image_comp * small, int xAxis, int yAxis) {
	int xCounter;
	int yCounter;

	int width = small->width;
	int height = small->height;

	int start = yAxis * big->stride + xAxis;

	float sum = 0;
	float product;

	for (yCounter = 0; yCounter < height; yCounter++) {
		for (xCounter = 0; xCounter < width; xCounter++) {
			product = big->handle[start + yCounter*big->stride + xCounter] * small->handle[yCounter * small->stride + xCounter];
			sum = sum + product;
		}
	}

	return sum;
}

float maxCorrelation(my_image_comp * big, my_image_comp * small, int * maxCorr) {
	int xCounter;
	int yCounter;

	float biggestCorr = 0;
	float Corr;

	for (yCounter = 0; yCounter < (big->height - small->height); yCounter++) {
		for (xCounter = 0; xCounter < (big->height - small->height); xCounter++) {
			Corr = correlation(big, small, xCounter, yCounter);
			if (Corr > biggestCorr) {
				maxCorr[0] = xCounter;
				maxCorr[1] = yCounter;
				biggestCorr = Corr;
			}
		}
	}
	return biggestCorr;
};

float correlationFloat(my_image_comp * big, my_image_comp * small, int xAxis, int yAxis, float * bigCorr, float * smallCorr) {
	int xCounter;
	int yCounter;

	int width = small->width;
	int height = small->height;

	int start = yAxis * big->stride + xAxis;

	float sum = 0;
	float product;

	for (yCounter = 0; yCounter < height; yCounter++) {
		for (xCounter = 0; xCounter < width; xCounter++) {
			product = bigCorr[start + yCounter*big->stride + xCounter] * smallCorr[yCounter * small->stride + xCounter];
			sum = sum + product;
		}
	}

	return sum;
}


float maxFloatCorrelation(my_image_comp * big, my_image_comp * small, int * maxCorr, float * bigCorr, float * smallCorr) {
	int xCounter;
	int yCounter;

	float biggestCorr = 0;
	float Corr;

	for (yCounter = 0; yCounter < (big->height - small->height); yCounter++) {
		for (xCounter = 0; xCounter < (big->height - small->height); xCounter++) {
			Corr = correlationFloat(big, small, xCounter, yCounter, bigCorr, smallCorr);
			if (Corr > biggestCorr) {
				maxCorr[0] = xCounter;
				maxCorr[1] = yCounter;
				biggestCorr = Corr;
			}
		}
	}
	return biggestCorr;
};