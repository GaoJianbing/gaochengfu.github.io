// mse.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <math.h>

#define MAX_SIZE 5000


int changeScale(double x[], double *y, int N, int scale)
{
	if(scale <= 0) return 0;

	double s = x[0];
	for(int i=1;i<N;i++)
	{
		if(i%scale == 0)
		{
			y[(i-1)/scale] = s/scale;
			s = 0;
		}
		s += x[i];
	}

	int len_s = N/scale;
	if(len_s*scale < N)
	{
		y[(N-1)/scale] = s/(N - len_s*scale);
		len_s += 1;
	}
	else
	{
		y[(N-1)/scale] = s/scale;
	}

	return len_s;
}


double MSE(double x[], int N, int m)
{
	int k = N - m + 1;
	double Cmr = 0.0;

	if(k > 0)
	{
		double x_min = x[0];
		double x_max = x[0];
		for(int i1=1;i1<N;i1++)
		{
			if(x_min > x[i1])x_min = x[i1];
			if(x_max < x[i1])x_max = x[i1];
		}
		double r = (x_max - x_min)*0.25;

		for(int j=0;j<(N-m-1);j++)
		{
			for(int h=j+1;h<(N-m);h++)
			{
				double d = 0.0;
				for(int i=0;i<m;i++)
				{
					if(d < abs(x[j+i]-x[h+i])) d = abs(x[j+i]-x[h+i]);
				}
				if(d < r) Cmr += 2;
			}
		}
		Cmr /= (double)(k*k);
	}

	return Cmr;
}


int main(int argc, char* argv[])
{
	double x[] = {60, 65, 70, 65, 66, 68, 69, 62, 65, 66, 68, 59, 60, 64,
		66, 68, 59, 60, 64, 60, 64, 66, 68, 59, 60, 64, 59, 60, 64, 66, 69, 67};
	double y[32];

	int N = sizeof(x)/sizeof(double);
	int m = 4;

	for(int scale=1;scale<20;scale++)
	{
		int len_y = changeScale(x, y, N, scale);
		if(len_y <= m) break;
		for(int i=0;i<len_y;i++)
			printf("%f ", y[i]);
		double mse = MSE(y, len_y, m);
		printf("\n--------------------------------------------------------\n");
		printf("%f\n\n\n", mse);
	}

	return 0;
}

