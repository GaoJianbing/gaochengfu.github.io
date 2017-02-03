// test_fft.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

// testdll.cpp : Defines the entry point for the DLL application.

#include "stdafx.h" 

///////////////////////////////////////////////
//
// 快速傅立叶变换 Fast Fourier Transform
// By rappizit@yahoo.com.cn
// 2007-07-20
// 版本 2.0
// 改进了《算法导论》的算法，旋转因子取 ωn-kj  (ωnkj 的共轭复数)
// 且只计算 n / 2 次，而未改进前需要计算 (n * lg n) / 2 次。
// 
////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const int N = 1024;
const double PI = 3.1416;

inline void swap (double &a, double &b)
{
    double t;
    t = a;
    a = b;
    b = t;
}

void bitrp (double xreal [], double ximag [], int n)
{
    // 位反转置换 Bit-reversal Permutation
    int i, j, a, b, p;

    for (i = 1, p = 0; i < n; i *= 2)
	{
        p ++;
	}
    for (i = 0; i < n; i ++)
	{
        a = i;
        b = 0;
        for (j = 0; j < p; j ++)
		{
            b = (b << 1) + (a & 1);    // b = b * 2 + a % 2;
            a >>= 1;        // a = a / 2;
		}
        if ( b > i)
		{
            swap (xreal [i], xreal [b]);
            swap (ximag [i], ximag [b]);
		}
	}
}

void FFT(double xreal [], double ximag [], int n)
{
    // 快速傅立叶变换，将复数 x 变换后仍保存在 x 中，xreal, ximag 分别是 x 的实部和虚部
    double wreal [N / 2], wimag [N / 2], treal, timag, ureal, uimag, arg;
    int m, k, j, t, index1, index2;
    
    bitrp (xreal, ximag, n);

    // 计算 1 的前 n / 2 个 n 次方根的共轭复数 W'j = wreal [j] + i * wimag [j] , j = 0, 1, ... , n / 2 - 1
    arg = - 2 * PI / n;
    treal = cos (arg);
    timag = sin (arg);
    wreal [0] = 1.0;
    wimag [0] = 0.0;
    for (j = 1; j < n / 2; j ++)
	{
        wreal [j] = wreal [j - 1] * treal - wimag [j - 1] * timag;
        wimag [j] = wreal [j - 1] * timag + wimag [j - 1] * treal;
	}

    for (m = 2; m <= n; m *= 2)
	{
        for (k = 0; k < n; k += m)
		{
            for (j = 0; j < m / 2; j ++)
			{
                index1 = k + j;
                index2 = index1 + m / 2;
                t = n * j / m;    // 旋转因子 w 的实部在 wreal [] 中的下标为 t
                treal = wreal [t] * xreal [index2] - wimag [t] * ximag [index2];
                timag = wreal [t] * ximag [index2] + wimag [t] * xreal [index2];
                ureal = xreal [index1];
                uimag = ximag [index1];
                xreal [index1] = ureal + treal;
                ximag [index1] = uimag + timag;
                xreal [index2] = ureal - treal;
                ximag [index2] = uimag - timag;
			}
		}
	}
}

void  IFFT (double xreal [], double ximag [], int n)
{
    // 快速傅立叶逆变换
    double wreal [N / 2], wimag [N / 2], treal, timag, ureal, uimag, arg;
    int m, k, j, t, index1, index2;
    
    bitrp (xreal, ximag, n);

    // 计算 1 的前 n / 2 个 n 次方根 Wj = wreal [j] + i * wimag [j] , j = 0, 1, ... , n / 2 - 1
    arg = 2 * PI / n;
    treal = cos (arg);
    timag = sin (arg);
    wreal [0] = 1.0;
    wimag [0] = 0.0;
    for (j = 1; j < n / 2; j ++)
	{
        wreal [j] = wreal [j - 1] * treal - wimag [j - 1] * timag;
        wimag [j] = wreal [j - 1] * timag + wimag [j - 1] * treal;
	}

    for (m = 2; m <= n; m *= 2)
	{
        for (k = 0; k < n; k += m)
		{
            for (j = 0; j < m / 2; j ++)
			{
                index1 = k + j;
                index2 = index1 + m / 2;
                t = n * j / m;    // 旋转因子 w 的实部在 wreal [] 中的下标为 t
                treal = wreal [t] * xreal [index2] - wimag [t] * ximag [index2];
                timag = wreal [t] * ximag [index2] + wimag [t] * xreal [index2];
                ureal = xreal [index1];
                uimag = ximag [index1];
                xreal [index1] = ureal + treal;
                ximag [index1] = uimag + timag;
                xreal [index2] = ureal - treal;
                ximag [index2] = uimag - timag;
			}
		}
	}

    for (j=0; j < n; j ++)
	{
        xreal [j] /= n;
        ximag [j] /= n;
	}
}

void FFT_test ()
{
    char inputfile [] = "d:/input.txt";    // 从文件 input.txt 中读入原始数据
    char outputfile [] = "d:/output.txt";    // 将结果输出到文件 output.txt 中
    double xreal [N], ximag [N];
    int i;
    FILE *input, *output;

    if (!(input = fopen (inputfile, "r")))
	{
        printf ("Cannot open file. ");
        exit (1);
	}
    if (!(output = fopen (outputfile, "w")))
	{
        printf ("Cannot open file. ");
        exit (1);
	}

	// 初始0填充
	for (i = 0; i < N; i ++)
		xreal[i] = ximag[i] = 0;
    i = 0;
	char real_tmp[50];
    while ((fscanf(input, "%s\n", real_tmp)) != EOF)
	{
		xreal[i] = atof(real_tmp);
        i++;
		if(i >= N) break;
	}

	// 镜面周期延拓
    if(i < N)
	{
		int k = 0;
		for(int j=i;j<N;j++)
		{
			k = j%i;
			if((j/i)%2 == 1) k = i - 1 - k;
			xreal[j] = xreal[k];
			ximag[j] = ximag[k];
		}
	}

    FFT(xreal, ximag, N);
	double d_module = 0;
    for (i = 0; i < N; i ++)
	{
		d_module = sqrt(xreal[i]*xreal[i] + ximag[i]*ximag[i]);
        fprintf (output, "%8.4f\n", d_module);
	}

    if ( fclose (input))
	{
        printf ("File close error. ");
        exit (1);
	}
    if ( fclose (output))
	{
        printf ("File close error. ");
        exit (1);
	}
}

int main(int argc, char* argv[])
{
    FFT_test ();

	return 0;
}
