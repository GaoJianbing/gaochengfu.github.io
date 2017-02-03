// test_pluse.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"


//#include <jni.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <string.h>
//#include "android/log.h"
static const char *TAG = "PulseTakingSO";
//#define LOGI(fmt, args...) __android_log_print(ANDROID_LOG_INFO,  TAG, fmt, ##args)
//#define LOGD(fmt, args...) __android_log_print(ANDROID_LOG_DEBUG, TAG, fmt, ##args)
//#define printf(fmt, args...) __android_log_print(ANDROID_LOG_ERROR, TAG, fmt, ##args)

#define fu_mai 1
#define chen_mai 2
#define chi_mai 3
#define shu_mai 4
#define ji_mai 5
#define jie_mai 6
#define dai_mai 7
#define xian_mai 8
#define hua_mai 9
#define se_mai 10
#define other_mai -1

#define MAX_PLUSE_COUNT 500
#define USED_PARA_COUNT 18

const int sample_rate = 50; // 采样率
const int N = 1024;
const int len_x = N;
const double pi = 3.14159265;
const double PI = pi;

const int len_dn_fc = 7;
const int len_lp_fc = 60;
const int len_little_wave = 13;
double dn_fc[len_dn_fc] = { 0.0597210210218931, -0.142474780383904,
		0.215573464294102, 0.734360590135818, 0.215573464294102,
		-0.142474780383904, 0.0597210210218931 };
double lp_fc[len_lp_fc] = { 0.00239475689068246, 0.00317225220096019,
		0.00400953515140809, 0.00490366592568990, 0.00585102449408378,
		0.00684732791564099, 0.00788765588444296, 0.00896648435932763,
		0.0100777270078738, 0.0112147840891489, 0.0123705982972578,
		0.0135377169905568, 0.0147083601409260, 0.0158744932550581,
		0.0170279044465399, 0.0181602847746731, 0.0192633109144647,
		0.0203287291828132, 0.0213484399192556, 0.0223145812061882,
		0.0232196109134842, 0.0240563860659984, 0.0248182385594457,
		0.0254990462902704, 0.0260932988178836, 0.0265961567423702,
		0.0270035040565964, 0.0273119928175867, 0.0275190795769211,
		0.0276230531124511, 0.0276230531124511, 0.0275190795769211,
		0.0273119928175867, 0.0270035040565964, 0.0265961567423702,
		0.0260932988178836, 0.0254990462902704, 0.0248182385594457,
		0.0240563860659984, 0.0232196109134842, 0.0223145812061882,
		0.0213484399192556, 0.0203287291828132, 0.0192633109144647,
		0.0181602847746731, 0.0170279044465399, 0.0158744932550581,
		0.0147083601409260, 0.0135377169905568, 0.0123705982972578,
		0.0112147840891489, 0.0100777270078738, 0.00896648435932763,
		0.00788765588444296, 0.00684732791564099, 0.00585102449408378,
		0.00490366592568990, 0.00400953515140809, 0.00317225220096019,
		0.00239475689068246 };
double little_wave[len_little_wave] = { -0.0190699938545502,
		-0.0412508126064754, -0.0674710839573549, -0.0942379809962896,
		-0.117436184372004, -0.133201262133710, 0.870574412894870,
		-0.133201262133710, -0.117436184372004, -0.0942379809962896,
		-0.0674710839573549, -0.0412508126064754, -0.0190699938545502 };

struct PulseParaList {
	double dPara[USED_PARA_COUNT];
};

//extern "C" {

//jstring strToJstring(JNIEnv* env, const char* pStr)
//{
//    int        strLen    = strlen(pStr);
//    jclass     jstrObj   = (*env)->FindClass(env, "java/lang/String");
//    jmethodID  methodId  = (*env)->GetMethodID(env, jstrObj, "<init>", "([BLjava/lang/String;)V");
//    jbyteArray byteArray = (*env)->NewByteArray(env, strLen);
//    jstring    encode    = (*env)->NewStringUTF(env, "utf-8");
//
//    (*env)->SetByteArrayRegion(env, byteArray, 0, strLen, (jbyte*)pStr);
//
//    return (jstring)(*env)->NewObject(env, jstrObj, methodId, byteArray, encode);
//}

int trans_to_int(const char* s0, int* x) {
	static char s[len_x * 20];
	strcpy(s, s0);
	char *delim = ",";
	char *tmp = strtok(s, delim);
	int len = 0;
	for (int i = 0; i < len_x; i++) {
		x[i] = 0;
		if (tmp) {
			x[i] = atoi(tmp);
			tmp = strtok(NULL, delim);
			len++;
		}
	}

	// 消除读入补零引入的偏移
	int kd = len-1;
	for (int j = len; j < len_x; j++)
		x[j] = x[kd];

	return len;
}

int trans_to_string(double x[], char* output, int len_data) {
	strcpy(output, "");
	char tmp[50] = "";
	for (int i = 0; i < len_data; i++) {
		sprintf(tmp, "%f,", x[i]);
		strcat(output, tmp);
	}
	output[strlen(output) - 1] = (char) 0;

	return 0;
}

// 读取数据
int readData(double *x, int len_x, char *filePath) {
	FILE *pFile = fopen(filePath, "r");
	char pBuf[100];
	int i = 0;
	int j = 0;
	while (!feof(pFile)) {
		fgets(pBuf, 100, pFile);
		x[i] = atof(pBuf);
		i++;
		if (i >= len_x)
			break;
	}
	fclose(pFile);

	return 0;
}


// 写数据
int writeData(double *x, int len_x, char *filePath) {
	FILE *pFileOut = fopen(filePath, "w");
	char pBuf[100];
	for (int i = 0; i < len_x; i++) {
		fputs(gcvt(x[i], 10, pBuf), pFileOut);
		fputs("\n", pFileOut);
	}
	fclose(pFileOut);

	return 0;
}


/*
 输入：脉诊数据
 输出：信号幅度
 原理：差分连续增监测
 */
int getAmplitude(int *y) {
	int a = 0;

	if (len_x >= sample_rate * 2) {
		int dx[len_x];
		dx[0] = 0;
		for (int i = 1; i < len_x; i++) {
			dx[i] = y[i] - y[i - 1];
		}

		int pulse_max = 0;
		int s = 0;
		int e = 0;
		int k = 0;
		while (k < len_x) {
			int pulse = 0;
			int sp = k;
			while (dx[k] > 0) {
				pulse += dx[k];
				k++;
				if (k >= len_x)
					break;
			}
			int ep = k;
			if (pulse_max < pulse) {
				pulse_max = pulse;
				s = sp;
				e = ep;
			}
			if (pulse == 0)
				k++;
		}
		a = pulse_max;
	}

	return a;
}

// 滤波器
void lowPassFilter(double x[], double *y, double a[], int k, int len_x) {
	double s = 0;
	for (int j = 0; j < k; j++) {
		s += a[j];
	}
	if (s > 0 || s < 0) {
		for (int j = 0; j < k; j++) {
			a[j] /= s;
		}
	}

	for (int i = 0; i < len_x; i++) {
		y[i] = x[i];
		if (i >= (k - 1)) {
			s = 0;
			for (int j = 0; j < k; j++) {
				s += x[i - j] * a[j];
			}
			y[i - k / 2] = s;
		}
	}
}

// 信号最小化平移
void stdMove(double *x, int len) {
	int j = 0;

	// 消除工频干扰
	int mk[len_x];
	double thr = 10;
	
	for (j = 1; j < (len-1); j++)
	{
		mk[j] = 0;
		if((x[j] - x[j-1]) > thr && (x[j+1] - x[j]) < -thr) mk[j] = 1;
		if((x[j] - x[j-1]) < -thr && (x[j+1] - x[j]) > thr) mk[j] = 1;
	}
	for (j = 0; j < len; j++)
	{
		if(mk[j] == 1)
		{
			x[j] = (x[j-1] + x[j+1]) / 2;
		}
	}
	x[0] = x[1];
	x[len-1] = x[len-2];

	// 尺度变换
	double a = 10.473;
	double b = 630.38;
	/*
	double a = 4.0081;
	double b = 567.736;*/
	for (j = 0; j < len; j++)
		x[j] = (x[j] - b) / a;
}

// 信号标准化平移并抑制两头畸变
void sigStdz(double *y, double y2[], int len_x, int len_lp_fc) {
	for (int i = 0; i < len_x; i++)
		y[i] = y[i] - y2[i];

	int r_win = len_lp_fc / 2;
	for (int j = 0; j < r_win; j++) {
		y[j] = y[r_win];
		y[len_x - 1 - j] = y[len_x - r_win];
	}
}

// 计算差分特征
void diff(double y[], double *dy2, double *dy, int len_x) {
	dy2[0] = 0;
	for (int i = 1; i < len_x; i++)
		dy2[i] = (y[i] - y[i - 1]) * (y[i] - y[i - 1]);

	double a[3] = { 1, 1, 1 };
	lowPassFilter(dy2, dy, a, 3, len_x);
}

// 计算潮波特征
void getFeature(double y[], double *z1, double little_wave[],
		int len_little_wave, int len_x, double *dyz) {
	lowPassFilter(y, z1, little_wave, len_little_wave, len_x);
	int i = 0;
	for (i = 1; i < len_x; i++)
		z1[i] = -z1[i] / 10;

	dyz[0] = 0;
	for (i = 1; i < len_x; i++)
		dyz[i] = (z1[i] - z1[i - 1]) - (y[i] - y[i - 1]);
}

// 探测波峰波谷
void detecPeck(double y[], double dy[], double *p_mark, int len) {
	// 寻找特征上的极大值
	const int max_peak_count = 500;
	double peak_dy[max_peak_count];
	int p[max_peak_count];

	int s = 0;
	int i = 2;
	int j = 0;

	for (i = 0; i < max_peak_count; i++) {
		peak_dy[i] = 0;
		p[i] = 0;
	}
	for (i = 2; i < (len - 2); i++) {
		if (dy[i] > dy[i - 2] && dy[i] >= dy[i - 1] && dy[i] >= dy[i + 1]
				&& dy[i] > dy[i + 2] && y[i] < y[i + 1] && y[i - 1] < y[i]) {
			p_mark[i] = 1;
			peak_dy[s] = dy[i];
			s += 1;
			if (s >= max_peak_count)
				break;
		}
	}
	if (s > 1) {
		double thr_dy = 0;
		for (i = 0; i < s; i++) {
			int k = i;
			for (j = i + 1; j < s; j++)
				if (peak_dy[j] < peak_dy[k])
					k = j;

			double tmp = peak_dy[i];
			peak_dy[i] = peak_dy[k];
			peak_dy[k] = tmp;
		}
		for (i = 1; i < s; i++)
			if (thr_dy < (peak_dy[i] - peak_dy[i - 1]))
				thr_dy = peak_dy[i] - peak_dy[i - 1];

		j = 0;
		for (i = 0; i < len; i++) {
			if (p_mark[i] == 1 && dy[i] <= thr_dy)
				p_mark[i] = 0;
			if (p_mark[i] == 1) {
				p[j] = i;
				j++;
			}
		}
		s = j;

		for (i = 0; i < s; i++) {
			int p_bottom = p[i];
			int p_peak = p[i];
			for (j = p[i]; j > 0; j--) {
				if (y[j - 1] >= y[j]) {
					p_bottom = j;
					break;
				}
			}
			for (j = p[i];
					j < (i < (s - 1) ? ((p[i + 1] + p[i]) / 2) : (len - 1));
					j++) {
				if (y[p_peak] < y[j] && y[j] >= y[j - 1] && y[j] >= y[j + 1])
					p_peak = j;
			}

			if (p_bottom != p[i] && p_peak != p[i] && (p_peak - p_bottom) >= 3
					&& (y[p_peak] - y[p_bottom]) >= 1) {
				p_mark[p_bottom] = 2;
				p_mark[p_peak] = 3;
			}
		}
	}
}

// 探测潮波
void detecTideWave(double y[], double dyz[], double *p_mark, int len) {
	int i, k, m, h;

	for (i = 0; i < len; i++) {
		if (p_mark[i] != 3)
			continue;

		int p_peak = i;
		int pre_flg = 0;
		for (m = p_peak; m > 0; m--) {
			if (p_mark[m] == 1) {
				for (h = m; h < p_peak; h++) {
					if ((y[h] > y[h - 1] && y[h] >= y[h + 1])
							|| ((y[h] - y[h - 1]) < (y[h - 1] - y[h - 2])
									&& (y[h] - y[h - 1]) <= (y[h + 1] - y[h]))) {
						pre_flg = 1;
						p_mark[h] = 4;
						break;
					}
				}
				break;
			}
		}
		if (pre_flg == 1)
			continue;

		int p_zero = p_peak + 1;
		for (k = (p_peak + 1); k < len; k++) {
			if (dyz[k] > 0) {
				p_zero = k;
				break;
			}
		}
		int p_zk = p_zero;
		for (k = p_zero; k < len; k++) {
			if (dyz[k] > dyz[k - 1] && dyz[k] >= dyz[k + 1]) {
				p_zk = k;
				break;
			}
		}
		int p_tide = (((p_zk - p_zero) <= 3) ? p_zk : p_zero);

		for (k = p_tide; k > p_peak; k--) {
			if (y[k] > y[k - 1] && y[k] >= y[k + 1]) {
				p_tide = k;
				break;
			}
			if ((y[k] - y[k - 1]) > (y[k - 1] - y[k - 2])
					&& (y[k] - y[k - 1]) >= (y[k + 1] - y[k])) {
				p_tide = k;
				break;
			}
		}

		// 增加潮波过低判断
		if (p_tide > p_peak) {
			int p_thr = ((p_tide > sample_rate) ? (p_tide - sample_rate) : 0);
			for (k = p_tide; k > p_thr; k--) {
				if (p_mark[k] == 1) {
					p_thr = k;
					break;
				}
			}
			if (y[p_tide] < y[p_thr]) {
				for (k = p_tide; k > p_peak; k--) {
					if (y[k] > y[p_thr]) {
						p_tide = k;
						break;
					}
				}
			}
		}

		p_mark[p_tide] = 4;
	}
}

//重搏波
void detecWave(double y[], double z[], double *p_mark, int len) {
	int p_tide, p_bottom2, flg;
	int i, j, k, v;
	i = 0;
	while (i < len) {
		if (p_mark[i] == 4) {
			flg = 0;
			p_tide = i;
			for (j = i; j < len; j++) {
				if (p_mark[j] == 2) {
					p_bottom2 = j;
					flg = 1;
					break;
				}
			}
			if (flg > 0) {
				// 寻找极小波谷
				int pv = p_tide + 1;
				for (v = (p_tide + 1); v < (p_bottom2 - 2); v++) {
					if (z[v] < z[v - 1] && z[v] <= z[v + 1] && z[v] <= z[pv])
						pv = v;
				}
				int pk = pv + 1;
				for (k = pv; k < (p_bottom2 - 2); k++) {
					if (z[k] > z[k - 1] && z[k] >= z[k + 1]) {
						pk = k;
						break;
					}
				}
				// 寻找次谷
				int pv2 = pk + 1;
				for (v = pv2; v < (p_bottom2 - 2); v++) {
					if (z[v] < z[v - 1] && z[v] <= z[v + 1] && z[v] <= z[pv2])
						pv2 = v;
				}
				int pk2 = pv2 + 1;
				for (k = pv2; k < p_bottom2; k++) {
					if (z[k] > z[k - 1] && z[k] >= z[k + 1]) {
						pk2 = k;
						break;
					}
				}

				// 极小次小波谷相邻
				int vv_adjacent = 1;
				for (v = (pv + 1); v <= (pv2 - 1); v++) {
					if (z[v] < z[v - 1] && z[v] <= z[v + 1]) {
						vv_adjacent = 0;
						break;
					}
				}
				// 比较次谷上升量，校正极小波谷
				double df = z[pk] - z[pv];
				double df2 = z[pk2] - z[pv2];
				if (df2 > df
						&& ((df2 - df) > (z[pv2] - z[pv]) || vv_adjacent == 1))
					pv = pv2;

				// 寻找波峰
				pk = pv;
				for (k = pv; k < p_bottom2; k++) {
					if (z[k] > z[k - 1] && z[k] >= z[k + 1]) {
						pk = k;
						// 二次探测避免小波扰动
						pv2 = pk + 1;
						pk2 = pv2 + 1;
						for (int k2 = pk2; k2 < p_bottom2; k2++) {
							if (z[k2] < z[k2 - 1] && z[k2] <= z[k2 + 1])
								pv2 = k2;
							if (z[k2] > z[k2 - 1] && z[k2] >= z[k2 + 1]) {
								pk2 = k2;
								break;
							}
						}
						if ((z[pk2] - z[pv2]) > (z[pk] - z[pv2]) * 2)
							pk = pk2;
						break;
					}
				}

				for (k = pv; k <= pk; k++) {
					if (y[k] < y[k - 1] && y[k] <= y[k + 1]) {
						pv = k;
						break;
					}
				}
				int out_pk = pk;
				for (k = pv; k <= pk; k++) {
					if (y[k] > y[k - 1] && y[k] >= y[k + 1])
						out_pk = k;
				}
				p_mark[pv] = 5;
				p_mark[out_pk] = 6;

				i = pk;
			}
		}
		i++;
	}
}

// 计算参数
int calcParaList(double y[], int p_list[], double *used_para_list) {
	int P0, P1, P2, P3, P4, P5, P6;
	P0 = p_list[0];
	P1 = p_list[1];
	P2 = p_list[2];
	P3 = p_list[3];
	P4 = p_list[4];
	P5 = p_list[5];
	P6 = p_list[6];
	double SamplingRate = (double) sample_rate;
	double h1 = y[P1] - y[P0];
	double h2 = y[P2] - y[P0];
	double h3 = y[P3] - y[P0];
	double h4 = y[P4] - y[P0];
	double h5 = y[P5] - y[P4];
	double h6 = y[P6] - y[P0];
	double t1 = (double) (P1 - P0) / SamplingRate;
	double t2 = (double) (P2 - P0) / SamplingRate;
	double t3 = (double) (P3 - P0) / SamplingRate;
	double t4 = (double) (P4 - P0) / SamplingRate;
	double t5 = (double) (P6 - P4) / SamplingRate;
	double t6 = (double) (P6 - P5) / SamplingRate;
	double t = (double) (P6 - P0) / SamplingRate;

	// 尺度转换:y(3.7mm=1g), x(25mm=1s)
	double uh = 3.7;
	double ut = 25;
	double alpha = atan(h1 * uh / (t1 * ut)) / (pi / 2) * 90;
	double theta = 0;
	double W = 0;
	double W1 = 0;
	double W2 = 0;
	double delta_h = h1 * uh / 3;
	double delta_h2 = h1 * uh / 5;
	for (int k1 = P1; k1 >= P0; k1--) {
		if ((y[P1] - y[k1]) * uh > delta_h) {
			for (int k2 = P1; k2 <= P6; k2++) {
				if ((y[P1] - y[k2]) * uh > delta_h2 && W2 == 0)
					W2 = (double) (k2 - k1) / SamplingRate;
				if ((y[P1] - y[k2]) * uh > delta_h) {
					W = (double) (k2 - k1) / SamplingRate;
					W1 = W;
					theta = atan(
							(double) (P1 - k1) / SamplingRate * ut / delta_h)
							+ atan(
									(double) (k2 - P1) / SamplingRate * ut
											/ delta_h);
					theta = theta / (pi / 2) * 90;
					break;
				}
			}
			break;
		}
	}

	double HR = 60 / (t > 0 ? t : 1);
	double A = 0;
	double As = 0;
	double Ad = 0;
	double dh = 0;
	int i = 0;
	for (i = P0; i < P6; i++) {
		dh = (y[i] - y[P0]) * uh / SamplingRate * ut;
		if (i < P4)
			As += dh;
		else
			Ad += dh;
	}
	A = As + Ad;

	double para_list[USED_PARA_COUNT] = { h1, h2, h3, h4, h5, t1, t2, t3, t4,
			t5, W1, W2, As, Ad, t, alpha, theta, HR };
	/*double para_list[USED_PARA_COUNT] = {h1, h2, h3, h4, h5, h6, t1, t2, t3, t4, t5, t6, t,
	 alpha, theta, W, h2/h1, h3/h1, h4/h1, h5/h1, h6/h1, (h1-h3)/h1,
	 t1/t, t1/t2, t6/t2, t1/t4, t5/t4, W/t, (t4-t1)/t, (t2-t1)/t, HR, A, As, Ad, W1, W2}; */
	for (i = 0; i < USED_PARA_COUNT; i++)
		used_para_list[i] = para_list[i];

	return 0;
}

// 整体判断脉律类型
int judgePluseFrequencyType(int p_mark[]) {
	int plus_pos[MAX_PLUSE_COUNT]; // 脉搏位置
	int plus_prebeat_mark[MAX_PLUSE_COUNT]; // 早搏标记: 0,正常; 1,早搏后
	int t[MAX_PLUSE_COUNT]; // 脉搏时间间隔
	int j = 0;
	int count_peak = 0;
	int i = 0;
	for (i = 0; i < MAX_PLUSE_COUNT; i++) {
		plus_pos[i] = 0;
		t[i] = 0;
		plus_prebeat_mark[i] = 0;
	}
	for (i = 0; i < len_x; i++) {
		if (p_mark[i] == 1 && j < MAX_PLUSE_COUNT) {
			plus_pos[j] = i;
			if (j > 0)
				t[j] = plus_pos[j] - plus_pos[j - 1];
			j++;
		}
	}
	t[0] = t[1];
	count_peak = j;

	// 脉律类型
	int pluse_frequency_type = 0;
	if (count_peak >= 5) {
		// 脉动周期不等，t值之差＞0.12s  t ≤ 0.66s
		float Fs = (float) sample_rate;
		int df_thr = (int) (Fs * 0.12);
		int fast_thr = (int) (Fs * 0.66);
		int df_count = 0;
		int flg = 0;
		for (i = 1; i < count_peak; i++) {
			if ((t[i] - t[i - 1]) > df_thr && t[i - 1] <= fast_thr) {
				df_count++;
				plus_prebeat_mark[i] = 1;
			}
		}
		if (df_count >= 3) {
			pluse_frequency_type = dai_mai;
			int df2_count = 0;
			int df3_count = 0;
			for (i = count_peak - 1; i > 2; i--) {
				if (plus_prebeat_mark[i] == 1 && plus_prebeat_mark[i - 2] == 1)
					df2_count++;
				if (plus_prebeat_mark[i] == 1 && plus_prebeat_mark[i - 3] == 1)
					df3_count++;
			}
			if (df2_count >= 3 || df3_count >= 3)
				pluse_frequency_type = dai_mai; // 代脉：脉搏二联脉/三联脉
			else
				pluse_frequency_type = jie_mai; // 结脉：脉图中有不规则的停搏
		}
	}

	return pluse_frequency_type;
}

// 计算整体参数
PulseParaList calcPluseSignleParaList(double y[], int p_mark[], int len_x) {
	int p_bottom[MAX_PLUSE_COUNT];
	int P0, P1, P2, P3, P4, P5, P6;
	int j = 0;
	int count_peak = 0;
	int i = 0;
	for (i = 0; i < len_x; i++) {
		if (p_mark[i] == 2 && j < MAX_PLUSE_COUNT) {
			p_bottom[j] = i;
			j++;
		}
	}
	count_peak = j;
	double used_para_list[USED_PARA_COUNT], out_para_list[USED_PARA_COUNT];
	for (i = 0; i < USED_PARA_COUNT; i++)
		used_para_list[i] = out_para_list[i] = 0;
	
	int last_p6 = (p_bottom[count_peak - 1] + sample_rate);
	if(last_p6 > len_x) last_p6 = (len_x-1);
	for (j = 0; j < count_peak; j++) {
		P0 = p_bottom[j];
		P6 = ((j + 1) < count_peak) ? p_bottom[j + 1] : last_p6;
		int p_peak, p_tide, pvs, pks;
		for (i = P0; i < P6; i++) {
			if (p_mark[i] == 3)
				p_peak = i;
			if (p_mark[i] == 4)
				p_tide = i;
			if (p_mark[i] == 5)
				pvs = i;
			if (p_mark[i] == 6)
				pks = i;
		}
		P1 = p_peak;
		P3 = p_tide;
		P4 = pvs;
		P5 = pks;
		P2 = (p_peak + p_tide) / 2;
		int p_list[7] = { P0, P1, P2, P3, P4, P5, P6 };
		calcParaList(y, p_list, used_para_list);
		for (i = 0; i < USED_PARA_COUNT; i++)
			out_para_list[i] += used_para_list[i];
	}

	PulseParaList ppl;
	for (i = 0; i < USED_PARA_COUNT; i++)
		ppl.dPara[i] = out_para_list[i] / (double) count_peak;

	return ppl;
}


/* 局部最大值差异判断信号杂乱程度：
 thr = 0.2
 disorder_factor > thr 则杂乱
 输出：0，规整；1，杂乱
 */
int calcDisorder(double y[], int len_x) {
	int sub_zoon_width = sample_rate * 2;
	int cnt = 0;
	double amplitude_list[MAX_PLUSE_COUNT];
	double mean_amplitude = 0;
	int i;
	for (i = 0; i < MAX_PLUSE_COUNT; i++) {
		amplitude_list[i] = 0;
		int j1 = (i - 1) * sub_zoon_width;
		int j2 = i * sub_zoon_width;
		if (j1 < 0)
			j1 = 0;
		if (j2 <= len_x) {
			double m1 = y[j1];
			double m2 = y[j1];
			for (int j = j1; j < j2; j++) {
				if (m1 < y[j])
					m1 = y[j];
				if (m2 > y[j])
					m2 = y[j];
			}
			amplitude_list[i] = m1 - m2;
			mean_amplitude += amplitude_list[i];
			cnt++;
		}
	}

	double disorder_factor = 0.0;
	if (cnt > 0) {
		mean_amplitude = mean_amplitude / (double) cnt;
		if (mean_amplitude <= 0)
			mean_amplitude = 0.01;

		for (i = 0; i < cnt; i++) {
			double diff = amplitude_list[i] - mean_amplitude;
			if (diff <= 0.0)
				diff = -diff;
			disorder_factor += diff / mean_amplitude;
		}

		disorder_factor = disorder_factor / (double) cnt;
	}

	int is_disorder = 0;
	double thr = 0.2;
	if (disorder_factor > thr)
		is_disorder = 1;

	return is_disorder;
}

// 信号滤波
int signalFilter(int input_data[], double *y, int len) {
	// 数据初始化
	double y2[len_x], x[len_x];
	for (int i = 0; i < len_x; i++) {
		x[i] = (double) input_data[i];
		y[i] = y2[i] = 0;
	}

	// 滤波处理
	stdMove(x, len_x);
	lowPassFilter(x, y, dn_fc, len_dn_fc, len_x);
	lowPassFilter(x, y2, lp_fc, len_lp_fc, len_x);
	sigStdz(y, y2, len_x, len_lp_fc);

	return 0;
}

// 将脉象信息及参数列表转化为字符串
const char* parameterlistToPlusetype(double para_list[],
		int pluse_frequence_type, double P, const char* user_info) {
	// 获得参数
	double h1 = para_list[0];
	double h2 = para_list[1];
	double h3 = para_list[2];
	double h4 = para_list[3];
	double h5 = para_list[4];
	double t1 = para_list[5];
	double t2 = para_list[6];
	double t3 = para_list[7];
	double t4 = para_list[8];
	double t5 = para_list[9];
	double W1 = para_list[10];
	double W2 = para_list[11];
	double As = para_list[12];
	double Ad = para_list[13];
	double t = para_list[14];
	double alpha = para_list[15];
	double theta = para_list[16];
	double HR = para_list[17];

	// 脉搏类型
	char *pulse_type = "平";
	if (h3 / h1 >= 0.7 && W1 / t > 0.2 && h4 / h1 > 0.5 && h5 / h1 <= 0.05)
		pulse_type = "弦";
	if (W1 / t < 0.2 && (t1 >= 0.07 && t1 <= 0.09)
			&& (theta >= 17 && theta <= 22) && h4 / h1 < 0.5 && h5 > 0.57)
		pulse_type = "滑";

	// 脉位
	char *pulse_position = "中";
	if (P <= 75.0 && h3 / h1 < 0.3)
		pulse_position = "浮";
	if (P >= 150.0)
		pulse_position = "沉";

	// 脉律
	char *str_pluse_frequence = "齐";
	if (pluse_frequence_type == jie_mai || pluse_frequence_type == dai_mai) {
		str_pluse_frequence = "不齐";
	}

	// 脉力: 脉力分有力无力，当以沉候为准。无论浮取脉力如何，只要沉取无力即为虚，沉取有力即为实。
	char *pulse_power = "中";
	if (P >= 150) {
		double thr1 = 15;
		double thr2 = 30;
		if (h1 < thr1)
			pulse_power = "无力";
		if (h1 > thr2)
			pulse_power = "有力";
	}

	// 脉形: abc, ab, ac
	char *pulse_shape = "abc";
	//printf("test------------------1---------------");
	// 脉势
	char pulse_tension[50] = "中"; // {低，中，高}
	if(strstr(pulse_type,  "弦")) strcpy(pulse_tension, "高");
	//printf("test------------------2---------------");
	char pulse_fluency[50] = "中"; // {涩，中，滑}
	if(strstr(pulse_type,  "滑")) strcpy(pulse_fluency, "滑");
	if(strstr(pulse_type,  "涩")) strcpy(pulse_fluency, "涩");
	char pulse_tension1[50] = "不滑不涩";
	char pulse_fluency1[50] = "不紧不弛";
	//printf("test------------------3---------------");
	if (!strstr(pulse_tension,  "中"))  strcpy(pulse_tension1, pulse_tension);
	if (!strstr(pulse_fluency,  "中"))  strcpy(pulse_fluency1, pulse_fluency);
	//printf("test------------------4---------------");

	char pulse_tendency[50];
	strcpy(pulse_tendency,  pulse_tension1);
	strcat (pulse_tendency,  pulse_fluency1);
	//printf("test------------------5---------------");
	// 脉率
	char *pulse_rate_str = "中";
	if (t > 1.0)
		pulse_rate_str = "迟";
	if (t <= 0.66 && t > 0.5)
		pulse_rate_str = "数";
	if (t < 0.49)
		pulse_rate_str = "疾";
	char pulse_regularity[30];
	sprintf(pulse_regularity, "%s(%d times/min)", pulse_rate_str, (int) HR);

	// 脉名
	char pulse_name[20];
	sprintf(pulse_name, "脉%s%s", pulse_type, pulse_rate_str);
	//printf("test---------------4------------------");
	// 脉名,脉位,脉数,脉力,脉形,脉势,脉律
	// "脉滑数,中,齐,中,ac,低平虚,数(94times/min),"
	static char str_pluse_info[2000];
	sprintf(str_pluse_info,
			"%s,%s,%s,%s,%s,%s,%s,%s,%s,%s, %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",
			user_info, pulse_name, pulse_position, str_pluse_frequence,
			pulse_power, pulse_shape, pulse_tendency, pulse_tension, pulse_fluency, pulse_regularity, h1, h2,
			h3, h4, h5, t1, t2, t3, t4, t5, W1, W2, As, Ad, t);

	return str_pluse_info;
}




// 频域分析模块

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

void FFT_test()
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


void spectrumAmplitude(double xn[], double* fx_db, int len_data)
{
	double xreal[len_x];
	double ximag[len_x];
	int i;

	// 数据初始化
	for(i=0;i<len_x;i++)
		ximag[i] = 0;
	for(i=0;i<len_data;i++)
		xreal[i] = xn[i];
	// 如果数据长度不足，则实部镜面周期延拓
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


	FFT(xreal, ximag, len_x);

	for(i=0;i<len_x;i++)
		fx_db[i] = sqrt(xreal[i]*xreal[i] + ximag[i]*ximag[i])*2;
}


// 计算谱能比与谐波峰值
void pluseSignalPSA(double xn[], double* param, int len_data, double mean_RR)
{
	double fxd[len_x], fx_db[len_x];
	spectrumAmplitude(xn, fxd, len_data);

	int i = 0;
	double a[sample_rate/5];
	for(i=0;i<sample_rate/5;i++)
		a[i] = 1;
	lowPassFilter(fxd, fx_db, a, sample_rate/5, len_x);
	
	double thr = 10;
	double x = 0;
	double s1 = 0;
	double s2 = 0;
	for(i=0;i<len_x/2;i++)
	{
		x = (double)i/(double)len_x*(double)sample_rate;
		if(x <= thr) s1 += fxd[i];
		if(x > thr && x < thr*2) s2 += fxd[i];
	}
	if(s2 <= 0) s2 = 1;
	
	double SER = s1/s2;
	

	// 计算各谐波幅值
	double f0 = 1/mean_RR;
	double delta = f0*0.2;
	double thr1,thr2;
	for(int j=1;j<7;j++)
	{
		param[j] = 0;
		thr1 = f0*j - delta;
		thr2 = f0*j + delta;
		for(i=0;i<len_x/2;i++)
		{
			x = (double)i/(double)len_x*(double)sample_rate;
			if(x < thr1 || x > thr2) continue;
			if(param[j] < fx_db[i]) param[j] = fx_db[i];
		}
	}
	param[0] = SER;
}



// 计算交感/副交感神经活性比
double nervousActivityRatio(int p_mark[], int len_data, double* param)
{
	// RR间期->心率->插值
	double Fs = (double)sample_rate;
	double yi[len_x];
	int R_pos[len_x];
	int R_count = 0;
	int i;
	for(i=0;i<len_x;i++)
	{
		yi[i] = 0;
		R_pos[i] = 0;
	}
	for(i=0;i<len_x;i++)
	{
		if(p_mark[i] == 1)
		{
			R_pos[R_count] = i;
			R_count++;
		}
	}
	param[0] = 1;
	param[1] = 0;
	if(R_count > 1){
		// 计算离散心率
		int p0 = R_pos[1];
		int p1 = R_pos[0];
		double RR = ((double)p0 - (double)p1)/Fs;
		double mean_RR = RR;
		yi[p1] = RR;
		for(i=1;i<R_count;i++)
		{
			p0 = R_pos[i];
			p1 = R_pos[i-1];
			RR = ((double)p0 - (double)p1)/Fs;
			mean_RR += RR;
			yi[p0] = RR;

			// 线性插值
			double stp = (yi[p0] - yi[p1])/(p0-p1);
			for(int j=p1;j<=p0;j++)
				yi[j] = yi[p1] + stp*(j-p1);
		}
		mean_RR = mean_RR/R_count;

		// 计算HRV
		double HRV = 0;
		for(i=1;i<R_count;i++)
		{
			p0 = R_pos[i];
			p1 = R_pos[i-1];
			RR = ((double)p0 - (double)p1)/Fs;
			HRV += (RR-mean_RR)*(RR-mean_RR);
		}
		HRV = sqrt(HRV/(R_count-1));

		// 边界邻近点插值
		for(i=len_x/2;i>=0;i--)
			if(yi[i] == 0) yi[i] = yi[i+1];
		for(i=len_x/2+1;i<len_x;i++)
			if(yi[i] == 0) yi[i] = yi[i-1];

		// 参数输出
		param[0] = mean_RR;
		param[1] = HRV;
	}

	// dB
	double fftR[len_x];
	spectrumAmplitude(yi, fftR, len_data);
	int len = len_x/2;
	double d_min = fftR[0];
	for(i=0;i<len;i++)
		if(d_min > fftR[i]) d_min = fftR[i];
	for(i=0;i<len;i++)
		fftR[i] = log10(fftR[i]-d_min+1)*10;

	// LF/HF
	double LF = 0;
	double HF = 0;
	double x = 0;
	for(i=0;i<len;i++)
	{
		x = (double)i/(double)len_x*Fs;
		if(x > 0.04 && x <= 0.15) LF += fftR[i];
		if(x > 0.15 && x <= 0.4) HF += fftR[i];
	}
	double TF = LF + HF;
	if(TF <= 0) TF = 1;

	double LFU = LF/TF;
	double HFU = HF/TF;
	if(HFU <= 0) HFU = 1;


	return LFU/HFU;
}



const char* getInfoPSA(double xn[], int p_mark[], int len_data) {
	// 脉搏波形频域分析参数：10Hz谱能比,[1:6]谐波振幅
	double param[7];
	double d_tmp[2]={0, 0};
	double LF_HF_Ratio = nervousActivityRatio(p_mark, len_data, d_tmp);
	double mean_RR = d_tmp[0];
	double HRV = d_tmp[1];
	pluseSignalPSA(xn, param, len_data, mean_RR);

	static char result[len_x * 20];
	sprintf(result, "%f,%f,%f,%f,%f,%f,%f,%f,%f",
		HRV,LF_HF_Ratio,param[0],param[1],param[2],param[3],param[4],param[5],param[6]);
	printf("%s\n", result);

	return result;
}



/*
 功能：计算参数并判断波形
 输入：
 input_data: 原始数据;
 len_x:数据长度;
 P: 脉诊仪传感器压力(浮脉75g;沉脉150g;中取100g)
 输出：
 脉象信息如下：
 脉名,脉位,脉数,脉力,脉形,脉势,脉律,H1,H2,H3,H4,H5,T1,T2,T3,T4,T5,W1,W2,ASS,AD,T
 脉滑数,中,齐,中,ac,低平虚,数(94times/min),19.085,9.042,7.538,5.408,2.039,0.104,0.19,0.218,0.278,0.368,0.094,0.068,68.873,30.702,0.647
 */
const char* judgeWave(const char* userinfo, int input_data[], double P, int len_data) {
	static char user_info[len_x*20];
	strcpy(user_info, userinfo);
	// 数据初始化
	double xn[len_x], y[len_x], dy[len_x], dy2[len_x], z1[len_x], dyz[len_x], double_mark[len_x];
	int i;
	for (i = 0; i < len_x; i++) {
		xn[i] = y[i] = (double) input_data[i];
		double_mark[i] = dy[i] = dy2[i] = z1[i] = dyz[i] = 0;
	}
	// 滤波
	signalFilter(input_data, y, len_x);

	// 判断信号是否杂乱，如果杂乱则提示杂乱并退出
	int is_disorder = calcDisorder(y, len_x);
	if (is_disorder == 1) {
		const char *bad_signal_message = "信号太乱!";
		//return bad_signal_message;
	}

	// 特征计算
	diff(y, dy2, dy, len_x); // 计算dy
	getFeature(y, z1, little_wave, len_little_wave, len_x, dyz); // 计算z1,dyz

	//波形识别
	detecPeck(y, dy, double_mark, len_x);
	detecTideWave(y, dyz, double_mark, len_x);
	detecWave(y, z1, double_mark, len_x);

	int p_mark[len_x];
	for (i = 0; i < len_x; i++)
		p_mark[i] = (int) double_mark[i];
	
	PulseParaList ppl = calcPluseSignleParaList(y, p_mark, len_x);
	double out_para_list[USED_PARA_COUNT];
	for (i = 0; i < USED_PARA_COUNT; i++){
		out_para_list[i] = ppl.dPara[i];
	}
	
	int pluse_frequence_type = judgePluseFrequencyType(p_mark);
	
	const char* PSA_info = getInfoPSA(xn, p_mark, len_data);

	const char* str_pluse_type = parameterlistToPlusetype(out_para_list,
			pluse_frequence_type, P, user_info);

	return str_pluse_type;
}


const char* test1(const char *s) {
	const char *user_info = "aaa";
	const char *pressure_info = "100";
	if (!s || strlen(s) <= 0)
		return "脉诊信号太短，或不合法";
	if (!user_info || strlen(user_info) <= 0)
		return "用户信息太短，或不合法";
	if (!pressure_info || strlen(pressure_info) <= 0)
		return "未输入脉诊压力，或不合法";

	double P = 80;
	P = atof(pressure_info);
	if (P <= 0)
		return "未输入脉诊压力，或不合法";

	char sDoub[256]={0};
	int x[len_x];
	int len_data = trans_to_int(s, x);
	const char* output_str = judgeWave(user_info, x, P, len_data);

	return output_str;
}


const char* test2(const char *s) {
	int x_src[len_x];
	int len_data = trans_to_int(s, x_src);
	// 数据初始化
	double x[len_x], y[len_x], y2[len_x];
	for (int i = 0; i < len_x; i++) {
		x[i] = (double) x_src[i];
		y[i] = y2[i] = 0;
	}

	// 滤波处理
	stdMove(x, len_x);
	lowPassFilter(x, y, dn_fc, len_dn_fc, len_x);
	lowPassFilter(x, y2, lp_fc, len_lp_fc, len_x);
	sigStdz(y, y2, len_x, len_lp_fc);
	static char result[len_x * 20];
	trans_to_string(y, result, len_data);
	int strLen = strlen(result);

	return result;
}




int main(int argc, char* argv[])
{
	const int len_str = 300;
	char file_name[len_str];
	char out_file[len_str];
	FILE *pFileIni = fopen("d:/ini.txt", "r");
	fgets(file_name, len_str, pFileIni);
	fgets(out_file, len_str, pFileIni);
	fclose(pFileIni);
	file_name[strlen(file_name)-1] = '\0';
	out_file[strlen(out_file)-1] = '\0';

	/*
	char file_name[len_str] = "D:/c_test_data/peak_position_test/stf_6.csv";
	char out_file[len_str] = "D:/out_info.csv";
	*/
	char s[10000];
	FILE *pFile = fopen(file_name, "r");
	fgets(s, 10000, pFile);
	fclose(pFile);

	// "文件名, 用户信息, 脉名, 脉位, 脉律, 脉力, 脉形, 脉势, 紧张度, 流利度, 脉率, H1, H2, H3, H4, H5, T1, T2, T3, T4, T5, W1, W2, AS, AD, T",
	// %s\n
	char output_info[10000];
	sprintf(
		output_info,
		"%s,%s\n",
		file_name,
		test1(s)
	);

	FILE *pFileOut = fopen(out_file, "a");
	fputs(output_info, pFileOut);
	fclose(pFileOut);


	return 0;
}

