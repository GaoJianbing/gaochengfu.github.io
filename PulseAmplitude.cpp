#include "stdafx.h"
#include <stdio.h>  
#include <stdlib.h>
#include <math.h>
#include"io.h"
#include <string.h>

/*
脉搏类型参数化定义：
浮脉：P≤75g h3/h1＜0.3
沉脉：P≥150g
迟脉：t＞1s
数脉：0.66s＞t＞0.5s
疾脉：t＜0.49s
结脉：脉动周期不等，t值之差＞0.12s  t ≤ 0.66s 脉图中有不规则的停搏
代脉：脉动周期不等，t值之差＞0.12s  t ≤ 0.66s 脉搏间有规律性停搏，歇止呈1:1比例（二联脉）、2:1比例（三联脉）等多种
弦脉：h3/h1≥0.7，W/t＞0.2，h4/h1＞0.5，h5/h1≤0.05?
滑脉：W/t＜0.2，t1：0.07～0.09s；主波夹角(θ)17～22°，h4/h1＜0.5，h5＞2mm
涩脉：升支、降支斜率小，脉图呈现低平土堡状。升支时值延长，t1可为0.09～0.16s，升支中可见顿挫。?主波峰顶圆钝，主波夹角28°～50°
*/
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

const double pi = 3.14159265;
const int len = 10; // 数据长度
const int sample_rate = 40; // 采样率

const int len_x = 1020;
const int len_dn_fc = 7;
const int len_lp_fc = 60;
const int len_little_wave = 13;
double dn_fc[len_dn_fc] = {0.0597210210218931, -0.142474780383904, 0.215573464294102, 0.734360590135818, 0.215573464294102, -0.142474780383904, 0.0597210210218931};
double lp_fc[len_lp_fc] = {0.00239475689068246, 0.00317225220096019, 0.00400953515140809, 0.00490366592568990, 0.00585102449408378, 0.00684732791564099, 0.00788765588444296, 0.00896648435932763, 0.0100777270078738, 0.0112147840891489, 0.0123705982972578, 0.0135377169905568, 0.0147083601409260, 0.0158744932550581, 0.0170279044465399, 0.0181602847746731, 0.0192633109144647, 0.0203287291828132, 0.0213484399192556, 0.0223145812061882, 0.0232196109134842, 0.0240563860659984, 0.0248182385594457, 0.0254990462902704, 0.0260932988178836, 0.0265961567423702, 0.0270035040565964, 0.0273119928175867, 0.0275190795769211, 0.0276230531124511, 0.0276230531124511, 0.0275190795769211, 0.0273119928175867, 0.0270035040565964, 0.0265961567423702, 0.0260932988178836, 0.0254990462902704, 0.0248182385594457, 0.0240563860659984, 0.0232196109134842, 0.0223145812061882, 0.0213484399192556, 0.0203287291828132, 0.0192633109144647, 0.0181602847746731, 0.0170279044465399, 0.0158744932550581, 0.0147083601409260, 0.0135377169905568, 0.0123705982972578, 0.0112147840891489, 0.0100777270078738, 0.00896648435932763, 0.00788765588444296, 0.00684732791564099, 0.00585102449408378, 0.00490366592568990, 0.00400953515140809, 0.00317225220096019, 0.00239475689068246};
double little_wave[len_little_wave] = {-0.0190699938545502, -0.0412508126064754, -0.0674710839573549, -0.0942379809962896, -0.117436184372004, -0.133201262133710, 0.870574412894870, -0.133201262133710, -0.117436184372004, -0.0942379809962896, -0.0674710839573549, -0.0412508126064754, -0.0190699938545502};

struct PulseParaList  
{
    double dPara[USED_PARA_COUNT];
};



int trans_to_int(const char* s0, int* x)
{
	static char s[len_x*20];
	strcpy(s, s0);
    char *delim = ",";
    char *tmp = strtok(s, delim);
	int len_data = 0;
	for(int i=0;i<len_x;i++)
	{
		if(tmp)
		{
			x[i] = atoi(tmp);
			tmp = strtok(NULL, delim);
			len_data++;
		}
		else
		{
			x[i] = 0;
		}
	}

	return len_data;
}


int trans_to_string(double x[], char* output, int len_data)
{
	char tmp[50] = "";
	for(int i=0;i<len_data;i++)
	{
		sprintf(tmp, "%f,", x[i]);
		strcat(output, tmp);
	}
	output[strlen(output)-1] = (char)0;

	return 0;
}


// 读取数据
int readData(double *x, int len_x, char *filePath)
{
	FILE *pFile = fopen(filePath, "r");
	char pBuf[100];
	int i = 0;
	int j = 0;
	while(!feof(pFile))
	{
		fgets(pBuf,100,pFile);
		x[i] = atof(pBuf);
		i++;
		if(i >= len_x) break;
	}
	fclose(pFile);

	return 0;
}


// 写数据
int writeData(double *x, int len_x, char *filePath)
{
	FILE *pFileOut = fopen(filePath, "w");
	char pBuf[100];
	for(int i=0;i<len_x;i++)
	{
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
int getAmplitude(int *y)
{
	int a = 0;

	if(len >= sample_rate*2)
	{
		int dx[len];
		dx[0] = 0;
		for(int i=1;i<len;i++)
		{
			dx[i] = y[i]-y[i-1];
		}

		int pulse_max = 0;
		int s = 0;
		int e = 0;
		int k = 0;
		while(k < len){
			int pulse = 0;
			int sp = k;
			while(dx[k]>0){
				pulse += dx[k];
				k++;
				if(k >= len) break;
			}
			int ep = k;
			if(pulse_max < pulse){
				pulse_max = pulse;
				s = sp;
				e = ep;
			}
			if(pulse == 0) k++;
		}
		a = pulse_max;
	}

	return a;
}


// 滤波器
void lowPassFilter(double x[], double *y, double a[], int k, int len_x)
{
	double s = 0;
	for(int j=0;j<k;j++)
	{
		s += a[j];
	}
	if(s>0 || s<0)
	{
		for(int j=0;j<k;j++)
		{
			a[j] /= s;
		}
	}
	
	
	for(int i=0;i<len_x;i++)
	{
		y[i] = x[i];
		if(i >= (k-1))
		{
			s = 0;
			for(int j=0;j<k;j++)
			{
				s += x[i-j]*a[j];
			}
			y[i-k/2] = s;
		}
	}
}


// 信号最小化平移
void stdMove(double *x, int len_x)
{
	for(int j=0;j<len_x;j++)
		x[j] = (x[j] - 567.736)/4.0081;
}


// 信号标准化平移并抑制两头畸变
void sigStdz(double *y, double y2[], int len_x, int len_lp_fc)
{
	for(int i=0;i<len_x;i++)
		y[i] = y[i] - y2[i];

	int r_win = len_lp_fc/2;
	for(int j=0;j<r_win;j++)
	{
		y[j] = y[r_win];
		y[len_x - 1 - j] = y[len_x - r_win];
	}
}


// 计算差分特征
void diff(double y[], double *dy2, double *dy, int len_x)
{
	dy2[0] = 0;
	for(int i=1;i<len_x;i++)
		dy2[i] = (y[i] - y[i-1])*(y[i] - y[i-1]);

	double a[3] = {1, 1, 1};
	lowPassFilter(dy2, dy, a, 3, len_x);
}


// 计算潮波特征
void getFeature(double y[], double *z1, double little_wave[], int len_little_wave, int len_x, double *dyz)
{
	lowPassFilter(y, z1, little_wave, len_little_wave, len_x);
	int i = 0;
	for(i=1;i<len_x;i++)
		z1[i] = -z1[i]/10;

	dyz[0] = 0;
	for(i=1;i<len_x;i++)
		dyz[i] = (z1[i] - z1[i-1]) - (y[i] - y[i-1]);
}


// 探测波峰波谷
void detecPeck(double y[], double dy[], double *p_mark, int len)
{
	// 寻找特征上的极大值
	const int max_peak_count = 500;
	double peak_dy[max_peak_count];
	int p[max_peak_count];
	
	int s = 0;
	int i = 2;
	int j = 0;

	for(i=0;i<max_peak_count;i++)
	{
		peak_dy[i] = 0;
		p[i] = 0;
	}
	for(i=2;i<(len - 2);i++)
	{
		if(dy[i] > dy[i-2] && dy[i] >= dy[i-1] && dy[i] >= dy[i+1] && dy[i] > dy[i+2] && y[i] < y[i+1] && y[i-1] < y[i])
		{
			p_mark[i] = 1;
			peak_dy[s] = dy[i];
			s += 1;
			if(s >= max_peak_count) break;
		}
	}
	if(s > 1)
	{
		double thr_dy = 0;
		for(i=0;i<s;i++)
		{
			int k = i;
			for(j=i+1;j<s;j++)
				if(peak_dy[j] < peak_dy[k]) k = j;

			double tmp = peak_dy[i];
			peak_dy[i] = peak_dy[k];
			peak_dy[k] = tmp;
		}
		for(i=1;i<s;i++)
			if(thr_dy < (peak_dy[i] - peak_dy[i - 1]))
				thr_dy = peak_dy[i] - peak_dy[i - 1];
		
		j = 0;
		for(i=0;i<len;i++)
		{
			if(p_mark[i] == 1 && dy[i] <= thr_dy) p_mark[i] = 0;
			if(p_mark[i] == 1)
			{
				p[j] = i;
				j++;
			}
		}
		s = j;

		for(i=0;i<s;i++)
		{
			int p_bottom = p[i];
			int p_peak = p[i];
			for(j=p[i];j>0;j--)
			{
				if(y[j-1] >= y[j])
				{
					p_bottom = j;
					break;
				}
			}
			for(j=p[i];j<(i<(s-1)?((p[i+1]+p[i])/2):(len - 1));j++)
			{
				if(y[p_peak] < y[j] && y[j] >= y[j-1] && y[j] >= y[j+1]) p_peak = j;
			}

			if (p_bottom != p[i] && p_peak != p[i] && (p_peak-p_bottom)>=3 && (y[p_peak] - y[p_bottom]) > 10)
			{
				p_mark[p_bottom] = 2;
				p_mark[p_peak] = 3;
			}
		}
	}
}


// 探测潮波
void detecTideWave(double y[], double dyz[], double *p_mark, int len)
{
	int i,k,m,h;

	for(i=0;i<len;i++)
	{
		if(p_mark[i] != 3) continue;

		int p_peak = i;
		int pre_flg = 0;
		for(m=p_peak;m>0;m--)
		{
			if(p_mark[m] == 1)
			{
				for(h=m;h<p_peak;h++)
				{
					if(
						(y[h]>y[h-1] && y[h]>=y[h+1])
						|| (
							(y[h]-y[h-1]) < (y[h-1]-y[h-2])
							&& (y[h]-y[h-1]) <= (y[h+1]-y[h])
						)
					)
					{
						pre_flg = 1;
						p_mark[h] =4;
						break;
					}
				}
				break;
			}
		}
		if(pre_flg == 1) continue;
    
		int p_zero = p_peak + 1;
		for(k=(p_peak + 1);k<len;k++)
		{
			if(dyz[k]>0)
			{
				p_zero = k;
				break;
			}
		}
		int p_zk = p_zero;
		for(k=p_zero;k<len;k++)
		{
			if(dyz[k]>dyz[k-1] && dyz[k]>=dyz[k+1])
			{
				p_zk = k;
				break;
			}
		}
		int p_tide = (((p_zk - p_zero) <= 3)?p_zk:p_zero);

		for(k=p_tide;k>p_peak;k--)
		{
			if(y[k]>y[k-1] && y[k]>=y[k+1])
			{
				p_tide = k;
				break;
			}
			if((y[k] - y[k-1]) > (y[k-1] - y[k-2]) && (y[k] - y[k-1]) >= (y[k+1] - y[k]))
			{
				p_tide = k;
				break;
			}
		}

		// 增加潮波过低判断
		if (p_tide > p_peak)
		{
			int p_thr = ((p_tide > sample_rate)? (p_tide - sample_rate):0);
			for(k=p_tide;k>p_thr;k--)
			{
				if(p_mark[k] == 1)
				{
					p_thr = k;
					break;
				}
			}
			if(y[p_tide] < y[p_thr])
			{
				for(k=p_tide;k>p_peak;k--)
				{
					if(y[k] > y[p_thr])
					{
						p_tide = k;
						break;
					}
				}
			}
		}

		p_mark[p_tide] =4;
	}
}


//重搏波
void detecWave(double y[], double z[], double *p_mark, int len)
{
	int p_tide, p_bottom2, flg;
	int i,j,k,v;
	i = 0;
	while(i<len)
	{
		if(p_mark[i] == 4)
		{
			flg = 0;
			p_tide = i;
			for(j=i;j<len;j++)
			{
				if(p_mark[j] == 2)
				{
					p_bottom2 = j;
					flg = 1;
					break;
				}
			}
			if(flg > 0)
			{
				// 寻找极小波谷
				int pv = p_tide+1;
				for(v=(p_tide+1);v<(p_bottom2-2);v++)
				{
					if(z[v]<z[v-1] && z[v]<=z[v+1] && z[v]<=z[pv]) pv = v;
				}
				int pk = pv+1;
				for(k=pv;k<(p_bottom2-2);k++)
				{
					if(z[k]>z[k-1] && z[k]>=z[k+1])
					{
						pk = k;
						break;
					}
				}
				// 寻找次谷
				int pv2 = pk + 1;
				for(v=pv2;v<(p_bottom2-2);v++)
				{
					if(z[v]<z[v-1] && z[v]<=z[v+1] && z[v]<=z[pv2]) pv2 = v;
				}
				int pk2 = pv2 + 1;
				for(k=pv2;k<p_bottom2;k++)
				{
					if(z[k]>z[k-1] && z[k]>=z[k+1])
					{
						pk2 = k;
						break;
					}
				}

				// 极小次小波谷相邻
				int vv_adjacent = 1;
				for(v=(pv + 1);v<=(pv2 - 1);v++)
				{
					if(z[v]<z[v-1] && z[v]<=z[v+1])
					{
					   vv_adjacent = 0;
					   break;
					}
				}
				// 比较次谷上升量，校正极小波谷
				double df = z[pk] - z[pv];
				double df2 = z[pk2] - z[pv2];
				if(df2>df && ((df2-df)>(z[pv2]-z[pv]) || vv_adjacent==1)) pv = pv2;
    
				// 寻找波峰
				pk = pv;
				for(k=pv;k<p_bottom2;k++)
				{
					if(z[k]>z[k-1] && z[k]>=z[k+1])
					{
						pk = k;
						// 二次探测避免小波扰动
						pv2 = pk + 1;
						pk2 = pv2 + 1;
						for(int k2=pk2;k2<p_bottom2;k2++)
						{
							if (z[k2]<z[k2-1] && z[k2]<=z[k2+1]) pv2 = k2;
							if (z[k2]>z[k2-1] && z[k2]>=z[k2+1])
							{
								pk2 = k2;
								break;
							}
						}
						if((z[pk2] - z[pv2]) > (z[pk] - z[pv2])*2) pk = pk2;
						break;
					}
				}

				for(k=pv;k<=pk;k++)
				{
					if(y[k] < y[k-1] && y[k] <= y[k+1])
					{
						pv = k;
						break;
					}
				}
				int out_pk = pk;
				for(k=pv;k<=pk;k++)
				{
					if(y[k] > y[k-1] && y[k] >= y[k+1]) out_pk = k;
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
int calcParaList(double y[], int p_list[], double *used_para_list)
{
	int P0, P1, P2, P3, P4, P5, P6;
	P0 = p_list[0];
	P1 = p_list[1];
	P2 = p_list[2];
	P3 = p_list[3];
	P4 = p_list[4];
	P5 = p_list[5];
	P6 = p_list[6];
	double SamplingRate = (double)sample_rate;
	double h1 = y[P1] - y[P0];
	double h2 = y[P2] - y[P0];
	double h3 = y[P3] - y[P0];
	double h4 = y[P4] - y[P0];
	double h5 = y[P5] - y[P4];
	double h6 = y[P6] - y[P0];
	double t1 = (double)(P1-P0)/SamplingRate;
	double t2 = (double)(P2-P0)/SamplingRate;
	double t3 = (double)(P3-P0)/SamplingRate;
	double t4 = (double)(P4-P0)/SamplingRate;
	double t5 = (double)(P6-P4)/SamplingRate;
	double t6 = (double)(P6-P5)/SamplingRate;
	double t = (double)(P6-P0)/SamplingRate;

	// 尺度转换:y(3.7mm=1g), x(25mm=1s)
	double uh = 3.7;
	double ut = 25;
	double alpha = atan(h1*uh/(t1*ut))/(pi/2)*90;
	double theta = 0;
	double W = 0;
	double W1 = 0;
	double W2 = 0;
	double delta_h = h1*uh/3;
	double delta_h2 = h1*uh/5;
	for(int k1=P1;k1>=P0;k1--)
	{
		if ((y[P1]-y[k1])*uh > delta_h)
		{
			for(int k2=P1;k2<=P6;k2++)
			{
				if ((y[P1]-y[k2])*uh > delta_h2 && W2 == 0)
					W2 = (double)(k2 - k1)/SamplingRate;
				if ((y[P1]-y[k2])*uh > delta_h)
				{
					W = (double)(k2 - k1)/SamplingRate;
					W1 = W;
					theta = atan((double)(P1-k1)/SamplingRate*ut/delta_h) + atan((double)(k2-P1)/SamplingRate*ut/delta_h);
					theta = theta/(pi/2)*90;
					break;
				}
			}
			break;
		}
	}

	double HR = 60/(t>0?t:1);
	double A = 0;
	double As = 0;
	double Ad = 0;
	double dh = 0;
	int i = P0;
	for(i=P0;i<P6;i++)
	{
		dh = (y[i] - y[P0])*uh/SamplingRate*ut;
		if(i < P4) As += dh;
		else Ad += dh;
	}
	A = As + Ad;

	double para_list[USED_PARA_COUNT] = {h1, h2, h3, h4, h5, t1, t2, t3, t4, t5, W1, W2, As, Ad, t, alpha, theta, HR};
	/*double para_list[USED_PARA_COUNT] = {h1, h2, h3, h4, h5, h6, t1, t2, t3, t4, t5, t6, t,
		alpha, theta, W, h2/h1, h3/h1, h4/h1, h5/h1, h6/h1, (h1-h3)/h1,
		t1/t, t1/t2, t6/t2, t1/t4, t5/t4, W/t, (t4-t1)/t, (t2-t1)/t, HR, A, As, Ad, W1, W2}; */
	for(i=0;i<USED_PARA_COUNT;i++)
		used_para_list[i] = para_list[i];

	return 0;
}


// 整体判断脉律类型
int judgePluseFrequencyType(int p_mark[])
{
	int plus_pos[MAX_PLUSE_COUNT]; // 脉搏位置
	int plus_prebeat_mark[MAX_PLUSE_COUNT]; // 早搏标记: 0,正常; 1,早搏后
	int t[MAX_PLUSE_COUNT]; // 脉搏时间间隔
	int j = 0;
	int count_peak = 0;
	int i = 0;
	for(i=0;i<MAX_PLUSE_COUNT;i++)
	{
		plus_pos[i] = 0;
		t[i] = 0;
		plus_prebeat_mark[i] = 0;
	}
	for(i=0;i<len_x;i++)
	{
		if(p_mark[i] == 1 && j < MAX_PLUSE_COUNT)
		{
			plus_pos[j] = i;
			if(j > 0) t[j] = plus_pos[j] - plus_pos[j-1];
			j++;
		}
	}
	t[0] = t[1];
	count_peak = j;

	// 脉律类型
	int pluse_frequency_type = 0;
	if(count_peak >= 5)
	{
		// 脉动周期不等，t值之差＞0.12s  t ≤ 0.66s
		float Fs = (float)sample_rate;
		int df_thr = (int)(Fs*0.12);
		int fast_thr = (int)(Fs*0.66);
		int df_count = 0;
		int flg = 0;
		for(i=1;i<count_peak;i++)
		{
			if((t[i] - t[i-1]) > df_thr && t[i-1] <= fast_thr)
			{
				df_count++;
				plus_prebeat_mark[i] = 1;
			}
		}
		if (df_count >= 3)
		{
			pluse_frequency_type = dai_mai;
			int df2_count = 0;
			int df3_count = 0;
			for(i=count_peak-1;i>2;i--)
			{
				if (plus_prebeat_mark[i] == 1 && plus_prebeat_mark[i-2] == 1) df2_count++;
				if (plus_prebeat_mark[i] == 1 && plus_prebeat_mark[i-3] == 1) df3_count++;
			}
			if (df2_count >= 3 || df3_count >= 3) pluse_frequency_type = dai_mai; // 代脉：脉搏二联脉/三联脉
			else  pluse_frequency_type = jie_mai; // 结脉：脉图中有不规则的停搏
		}
	}

	return pluse_frequency_type;
}


// 计算整体参数
PulseParaList calcPluseSignleParaList(double y[], int p_mark[], int len_x)
{
	int p_bottom[MAX_PLUSE_COUNT];
	int P0, P1, P2, P3, P4, P5, P6;
	int j = 0;
	int count_peak = 0;
	int i = 0;
	for(i=0;i<len_x;i++)
	{
		if(p_mark[i] == 2 && j < MAX_PLUSE_COUNT)
		{
			p_bottom[j] = i;
			j++;
		}
	}
	count_peak = j;
	double used_para_list[USED_PARA_COUNT], out_para_list[USED_PARA_COUNT];
	for(i=0;i<USED_PARA_COUNT;i++)
		used_para_list[i] = out_para_list[i] = 0;

	
	for(j=0;j<count_peak;j++)
	{
		P0 = p_bottom[j];
		P6 = ((j+1) < count_peak) ? p_bottom[j+1] : len_x;
		int p_peak, p_tide, pvs, pks;
		for(i=P0;i<P6;i++)
		{
			if(p_mark[i] == 3)p_peak = i;
			if(p_mark[i] == 4)p_tide = i;
			if(p_mark[i] == 5)pvs = i;
			if(p_mark[i] == 6)pks = i;
		}
		P1 = p_peak;
		P3 = p_tide;
		P4 = pvs;
		P5 = pks;
		P2 = (p_peak+p_tide)/2;
		int p_list[7] = {P0, P1, P2, P3, P4, P5, P6};
		calcParaList(y, p_list, used_para_list);
		for(i=0;i<USED_PARA_COUNT;i++)
			out_para_list[i] += used_para_list[i];
	}
	
	PulseParaList ppl;
	for(i=0;i<USED_PARA_COUNT;i++)
		ppl.dPara[i]= out_para_list[i]/(double)count_peak;

	return ppl;
}


/* 局部最大值差异判断信号杂乱程度：
	thr = 0.2
	disorder_factor > thr 则杂乱
	输出：0，规整；1，杂乱
*/
int calcDisorder(double y[], int len_x)
{
	int sub_zoon_width = sample_rate*2;
	int cnt = 0;
	double amplitude_list[MAX_PLUSE_COUNT];
	double mean_amplitude = 0;
	int i;
	for(i=0;i<MAX_PLUSE_COUNT;i++)
	{
		amplitude_list[i] = 0;
		int j1 = (i-1)*sub_zoon_width;
		int j2 = i*sub_zoon_width;
		if(j1<0) j1 = 0;
		if(j2 <= len_x)
		{
			double m1 = y[j1];
			double m2 = y[j1];
			for(int j=j1;j<j2;j++)
			{
				if(m1 < y[j]) m1 = y[j];
				if(m2 > y[j]) m2 = y[j];
			}
			amplitude_list[i] = m1 - m2;
			mean_amplitude += amplitude_list[i];
			cnt++;
		}
	}
	
	double disorder_factor = 0.0;
	if(cnt > 0)
	{
		mean_amplitude = mean_amplitude/(double)cnt;
		if(mean_amplitude <=0) mean_amplitude = 0.01;
		
		for(i=0;i<cnt;i++)
		{
			double diff = amplitude_list[i] - mean_amplitude;
			if(diff <= 0.0) diff = -diff;
			disorder_factor += diff/mean_amplitude;
		}

		disorder_factor = disorder_factor/(double)cnt;
	}

	int is_disorder = 0;
	double thr = 0.2;
	if(disorder_factor > thr) is_disorder = 1;

	return is_disorder;
}


// 信号滤波
int signalFilter(int input_data[], double *y, int len)
{
	// 数据初始化
	double y2[len_x], x[len_x];
	for(int i=0;i<len_x;i++)
	{
		x[i] = (double)input_data[i];
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
const char* parameterlistToPlusetype(double para_list[], int pluse_frequence_type, double P, const char* user_info)
{
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
	if(h3/h1>=0.7 && W1/t>0.2 && h4/h1>0.5 && h5/h1<=0.05) pulse_type = "弦";
	if(W1/t<0.2 && (t1>=0.07 && t1<=0.09) &&(theta >= 17 && theta <= 22) && h4/h1<0.5 && h5>0.57) pulse_type = "滑";
	
	// 脉位
	char *pulse_position = "中";
	if(P<=75.0 && h3/h1<0.3) pulse_position = "浮";
	if(P>=150.0) pulse_position = "沉";
	
	// 脉数
	char *str_pluse_frequence = "齐";
	if(pluse_frequence_type == jie_mai || pluse_frequence_type == dai_mai)
	{
		str_pluse_frequence = "不齐";
	}
	
	// 脉力: 脉力分有力无力，当以沉候为准。无论浮取脉力如何，只要沉取无力即为虚，沉取有力即为实。
	char *pulse_power = "unknown";
	if(P >= 150)
	{
		double thr1 = 15;
		double thr2 = 30;
		pulse_power = "中";
		if(h1 < thr1) pulse_power = "无力";
		if(h1 > thr2) pulse_power = "有力";
	}
	
	// 脉形: abc, ab, ac
	char *pulse_shape = "unknown";
	
	// 脉势: 正常，低平虚
	char *pulse_tendency = "unknown";
	if(P >= 75 && P <= 150)
	{
		double thr1 = 15;
		double thr2 = 30;
		pulse_power = "正常";
		if(h1 < thr1) pulse_power = "低平虚";
		if(h1 > thr2) pulse_power = "亢实";
	}

	// 脉律
	char *pulse_rate_str = "平";
	if(t>1.0) pulse_rate_str = "迟";
	if(t<=0.66 && t>0.5) pulse_rate_str = "数";
	if(t<0.49) pulse_rate_str = "疾";
	char pulse_regularity[30];
	sprintf(pulse_regularity, "%s(%d times/min)", pulse_rate_str, (int)HR);

	// 脉名
	char pulse_name[20];
	sprintf(pulse_name, "脉%s%s", pulse_type, pulse_rate_str);
	
	// 脉名,脉位,脉数,脉力,脉形,脉势,脉律
	// "脉滑数,中,齐,中,ac,低平虚,数(94times/min),"
	static char str_pluse_info[2000];
	sprintf(str_pluse_info, "%s,%s,%s,%s,%s,%s,%s,%s, %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",
		user_info,pulse_name, pulse_position, str_pluse_frequence, pulse_power, pulse_shape, pulse_tendency, pulse_regularity,
		h1, h2, h3, h4, h5, t1, t2, t3, t4, t5, W1, W2, As, Ad, t);

	return str_pluse_info;
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
const char* judgeWave(const char* user_info, int input_data[], double P)
{
	// 数据初始化
	double y[len_x], dy[len_x], dy2[len_x], z1[len_x], dyz[len_x], double_mark[len_x];
	int i;
	for(i=0;i<len_x;i++)
	{
		y[i] = (double)input_data[i];
		double_mark[i] = dy[i] = dy2[i] = z1[i] = dyz[i] = 0;
	}
	// 滤波
	signalFilter(input_data, y, len_x);

	// 判断信号是否杂乱，如果杂乱则提示杂乱并退出
	int is_disorder = calcDisorder(y, len_x);
	if(is_disorder == 1)
	{
		const char *bad_signal_message = "信号太乱!";
		return bad_signal_message;
	}

	// 特征计算
	diff(y, dy2, dy, len_x);// 计算dy
	getFeature(y, z1, little_wave, len_little_wave, len_x, dyz); // 计算z1,dyz
	
	//波形识别
	detecPeck(y, dy, double_mark, len_x);
	detecTideWave(y, dyz, double_mark, len_x);
	detecWave(y, z1, double_mark, len_x);

	int p_mark[len_x];
	for(i=0;i<len_x;i++)
		p_mark[i] = (int)double_mark[i];

	PulseParaList ppl = calcPluseSignleParaList(y, p_mark, len_x);
	double out_para_list[USED_PARA_COUNT];
	for(i=0;i<USED_PARA_COUNT;i++)
		out_para_list[i] = ppl.dPara[i];
	int pluse_frequence_type = judgePluseFrequencyType(p_mark);
	const char* str_pluse_type = parameterlistToPlusetype(out_para_list, pluse_frequence_type, P, user_info);

	return str_pluse_type;
}



// 滤波器测试
char* outStringSignalFilter(char *s)
{
	int x_src[len_x];
	int len_data = trans_to_int(s, x_src);

	// 数据初始化
	double x[len_x], y[len_x], y2[len_x];
	for(int i=0;i<len_x;i++)
	{
		x[i] = (double)x_src[i];
		y[i] = y2[i] = 0;
	}

	// 滤波处理
	stdMove(x, len_x);
	lowPassFilter(x, y, dn_fc, len_dn_fc, len_x);
	lowPassFilter(x, y2, lp_fc, len_lp_fc, len_x);
	sigStdz(y, y2, len_x, len_lp_fc);

	static char result[len_x*20];
	trans_to_string(y, result, len_data);

	return result;
}


int test_recognize(char *input_file, char *output_file)
{
	// 数据初始化
	double y[len_x], dy[len_x], dy2[len_x], z1[len_x], dyz[len_x], p_mark[len_x];
	
	readData(y, len_x, input_file);
	for(int i=0;i<len_x;i++)
		p_mark[i] = dy[i] = dy2[i] = z1[i] = dyz[i] = 0;

	// 特征计算
	diff(y, dy2, dy, len_x);// 计算dy
	getFeature(y, z1, little_wave, len_little_wave, len_x, dyz); // 计算z1,dyz
	
	//波形识别
	detecPeck(y, dy, p_mark, len_x);
	detecTideWave(y, dyz, p_mark, len_x);
	detecWave(y, z1, p_mark, len_x);

	// 存储结果
	writeData(p_mark, len_x, output_file);
	

	return 0;
}


const char* outCalcPluseType(char* s, char* user_info, const char* pressure_info)
{
	if(!s || strlen(s) <= 0) return "脉诊信号太短，或不合法";
	if(!user_info || strlen(user_info) <= 0) return "用户信息太短，或不合法";
	if(!pressure_info || strlen(pressure_info) <= 0) return "未输入脉诊压力，或不合法";

	double P = 80;
	P = atof(pressure_info);
	if(P <= 0) return "脉诊压力不合法";

	int x[len_x];
	int len_data = trans_to_int(s, x);

	const char* output_str = judgeWave(user_info, x, P);
	
	return output_str;
}


int test_daosheng()
{
	struct _finddata_t files;
	int File_Handle;
	int i=0;
	
	char dir_str[500] = "D:/matlab_data/pic_data/data/out/";
	char *file_name = strcat(dir_str, "*.csv");
	File_Handle = _findfirst(file_name,&files);
	if(File_Handle==-1)
	{
		printf("error\n");
		return 0;
	}
	do
	{
		sprintf(file_name, "%s%s", dir_str, files.name);
		printf("%s \n", file_name);
		i++;
	}while(0==_findnext(File_Handle,&files));
	_findclose(File_Handle);
	printf("Find %d files\n",i);

	return 0;
}


int main(int argc, char* argv[])
{
	char f1[] = "D:/c_test_data/src_0.csv";
	char f2[] = "D:/c_test_data/filter_0.csv";
	char f3[] = "D:/c_test_data/p_mark_0.csv";
	char f4[] = "D:/c_test_data/p_type_0.csv";
	for(int i=1;i<=8;i++)
	{
		f1[sizeof(f1)-6] = f2[sizeof(f2)-6] = f3[sizeof(f3)-6] = f4[sizeof(f4)-6] = (char)(i+(int)'0');
		//test_filter(f1, f2);
		//test_recognize(f2, f3);
		
		double y[len_x];
		readData(y, len_x, f1);

		static char s0[len_x*20], s[len_x*20];
		trans_to_string(y, s0, len_x);
		strcpy(s, s0);
		
		/*
		char* filter_result = outStringSignalFilter(s0);
		printf("file: %s\n filtered data: %s\n\n", f1, filter_result);
		printf("%s\n", s);
		*/
		char* user_info = "e150618092841,9,李志星,男,0,0,2015/6/18 15:47,左手";
		const char* pressure_info = "80";
		const char* output_str = outCalcPluseType(s, user_info, pressure_info);
		printf("file: %s\n pressure_info: %s\n result: %s\n\n\n\n", f1, pressure_info, output_str);

		break;
	}
	return 0;
}

