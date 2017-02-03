#include "stdafx.h"
#include "stdio.h" 

const int len_x = 1020;

int trans_to_int(const char* s0, int* x) {
	static char s[len_x * 20];
	strcpy(s, s0);
	char *delim = ",";
	char *tmp = strtok(s, delim);
	int len_data = 0;
	for (int i = 0; i < len_x; i++) {
		x[i] = 0;
		if (tmp) {
			x[i] = atoi(tmp);
			tmp = strtok(NULL, delim);
			len_data++;
		}
	}

	return len_data;
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

// ¶ÁÈ¡Êý¾Ý
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


int main(int argc, char* argv[])
{
	char filepath[MAX_PATH]="D:/tough/crack_test/test_no_crack/resized_pic/"; 
	find(filepath);

	
	return 0;
}

