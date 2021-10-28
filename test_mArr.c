#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include "src/statNS.h"

int main(int argc, char *argv[]){
	int pf;

	/*这里给定感兴趣的质量序列*/
	double massArray[]={1.25, 1.325, 1.4, 1.44 ,1.46, 1.48, 1.5, 1.55, 1.62, 1.64, 1.67, 1.7123, 1.7456, 1.8257, 1.90123, 1.95};

	double rhocArray[16], massInspect[16];
	CompactStar_t Results[16];
	EoS_t myEoS;

	loadEoS(&myEoS, "APR.txt");

	/*这个函数能以多线程+插值的方式快速确定整个质量序列对应的中心密度，序列越长平均效率越高，但不建议太长，现在内部设置了的质量-密度存储空间为2048，如果超过可能会出错*/
	M2Rhoc_Arr_fm(rhocArray, &myEoS, massArray, 16, 8);

	/*获取了密度序列rhocArray之后可以利用之前的多线程程序，快速获得对应的f-mode频率*/
	fmode_mt(Results, &myEoS, rhocArray, 16, 8);

	/*这里的做法可以有效避免同时加载多个物态方程带来的内存开销问题*/

	for(pf=0;pf<16;pf++){
		fprintf(stderr,"Rhoc=%e -->\tM=%lf freq=%lf dTime=%lf\n",rhocArray[pf],Results[pf].M,Results[pf].freq,Results[pf].dampTime);
	}

	return 0;
}
