#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "src/statNS.h"

int main(){
    int EoSlistLen=0,pf;
    
    struct EoSpath_t EoSlist[1000];/*用来存放EoS的文件路径，这里能存1000个*/
    
    /*用你喜欢的办法将EoS文件的路径写进EoSlist[].FilePath上*/
    system("ls EoS_lib/*.txt > tmpEoS.list");
    FILE *fp=fopen("tmpEoS.list","r");
    while(fgets(EoSlist[EoSlistLen].FilePath,128,fp)){
        EoSlist[EoSlistLen].FilePath[strlen(EoSlist[EoSlistLen].FilePath)-1]=0x00;
        EoSlistLen++;
    }
    fclose(fp);

/*初始化一个数组来存放计算结果，它未必要跟EoS的总数量一样长，因为你可以分多次进行计算*/
    struct CompactStar_t Results[EoSlistLen];

/*这里指定用6个线程进行计算，结果存放到Results[]上，目标质量为1.44倍太阳质量，EoS文件名列表入口，总共要计算的EoS数量*/
    m2fmode_fm(6,Results,1.44,EoSlist,EoSlistLen);
    
/*为了能操控中间结果，我推荐每次不要算太多，这样数组Results[]可以短一些，
* 而且计算虽然是多线程，但是在没有全部完成之前你是无法查看中间任何结果的，
* 为了防止出现计算崩溃导致结果全部丢失的问题，我建议手动对EoSlist进行分割，
* 比如说，10000个物态方程，我可以先将它们的路径全部保存至EoSlist[]，此时EoSlistLen=10000，
* 我们可以通过调用50次 m2fmode_fm() 函数来达到计算目的（这样每次要计算200个）：
* struct CompactStar_t Results[200];
* int pf;
* for(pf=0;pf<50;pf++){
*     m2fmode_fm(10,Results,1.44,EoSlist+200*pf,200);
*     ...
* }
* 
*/

/*一种直接打印结果的方法，当然你可以对这个结果做其他处理*/
    for(pf=0;pf<EoSlistLen;pf++){
        fprintf(stderr,"%s Rho=%e M=%lf freq=%lf damptime=%lf\n",EoSlist[pf].FilePath,Results[pf].Rho,Results[pf].M,Results[pf].freq,Results[pf].dampTime);
    }
	
    return 0;
}
