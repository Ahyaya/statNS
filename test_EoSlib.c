/*
* 本测试示例为m2fmode_mt()使用方法，
* 用于利用多线程同时计算多个物态方程在某个质量上的fmode频率与衰减时间。
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>

#include "src/statNS.h"

int main(){

    char EoS_lib_dir[64]="EoS_lib/";

    /*判断路径是否以斜杠结尾，没有的话自动加上*/
    EoS_lib_dir[strlen(EoS_lib_dir)]=(EoS_lib_dir[strlen(EoS_lib_dir)-1]==0x2F)?0x00:0x2F;

    int EoSlistLen=0,pf;
    DIR *dirPt;
    struct dirent *fileBlock;

    if((dirPt=opendir(EoS_lib_dir))==NULL){
        fprintf(stderr,"unable to access directory\n");
        return -1;
    }
    
    Path_t EoSlist[8];/*用来存放EoS的文件路径，这里能存8个*/
    CompactStar_t Results[8];/*用来存放计算结果*/
    
    /*用你喜欢的办法将EoS文件的路径写进EoSlist[].FilePath上*/
    while(fileBlock=readdir(dirPt)){
        if(EoSlistLen<8){
            /*判断是否为current dir和parent dir标识符*/
            if(fileBlock->d_name[0]==0x2E){
                continue;
            }
            sprintf(EoSlist[EoSlistLen].FilePath,"%s%s",EoS_lib_dir,fileBlock->d_name);
            EoSlistLen++;
        }else{
            m2fmode_mt(Results,1.4,EoSlist,EoSlistLen,8);
            for(pf=0;pf<EoSlistLen;pf++){
                fprintf(stderr,"%s Rho=%.3e M=%.3lf freq=%lf dptime=%.3lf\n",EoSlist[pf].FilePath,Results[pf].Rho,Results[pf].M,Results[pf].freq,Results[pf].dampTime);
            }
            EoSlistLen=0;
        }
    }
    closedir(dirPt);

    /*check if all EoS have been computed*/
    if(EoSlistLen){
        m2fmode_mt(Results,1.4,EoSlist,EoSlistLen,8);
        for(pf=0;pf<EoSlistLen;pf++){
            fprintf(stderr,"%s Rho=%.3e M=%.3lf freq=%lf dptime=%.3lf\n",EoSlist[pf].FilePath,Results[pf].Rho,Results[pf].M,Results[pf].freq,Results[pf].dampTime);
        }
        EoSlistLen=0;
    }
    
    return 0;
}
