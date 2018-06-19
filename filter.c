#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"

void filterData(double ***fieldR,double ***fieldI,int numMode,int istart,int iend,int jstart,int jend);

void filterFieldC(Domain *D)
{
    int i,j,m,istart,iend,jstart,jend,numMode; 
    int myrank,nTasks;  

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

/*
    if(D->L>1)  {
      MPI_Transfer6F_Xminus(D,D->JzR,D->JrR,D->JpR,D->JzI,D->JrI,D->JpI,D->nySub+5,3);
      MPI_Transfer6F_Xplus(D,D->JzR,D->JrR,D->JpR,D->JzI,D->JrI,D->JpI,D->nySub+5,3);
    } else  ;
*/

//    filterData(D->EzCR,D->EzCI,numMode,istart,iend,jstart,jend);
//    filterData(D->PlCR,D->PlCI,numMode,istart,iend,jstart,jend);
//    filterData(D->SlCR,D->SlCI,numMode,istart,iend,jstart,jend);

}

void filterField(Domain *D)
{
    int i,j,m,istart,iend,jstart,jend,numMode; 
    int myrank,nTasks;  

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

/*
    if(D->L>1)  {
      MPI_Transfer6F_Xminus(D,D->JzR,D->JrR,D->JpR,D->JzI,D->JrI,D->JpI,D->nySub+5,3);
      MPI_Transfer6F_Xplus(D,D->JzR,D->JrR,D->JpR,D->JzI,D->JrI,D->JpI,D->nySub+5,3);
    } else  ;
*/

    filterData(D->JzR,D->JzI,numMode,istart,iend,jstart,jend);
    filterData(D->JrR,D->JrI,numMode,istart,iend,jstart,jend);
    filterData(D->JpR,D->JpI,numMode,istart,iend,jstart,jend);
//    filterData(D->PlR,D->PlI,numMode,istart,iend,jstart,jend);
//    filterData(D->SlR,D->SlI,numMode,istart,iend,jstart,jend);

}

void filterData(double ***fieldR,double ***fieldI,int numMode,int istart,int iend,int jstart,int jend)
{
  int i,j,m,a,alpha[2];
  double nowR,nowI,beforeR,beforeI;

  alpha[0]=0.5;
  alpha[1]=1.5;
         
  for(a=0; a<2; a++)
  {
    for(m=0; m<numMode; m++)
    {
      for(j=jstart+1; j<jstart+3; j++)
      {
        i=istart-1;
        nowR=fieldR[m][i][j];
        nowI=fieldI[m][i][j];
        for(i=istart; i<iend; i++)
        {
          beforeR=nowR;
          beforeI=nowI;
          fieldR[m][i][j]=alpha[a]*nowR+(1.0-alpha[a])*(beforeR+fieldR[m][i+1][j])*0.5;       
          fieldI[m][i][j]=alpha[a]*nowI+(1.0-alpha[a])*(beforeI+fieldI[m][i+1][j])*0.5;       
        }
      } 
    }
  }
}    
