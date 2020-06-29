#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"

void MPI_filterNew_Xminus(Domain *D,double ***f1,double ***f2,int ny,int share,int m);
void MPI_filterNew_Xplus(Domain *D,double ***f1,double ***f2,int ny,int share,int m);

void filter(Domain *D,double ***dataR,double ***dataI)
{
    int i,j,m,istart,iend,jstart,jend,numMode,iter;
	 int nxSub,nySub,n,ii,jj; 
	 double ***valR,***valI;
    double sumR,sumI,wz[3],wr[3];
    int myrank,nTasks;  

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;
	 iter=D->filterIter;
	 nxSub=D->nxSub+5;
	 nySub=D->nySub+5;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    valR=(double ***)malloc(1*sizeof(double ** ));
    valI=(double ***)malloc(1*sizeof(double ** ));
    for(m=0; m<1; m++) {
      valR[m]=(double **)malloc(nxSub*sizeof(double * ));
      valI[m]=(double **)malloc(nxSub*sizeof(double * ));
      for(i=0; i<nxSub; i++) {
        valR[m][i]=(double *)malloc(nySub*sizeof(double  ));
        valI[m][i]=(double *)malloc(nySub*sizeof(double  ));
		}
    }
    for(i=0; i<nxSub; i++)
      for(j=0; j<nySub; j++) {
         valR[0][i][j]=0.0;
         valI[0][i][j]=0.0;
		}

    wz[0]=0.25; wz[1]=0.5; wz[2]=0.25;

    m=0;
      n=0;
      while(n<iter)  {
        //first
        j=jstart;
        wr[0]=0.5; wr[1]=0.5;
        for(i=istart; i<iend; i++) {
          sumR=0.0;
          for(ii=0; ii<3; ii++)
            for(jj=0; jj<2; jj++)
              sumR+=wz[ii]*wr[jj]*dataR[m][i-1+ii][j+jj];
          valR[0][i][j]=sumR;
        }

        wr[0]=0.25; wr[1]=0.5; wr[2]=0.25;
        for(i=istart; i<iend; i++)
          for(j=jstart+1; j<jend; j++) {
            sumR=0.0;
            for(ii=0; ii<3; ii++)
              for(jj=0; jj<3; jj++)
                sumR+=wz[ii]*wr[jj]*dataR[m][i-1+ii][j-1+jj];
            valR[0][i][j]=sumR;
          }

        if(D->L>1)  {
          MPI_filterNew_Xplus(D,valR,valI,nySub,3,0);
          MPI_filterNew_Xminus(D,valR,valI,nySub,3,0);
        } else ;

        //second
        j=jstart;
        wr[0]=0.5; wr[1]=0.5;
        for(i=istart; i<iend; i++) {
          sumR=0.0;
          for(ii=0; ii<3; ii++)
            for(jj=0; jj<2; jj++)
              sumR+=wz[ii]*wr[jj]*valR[0][i-1+ii][j+jj];
          dataR[m][i][j]=sumR;
        }

        wr[0]=0.25; wr[1]=0.5; wr[2]=0.25;
        for(i=istart; i<iend; i++)
          for(j=jstart+1; j<jend; j++) {
            sumR=0.0;
            for(ii=0; ii<3; ii++)
              for(jj=0; jj<3; jj++)
                sumR+=wz[ii]*wr[jj]*valR[0][i-1+ii][j-1+jj];
            dataR[m][i][j]=sumR;
          }

        if(D->L>1)  {
          MPI_filterNew_Xplus(D,dataR,dataI,nySub,3,m);
          MPI_filterNew_Xminus(D,dataR,dataI,nySub,3,m);
        } else ;
		  
        n++;
      }

    for(m=1; m<numMode; m++) {
      n=0;
      while(n<iter)  {
        //first
        j=jstart;
        wr[0]=0.5; wr[1]=0.0;
        for(i=istart; i<iend; i++) {
          sumR=0.0;
          sumI=0.0;
          for(ii=0; ii<3; ii++)
            for(jj=0; jj<2; jj++) {
              sumR+=wz[ii]*wr[jj]*dataR[m][i-1+ii][j+jj];
              sumI+=wz[ii]*wr[jj]*dataI[m][i-1+ii][j+jj];
				}
          valR[0][i][j]=sumR;
          valI[0][i][j]=sumI;
        }

        wr[0]=0.25; wr[1]=0.5; wr[2]=0.25;
        for(i=istart; i<iend; i++)
          for(j=jstart+1; j<jend; j++) {
            sumR=0.0;
            sumI=0.0;
            for(ii=0; ii<3; ii++)
              for(jj=0; jj<3; jj++) {
                sumR+=wz[ii]*wr[jj]*dataR[m][i-1+ii][j-1+jj];
                sumI+=wz[ii]*wr[jj]*dataI[m][i-1+ii][j-1+jj];
				  }
            valR[0][i][j]=sumR;
            valI[0][i][j]=sumI;
          }

        if(D->L>1)  {
          MPI_filterNew_Xplus(D,valR,valI,nySub,3,0);
          MPI_filterNew_Xminus(D,valR,valI,nySub,3,0);
        } else ;

        //second
        j=jstart;
        wr[0]=0.5; wr[1]=0.0;
        for(i=istart; i<iend; i++) {
          sumR=0.0;
          sumI=0.0;
          for(ii=0; ii<3; ii++)
            for(jj=0; jj<2; jj++) {
              sumR+=wz[ii]*wr[jj]*valR[0][i-1+ii][j+jj];
              sumI+=wz[ii]*wr[jj]*valI[0][i-1+ii][j+jj];
				}
          dataR[m][i][j]=sumR;
          dataI[m][i][j]=sumI;
        }

        wr[0]=0.25; wr[1]=0.5; wr[2]=0.25;
        for(i=istart; i<iend; i++)
          for(j=jstart+1; j<jend; j++) {
            sumR=0.0;
            sumI=0.0;
            for(ii=0; ii<3; ii++)
              for(jj=0; jj<3; jj++) {
                sumR+=wz[ii]*wr[jj]*valR[0][i-1+ii][j-1+jj];
                sumI+=wz[ii]*wr[jj]*valI[0][i-1+ii][j-1+jj];
				  }
            dataR[m][i][j]=sumR;
            dataI[m][i][j]=sumI;
          }

        if(D->L>1)  {
          MPI_filterNew_Xplus(D,dataR,dataI,nySub,3,m);
          MPI_filterNew_Xminus(D,dataR,dataI,nySub,3,m);
        } else ;

        n++;
		}
    }
	
    for(i=0; i<nxSub; i++) { 
      free(valR[0][i]);
      free(valI[0][i]);
    }
    free(valR[0]); free(valR); 
    free(valI[0]); free(valI); 
}

void MPI_filterNew_Xminus(Domain *D,double ***f1,double ***f2,int ny,int share,int m)
{
    int i,j,num,start,end,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank;
    double *data;

    MPI_Status status;

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    rank=myrank/D->M;

    num=2*ny*3;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0;
    for(i=0; i<share; i++) {
      for(j=0; j<ny; j++)  data[start+j]=f1[m][i+istart][j];  start+=ny;
      for(j=0; j<ny; j++)  data[start+j]=f2[m][i+istart][j];  start+=ny;
	 }

    if(rank%2==0 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);
      start=0;
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==1)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0;
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][i+istart][j]; start+=ny;
      }

    if(rank%2==1 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);
      start=0;
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_filterNew_Xplus(Domain *D,double ***f1,double ***f2,int ny,int share,int m)
{
    int i,j,num,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank,start;
    double *data;

    MPI_Status status;

    istart=D->istart;   iend=D->iend;
    jstart=D->jstart;   jend=D->jend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    rank=myrank/D->M;

    num=2*ny*2;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0;
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
      }

    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);
      start=0;
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==0 && rank!=D->L-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0;
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
      }

    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);
      start=0;
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}


