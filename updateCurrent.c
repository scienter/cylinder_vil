#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

double maximum(double x1,double x2)  {
   double result;
   if(x1>=x2)      result=x1;
   else           result=x2;
   return result;
}

double minimum(double x1,double x2)  {
   double result;
   if(x1>=x2)      result=x2;
   else            result=x1; 
   return result;
}


void updateCurrent_1st(Domain *D,int nSpecies,int iteration);
void updateCurrent_villasenor(Domain *D,int nSpecies,int iteration);
void updateCurrent_umeda(Domain *D,int nSpecies,int iteration);
void calculaionRally(double *xr,double *yr,double rr,double x1,double x2,double y1,double y2,double r1,double r2,int iteration);

void updateCurrent(Domain D,int iteration)
{
  void MPI_TransferJ_Xplus();
  void MPI_TransferJ_Xminus();
  void MPI_TransferJ_Yplus();
  void MPI_TransferJ_Yminus();

  double ***val;
  int nxSub,nySub,i,j,n,iter;
  int nSpecies;

  nxSub=D.nxSub+5;
  nySub=D.nySub+5;
  nSpecies=D.nSpecies;
  iter=D.filterIter;
  if(D.boostOn==ON && D.boostIon==OFF)
    nSpecies=1;    

  switch(D.currentType)  {
  case 1 :    
//    updateCurrent_1st(&D,nSpecies,iteration);
//    updateCurrent_villasenor(&D,nSpecies,iteration);
    updateCurrent_umeda(&D,nSpecies,iteration);

    if(D.L>1)  {
      MPI_TransferJ_Xplus(&D,D.JzR,D.JrR,D.JpR,D.JzI,D.JrI,D.JpI,D.nySub+5,3);
      MPI_TransferJ_Xminus(&D,D.JzR,D.JrR,D.JpR,D.JzI,D.JrI,D.JpI,D.nySub+5,3);
      if(D.fieldType==NDFX) {
        MPI_Transfer12F_Xplus(&D,D.JrCR,D.JrCI,D.JrR,D.JrI,D.JpCR,D.JpCI,D.JpR,D.JpI,D.JzCR,D.JzCI,D.JzR,D.JzI,D.nySub+5,3);
        MPI_Transfer12F_Xminus(&D,D.JrCR,D.JrCI,D.JrR,D.JrI,D.JpCR,D.JpCI,D.JpR,D.JpI,D.JzCR,D.JzCI,D.JzR,D.JzI,D.nySub+5,3);
      } else;
    }  else	;

    filter(&D,D.JzR,D.JzI);
    filter(&D,D.JrR,D.JrI);
    filter(&D,D.JpR,D.JpI);
    break;
  }
}



void updateCurrent_umeda(Domain *D,int nSpecies,int iteration)
{
    int i,j,m,s,i1,i2,j1,j2,ii,jj,numMode;
    int istart,iend,jstart,jend,minRSub;
    int nxSub,nySub;
    double oldZ,oldR,weight,gamma;
    double x1,x2,y1,y2,z1,z2,r1,r2,alpha;
    double Fz,Fr,factor,factM;
    double Wz[2],Wr[2],Wz1[2],Wz2[2],Wr1[2],Wr2[2];
    double tmpZ[2],tmpR[2],tmpP[2][2];
    double vp,xr,yr,zr,rr,xc,yc,zc,rc;
    double coss[D->numMode],sins[D->numMode];
    double cos1[D->numMode],sin1[D->numMode];
    double cos2[D->numMode],sin2[D->numMode];
    double dz,dr,dt,inverDt,drBydt,dzBydt;

    ptclList *p;
    LoadList *LL;
    Particle **particle;
    particle=D->particle;

//    double maximum();
//    double minimum();
    double coeff[nSpecies];

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    nxSub=D->nxSub;      nySub=D->nySub;
    minRSub=D->minYSub;

    dt=D->dt; dz=D->dz; dr=D->dr; inverDt=1.0/dt;
    numMode=D->numMode;
    drBydt=dr/dt;  dzBydt=dz/dt;

    s=0;
    LL=D->loadList;
    while(LL->next) {
       coeff[s]=LL->density/LL->criticalDensity;
       LL=LL->next;
       s++;
    }
    int myrank,nTasks;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//    if(myrank==0) { istart=D->istart+2; } else ;
//    if(myrank==D->L-1) { iend=D->iend-3; } else ;

    //initialize J
    if(D->fieldIonization==OFF) {
      if(D->fieldType==NDFX)  {
        for(m=0; m<numMode; m++)
         for(i=0; i<nxSub+5; i++)
           for(j=0; j<nySub+5; j++)   {
             D->JzCR[m][i][j]=D->JzR[m][i][j];
             D->JpCR[m][i][j]=D->JpR[m][i][j];
             D->JrCR[m][i][j]=D->JrR[m][i][j];
             D->JzCI[m][i][j]=D->JzI[m][i][j];
             D->JpCI[m][i][j]=D->JpI[m][i][j];
             D->JrCI[m][i][j]=D->JrI[m][i][j];
           }
      } else ;
      for(m=0; m<numMode; m++)
       for(i=0; i<nxSub+5; i++)
         for(j=0; j<nySub+5; j++)   {
           D->JzR[m][i][j]=0.0;
           D->JrR[m][i][j]=0.0;
           D->JpR[m][i][j]=0.0;
           D->JzI[m][i][j]=0.0;
           D->JrI[m][i][j]=0.0;
           D->JpI[m][i][j]=0.0;
         }
    }
    alpha=2.0;

    for(i=istart; i<iend; i++)
      for(j=jstart+2; j<jend; j++)
      {
        for(s=0; s<nSpecies; s++)
        {
          p=particle[i][j].head[s]->pt;
          while(p)
          {
            weight=p->weight*p->charge;
            gamma=sqrt(1.0+p->pz*p->pz+p->px*p->px+p->py*p->py);
            x1=p->oldX; y1=p->oldY; r1=sqrt(x1*x1+y1*y1);
            x2=p->x;    y2=p->y;    r2=sqrt(x2*x2+y2*y2);
            z2=p->z+i;  z1=p->oldZ;

            i1=(int)z1;           j1=(int)r1;
            i2=(int)z2;           j2=(int)r2;

            //rally calculation
            if(i1==i2)       zr=0.5*(z1+z2);
            else             zr=maximum(i1*1.0,i2*1.0);
            if(j1==j2)       rr=0.5*(r1+r2);
            else             rr=maximum(j1*1.0,j2*1.0);
            xr=x1; yr=y1;

            calculaionRally(&xr,&yr,rr,x1,x2,y1,y2,r1,r2,iteration);

//step 1 --------------------------------------------------------------
            rc=0.5*(r1+rr);     zc=0.5*(z1+zr);
            Wr[0]=((j1+1.0)*(j1+1.0)-rc*rc)/(2.0*j1+1.0); Wr[1]=1.0-Wr[0];
            Wz[1]=zc-i1;        Wz[0]=1.0-Wz[1];
            xc=x1; yc=y1;            
            calculaionRally(&xc,&yc,rc,x1,xr,y1,yr,r1,rr,iteration);
            coss[1]=xc/rc; sins[1]=yc/rc;
            for(m=2; m<numMode; m++) {
              coss[m]=coss[m-1]*coss[1]-sins[m-1]*sins[1];
              sins[m]=sins[m-1]*coss[1]+coss[m-1]*sins[1];
            }

            Fz=(zr-z1)*dzBydt; Fr=(rr-r1)*drBydt;
            factor=weight*coeff[s]/(2.0*j1+1.0);

            tmpZ[0]=Fz*Wr[0]*factor;
            tmpZ[1]=Fz*Wr[1]*factor;
            for(jj=0; jj<2; jj++)
              D->JzR[0][i1][jj+j1+jstart]+=tmpZ[jj];
            for(m=1; m<numMode; m++) 
              for(jj=0; jj<2; jj++) {
                D->JzR[m][i1][jj+j1+jstart]+=tmpZ[jj]*coss[m]*alpha;
                D->JzI[m][i1][jj+j1+jstart]-=tmpZ[jj]*sins[m]*alpha;
              }
            
            tmpR[0]=Fr*Wz[0]*factor;
            tmpR[1]=Fr*Wz[1]*factor;
            for(ii=0; ii<2; ii++)
              D->JrR[0][ii+i1][j1+jstart]+=tmpR[ii];
            for(m=1; m<numMode; m++)  
              for(ii=0; ii<2; ii++) {
                D->JrR[m][ii+i1][j1+jstart]+=tmpR[ii]*coss[m]*alpha;
                D->JrI[m][ii+i1][j1+jstart]-=tmpR[ii]*sins[m]*alpha;
              }

              //Jp
            vp=coss[1]*(yr-y1)*drBydt-sins[1]*(xr-x1)*drBydt;
            factor=vp*weight*coeff[s]/(2.0*j1+1.0);
            tmpP[0][0]=Wz[0]*Wr[0]*factor;
            tmpP[1][0]=Wz[1]*Wr[0]*factor;
            tmpP[0][1]=Wz[0]*Wr[1]*factor;
            tmpP[1][1]=Wz[1]*Wr[1]*factor;
            for(ii=0; ii<2; ii++)
              for(jj=0; jj<2; jj++)
                D->JpR[0][i1+ii][jj+j1+jstart]+=tmpP[ii][jj];

			   if(D->currentCons==Lifschitz) {
              for(m=1; m<numMode; m++)  
                for(ii=0; ii<2; ii++) 
                  for(jj=0; jj<2; jj++) {
                    D->JpR[m][i1+ii][jj+j1+jstart]+=tmpP[ii][jj]*coss[m]*alpha;
                    D->JpI[m][i1+ii][jj+j1+jstart]-=tmpP[ii][jj]*sins[m]*alpha;
			         }
				} else if(D->currentCons==Davidson) {
               //delta function and weights
              cos1[1]=x1/r1; sin1[1]=y1/r1;
              for(m=2; m<numMode; m++) {
                cos1[m]=cos1[m-1]*cos1[1]-sin1[m-1]*sin1[1];
                sin1[m]=sin1[m-1]*cos1[1]+cos1[m-1]*sin1[1];
              }
              cos2[1]=xr/rr; sin2[1]=yr/rr;
              for(m=2; m<numMode; m++) {
                cos2[m]=cos2[m-1]*cos2[1]-sin2[m-1]*sin2[1];
                sin2[m]=sin2[m-1]*cos2[1]+cos2[m-1]*sin2[1];
              }
              Wz1[1]=z1-i1; Wz1[0]=1.0-Wz1[1];
              Wz2[1]=zr-i1; Wz2[0]=1.0-Wz2[1];
              Wr1[0]=((j1+1.0)*(j1+1.0)-r1*r1)/(2.0*j1+1.0); Wr1[1]=1.0-Wr1[0];
              Wr2[0]=((j1+1.0)*(j1+1.0)-rr*rr)/(2.0*j1+1.0); Wr2[1]=1.0-Wr2[0];

              factor=weight*coeff[s]/(2.0*j1+1.0); 
              for(m=1; m<numMode; m++)  
                for(ii=0; ii<2; ii++) 
                  for(jj=0; jj<2; jj++)  {
                    factM=alpha*(j1+jj)/(m*1.0)*drBydt;
                    D->JpR[m][ii+i1][jj+j1+jstart]+=factM*factor*(Wz2[ii]*Wr2[jj]*(sin2[m]-sins[m])-Wz1[ii]*Wr1[jj]*(sin1[m]-sins[m]));
                    D->JpI[m][ii+i1][jj+j1+jstart]+=factM*factor*(Wz2[ii]*Wr2[jj]*(cos2[m]-coss[m])-Wz1[ii]*Wr1[jj]*(cos1[m]-coss[m]));
                  }
				}

//step 2 --------------------------------------------------------------
            rc=0.5*(rr+r2);     zc=0.5*(zr+z2);
            Wr[0]=((j2+1.0)*(j2+1.0)-rc*rc)/(2.0*j2+1.0); Wr[1]=1.0-Wr[0];
            Wz[1]=zc-i2;        Wz[0]=1.0-Wz[1];
            xc=xr; yc=yr;            
            calculaionRally(&xc,&yc,rc,xr,x2,yr,y2,rr,r2,iteration);
            coss[1]=xc/rc; sins[1]=yc/rc;
            for(m=2; m<numMode; m++) {
              coss[m]=coss[m-1]*coss[1]-sins[m-1]*sins[1];
              sins[m]=sins[m-1]*coss[1]+coss[m-1]*sins[1];
            }

            Fz=(z2-zr)*dzBydt; Fr=(r2-rr)*drBydt;
            factor=weight*coeff[s]/(2.0*j2+1.0);

            tmpZ[0]=Fz*Wr[0]*factor;
            tmpZ[1]=Fz*Wr[1]*factor;
            for(jj=0; jj<2; jj++)
              D->JzR[0][i2][jj+j2+jstart]+=tmpZ[jj];
            for(m=1; m<numMode; m++)  
              for(jj=0; jj<2; jj++) {
                D->JzR[m][i2][jj+j2+jstart]+=tmpZ[jj]*coss[m]*alpha;
                D->JzI[m][i2][jj+j2+jstart]-=tmpZ[jj]*sins[m]*alpha;
              }
            
            tmpR[0]=Fr*Wz[0]*factor;
            tmpR[1]=Fr*Wz[1]*factor;
            for(ii=0; ii<2; ii++)
              D->JrR[0][ii+i2][j2+jstart]+=tmpR[ii];
            for(m=1; m<numMode; m++)  {
              for(ii=0; ii<2; ii++) {
                D->JrR[m][ii+i2][j2+jstart]+=tmpR[ii]*coss[m]*alpha;
                D->JrI[m][ii+i2][j2+jstart]-=tmpR[ii]*sins[m]*alpha;
              }
            }

              //Jp
            vp=coss[1]*(y2-yr)*drBydt-sins[1]*(x2-xr)*drBydt;
            factor=vp*weight*coeff[s]/(2.0*j2+1.0);
            tmpP[0][0]=Wz[0]*Wr[0]*factor;
            tmpP[1][0]=Wz[1]*Wr[0]*factor;
            tmpP[0][1]=Wz[0]*Wr[1]*factor;
            tmpP[1][1]=Wz[1]*Wr[1]*factor;
            for(ii=0; ii<2; ii++)
              for(jj=0; jj<2; jj++)
                D->JpR[0][i2+ii][jj+j2+jstart]+=tmpP[ii][jj];
				if(D->currentCons==Lifschitz) {
              for(m=1; m<numMode; m++)  
                for(ii=0; ii<2; ii++) 
                  for(jj=0; jj<2; jj++) {
                    D->JpR[m][i2+ii][jj+j2+jstart]+=tmpP[ii][jj]*coss[m]*alpha;
                    D->JpI[m][i2+ii][jj+j2+jstart]-=tmpP[ii][jj]*sins[m]*alpha;
			         }
				} else if(D->currentCons==Davidson) {
               //delta function and weights
              cos1[1]=xr/rr; sin1[1]=yr/rr;
              for(m=2; m<numMode; m++) {
                cos1[m]=cos1[m-1]*cos1[1]-sin1[m-1]*sin1[1];
                sin1[m]=sin1[m-1]*cos1[1]+cos1[m-1]*sin1[1];
              }
              cos2[1]=x2/r2; sin2[1]=y2/r2;
              for(m=2; m<numMode; m++) {
                cos2[m]=cos2[m-1]*cos2[1]-sin2[m-1]*sin2[1];
                sin2[m]=sin2[m-1]*cos2[1]+cos2[m-1]*sin2[1];
              }
              Wz1[1]=zr-i2; Wz1[0]=1.0-Wz1[1];
              Wz2[1]=z2-i2; Wz2[0]=1.0-Wz2[1];
              Wr1[0]=((j2+1.0)*(j2+1.0)-rr*rr)/(2.0*j2+1.0); Wr1[1]=1.0-Wr1[0];
              Wr2[0]=((j2+1.0)*(j2+1.0)-r2*r2)/(2.0*j2+1.0); Wr2[1]=1.0-Wr2[0];

              factor=weight*coeff[s]/(2.0*j2+1.0); 
              for(m=1; m<numMode; m++)  
                for(ii=0; ii<2; ii++) 
                  for(jj=0; jj<2; jj++)  {
                    factM=alpha*(j2+jj)/(m*1.0)*drBydt;
                    D->JpR[m][ii+i2][jj+j2+jstart]+=factM*factor*(Wz2[ii]*Wr2[jj]*(sin2[m]-sins[m])-Wz1[ii]*Wr1[jj]*(sin1[m]-sins[m]));
                    D->JpI[m][ii+i2][jj+j2+jstart]+=factM*factor*(Wz2[ii]*Wr2[jj]*(cos2[m]-coss[m])-Wz1[ii]*Wr1[jj]*(cos1[m]-coss[m]));
                  }
				}
           
            p=p->next;
          }    //End of while(p)

        }    //End of for(s)     
      }      //End of for(i,j)

    // for Axis
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jstart+2; j++)
      {
        for(s=0; s<nSpecies; s++)
        {
          p=particle[i][j].head[s]->pt;
          while(p)
          {
            weight=p->weight*p->charge;
            gamma=sqrt(1.0+p->pz*p->pz+p->px*p->px+p->py*p->py);
            x1=p->oldX; y1=p->oldY; r1=sqrt(x1*x1+y1*y1);
            x2=p->x;    y2=p->y;    r2=sqrt(x2*x2+y2*y2);
            z2=p->z+i;  z1=p->oldZ;

            i1=(int)z1;           j1=(int)r1;
            i2=(int)z2;           j2=(int)r2;

            //rally calculation
            if(i1==i2)       zr=0.5*(z1+z2);
            else             zr=maximum(i1*1.0,i2*1.0);
            if(j1==j2)       rr=0.5*(r1+r2);
            else             rr=maximum(j1*1.0,j2*1.0);
            xr=x1; yr=y1;

            calculaionRally(&xr,&yr,rr,x1,x2,y1,y2,r1,r2,iteration);

//step 1 --------------------------------------------------------------
            rc=0.5*(r1+rr);     zc=0.5*(z1+zr);
            Wr[0]=((j1+1.0)*(j1+1.0)-rc*rc)/(2.0*j1+1.0); Wr[1]=1.0-Wr[0];
            Wz[1]=zc-i1;        Wz[0]=1.0-Wz[1];
            xc=x1; yc=y1;            
            calculaionRally(&xc,&yc,rc,x1,xr,y1,yr,r1,rr,iteration);
            coss[1]=xc/rc; sins[1]=yc/rc;
            for(m=2; m<numMode; m++) {
              coss[m]=coss[m-1]*coss[1]-sins[m-1]*sins[1];
              sins[m]=sins[m-1]*coss[1]+coss[m-1]*sins[1];
            }

            Fz=(zr-z1)*dzBydt; Fr=(rr-r1)*drBydt;
            factor=weight*coeff[s]/(2.0*j1+1.0);

            if(j1==0) {
              tmpZ[0]=Fz*Wr[0]*factor*2.0;
              tmpZ[1]=Fz*Wr[1]*factor;
              for(jj=0; jj<2; jj++)
                D->JzR[0][i1][jj+j1+jstart]+=tmpZ[jj];
              for(m=1; m<numMode; m++)  
                for(jj=1; jj<2; jj++) {
                  D->JzR[m][i1][jj+j1+jstart]+=tmpZ[jj]*coss[m]*alpha;
                  D->JzI[m][i1][jj+j1+jstart]-=tmpZ[jj]*sins[m]*alpha;
                }
				} else {
              tmpZ[0]=Fz*Wr[0]*factor;
              tmpZ[1]=Fz*Wr[1]*factor;
              for(jj=0; jj<2; jj++)
                D->JzR[0][i1][jj+j1+jstart]+=tmpZ[jj];
              for(m=1; m<numMode; m++)  
                for(jj=0; jj<2; jj++) {
                  D->JzR[m][i1][jj+j1+jstart]+=tmpZ[jj]*coss[m]*alpha;
                  D->JzI[m][i1][jj+j1+jstart]-=tmpZ[jj]*sins[m]*alpha;
                }
				}
            
            tmpR[0]=Fr*Wz[0]*factor;
            tmpR[1]=Fr*Wz[1]*factor;
            for(ii=0; ii<2; ii++)
              D->JrR[0][ii+i1][j1+jstart]+=tmpR[ii];
            for(m=1; m<numMode; m++)  
              for(ii=0; ii<2; ii++) {
                D->JrR[m][ii+i1][j1+jstart]+=tmpR[ii]*coss[m]*alpha;
                D->JrI[m][ii+i1][j1+jstart]-=tmpR[ii]*sins[m]*alpha;
              }

               //delta function and weights

              //Jp
            vp=coss[1]*(yr-y1)*drBydt-sins[1]*(xr-x1)*drBydt;
            factor=vp*weight*coeff[s]/(2.0*j1+1.0);
				if(j1==0) {
              tmpP[0][0]=Wz[0]*Wr[0]*factor*2.0;
              tmpP[1][0]=Wz[1]*Wr[0]*factor*2.0;
              tmpP[0][1]=Wz[0]*Wr[1]*factor;
              tmpP[1][1]=Wz[1]*Wr[1]*factor;
              for(ii=0; ii<2; ii++)
                for(jj=0; jj<2; jj++)
                  D->JpR[0][i1+ii][jj+j1+jstart]+=tmpP[ii][jj];
				  if(D->currentCons==Lifschitz) {
                m=1;
                for(ii=0; ii<2; ii++) 
                  for(jj=0; jj<2; jj++) {
                    D->JpR[m][i1+ii][jj+j1+jstart]+=tmpP[ii][jj]*coss[m]*alpha;
                    D->JpI[m][i1+ii][jj+j1+jstart]-=tmpP[ii][jj]*sins[m]*alpha;
			       }
				    jj=1;
                for(m=2; m<numMode; m++)  
                  for(ii=0; ii<2; ii++) {
                    D->JpR[m][i1+ii][jj+j1+jstart]+=tmpP[ii][jj]*coss[m]*alpha;
                    D->JpI[m][i1+ii][jj+j1+jstart]-=tmpP[ii][jj]*sins[m]*alpha;
                  }
				  } else if(D->currentCons==Davidson) {
                cos1[1]=x1/r1; sin1[1]=y1/r1;
                for(m=2; m<numMode; m++) {
                  cos1[m]=cos1[m-1]*cos1[1]-sin1[m-1]*sin1[1];
                  sin1[m]=sin1[m-1]*cos1[1]+cos1[m-1]*sin1[1];
                }
                cos2[1]=xr/rr; sin2[1]=yr/rr;
                for(m=2; m<numMode; m++) {
                  cos2[m]=cos2[m-1]*cos2[1]-sin2[m-1]*sin2[1];
                  sin2[m]=sin2[m-1]*cos2[1]+cos2[m-1]*sin2[1];
                }
                Wz1[1]=z1-i1; Wz1[0]=1.0-Wz1[1];
                Wz2[1]=zr-i1; Wz2[0]=1.0-Wz2[1];
                Wr1[0]=((j1+1.0)*(j1+1.0)-r1*r1)/(2.0*j1+1.0); Wr1[1]=1.0-Wr1[0];
                Wr2[0]=((j1+1.0)*(j1+1.0)-rr*rr)/(2.0*j1+1.0); Wr2[1]=1.0-Wr2[0];

                factor=weight*coeff[s]/(2.0*j1+1.0);
                m=1;
                for(ii=0; ii<2; ii++)
                  for(jj=0; jj<1; jj++)  {
                    factM=2*alpha*(j1+jj)/(m*1.0)*drBydt;
                    D->JpR[m][ii+i1][jj+j1+jstart]+=factM*factor*(Wz2[ii]*Wr2[jj]*(sin2[m]-sins[m])-Wz1[ii]*Wr1[jj]*(sin1[m]-sins[m]));
                    D->JpI[m][ii+i1][jj+j1+jstart]+=factM*factor*(Wz2[ii]*Wr2[jj]*(cos2[m]-coss[m])-Wz1[ii]*Wr1[jj]*(cos1[m]-coss[m]));
                  }
                for(ii=0; ii<2; ii++)
                  for(jj=1; jj<2; jj++)  {
                    factM=alpha*(j1+jj)/(m*1.0)*drBydt;
                    D->JpR[m][ii+i1][jj+j1+jstart]+=factM*factor*(Wz2[ii]*Wr2[jj]*(sin2[m]-sins[m])-Wz1[ii]*Wr1[jj]*(sin1[m]-sins[m]));
                    D->JpI[m][ii+i1][jj+j1+jstart]+=factM*factor*(Wz2[ii]*Wr2[jj]*(cos2[m]-coss[m])-Wz1[ii]*Wr1[jj]*(cos1[m]-coss[m]));
  					   }
						for(m=2; m<numMode; m++)  {
						  for(ii=0; ii<2; ii++)
							 for(jj=1; jj<2; jj++)  {	
								factM=alpha*(j1+jj)/(m*1.0)*drBydt;
								D->JpR[m][ii+i1][jj+j1+jstart]+=factM*factor*(Wz2[ii]*Wr2[jj]*(sin2[m]-sins[m])-Wz1[ii]*Wr1[jj]*(sin1[m]-sins[m]));
								D->JpI[m][ii+i1][jj+j1+jstart]+=factM*factor*(Wz2[ii]*Wr2[jj]*(cos2[m]-coss[m])-Wz1[ii]*Wr1[jj]*(cos1[m]-coss[m]));
						    }
						}
				  }
				} else {
					
              tmpP[0][0]=Wz[0]*Wr[0]*factor;
              tmpP[1][0]=Wz[1]*Wr[0]*factor;
              tmpP[0][1]=Wz[0]*Wr[1]*factor;
              tmpP[1][1]=Wz[1]*Wr[1]*factor;
              for(ii=0; ii<2; ii++)
                for(jj=0; jj<2; jj++)
                  D->JpR[0][i1+ii][jj+j1+jstart]+=tmpP[ii][jj];
				  if(D->currentCons==Lifschitz) {
                for(m=1; m<numMode; m++)  
                  for(ii=0; ii<2; ii++) 
                    for(jj=0; jj<2; jj++) {
                      D->JpR[m][i1+ii][jj+j1+jstart]+=tmpP[ii][jj]*coss[m]*alpha;
                      D->JpI[m][i1+ii][jj+j1+jstart]-=tmpP[ii][jj]*sins[m]*alpha;
                    }
				  } else if(D->currentCons==Davidson) {					  
                cos1[1]=x1/r1; sin1[1]=y1/r1;
                for(m=2; m<numMode; m++) {
                  cos1[m]=cos1[m-1]*cos1[1]-sin1[m-1]*sin1[1];
                  sin1[m]=sin1[m-1]*cos1[1]+cos1[m-1]*sin1[1];
                }
                cos2[1]=xr/rr; sin2[1]=yr/rr;
                for(m=2; m<numMode; m++) {
                  cos2[m]=cos2[m-1]*cos2[1]-sin2[m-1]*sin2[1];
                  sin2[m]=sin2[m-1]*cos2[1]+cos2[m-1]*sin2[1];
                }
                Wz1[1]=z1-i1; Wz1[0]=1.0-Wz1[1];
                Wz2[1]=zr-i1; Wz2[0]=1.0-Wz2[1];
                Wr1[0]=((j1+1.0)*(j1+1.0)-r1*r1)/(2.0*j1+1.0); Wr1[1]=1.0-Wr1[0];
                Wr2[0]=((j1+1.0)*(j1+1.0)-rr*rr)/(2.0*j1+1.0); Wr2[1]=1.0-Wr2[0];

                factor=weight*coeff[s]/(2.0*j1+1.0); 
                for(m=1; m<numMode; m++)  
                  for(ii=0; ii<2; ii++) 
                    for(jj=0; jj<2; jj++)  {
                      factM=alpha*(j1+jj)/(m*1.0)*drBydt;
                      D->JpR[m][ii+i1][jj+j1+jstart]+=factM*factor*(Wz2[ii]*Wr2[jj]*(sin2[m]-sins[m])-Wz1[ii]*Wr1[jj]*(sin1[m]-sins[m]));
                      D->JpI[m][ii+i1][jj+j1+jstart]+=factM*factor*(Wz2[ii]*Wr2[jj]*(cos2[m]-coss[m])-Wz1[ii]*Wr1[jj]*(cos1[m]-coss[m]));
                    }
				    }
				  }

//step 2 --------------------------------------------------------------
            rc=0.5*(rr+r2);     zc=0.5*(zr+z2);
            Wr[0]=((j2+1.0)*(j2+1.0)-rc*rc)/(2.0*j2+1.0); Wr[1]=1.0-Wr[0];
            Wz[1]=zc-i2;        Wz[0]=1.0-Wz[1];
            xc=xr; yc=yr;            
            calculaionRally(&xc,&yc,rc,xr,x2,yr,y2,rr,r2,iteration);
            coss[1]=xc/rc; sins[1]=yc/rc;
            for(m=2; m<numMode; m++) {
              coss[m]=coss[m-1]*coss[1]-sins[m-1]*sins[1];
              sins[m]=sins[m-1]*coss[1]+coss[m-1]*sins[1];
            }

            Fz=(z2-zr)*dzBydt; Fr=(r2-rr)*drBydt;
            factor=weight*coeff[s]/(2.0*j2+1.0);

            if(j2==0) {
              tmpZ[0]=Fz*Wr[0]*factor*2.0;
              tmpZ[1]=Fz*Wr[1]*factor;
              for(jj=0; jj<2; jj++)
                D->JzR[0][i2][jj+j2+jstart]+=tmpZ[jj];
              for(m=1; m<numMode; m++)  {
                for(jj=1; jj<2; jj++) {
                  D->JzR[m][i2][jj+j2+jstart]+=tmpZ[jj]*coss[m]*alpha;
                  D->JzI[m][i2][jj+j2+jstart]-=tmpZ[jj]*sins[m]*alpha;
                }
              }
				} else {
              tmpZ[0]=Fz*Wr[0]*factor;
              tmpZ[1]=Fz*Wr[1]*factor;
              for(jj=0; jj<2; jj++)
                D->JzR[0][i2][jj+j2+jstart]+=tmpZ[jj];
              for(m=1; m<numMode; m++)  {
                for(jj=0; jj<2; jj++) {
                  D->JzR[m][i2][jj+j2+jstart]+=tmpZ[jj]*coss[m]*alpha;
                  D->JzI[m][i2][jj+j2+jstart]-=tmpZ[jj]*sins[m]*alpha;
                }
              }
				}

            
            tmpR[0]=Fr*Wz[0]*factor;
            tmpR[1]=Fr*Wz[1]*factor;
            for(ii=0; ii<2; ii++)
              D->JrR[0][ii+i2][j2+jstart]+=tmpR[ii];
            for(m=1; m<numMode; m++)  
              for(ii=0; ii<2; ii++) {
                D->JrR[m][ii+i2][j2+jstart]+=tmpR[ii]*coss[m]*alpha;
                D->JrI[m][ii+i2][j2+jstart]-=tmpR[ii]*sins[m]*alpha;
              }
            

            //Jp
            vp=coss[1]*(y2-yr)*drBydt-sins[1]*(x2-xr)*drBydt;
            factor=vp*weight*coeff[s]/(2.0*j2+1.0);
				if(j2==0) {
              tmpP[0][0]=Wz[0]*Wr[0]*factor*2.0;
              tmpP[1][0]=Wz[1]*Wr[0]*factor*2.0;
              tmpP[0][1]=Wz[0]*Wr[1]*factor;
              tmpP[1][1]=Wz[1]*Wr[1]*factor;
              for(ii=0; ii<2; ii++)
                for(jj=1; jj<2; jj++)
                  D->JpR[0][i2+ii][jj+j2+jstart]+=tmpP[ii][jj];

				  if(D->currentCons==Lifschitz) {
                m=1;
                for(ii=0; ii<2; ii++) 
                  for(jj=0; jj<2; jj++) {
                    D->JpR[m][i2+ii][jj+j2+jstart]+=tmpP[ii][jj]*coss[m]*alpha;
                    D->JpI[m][i2+ii][jj+j2+jstart]-=tmpP[ii][jj]*sins[m]*alpha;
			         }
				    jj=1;
                for(m=2; m<numMode; m++)  
                  for(ii=0; ii<2; ii++) {
                    D->JpR[m][i2+ii][jj+j2+jstart]+=tmpP[ii][jj]*coss[m]*alpha;
                    D->JpI[m][i2+ii][jj+j2+jstart]-=tmpP[ii][jj]*sins[m]*alpha;
					   }
              } else if(D->currentCons==Davidson) {
				    //delta function and weights
                cos1[1]=xr/rr; sin1[1]=yr/rr;
                for(m=2; m<numMode; m++) {
                  cos1[m]=cos1[m-1]*cos1[1]-sin1[m-1]*sin1[1];
                  sin1[m]=sin1[m-1]*cos1[1]+cos1[m-1]*sin1[1];
                }
                cos2[1]=x2/r2; sin2[1]=y2/r2;
                for(m=2; m<numMode; m++) {
                  cos2[m]=cos2[m-1]*cos2[1]-sin2[m-1]*sin2[1];
                  sin2[m]=sin2[m-1]*cos2[1]+cos2[m-1]*sin2[1];
                }
                Wz1[1]=zr-i2; Wz1[0]=1.0-Wz1[1];
                Wz2[1]=z2-i2; Wz2[0]=1.0-Wz2[1];
                Wr1[0]=((j2+1.0)*(j2+1.0)-rr*rr)/(2.0*j2+1.0); Wr1[1]=1.0-Wr1[0];
                Wr2[0]=((j2+1.0)*(j2+1.0)-r2*r2)/(2.0*j2+1.0); Wr2[1]=1.0-Wr2[0];

                factor=weight*coeff[s]/(2.0*j2+1.0);
                m=1;
                for(ii=0; ii<2; ii++)
                  for(jj=0; jj<1; jj++)  {
                    factM=2*alpha*(j2+jj)/(m*1.0)*drBydt;
                    D->JpR[m][ii+i2][jj+j2+jstart]+=factM*factor*(Wz2[ii]*Wr2[jj]*(sin2[m]-sins[m])-Wz1[ii]*Wr1[jj]*(sin1[m]-sins[m]));
                    D->JpI[m][ii+i2][jj+j2+jstart]+=factM*factor*(Wz2[ii]*Wr2[jj]*(cos2[m]-coss[m])-Wz1[ii]*Wr1[jj]*(cos1[m]-coss[m]));
                  }
                for(ii=0; ii<2; ii++)
                  for(jj=1; jj<2; jj++)  {
                    factM=alpha*(j2+jj)/(m*1.0)*drBydt;
                    D->JpR[m][ii+i2][jj+j2+jstart]+=factM*factor*(Wz2[ii]*Wr2[jj]*(sin2[m]-sins[m])-Wz1[ii]*Wr1[jj]*(sin1[m]-sins[m]));
                    D->JpI[m][ii+i2][jj+j2+jstart]+=factM*factor*(Wz2[ii]*Wr2[jj]*(cos2[m]-coss[m])-Wz1[ii]*Wr1[jj]*(cos1[m]-coss[m]));
  					   }
						for(m=2; m<numMode; m++)  {
						  for(ii=0; ii<2; ii++)
							 for(jj=1; jj<2; jj++)  {
								factM=alpha*(j2+jj)/(m*1.0)*drBydt;
								D->JpR[m][ii+i2][jj+j2+jstart]+=factM*factor*(Wz2[ii]*Wr2[jj]*(sin2[m]-sins[m])-Wz1[ii]*Wr1[jj]*(sin1[m]-sins[m]));
								D->JpI[m][ii+i2][jj+j2+jstart]+=factM*factor*(Wz2[ii]*Wr2[jj]*(cos2[m]-coss[m])-Wz1[ii]*Wr1[jj]*(cos1[m]-coss[m]));
						    }
						}
				  }
				} else {
              tmpP[0][0]=Wz[0]*Wr[0]*factor;
              tmpP[1][0]=Wz[1]*Wr[0]*factor;
              tmpP[0][1]=Wz[0]*Wr[1]*factor;
              tmpP[1][1]=Wz[1]*Wr[1]*factor;
              for(ii=0; ii<2; ii++)
                for(jj=0; jj<2; jj++)
                  D->JpR[0][i2+ii][jj+j2+jstart]+=tmpP[ii][jj];
				  if(D->currentCons==Lifschitz) {
                for(m=1; m<numMode; m++)  
                  for(ii=0; ii<2; ii++) 
                    for(jj=0; jj<2; jj++) {
                      D->JpR[m][i2+ii][jj+j2+jstart]+=tmpP[ii][jj]*coss[m]*alpha;
                      D->JpI[m][i2+ii][jj+j2+jstart]-=tmpP[ii][jj]*sins[m]*alpha;
						  }
				  } else if(D->currentCons==Davidson) {
				    //delta function and weights
                cos1[1]=xr/rr; sin1[1]=yr/rr;
                for(m=2; m<numMode; m++) {
                  cos1[m]=cos1[m-1]*cos1[1]-sin1[m-1]*sin1[1];
                  sin1[m]=sin1[m-1]*cos1[1]+cos1[m-1]*sin1[1];
                }
                cos2[1]=x2/r2; sin2[1]=y2/r2;
                for(m=2; m<numMode; m++) {
                  cos2[m]=cos2[m-1]*cos2[1]-sin2[m-1]*sin2[1];
                  sin2[m]=sin2[m-1]*cos2[1]+cos2[m-1]*sin2[1];
                }
                Wz1[1]=zr-i2; Wz1[0]=1.0-Wz1[1];
                Wz2[1]=z2-i2; Wz2[0]=1.0-Wz2[1];
                Wr1[0]=((j2+1.0)*(j2+1.0)-rr*rr)/(2.0*j2+1.0); Wr1[1]=1.0-Wr1[0];
                Wr2[0]=((j2+1.0)*(j2+1.0)-r2*r2)/(2.0*j2+1.0); Wr2[1]=1.0-Wr2[0];

                factor=weight*coeff[s]/(2.0*j2+1.0); 
                for(m=1; m<numMode; m++)  {
                  for(ii=0; ii<2; ii++) 
                    for(jj=0; jj<2; jj++)  {
                      factM=alpha*(j2+jj)/(m*1.0)*drBydt;
                      D->JpR[m][ii+i2][jj+j2+jstart]+=factM*factor*(Wz2[ii]*Wr2[jj]*(sin2[m]-sins[m])-Wz1[ii]*Wr1[jj]*(sin1[m]-sins[m]));
                      D->JpI[m][ii+i2][jj+j2+jstart]+=factM*factor*(Wz2[ii]*Wr2[jj]*(cos2[m]-coss[m])-Wz1[ii]*Wr1[jj]*(cos1[m]-coss[m]));
                    }
                }
				  }
				  //lala
				}


            p=p->next;
          }    //End of while(p)

        }    //End of for(s)     
      }      //End of for(i,j)

}

void calculaionRally(double *xr,double *yr,double rr,double x1,double x2,double y1,double y2,double r1,double r2,int iteration)
{
   int mode;
   double alpha,xc,yc,sign;

   xc=0.5*(x1+x2);
   yc=0.5*(y1+y2);
   if(x1!=x2)  {
     alpha=yc/xc;
     sign=xc/fabs(xc);
     *xr=sign*sqrt(rr*rr/(1.0+alpha*alpha));
     *yr=alpha*(*xr);
     mode=1;
   } else { 
     if(yc!=0.0) {
       *xr=x1; 
       sign=yc/fabs(yc);
//if(isnan(sign)) printf("sign=%g, yc=%g\n",sign,yc);
       *yr=sign*sqrt(rr*rr-x1*x1);
       mode=2;
     } else {
       sign=x1/fabs(x1);
       *xr=sign*rr;
//if(isnan(sign)) printf("sign=%g, yc=%g\n",sign,yc);
       *yr=0.0;
       mode=3;
     }
   }
}




void updateCurrent_villasenor(Domain *D,int nSpecies,int iteration)
{
    int n,i,j,m,s,i1,j1,i2,j2,idI,idJ,numMode,ii[3],jj[3],mode;
    int istart,iend,jstart,jend,index,minRSub;
    int nxSub,nySub,difI,difJ,cnt,maxI,maxJ;
    int indexIC,indexJC,indexI1,indexJ1,indexI2,indexJ2;
    double dz,dr,dt,dzBydt,drBydt,factor,tmp;
    double WzC[2],WrC[2],Wz1[2],Wr1[2],portion[3];
    double x[4],y[4],z[4],r[4],xc,yc,zc,rc,ddZ,ddR,ddX,ddY;
    double z1,x1,y1,r1,z2,x2,y2,r2,rr;
    double tmpP1,tmpP2,tmpP3,tmpP4;
    double J0,vp,alpha,weight,gamma,capR[3],capZ[3],delZ[3],delR[3];
    double cossC[D->numMode],sinsC[D->numMode];
    double coss[D->numMode],sins[D->numMode];
    double coss1[D->numMode],sins1[D->numMode];
    double coss2[D->numMode],sins2[D->numMode];
    
    ptclList *p;
    LoadList *LL;   
    Particle **particle;
    particle=D->particle;

//    double maximum();
//    double minimum();

    double coeff[nSpecies];

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    nxSub=D->nxSub;      nySub=D->nySub;
    minRSub=D->minYSub;
   
    dt=D->dt; dz=D->dz; dr=D->dr; dzBydt=dz/dt; drBydt=dr/dt;
    numMode=D->numMode;

    s=0;
    LL=D->loadList;
    while(LL->next) {
       coeff[s]=LL->density/LL->criticalDensity;
       LL=LL->next;
       s++;
    }

    int myrank,nTasks;
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//    if(myrank==0) { istart=D->istart+2; } else ;
//    if(myrank==D->L-1) { iend=D->iend-3; } else ;

    //initialize J
    if(D->fieldIonization==OFF) {
      if(D->fieldType==NDFX)  {
        for(m=0; m<numMode; m++)  
         for(i=0; i<nxSub+5; i++)
           for(j=0; j<nySub+5; j++)   {
             D->JzCR[m][i][j]=D->JzR[m][i][j];
             D->JpCR[m][i][j]=D->JpR[m][i][j];
             D->JrCR[m][i][j]=D->JrR[m][i][j];
             D->JzCI[m][i][j]=D->JzI[m][i][j];
             D->JpCI[m][i][j]=D->JpI[m][i][j];
             D->JrCI[m][i][j]=D->JrI[m][i][j];
           }
      } else ;
      for(m=0; m<numMode; m++)  
       for(i=0; i<nxSub+5; i++)
         for(j=0; j<nySub+5; j++)   {
           D->JzR[m][i][j]=0.0;
           D->JrR[m][i][j]=0.0;
           D->JpR[m][i][j]=0.0;
           D->JzI[m][i][j]=0.0;
           D->JrI[m][i][j]=0.0;
           D->JpI[m][i][j]=0.0;
         }
    }
    alpha=2.0;

    for(i=istart; i<iend; i++)
      for(j=jstart+2; j<jend; j++)
      {
        for(s=0; s<nSpecies; s++)
        {
          p=particle[i][j].head[s]->pt;     
          while(p) 
          {
            cnt=0;
            weight=p->weight*p->charge;
            gamma=sqrt(1.0+p->pz*p->pz+p->px*p->px+p->py*p->py);
            x1=p->oldX; y1=p->oldY; r1=sqrt(x1*x1+y1*y1);
            x2=p->x;    y2=p->y;    r2=sqrt(x2*x2+y2*y2);
            z2=p->z+i;  z1=p->oldZ;
            z[0]=z1; x[0]=x1; y[0]=y1; r[0]=r1;

            i1=(int)(z1); j1=(int)(r1)+jstart;
            i2=(int)(z2); j2=(int)(r2)+jstart;    
            ii[0]=i1;       jj[0]=j1;
            difI=i2-i1; difJ=j2-j1;
            ddZ=z2-z1; ddR=r2-r1; ddX=x2-x1; ddY=y2-y1;
            maxI=MAX(i1,i2); maxJ=MAX(j1,j2);
            if(difI==0 && difJ==0) { 
              cnt=1; mode=1; 
              z[1]=z2; x[1]=x2; y[1]=y2; r[1]=r2;
              ii[1]=i2;       jj[1]=j2;
                
              capZ[0]=z[0]-(ii[0]+0.5);
              capR[0]=r[0]+jstart-(jj[0]+0.5);
              delZ[0]=ddZ;   delR[0]=ddR;
              portion[0]=1.0;
            } else if(difI!=0 && difJ==0) { 
              cnt=2; mode=2;
              z[2]=z2; x[2]=x2; y[2]=y2; r[2]=r2;
              ii[1]=i2;       jj[1]=j2;

              capZ[0]=z[0]-(ii[0]+0.5);
              capR[0]=r[0]+jstart-(jj[0]+0.5);
              delZ[0]=maxI-z[0]; 
              delR[0]=ddR*(delZ[0]/ddZ);
              capZ[1]=0.5*(ii[0]+ii[1])-ii[1]; 
              capR[1]=capR[0]+delR[0];
              delZ[1]=ddZ-delZ[0];    delR[1]=ddR-delR[0];

              z[1]=maxI       ; r[1]=r[0]+delR[0];
              xc=ddX/ddZ*(z[1]-z[0])+x[0];
              yc=ddY/ddZ*(z[1]-z[0])+y[0];
              rc=sqrt(xc*xc+yc*yc);
              x[1]=xc/rc*r[1]; y[1]=yc/rc*r[1]; 
              for(n=0; n<cnt; n++)
                portion[n]=sqrt((delZ[n]*delZ[n]+delR[n]*delR[n])/(ddZ*ddZ+ddR*ddR));
            } else if(difI==0 && difJ!=0) { 
              cnt=2; mode=3;
              z[2]=z2; x[2]=x2; y[2]=y2; r[2]=r2;
              ii[1]=i2;       jj[1]=j2;

              capZ[0]=z[0]-(ii[0]+0.5);
              capR[0]=r[0]+jstart-(jj[0]+0.5);
              delR[0]=maxJ-jstart-r[0]; 
              delZ[0]=ddZ*(delR[0]/ddR);
              capR[1]=0.5*(jj[0]+jj[1])-jj[1]; 
              capZ[1]=capZ[0]+delZ[0];
              delR[1]=ddR-delR[0];    delZ[1]=ddZ-delZ[0];

              z[1]=z[0]+delZ[0]; r[1]=maxJ-jstart;
              xc=ddX/ddZ*(z[1]-z[0])+x[0];
              yc=ddY/ddZ*(z[1]-z[0])+y[0];
              rc=sqrt(xc*xc+yc*yc);
              x[1]=xc/rc*r[1]; y[1]=yc/rc*r[1]; 
              if(ddZ==0.0) cnt=0; else ;
              for(n=0; n<cnt; n++)
                portion[n]=sqrt((delZ[n]*delZ[n]+delR[n]*delR[n])/(ddZ*ddZ+ddR*ddR));
            } else { 
              if(fabs(ddR/ddZ)<fabs((maxJ-jstart-r[0])/(maxI-z[0]))) {
                cnt=3; mode=4;
                z[3]=z2; x[3]=x2; y[3]=y2; r[3]=r2;
                ii[2]=i2;       jj[2]=j2;

                capZ[0]=z[0]-(ii[0]+0.5);
                capR[0]=r[0]+jstart-(jj[0]+0.5);
                z[1]=maxI;       delZ[0]=z[1]-z[0]; 
                ii[1]=ii[2];     capZ[1]=0.5*(ii[0]+ii[1])-ii[1]; 
                delR[0]=ddR*delZ[0]/ddZ; r[1]=r[0]+delR[0];
                jj[1]=jj[0];     capR[1]=r[1]+jstart-maxJ;
                delR[1]=maxJ-jstart-r[1]; delZ[1]=ddZ*(delR[1]/ddR);
                z[2]=z[1]+delZ[1];    capZ[2]=z[2]-(ii[2]+0.5); 
                r[2]=maxJ-jstart; capR[2]=0.5*(jj[0]+jj[2])-jj[2];
                delZ[2]=ddZ-delZ[0]-delZ[1];
                delR[2]=ddR-delR[0]-delR[1];

                xc=ddX/ddZ*(z[1]-z[0])+x[0];                
                yc=ddY/ddZ*(z[1]-z[0])+y[0]; 
                rc=sqrt(xc*xc+yc*yc);
                x[1]=xc/rc*r[1]; y[1]=yc/rc*r[1]; 
                xc=ddX/ddZ*(z[2]-z[0])+x[0];                
                yc=ddY/ddZ*(z[2]-z[0])+y[0];                
                rc=sqrt(xc*xc+yc*yc);
                x[2]=xc/rc*r[2]; y[2]=yc/rc*r[2]; 
                for(n=0; n<cnt; n++)
                  portion[n]=sqrt((delZ[n]*delZ[n]+delR[n]*delR[n])/(ddZ*ddZ+ddR*ddR));
              } else if(fabs(ddR/ddZ)>fabs((maxJ-jstart-r[0])/(maxI-z[0]))) {
                cnt=3; mode=5;
                z[3]=z2; x[3]=x2; y[3]=y2; r[3]=r2;
                ii[2]=i2;       jj[2]=j2;

                capZ[0]=z[0]-(ii[0]+0.5);
                capR[0]=r[0]+jstart-(jj[0]+0.5);
                r[1]=maxJ-jstart; delR[0]=r[1]-r[0]; 
                jj[1]=jj[2];      capR[1]=0.5*(jj[1]+jj[1])-jj[1];
                delZ[0]=ddZ*delR[0]/ddR; z[1]=z[0]+delZ[0];
                ii[1]=ii[0];      capZ[1]=z[1]-maxI;
                delZ[1]=maxI-z[1]; delR[1]=ddR*(delZ[1]/ddZ);
                r[2]=r[1]+delR[1]; capR[2]=r[2]+jstart-(jj[2]+0.5); 
                z[2]=maxI;         capZ[2]=0.5*(ii[0]+ii[2])-ii[2];
                delZ[2]=ddZ-delZ[0]-delZ[1];
                delR[2]=ddR-delR[0]-delR[1];

                xc=ddX/ddZ*(z[1]-z[0])+x[0];                
                yc=ddY/ddZ*(z[1]-z[0])+y[0];                
                rc=sqrt(xc*xc+yc*yc);
                x[1]=xc/rc*r[1]; y[1]=yc/rc*r[1]; 
                xc=ddX/ddZ*(z[1]-z[0])+x[0];                
                yc=ddY/ddZ*(z[1]-z[0])+y[0];                
                rc=sqrt(xc*xc+yc*yc);
                x[2]=xc/rc*r[2]; y[2]=yc/rc*r[2]; 
                for(n=0; n<cnt; n++)
                  portion[n]=sqrt((delZ[n]*delZ[n]+delR[n]*delR[n])/(ddZ*ddZ+ddR*ddR));
              } else { cnt=0; }
            }
            
            for(n=0; n<cnt; n++) {
              zc=0.5*(z[n]+z[n+1]);
              xc=0.5*(x[n]+x[n+1]);
              yc=0.5*(y[n]+y[n+1]);
              rc=0.5*(r[n]+r[n+1]);
              tmp=sqrt(xc*xc+yc*yc);
              xc=xc/tmp*rc; yc=yc/tmp*rc; 
              indexI1=ii[n]; indexJ1=jj[n];

              cossC[1]=xc/rc; sinsC[1]=yc/rc;
              for(m=2; m<numMode; m++) {
                cossC[m]=cossC[m-1]*cossC[1]-sinsC[m-1]*sinsC[1];
                sinsC[m]=sinsC[m-1]*cossC[1]+cossC[m-1]*sinsC[1];
              }

              Wr1[0]=0.5-capR[n]-0.5*delR[n]; Wr1[1]=1.0-Wr1[0];
              Wz1[0]=0.5-capZ[n]-0.5*delZ[n]; Wz1[1]=1.0-Wz1[0];
//if(Wr1[0]>1.0 || Wr1[0]<0 || Wz1[0]>1.0 || Wz1[0]<0) {
//  printf("mode=%d,cnt=%d,n=%d,Wr1[0]=%g, Wr1[1]=%g,capR[0]=%g,capR[1]=%g,capR[2]=%g,ddR=%g,ddZ=%g,delR[0]=%g,delR[1]=%g,delR[2]=%g\n",mode,cnt,n,Wr1[0],Wr1[1],capR[0],capR[1],capR[2],ddR,ddZ,delR[0],delR[1],delR[2]);
//  printf("mode=%d,cnt=%d,n=%d,Wz1[0]=%g, Wz1[1]=%g,capZ[0]=%g,capZ[1]=%g,capZ[2]=%g,ddR=%g,ddZ=%g,delZ[0]=%g,delZ[1]=%g,delZ[2]=%g\n",mode,cnt,n,Wz1[0],Wz1[1],capZ[0],capZ[1],capZ[2],ddR,ddZ,delZ[0],delZ[1],delZ[2]);
//}
              rr=indexJ1-jstart+0.5;
              factor=weight*coeff[s]/rr;

              for(idI=0; idI<2; idI++) {
//                factor=weight*coeff[s]/(2.0*(indexJ1+idI-jstart));
                J0=delZ[n]*Wr1[idI]*factor*dzBydt;
                D->JzR[0][indexI1][indexJ1+idI]+=J0;
                for(m=1; m<numMode; m++)  {
                  D->JzR[m][indexI1][indexJ1+idI]+=J0*cossC[m]*alpha;
                  D->JzI[m][indexI1][indexJ1+idI]-=J0*sinsC[m]*alpha;
                }
//                factor=weight*coeff[s]/(2.0*(indexJ1-jstart+0.5));
                J0=delR[n]*Wz1[idI]*factor*drBydt;
                D->JrR[0][indexI1+idI][indexJ1]+=J0;
                for(m=1; m<numMode; m++)  {
                  D->JrR[m][indexI1+idI][indexJ1]+=J0*cossC[m]*alpha;
                  D->JrI[m][indexI1+idI][indexJ1]-=J0*sinsC[m]*alpha;
                }
              }
/*
              //J_phi
              indexIC=(int)(zc);  indexJC=(int)(rc)+jstart;
//if(indexIC!=indexI1 || indexJC!=indexJ1) printf("indexIC=%d,indexJC=%d,indexI1=%d,indexJ1=%d,mode=%d,zc=%g,z[%d]=%g,z[%d]=%g\n",indexIC,indexJC,indexI1,indexJ1,mode,zc,n,z[n],n+1,z[n+1]);


              vp=(cossC[1]*p->py-sinsC[1]*p->px)/gamma*portion[n];
              WzC[1]=zc-(int)(zc); WzC[0]=1.0-WzC[1];
              index=indexJC-jstart+1.0;
              WrC[0]=(index*index-rc*rc)/(2.0*index-1);
              WrC[1]=1.0-WrC[0];
              
              for(idI=0; idI<2; idI++)
                for(idJ=0; idJ<2; idJ++) {
//                  factor=weight*coeff[s]/(2.0*(indexJC+idJ-jstart));
                  J0=WzC[idI]*WrC[idJ]*factor*vp;
                  D->JpR[0][indexIC+idI][indexJC+idJ]+=J0;
                  for(m=1; m<numMode; m++)  {
                    D->JpR[m][indexIC+idI][indexJC+idJ]+=J0*cossC[m]*alpha;
                    D->JpI[m][indexIC+idI][indexJC+idJ]-=J0*sinsC[m]*alpha;
                  }
                }
*/
            }  //End of n

            //Jp calculation
            rc=0.5*(r1+r2); zc=0.5*(z1+z2); index=(int)rc;
            xc=x1; yc=y1;           
            idI=(int)zc; idJ=(int)rc+jstart; 
            WzC[1]=zc-idI;     WzC[0]=1.0-WzC[1]; 
            WrC[1]=rc-(int)rc; WrC[0]=1.0-WrC[1]; 
            calculaionRally(&xc,&yc,rc,x1,x2,y1,y2,r1,r2,iteration);
            coss[1]=xc/rc; sins[1]=yc/rc;
            for(m=2; m<numMode; m++) {
              coss[m]=coss[m-1]*coss[1]-sins[m-1]*sins[1];
              sins[m]=sins[m-1]*coss[1]+coss[m-1]*sins[1];
            }
            vp=(coss[1]*p->py-sins[1]*p->px)/gamma;
            factor=vp*weight*coeff[s]/(index+0.5);
            tmpP1=WzC[0]*WrC[0]*factor;
            tmpP2=WzC[1]*WrC[0]*factor;
            tmpP3=WzC[0]*WrC[1]*factor;
            tmpP4=WzC[1]*WrC[1]*factor;
              D->JpR[0][idI][idJ]+=tmpP1;
              D->JpR[0][idI+1][idJ]+=tmpP2;
              D->JpR[0][idI][idJ+1]+=tmpP3;
              D->JpR[0][idI+1][idJ+1]+=tmpP4;
              for(m=1; m<numMode; m++)  {
                D->JpR[m][idI][idJ]+=tmpP1*coss[m]*alpha;
                D->JpR[m][idI+1][idJ]+=tmpP2*coss[m]*alpha;
                D->JpR[m][idI][idJ+1]+=tmpP3*coss[m]*alpha;
                D->JpR[m][idI+1][idJ+1]+=tmpP4*coss[m]*alpha;
                D->JpI[m][idI][idJ]-=tmpP1*sins[m]*alpha;
                D->JpI[m][idI+1][idJ]-=tmpP2*sins[m]*alpha;
                D->JpI[m][idI][idJ+1]-=tmpP3*sins[m]*alpha;
                D->JpI[m][idI+1][idJ+1]-=tmpP4*sins[m]*alpha;
              }

            p=p->next;
          }    //End of while(p)

        }  	 //End of for(s)     
      }  	 //End of for(i,j)

    //At axis
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jstart+2; j++)
      {
        for(s=0; s<nSpecies; s++)
        {
          p=particle[i][j].head[s]->pt;     
          while(p) 
          {
            cnt=0;
            weight=p->weight*p->charge;
            gamma=sqrt(1.0+p->pz*p->pz+p->px*p->px+p->py*p->py);
            x1=p->oldX; y1=p->oldY; r1=sqrt(x1*x1+y1*y1);
            x2=p->x;    y2=p->y;    r2=sqrt(x2*x2+y2*y2);
            z2=p->z+i;  z1=p->oldZ;
            z[0]=z1; x[0]=x1; y[0]=y1; r[0]=r1;

            i1=(int)(z1); j1=(int)(r1)+jstart;
            i2=(int)(z2); j2=(int)(r2)+jstart;    
            ii[0]=i1;       jj[0]=j1;
            difI=i2-i1; difJ=j2-j1;
            ddZ=z2-z1; ddR=r2-r1; ddX=x2-x1; ddY=y2-y1;
            maxI=MAX(i1,i2); maxJ=MAX(j1,j2);
            if(difI==0 && difJ==0) { 
              cnt=1; mode=1; 
              z[1]=z2; x[1]=x2; y[1]=y2; r[1]=r2;
                
              capZ[0]=z[0]-(ii[0]+0.5);
              capR[0]=r[0]+jstart-(jj[0]+0.5);
              delZ[0]=ddZ;   delR[0]=ddR;
              portion[0]=1.0;
            } else if(difI!=0 && difJ==0) { 
              cnt=2; mode=2;
              z[2]=z2; x[2]=x2; y[2]=y2; r[2]=r2;
              ii[1]=i2;       jj[1]=j2;

              capZ[0]=z[0]-(ii[0]+0.5);
              capR[0]=r[0]+jstart-(jj[0]+0.5);
              delZ[0]=maxI-z[0]; 
              delR[0]=ddR*(delZ[0]/ddZ);
              capZ[1]=0.5*(ii[0]+ii[1])-ii[1]; 
              capR[1]=capR[0]+delR[0];
              delZ[1]=ddZ-delZ[0];    delR[1]=ddR-delR[0];

              z[1]=z[0]+delZ[0]; r[1]=r[0]+delR[0];
              xc=ddX/ddZ*(z[1]-z[0])+x[0];
              yc=ddY/ddZ*(z[1]-z[0])+y[0];
              rc=sqrt(xc*xc+yc*yc);
              x[1]=xc/rc*r[1]; y[1]=yc/rc*r[1]; 
              for(n=0; n<cnt; n++)
                portion[n]=sqrt((delZ[n]*delZ[n]+delR[n]*delR[n])/(ddZ*ddZ+ddR*ddR));
            } else if(difI==0 && difJ!=0) { 
              cnt=2; mode=3;
              z[2]=z2; x[2]=x2; y[2]=y2; r[2]=r2;
              ii[1]=i2;       jj[1]=j2;

              capZ[0]=z[0]-(ii[0]+0.5);
              capR[0]=r[0]+jstart-(jj[0]+0.5);
              delR[0]=maxJ-jstart-r[0]; 
              delZ[0]=ddZ*(delR[0]/ddR);
              capR[1]=0.5*(jj[0]+jj[1])-jj[1]; 
              capZ[1]=capZ[0]+delZ[0];
              delR[1]=ddR-delR[0];    delZ[1]=ddZ-delZ[0];

              z[1]=z[0]+delZ[0]; r[1]=r[0]+delR[0];
              xc=ddX/ddZ*(z[1]-z[0])+x[0];
              yc=ddY/ddZ*(z[1]-z[0])+y[0];
              rc=sqrt(xc*xc+yc*yc);
              x[1]=xc/rc*r[1]; y[1]=yc/rc*r[1]; 
              if(ddZ==0.0) cnt=0; else ;
              for(n=0; n<cnt; n++)
                portion[n]=sqrt((delZ[n]*delZ[n]+delR[n]*delR[n])/(ddZ*ddZ+ddR*ddR));
            } else { 
              if(fabs(ddR/ddZ)<fabs((maxJ-jstart-r[0])/(maxI-z[0]))) {
                cnt=3; mode=4;
                z[3]=z2; x[3]=x2; y[3]=y2; r[3]=r2;
                ii[2]=i2;       jj[2]=j2;

                capZ[0]=z[0]-(ii[0]+0.5);
                capR[0]=r[0]+jstart-(jj[0]+0.5);
                z[1]=maxI;       delZ[0]=z[1]-z[0]; 
                ii[1]=ii[2];     capZ[1]=0.5*(ii[0]+ii[1])-ii[1]; 
                delR[0]=ddR*delZ[0]/ddZ; r[1]=r[0]+delR[0];
                jj[1]=jj[0];     capR[1]=r[1]+jstart-maxJ;
                delR[1]=maxJ-jstart-r[1]; delZ[1]=ddZ*(delR[1]/ddR);
                z[2]=z[1]+delZ[1];    capZ[2]=z[2]-(ii[2]+0.5); 
                r[2]=maxJ-jstart; capR[2]=0.5*(jj[0]+jj[2])-jj[2];
                delZ[2]=ddZ-delZ[0]-delZ[1];
                delR[2]=ddR-delR[0]-delR[1];

                xc=ddX/ddZ*(z[1]-z[0])+x[0];                
                yc=ddY/ddZ*(z[1]-z[0])+y[0]; 
                rc=sqrt(xc*xc+yc*yc);
                x[1]=xc/rc*r[1]; y[1]=yc/rc*r[1]; 
                xc=ddX/ddZ*(z[2]-z[0])+x[0];                
                yc=ddY/ddZ*(z[2]-z[0])+y[0];                
                rc=sqrt(xc*xc+yc*yc);
                x[2]=xc/rc*r[2]; y[2]=yc/rc*r[2]; 
                for(n=0; n<cnt; n++)
                  portion[n]=sqrt((delZ[n]*delZ[n]+delR[n]*delR[n])/(ddZ*ddZ+ddR*ddR));
              } else if(fabs(ddR/ddZ)>fabs((maxJ-jstart-r[0])/(maxI-z[0]))) {
                cnt=3; mode=5;
                z[3]=z2; x[3]=x2; y[3]=y2; r[3]=r2;
                ii[2]=i2;       jj[2]=j2;

                capZ[0]=z[0]-(ii[0]+0.5);
                capR[0]=r[0]+jstart-(jj[0]+0.5);
                r[1]=maxJ-jstart; delR[0]=r[1]-r[0]; 
                jj[1]=jj[2];      capR[1]=0.5*(jj[1]+jj[1])-jj[1];
                delZ[0]=ddZ*delR[0]/ddR; z[1]=z[0]+delZ[0];
                ii[1]=ii[0];      capZ[1]=z[1]-maxI;
                delZ[1]=maxI-z[1]; delR[1]=ddR*(delZ[1]/ddZ);
                r[2]=r[1]+delR[1]; capR[2]=r[2]+jstart-(jj[2]+0.5); 
                z[2]=maxI;         capZ[2]=0.5*(ii[0]+ii[2])-ii[2];
                delZ[2]=ddZ-delZ[0]-delZ[1];
                delR[2]=ddR-delR[0]-delR[1];

                xc=ddX/ddZ*(z[1]-z[0])+x[0];                
                yc=ddY/ddZ*(z[1]-z[0])+y[0];                
                rc=sqrt(xc*xc+yc*yc);
                x[1]=xc/rc*r[1]; y[1]=yc/rc*r[1]; 
                xc=ddX/ddZ*(z[1]-z[0])+x[0];                
                yc=ddY/ddZ*(z[1]-z[0])+y[0];                
                rc=sqrt(xc*xc+yc*yc);
                x[2]=xc/rc*r[2]; y[2]=yc/rc*r[2]; 
                for(n=0; n<cnt; n++)
                  portion[n]=sqrt((delZ[n]*delZ[n]+delR[n]*delR[n])/(ddZ*ddZ+ddR*ddR));
              } else { cnt=0; }
            }

            for(n=0; n<cnt; n++) {
              zc=0.5*(z[n]+z[n+1]);
              xc=0.5*(x[n]+x[n+1]);
              yc=0.5*(y[n]+y[n+1]);
              rc=0.5*(r[n]+r[n+1]);
              tmp=sqrt(xc*xc+yc*yc);
              xc=xc/tmp*rc; yc=yc/tmp*rc; 
              indexI1=ii[n]; indexJ1=jj[n];

              cossC[1]=xc/rc; sinsC[1]=yc/rc;
              for(m=2; m<numMode; m++) {
                cossC[m]=cossC[m-1]*cossC[1]-sinsC[m-1]*sinsC[1];
                sinsC[m]=sinsC[m-1]*cossC[1]+cossC[m-1]*sinsC[1];
              }

              Wr1[0]=0.5-capR[n]-0.5*delR[n]; Wr1[1]=1.0-Wr1[0];
              Wz1[0]=0.5-capZ[n]-0.5*delZ[n]; Wz1[1]=1.0-Wz1[0];
//if(Wr1[0]>1.0 || Wr1[0]<0 || Wz1[0]>1.0 || Wz1[0]<0) {
//  printf("mode=%d,cnt=%d,n=%d,Wr1[0]=%g, Wr1[1]=%g,capR[0]=%g,capR[1]=%g,capR[2]=%g,ddR=%g,ddZ=%g,delR[0]=%g,delR[1]=%g,delR[2]=%g\n",mode,cnt,n,Wr1[0],Wr1[1],capR[0],capR[1],capR[2],ddR,ddZ,delR[0],delR[1],delR[2]);
//  printf("mode=%d,cnt=%d,n=%d,Wz1[0]=%g, Wz1[1]=%g,capZ[0]=%g,capZ[1]=%g,capZ[2]=%g,ddR=%g,ddZ=%g,delZ[0]=%g,delZ[1]=%g,delZ[2]=%g\n",mode,cnt,n,Wz1[0],Wz1[1],capZ[0],capZ[1],capZ[2],ddR,ddZ,delZ[0],delZ[1],delZ[2]);
//}
              rr=indexJ1-jstart+0.5;
              factor=weight*coeff[s]/rr;

              if(indexJ1==jstart) {
//                if(rc<0.5)  factor=weight*coeff[s]/(rc*rc+rc+0.25); 
//                else        factor=weight*coeff[s];
                J0=delZ[n]*Wr1[0]*factor*dzBydt;
                D->JzR[0][indexI1][indexJ1]+=J0;

//                factor=weight*coeff[s]/(2.0*(indexJ1+1-jstart));
                J0=delZ[n]*Wr1[1]*factor*dzBydt;
                D->JzR[0][indexI1][indexJ1+1]+=J0;
                for(m=1; m<numMode; m++)  {
                  D->JzR[m][indexI1][indexJ1+1]+=J0*cossC[m]*alpha;
                  D->JzI[m][indexI1][indexJ1+1]-=J0*sinsC[m]*alpha;
                }

//                factor=weight*coeff[s]/(2.0*(indexJ1-jstart+0.5));
                for(idI=0; idI<2; idI++) {
                  J0=delR[n]*Wz1[idI]*factor*drBydt;
                  D->JrR[0][indexI1+idI][indexJ1]+=J0;
                  for(m=1; m<numMode; m++)  {
                    D->JrR[m][indexI1+idI][indexJ1]+=J0*cossC[m]*alpha;
                    D->JrI[m][indexI1+idI][indexJ1]-=J0*sinsC[m]*alpha;
                  }
                }
              } else {
//                factor=weight*coeff[s]/(2.0*(indexJ1-jstart+0.5));
                for(idI=0; idI<2; idI++) {
//                  factor=weight*coeff[s]/(2.0*(indexJ1+idI-jstart));
                  J0=delZ[n]*Wr1[idI]*factor*dzBydt;
                  D->JzR[0][indexI1][indexJ1+idI]+=J0;
                  for(m=1; m<numMode; m++)  {
                    D->JzR[m][indexI1][indexJ1+idI]+=J0*cossC[m]*alpha;
                    D->JzI[m][indexI1][indexJ1+idI]-=J0*sinsC[m]*alpha;
                  }

//                  factor=weight*coeff[s]/(2.0*(indexJ1-jstart+0.5));
                  J0=delR[n]*Wz1[idI]*factor*drBydt;
                  D->JrR[0][indexI1+idI][indexJ1]+=J0;
                  for(m=1; m<numMode; m++)  {
                    D->JrR[m][indexI1+idI][indexJ1]+=J0*cossC[m]*alpha;
                    D->JrI[m][indexI1+idI][indexJ1]-=J0*sinsC[m]*alpha;
                  }
                }
              }
/*
              //J_phi
              indexIC=(int)(zc);  indexJC=(int)(rc)+jstart;
//if(indexIC!=indexI1 || indexJC!=indexJ1) printf("indexIC=%d,indexJC=%d,indexI1=%d,indexJ1=%d,mode=%d,zc=%g,z[%d]=%g,z[%d]=%g\n",indexIC,indexJC,indexI1,indexJ1,mode,zc,n,z[n],n+1,z[n+1]);
              vp=(cossC[1]*p->py-sinsC[1]*p->px)/gamma*portion[n];
              WzC[1]=zc-(int)(zc); WzC[0]=1.0-WzC[1];
              index=indexJC-jstart+1.0;
              WrC[0]=(index*index-rc*rc)/(2.0*index-1);
              WrC[1]=1.0-WrC[0];

              if(indexJC==jstart) { 
//                if(rc<0.5)  factor=weight*coeff[s]/(rc*rc+rc+0.25); else;
//                factor=weight*coeff[s]/(2.0*(indexJC-jstart+0.5));
                for(idI=0; idI<2; idI++) {
                  //bottom
//                  factor=weight*coeff[s]*4.0;
                  J0=WzC[idI]*WrC[0]*factor*vp;
                  m=1;
                  D->JpR[m][indexIC+idI][indexJC]+=J0*cossC[m]*alpha;
                  D->JpI[m][indexIC+idI][indexJC]-=J0*sinsC[m]*alpha;
                  //up
//                  factor=weight*coeff[s]/(2.0*(indexJC+1-jstart));
                  J0=WzC[idI]*WrC[1]*factor*vp;
                  m=0;
                  D->JpR[m][indexIC+idI][indexJC+1]+=J0;
                  for(m=1; m<numMode; m++)  {
                    D->JpR[m][indexIC+idI][indexJC+1]+=J0*cossC[m]*alpha;
                    D->JpI[m][indexIC+idI][indexJC+1]-=J0*sinsC[m]*alpha;
                  }
                }
              } else {
//                factor=weight*coeff[s]/(2.0*(indexJC-jstart+0.5));
                for(idI=0; idI<2; idI++)
                  for(idJ=0; idJ<2; idJ++) {
//                    factor=weight*coeff[s]/(2.0*(indexJC+idJ-jstart));
                    J0=WzC[idI]*WrC[idJ]*factor*vp;
                    D->JpR[0][indexIC+idI][indexJC+idJ]+=J0;
                    for(m=1; m<numMode; m++)  {
                      D->JpR[m][indexIC+idI][indexJC+idJ]+=J0*cossC[m]*alpha;
                      D->JpI[m][indexIC+idI][indexJC+idJ]-=J0*sinsC[m]*alpha;
                    }
                  }
              }
*/
            }  //End of n

            //Jp calculation
            rc=0.5*(r1+r2); zc=0.5*(z1+z2); index=(int)rc;
            xc=x1; yc=y1;            
            idI=(int)zc; idJ=(int)rc+jstart; 
            WzC[1]=zc-idI;     WzC[0]=1.0-WzC[1]; 
            WrC[1]=rc-(int)rc; WrC[0]=1.0-WrC[1]; 
            calculaionRally(&xc,&yc,rc,x1,x2,y1,y2,r1,r2,iteration);
            coss[1]=xc/rc; sins[1]=yc/rc;
            for(m=2; m<numMode; m++) {
              coss[m]=coss[m-1]*coss[1]-sins[m-1]*sins[1];
              sins[m]=sins[m-1]*coss[1]+coss[m-1]*sins[1];
            }
            vp=(coss[1]*p->py-sins[1]*p->px)/gamma;
            factor=vp*weight*coeff[s]/(index+0.5);
            tmpP1=WzC[0]*WrC[0]*factor;
            tmpP2=WzC[1]*WrC[0]*factor;
            tmpP3=WzC[0]*WrC[1]*factor;
            tmpP4=WzC[1]*WrC[1]*factor;
            if(index==0) {
              D->JpR[0][idI][idJ+1]+=tmpP3;
              D->JpR[0][idI+1][idJ+1]+=tmpP4;
              m=1;
                D->JpR[m][idI][idJ]+=tmpP1*coss[m]*alpha;
                D->JpR[m][idI+1][idJ]+=tmpP2*coss[m]*alpha;
                D->JpR[m][idI][idJ+1]+=tmpP3*coss[m]*alpha;
                D->JpR[m][idI+1][idJ+1]+=tmpP4*coss[m]*alpha;
                D->JpI[m][idI][idJ]-=tmpP1*sins[m]*alpha;
                D->JpI[m][idI+1][idJ]-=tmpP2*sins[m]*alpha;
                D->JpI[m][idI][idJ+1]-=tmpP3*sins[m]*alpha;
                D->JpI[m][idI+1][idJ+1]-=tmpP4*sins[m]*alpha;
              for(m=2; m<numMode; m++)  {
                D->JpR[m][idI][idJ+1]+=tmpP3*coss[m]*alpha;
                D->JpR[m][idI+1][idJ+1]+=tmpP4*coss[m]*alpha;
                D->JpI[m][idI][idJ+1]-=tmpP3*sins[m]*alpha;
                D->JpI[m][idI+1][idJ+1]-=tmpP4*sins[m]*alpha;
              }
            } else {
              D->JpR[0][idI][idJ]+=tmpP1;
              D->JpR[0][idI+1][idJ]+=tmpP2;
              D->JpR[0][idI][idJ+1]+=tmpP3;
              D->JpR[0][idI+1][idJ+1]+=tmpP4;
              for(m=1; m<numMode; m++)  {
                D->JpR[m][idI][idJ]+=tmpP1*coss[m]*alpha;
                D->JpR[m][idI+1][idJ]+=tmpP2*coss[m]*alpha;
                D->JpR[m][idI][idJ+1]+=tmpP3*coss[m]*alpha;
                D->JpR[m][idI+1][idJ+1]+=tmpP4*coss[m]*alpha;
                D->JpI[m][idI][idJ]-=tmpP1*sins[m]*alpha;
                D->JpI[m][idI+1][idJ]-=tmpP2*sins[m]*alpha;
                D->JpI[m][idI][idJ+1]-=tmpP3*sins[m]*alpha;
                D->JpI[m][idI+1][idJ+1]-=tmpP4*sins[m]*alpha;
              }
            }

            p=p->next;
          }    //End of while(p)

        }  	 //End of for(s)     
      }  	 //End of for(i,j)

}



void findPhaseDif(int m,double  cosP1,double cosP2,double sinP1,double sinP2,double s1,double s2,double *A,double *B,double *X,double *Y,double *coss,double *sins)
{
  switch (m) {
    case 1 :
      *A=0.5*(sqrt((1.0+cosP2)*(1.0+cosP1))+s1*s2*sqrt((1.0-cosP2)*(1.0-cosP1)));
      *B=0.5*(s2*sqrt((1.0-cosP2)*(1.0+cosP1))-s1*sqrt((1.0+cosP2)*(1.0-cosP1)));
      *coss=*A;
      *sins=*B;
      break;
    case 2 :
      *X=cosP2*cosP1+sinP2*sinP1;
      *Y=sinP2*cosP1-cosP2*sinP1;
      *coss=*X;
      *sins=*Y;
      break;
    case 3 :
      *coss=(*A)*(*X)-(*B)*(*Y);
      *sins=(*B)*(*X)+(*A)*(*Y);
      break;
    case 4 :
      *coss=(*X)*(*X)-(*Y)*(*Y);
      *sins=2.0*(*X)*(*Y);
      break;
    default :
      printf("out of mode. Check number_mode.\n");
  }

}

void findPhaseAdd(int m,double  cosP1,double cosP2,double sinP1,double sinP2,double s1,double s2,double *A,double *B,double *X,double *Y,double *coss,double *sins)
{
  switch (m) {
    case 1 :
      *A=0.5*(sqrt((1.0+cosP2)*(1.0+cosP1))-s1*s2*sqrt((1.0-cosP2)*(1.0-cosP1)));
      *B=0.5*(s2*sqrt((1.0-cosP2)*(1.0+cosP1))+s1*sqrt((1.0+cosP2)*(1.0-cosP1)));
      *coss=*A;
      *sins=*B;
      break;
    case 2 :
      *X=cosP2*cosP1-sinP2*sinP1;
      *Y=sinP2*cosP1+cosP2*sinP1;
      *coss=*X;
      *sins=*Y;
      break;
    case 3 :
      *coss=(*A)*(*X)-(*B)*(*Y);
      *sins=(*B)*(*X)+(*A)*(*Y);
      break;
    case 4 :
      *coss=(*X)*(*X)-(*Y)*(*Y);
      *sins=2.0*(*X)*(*Y);
      break;
    default :
      printf("out of mode. Check number_mode.\n");
  }

}

double findR(double x1, double x2,double x3, double x4)
{
//  double minimum();
//  double maximum();
  double result,result1,result2,result3;

  result1=MIN(x1-0.5,x2-0.5);
  result2=MAX(x1-1.5,x2-1.5);
  result3=MAX(result2,(x3+x4)*0.5);
  result=MIN(result1,result3);

  return result;
}


int intmaximum(int x1,int x2)
{
   int result;

   if(x1>=x2)
      result=x1;
   else
      result=x2;
  
   return result;
}

int intminimum(int x1,int x2)
{
   int result;

   if(x1>=x2)
      result=x2;
   else
      result=x1;
  
   return result;
}
