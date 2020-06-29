#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_qrng.h>


double maxwellianVelocity(double temperature);
int randomInt(int range);
void loadPlasma_crystal(Domain *D,LoadList *LL,int s);
void random1D_sobol(double *x,gsl_qrng *q);
void random2D_sobol(double *x,double *y,gsl_qrng *q);
void random3D_sobol(double *z,double *x,double *y,gsl_qrng *q);
double randomValue(double beta);
void loadPolygonPlasma(Domain *D,LoadList *LL,int s,int iteration,int istart,int iend,int jstart,int jend);

double applyFunctionX(int mode,double centerX,double x,double gaussCoefX,double polyCoefX)
{
  double result;

  switch (mode)  {
  case 0 :	//Costant
    result=1.0;
    break;
  case 1 :	//Gaussian
    result=exp(-(x-centerX)*(x-centerX)/gaussCoefX/gaussCoefX);
    break;
  case 2 :	//2nd polynomial
    result=1.0+polyCoefX*(x-centerX)*(x-centerX);
    break;
  }
  return result;
}

double applyFunctionYZ(int mode,double centerY,double y,double centerZ,double z,double gaussCoefYZ,double polyCoefYZ)
{
  double result;

  switch (mode)  {
  case 0 :	//Constant
    result=1.0;
    break;
  case 1 :	//Gaussian
    result=exp(-((y-centerY)*(y-centerY)+(z-centerZ)*(z-centerZ))/gaussCoefYZ/gaussCoefYZ);
    break;
  case 2 :	//2nd polynomial
    result=1.0+polyCoefYZ*((y-centerY)*(y-centerY)+(z-centerZ)*(z-centerZ));
    break;
  }
  return result;
}

void loadPlasma(Domain *D,LoadList *LL,int s,int iteration,int istart,int iend,int jstart,int jend)
{
  void loadChannelPlasma();

  switch(LL->type)  {
  case Polygon:
    loadPolygonPlasma(D,LL,s,iteration,istart,iend,jstart,jend); 
    break;
//  case ((Polygon-1)*3+3):
//    loadPolygonPlasma3D(D,LL,s,iteration); 
    break;

  default:
    ;
  }
}


void loadPolygonPlasma(Domain *D,LoadList *LL,int s,int iteration,int istart,int iend,int jstart,int jend)
{
   int n,i,j,numberRZ,numPhi,cnt,l,t,index;
   int modeX,modeYZ,minRSub,minZSub;
   double z,r,R,posX,posY,posZ,v1,v2,v3,centerX,centerY,centerZ,tmp,weight,charge;
   double gaussCoefX,polyCoefX,gaussCoefYZ,polyCoefYZ,minZ,maxZ;
   double tz,tx,ty,positionR;
   double ne,randTest,dPhi,phi,z0,pz0,alpha,cosPhi,sinPhi;
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   Particle **particle;
   particle=D->particle;

   ptclList *New,*p;   

   minZSub=D->minXSub;   minRSub=D->minYSub;
   numPhi=LL->numberPhi; numberRZ=LL->numberRZ;
   pz0=LL->pz0;
   z0=LL->z0;
   alpha=LL->delPz/LL->delZ;
   minZ=z0-LL->delZ*0.5;
   maxZ=z0+LL->delZ*0.5;

   centerX=LL->centerX;
   centerY=LL->centerY;
   centerZ=0.0;
   gaussCoefX=LL->gaussCoefX;
   polyCoefX=LL->polyCoefX;
   gaussCoefYZ=LL->gaussCoefYZ;
   polyCoefYZ=LL->polyCoefYZ;
   modeX=LL->modeX;
   modeYZ=LL->modeYZ;
   charge=LL->charge;


   srand((iteration+1)*myrank);

/*
   //position define      
   gsl_qrng *q = gsl_qrng_alloc (gsl_qrng_sobol,2);
   for(i=istart; i<iend; i++)
     for(j=jstart+2; j<jend; j++)
     {
       for(l=0; l<LL->xnodes-1; l++)
         for(t=0; t<LL->ynodes-1; t++)
         {
           posX=(double)(i+D->minXSub-istart);
           posY=(double)(j+D->minYSub-jstart);
           posZ=0.0;
           if(posX>=LL->xpoint[l] && posX<LL->xpoint[l+1] &&
              posY>=LL->ypoint[t] && posY<LL->ypoint[t+1])
           {
             ne=((LL->xn[l+1]-LL->xn[l])/(LL->xpoint[l+1]-LL->xpoint[l])*(posX-LL->xpoint[l])+LL->xn[l]);
             ne*=((LL->yn[t+1]-LL->yn[t])/(LL->ypoint[t+1]-LL->ypoint[t])*(posY-LL->ypoint[t])+LL->yn[t]);
//             tmp=applyFunctionX(modeX,centerX,posX,gaussCoefX,polyCoefX);
//             ne*=tmp;
             tmp=applyFunctionYZ(modeYZ,centerY,posY,centerZ,posZ,gaussCoefYZ,polyCoefYZ);
             ne*=tmp;
//             ne*=LL->numberRZ;	//it is the double number of superparticles.
//             intNum=(int)ne;
//             randTest=ne-intNum;
//             if(randTest>randomValue(1.0)) intNum=intNum+1;
             weight=ne/(double)(numPhi*numberRZ);
             dPhi=2.0*pi/((double)numPhi);
 
             for(n=0; n<numberRZ; n++) {
//					r=2.0;
//               while(r>=1.0 || r==0.0) {
                 random2D_sobol(&tz,&positionR,q);
//                 tz=randomValue(1.0);
//                 tx=randomValue(1.0);
//                 ty=randomValue(1.0);
//                 tx=2.0*tx-1.0;
//                 ty=2.0*ty-1.0;
//                 r=sqrt(tx*tx+ty*ty);
//               }

               cnt=0;
               phi=randomValue(1.0)*2.0*pi;             
//               phi=0.5*dPhi;
               while(cnt<numPhi)  {      
//                 positionZ=randomValue(1.0);
//                 positionR=randomValue(1.0);

                 New = (ptclList *)malloc(sizeof(ptclList)); 
                 New->next = particle[i][j].head[s]->pt;
                 particle[i][j].head[s]->pt = New;
 
                 New->z = tz;
                 New->oldZ= i+tz;
                 index=j-jstart+minRSub;
                 cosPhi=cos(phi);
                 sinPhi=sin(phi);
                   r=positionR+j-jstart;
                   New->x=r*cosPhi;
                   New->y=r*sinPhi;
//                   New->x=(index+r)*tx/r;
//                   New->y=(index+r)*ty/r;
                   New->oldX=New->x;
                   New->oldY=New->y;
                   New->weight=weight*(2.0*index+1.0);
                 New->charge=charge;

                 New->Ez=New->Ex=New->Ey=0.0;
                 New->Bz=New->Bx=New->By=0.0;
                 
                 z=i+tz-istart+minZSub;
                 if(z>minZ && z<maxZ) {
                   v1=alpha*(z-z0)+pz0;
                   v1=maxwellianVelocity(LL->temperatureZ)/velocityC+v1;
                   v2=maxwellianVelocity(LL->temperatureR)/velocityC;
                   v3=maxwellianVelocity(LL->temperatureR)/velocityC;
                 } else {
                   v1=maxwellianVelocity(LL->temperatureZ)/velocityC;
                   v2=maxwellianVelocity(LL->temperatureR)/velocityC;
                   v3=maxwellianVelocity(LL->temperatureR)/velocityC;
                 }
                 New->pz=v1;      New->px=v2;       New->py=v3;
                 LL->index+=1;
                 New->index=LL->index;            
                 New->core=myrank;            

                 cnt++; 
                 phi+=dPhi;
               }		//end of while(cnt)
             }			//end of for(n)
           }	
         } 		//end of for(lnodes)  

     }			//End of for(i,j)
     gsl_qrng_free(q);
*/

   gsl_qrng *q = gsl_qrng_alloc (gsl_qrng_sobol,3);
   for(i=istart; i<iend; i++)
     for(j=jstart; j<jend; j++)
     {
       for(l=0; l<LL->xnodes-1; l++)
         for(t=0; t<LL->ynodes-1; t++)
         {
           posX=(double)(i+D->minXSub-istart);
           posY=(double)(j+D->minYSub-jstart);
           posZ=0.0;
           if(posX>=LL->xpoint[l] && posX<LL->xpoint[l+1] &&
              posY>=LL->ypoint[t] && posY<LL->ypoint[t+1])
           {
             ne=((LL->xn[l+1]-LL->xn[l])/(LL->xpoint[l+1]-LL->xpoint[l])*(posX-LL->xpoint[l])+LL->xn[l]);
             ne*=((LL->yn[t+1]-LL->yn[t])/(LL->ypoint[t+1]-LL->ypoint[t])*(posY-LL->ypoint[t])+LL->yn[t]);
//             tmp=applyFunctionX(modeX,centerX,posX,gaussCoefX,polyCoefX);
//             ne*=tmp;
             tmp=applyFunctionYZ(modeYZ,centerY,posY,centerZ,posZ,gaussCoefYZ,polyCoefYZ);
             ne*=tmp;
//             ne*=LL->numberRZ;	//it is the double number of superparticles.
//             intNum=(int)ne;
//             randTest=ne-intNum;
//             if(randTest>randomValue(1.0)) intNum=intNum+1;
             weight=ne/(double)(numPhi*numberRZ);
 
             for(n=0; n<numberRZ; n++) {
					r=2.0;
               while(r>=1.0 || r==0.0) {
                 tz=randomValue(1.0);
                 tx=randomValue(1.0);
                 ty=randomValue(1.0);
                 tx=2.0*tx-1.0;
                 ty=2.0*ty-1.0;
                 r=sqrt(tx*tx+ty*ty);
               }

               cnt=0;
               while(cnt<numPhi)  {      
                 New = (ptclList *)malloc(sizeof(ptclList)); 
                 New->next = particle[i][j].head[s]->pt;
                 particle[i][j].head[s]->pt = New;
 
                 New->z = tz;
                 New->oldZ= i+tz;
                 index=j-jstart+minRSub;
                 New->x=(index+r)*tx/r;
                 New->y=(index+r)*ty/r;
                 New->oldX=New->x;
                 New->oldY=New->y;
                 New->weight=weight*(2.0*index+1.0);
                 New->charge=charge;

                 New->Ez=New->Ex=New->Ey=0.0;
                 New->Bz=New->Bx=New->By=0.0;
                 
                 z=i+tz-istart+minZSub;
                 if(z>minZ && z<maxZ) {
                   v1=alpha*(z-z0)+pz0;
                   v1=maxwellianVelocity(LL->temperatureZ)/velocityC+v1;
                   v2=maxwellianVelocity(LL->temperatureR)/velocityC;
                   v3=maxwellianVelocity(LL->temperatureR)/velocityC;
                 } else {
                   v1=maxwellianVelocity(LL->temperatureZ)/velocityC;
                   v2=maxwellianVelocity(LL->temperatureR)/velocityC;
                   v3=maxwellianVelocity(LL->temperatureR)/velocityC;
                 }
                 New->pz=v1;      New->px=v2;       New->py=v3;
                 LL->index+=1;
                 New->index=LL->index;            
                 New->core=myrank;            

                 cnt++; 
               }		//end of while(cnt)
             }			//end of for(n)
           }	
         } 		//end of for(lnodes)  

     }			//End of for(i,j)
     gsl_qrng_free(q);
}

double maxwellianVelocity(double temperature)
{
   double vth,r,prob,v,random;
   int intRand,randRange=1e5;

   vth=sqrt(2.0*eCharge*temperature/eMass);
   
   r=1.0;
   prob=0.0;
   while (r>prob)  {
      intRand = rand() % randRange;
      r = ((double)intRand)/randRange;
      intRand = rand() % randRange;
      random = ((double)intRand)/randRange;
      v = 6.0*(random-0.5);
      prob=exp(-v*v);
   }
   return vth*v;
}

int randomInt(int range)
{
   int r;

   r = rand() % range;

   return r;
}

void random1D_sobol(double *x,gsl_qrng *q)
{
   double v[1];

   gsl_qrng_get(q,v);
   *x=v[0];
}

void random2D_sobol(double *x,double *y,gsl_qrng *q)
{
   double v[2];

   gsl_qrng_get(q,v);
   *x=v[0];
   *y=v[1];
}

void random3D_sobol(double *z,double *x,double *y,gsl_qrng *q)
{
   double v[3];

   gsl_qrng_get(q,v);
   *z=v[0];
   *x=v[1];
   *y=v[2];
}

