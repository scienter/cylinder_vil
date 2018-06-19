#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "mesh.h"
#include "constants.h"
#include "math.h"

void absorb_R(Domain *D,double *rtr,double *rtd,double x,double rL,double LdR,double rr,double rd);
void absorb_L(Domain *D,double *lftr,double *lftd,double x,double lL,double LdL,double rr,double rd);
void absorb_U(Domain *D,double *upr,double *upd,double y,double upL,double LdU,double rr,double rd);
void Bsolve2D_Yee(Domain *D);
void Esolve2D_Yee(Domain *D,double dF,int iteration);
void Bsolve2D_NoCherenkov(Domain *D,int iteration);
void solve_NDFX_C(Domain *D);
void solve_NDFX(Domain *D);

void fieldSolve1(Domain D,double t,int iteration,double dF)
{
  int rankX,rankY;
  float limit;
  LaserList *L;
  int myrank, nTasks,rank,rankM,rankN;

  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  switch(D.fieldType) {

  case Yee :
    //load laser
    if(D.boostOn==OFF)   {
      L=D.laserList;
      while(L->next)  {
        limit=((L->rU+L->rD)/D.dtRatio*2.1+L->retard)*1.0/D.divisionLambda;
        if(iteration<=limit) loadLaser(&D,L,t); else ;
        L=L->next;
      }
    } else ;

    Bsolve2D_Yee(&D);
    if(D.L>1)  {
      MPI_Transfer12F_Xminus(&D,D.BzR,D.BrR,D.BpR,D.BzI,D.BrI,D.BpI,D.BzNowR,D.BrNowR,D.BpNowR,D.BzNowI,D.BrNowI,D.BpNowI,D.nySub+5,3);
      MPI_Transfer12F_Xplus(&D,D.BzR,D.BrR,D.BpR,D.BzI,D.BrI,D.BpI,D.BzNowR,D.BrNowR,D.BpNowR,D.BzNowI,D.BrNowI,D.BpNowI,D.nySub+5,3);
    } else	;

/*
    if(D.pmlLeft==ON && iteration>D.pmlStart) {
      if(myrank==0) PML_Bsolve2D_Left(&D); else ;
      MPI_Barrier(MPI_COMM_WORLD);
      if(myrank==0) copyPML_Left_B(&D,D.pmlCellLeft); else ;
    } else ;
*/

    break ;

  case NoCherenkov :
    Bsolve2D_NoCherenkov(&D,iteration);
    if(D.L>1)  {
      MPI_Transfer12F_Xminus(&D,D.BzR,D.BrR,D.BpR,D.BzI,D.BrI,D.BpI,D.BzNowR,D.BrNowR,D.BpNowR,D.BzNowI,D.BrNowI,D.BpNowI,D.nySub+5,3);
      MPI_Transfer12F_Xplus(&D,D.BzR,D.BrR,D.BpR,D.BzI,D.BrI,D.BpI,D.BzNowR,D.BrNowR,D.BpNowR,D.BzNowI,D.BrNowI,D.BpNowI,D.nySub+5,3);
    } else	;
    if(D.M>1)  {
//    MPI_Transfer6F_Yminus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nxSub+5,1,3);
//    MPI_Transfer6F_Yplus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nxSub+5,1,3);
    } else	;

    //load laser
    if(D.boostOn==OFF)   {
      L=D.laserList;
      while(L->next)  {
        limit=((L->rU+L->rD)*2.1+L->retard)*1.0/D.divisionLambda;
        if(iteration<=limit) loadLaser(&D,L,t); else ;
        L=L->next;
      }
    } else ;
    break ;

  case NDFX :
     ;
    break ;
  }
}

void fieldSolve2(Domain D,double t,int iteration,double dF)
{
  int rankX,rankY;
  float limit;
  LaserList *L;
  int myrank, nTasks,rank,rankM,rankN;

  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  switch(D.fieldType) {
  case Yee :
    Esolve2D_Yee(&D,dF,iteration);
    if(D.L>1)  {
      MPI_Transfer6F_Xminus(&D,D.EzR,D.ErR,D.EpR,D.EzI,D.ErI,D.EpI,D.nySub+5,3);
      MPI_Transfer6F_Xplus(&D,D.EzR,D.ErR,D.EpR,D.EzI,D.ErI,D.EpI,D.nySub+5,3);
    } else	;

    break ;

  case NoCherenkov :
    Esolve2D_Yee(&D,dF,iteration);
    if(D.L>1)  {
      MPI_Transfer6F_Xminus(&D,D.EzR,D.ErR,D.EpR,D.EzI,D.ErI,D.EpI,D.nySub+5,3);
      MPI_Transfer6F_Xplus(&D,D.EzR,D.ErR,D.EpR,D.EzI,D.ErI,D.EpI,D.nySub+5,3);
    } else	;
    break ;

  case NDFX :
    //load laser
    if(D.boostOn==OFF)   {
      L=D.laserList;
      while(L->next)  {
        limit=((L->rU+L->rD)/D.dtRatio*2.1+L->retard)*1.0/D.divisionLambda;
        if(iteration<=limit) loadLaser(&D,L,t); else ;
        L=L->next;
      }
    } else ;

    solve_NDFX_C(&D);
    if(D.L>1)  {
      MPI_Transfer12F_Xplus(&D,D.PrCR,D.PrCI,D.SrCR,D.SrCI,D.PlCR,D.PlCI,D.SlCR,D.SlCI,D.EzCR,D.EzCI,D.BzCR,D.BzCI,D.nySub+5,3);
      MPI_Transfer12F_Xminus(&D,D.PrCR,D.PrCI,D.SrCR,D.SrCI,D.PlCR,D.PlCI,D.SlCR,D.SlCI,D.EzCR,D.EzCI,D.BzCR,D.BzCI,D.nySub+5,3);
    } else	;

    solve_NDFX(&D);
    if(D.L>1)  {
      MPI_Transfer12F_Xminus(&D,D.PrR,D.PrI,D.SrR,D.SrI,D.PlR,D.PlI,D.SlR,D.SlI,D.EzR,D.BzR,D.EzI,D.BzI,D.nySub+5,3);
      MPI_Transfer12F_Xplus(&D,D.PrR,D.PrI,D.SrR,D.SrI,D.PlR,D.PlI,D.SlR,D.SlI,D.EzR,D.BzR,D.EzI,D.BzI,D.nySub+5,3);
    } else	;
    break ;
  }
}

void Bsolve2D_NoCherenkov(Domain *D,int iteration)
{
  int i,j,m,numMode,istart,iend,jstart,jend;  
  double x,minZSub,dtBydr,dtBydz,r,dr,dz,dt,BzR,BrR,BpR,BzI,BrI,BpI;
  double oldBzR,oldBrR,oldBpR,oldBzI,oldBrI,oldBpI;
  double deltaZ,betaPZ,betaRZ,alphaZ,alphaP,alphaR;
  double upr,upd,lftr,lftd,upL,leftL,LdU,LdL,rr,rd,tmp,tmpr,tmpd;

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  numMode=D->numMode;
  minZSub=D->minXSub;

  dtBydr=D->dt/D->dr;
  dtBydz=D->dt/D->dz;
  dr=D->dr; dz=D->dz; dt=D->dt;

  leftL=(double)(D->minXDomain+D->pmlCellLeft);
  upL=(double)(D->ny+D->minYDomain-D->pmlCellUp);
  LdU=D->pmlCellUp;
  LdL=D->pmlCellLeft;
  rr=D->pmlr;    rd=D->pmld;
  lftr=lftd=upr=upd=1.0;

  tmp=sin(0.5*pi*dtBydz);
  deltaZ=0.25*(1.0-tmp*tmp/dtBydz/dtBydz);
  betaPZ=0.25;
  betaRZ=0.25;
  alphaZ=1.0-3.0*deltaZ;
  alphaP=1.0-2.0*betaPZ;
  alphaR=1.0-2.0*betaRZ;

  for(m=0; m<numMode; m++)
    for(i=istart; i<iend; i++)
    {
      x=(i-istart)+minZSub;
      if(D->pmlLeft==ON && iteration>D->pmlStart)
          absorb_L(D,&lftr,&lftd,x,leftL,LdL,rr,rd);
      else ;
      for(j=jstart+1; j<jend; j++)  
      {
        r=(double)(j-jstart);
        if(D->pmlUp==ON)
          absorb_U(D,&upr,&upd,r,upL,LdU,rr,rd);
        else    ;
        oldBzR=D->BzR[m][i][j];
        oldBrR=D->BrR[m][i][j];
        oldBpR=D->BpR[m][i][j];
        oldBzI=D->BzI[m][i][j];
        oldBrI=D->BrI[m][i][j];
        oldBpI=D->BpI[m][i][j];

        tmpr=upr*lftr;
        tmpd=upd*lftd;
        tmp=tmpr*dtBydr/r*(-alphaR*((r+0.5)*D->EpR[m][i][j+1]-(r-0.5)*D->EpR[m][i][j])-betaRZ*((r+0.5)*D->EpR[m][i+1][j+1]-(r-0.5)*D->EpR[m][i+1][j]+(r+0.5)*D->EpR[m][i-1][j+1]-(r-0.5)*D->EpR[m][i-1][j])-m*(alphaP*D->ErI[m][i][j]+betaPZ*(D->ErI[m][i+1][j]+D->ErI[m][i-1][j])));
        D->BzR[m][i][j]=tmpd*(oldBzR+tmp);
        tmp=tmpr*dtBydr/r*(-alphaR*((r+0.5)*D->EpI[m][i][j+1]-(r-0.5)*D->EpI[m][i][j])-betaRZ*((r+0.5)*D->EpI[m][i+1][j+1]-(r-0.5)*D->EpI[m][i+1][j]+(r+0.5)*D->EpI[m][i-1][j+1]-(r-0.5)*D->EpI[m][i-1][j])+m*(alphaP*D->ErR[m][i][j]+betaPZ*(D->ErR[m][i+1][j]+D->ErR[m][i-1][j])));
        D->BzI[m][i][j]=tmpd*(oldBzI+tmp);

        tmp=tmpr*(dtBydz*(alphaZ*(D->EpR[m][i+1][j]-D->EpR[m][i][j])+deltaZ*(D->EpR[m][i+2][j]-D->EpR[m][i-1][j]))+m*dtBydr/(r-0.5)*(alphaP*D->EzI[m][i][j]+betaPZ*(D->EzI[m][i+1][j]+D->EzI[m][i-1][j])));
        D->BrR[m][i][j]=tmpd*(oldBrR+tmp);
        tmp=tmpr*(dtBydz*(alphaZ*(D->EpI[m][i+1][j]-D->EpI[m][i][j])+deltaZ*(D->EpI[m][i+2][j]-D->EpI[m][i-1][j]))-m*dtBydr/(r-0.5)*(alphaP*D->EzR[m][i][j]+betaPZ*(D->EzR[m][i+1][j]+D->EzR[m][i-1][j])));
        D->BrI[m][i][j]=tmpd*(oldBrI+tmp);

        tmp=tmpr*(dtBydr*(alphaR*(D->EzR[m][i][j+1]-D->EzR[m][i][j])+betaRZ*(D->EzR[m][i+1][j+1]-D->EzR[m][i+1][j]+D->EzR[m][i-1][j+1]-D->EzR[m][i-1][j]))-dtBydz*(alphaZ*(D->ErR[m][i+1][j]-D->ErR[m][i][j])+deltaZ*(D->ErR[m][i+2][j]-D->ErR[m][i-1][j])));
        D->BpR[m][i][j]=tmpd*(oldBpR+tmp);
        tmp=tmpr*(dtBydr*(alphaR*(D->EzI[m][i][j+1]-D->EzI[m][i][j])+betaRZ*(D->EzI[m][i+1][j+1]-D->EzI[m][i+1][j]+D->EzI[m][i-1][j+1]-D->EzI[m][i-1][j]))-dtBydz*(alphaZ*(D->ErI[m][i+1][j]-D->ErI[m][i][j])+deltaZ*(D->ErI[m][i+2][j]-D->ErI[m][i-1][j])));
        D->BpI[m][i][j]=tmpd*(oldBpI+tmp);

        D->BzNowR[m][i][j]=0.5*(D->BzR[m][i][j]+oldBzR);
        D->BrNowR[m][i][j]=0.5*(D->BrR[m][i][j]+oldBrR);
        D->BpNowR[m][i][j]=0.5*(D->BpR[m][i][j]+oldBpR);
        D->BzNowI[m][i][j]=0.5*(D->BzI[m][i][j]+oldBzI);
        D->BrNowI[m][i][j]=0.5*(D->BrI[m][i][j]+oldBrI);
        D->BpNowI[m][i][j]=0.5*(D->BpI[m][i][j]+oldBpI);
      }
    }

  lftr=lftd=upr=upd=1.0;
  j=jstart;
    for(i=istart; i<iend; i++) {
      x=(i-istart)+minZSub;
      if(D->pmlLeft==ON && iteration>D->pmlStart)
          absorb_L(D,&lftr,&lftd,x,leftL,LdL,rr,rd);
      else ;

      oldBzR=D->BzR[0][i][j];
      oldBpR=D->BpR[1][i][j];
      oldBpI=D->BpI[1][i][j];

      tmpr=upr*lftr;
      tmpd=upd*lftd;

      D->BzR[0][i][j]+=-4.0*dtBydr*D->EpR[0][i][j+1];

      tmp=tmpr*(dt*(4.0*D->EzR[1][i][j+1]-D->EzR[1][i][j+2])-dtBydz*(alphaZ*(D->ErR[1][i+1][j]-D->ErR[1][i][j])+deltaZ*(D->ErR[1][i+2][j]-D->ErR[1][i-1][j])));
      D->BpR[1][i][j]=tmpd*(oldBpR+tmp);
      tmp=tmpr*(dt*(4.0*D->EzI[1][i][j+1]-D->EzI[1][i][j+2])-dtBydz*(alphaZ*(D->ErI[1][i+1][j]-D->ErI[1][i][j])+deltaZ*(D->ErI[1][i+2][j]-D->ErI[1][i-1][j])));
      D->BpI[1][i][j]=tmpd*(oldBpI+tmp);
      D->BzNowR[0][i][j]=0.5*(D->BzR[0][i][j]+oldBzR);
      D->BpNowR[1][i][j]=0.5*(D->BpR[1][i][j]+oldBpR);
      D->BpNowI[1][i][j]=0.5*(D->BpI[1][i][j]+oldBpI);
    }
}

void Bsolve2D_Yee(Domain *D)
{
  int i,j,m,numMode,istart,iend,jstart,jend;  
  double dtBydr,dtBydz,r,dr,BzR,BrR,BpR,BzI,BrI,BpI;
  double oldBzR,oldBrR,oldBpR,oldBzI,oldBrI,oldBpI;
  double upr,upd,upL,LdU,rr,rd,tmp,tmpr,tmpd,coef1,coef2,coef3,coef4;

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  numMode=D->numMode;

  upL=(double)(D->ny+D->minYDomain-D->pmlCellUp);
  LdU=D->pmlCellUp;
  rr=D->pmlr;    rd=D->pmld;
  upr=upd=1.0;

  dtBydr=D->dt/D->dr;
  dtBydz=D->dt/D->dz;
  dr=D->dr;
  for(m=0; m<numMode; m++)
    for(i=istart; i<iend; i++)
    {
      for(j=jstart+1; j<jend; j++)  
      {
        r=(double)(j-jstart);
        if(D->pmlUp==ON)
          absorb_U(D,&upr,&upd,r,upL,LdU,rr,rd);
        else    ;

        oldBzR=D->BzR[m][i][j];
        oldBrR=D->BrR[m][i][j];
        oldBpR=D->BpR[m][i][j];
        oldBzI=D->BzI[m][i][j];
        oldBrI=D->BrI[m][i][j];
        oldBpI=D->BpI[m][i][j];

        tmpr=upr;    tmpd=upd;

        tmp=tmpr*(-dtBydr*m/(r+0.5)*D->ErI[m][i][j]
             -0.5*dtBydr/(r+0.5)*(D->EpR[m][i][j+1]+D->EpR[m][i][j])
             -dtBydr*(D->EpR[m][i][j+1]-D->EpR[m][i][j]));
        D->BzR[m][i][j]=tmpd*(oldBzR+tmp);

        tmp=tmpr*(dtBydr*m/(r+0.5)*D->ErR[m][i][j]
             -0.5*dtBydr/(r+0.5)*(D->EpI[m][i][j+1]+D->EpI[m][i][j])
             -dtBydr*(D->EpI[m][i][j+1]-D->EpI[m][i][j]));
        D->BzI[m][i][j]=tmpd*(oldBzI+tmp);

        tmp=tmpr*(dtBydz*(D->EpR[m][i+1][j]-D->EpR[m][i][j])
              +dtBydr*m/r*D->EzI[m][i][j]);
        D->BrR[m][i][j]=tmpd*(oldBrR+tmp);
        tmp=tmpr*(dtBydz*(D->EpI[m][i+1][j]-D->EpI[m][i][j])
              -dtBydr*m/r*D->EzR[m][i][j]);
        D->BrI[m][i][j]=tmpd*(oldBrI+tmp);

        tmp=tmpr*(dtBydr*(D->EzR[m][i][j+1]-D->EzR[m][i][j])
             -dtBydz*(D->ErR[m][i+1][j]-D->ErR[m][i][j]));
        D->BpR[m][i][j]=tmpd*(oldBpR+tmp);
        tmp=tmpr*(dtBydr*(D->EzI[m][i][j+1]-D->EzI[m][i][j])
             -dtBydz*(D->ErI[m][i+1][j]-D->ErI[m][i][j]));
        D->BpI[m][i][j]=tmpd*(oldBpI+tmp);

        D->BzNowR[m][i][j]=0.5*(D->BzR[m][i][j]+oldBzR);
        D->BrNowR[m][i][j]=0.5*(D->BrR[m][i][j]+oldBrR);
        D->BpNowR[m][i][j]=0.5*(D->BpR[m][i][j]+oldBpR);
        D->BzNowI[m][i][j]=0.5*(D->BzI[m][i][j]+oldBzI);
        D->BrNowI[m][i][j]=0.5*(D->BrI[m][i][j]+oldBrI);
        D->BpNowI[m][i][j]=0.5*(D->BpI[m][i][j]+oldBpI);
      }
    }

  j=jstart; m=1;
    for(i=istart; i<iend; i++) {
      oldBrR=D->BrR[m][i][j];
      oldBrI=D->BrI[m][i][j];
      D->BrR[m][i][j]+=dtBydz*(D->EpR[m][i+1][j]-D->EpR[m][i][j])
        +dtBydr*D->EzI[m][i][j+1];
      D->BrI[m][i][j]+=dtBydz*(D->EpI[m][i+1][j]-D->EpI[m][i][j])
        -dtBydr*D->EzR[m][i][j+1];
      D->BrNowR[m][i][j]=0.5*(D->BrR[m][i][j]+oldBrR);
      D->BrNowI[m][i][j]=0.5*(D->BrI[m][i][j]+oldBrI);
    }
  r=0.0;
  for(m=0; m<numMode; m++)
    for(i=istart; i<iend; i++)
    {
//        if(D->pmlUp==ON)
//          absorb_U(D,&upr,&upd,r,upL,LdU,rr,rd);
//        else    ;
        oldBzR=D->BzR[m][i][j];
        oldBpR=D->BpR[m][i][j];
        oldBzI=D->BzI[m][i][j];
        oldBpI=D->BpI[m][i][j];

        D->BzR[m][i][j]+=-dtBydr*m/(r+0.5)*D->ErI[m][i][j]
          -0.5*dtBydr/(r+0.5)*(D->EpR[m][i][j+1]+D->EpR[m][i][j])
          -dtBydr*(D->EpR[m][i][j+1]-D->EpR[m][i][j]);
        D->BzI[m][i][j]+=dtBydr*m/(r+0.5)*D->ErR[m][i][j]
          -0.5*dtBydr/(r+0.5)*(D->EpI[m][i][j+1]+D->EpI[m][i][j])
          -dtBydr*(D->EpI[m][i][j+1]-D->EpI[m][i][j]);

        D->BpR[m][i][j]+=dtBydr*(D->EzR[m][i][j+1]-D->EzR[m][i][j])
          -dtBydz*(D->ErR[m][i+1][j]-D->ErR[m][i][j]);
        D->BpI[m][i][j]+=dtBydr*(D->EzI[m][i][j+1]-D->EzI[m][i][j])
          -dtBydz*(D->ErI[m][i+1][j]-D->ErI[m][i][j]);

        D->BzNowR[m][i][j]=0.5*(D->BzR[m][i][j]+oldBzR);
        D->BpNowR[m][i][j]=0.5*(D->BpR[m][i][j]+oldBpR);
        D->BzNowI[m][i][j]=0.5*(D->BzI[m][i][j]+oldBzI);
        D->BpNowI[m][i][j]=0.5*(D->BpI[m][i][j]+oldBpI);
    }
}

//lala
void Esolve2D_Yee(Domain *D,double dF,int iteration)
{
  int i,j,m,numMode,istart,iend,jstart,jend,a;  
  double dtBydr,dtBydz,r,dt,dr,dz;
  double EzR,ErR,EpR,EzI,ErI,EpI,x,minZSub;
  double upr,upd,lftr,lftd,upL,leftL,LdU,LdL,rr,rd,tmp,tmpr,tmpd,alpha[2];
  double oldEzR,oldEzI,oldErR,oldErI,oldEpR,oldEpI;
  double beforeEpR,beforeEpI,beforeEzR,beforeEzI,beforeErR,beforeErI;
  double nowEpR,nowEpI,nowEzR,nowEzI,nowErR,nowErI;
  double beforeBpR,beforeBpI,beforeBzR,beforeBzI,beforeBrR,beforeBrI;
  double nowBpR,nowBpI,nowBzR,nowBzI,nowBrR,nowBrI;

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  numMode=D->numMode;
  minZSub=D->minXSub;

  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  leftL=(double)(D->minXDomain+D->pmlCellLeft);
  upL=(double)(D->ny+D->minYDomain-D->pmlCellUp);
  LdU=D->pmlCellUp;
  LdL=D->pmlCellLeft;
  rr=D->pmlr;    rd=D->pmld;
  lftr=lftd=upr=upd=1.0;

  dtBydr=D->dt/D->dr; dr=D->dr; dz=D->dz;
  dtBydz=D->dt/D->dz; dt=D->dt;
  m=0;
    for(i=istart; i<iend; i++) 
    {
//      x=(i-istart)+minZSub;
//      if(D->pmlLeft==ON && iteration>D->pmlStart)  
//        absorb_L(D,&lftr,&lftd,x,leftL,LdL,rr,rd);
//      else   ;
      for(j=jstart+1; j<jend; j++) 
      {
        r=(double)(j-jstart);
        if(D->pmlUp==ON) absorb_U(D,&upr,&upd,r,upL,LdU,rr,rd);
        else    ;

        oldEzR=D->EzR[m][i][j];
        oldEzI=D->EzI[m][i][j];
        oldErR=D->ErR[m][i][j];
        oldErI=D->ErI[m][i][j];
        oldEpR=D->EpR[m][i][j];
        oldEpI=D->EpI[m][i][j];

        tmpr=upr*lftr;     tmpd=upd*lftd;

        tmp=tmpr*(dtBydr*m/r*D->BrI[m][i][j]
              +0.5*dtBydr/r*(D->BpR[m][i][j]+D->BpR[m][i][j-1])
              +dtBydr*(D->BpR[m][i][j]-D->BpR[m][i][j-1])
              -2.0*pi*dt*D->JzR[m][i][j]
              +dF*dtBydz*(D->FR[m][i+1][j]-D->FR[m][i][j]));
        D->EzR[m][i][j]=tmpd*(oldEzR+tmp);
        tmp=tmpr*(-dtBydr*m/r*D->BrR[m][i][j]
              +0.5*dtBydr/r*(D->BpI[m][i][j]+D->BpI[m][i][j-1])
              +dtBydr*(D->BpI[m][i][j]-D->BpI[m][i][j-1])
              -2.0*pi*dt*D->JzI[m][i][j]
              +dF*dtBydz*(D->FI[m][i+1][j]-D->FI[m][i][j]));
        D->EzI[m][i][j]=tmpd*(oldEzI+tmp);

        tmp=tmpr*(-2.0*pi*dt*D->JrR[m][i][j]
              -dtBydz*(D->BpR[m][i][j]-D->BpR[m][i-1][j])
              -dtBydr*m/(r+0.5)*D->BzI[m][i][j]
              +dF*dtBydr*(D->FR[m][i][j+1]-D->FR[m][i][j]));
        D->ErR[m][i][j]=tmpd*(oldErR+tmp);
        tmp=tmpr*(-2.0*pi*dt*D->JrI[m][i][j]
              -dtBydz*(D->BpI[m][i][j]-D->BpI[m][i-1][j])
              +dtBydr*m/(r+0.5)*D->BzR[m][i][j]
              +dF*dtBydr*(D->FI[m][i][j+1]-D->FI[m][i][j]));
        D->ErI[m][i][j]=tmpd*(oldErI+tmp);
        
        tmp=tmpr*(-2.0*pi*dt*D->JpR[m][i][j]
              -dtBydr*(D->BzR[m][i][j]-D->BzR[m][i][j-1])
              +dtBydz*(D->BrR[m][i][j]-D->BrR[m][i-1][j]));
        D->EpR[m][i][j]=tmpd*(oldEpR+tmp);
        tmp=tmpr*(-2.0*pi*dt*D->JpI[m][i][j]
              -dtBydr*(D->BzI[m][i][j]-D->BzI[m][i][j-1])
              +dtBydz*(D->BrI[m][i][j]-D->BrI[m][i-1][j]));
        D->EpI[m][i][j]=tmpd*(oldEpI+tmp);
      }
    }
  for(m=1; m<numMode; m++) 
    for(i=istart; i<iend; i++) 
    {
//      x=(i-istart)+minZSub;
//      if(D->pmlLeft==ON && iteration>D->pmlStart)  
//        absorb_L(D,&lftr,&lftd,x,leftL,LdL,rr,rd);
//      else   ;
      for(j=jstart+1; j<jend; j++) 
      {
        r=(double)(j-jstart);
        if(D->pmlUp==ON) absorb_U(D,&upr,&upd,r,upL,LdU,rr,rd);
        else    ;

        oldEzR=D->EzR[m][i][j];
        oldEzI=D->EzI[m][i][j];
        oldErR=D->ErR[m][i][j];
        oldErI=D->ErI[m][i][j];
        oldEpR=D->EpR[m][i][j];
        oldEpI=D->EpI[m][i][j];

        tmpr=upr*lftr;     tmpd=upd*lftd;

        tmp=tmpr*(dtBydr*m/r*D->BrI[m][i][j]
              +0.5*dtBydr/r*(D->BpR[m][i][j]+D->BpR[m][i][j-1])
              +dtBydr*(D->BpR[m][i][j]-D->BpR[m][i][j-1])
              -2.0*pi*dt*D->JzR[m][i][j]
              +dF*dtBydz*(D->FR[m][i+1][j]-D->FR[m][i][j]));
        D->EzR[m][i][j]=tmpd*(oldEzR+tmp);
        tmp=tmpr*(-dtBydr*m/r*D->BrR[m][i][j]
              +0.5*dtBydr/r*(D->BpI[m][i][j]+D->BpI[m][i][j-1])
              +dtBydr*(D->BpI[m][i][j]-D->BpI[m][i][j-1])
              -2.0*pi*dt*D->JzI[m][i][j]
              +dF*dtBydz*(D->FI[m][i+1][j]-D->FI[m][i][j]));
        D->EzI[m][i][j]=tmpd*(oldEzI+tmp);
       
        tmp=tmpr*(-2.0*pi*dt*D->JrR[m][i][j]
              -dtBydz*(D->BpR[m][i][j]-D->BpR[m][i-1][j])
              -dtBydr*m/(r+0.5)*D->BzI[m][i][j]
              +dF*dtBydr*(D->FR[m][i][j+1]-D->FR[m][i][j]));
        D->ErR[m][i][j]=tmpd*(oldErR+tmp);
        tmp=tmpr*(-2.0*pi*dt*D->JrI[m][i][j]
              -dtBydz*(D->BpI[m][i][j]-D->BpI[m][i-1][j])
              +dtBydr*m/(r+0.5)*D->BzR[m][i][j]
              +dF*dtBydr*(D->FI[m][i][j+1]-D->FI[m][i][j]));
        D->ErI[m][i][j]=tmpd*(oldErI+tmp);
        
        tmp=tmpr*(-2.0*pi*dt*D->JpR[m][i][j]
              -dtBydr*(D->BzR[m][i][j]-D->BzR[m][i][j-1])
              +dtBydz*(D->BrR[m][i][j]-D->BrR[m][i-1][j])
              -dF*dtBydr/r*m*D->FI[m][i][j]);
        D->EpR[m][i][j]=tmpd*(oldEpR+tmp);
        tmp=tmpr*(-2.0*pi*dt*D->JpI[m][i][j]
              -dtBydr*(D->BzI[m][i][j]-D->BzI[m][i][j-1])
              +dtBydz*(D->BrI[m][i][j]-D->BrI[m][i-1][j])
              +dF*dtBydr/r*m*D->FR[m][i][j]);
        D->EpI[m][i][j]=tmpd*(oldEpI+tmp);
      }
    }

  j=jstart;
    for(i=istart; i<iend; i++) {
      m=1;
      D->EpR[m][i][j]+=-2.0*pi*dt*D->JpR[m][i][j]
        -2.0*dtBydr*(D->BzR[m][i][j])
        +dtBydz*(D->BrR[m][i][j]-D->BrR[m][i-1][j])
        -dF*dtBydr*m*D->FI[m][i][j+1];
      D->EpI[m][i][j]+=-2.0*pi*dt*D->JpI[m][i][j]
        -2.0*dtBydr*(D->BzI[m][i][j])
        +dtBydz*(D->BrI[m][i][j]-D->BrI[m][i-1][j])
        +dF*dtBydr*m*D->FR[m][i][j+1];
      m=0;
      D->EzR[m][i][j]+=4.0*dtBydr*D->BpR[m][i][j]
        -2.0*pi*dt*D->JzR[m][i][j]
        +dF*dtBydz*(D->FR[m][i+1][j]-D->FR[m][i][j]);
    }
  for(m=0; m<numMode; m++) 
    for(i=istart; i<iend; i++) 
    {
        r=(double)(j-jstart);
        oldErR=D->ErR[m][i][j];
        oldErI=D->ErI[m][i][j];

        D->ErR[m][i][j]+=-2.0*pi*dt*D->JrR[m][i][j]
          -dtBydz*(D->BpR[m][i][j]-D->BpR[m][i-1][j])
          -dtBydr*m/(r+0.5)*D->BzI[m][i][j]
          +dF*dtBydr*(D->FR[m][i][j+1]-D->FR[m][i][j]);
        D->ErI[m][i][j]+=-2.0*pi*dt*D->JrI[m][i][j]
          -dtBydz*(D->BpI[m][i][j]-D->BpI[m][i-1][j])
          +dtBydr*m/(r+0.5)*D->BzR[m][i][j] 
          +dF*dtBydr*(D->FI[m][i][j+1]-D->FI[m][i][j]);
    }
}

void solve_NDFX_C(Domain *D)
{
  int i,j,m,numMode,istart,iend,jstart,jend;  
  double dtBydr,dtBydz,r,dr,dz,dt,coefEB1,coefEB2,coefEB3,coefRL;
  double nowPrCR,nowSrCR,prevPrCR,prevSrCR;
  double nowPrCI,nowSrCI,prevPrCI,prevSrCI;
  double upPrR,upPlR,upSrR,upSlR;
  double dnPrR,dnPlR,dnSrR,dnSlR,PrR,PlR,SrR,SlR;
  double upPrI,upPlI,upSrI,upSlI;
  double dnPrI,dnPlI,dnSrI,dnSlI,PrI,PlI,SrI,SlI;
  double rr,rd,upL,LdU,tmpr,tmpd,tmp,upr,upd;
  double oldEzCR,oldBzCR,oldPrCR,oldPlCR,oldSrCR,oldSlCR;
  double oldEzCI,oldBzCI,oldPrCI,oldPlCI,oldSrCI,oldSlCI;

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  numMode=D->numMode;

  dtBydr=D->dt/D->dr;
  dtBydz=D->dt/D->dz;
  dr=D->dr; dz=D->dz; dt=D->dt;
  upL=(double)(D->ny+D->minYDomain-D->pmlCellUp);
  LdU=D->pmlCellUp;
  rr=D->pmlr;    rd=D->pmld;

  upr=upd=1.0;

  for(m=0; m<numMode; m++) 
    for(j=jstart+1; j<jend; j++) {
      nowPrCR=D->PrCR[m][istart-1][j];
      nowSrCR=D->SrCR[m][istart-1][j];
      nowPrCI=D->PrCI[m][istart-1][j];
      nowSrCI=D->SrCI[m][istart-1][j];
      r=(double)(j-jstart);
      if(D->pmlUp==ON) absorb_U(D,&upr,&upd,r,upL,LdU,rr,rd);
      else    ;
      for(i=istart; i<iend; i++) {
        upPrR=0.5*(D->PrR[m][i-1][j]  +D->PrR[m][i][j]);
        dnPrR=0.5*(D->PrR[m][i-1][j-1]+D->PrR[m][i][j-1]);
        upPrI=0.5*(D->PrI[m][i-1][j]  +D->PrI[m][i][j]);
        dnPrI=0.5*(D->PrI[m][i-1][j-1]+D->PrI[m][i][j-1]);
        upSrR=0.5*(D->SrR[m][i-1][j]  +D->SrR[m][i][j]);
        dnSrR=0.5*(D->SrR[m][i-1][j-1]+D->SrR[m][i][j-1]);
        upSrI=0.5*(D->SrI[m][i-1][j]  +D->SrI[m][i][j]);
        dnSrI=0.5*(D->SrI[m][i-1][j-1]+D->SrI[m][i][j-1]);
        upPlR=0.5*(D->PlR[m][i-1][j]  +D->PlR[m][i][j]);
        dnPlR=0.5*(D->PlR[m][i-1][j-1]+D->PlR[m][i][j-1]);
        upPlI=0.5*(D->PlI[m][i-1][j]  +D->PlI[m][i][j]);
        dnPlI=0.5*(D->PlI[m][i-1][j-1]+D->PlI[m][i][j-1]);
        upSlR=0.5*(D->SlR[m][i-1][j]  +D->SlR[m][i][j]);
        dnSlR=0.5*(D->SlR[m][i-1][j-1]+D->SlR[m][i][j-1]);
        upSlI=0.5*(D->SlI[m][i-1][j]  +D->SlI[m][i][j]);
        dnSlI=0.5*(D->SlI[m][i-1][j-1]+D->SlI[m][i][j-1]);
        PrR=0.5*(upPrR+dnPrR); PlR=0.5*(upPlR+dnPlR); 
        PrI=0.5*(upPrI+dnPrI); PlI=0.5*(upPlI+dnPlI); 
        SrR=0.5*(upSrR+dnSrR); SlR=0.5*(upSlR+dnSlR); 
        SrI=0.5*(upSrI+dnSrI); SlI=0.5*(upSlI+dnSlI); 

        oldEzCR=D->EzCR[m][i][j]; oldEzCI=D->EzCI[m][i][j];
        oldBzCR=D->BzCR[m][i][j]; oldBzCI=D->BzCI[m][i][j];
        oldPrCR=D->PrCR[m][i][j]; oldPrCI=D->PrCI[m][i][j];
        oldPlCR=D->PlCR[m][i][j]; oldPlCI=D->PlCI[m][i][j];
        oldSrCR=D->SrCR[m][i][j]; oldSrCI=D->SrCI[m][i][j];
        oldSlCR=D->SlCR[m][i][j]; oldSlCI=D->SlCI[m][i][j];

        tmpr=upr;     tmpd=upd;

        tmp=tmpr*(0.5/r*dtBydr*(PrR-PlR)
              +0.5/r*m*dtBydr*(SlI-SrI)
              +0.5*dtBydr*((upPrR-dnPrR)-(upPlR-dnPlR))
              -0.5*pi*dt*(D->JzR[m][i-1][j]+D->JzR[m][i][j]+D->JzCR[m][i-1][j]+D->JzCR[m][i][j]));
        D->EzCR[m][i][j]=tmpd*(oldEzCR+tmp);
        tmp=tmpr*(0.5/r*dtBydr*(PrI-PlI)
              -0.5/r*m*dtBydr*(SlR-SrR)
              +0.5*dtBydr*((upPrI-dnPrI)-(upPlI-dnPlI))
              -0.5*pi*dt*(D->JzI[m][i-1][j]+D->JzI[m][i][j]+D->JzCI[m][i-1][j]+D->JzCI[m][i][j]));
        D->EzCI[m][i][j]=tmpd*(oldEzCI+tmp);

        tmp=tmpr*(-0.5/r*dtBydr*(SrR+SlR)
              -0.5/r*m*dtBydr*(PlI+PrI)
              -0.5*dtBydr*((upSrR-dnSrR)+(upSlR-dnSlR)));
        D->BzCR[m][i][j]=tmpd*(oldBzCR+tmp);
        tmp=tmpr*(-0.5/r*dtBydr*(SrI+SlI)
              +0.5/r*m*dtBydr*(PlR+PrR)
              -0.5*dtBydr*((upSrI-dnSrI)+(upSlI-dnSlI)));
        D->BzCI[m][i][j]=tmpd*(oldBzCI+tmp);

        prevPrCR=nowPrCR;
        nowPrCR=D->PrCR[m][i][j];
        tmp=tmpr*(-0.5/(r+0.5)*m*dtBydr*(D->BzI[m][i][j]+D->BzI[m][i][j+1])
              +dtBydr*(D->EzR[m][i][j+1]-D->EzR[m][i][j])
              -pi*dt*(D->JrR[m][i][j]+D->JrCR[m][i][j]));
        D->PrCR[m][i][j]=tmpd*(prevPrCR+tmp);
        prevPrCI=nowPrCI;
        nowPrCI=D->PrCI[m][i][j];
        tmp=tmpr*(0.5/(r+0.5)*m*dtBydr*(D->BzR[m][i][j]+D->BzR[m][i][j+1])
              +dtBydr*(D->EzI[m][i][j+1]-D->EzI[m][i][j])
              -pi*dt*(D->JrI[m][i][j]+D->JrCI[m][i][j]));
        D->PrCI[m][i][j]=tmpd*(prevPrCI+tmp);

        prevSrCR=nowSrCR;
        nowSrCR=D->SrCR[m][i][j];
        tmp=tmpr*(-0.5/(r+0.5)*m*dtBydr*(D->EzI[m][i][j]+D->EzI[m][i][j+1])
              -dtBydr*(D->BzR[m][i][j+1]-D->BzR[m][i][j])
              -0.5*pi*dt*(D->JpR[m][i][j]+D->JpR[m][i][j+1]+D->JpCR[m][i][j]+D->JpCR[m][i][j+1]));
        D->SrCR[m][i][j]=tmpd*(prevSrCR+tmp);
        prevSrCI=nowSrCI;
        nowSrCI=D->SrCI[m][i][j];
        tmp=tmpr*(0.5/(r+0.5)*m*dtBydr*(D->EzR[m][i][j]+D->EzR[m][i][j+1])
              -dtBydr*(D->BzI[m][i][j+1]-D->BzI[m][i][j])
              -0.5*pi*dt*(D->JpI[m][i][j]+D->JpI[m][i][j+1]+D->JpCI[m][i][j]+D->JpCI[m][i][j+1]));
        D->SrCI[m][i][j]=tmpd*(prevSrCI+tmp);

        tmp=tmpr*(-0.5/(r+0.5)*m*dtBydr*(D->BzI[m][i+1][j+1]+D->BzI[m][i+1][j])
              -dtBydr*(D->EzR[m][i+1][j+1]-D->EzR[m][i+1][j])
              -pi*dt*(D->JrR[m][i+1][j]+D->JrCR[m][i+1][j]));
        D->PlCR[m][i][j]=tmpd*(D->PlCR[m][i+1][j]+tmp);
        tmp=tmpr*(0.5/(r+0.5)*m*dtBydr*(D->BzR[m][i+1][j+1]+D->BzR[m][i+1][j])
              -dtBydr*(D->EzI[m][i+1][j+1]-D->EzI[m][i+1][j])
              -pi*dt*(D->JrI[m][i+1][j]+D->JrCI[m][i+1][j]));
        D->PlCI[m][i][j]=tmpd*(D->PlCI[m][i+1][j]+tmp);

        tmp=tmpr*(0.5/(r+0.5)*m*dtBydr*(D->EzI[m][i+1][j+1]+D->EzI[m][i+1][j])
              -dtBydr*(D->BzR[m][i+1][j+1]-D->BzR[m][i+1][j])
              -0.5*pi*dt*(D->JpR[m][i+1][j]+D->JpR[m][i+1][j+1]+D->JpCR[m][i+1][j]+D->JpCR[m][i+1][j+1]));
        D->SlCR[m][i][j]=tmpd*(D->SlCR[m][i+1][j]+tmp);
        tmp=tmpr*(-0.5/(r+0.5)*m*dtBydr*(D->EzR[m][i+1][j+1]+D->EzR[m][i+1][j])
              -dtBydr*(D->BzI[m][i+1][j+1]-D->BzI[m][i+1][j])
              -0.5*pi*dt*(D->JpI[m][i+1][j]+D->JpI[m][i+1][j+1]+D->JpCI[m][i+1][j]+D->JpCI[m][i+1][j+1]));
        D->SlCI[m][i][j]=tmpd*(D->SlCR[m][i+1][j]+tmp);
      }       //End of i
    }         //End of j
  
  j=jstart; m=0;
  for(i=istart; i<iend; i++) {
    D->BzCR[m][i][j]+=
      -dtBydr*(D->SrR[m][i-1][j]+D->SrR[m][i][j]+D->SlR[m][i-1][j]+D->SlR[m][i][j]);
    D->EzCR[m][i][j]+=
       dtBydr*(D->PrR[m][i-1][j]+D->PrR[m][i][j]-D->PlR[m][i-1][j]-D->PlR[m][i][j])
      -0.5*pi*dt*(D->JzCR[m][i-1][j]+D->JzCR[m][i][j]+D->JzR[m][i-1][j]+D->JzR[m][i][j]);
  }
  for(m=0; m<numMode; m++) { 
      nowPrCR=D->PrCR[m][istart-1][j];
      nowSrCR=D->SrCR[m][istart-1][j];
      nowPrCI=D->PrCI[m][istart-1][j];
      nowSrCI=D->SrCI[m][istart-1][j];
      r=(double)(j-jstart);
      for(i=istart; i<iend; i++) {
        prevPrCR=nowPrCR;
        nowPrCR=D->PrCR[m][i][j];
        D->PrCR[m][i][j]=prevPrCR
          -0.5/(r+0.5)*m*dtBydr*(D->BzI[m][i][j]+D->BzI[m][i][j+1])
          +dtBydr*(D->EzR[m][i][j+1]-D->EzR[m][i][j])
          -pi*dt*(D->JrR[m][i][j]+D->JrCR[m][i][j]);
        prevPrCI=nowPrCI;
        nowPrCI=D->PrCI[m][i][j];
        D->PrCI[m][i][j]=prevPrCI
          +0.5/(r+0.5)*m*dtBydr*(D->BzR[m][i][j]+D->BzR[m][i][j+1])
          +dtBydr*(D->EzI[m][i][j+1]-D->EzI[m][i][j])
          -pi*dt*(D->JrI[m][i][j]+D->JrCI[m][i][j]);
        prevSrCR=nowSrCR;
        nowSrCR=D->SrCR[m][i][j];
        D->SrCR[m][i][j]=prevSrCR
          -0.5/(r+0.5)*m*dtBydr*(D->EzI[m][i][j]+D->EzI[m][i][j+1])
          -dtBydr*(D->BzR[m][i][j+1]-D->BzR[m][i][j])
          -0.5*pi*dt*(D->JpR[m][i][j]+D->JpR[m][i][j+1]+D->JpCR[m][i][j]+D->JpCR[m][i][j+1]);
        prevSrCI=nowSrCI;
        nowSrCI=D->SrCI[m][i][j];
        D->SrCI[m][i][j]=prevSrCI
          +0.5/(r+0.5)*m*dtBydr*(D->EzR[m][i][j]+D->EzR[m][i][j+1])
          -dtBydr*(D->BzI[m][i][j+1]-D->BzI[m][i][j])
          -0.5*pi*dt*(D->JpI[m][i][j]+D->JpI[m][i][j+1]+D->JpCI[m][i][j]+D->JpCI[m][i][j+1]);

        D->PlCR[m][i][j]=D->PlCR[m][i+1][j]
          -0.5/(r+0.5)*m*dtBydr*(D->BzI[m][i+1][j+1]+D->BzI[m][i+1][j])
          -dtBydr*(D->EzR[m][i+1][j+1]-D->EzR[m][i+1][j])
          -pi*dt*(D->JrR[m][i+1][j]+D->JrCR[m][i+1][j]);
        D->PlCI[m][i][j]=D->PlCI[m][i+1][j]
          +0.5/(r+0.5)*m*dtBydr*(D->BzR[m][i+1][j+1]+D->BzR[m][i+1][j])
          -dtBydr*(D->EzI[m][i+1][j+1]-D->EzI[m][i+1][j])
          -pi*dt*(D->JrI[m][i+1][j]+D->JrCI[m][i+1][j]);
        D->SlCR[m][i][j]=D->SlCR[m][i+1][j]
          +0.5/(r+0.5)*m*dtBydr*(D->EzI[m][i+1][j+1]+D->EzI[m][i+1][j])
          -dtBydr*(D->BzR[m][i+1][j+1]-D->BzR[m][i+1][j])
          -0.5*pi*dt*(D->JpR[m][i+1][j]+D->JpR[m][i+1][j+1]+D->JpCR[m][i+1][j]+D->JpCR[m][i+1][j+1]);
        D->SlCI[m][i][j]=D->SlCI[m][i+1][j]
          -0.5/(r+0.5)*m*dtBydr*(D->EzR[m][i+1][j+1]+D->EzR[m][i+1][j])
          -dtBydr*(D->BzI[m][i+1][j+1]-D->BzI[m][i+1][j])
          -0.5*pi*dt*(D->JpI[m][i+1][j]+D->JpI[m][i+1][j+1]+D->JpCI[m][i+1][j]+D->JpCI[m][i+1][j+1]);
      }       //End of i
  }         //End of m
}

//lala
void solve_NDFX(Domain *D)
{
  int i,j,m,numMode,istart,iend,jstart,jend;  
  double dtBydr,dtBydz,r,dr,dz,dt,coefEB1,coefEB2,coefEB3,coefRL,coefF;
  double nowPrR,nowSrR,prevPrR,prevSrR;
  double nowPrI,nowSrI,prevPrI,prevSrI;
  double upPrCR,upPlCR,upSrCR,upSlCR;
  double dnPrCR,dnPlCR,dnSrCR,dnSlCR,PrCR,PlCR,SrCR,SlCR;
  double upPrCI,upPlCI,upSrCI,upSlCI;
  double dnPrCI,dnPlCI,dnSrCI,dnSlCI,PrCI,PlCI,SrCI,SlCI;
  double upFR,dnFR,upFI,dnFI,dF;
  double rr,rd,upL,LdU,tmpr,tmpd,tmp,upr,upd;
  double oldEzR,oldBzR,oldPrR,oldPlR,oldSrR,oldSlR;
  double oldEzI,oldBzI,oldPrI,oldPlI,oldSrI,oldSlI;

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  numMode=D->numMode;

  dtBydr=D->dt/D->dr;
  dtBydz=D->dt/D->dz;
  dr=D->dr; dz=D->dz; dt=D->dt; dF=D->dF;
  upL=(double)(D->ny+D->minYDomain-D->pmlCellUp);
  LdU=D->pmlCellUp;
  rr=D->pmlr;    rd=D->pmld;

  upr=upd=1.0;

  for(m=0; m<numMode; m++) 
    for(j=jstart+1; j<jend; j++) {
      r=(double)(j-jstart);
      if(D->pmlUp==ON) absorb_U(D,&upr,&upd,r,upL,LdU,rr,rd);
      else ;
      nowPrR=D->PrR[m][istart-1][j];
      nowSrR=D->SrR[m][istart-1][j];
      nowPrI=D->PrI[m][istart-1][j];
      nowSrI=D->SrI[m][istart-1][j];
      coefF=dF*m/(r-0.5)*dtBydr;
      for(i=istart; i<iend; i++) {
        upPrCR=0.5*(D->PrCR[m][i-1][j]  +D->PrCR[m][i][j]);
        dnPrCR=0.5*(D->PrCR[m][i-1][j-1]+D->PrCR[m][i][j-1]);
        upPrCI=0.5*(D->PrCI[m][i-1][j]  +D->PrCI[m][i][j]);
        dnPrCI=0.5*(D->PrCI[m][i-1][j-1]+D->PrCI[m][i][j-1]);
        upSrCR=0.5*(D->SrCR[m][i-1][j]  +D->SrCR[m][i][j]);
        dnSrCR=0.5*(D->SrCR[m][i-1][j-1]+D->SrCR[m][i][j-1]);
        upSrCI=0.5*(D->SrCI[m][i-1][j]  +D->SrCI[m][i][j]);
        dnSrCI=0.5*(D->SrCI[m][i-1][j-1]+D->SrCI[m][i][j-1]);
        upPlCR=0.5*(D->PlCR[m][i-1][j]  +D->PlCR[m][i][j]);
        dnPlCR=0.5*(D->PlCR[m][i-1][j-1]+D->PlCR[m][i][j-1]);
        upPlCI=0.5*(D->PlCI[m][i-1][j]  +D->PlCI[m][i][j]);
        dnPlCI=0.5*(D->PlCI[m][i-1][j-1]+D->PlCI[m][i][j-1]);
        upSlCR=0.5*(D->SlCR[m][i-1][j]  +D->SlCR[m][i][j]);
        dnSlCR=0.5*(D->SlCR[m][i-1][j-1]+D->SlCR[m][i][j-1]);
        upSlCI=0.5*(D->SlCI[m][i-1][j]  +D->SlCI[m][i][j]);
        dnSlCI=0.5*(D->SlCI[m][i-1][j-1]+D->SlCI[m][i][j-1]);
        PrCR=0.5*(upPrCR+dnPrCR); PlCR=0.5*(upPlCR+dnPlCR); 
        PrCI=0.5*(upPrCI+dnPrCI); PlCI=0.5*(upPlCI+dnPlCI); 
        SrCR=0.5*(upSrCR+dnSrCR); SlCR=0.5*(upSlCR+dnSlCR); 
        SrCI=0.5*(upSrCI+dnSrCI); SlCI=0.5*(upSlCI+dnSlCI); 
        upFR=0.5*(D->FR[m][i][j+1]+D->FR[m][i-1][j+1]);
        dnFR=0.5*(D->FR[m][i][j]+D->FR[m][i-1][j]);
        upFI=0.5*(D->FI[m][i][j+1]+D->FI[m][i-1][j+1]);
        dnFI=0.5*(D->FI[m][i][j]+D->FI[m][i-1][j]);

        oldEzR=D->EzR[m][i][j]; oldEzI=D->EzI[m][i][j];
        oldBzR=D->BzR[m][i][j]; oldBzI=D->BzI[m][i][j];
        oldPrR=D->PrR[m][i][j]; oldPrI=D->PrI[m][i][j];
        oldPlR=D->PlR[m][i][j]; oldPlI=D->PlI[m][i][j];
        oldSrR=D->SrR[m][i][j]; oldSrI=D->SrI[m][i][j];
        oldSlR=D->SlR[m][i][j]; oldSlI=D->SlI[m][i][j];

        tmpr=upr;     tmpd=upd;

        tmp=tmpr*(0.5/r*dtBydr*(PrCR-PlCR)
              +0.5/r*m*dtBydr*(SlCI-SrCI)
              +0.5*dtBydr*((upPrCR-dnPrCR)-(upPlCR-dnPlCR))
              -pi*dt*(D->JzR[m][i-1][j]+D->JzR[m][i][j])
              +dF*dtBydz*(D->FR[m][i][j]-D->FR[m][i-1][j]));
        D->EzR[m][i][j]=tmpd*(oldEzR+tmp);
        tmp=tmpr*(0.5/r*dtBydr*(PrCI-PlCI)
              -0.5/r*m*dtBydr*(SlCR-SrCR)
              +0.5*dtBydr*((upPrCI-dnPrCI)-(upPlCI-dnPlCI))
              -pi*dt*(D->JzI[m][i-1][j]+D->JzI[m][i][j])
              +dF*dtBydz*(D->FI[m][i][j]-D->FI[m][i-1][j]));
        D->EzI[m][i][j]=tmpd*(oldEzI+tmp);

        tmp=tmpr*(-0.5/r*dtBydr*(SrCR+SlCR)
              -0.5/r*m*dtBydr*(PlCI+PrCI)
              -0.5*dtBydr*((upSrCR-dnSrCR)+(upSlCR-dnSlCR)));
        D->BzR[m][i][j]=tmpd*(oldBzR+tmp);
        tmp=tmpr*(-0.5/r*dtBydr*(SrCI+SlCI)
              +0.5/r*m*dtBydr*(PlCR+PrCR)
              -0.5*dtBydr*((upSrCI-dnSrCI)+(upSlCI-dnSlCI)));
        D->BzI[m][i][j]=tmpd*(oldBzI+tmp);

        prevPrR=nowPrR;
        nowPrR=D->PrR[m][i][j];
        tmp=tmpr*(-0.5/(r+0.5)*m*dtBydr*(D->BzCI[m][i][j]+D->BzCI[m][i][j+1])
              +dtBydr*(D->EzCR[m][i][j+1]-D->EzCR[m][i][j])
              -2.0*pi*dt*D->JrR[m][i][j]
              +dF*dtBydr*(upFR-dnFR));
        D->PrR[m][i][j]=tmpd*(prevPrR+tmp);
        prevPrI=nowPrI;
        nowPrI=D->PrI[m][i][j];
        tmp=tmpr*(0.5/(r+0.5)*m*dtBydr*(D->BzCR[m][i][j]+D->BzCR[m][i][j+1])
              +dtBydr*(D->EzCI[m][i][j+1]-D->EzCI[m][i][j])
              -2.0*pi*dt*D->JrI[m][i][j]
              +dF*dtBydr*(upFI-dnFI));
        D->PrI[m][i][j]=tmpd*(prevPrI+tmp);

        prevSrR=nowSrR;
        nowSrR=D->SrR[m][i][j];
        tmp=tmpr*(-0.5/(r+0.5)*m*dtBydr*(D->EzCI[m][i][j]+D->EzCI[m][i][j+1])
              -dtBydr*(D->BzCR[m][i][j+1]-D->BzCR[m][i][j])
              -pi*dt*(D->JpR[m][i][j]+D->JpR[m][i][j+1])
              -dF*m/(r+0.5)*dtBydr*0.5*(upFI+dnFI));
        D->SrR[m][i][j]=tmpd*(prevSrR+tmp);
        prevSrI=nowSrI;
        nowSrI=D->SrI[m][i][j];
        tmp=tmpr*(0.5/(r+0.5)*m*dtBydr*(D->EzCR[m][i][j]+D->EzCR[m][i][j+1])
              -dtBydr*(D->BzCI[m][i][j+1]-D->BzCI[m][i][j])
              -pi*dt*(D->JpI[m][i][j]+D->JpI[m][i][j+1])
              +dF*m/(r+0.5)*dtBydr*0.5*(upFR+dnFR));
        D->SrI[m][i][j]=tmpd*(prevSrI+tmp);

        tmp=tmpr*(-0.5/(r+0.5)*m*dtBydr*(D->BzCI[m][i+1][j+1]+D->BzCI[m][i+1][j])
              -dtBydr*(D->EzCR[m][i+1][j+1]-D->EzCR[m][i+1][j])
              -2.0*pi*dt*D->JrR[m][i+1][j]
              +dF*dtBydr*(upFR-dnFR));
        D->PlR[m][i][j]=tmpd*(D->PlR[m][i+1][j]+tmp);
        tmp=tmpr*(0.5/(r+0.5)*m*dtBydr*(D->BzCR[m][i+1][j+1]+D->BzCR[m][i+1][j])
              -dtBydr*(D->EzCI[m][i+1][j+1]-D->EzCI[m][i+1][j])
              -2.0*pi*dt*D->JrI[m][i+1][j]
              +dF*dtBydr*(upFI-dnFI));
        D->PlI[m][i][j]=tmpd*(D->PlI[m][i+1][j]+tmp);

        tmp=tmpr*(0.5/(r+0.5)*m*dtBydr*(D->EzCI[m][i+1][j+1]+D->EzCI[m][i+1][j])
              -dtBydr*(D->BzCR[m][i+1][j+1]-D->BzCR[m][i+1][j])
              -pi*dt*(D->JpR[m][i+1][j]+D->JpR[m][i+1][j+1])
              -dF*m/(r+0.5)*dtBydr*0.5*(upFI+dnFI));
        D->SlR[m][i][j]=tmpd*(D->SlR[m][i+1][j]+tmp);
        tmp=tmpr*(-0.5/(r+0.5)*m*dtBydr*(D->EzCR[m][i+1][j+1]+D->EzCR[m][i+1][j])
              -dtBydr*(D->BzCI[m][i+1][j+1]-D->BzCI[m][i+1][j])
              -pi*dt*(D->JpI[m][i+1][j]+D->JpI[m][i+1][j+1])
              +dF*m/(r+0.5)*dtBydr*0.5*(upFR+dnFR));
        D->SlI[m][i][j]=tmpd*(D->SlI[m][i+1][j]+tmp);
      }       //End of i
    }         //End of j

  j=jstart; m=0;
  for(i=istart; i<iend; i++) {
    D->BzR[m][i][j]+=
      -dtBydr*(D->SrCR[m][i-1][j]+D->SrCR[m][i][j]+D->SlCR[m][i-1][j]+D->SlCR[m][i][j]);
    D->EzR[m][i][j]+=
       dtBydr*(D->PrCR[m][i-1][j]+D->PrCR[m][i][j]-D->PlCR[m][i-1][j]-D->PlCR[m][i][j])
      -pi*dt*(D->JzR[m][i-1][j]+D->JzR[m][i][j])
      +dF*dtBydz*(D->FR[m][i][j]-D->FR[m][i-1][j]);
  }

  for(m=0; m<numMode; m++) { 
      r=(double)(j-jstart);
      nowPrR=D->PrR[m][istart-1][j];
      nowSrR=D->SrR[m][istart-1][j];
      nowPrI=D->PrI[m][istart-1][j];
      nowSrI=D->SrI[m][istart-1][j];
      coefF=dF*m/(r-0.5)*dtBydr;
      for(i=istart; i<iend; i++) {
        upFR=0.5*(D->FR[m][i][j]+D->FR[m][i-1][j]);
        dnFR=0.5*(D->FR[m][i][j-1]+D->FR[m][i-1][j-1]);
        upFI=0.5*(D->FI[m][i][j]+D->FI[m][i-1][j]);
        dnFI=0.5*(D->FI[m][i][j-1]+D->FI[m][i-1][j-1]);

        prevPrR=nowPrR;
        nowPrR=D->PrR[m][i][j];
        D->PrR[m][i][j]=prevPrR
          -0.5/(r+0.5)*m*dtBydr*(D->BzCI[m][i][j]+D->BzCI[m][i][j+1])
          +dtBydr*(D->EzCR[m][i][j+1]-D->EzCR[m][i][j])
          -2.0*pi*dt*D->JrR[m][i][j];
          +dF*dtBydr*(upFR-dnFR);
        prevPrI=nowPrI;
        nowPrI=D->PrI[m][i][j];
        D->PrI[m][i][j]=prevPrI
          +0.5/(r+0.5)*m*dtBydr*(D->BzCR[m][i][j]+D->BzCR[m][i][j+1])
          +dtBydr*(D->EzCI[m][i][j+1]-D->EzCI[m][i][j])
          -2.0*pi*dt*D->JrI[m][i][j]
          +dF*dtBydr*(upFI-dnFI);
        prevSrR=nowSrR;
        nowSrR=D->SrR[m][i][j];
        D->SrR[m][i][j]=prevSrR
          -0.5/(r+0.5)*m*dtBydr*(D->EzCI[m][i][j]+D->EzCI[m][i][j+1])
          -dtBydr*(D->BzCR[m][i][j+1]-D->BzCR[m][i][j])
          -pi*dt*(D->JpR[m][i][j]+D->JpR[m][i][j+1])
          -dF*m/(r+0.5)*dtBydr*0.5*(upFI+dnFI);

        prevSrI=nowSrI;
        nowSrI=D->SrI[m][i][j];
        D->SrI[m][i][j]=prevSrI
          +0.5/(r+0.5)*m*dtBydr*(D->EzCR[m][i][j]+D->EzCR[m][i][j+1])
          -dtBydr*(D->BzCI[m][i][j+1]-D->BzCI[m][i][j])
          -pi*dt*(D->JpI[m][i][j]+D->JpI[m][i][j+1]);
          +dF*m/(r+0.5)*dtBydr*0.5*(upFR+dnFR);

        D->PlR[m][i][j]=D->PlR[m][i+1][j]
          -0.5/(r+0.5)*m*dtBydr*(D->BzCI[m][i+1][j+1]+D->BzCI[m][i+1][j])
          -dtBydr*(D->EzCR[m][i+1][j+1]-D->EzCR[m][i+1][j])
          -2.0*pi*dt*D->JrR[m][i+1][j]
          +dF*dtBydr*(upFR-dnFR);
        D->PlI[m][i][j]=D->PlI[m][i+1][j]
          +0.5/(r+0.5)*m*dtBydr*(D->BzCR[m][i+1][j+1]+D->BzCR[m][i+1][j])
          -dtBydr*(D->EzCI[m][i+1][j+1]-D->EzCI[m][i+1][j])
          -2.0*pi*dt*D->JrI[m][i+1][j]
          +dF*dtBydr*(upFI-dnFI);
        D->SlR[m][i][j]=D->SlR[m][i+1][j]
          +0.5/(r+0.5)*m*dtBydr*(D->EzCI[m][i+1][j+1]+D->EzCI[m][i+1][j])
          -dtBydr*(D->BzCR[m][i+1][j+1]-D->BzCR[m][i+1][j])
          -pi*dt*(D->JpR[m][i+1][j]+D->JpR[m][i+1][j+1])
          -dF*m/(r+0.5)*dtBydr*0.5*(upFI+dnFI);
        D->SlI[m][i][j]=D->SlI[m][i+1][j]
          -0.5/(r+0.5)*m*dtBydr*(D->EzCR[m][i+1][j+1]+D->EzCR[m][i+1][j])
          -dtBydr*(D->BzCI[m][i+1][j+1]-D->BzCI[m][i+1][j])
          -pi*dt*(D->JpI[m][i+1][j]+D->JpI[m][i+1][j+1])
          +dF*m/(r+0.5)*dtBydr*0.5*(upFR+dnFR);
      }       //End of i
    }         //End of j

}




void filter(Domain *D)
{
  int a,i,j,m,istart,iend,jstart,jend,numMode;
  double alpha[2];
  double nowEpR,nowEzR,nowErR,nowBpR,nowBzR,nowBrR;
  double nowEpI,nowEzI,nowErI,nowBpI,nowBzI,nowBrI;
  double beforeEpR,beforeEzR,beforeErR,beforeBpR,beforeBzR,beforeBrR;
  double beforeEpI,beforeEzI,beforeErI,beforeBpI,beforeBzI,beforeBrI;

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  numMode=D->numMode;

  alpha[0]=0.5;
  alpha[1]=1.5;
  j=jstart;
  for(a=0; a<2; a++) 
    for(m=0; m<numMode; m++)
    {
      nowEpR=D->EpR[m][0][j];  nowEpI=D->EpI[m][0][j];
      nowEzR=D->EzR[m][0][j];  nowEzI=D->EzI[m][0][j];
      nowErR=D->ErR[m][0][j];  nowErI=D->ErI[m][0][j];
      nowBpR=D->BpR[m][0][j];  nowBpI=D->BpI[m][0][j];
      nowBzR=D->BzR[m][0][j];  nowBzI=D->BzI[m][0][j];
      nowBrR=D->BrR[m][0][j];  nowBrI=D->BrI[m][0][j];
      for(i=istart-1; i<iend+2; i++) {
        beforeEpR=nowEpR;      beforeEpI=nowEpI;
        beforeEzR=nowEzR;      beforeEzI=nowEzI;
        beforeErR=nowErR;      beforeErI=nowErI;
        beforeBpR=nowBpR;      beforeBpI=nowBpI;
        beforeBzR=nowBzR;      beforeBzI=nowBzI;
        beforeBrR=nowBrR;      beforeBrI=nowBrI;
        nowEpR=D->EpR[m][i][j];    nowEpI=D->EpI[m][i][j];
        nowEzR=D->EzR[m][i][j];    nowEzI=D->EzI[m][i][j];
        nowErR=D->ErR[m][i][j];    nowErI=D->ErI[m][i][j];

        D->EpR[m][i][j+1]=alpha[a]*nowEpR+(1.0-alpha[a])*(beforeEpR+D->EpR[m][i+1][j+1])*0.5;
        D->EpI[m][i][j+1]=alpha[a]*nowEpI+(1.0-alpha[a])*(beforeEpI+D->EpI[m][i+1][j+1])*0.5;
        D->EzR[m][i][j+1]=alpha[a]*nowEzR+(1.0-alpha[a])*(beforeEzR+D->EzR[m][i+1][j+1])*0.5;
        D->EzI[m][i][j+1]=alpha[a]*nowEzI+(1.0-alpha[a])*(beforeEzI+D->EzI[m][i+1][j+1])*0.5;
        D->ErR[m][i][j]=alpha[a]*nowErR+(1.0-alpha[a])*(beforeErR+D->ErR[m][i+1][j])*0.5;
        D->ErI[m][i][j]=alpha[a]*nowErI+(1.0-alpha[a])*(beforeErI+D->ErI[m][i+1][j])*0.5;
        D->BpR[m][i][j]=alpha[a]*nowBpR+(1.0-alpha[a])*(beforeBpR+D->BpR[m][i+1][j])*0.5;
        D->BpI[m][i][j]=alpha[a]*nowBpI+(1.0-alpha[a])*(beforeBpI+D->BpI[m][i+1][j])*0.5;
        D->BzR[m][i][j]=alpha[a]*nowBzR+(1.0-alpha[a])*(beforeBzR+D->BzR[m][i+1][j])*0.5;
        D->BzI[m][i][j]=alpha[a]*nowBzI+(1.0-alpha[a])*(beforeBzI+D->BzI[m][i+1][j])*0.5;
        D->BrR[m][i][j+1]=alpha[a]*nowBrR+(1.0-alpha[a])*(beforeBrR+D->BrR[m][i+1][j+1])*0.5;
        D->BrI[m][i][j+1]=alpha[a]*nowBrI+(1.0-alpha[a])*(beforeBrI+D->BrI[m][i+1][j+1])*0.5;
      }
    }
}
