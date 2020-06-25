#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include "mpi.h"
#include <time.h>

int main(int argc, char *argv[])
{
    int i,j,k,n,s,iteration=0,boost,filterStep,labSaveStep;
    int rnk,suddenDump=OFF,shiftIteration,iter=1;
    double factor,time_spent,t,dF=0.0,x,***val;
    clock_t begin,end;
    struct tm *t_now;
    time_t timer; 	//measure time
    char name[100];
    FILE *out;
    Domain D;  
    LaserList *L;
    LoadList *LL;
    External Ext;
    int myrank, nTasks;
    MPI_Status status; 

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


    begin=clock();

    if(argc < 2) 
    {  
      printf("mpirun -np N show [inputFile] [dumpNum]\n"); 
      exit(0); 
    }

    timer=time(NULL);
    t_now=localtime(&timer);
    if(myrank==0) {  
      sprintf(name,"report");
      out = fopen(name,"a");
      fprintf(out,"simulation start.\n");
      fprintf(out,"%d-%d-%d %d:%d:%d\n",t_now->tm_year+1900,t_now->tm_mon+1,t_now->tm_mday,t_now->tm_hour,t_now->tm_min,t_now->tm_sec);
    } else ;

    //parameter setting
    parameterSetting(&D,&Ext,argv[1]);
    if(argc >= 3) { 
      D.dumpStep = atoi(argv[2]); 
      if(D.dumpStart==D.dumpStep) D.dumpStart+=1; else;
    } else;

    //create mesh
    boundary(&D,&Ext);
    MPI_Barrier(MPI_COMM_WORLD);

    //load plasma or load dump file
    if(argc >= 3)  {   
      iteration=D.dumpStep;
      restoreDump(D,iteration);
      t=D.dt*iteration;
      sprintf(name,"dumpField%d.h5",iteration);
      if(myrank==0) restoreIntMeta(name,"/minXDomain",&(D.minXDomain),1);
      else ;
      MPI_Bcast(&(D.minXDomain),1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
      D.maxXDomain=D.minXDomain+D.nx;
      D.minXSub+=D.minXDomain;
    }  else   {
      LL=D.loadList;
      s=0;
      while(LL->next)      {
        loadPlasma(&D,LL,s,iteration,D.istart,D.iend,D.jstart,D.jend);
        LL=LL->next;
        s++;
      }
      t=0;

    }
    MPI_Barrier(MPI_COMM_WORLD);

    //pair charge
    if(iteration==0) {
      LL=D.loadList; s=0;
      while(LL->next)      {
        if(LL->pair==ON) {
          solveCharge(&D,LL,D.RhoPairR,D.RhoPairI,D.istart,D.iend,D.jstart,D.jend,s,-1.0);
        }  else ;
        LL=LL->next; s++;
      }
      if(D.L>1)  {
        MPI_TransferDen_Xplus(&D,D.RhoPairR,D.RhoPairI,D.nySub+5,3);
        MPI_TransferDen_Xminus(&D,D.RhoPairR,D.RhoPairI,D.nySub+5,3);
      }  else ;

      val=(double ***)malloc(1*sizeof(double ** ));
      for(n=0; n<1; n++) {
        val[n]=(double **)malloc((D.nxSub+5)*sizeof(double * ));
        for(i=0; i<D.nxSub+5; i++)
          val[n][i]=(double *)malloc((D.nySub+5)*sizeof(double  ));
      }
      for(i=0; i<D.nxSub+5; i++)
        for(j=0; j<D.nySub+5; j++)
          val[0][i][j]=0.0;

      filter_current(&D,val,D.RhoPairR,iter);
      filter_current(&D,val,D.RhoPairI,iter);

      for(i=0; i<D.nxSub+5; i++) free(val[0][i]);
      free(val[0]); free(val);

    } else ;


    //rooping time 
    while(iteration<=D.maxStep)
    {
       //save dump File
       if(D.dumpSave==ON && iteration>=D.dumpStart && iteration%D.dumpSaveStep==0) {
         saveDump(D,iteration); 
         end=clock();
         time_spent=(end-begin)/CLOCKS_PER_SEC/60.0;
         if(myrank==0) {
           fprintf(out,"Time duration at %ddump:%.4gmin\n",iteration,time_spent);
           printf("Time duration at %ddump:%.4gmin\n",iteration,time_spent);
         } else ;
       } else	;


       //save File      
       if(iteration%D.saveStep==0 && iteration>=D.saveStart)
         saveFile(D,iteration);
       else	;
       //save center field      
       if(iteration%D.centerStep==0)  {
         saveCenterField(&D,iteration);
//         saveCenterDensity(&D,iteration);
       }  else	;
       MPI_Barrier(MPI_COMM_WORLD);

       //redistributing particles lala
//       if(D.redist==ON && iteration>0) particle_redist(&D,iteration,&Ext); else;
       fieldSolve1(D,t,iteration,dF);
//printf("fieldSolve1,iteration=%d\n",iteration);

       //filter
       interpolation(&D,&Ext,iteration);
//printf("interpolation,iteration=%d\n",iteration);

       particlePush(&D,iteration);
//printf("particle Push,iteration=%d\n",iteration);

       if(D.fieldIonization==ON) fieldIonization(&D); else;

       updateCurrent(D,iteration);
//printf("current,iteration=%d\n",iteration);
//       calConservation(D,iteration);

       if(D.moveIt==ON) {
         if(iteration%D.shiftDuration!=0)    {
           dF=D.dF;
           movingDomain(&D,iteration);
           if(myrank==D.L-1) {
             LL=D.loadList; s=0;
             while(LL->next)      {
               loadPlasma(&D,LL,s,iteration,D.iend-1,D.iend,D.jstart,D.jend);
               LL=LL->next; s++;
             }
           } else ;
           rearrangeParticles(&D);
           if(D.L>1)   particleShareX(D);   else	;
//           if(D.M>1)   particleShareY(D);   else	;
           if(myrank==D.L-1) {
             LL=D.loadList; s=0;
             while(LL->next)      {
               if(LL->pair==ON)
                 solveCharge(&D,LL,D.RhoPairR,D.RhoPairI,D.iend-1,D.iend,D.jstart,D.jend,s,-1.0);
               else ;
               LL=LL->next;  s++;
             }
           } else ;
           removeEdge(&D);
         } else  {
           rearrangeParticles(&D);
           if(D.L>1)  particleShareX(D);   else	;
    //       if(D.M>1)  particleShareY(D);   else	;
         } 
       } 
       else {
//         if(iteration>=D.nx && D.moving==ON && D.boostOn==OFF 
//            && (iteration-D.nx)%D.shiftDuration!=0)    {
         x=D.movingV*iteration;
         if(D.moving==ON && x>D.maxXDomain-1)    {

           dF=D.dF;
           movingDomain(&D,iteration);

           if(myrank==D.L-1) {
             LL=D.loadList; s=0;
             while(LL->next)      {
               loadPlasma(&D,LL,s,iteration,D.iend-1,D.iend,D.jstart,D.jend);
               LL=LL->next; s++;
             }
           } else ; 
           rearrangeParticles(&D);
           if(D.L>1)   particleShareX(D);   else	;

           if(myrank==D.L-1) {
             LL=D.loadList; s=0;
             while(LL->next)      {
               if(LL->pair==ON)
                 solveCharge(&D,LL,D.RhoPairR,D.RhoPairI,D.iend-1,D.iend,D.jstart,D.jend,s,-1.0);
               else ;
               LL=LL->next;  s++;
             }
           } else ;
           removeEdge(&D);

         } else       {
           rearrangeParticles(&D);
           if(D.L>1)  particleShareX(D);   else	;
           removeEdge(&D);
         }  
       }

//       if(D.filter==ON) filterField(&D); else ;
       solveF(D);
       fieldSolve2(D,t,iteration,dF);
//printf("fieldSolve2,iteration=%d\n",iteration);

       //time update
       if(iteration%10==0 && myrank==0)  
          printf("iteration = %d\n",iteration);           
       iteration+=1;
       t=D.dt*iteration;  

    }     //end of time roop                  
//    if(D.tracking==ON)
//      saveTracking(&D);

    end=clock();
    time_spent=(end-begin)/CLOCKS_PER_SEC;

    //make 'report' file
    if(myrank==0) {
      fprintf(out,"nx=%d, ",D.nx);
      fprintf(out,"ny=%d, ",D.ny);
      fprintf(out,"cores=%d, \n",nTasks);
      fprintf(out,"nSpecies=%d\n",D.nSpecies);
      fprintf(out,"running time=%.4gm\n",time_spent/60.0);
      fprintf(out,"\n");
      fclose(out);
    }
    else	;

    cleanMemory(&D);
  
    MPI_Finalize();

    return 0;
}
