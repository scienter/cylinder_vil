#define DEFAULT		0
#define ADDITION	1

typedef struct _LaserList  {
   int polarity;
   int mode;
   int add;
   double lambda;
   double omega;
   double amplitude;   //unit is a0.
   double rU;
   double rD;
   double retard;
   double flat;
   int loadPointX; 
   int loadPointY; 
   int loadPointZ; 

   double rayleighLength;
   double beamWaist;
   double focus;
   int direction;
   double gdd;
  
   struct _LaserList *next;
} LaserList;
