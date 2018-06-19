#define Electron 	1
#define HPlus0 	 	100
#define HPlus1 	 	101
#define HePlus0 	200
#define HePlus1 	201
#define HePlus2 	202
#define CPlus0          600
#define CPlus1          601
#define CPlus2          602
#define CPlus3          603
#define CPlus4          604
#define CPlus5          605
#define CPlus6          606
#define NPlus0          700
#define NPlus1          701
#define NPlus2          702
#define NPlus3          703
#define NPlus4          704
#define NPlus5          705
#define NPlus6          706
#define NPlus7          707
#define NePlus0         1000
#define NePlus1         1001
#define NePlus2         1002
#define NePlus3         1003
#define NePlus4         1004
#define NePlus5         1005
#define NePlus6         1006
#define NePlus7         1007
#define NePlus8         1008
#define NePlus9         1009
#define NePlus10        1010
#define AlPlus4         1304
#define userDefined   	9999999

#define Polygon    	1
#define Defined    	2
#define Channel    	3
#define BoostFrame	4
#define Circle    	5
#define Exp    		6

#define Constant   	0
#define Gaussian   	1
#define Polynomial   	2

#define byNumber	0
#define byDensity	1

typedef struct _LoadList  {
   int type;
   int species;
   double superP;
   double density;
   double numberRZ;
   double numberPhi;
   double criticalDensity;
   double targetW;
   int index;  
   double num;      //exceeded number of particle which is less than 1
   int xnodes;     //longitudinal point number
   double *xn;      //longitudinal density (1 is P->density)
   double *xpoint;    //longitudinal point
   int ynodes;     //transverse point number
   double *yn;      //transverse density (1 is P->density)
   double *ypoint;    //transverse point
   int znodes;     //transverse point number
   double *zn;      //transverse density (1 is P->density)
   double *zpoint;    //transverse point
   double givenMinPx;	//for saveParticle option

   //initial momentum distribution
   double z0;
   double pz0;
   double delPz;
   double delZ;

   //defined plasma
   int defineMode;
   double *xPosition;
   double *yPosition;
   double *zPosition;
   int numDefined;
   double xLengDef;
   double yLengDef;
   double zLengDef;
   int numDefPtcls;
   double **define;
   int maxLoadTime;
   int minLoadTime;
   double minX;
   double maxX;
   double minY;
   double maxY;
   double minZ;
   double maxZ;

   int pointPosition;   
   double p1;
   double p2;
   double p3;

   double mass;
   int charge;
   
   double temperatureR, temperatureZ;
   
   //applying function
   double centerX;
   double centerY;
   double centerZ;
   double gaussCoefX;
   double polyCoefX;
   double gaussCoefYZ;
   double polyCoefYZ;
   int modeX;
   int modeYZ;

   //ionization
   int levels;
   int ionFinal;
   double givenMinA0;
   double *ionEnergy;
   double *ionW;

   //pair current
   int pair;

   struct _LoadList *next;
} LoadList;
