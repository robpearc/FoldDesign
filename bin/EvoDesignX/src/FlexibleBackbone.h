///////////////////////////////////////////////////////////////////////////////////////
//Copyright (C) Yang Zhang Lab, University of Michigan - All Rights Reserved
//Unauthorized copying of this file, via any medium is strictly prohibited
//Written by Robin Pearce <robpearc@umich.edu>
///////////////////////////////////////////////////////////////////////////////////////

#ifndef FLEXIBLEBACKBONE_H
#define FLEXIBLEBACKBONE_H

#include "EnergyFunction.h"
#include "Structure.h"

#define MODULUS    2147483647 /* DON'T CHANGE THIS VALUE                  */
#define MULTIPLIER 48271      /* DON'T CHANGE THIS VALUE                  */
#define CHECK      399268537  /* DON'T CHANGE THIS VALUE                  */
#define STREAMS    256        /* # of streams, DON'T CHANGE THIS VALUE    */
#define A256       22925      /* jump multiplier, DON'T CHANGE THIS VALUE */
#define DEFAULT    123456789  /* initial seed, use 0 < DEFAULT < MODULUS  */
#define PI         3.14159
#define err        1e-7

typedef struct _sssegment{
  int init;
  int term;
  char ss;
} sssegment;

typedef struct _segment{
  int init;
  int term;
} segment;

typedef struct _paircont{
  int init;
  int term;
  double dist;
  double prob;
} paircont;

typedef struct _pairaa{
  double dist;
  double dstd;
  double weight;
} pairaa;


int LoadRamaOutlier(double*** rlogduke,int*** ramaduke,char* ramafile);
int LoadTorDistribution(double** tordis,double* accdis,char* torDisFile);
int LoadPhiPsiProb(double ***phipsiprob,char *torFileName);
int Loadca2ncbins(float ****cancbins,char *filename);
int Loadbbord(double *bborder,char *filename);
int LoadSecondaryStructureDesign(ChainFlex *pChainFlex,sssegment *sse,int &numsse,char *ssFile);
void calcbbptable(ChainFlex *pChainFlex,sssegment *sse,int &numbbptb,int numsse,double *bborder,paircont *bbptb,int *alphasind,int &numsalpha);
int loadContact(double **contactCA,double **contactCB,int numseq,char *filename);
int loadDistance(double **distanceCB,double **distanceWeight,int numseq,char *filename);
void calcContactEnergy(Chain *pChain,double **contactCA,double **contactCB,double &caEnergy,double &cbEnergy);
double calcDistanceEnergy(Chain *pChain,double **distanceCB,double **distanceWeights);
double calcConstraintEnergyWeighted(Chain *pChain,pairaa *natDistRestr,double globalWeight);
double Random(void);
XYZ rana(double theta);
XYZ ranv_ergodic(double fac);
void swap(int *a, int *b);
void randomizeArray(int arr[], int numDesignSites);
int OptimizeSideChains(Structure* pStructure, int numCycles, int deschn);
int getResidueIndex(Residue *newResi);
void rotate_trans_chain(Chain* pChainNew, Chain* pChainOrig,int chain1_len,int chain2_len, int start, int end, int type);

float calcTorsion(float xi,float yi,float zi,
                  float xj,float yj,float zj,
                  float xk,float yk,float zk,
                  float xl,float yl,float zl);

int str2torp(Chain *pChain,ChainFlex *pChainFlex,int istart,int iend);
int str2tor(Chain *pChain,ChainFlex *pChainFlex);
BOOL tor2strp2p2(Chain *pChain,ChainFlex *pChainFlex,int istart,int iend);


BOOL NewMainChainPos(Chain* pChain, ChainFlex* pChainFlex, int istart);
BOOL convertTorToPos(double xi,double yi,double zi,
                     double xj,double yj,double zj,
                     double xk,double yk,double zk,
                     double tang,double tleng,double tinner,
                     double *xl,double *yl,double *zl);

void extractCaDistances(Chain *pChain,pairaa *natDistRestr);
double calcConstraintEnergy(Chain *pChain,pairaa *natDistRestr);


int StructureCalcEnergyFlex(Structure* pStructure,BBdepRotamerLib* pBBdepRotLib,double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]);
int getAngleOutliers(Chain *pChain,ChainFlex *pChainFlex);
int getBondLengthOutliers(Chain *pChain,ChainFlex *pChainFlex);
int getRamaOutliers(Chain *pChain,ChainFlex *pChainFlex,int ***ramaduke);
int getClashes(Chain *origChain,ChainFlex *origChainFlex);

BOOL makeConformationalMove(Structure *decstr,StructureFlex *pStructureFlex,int chains[],int numchains,double &oldenergy,double &newenergy,
                            BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,
                            double *accdis,double ***phipsiprob,float ****cancbins,sssegment *sse,int numsse,paircont *bbptb,
                            int numbbptb,int *alphasind,int numsalpha,pairaa *natDistRestr,double **contactCA,double **contactCB);


BOOL moveAngle(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
               bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,pairaa *natDistRestr,
               double **contactCA,double **contactCB);


BOOL moveLength(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
               bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,pairaa *natDistRestr,
               double **contactCA,double **contactCB);

BOOL moveOmega(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
               bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,pairaa *natDistRestr,
               double **contactCA,double **contactCB);
  
BOOL moveTor(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
             bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,
             double *accdis,pairaa *natDistRestr,double **contactCA,double **contactCB);

BOOL movePhi(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
             bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,
             double *accdis,double ***phipsiprob,pairaa *natDistRestr,double **contactCA,double **contactCB);

BOOL movePsi(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
             bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,
             double *accdis,double ***phipsiprob,pairaa *natDistRestr,double **contactCA,double **contactCB);

BOOL moveAnc(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
             bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,
             double *accdis,double ***phipsiprob,float ****cancbins,pairaa *natDistRestr,double **contactCA,double **contactCB);

BOOL moveRotation(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
                  bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,pairaa *natDistRestr,
                  double **contactCA,double **contactCB);

BOOL moveBackrub(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
                 bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,pairaa *natDistRestr,
                 double **contactCA,double **contactCB);

BOOL moveLMP(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
             bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,pairaa *natDistRestr,
             double **contactCA,double **contactCB);

void singlemoveLMP(Chain *pChain,int k,int mtype,int flagdir);

BOOL mcfragsweepLMP2(Chain *pChain,int poss,int pose);


BOOL moveShift(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
             bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,pairaa *natDistRestr,
             double **contactCA,double **contactCB);

BOOL movePos(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
             bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,pairaa *natDistRestr,
             double **contactCA,double **contactCB);

BOOL moveBetaPack(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
                  bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,
                  paircont *bbptb,int numbbptb,pairaa *natDistRestr,double **contactCA,double **contactCB);

BOOL moveHelixPack(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
                   bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,
                   sssegment *sse,int *alphasind,int numsalpha,pairaa *natDistRestr,double **contactCA,double **contactCB);

BOOL moveHelix(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
               bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,pairaa *natDistRestr,
               double **contactCA,double **contactCB);

BOOL genesse(ChainFlex *pChainFlex,sssegment *sse,int &numsse);

void getfrombbptable(double tval,segment *sgval,paircont* tbbptb,int tnumbbptb);

#endif
