///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef FOLDDESIGNMOVEMENT_H
#define FOLDDESIGNMOVEMENT_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include "ParsePDB.h"
#include "BasicFunc.h"
#include "FoldDesignEnergyFunction.h"
#include "GeometryCalc.h"
#include "CommonPara.h"
#include "InputData.h"
#include "ParseFoldDesignMovementFiles.h"
#define MAX_ITERATION 1000

class FoldDesignMovement 
{
  public:
    FoldDesignMovement();
    virtual ~FoldDesignMovement();
    FoldDesignMovement & operator=(const FoldDesignMovement &inputMovement);
    void configureSimulation(char *libDir,char *dataDir,
                             const ParseFoldDesignMovementFiles &movementParameters,
                             FoldDesignEnergyFunction inputEnergyFunction,
                             InputData input,GeometryCalc inputGeometry,
                             double tempLow,double tempHigh,double scale1,double scale2,
                             int numRep);

    //---Movement acceptance statistics--->
    int summcsub[3];
    int summctop[3];
    int summcphi[3];
    int summcpsi[3];
    int summcome[3];
    int summclen[3];
    int summcang[3];
    int summcccd[3];
    int summcLMP[3];
    int summcmid[3];
    int summcbbp[3];
    int summcaaa[3];
    int summcbtn[3];
    int summcsft[3];
    int summctra[3];
    int summctot[20];

    void setHbondRamp(int numcycles);
    int getMoveType();
    void setpara8(int seqLength);
    void setInitDecoy(point3f *decstr,int numseq);
    void setInitDecoyFixed(point3f *decstr,int numseq);
    bool attemptConformationalMove(point3f *decstr,int numseq,double tbeta,
                                   double oldenergy,double *newenergy,
                                   int *movtype);//full composite 
    bool mcfragsweepcom_noene(point3f *decstr,int numseq, const int ntype);
    int getmovetype(double tmov[],int totmov,double trandnum);
    void settmstr(int seqLength);
    void calcabind();
    void calcbbptable(int numseq);
    double enelist[20];//used for output
    double enelistbk[20];
    double *temparray,*ttemparray,*tttemparray;
    int indCyc;
    vector<int> fixed_pos;
    bool fixed_flag;
  
  private:
    //---Object instantiation--->
    BasicFunc bf;
    ParsePDB pp;
    FoldDesignEnergyFunction energyFunction;
    GeometryCalc geometry; 
    InputData inputInfo;

    //------------------Simulation parameters-------------------->
    int hbondWeightRamp;       //ramp up hbond weight as simulation progresses 
    int n_rep;                 //number of REMC replicas
    double T1,T2,S1,S2;        //temperatures for REMC
    int i_move_type;
    int numbeta,numalpha,numsalpha;
    int indcut;
    int ntopfrag;
    int numbbptb;
    int numeturn;
    int *betaind,*alphaind,*alphasind;
    int *listeturn;
    int flagres;//0no 1 partrmsd 2 distrest 3 both
    int flagMove;//0 normal 1 only local
    vector<int> start_fixed,end_fixed;
    bool flagcaca;//false intotal true notintot
    paircont *bbptb;
    double *bborder;
    double *probeturn;
    double **phiprob;
    double **psiprob;
    point3f *tmstr;
    point3f *tmstr2;
    point3f *fixed_decstr;
    point3f **fragcont[20];
    point3f *topdhlist[nosegdh];
    int *numtopdh;	
    bool flagLocalMove;
    bool flagVerbose;

    //---------------Functions for conformational movement--------------------->
    bool mcfragsweepfragsub(point3f *decstr,int numseq,double tbeta,
                            double oldenergy,double *newenergy);//dihedral and substitution	
    bool mcfragsweepmid(point3f *decstr,int numseq,double tbeta,
                        double oldenergy,double *newenergy);//rot middle
    bool mcfragsweeptra(point3f *decstr,int numseq,double tbeta,
                        double oldenergy,double *newenergy);//translate middle
    bool mcfragsweeptopdh(point3f *decstr,int numseq,double tbeta,double oldenergy,double *newenergy);
    bool mcfragsweepphi(point3f *decstr,int numseq,double tbeta,double oldenergy,double *newenergy);//phi 
    bool mcfragsweeppsi(point3f *decstr,int numseq,double tbeta,double oldenergy,double *newenergy);//psi 
    bool mcfragsweepome(point3f *decstr,int numseq,double tbeta,double oldenergy,double *newenergy);//omega 
    bool mcfragsweeplen(point3f *decstr,int numseq,double tbeta,double oldenergy,double *newenergy);//length
    bool mcfragsweepang(point3f *decstr,int numseq,double tbeta,double oldenergy,double *newenergy);//angle
    bool mcfragsweepLMP2(point3f *decstr,int numseq,int poss,int pose);//fixed positions
    bool mcfragsweepLMP(point3f *decstr,int numseq,double tbeta,double oldenergy,double *newenergy);//LMP
    void singlemoveLMPf(point3f *tmstr,int k,int mtype);
    void singlemoveLMPb(point3f *tmstr,int k,int mtype);
    bool mcfragsweepsft(point3f *decstr,int numseq,double tbeta,double oldenergy,double *newenergy);
    bool mcfragsweepCCD3(point3f *decstr,int numseq,int lps,int lpe,int lpt,
                         point3d pt1,point3d pt2,point3d pt3,double oldenergy,
                         double *newenergy);//bsheet 3 points
    bool mcfragsweepCCD4(point3f *decstr,int numseq,int lps,int lpe,point3d pt1,
                         int ind1,point3d pt2,int ind2,double oldenergy,
                         double *newenergy);//aaa 2 points
    bool mcfragsweepCCD5(point3f *decstr,int numseq,int lps,int lpe,point3d pt1,
                         int ind1,point3d pt2,int ind2,double oldenergy,
                         double *newenergy);//bturn 2points
    bool mcfragsweepCCD6(point3f *decstr,int numseq,int lps,int lpe,point3d pt1,
                         int ind1,point3d pt2,int ind2,double oldenergy,
                         double *newenergy);//sub 2 end points
    bool mcfragsweepbbp(point3f *decstr,int numseq,double tbeta,double oldenergy,double *newenergy);//strand strand pair
    bool mcfragsweepaaa(point3f *decstr,int numseq,double tbeta,double oldenergy,double *newenergy);//helix helix angle
    bool mcfragsweepbtn(point3f *decstr,int numseq,double tbeta,double oldenergy,double *newenergy);//beta turn
    void getfrombbptable(double tval,segment *sgval,paircont* tbbptb,int tnumbbptb);
    void getfrombtntable(double tval,int *tlisteturn,double* tprobeturn,int *tindet,int *tindpos,int tnumeturn);

};

#endif 
