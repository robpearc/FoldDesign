///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef FOLDDESIGNENERGYFUNCTION_H
#define FOLDDESIGNENERGYFUNCTION_H

#include "CommonPara.h"
#include "BasicFunc.h"
#include "GeometryCalc.h"
#include "InputData.h"
#include "ParseFoldDesignEnergyFiles.h"
#define MAX_ENERGY_TERM_WEIGHTS 30
#define MAX_ENERGY_TERMS 20

class FoldDesignEnergyFunction
{
  public:
    FoldDesignEnergyFunction();
    virtual ~FoldDesignEnergyFunction();

    void setEnergyFunctionParameters(const ParseFoldDesignEnergyFiles &energyParameters,
                                     GeometryCalc inputGeometry,InputData input);
    void estimateContactNum(int seqLength);

    void setContactParameters(int seqLength);
    void setLongestH();

    void calcAllAtomEnergy(point3f *decstr,int numseq,double *enelist);
    double calcHbondBackboneEnergy(point3f *decstr,int numseq);
    double calcSSConstrEnergy(point3f *decstr,int numseq);
    double calcRamaEnergy(point3f *decstr,int numseq);
    double calcBABMotifPenalty(point3f *decstr,int numseq);
    vector<double> calcStrandStrandEnergy(point3f *decstr,int numseq);
    vector<double> calcHelixHelixEnergy(point3f *decstr,int numseq);
    vector<double> calcHelixStrandEnergy(point3f *decstr,int numseq);
    double calcTotalEnergy(point3f *decstr,int numseq);
    void setHbondRamp(int curCycle,int numCycles);

  private:
    //---Object instantiation--->
    GeometryCalc geometry;
    InputData inputInfo;    
    BasicFunc bf;


    int expectedShort,expectedMed,expectedLong;
    double weights[MAX_ENERGY_TERM_WEIGHTS];
    double enelist[MAX_ENERGY_TERMS];
    bool hbondRamp;

    //---Ramachandran energy tables--->
    double **H_general_rama;
    double **E_general_rama;
    double **C_general_rama;
    double **H_specific_rama;
    double **E_specific_rama;
    double **C_specific_rama;
    double **G_specific_rama;
    double **I_specific_rama;
    double **B_specific_rama;
    double **S_specific_rama;
    double **T_specific_rama;

    //-----SSE packing energy tables---->
    double ***energyHHPackAngle2;
    double ***energyHHPackAngle3;
    double ***energyHHPackAngle4;
    double **energyHHPackDist2;
    double **energyHHPackDist3;
    double ***energySSPackAngle2;
    double ***energySSPackAngle3;
    double **energySSPackDist1;
    double **energySSPackDist2;
    double **energySSPackDist3;
    double **energySSPackDist4;
    double ***energyHSPackAngle2;
    double ***energyHSPackAngle3;
    double ***energyHSPackAngle4;
    double **energyHSPackDist2;
    double **energyHSPackDist3;
    double **energyHSPackDist4;

    double da,db,dc,dd,d8,dwell,d10;
    double longestdist;

    double ***kbp;
    double *solweight;
    float ***contactConstr;
    float ***distanceConstr; 
    float ***distanceWeight;
    int ***distRestrType;
    pairaa *paa;
};

#endif

