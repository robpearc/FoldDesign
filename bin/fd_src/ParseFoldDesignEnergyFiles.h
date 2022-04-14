///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PARSEFOLDDESIGNENERGYFILES_H
#define PARSEFOLDDESIGNENERGYFILES_H

#include "CommonPara.h"
#include "BasicFunc.h"
#include "GeometryCalc.h"
#include "InputData.h"
#define MAX_ENERGY_TERM_WEIGHTS 30
#define MAX_ENERGY_TERMS 20

class ParseFoldDesignEnergyFiles
{
  public:
    ParseFoldDesignEnergyFiles();
    virtual ~ParseFoldDesignEnergyFiles();

    bool loadFiles(char *libDir,char *dataDir,char *energyWeightFile,char *contactConstrFile,
                   char *distConstrFile,InputData inputInfo);

    bool loadEnergyWeights(char *libDir,char *weightFile);
    bool loadAllRamachandranFiles(char *libDir);
    bool loadAllHelixHelixPackingFiles(char *libDir);
    bool loadAllStrandStrandPackingFiles(char *libDir);
    bool loadAllHelixStrandPackingFiles(char *libDir);
    bool loadRamaEnergyFile(char *ramaFile,double **ramaEnergyArray);
    bool loadSSEPackingFilePhiPsiTheta(char *sseFile,double ***phiPsiThetaEnergyArray);
    bool loadSSEPackingFileDistTheta(char *sseFile,double **distThetaEnergyArray);
    bool loadSolSeq(char *solFile);
    bool loadDfire(char *filename);
    bool loadContactRestr(char *contactConstrFile,int seqLength);
    bool loadDistanceRestr(char *distancename,int seqnum);
    void getFragDistConstr(char *seqdat,InputData inputInfo);

  private:
    double weights[MAX_ENERGY_TERM_WEIGHTS];
    double enelist[MAX_ENERGY_TERMS];

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

    int ***countptr;
    double ***kbp;
    double *solweight;
    float ***contactConstr;
    float ***distanceConstr; 
    float ***distanceWeight;
    int ***distRestrType;
    pairaa *paa;

    friend class FoldDesignEnergyFunction;
};

#endif

