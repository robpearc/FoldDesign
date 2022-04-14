///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PARSEFOLDDESIGNMOVEMENTFILES_H
#define PARSEFOLDDESIGNMOVEMENTFILES_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "ParsePDB.h"
#include "BasicFunc.h"
#include "FoldDesignEnergyFunction.h"
#include "GeometryCalc.h"
#include "CommonPara.h"
#include "InputData.h"

class ParseFoldDesignMovementFiles
{
  public:
    ParseFoldDesignMovementFiles();
    virtual ~ParseFoldDesignMovementFiles();

    void loadFiles(char *libDir,char *dataDir,InputData input);
    bool loadPhiProb(char *fileName);
    bool loadPsiProb(char *fileName);
    bool loadTopdh(char *filename,int numseq);
    bool loadbbord(char *filename);
    bool loadeturn(char *namefile);

  private:
    BasicFunc bf;
    InputData inputInfo;
    point3f *topdhlist[nosegdh];
    double **phiprob;
    double **psiprob;
    double *bborder;
    double *probeturn;
    int *listeturn;
    int numeturn;
    int *numtopdh;

    friend class FoldDesignMovement;
};

#endif 
