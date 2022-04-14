///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PROGRAMFUNCTIONS_H
#define PROGRAMFUNCTIONS_H

class ProgramFunctions
{
  public:
    void FoldDesignScore(char *libDir,char *dataDir,char *inputFile,char *energyWeightFile,
                         char *contactConstrFile,char *distConstrFile,char *pdbFile);
    void generateFragments(char *libDir,char *dataDir,char *inputFile,char *templateList,
                           bool zipFlag);
    void FoldDesignREMC(char *libDir,char *dataDir,char *inputFile,char *energyWeightFile,
                        char *contactConstrFile,char *distConstrFile,
                        char *fixedFileName,int cycles,int randomNum,bool fixFlag);
    void refineMainChain(char *dataDir,char *libDir,char *refname,char *fixname,int rannum,
                         double cuttime);
};

#endif
