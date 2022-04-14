///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef INPUTDATA_H
#define INPUTDATA_H

#include "CommonPara.h"
#include "ParseInput.h"

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class InputData
{
  public:
    InputData();
    virtual ~InputData();

    void configureInput(const ParseInput &inputInfo);
    bool geneSSE(point3f *decstr,int numSeq);
    int getProteinType();
    sssegment *getSSE();
    int getNumSSE();
    ssef *get3StateSS();
    ssef *get8StateSS();
    char get3StateSSatPos(int i);
    char get8StateSSatPos(int i);
    int getNumBAB();
    abssinfo *getABSS();
    triplebab *getBABInfo();
    point3f **getFragConformation(int i);
    int getSeqLength();
    char *getSequence();
    char getResAtPos(int i);

  private:
    bool flagVerbose;
    int oldSegLength;
    char *sequence;
    ssef *ss3;
    ssef *ss8;   
    sssegment *sse;
    int numsse;
    int seqLength;
    double *dispos;
    double *disposcp;
    betastrand *bstrand;
    alphahelix *ahelix;
    abssinfo *abss;
    triplebab babinfo[30];
    int numbab;
    int numbstrand;
    int numahelix;
    int numabss;
    point3f **fragcont[20];
    int protType;
};

#endif // !defined(INPUTDATA_H)
