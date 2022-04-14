///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PARSEINPUT_H
#define PARSEINPUT_H

#include "CommonPara.h"

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class ParseInput
{
  public:
    ParseInput();
    virtual ~ParseInput();


    bool configureInputFoldDesign(char *dataDir,char *inputFile);
    bool configureInputFragments(char *dataDir,char *inputFile);

    bool loadInput(char *dataDir,char *inputFile,bool &flagUnassigned);
    bool loadFrags(char *dataDir); //Load seq.txt
    bool assignSecStructFromFrags();

    void pdb2seq(boneinfo *bb,int numbb);
    void modifySS();
    bool geneSSE();
    void removeSingletonHelices(point3f *decstr,int numSeq);
    void getAhelix();
    void getBstrand();
    void getABseq();
    void calcBAB();
    void determineProteinType();
    int getSeqLength();
    char getSeq(int i);

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
    int *fixPos;
    int numFixed;    

    friend class InputData;
};

#endif // !defined(PARSEINPUT_H)
