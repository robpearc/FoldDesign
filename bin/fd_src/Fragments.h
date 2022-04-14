///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef FRAGMENTS_H
#define FRAGMENTS_H

#include "CommonPara.h"
#include "InputData.h"

class Fragments
{
  public:
    Fragments();
    virtual ~Fragments();   

    bool loadPISCES(char *selfile);
    bool loadSS2(char *ss2file,ssef *sssef,int seqlength);
    void normpSS2();
    void calcfragdh(int seqlength,char *outname,int segleng,int cuttop);
    void scoreFragments(seginfo **seginfall,int seqlength,InputData inputInfo,
                        char *dbfile,int segleng,int topleng,double *wtterm);
    bool generateFragments(char *dataDir,char *listName,InputData inputInfo,char *mtxDbName,
                           char *dbFeatName,char *outName,int fragSize);
    void calcfragdh(char *dataDir,int seqlength,char *outname,
                    int segleng,int cuttop);
    bool loadfragd(int numseq,int seglength,int ntop,char *filename);

  private:
    int numsel;
    int numhomo;
    int *homolist;
    int *dhcennum;
    point3f **fragcont[20];
    dihedral *fragdh[nosegdh];
    namep *namepdb;
    int oldseglength;
};

#endif
