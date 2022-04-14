///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PRINTFUNC_H
#define PRINTFUNC_H

#include "CommonPara.h"

class PrintFunc
{
  public:
    PrintFunc();
    void printSS(const point3f *decstr,const int resnum);
    bool writetradesign(point3f *decstr, int numseq, char *outname, int iind,
                        double denergy,int ss_score,int ftype);
    bool writetra(point3f *decstr, int numseq, char *outname, int iind,
                  double denergy,int ftype);
    bool writepdb(point3f *decstr, int numseq, int *bbres,
                  const char *outname, int iind, double denergy, double *eneterms,
                  int numene, int atype, int ftype);

    bool writepdb(point3f *decstr,int numseq,char *outname,int iind,double denergy,
                  double *eneterms,int numene,int atype,int ftype);

    bool writepdb(point3f *decstr,int numseq,int *bbres,char *outname,int iind,
                  double denergy,double *eneterms,int numene,int atype,int ftype);

    virtual ~PrintFunc();
};

#endif 
