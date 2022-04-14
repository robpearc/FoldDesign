///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef GEOMETRYCALC_H
#define GEOMETRYCALC_H

#include "CommonPara.h"
#include "ParseGeometryFiles.h"

class GeometryCalc
{
  public:
    GeometryCalc();
    virtual ~GeometryCalc();
    bool configureGeometry(const ParseGeometryFiles &geometryParameters);
    bool tor2str(point3f *decstr,int seqnum,int type);
    bool tor2strp(point3f *decstr,int seqnum,int istart);
    void str2torp(point3f *decstr,int seqnum,int iStart,int iEnd);
    void str2tor(point3f *decstr,int seqnum,int type);
    bool tor2strsg(point3f *decstr,int seqnum);
    bool tor2strsg2(point3f *decstr,int seqnum);

  private:
    double ****sgpos;
    double ****sgpos2;
    bool flagVerbose;
};

#endif
