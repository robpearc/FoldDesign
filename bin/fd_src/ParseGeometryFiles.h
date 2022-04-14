///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PARSEGEOMETRYFILES_H
#define PARSEGEOMETRYFILES_H

class ParseGeometryFiles
{
  public:
    ParseGeometryFiles();
    virtual ~ParseGeometryFiles();

    bool loadFiles(char *libDir);
    bool loadSGpos(char *fileName);

  private:
    double ****sgpos;
    double ****sgpos2;

    friend class GeometryCalc;
};

#endif
