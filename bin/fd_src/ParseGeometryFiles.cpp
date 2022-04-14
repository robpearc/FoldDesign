///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "CommonPara.h"
#include "ParseGeometryFiles.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////lo////////////////

ParseGeometryFiles::ParseGeometryFiles()
{
  sgpos=NULL;
  sgpos2=NULL;
}

ParseGeometryFiles::~ParseGeometryFiles()
{
  int i,j,k;

  if(sgpos)
  {  
    for(i=0;i<6;i++)
    {
      for(j=0;j<20;j++)
      {
        for(k=0;k<180;k++)
        {
          delete[]sgpos[i][j][k];
        }
        delete[]sgpos[i][j];
      }
      delete[]sgpos[i];
    }
    delete[]sgpos;
    sgpos=NULL;
  }
}

bool ParseGeometryFiles::loadFiles(
  char *libDir
)
{
  char fileName[STD_FILE_NAME_LENGTH+1];
  sprintf(fileName,"%s/newsgdistriaanew72.txt",libDir);
  loadSGpos(fileName); 
}

bool ParseGeometryFiles::loadSGpos(
  char *fileName
)
{
  FILE *file;
  file=fopen(fileName,"rt");
  if(!file)
  {
    printf("no sgpos2 file %s\n",fileName);
    return false;
  }

  int i,j,k;
  int ndim=72;
  char line[STD_ARRAY_SIZE];
  if(!sgpos)
  {
    sgpos=new double***[6];
    for(i=0;i<6;i++)
    {
      sgpos[i]=new double**[20];
      for(j=0;j<20;j++)
      {
        sgpos[i][j]=new double*[180];
        for(k=0;k<180;k++)
        {
          sgpos[i][j][k]=new double[180];
        }		
      }		
    }	
  }
  for(i=0;i<20;i++)
  {
    for(j=0;j<ndim;j++)
    {
      for(k=0;k<ndim;k++)
      {
        fgets(line,STD_ARRAY_SIZE,file);
        sscanf(line,"%lf %lf %lf %lf %lf %lf",&sgpos[0][i][j][k],&sgpos[1][i][j][k],
               &sgpos[2][i][j][k],&sgpos[3][i][j][k],&sgpos[4][i][j][k],
               &sgpos[5][i][j][k]);
      }
    }
    //if(i==19)
    //printf("sgposition %d %f %f %f %f %f %f\n",i,sgpos2[0][i][ndim-1][ndim-1],
    //       sgpos2[1][i][ndim-1][ndim-1],sgpos2[2][i][ndim-1][ndim-1],
    //       sgpos2[3][i][ndim-1][ndim-1],sgpos2[4][i][ndim-1][ndim-1],
    //       sgpos2[5][i][ndim-1][ndim-1]);
  }
  fclose(file);

  return true;
}

