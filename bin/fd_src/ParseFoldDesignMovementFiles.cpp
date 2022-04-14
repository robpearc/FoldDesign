///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "PrintFunc.h"
#include "ParseFoldDesignMovementFiles.h"

///////////////////////////////////////////////////////////////////////////////////////////////
// Construction/Destruction
///////////////////////////////////////////////////////////////////////////////////////////////

ParseFoldDesignMovementFiles::ParseFoldDesignMovementFiles()  // initialization of parameters/options
{
  int i;
  bborder=NULL;
  probeturn=NULL;
  listeturn=NULL;
  for(i=0;i<nosegdh;i++) topdhlist[i]=NULL;
  numtopdh=NULL;
  phiprob=NULL;
  psiprob=NULL;

}

ParseFoldDesignMovementFiles::~ParseFoldDesignMovementFiles()
{
  int i;
  if(phiprob)
  {
    for(i=0;i<11;i++)
    {
      delete[]phiprob[i];
    }
    delete[]phiprob;
    phiprob=NULL;
  }
  if(psiprob)
  {
    for(i=0;i<11;i++)
    {
      delete[]psiprob[i];
    }
    delete[]psiprob;
    psiprob=NULL;
  }
  if(bborder)
  {
    delete[]bborder;
    bborder=NULL;
  }
  if(probeturn)
  {
    delete[]probeturn;
    probeturn=NULL;
  }
  if(listeturn)
  {
    delete[]listeturn;
    listeturn=NULL;
  }
  for(i=0;i<nosegdh;i++)
  {
    if(topdhlist[i]) delete[]topdhlist[i];
  }
  if(numtopdh) 
  {
    delete[]numtopdh;
    numtopdh=NULL;
  }
}

void ParseFoldDesignMovementFiles::loadFiles(
  char *libDir,
  char *dataDir,
  InputData input
)
{
  inputInfo=input;
  int seqLength=inputInfo.getSeqLength();
  char fileName[STD_FILE_NAME_LENGTH+1];
  sprintf(fileName,"%s/phiprob.txt",libDir);
  loadPhiProb(fileName);
  sprintf(fileName,"%s/psiprob.txt",libDir);
  loadPsiProb(fileName);
  sprintf(fileName,"%s/betaorder.txt",libDir); //keep
  loadbbord(fileName);
  sprintf(fileName,"%s/turn.txt",dataDir); //keep
  loadeturn(fileName);
  sprintf(fileName,"%s/topdh.topdh",dataDir);
  loadTopdh(fileName,seqLength);
} 

bool ParseFoldDesignMovementFiles::loadPhiProb(
  char *fileName
)
{
  FILE *file;
  if((file=fopen(fileName,"rt"))==NULL)
  {
    printf("Unable to open phiprob %s\n",fileName);
    return false;
  }

  const int NUM_SS=11,NUM_PHI_BINS=360;
  int i,j,phi;
  char oneline[300],label[300];

  phiprob=new double*[NUM_SS];
  for(i=0;i<NUM_SS;i++)
  {
    phiprob[i]=new double[NUM_PHI_BINS];
  }

  for(i=0;i<NUM_SS;i++)
  {
    for(j=0;j<NUM_PHI_BINS;j++)
    {
      fgets(oneline,300,file);
      sscanf(oneline,"%s %d %lf",&label,&phi,&phiprob[i][j]);
    }
  }
  fclose(file);

  for(i=0;i<NUM_SS;i++)
  {
    for(j=1;j<NUM_PHI_BINS;j++)
    {
      phiprob[i][j]+=phiprob[i][j-1];
    }
  }

  return true;
}

bool ParseFoldDesignMovementFiles::loadPsiProb(
  char *fileName
)
{
  FILE *file;
  if((file=fopen(fileName,"rt"))==NULL)
  {
    printf("Unable to open psiprob %s\n",fileName);
    return false;
  }

  const int NUM_SS=11,NUM_PSI_BINS=360;
  int i,j,psi;
  char oneline[300],label[300];
  psiprob=new double*[NUM_SS];
  for(i=0;i<NUM_SS;i++)
  {
    psiprob[i]=new double[NUM_PSI_BINS];
  }

  for(i=0;i<NUM_SS;i++)
  {
    for(j=0;j<NUM_PSI_BINS;j++)
    {
      fgets(oneline,300,file);
      sscanf(oneline,"%s %d %lf",&label,&psi,&psiprob[i][j]);
    }
  }
  fclose(file);

  for(i=0;i<NUM_SS;i++){
    for(j=1;j<NUM_PSI_BINS;j++){
      psiprob[i][j]+=psiprob[i][j-1];
    }
  }

  return true;
}

bool ParseFoldDesignMovementFiles::loadbbord(
  char *filename
)
{
  FILE *file;
  file=fopen(filename,"rt");
  if(!file)
  {
    printf("cannot load border statics %s\n",filename);
    return false;
  }
  if(bborder)
  {
    delete[]bborder;
    bborder=NULL;
  }
  bborder=new double[1000];
  int i,j;
  float tval;
  char oneline[100];
  fgets(oneline,100,file);
  for(i=0;i<1000;i++)
  {
    fgets(oneline,100,file);
    sscanf(oneline,"%d %f",&j,&tval);
    bborder[i]=tval*30.0;
  }
  fclose(file);
  return true;
}

bool ParseFoldDesignMovementFiles::loadeturn(
  char *namefile
)
{
  FILE *file2;
  int i,tot;
  double totsum=0;
  file2=fopen(namefile,"rt");
  if(file2==NULL)
  {
    printf("Error when load turn file %s\n",namefile);
    return false;
  }
  if(probeturn)
  {
    delete[]probeturn;
    probeturn=NULL;
  }
  if(listeturn)
  {
    delete[]listeturn;
    listeturn=NULL;
  }
  char oneline[200];
  numeturn=0;
  fgets(oneline,200,file2);
  sscanf(oneline,"%d",&tot);
  probeturn=new double[tot];
  listeturn=new int[tot];
  int ti,tj;
  float tk;
  for(i=0;i<tot-3;i++)
  {
    fgets(oneline,200,file2);
    sscanf(oneline,"%d %d %f",&ti,&tj,&tk);
    if(tk>0)
    {
      listeturn[numeturn]=i;
      probeturn[numeturn]=tk;
      totsum+=tk;
      numeturn++;
    }
  }
  fclose(file2);
  if(numeturn==0) return true;
  for(i=0;i<numeturn;i++) probeturn[i]/=totsum;
  for(i=1;i<numeturn;i++) probeturn[i]+=probeturn[i-1];
  return true;
}

bool ParseFoldDesignMovementFiles::loadTopdh(
  char *filename,
  int numseq
)
{
  FILE *file;
  file=fopen(filename,"rt");
  if(!file)
  {
    printf("cannot load topdh %s\n",filename);
    return false;
  }
  ssef *ss3=inputInfo.get3StateSS();
  int i,j,k,tmpind;
  char tmpname[300];
  if(numtopdh) delete[]numtopdh;
  numtopdh=new int[numseq];
  for(i=0;i<nosegdh;i++)
  {
    if(topdhlist[i]) delete[]topdhlist[i];
    topdhlist[i]=new point3f[numseq];
  }

  for(i=0;i<numseq;i++)
  {
    fgets(tmpname,300,file);
    sscanf(tmpname,"%d %d",&j,&tmpind);
    numtopdh[i]=tmpind;
    for(k=0;k<numtopdh[i];k++)
    {
      fgets(tmpname,300,file);
      sscanf(tmpname,"%c %f %f %f %c %f %f %f %f %f %f %f %f %f %f %f %f %d %s %f",
             &topdhlist[k][i].residueid,
             &topdhlist[k][i].x,
             &topdhlist[k][i].y,
             &topdhlist[k][i].z,
             &topdhlist[k][i].stype,
             &topdhlist[k][i].phi,
             &topdhlist[k][i].leng,
             &topdhlist[k][i].angl,
             &topdhlist[k][i].tor[0],
             &topdhlist[k][i].len[0],
             &topdhlist[k][i].ang[0],
             &topdhlist[k][i].tor[1],
             &topdhlist[k][i].len[1],
             &topdhlist[k][i].ang[1],
             &topdhlist[k][i].tor[2],
             &topdhlist[k][i].len[2],
             &topdhlist[k][i].ang[2],
             &topdhlist[k][i].resind,
             &topdhlist[k][i].name,
             &topdhlist[k][i].prob);
      topdhlist[k][i].ss3=ss3[i].ss;
      topdhlist[k][i].aaa=ss3[i].res;
      topdhlist[k][i].iaa=(unsigned char)(bf.aminoid(topdhlist[k][i].aaa));
      if(topdhlist[k][i].iaa>19) topdhlist[k][i].iaa=5;
    }
  }
  fclose(file);
    
  return true;
}


