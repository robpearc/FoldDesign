///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "InputData.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

InputData::InputData()
{
  int i;
  sequence=NULL;
  ss3=NULL;
  ss8=NULL;
  sse=NULL;
  numsse=0;
  seqLength=0;
  dispos=NULL;
  disposcp=NULL;
  bstrand=NULL;
  ahelix=NULL;
  abss=NULL;
  numbstrand=0;
  numahelix=0;
  numabss=0;
  for(i=0;i<20;i++)
  {
    fragcont[i]=NULL;
  }
  oldSegLength=0;
  protType=2; 
}

InputData::~InputData()
{
}

void InputData::configureInput(
  const ParseInput &inputInfo
)
{
  int i;
  seqLength=inputInfo.seqLength;
  sequence=inputInfo.sequence;
  ss3=inputInfo.ss3;
  ss8=inputInfo.ss8;
  sse=inputInfo.sse;
  dispos=inputInfo.dispos;
  disposcp=inputInfo.disposcp;
  bstrand=inputInfo.bstrand;
  ahelix=inputInfo.ahelix;

  numsse=inputInfo.numsse;
  numbstrand=inputInfo.numbstrand;
  numahelix=inputInfo.numahelix;
  numabss=inputInfo.numabss; 
  for(i=0;i<20;i++)
  {
    fragcont[i]=inputInfo.fragcont[i];
  }
  for(i=0;i<30;i++)
  {
    babinfo[i]=inputInfo.babinfo[i];

  }
  oldSegLength=inputInfo.oldSegLength;
  protType=inputInfo.protType;
}

bool InputData::geneSSE(
  point3f *decstr,
  int numSeq
)
{
  int i,j;
  if(sse)
  {
    delete[]sse;
    sse=NULL;
  }
  sse=new sssegment[numSeq];
  numsse=0;
  for(i=0;i<numSeq;i++)
  {
    j=i;
    while(j<numSeq && decstr[j].ssm==decstr[i].ssm)
    {
      j++;
    }
    sse[numsse].init=i;
    sse[numsse].term=j-1;
    sse[numsse].ss=decstr[i].ssm;
    numsse++;
    i=j-1;
  }
  if(numsse!=0)
  sse=(sssegment *)realloc(sse,numsse*sizeof(sssegment));

  return true;
}

int InputData::getProteinType()
{
  return protType;
}

sssegment *InputData::getSSE()
{
  return sse;
}

int InputData::getSeqLength()
{
  return seqLength;
}

int InputData::getNumSSE()
{
  return numsse;
}

ssef *InputData::get3StateSS()
{
  return ss3;
}

ssef *InputData::get8StateSS()
{
  return ss8;
}

char InputData::get3StateSSatPos(
  int i
)
{
  return ss3[i].ss;
}

char InputData::get8StateSSatPos(
  int i
)
{
  return ss8[i].ss;
}

int InputData::getNumBAB()
{
  return numbab;
}

abssinfo *InputData::getABSS()
{
  return abss;
}

triplebab *InputData::getBABInfo()
{
  return babinfo;
}

point3f **InputData::getFragConformation(
  int i
)
{
  return fragcont[i];
}

char *InputData::getSequence()
{
  return sequence;
}

char InputData::getResAtPos(
  int i
)
{
  return sequence[i];
}

