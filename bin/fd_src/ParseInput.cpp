///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "ParseInput.h"
#include "BasicFunc.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ParseInput::ParseInput()
{
  int i;

  seqLength=0;
  sequence=NULL;
  ss3=NULL;
  ss8=NULL;
  sse=NULL;
  //srand((int)time(0));

  for(i=0;i<20;i++)
  {
    fragcont[i]=NULL;
  }	
  dispos=NULL;
  disposcp=NULL;
  bstrand=NULL;
  ahelix=NULL;
  abss=NULL;
  numbstrand=0;
  numahelix=0;
  numabss=0;
  oldSegLength=0;
  protType=2; 
  fixPos=NULL;
  numFixed=0;
  flagVerbose=true;
}

ParseInput::~ParseInput()
{
  int i,j,k;
  if(sequence)
  {
    delete[]sequence;
    sequence=NULL;
  }
  if(ss3)
  {
    delete[]ss3;
    ss3=NULL;
  }
  if(ss8)
  {
    delete[]ss8;
    ss8=NULL;
  }	
  if(sse)
  {
    delete[]sse;
    sse=NULL;
  }

  if(dispos)
  {
    delete[]dispos;
    dispos=NULL;
  }
  if(disposcp)
  {
    delete[]disposcp;
    disposcp=NULL;
  }
  if(bstrand)
  {
    delete[]bstrand;
    bstrand=NULL;
  }
  if(ahelix)
  {
    delete[]ahelix;  
    ahelix=NULL;
  }
  if(abss)
  {
    delete[]abss;
    abss=NULL;
  }

  if(oldSegLength>0)
  {
    if(fragcont[0])
    {
      for(j=0;j<oldSegLength*topno;j++)
      {
        delete[]fragcont[0][j];
      }
      delete[]fragcont[0];
      fragcont[0]=NULL;
    }
  }
  for(i=0;i<20;i++)
  {
    if(fragcont[i])
    {
      for(j=0;j<(i+1)*topno;j++)
      {
        delete[]fragcont[i][j];
      }
      delete[]fragcont[i];
    }
  }
}

bool ParseInput::configureInputFoldDesign(
  char *dataDir,
  char *inputFile
)
{
  bool flagUnassigned=false;
  bool flagSuccess;

  flagSuccess=loadInput(dataDir,inputFile,flagUnassigned);
  loadFrags(dataDir);
  if(flagUnassigned)
  {
    assignSecStructFromFrags();
  }

  if(!flagSuccess)
  {
    printf("Fatal Error: There were too many undefined secondary structure values in the input. At most 18 consecutive secondary structure values may be assigned as X!\n");
    return false;
  }

  modifySS();
  geneSSE();
  getAhelix();
  getBstrand();     
  getABseq();
  calcBAB();
  determineProteinType();

  return true;
}

bool ParseInput::configureInputFragments(
  char *dataDir,
  char *inputFile
)
{
  bool flagUnassigned=false;
  bool flagSuccess;

  flagSuccess=loadInput(dataDir,inputFile,flagUnassigned);

  modifySS();
  geneSSE();
  getAhelix();
  getBstrand();
  getABseq();
  calcBAB();
  determineProteinType();

  return flagSuccess;
  //normpSS1();
}


bool ParseInput::loadInput(
  char *dataDir,
  char *inputFile,
  bool &flagUnassigned
)
{
  char fileName[STD_ARRAY_SIZE];
  sprintf(fileName,"%s/%s",dataDir,inputFile);
  FILE* file=fopen(fileName,"r");
  if(!file)
  {
    fprintf(stderr,"Error in file %s function %s line %d, cannot open file %s\n",__FILE__,
            __FUNCTION__,__LINE__,inputFile);
    return false;
  }

  char line[STD_ARRAY_SIZE];
  char residue,secStruct,fixChar;
  seqLength=0;
  numFixed=0;
  int allocation=100;
  int allocationFixed=100;
  int seqAllocation=65525;
  fixChar='R';

  if(ss3)
  {
    delete[]ss3;
    ss3=NULL;
  }

  if(ss8)
  {
    delete[]ss8;
    ss8=NULL;
  }
  if(fixPos)
  {
    delete[]fixPos;
    fixPos=NULL;
  }
  if(sequence)
  {
    delete[]sequence;
    sequence=NULL;
  }

  ss3=new ssef[allocation];
  ss8=new ssef[allocation];
  fixPos=new int[allocationFixed];
  sequence=new char[seqAllocation];

  while(fgets(line,STD_ARRAY_SIZE,file))
  {
    sscanf(line,"%c %c",&residue,&secStruct);
    
    if(0)//fixChar=='F')
    {
      fixPos[numFixed]=seqLength;
      numFixed++;
    }

    if(!(residue>='A' && residue<='Z'))
    {
      sequence[seqLength]='V';
      ss3[seqLength].res='V';
      ss8[seqLength].res='V';
      //fprintf(stderr,"Error in sequence code\n");
      //delete[]seqdata;
      //seqdata=NULL;
      //return false;
    }
    else if(residue>='A' && residue<='Z')
    {
      sequence[seqLength]=residue;
      ss3[seqLength].res=residue;
      ss8[seqLength].res=residue;
    }

    if(secStruct=='H' || secStruct=='h' || secStruct=='G' || secStruct=='I')
    {
      ss3[seqLength].ss='H';
      ss3[seqLength].a[0]=0.000;
      ss3[seqLength].a[1]=0.999;
      ss3[seqLength].a[2]=0.000;
    }
    else if(secStruct=='E' || secStruct=='e' || secStruct=='B')
    {
      ss3[seqLength].ss='E';
      ss3[seqLength].a[0]=0.000;
      ss3[seqLength].a[1]=0.000;
      ss3[seqLength].a[2]=0.999;

    }
    else if(secStruct=='C' || secStruct=='c' || secStruct=='T' || secStruct=='S')
    {
      ss3[seqLength].ss='C';
      ss3[seqLength].a[0]=0.999;
      ss3[seqLength].a[1]=0.000;
      ss3[seqLength].a[2]=0.000;
    }
    else
    {
      ss3[seqLength].ss='X';
      ss3[seqLength].a[0]=0.000;
      ss3[seqLength].a[1]=0.000;
      ss3[seqLength].a[2]=0.000;
      flagUnassigned=true;
    }

    if(secStruct=='H')
    {
      ss8[seqLength].ss='H';
      ss8[seqLength].a[0]=0.000;
      ss8[seqLength].a[1]=0.999;
      ss8[seqLength].a[2]=0.000;
    }
    else if(secStruct=='h')
    {
      ss8[seqLength].ss='h';
      ss8[seqLength].a[0]=0.000;
      ss8[seqLength].a[1]=0.500;
      ss8[seqLength].a[2]=0.000;
    }
    else if(secStruct=='G')
    {
      ss8[seqLength].ss='G';
      ss8[seqLength].a[0]=0.000;
      ss8[seqLength].a[1]=0.500;
      ss8[seqLength].a[2]=0.000;
    }
    else if(secStruct=='I')
    {
      ss8[seqLength].ss='I';
      ss8[seqLength].a[0]=0.000;
      ss8[seqLength].a[1]=0.500;
      ss8[seqLength].a[2]=0.000;
    }
    else if(secStruct=='E')
    {
      ss8[seqLength].ss='E';
      ss8[seqLength].a[0]=0.000;
      ss8[seqLength].a[1]=0.000;
      ss8[seqLength].a[2]=0.999;
    }
    else if(secStruct=='e')
    {
      ss8[seqLength].ss='e';
      ss8[seqLength].a[0]=0.000;
      ss8[seqLength].a[1]=0.000;
      ss8[seqLength].a[2]=0.500;
    }
    else if(secStruct=='B')
    {
      ss8[seqLength].ss='B';
      ss8[seqLength].a[0]=0.000;
      ss8[seqLength].a[1]=0.000;
      ss8[seqLength].a[2]=0.500;
    }
    else if(secStruct=='C')
    {
      ss8[seqLength].ss='C';
      ss8[seqLength].a[0]=0.000;
      ss8[seqLength].a[1]=0.000;
      ss8[seqLength].a[2]=0.999;
    }
    else if(secStruct=='c')
    {
      ss8[seqLength].ss='c';
      ss8[seqLength].a[0]=0.000;
      ss8[seqLength].a[1]=0.000;
      ss8[seqLength].a[2]=0.500;
    }
    else if(secStruct=='T')
    {
      ss8[seqLength].ss='T';
      ss8[seqLength].a[0]=0.000;
      ss8[seqLength].a[1]=0.000;
      ss8[seqLength].a[2]=0.500;
    }
    else if(secStruct=='S')
    {
      ss8[seqLength].ss='S';
      ss8[seqLength].a[0]=0.000;
      ss8[seqLength].a[1]=0.000;
      ss8[seqLength].a[2]=0.500;
    }
    else
    {
      ss8[seqLength].ss='X';
      ss8[seqLength].a[0]=0.000;
      ss8[seqLength].a[1]=0.000;
      ss8[seqLength].a[2]=0.000;
    }
    seqLength++;

    if(seqLength>=allocation)
    {
      allocation*=2;
      ss3=(ssef *)realloc(ss3,allocation*sizeof(ssef));
      ss8=(ssef *)realloc(ss8,allocation*sizeof(ssef));
    }
    if(numFixed>=allocationFixed)
    {
      allocationFixed*=2;    
      fixPos=(int *)realloc(fixPos,allocationFixed*sizeof(int));
    }
  }
  fclose(file);

  if(seqLength!=0)
  { 
    ss3=(ssef *)realloc(ss3,seqLength*sizeof(ssef));
    ss8=(ssef *)realloc(ss8,seqLength*sizeof(ssef));
  }
  if(numFixed!=0)
  {
    fixPos=(int *)realloc(fixPos,numFixed*sizeof(int));
  }

  sequence[seqLength]='\0';

  return true;
}

bool ParseInput::loadFrags(
  char *dataDir
)
{
  int i,j,k,l,tmpind,istart;
  FILE *file;
  char tmpname[300],tmpstr[8];
  int totfrag=20;

  for(i=0;i<totfrag;i++)
  {
    if(fragcont[i])
    {
      for(j=0;j<(i+1)*topno;j++) delete[]fragcont[i][j];
      delete[]fragcont[i];
    }
  }
  for(i=0;i<totfrag;i++)
  {
    fragcont[i]=new point3f*[(i+1)*topno];
    for(j=0;j<(i+1)*topno;j++) fragcont[i][j]=new point3f[seqLength-(i+1)+1];
  }
  for(i=0;i<totfrag;i++)
  {
    sprintf(tmpname,"%s/%dseqfra.topse",dataDir,i+1);
    file=fopen(tmpname,"rt");
    if(!file)
    {
      printf("Cannot load fragments %s\n",tmpname);
      continue;
    }
    for(j=0;j<seqLength-(i+1)+1;j++)
    {
      for(l=0;l<topno;l++)
      {
        istart=l*(i+1);
        fgets(tmpname,300,file);
        sscanf(tmpname,"%s %d %f",tmpstr,&tmpind,&fragcont[i][istart][j].ptsg.x);
        for(k=0;k<(i+1);k++)
        {
          fgets(tmpname,300,file);
          sscanf(tmpname,
                 "%c %f %f %f %c %f %f %f %f %f %f %f %f %f %f %f %f",
                 &fragcont[i][istart+k][j].residueid,
                 &fragcont[i][istart+k][j].x,
                 &fragcont[i][istart+k][j].y,
                 &fragcont[i][istart+k][j].z,
                 &fragcont[i][istart+k][j].stype,
                 &fragcont[i][istart+k][j].phi,
                 &fragcont[i][istart+k][j].leng,
                 &fragcont[i][istart+k][j].angl,
                 &fragcont[i][istart+k][j].tor[0],
                 &fragcont[i][istart+k][j].len[0],
                 &fragcont[i][istart+k][j].ang[0],
                 &fragcont[i][istart+k][j].tor[1],
                 &fragcont[i][istart+k][j].len[1],
                 &fragcont[i][istart+k][j].ang[1],
                 &fragcont[i][istart+k][j].tor[2],
                 &fragcont[i][istart+k][j].len[2],
                 &fragcont[i][istart+k][j].ang[2]);
          fragcont[i][istart+k][j].resind=tmpind+k;
          fragcont[i][istart+k][j].ss3=ss3[j+k].ss;
          strcpy(fragcont[i][istart+k][j].name,tmpstr);
        }
        if(j==seqLength-(i+1) && l==topno-1)
        {
          printf("loading fragments %2d %d %d %d %f %f %f\n",
                 i+1,j,l,i+1,
                 fragcont[i][istart+i][j].tor[2],
                 fragcont[i][istart+i][j].len[2],
                 fragcont[i][istart+i][j].ang[2]);
        }
      }
    }
    fclose(file);
  }

  return true;
}


bool ParseInput::assignSecStructFromFrags()
{
  int i,j,k;
  int largestFrag=20;
  int istart;
  for(i=0;i<seqLength;i++)
  {
    if(ss3[i].ss=='X')
    {
      int numX=1;
      for(j=1;j<seqLength-i;j++)
      {
        if(ss3[i+j].ss=='X')
        {
          numX++;
        }
        else
        {
          break;
        }
      }
      if(numX>largestFrag-2)
      {
        return false;
      }
      int E_count=0;
      int H_count=0;
      int C_count=0;
      for(int num=numX+1; num<largestFrag; num++)
      {
        for(int pos=i-num;pos<=i;pos++)
        {
          if(i>=0 && i<seqLength-num)
          {
            for(int index=0;index<topno;index++)
            {
              istart=index*(num+1)+i-pos;
              if(fragcont[num][istart][pos].stype=='E')
              {
                E_count++;
              }
              else if(fragcont[num][istart][pos].stype=='H')
              {
                H_count++;
              }
              else if(fragcont[num][istart][pos].stype=='C')
              {
                C_count++;
              }
            }
          }
        }
      }
      double probE = E_count/(E_count+H_count+C_count);
      double probH = H_count/(E_count+H_count+C_count);
      double probC = C_count/(E_count+H_count+C_count);
      ss3[i].a[0] = probC;
      ss3[i].a[1] = probH;
      ss3[i].a[2] = probE;
      ss8[i].a[0] = probC;
      ss8[i].a[1] = probH;
      ss8[i].a[2] = probE;

      if(probE>probH && probE>probC)
      {
        ss3[i].ss='E';
        ss8[i].ss='E';
      }
      else if(probH>probE && probH>probC)
      {
        ss3[i].ss='H';
        ss8[i].ss='H';
      }
      else
      {
        ss3[i].ss='C';
        ss8[i].ss='C';
      }
    }
  }

  return true;
}

void ParseInput::pdb2seq(
  boneinfo *bb,
  int numbb
)
{
  int seqAllocation=65525;
  int i;
  int indres;
  BasicFunc bf;

  if(!sequence)
  {
    sequence=new char[seqAllocation];
  }
  for(i=0;i<numbb;i++)
  {
    indres=bf.aminoid(bb[i].resid);
    if(indres>19) indres=5;
    sequence[i]=aad1[indres];
  }
  sequence[numbb]='\0';
  seqLength=numbb;
}

/*
void ParseInput::normpSS8()
{
  if(!ss1) return;
  int i;
  float t;
  for(i=0;i<numss1;i++)
  {
    t=ss1[i].a[0]+ss1[i].a[1]+ss1[i].a[2];
    if(t<epsilon) t=1;
    ss1[i].a[0]/=t;
    ss1[i].a[1]/=t;
    ss1[i].a[2]/=t;
  }
}
*/

/*
void ParseInput::normpSS3()
{
  if(!ss2) return;
  int i;
  float t;
  for(i=0;i<numss2;i++)
  {
    t=ss3[i].a[0]+ss3[i].a[1]+ss3[i].a[2];
    if(t<epsilon) t=1;
    ss3[i].a[0]/=t;
    ss3[i].a[1]/=t;
    ss3[i].a[2]/=t;
  }
}
*/
void ParseInput::modifySS()
{
  if(!ss3) return;
  int i,j,k;
  for(i=1;i<seqLength-1;i++)
  {
    if(ss3[i].ss=='H' && ss3[i-1].ss=='E' && ss3[i+1].ss=='E' && ss3[i].a[1]!=0.999 &&
       ss3[i].a[0]<ss3[i].a[2])
    {
      ss3[i].ss='E';
      ss8[i].ss='E';
    }
    else if(ss3[i].ss=='H' && ss3[i-1].ss=='E' && ss3[i+1].ss=='E' && ss3[i].a[1]!=0.999)
    {
      ss3[i].ss='C';
      ss8[i].ss='C';
    }
  }

  for(i=0;i<seqLength;i++)
  {
    j=i;
    int index;
    if(ss3[j].ss=='E')
    {
      index=2;
    }
    else if(ss3[j].ss=='H')
    {
      index=1;
    }
    else if(ss3[j].ss=='C')
    {
      index=0;
    }
    bool flagIter = false;
    while(j<seqLength && ss3[j].ss==ss3[i].ss && ss3[j].a[index]!=0.999 && 
          ss3[i].a[index]!=0.999)
    {
      flagIter = true;
      j++;
    }
    if(flagIter==false)
    {
      j++;
    }
    if((j-i==2 || j-i==1) && ss3[i].ss=='H')
    {
      for(k=i;k<j;k++)
      { 
        ss3[k].ss='C';
        ss8[k].ss='C';
      }
    }
    i=j-1;
  }
}

bool ParseInput::geneSSE()
{
  if(!ss3)
    return false;
  int i,j;
  if(sse)
  {
    delete[]sse;
    sse=NULL;
  }
  sse=new sssegment[seqLength];
  numsse=0;
  for(i=0;i<seqLength;i++)
  {
    char s1=ss3[i].ss;
    char s2=ss3[i].ss;

    j=i;
    while(j<seqLength && s1==s2)
    {
      s1=ss3[j].ss;
      j++;
    }
    sse[numsse].init=i;
    if(j>=seqLength-1)
    {
      sse[numsse].term=j-1;
    }
    else
    {
      sse[numsse].term=j-2;
    }
    sse[numsse].ss=s2;
    numsse++;
    if(j>=seqLength-1)
    {
      break;
    }
    else
    {
      i=j-2;
    }
  }
  if(numsse!=0)
    sse=(sssegment *)realloc(sse,numsse*sizeof(sssegment));

  ////////know the proba of disconnect
  if(dispos)
  {
    delete[]dispos;
    dispos=NULL;
  }
  dispos=new double[seqLength];
  for(i=0;i<seqLength-1;i++)
  {
    dispos[i]=ss3[i].a[0]+ss3[i+1].a[0];
  }
  for(i=0;i<numsse-1;i++)
  {
    dispos[sse[i].term]+=2.0;
  }
  dispos[seqLength-1]=0;
  if(disposcp)
  {
    delete[]disposcp;
    disposcp=NULL;
  }
  disposcp=new double[seqLength];
  for(i=0;i<seqLength;i++)
  {
    disposcp[i]=dispos[i];
  }

  return true;
}

void ParseInput::removeSingletonHelices(
  point3f *decstr,
  int numSeq
)
{
  int i,j;
  if(!dispos)
  {
    dispos=new double[numSeq];
  }
  char *tmpss;
  tmpss=new char[numSeq];
  for(i=0;i<numSeq;i++)
  {
    if(decstr[i].stype=='H')      tmpss[i]='H';
    else if(decstr[i].stype=='E') tmpss[i]='E';
    else                          tmpss[i]='C';
  }

  //remove single H or form at least three H, don't modify two H
  i=0;
  while(i<numSeq)
  {
    j=i;
    while(j<numSeq && tmpss[j]==tmpss[i]) j++;
    if(j-i>=2) i=j;
    else if(j-i==1 && tmpss[i]=='H')
    {
      if(i<numSeq-2 && tmpss[i+2]=='H') tmpss[i+1]='H';
      else if(i>1 && tmpss[i-2]=='H')
      {
        tmpss[i-1]='H';
        i-=2;
      }
      else tmpss[i]='C';
    }
    else i=j;
  }
  for(i=0;i<numSeq;i++)
  {
    dispos[i]=0;
    if(tmpss[i]=='C') dispos[i]=2;
    else if(i+1<numSeq && tmpss[i]!='C' && tmpss[i+1]=='C') dispos[i]=2;
  }
  for(i=0;i<numSeq;i++) decstr[i].stype=tmpss[i];
  delete[]tmpss;
  tmpss=NULL;
}

void ParseInput::getAhelix()
{
  int i;

  if(ahelix)
  {
    delete[]ahelix;
    ahelix=NULL;
  }
  numahelix=0;
  ahelix=new alphahelix[numsse];
  for(i=0;i<numsse;i++)
  {
    if(sse[i].ss=='H')
    {
      ahelix[numahelix].seg.indss2=i;
      ahelix[numahelix].seg.init=sse[i].init;
      ahelix[numahelix].seg.term=sse[i].term;
      numahelix++;
    }
  }
  if(numahelix!=0)
  ahelix=(alphahelix *)realloc(ahelix,numahelix*sizeof(alphahelix));
}

void ParseInput::getBstrand()
{
  if(bstrand)
  {
    delete[]bstrand;
    bstrand=NULL;
  }
  numbstrand=0;
  bstrand=new betastrand[numsse];
  int i;
  for(i=0;i<numsse;i++)
  {
    if(sse[i].ss=='E')
    {
      bstrand[numbstrand].seg.indss2=i;
      bstrand[numbstrand].seg.init=sse[i].init;
      bstrand[numbstrand].seg.term=sse[i].term;
      numbstrand++;
    }
  }
  if(numbstrand!=0)
  bstrand=(betastrand *)realloc(bstrand,numbstrand*sizeof(betastrand));
}

void ParseInput::getABseq()
{
  if(abss)
  {
    delete[]abss;
    abss=NULL;
  }
  numabss=0;
  if(numahelix+numbstrand==0) return;
  abss=new abssinfo[numahelix+numbstrand];
  int indi,indj;
  indi=0;indj=0;
  while(indi<numahelix || indj<numbstrand)
  {
    if(indj==numbstrand && indi<numahelix)
    {
      abss[numabss].ss='H';
      abss[numabss].inds=indi;
      abss[numabss].seg=ahelix[indi].seg;
      indi++;
      numabss++;
    }
    else if(indi==numahelix && indj<numbstrand)
    {
      abss[numabss].ss='E';
      abss[numabss].inds=indj;
      abss[numabss].seg=bstrand[indj].seg;
      indj++;
      numabss++;
    }
    else if(ahelix[indi].seg.init<bstrand[indj].seg.init)
    {
      abss[numabss].ss='H';
      abss[numabss].inds=indi;
      abss[numabss].seg=ahelix[indi].seg;
      indi++;
      numabss++;
    }
    else if(ahelix[indi].seg.init>bstrand[indj].seg.init)
    {
      abss[numabss].ss='E';
      abss[numabss].inds=indj;
      abss[numabss].seg=bstrand[indj].seg;
      indj++;
      numabss++;
    }
  }
}

void ParseInput::calcBAB()
{
  int i,j,k;
  int dellen;
  int cutlen=10;
  numbab=0;
  for(i=0;i<numabss;i++)
  {
    if(abss[i].ss=='E')
    {
      j=i+1;
      k=i+2;
      if(k>=numabss) continue;
      dellen=abs(abss[i].seg.term-abss[i].seg.init-abss[k].seg.term+abss[k].seg.init);
      if(abss[j].ss=='H' && abss[k].ss=='E' && dellen<=cutlen)
      {
        babinfo[numbab].a=i;
        babinfo[numbab].b=j;
        babinfo[numbab].c=k;
        numbab++;
      }
      k=i+3;
      if(k>=numabss) continue;
      dellen=abs(abss[i].seg.term-abss[i].seg.init-abss[k].seg.term+abss[k].seg.init);
      if(abss[j].ss=='H' && abss[k].ss=='E' && dellen<=cutlen)
      {
        babinfo[numbab].a=i;
        babinfo[numbab].b=j;
        babinfo[numbab].c=k;
        numbab++;
      }
      j=i+2;
      dellen=abs(abss[i].seg.term-abss[i].seg.init-abss[k].seg.term+abss[k].seg.init);
      if(abss[j].ss=='H' && abss[k].ss=='E' && dellen<=cutlen)
      {
        babinfo[numbab].a=i;
        babinfo[numbab].b=j;
        babinfo[numbab].c=k;
        numbab++;
      }
      if(numbab>26) break;
    }
  }
  i=-1;
  j=-1;
  if(numabss>=4)
  {
    if(abss[0].ss=='E') i=0;
    else if(abss[1].ss=='E') i=1;
    if(abss[numabss-1].ss=='E') j=numabss-1;
    else if(abss[numabss-2].ss=='E') j=numabss-2;
  }
  babinfo[numbab].a=i;
  babinfo[numbab].b=j;
}

void ParseInput::determineProteinType()
{
  const int ALPHA_PROTEIN=1;
  const int NOT_ALPHA_PROTEIN=2;

  if(numbstrand<2)
  {
    protType=ALPHA_PROTEIN;
  }
  else
  {
    protType=NOT_ALPHA_PROTEIN;
  }

}

int ParseInput::getSeqLength()
{
  return seqLength;
}

char ParseInput::getSeq(
  int i
)
{
  return sequence[i];
}
