///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "PrintFunc.h"
#include "FoldDesignMovement.h"

///////////////////////////////////////////////////////////////////////////////////////////////
// Construction/Destruction
///////////////////////////////////////////////////////////////////////////////////////////////

FoldDesignMovement::FoldDesignMovement()  // initialization of parameters/options
{
  int i,j;
 
  for(i=0;i<3;i++)
  {
    summcsub[i]=0;
    summctop[i]=0;
    summcphi[i]=0;
    summcpsi[i]=0;
    summcome[i]=0;
    summclen[i]=0;
    summcang[i]=0;
    summcccd[i]=0;
    summcLMP[i]=0;
    summcmid[i]=0;
    summcbbp[i]=0;
    summcaaa[i]=0;
    summcbtn[i]=0;
    summcsft[i]=0;
    summctra[i]=0;
  }

  for(i=0;i<20;i++) summctot[i]=0;
  temparray=NULL;
  ttemparray=NULL;
  tttemparray=NULL;
  tmstr=NULL;
  tmstr2=NULL;
  fixed_decstr=NULL;
  alphaind=NULL;
  betaind=NULL;
  alphasind=NULL;
  numbbptb=0;
  numeturn=0;
  bbptb=NULL;
  bborder=NULL;
  probeturn=NULL;
  listeturn=NULL;
  
  for(i=0;i<nosegdh;i++) topdhlist[i]=NULL;
  numtopdh=NULL;

  phiprob=NULL;
  psiprob=NULL;

  for(i=0;i<20;i++) fragcont[i]=NULL;
  

  indCyc=0;
  ntopfrag=1;
  flagres=0;
  flagLocalMove=false;
  flagcaca=false;
  flagVerbose=true;
  hbondWeightRamp=1;
}

FoldDesignMovement::~FoldDesignMovement()
{
  int i,j;
  if(temparray)
  {
    delete[]temparray;
    temparray=NULL;
  }
  if(ttemparray)
  {
    delete[]ttemparray;
    ttemparray=NULL;
  }
  if(tttemparray)
  {
    delete[]tttemparray;
    tttemparray=NULL;
  }
  if(tmstr)
  {
    delete[]tmstr;
    tmstr=NULL;
  }
  if(tmstr2)
  {
    delete[]tmstr2;
    tmstr2=NULL;
  }
  if(alphaind)
  {
    delete[]alphaind;
    alphaind=NULL;
  }
  if(betaind)
  {
    delete[]betaind;
    betaind=NULL;
  }
  if(alphasind)
  {
    delete[]alphasind;
    alphasind=NULL;
  }
  if(bbptb)
  {
    delete[]bbptb;
    bbptb=NULL;
  }
  flagres=0;
}

FoldDesignMovement &FoldDesignMovement::operator=(
  const FoldDesignMovement &inputMovement
)
{
  if(this==&inputMovement)
  {
    return *this;
  }

  if(tmstr)
  {
    delete[]tmstr;
    tmstr=NULL;
  }
  if(tmstr2)
  {
    delete[]tmstr2;
    tmstr2=NULL;
  }
  if(fixed_decstr)
  {
    delete[]fixed_decstr;
    fixed_decstr=NULL;
  }
  if(alphaind)
  {
    delete[]alphaind;
    alphaind=NULL;
  }
  if(betaind)
  {
    delete[]betaind;
    betaind=NULL;
  }
  if(alphasind)
  {
    delete[]alphasind;
    alphasind=NULL;
  }
  /*if(phiprob)
  {
    delete[]phiprob;
    phiprob=NULL;
  }
  if(psiprob)
  {
    delete[]psiprob;
    psiprob=NULL;
  }
  */
  if(temparray)
  {
    delete[]temparray;
    temparray=NULL;
  }
  if(ttemparray)
  {
    delete[]ttemparray;
    ttemparray=NULL;
  }
  if(tttemparray)
  {
    delete[]tttemparray;
    tttemparray=NULL;
  }
  if(bbptb)
  {
    delete[]bbptb;
    bbptb=NULL;
  }

  T1=inputMovement.T1;
  T2=inputMovement.T2;
  S1=inputMovement.S1;
  S2=inputMovement.S2;
  n_rep=inputMovement.n_rep;
  inputInfo=inputMovement.inputInfo;
  geometry=inputMovement.geometry;
  energyFunction=inputMovement.energyFunction;

  int i;
  int seqLength=inputInfo.getSeqLength();
  int numSSE=inputInfo.getNumSSE();
  for(i=0;i<20;i++)
  {
    fragcont[i]=inputMovement.fragcont[i];
  }
  for(i=0;i<nosegdh;i++)
  {
     topdhlist[i]=inputMovement.topdhlist[i];
  }
  phiprob=inputMovement.phiprob;
  psiprob=inputMovement.psiprob;
  numtopdh=inputMovement.numtopdh;
  probeturn=inputMovement.probeturn;
  listeturn=inputMovement.listeturn;
  bborder=inputMovement.bborder;
  numsalpha=inputMovement.numsalpha;
  numalpha=inputMovement.numalpha;
  numbeta=inputMovement.numbeta;
  numeturn=inputMovement.numeturn;
  
  numbbptb=inputMovement.numbbptb;
  bbptb=new paircont[numbbptb];
  for(i=0;i<numbbptb;i++)
  { 
    bbptb[i].init=inputMovement.bbptb[i].init;
    bbptb[i].term=inputMovement.bbptb[i].term;
    bbptb[i].dist=inputMovement.bbptb[i].dist;
  }

  tmstr=new point3f[seqLength];
  tmstr2=new point3f[seqLength];
  fixed_decstr=new point3f[seqLength];
  alphaind=new int[numSSE];
  alphasind=new int[numSSE];
  betaind=new int[numSSE];
  for(i=0;i<numSSE;i++)
  {
    alphaind[i]=inputMovement.alphaind[i];
    betaind[i]=inputMovement.betaind[i];
    alphasind[i]=inputMovement.alphasind[i];
  }

  temparray=new double[n_rep];
  ttemparray=new double[n_rep];
  tttemparray=new double[n_rep];
  for(i=0;i<n_rep;i++)
  {
    temparray[i]=inputMovement.temparray[i];
    ttemparray[i]=inputMovement.ttemparray[i];
    tttemparray[i]=inputMovement.tttemparray[i];
  }

  return *this;
}

void FoldDesignMovement::configureSimulation(
  char *libDir,
  char *dataDir,
  const ParseFoldDesignMovementFiles &movementParameters,
  FoldDesignEnergyFunction inputEnergyFunction,
  InputData input,
  GeometryCalc inputGeometry,
  double tempLow,
  double tempHigh,
  double scale1,
  double scale2,
  int numRep

)
{
  int i;
  int seqLength;
  geometry=inputGeometry;
  inputInfo=input;
  energyFunction=inputEnergyFunction;

  phiprob=movementParameters.phiprob;
  psiprob=movementParameters.psiprob;
  bborder=movementParameters.bborder;
  probeturn=movementParameters.probeturn;
  listeturn=movementParameters.listeturn;
  numtopdh=movementParameters.numtopdh;
  numeturn=movementParameters.numeturn;
  for(i=0;i<nosegdh;i++)
  {
     topdhlist[i]=movementParameters.topdhlist[i];
  }
  for(i=0;i<20;i++)
  {
    fragcont[i]=inputInfo.getFragConformation(i);
  }
  seqLength=inputInfo.getSeqLength();
  T1=tempLow; 
  T2=tempHigh;
  S1=scale1;
  S2=scale2;
  n_rep=numRep;

  calcabind();
  calcbbptable(seqLength);
  setpara8(seqLength);
  settmstr(seqLength);

}

void FoldDesignMovement::setHbondRamp(
  int numcycles
)
{
  energyFunction.setHbondRamp(indCyc,numcycles);
}

int FoldDesignMovement::getMoveType()
{
  return i_move_type;
}

bool FoldDesignMovement::mcfragsweepCCD6(
  point3f *decstr,
  int numseq,
  int lps,
  int lpe,
  point3d pt1,
  int ind1,
  point3d pt2,
  int ind2,
  double oldenergy,
  double *newenergy
)
{
  int trandpos;
  bool flagphi,flagtor;
  int flagpt,indp[2];
  indp[0]=ind1;indp[1]=ind2;
  point3d tp[2]; 
  double cutdist=5.0;
  tp[0].x=pt1.x;tp[0].y=pt1.y;tp[0].z=pt1.z;
  tp[1].x=pt2.x;tp[1].y=pt2.y;tp[1].z=pt2.z;
  point3d pcur,p12,p13;
  int numiter=0;
  double tdist[2],ttheta,tphi,tpsi,tinner;
  memcpy(tmstr2,decstr,numseq*sizeof(point3f));
  flagpt=0;
  int threshiter=15*(lpe-lps+1);
  if(threshiter>100) threshiter=100;
  
  do
  {
    tdist[0]=10000;
    tdist[1]=10000;
    trandpos=lps+int((lpe-lps+1)*Random());
    if(Random()<0.5) flagphi=true; //Change phi n ca
    else flagphi=false; //Change psi ca c in pos+1
        
    flagpt=(flagpt+1)%2;
    if(flagphi)
    {
      p12=bf.setv(tmstr2[trandpos].x-tmstr2[trandpos].ptn.x,
                  tmstr2[trandpos].y-tmstr2[trandpos].ptn.y,
                  tmstr2[trandpos].z-tmstr2[trandpos].ptn.z);
      p13=bf.setv(tmstr2[indp[flagpt]].x-tmstr2[trandpos].ptn.x,
                  tmstr2[indp[flagpt]].y-tmstr2[trandpos].ptn.y,
                  tmstr2[indp[flagpt]].z-tmstr2[trandpos].ptn.z);
      tinner=bf.angv(p12,p13);
      //In one line no affect when rotate
      if(tinner<epsilon || PI-tinner<epsilon || bf.norm(p12)<epsilon) 
      {
        numiter++;
        continue;
      }
      pcur.x=tmstr2[indp[flagpt]].x;
      pcur.y=tmstr2[indp[flagpt]].y;
      pcur.z=tmstr2[indp[flagpt]].z;
      p13.x=tmstr2[trandpos].x;
      p13.y=tmstr2[trandpos].y;
      p13.z=tmstr2[trandpos].z;
      p12.x=tmstr2[trandpos].ptn.x;
      p12.y=tmstr2[trandpos].ptn.y;
      p12.z=tmstr2[trandpos].ptn.z;
      ttheta=bf.phi(pcur.x,pcur.y,pcur.z,
                    p12.x,p12.y,p12.z,
                    p13.x,p13.y,p13.z,
                    tp[flagpt].x,tp[flagpt].y,tp[flagpt].z);
      tphi=tmstr2[trandpos].tor[2]+ttheta;
      //if(tphi>360.0) tphi-=360.0;
      while(tphi>=360.0) tphi-=360.0;
      tmstr2[trandpos].tor[2]=tphi;   
    }
    else
    {
      p12=bf.setv(tmstr2[trandpos].ptc.x-tmstr2[trandpos].x,
                  tmstr2[trandpos].ptc.y-tmstr2[trandpos].y,
                  tmstr2[trandpos].ptc.z-tmstr2[trandpos].z);
      p13=bf.setv(tmstr2[indp[flagpt]].x-tmstr2[trandpos].x,
                  tmstr2[indp[flagpt]].y-tmstr2[trandpos].y,
                  tmstr2[indp[flagpt]].z-tmstr2[trandpos].z);
      tinner=bf.angv(p12,p13);

      //In one line no affect when rotate
      if(tinner<epsilon || PI-tinner<epsilon || bf.norm(p12)<epsilon) 
      {
        numiter++;
        continue;
      }
      pcur.x=tmstr2[indp[flagpt]].x;
      pcur.y=tmstr2[indp[flagpt]].y;
      pcur.z=tmstr2[indp[flagpt]].z;
      p12.x=tmstr2[trandpos].x;
      p12.y=tmstr2[trandpos].y;
      p12.z=tmstr2[trandpos].z;
      p13.x=tmstr2[trandpos].ptc.x;
      p13.y=tmstr2[trandpos].ptc.y;
      p13.z=tmstr2[trandpos].ptc.z;
      ttheta=bf.phi(pcur.x,pcur.y,pcur.z,
                    p12.x,p12.y,p12.z,
                    p13.x,p13.y,p13.z,
                    tp[flagpt].x,tp[flagpt].y,tp[flagpt].z);
      tpsi=tmstr2[trandpos+1].tor[0]+ttheta;
      //if(tpsi>360.0) tpsi-=360.0;
      while(tpsi>=360.0) tpsi-=360.0;
      tmstr2[trandpos+1].tor[0]=tpsi;
    }       
    numiter++;
    if(trandpos==0)
    {
      flagtor=geometry.tor2str(tmstr2,numseq,3);
    }
    else
    {
      flagtor=geometry.tor2strp(tmstr2,numseq,trandpos);
    }
    if(!flagtor && flagVerbose)
    {
      printf("tor2str wrong in CCD6 %d\n",trandpos);
    }
    pcur.x=tmstr2[ind1].x-tp[0].x;
    pcur.y=tmstr2[ind1].y-tp[0].y;
    pcur.z=tmstr2[ind1].z-tp[0].z;
    tdist[0]=pcur.x*pcur.x+pcur.y*pcur.y+pcur.z*pcur.z;
    pcur.x=tmstr2[ind2].x-tp[1].x;
    pcur.y=tmstr2[ind2].y-tp[1].y;
    pcur.z=tmstr2[ind2].z-tp[1].z;
    tdist[1]=pcur.x*pcur.x+pcur.y*pcur.y+pcur.z*pcur.z;
  } while(numiter<threshiter && (tdist[0]>cutdist || tdist[1]>cutdist));
  
  if(numiter<threshiter)
  {
    memcpy(decstr,tmstr2,numseq*sizeof(point3f));
    return true;
  }
  else
  {
    memcpy(decstr,tmstr2,numseq*sizeof(point3f));
    return false;
  }
}

/* Movement M5 substitutes one fragment in the decoy by another one randomly
 * selected from the position-specific fragment structures. It is one of the
 * important local movements in FoldDesign that can help reduce the conformational
 * search space and increase the quality of the local structures. To minimize
 * the conformational change of the surrounding residues and to increase the
 * acceptance rate, a Cyclic Coordinate Descent movement (CCD gap closure,
 * implemented by FoldDesignMovement::mcfragsweepCCD6) movement is followed, which
 * tries to adjust the segment conformation and makes it connect with the
 * anchor points of the surrounding chain. However, M5 still can be hardly
 * accepted when the fragment is long and the decoy structure becomes compact
 * after a number of cycles of simulations. Hence, we try to use more long
 * fragments at the commencement of the simulation where the decoy structure
 * has not been well packed, as implemented by:
 *    for(i=0;i<segcycle;i++) tprob[i]=subratio[indline][i];
 *    indsub=getmovetype(tprob,segcycle-1,rfrag)+1; */
bool FoldDesignMovement::mcfragsweepfragsub(
  point3f *decstr, 
  int numseq, 
  double tbeta,
  double oldenergy,
  double *newenergy
)
{
  int i,j,k;
  int attemptNum=0;
  int trandpos,trand2;
  bool flagclash=true;
  bool flagtor,swap_fixed;
  double threshclash=0.01;
  int numsub;
  int indsub;
  int startpos;
  double tprob[40];
  int indline=indCyc/10;
  int segcycle=21;
  if(indline>=segcycle) indline=segcycle-1;
  else if(indline<0) indline=0;
  for(i=0;i<segcycle;i++) tprob[i]=subratio[indline][i];

  do
  {
    flagclash=false;
    memcpy(tmstr,decstr,numseq*sizeof(point3f));
    int startinseg;
    int numfragsub;
    double rfrag=Random();
    indsub=getmovetype(tprob,segcycle-1,rfrag)+1;
    numsub=indsub;
    numfragsub=numsub;
    if(indsub<5)
    {
      indsub=7+int(6.0*Random());
      startinseg=(indsub-numsub)/2;
      numfragsub=indsub;
    }
    else startinseg=0;

    //Use fragment in pos trand2, starting seq from trandpos-startinseg, 
    //seglength=numfragsub, only substitute numsub 
    trandpos=startinseg+int((numseq-numfragsub-startinseg+1)*Random());

    swap_fixed=false;
    if(fixed_flag)
    {
      if(Random()<0.30)
      {
        for(i=0;i<start_fixed.size();i++)
        {
          if(trandpos>=start_fixed[i] && trandpos<=end_fixed[i])
          {
            numsub=int(Random()*(end_fixed[i]-trandpos))+1;
            if(numsub>=end_fixed[i]-trandpos+1) numsub--;
            swap_fixed=true;
            break;
          }
        }
      }
    }

    if(!flagclash) 
    {
      if(!swap_fixed)
      {
        trand2=int((topno)*Random());
        if(trandpos>=numseq-numfragsub+1) trandpos=numseq-numfragsub;
        if(trand2>=topno) trand2=topno-1;

        startpos=trand2*numfragsub+startinseg;
        for(k=0;k<numsub;k++)
        {
          tmstr[trandpos+k]=fragcont[indsub-1][startpos+k][trandpos-startinseg];
          tmstr[trandpos+k].ss3=decstr[trandpos+k].ss3;
          tmstr[trandpos+k].ss8=decstr[trandpos+k].ss8;
          tmstr[trandpos+k].ssm=decstr[trandpos+k].ssm;
          tmstr[trandpos+k].aaa=decstr[trandpos+k].aaa;
          tmstr[trandpos+k].iaa=decstr[trandpos+k].iaa;
        }
      }
      else
      {
        for(k=0;k<numsub;k++)
        {
          tmstr[trandpos+k]=fixed_decstr[trandpos+k];
          tmstr[trandpos+k].ss3=decstr[trandpos+k].ss3;
          tmstr[trandpos+k].ss8=decstr[trandpos+k].ss8;
          tmstr[trandpos+k].ssm=decstr[trandpos+k].ssm;
          tmstr[trandpos+k].aaa=decstr[trandpos+k].aaa;
          tmstr[trandpos+k].iaa=decstr[trandpos+k].iaa;
        }
      }

      //////////////////don't change
      if(trandpos==0) flagtor=geometry.tor2str(tmstr,numseq,3);
      else flagtor=geometry.tor2strp(tmstr,numseq,trandpos);
      if(!flagtor && flagVerbose) printf("tor2str wrong in sub %d\n",trandpos);
      flagclash=false;
    }
    else if(attemptNum>=MAX_ITERATION)
    {
      summcsub[2]++;
      return false;
    }
  } while(flagclash);

  if(tbeta<=0 && (trandpos+numsub>=numseq || trandpos<=0 || tbeta>-1.5))
  {   // single anchored
    summcsub[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }
  else if (tbeta>0) *newenergy=energyFunction.calcTotalEnergy(tmstr,numseq);

  if(tbeta<=0 || (*newenergy>oldenergy && trandpos+numsub<numseq && trandpos>0))
  {
    int ips,ipe;
    ipe=trandpos+numsub-1;
    if(ipe>=numseq-1) ipe=numseq-2;
    while(ipe>=0 && tmstr[ipe].ss3!='C') ipe--;
    ips=ipe;
    while(ips>=0 && tmstr[ips].ss3=='C') ips--;
    if(ipe>=0) ips++;
    int ips2,ipe2;
    ips2=trandpos+numsub;
    while(ips2<numseq && tmstr[ips2].ss3!='C') ips2++;
    ipe2=ips2;
    while(ipe2<numseq && tmstr[ipe2].ss3=='C') ipe2++;
    if(ips2<numseq) ipe2--;
    double toldener=oldenergy;
    double tnewener;
    point3d ptran;
    ptran.x=decstr[numseq-1].x-tmstr[numseq-1].x;
    ptran.y=decstr[numseq-1].y-tmstr[numseq-1].y;
    ptran.z=decstr[numseq-1].z-tmstr[numseq-1].z;
    if(ipe>=0 && ipe<numseq-1 && trandpos+numsub-ipe<ips2-trandpos+numsub)
    {
      if(ipe-ips>10) ips=ipe-10;
      point3d pt1,pt2;
      pt1=bf.setv(decstr[numseq-2].x,decstr[numseq-2].y,decstr[numseq-2].z);
      pt2=bf.setv(decstr[numseq-1].x,decstr[numseq-1].y,decstr[numseq-1].z);
      mcfragsweepCCD6(tmstr,numseq,ips,ipe,pt1,numseq-2,pt2,numseq-1,toldener,&tnewener);
    }
    else if(ips2<numseq-1 && ipe2<numseq-1)
    {
      if(ipe2-ips2>10) ipe2=ips2+10;
      point3d pt1,pt2;
      pt1=bf.setv(decstr[numseq-2].x,decstr[numseq-2].y,decstr[numseq-2].z);
      pt2=bf.setv(decstr[numseq-1].x,decstr[numseq-1].y,decstr[numseq-1].z);
      mcfragsweepCCD6(tmstr,numseq,ips2,ipe2,pt1,numseq-2,pt2,numseq-1,toldener,&tnewener);
    }
    else 
    {
      point3d p12,p13,pcur,pnew;
      int ppos=trandpos+numsub;
      if(ppos==numseq-1) ppos--;
      pnew.x=decstr[numseq-1].x;
      pnew.y=decstr[numseq-1].y;
      pnew.z=decstr[numseq-1].z;
      pcur.x=tmstr[numseq-1].x;pcur.y=tmstr[numseq-1].y;pcur.z=tmstr[numseq-1].z;
      p13.x=tmstr[ppos].x;p13.y=tmstr[ppos].y;p13.z=tmstr[ppos].z;
      p12.x=tmstr[ppos].ptn.x;p12.y=tmstr[ppos].ptn.y;p12.z=tmstr[ppos].ptn.z;
      double ttheta=bf.phi(pcur.x,pcur.y,pcur.z,
                           p12.x,p12.y,p12.z,
                           p13.x,p13.y,p13.z,
                           pnew.x,pnew.y,pnew.z);
      tmstr[trandpos].tor[2]+=ttheta;
      while(tmstr[trandpos].tor[2]>=360.0) tmstr[trandpos].tor[2]-=360.0;
      if(trandpos-1==0) flagtor=geometry.tor2str(tmstr,numseq,3);
      else flagtor=geometry.tor2strp(tmstr,numseq,trandpos-1);
      if(!flagtor && flagVerbose) printf("tor2str wrong in suba\n");
      //////
      p12.x=tmstr[ppos].x;p12.y=tmstr[ppos].y;p12.z=tmstr[ppos].z;
      p13.x=tmstr[ppos].ptc.x;p13.y=tmstr[ppos].ptc.y;p13.z=tmstr[ppos].ptc.z;
      ttheta=bf.phi(pcur.x,pcur.y,pcur.z,
                    p12.x,p12.y,p12.z,
                    p13.x,p13.y,p13.z,
                    pnew.x,pnew.y,pnew.z);
      tmstr[ppos+1].tor[0]+=ttheta;
      while(tmstr[ppos+1].tor[0]>=360.0) tmstr[ppos+1].tor[0]-=360.0;
      if(ppos==0) flagtor=geometry.tor2str(tmstr,numseq,3);
      else flagtor=geometry.tor2strp(tmstr,numseq,ppos);
      if(!flagtor && flagVerbose) printf("tor2str wrong in subb\n");
      //////
      ppos=trandpos-1;
      if(ppos==0) ppos++;
      pnew.x=decstr[numseq-1].x;
      pnew.y=decstr[numseq-1].y;
      pnew.z=decstr[numseq-1].z;
      pcur.x=tmstr[numseq-1].x;pcur.y=tmstr[numseq-1].y;pcur.z=tmstr[numseq-1].z;
      //////
      p13.x=tmstr[ppos].x;p13.y=tmstr[ppos].y;p13.z=tmstr[ppos].z;
      p12.x=tmstr[ppos].ptn.x;p12.y=tmstr[ppos].ptn.y;p12.z=tmstr[ppos].ptn.z;
      ttheta=bf.phi(pcur.x,pcur.y,pcur.z,
                    p12.x,p12.y,p12.z,
                    p13.x,p13.y,p13.z,
                    pnew.x,pnew.y,pnew.z);
      tmstr[trandpos].tor[2]+=ttheta;
      while(tmstr[trandpos].tor[2]>=360.0) tmstr[trandpos].tor[2]-=360.0;
      if(trandpos-1==0) flagtor=geometry.tor2str(tmstr,numseq,3);
      else flagtor=geometry.tor2strp(tmstr,numseq,trandpos-1);
      if(!flagtor && flagVerbose) printf("tor2str wrong in subc\n");
      //////
      p12.x=tmstr[ppos].x;p12.y=tmstr[ppos].y;p12.z=tmstr[ppos].z;
      p13.x=tmstr[ppos].ptc.x;p13.y=tmstr[ppos].ptc.y;p13.z=tmstr[ppos].ptc.z;
      ttheta=bf.phi(pcur.x,pcur.y,pcur.z,
                    p12.x,p12.y,p12.z,
                    p13.x,p13.y,p13.z,
                    pnew.x,pnew.y,pnew.z);
      tmstr[ppos+1].tor[0]+=ttheta;
      while(tmstr[ppos+1].tor[0]>=360.0) tmstr[ppos+1].tor[0]-=360.0;
      if(ppos==0) flagtor=geometry.tor2str(tmstr,numseq,3);
      else flagtor=geometry.tor2strp(tmstr,numseq,ppos);
      if(!flagtor && flagVerbose) printf("tor2str wrong in subd\n");
    }

    if(tbeta<=0) //Double anchored
    {
      summcsub[0]++;
      memcpy(decstr,tmstr,numseq*sizeof(point3f));
      return true;
    }
    *newenergy=energyFunction.calcTotalEnergy(tmstr,numseq);
  }

  //////////Make decision
  if(*newenergy<oldenergy)
  {
    summcsub[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }
  else
  {
    double trand;
    trand=Random();
    if(trand<exp(-tbeta*(*newenergy-oldenergy)))
    {
      summcsub[1]++;
      memcpy(decstr,tmstr,numseq*sizeof(point3f));
      return true;
    }
    else
    {
      summcsub[2]++;
     return false;
    }
  }
}

/* M3: change phi (C[-1]-N-CA-C) backbone torsion angle. */
bool FoldDesignMovement::mcfragsweepphi(
  point3f *decstr, 
  int numseq, 
  double tbeta,
  double oldenergy, 
  double *newenergy
)
{
  int i;
  int attemptNum=0; 
  int trandpos;  //Random residue
  double trand2; //Random angle
  bool flagclash=true;
  bool flagtor;
  double threshclash=0.01;
  int indp;
  int indss; //Secondary structure index
  int psi;

  do
  {
    flagclash=false;
    memcpy(tmstr,decstr,numseq*sizeof(point3f));
    trandpos=int((numseq)*Random());
    if(tmstr[trandpos].ssm!='C' && Random()>0.2)
    {
      summcphi[2]++;
      return false;
    }
    trand2=Random();

    if(!flagclash) 
    {
      if(tmstr[trandpos].ss8=='H') 
      {
        indp=bf.findpos2(phiprob[0],0,359,trand2);
      }
      else if(tmstr[trandpos].ss8=='E') 
      {
        indp=bf.findpos2(phiprob[1],0,359,trand2);
      }
      else if(tmstr[trandpos].ss8=='C') 
      {
        indp=bf.findpos2(phiprob[2],0,359,trand2);
      }
      else if(tmstr[trandpos].ss8=='h') 
      {
        indp=bf.findpos2(phiprob[3],0,359,trand2);
      }
      else if(tmstr[trandpos].ss8=='e') 
      {
        indp=bf.findpos2(phiprob[4],0,359,trand2);
      }
      else if(tmstr[trandpos].ss8=='c') 
      {
        indp=bf.findpos2(phiprob[5],0,359,trand2);
      }
      else if(tmstr[trandpos].ss8=='G') 
      {
        indp=bf.findpos2(phiprob[6],0,359,trand2);
      }
      else if(tmstr[trandpos].ss8=='I') 
      {
        indp=bf.findpos2(phiprob[7],0,359,trand2);
      }
      else if(tmstr[trandpos].ss8=='B') 
      {
        indp=bf.findpos2(phiprob[8],0,359,trand2);
      }
      else if(tmstr[trandpos].ss8=='S') 
      {
        indp=bf.findpos2(phiprob[9],0,359,trand2);
      }
      else if(tmstr[trandpos].ss8=='T') 
      {
        indp=bf.findpos2(phiprob[10],0,359,trand2);
      }
	
      if(indp>360 || indp<-360) 
      {
        summcphi[2]++;
        return false;
      }

      tmstr[trandpos].tor[2]=indp+0.5;
      //Correction for periodicity of torsion angles
      if(tmstr[trandpos].tor[2]<0) tmstr[trandpos].tor[2]+=360;
      else if(tmstr[trandpos].tor[2]>=360) tmstr[trandpos].tor[2]-=360;
      //Make conformation change
      if(trandpos==0) flagtor=geometry.tor2str(tmstr,numseq,3);
      else flagtor=geometry.tor2strp(tmstr,numseq,trandpos);
      if(!flagtor && flagVerbose) printf("tor2str wrong in phi %d\n",trandpos);
      flagclash=false;
    }        
    else if(attemptNum>=MAX_ITERATION)
    {
      summcphi[2]++;
      return false;
    }
  } while(flagclash);

  if(tbeta<0) //Skip energy calculation
  {
    summcphi[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f)); //Should be conformation copying
    return true;
  }
  *newenergy=energyFunction.calcTotalEnergy(tmstr,numseq); 
  
  //printf("-------- phi--------> beta=%8.3f, temp=%8.3f\n",tbeta,1/tbeta);
  if(*newenergy<oldenergy)
  {
    summcphi[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f)); //Should be conformation copying
    return true;
  }
  else
  {
    double trand;
    trand=Random();
    if(trand<exp(-tbeta*(*newenergy-oldenergy))) //Metropolis criterion
    {       
      summcphi[1]++;
      memcpy(decstr,tmstr,numseq*sizeof(point3f));
      return true;
    }
    else
    {
      summcphi[2]++;
      return false;
    }
  }
}

/* M3: change psi (CA-C-N-CA[+1]) backbone torsion angle. */
bool FoldDesignMovement::mcfragsweeppsi(
  point3f *decstr,
  int numseq,
  double tbeta,
  double oldenergy,
  double *newenergy
)
{
  int i;
  int attemptNum=0;
  int trandpos;
  double trand2;
  bool flagclash=true;
  bool flagtor;
  double threshclash=0.01;
  int indp;
  int phi;

  do
  {
    flagclash=false;
    memcpy(tmstr,decstr,numseq*sizeof(point3f));
    trandpos=int((numseq)*Random());
    if(tmstr[trandpos].ssm!='C' && Random()>0.2)
    {
      summcpsi[2]++;
      return false;
    }

    if(trandpos==0)
    {
      summcpsi[2]++;
      return false;
    }
    trand2=Random();
    phi=tmstr[trandpos].tor[2];
   
    if(!flagclash) 
    {
      if(tmstr[trandpos].ss8=='H') 
      {
        indp=bf.findpos2(psiprob[0],0,359,trand2);
      }
      else if(tmstr[trandpos].ss8=='E') 
      {
        indp=bf.findpos2(psiprob[1],0,359,trand2);
      }
      else if(tmstr[trandpos].ss8=='C') 
      {
        indp=bf.findpos2(psiprob[2],0,359,trand2);
      }
      else if(tmstr[trandpos].ss8=='h') 
      {
        indp=bf.findpos2(psiprob[3],0,359,trand2);
      }
      else if(tmstr[trandpos].ss8=='e') 
      {
        indp=bf.findpos2(psiprob[4],0,359,trand2);
      }
      else if(tmstr[trandpos].ss8=='c') 
      {
        indp=bf.findpos2(psiprob[5],0,359,trand2);
      }
      else if(tmstr[trandpos].ss8=='G') 
      {
        indp=bf.findpos2(psiprob[6],0,359,trand2);
      }
      else if(tmstr[trandpos].ss8=='I') 
      {
        indp=bf.findpos2(psiprob[7],0,359,trand2);
      }
      else if(tmstr[trandpos].ss8=='B') 
      {
        indp=bf.findpos2(psiprob[8],0,359,trand2);
      }
      else if(tmstr[trandpos].ss8=='S') 
      {
        indp=bf.findpos2(psiprob[9],0,359,trand2);
      }
      else if(tmstr[trandpos].ss8=='T') 
      {
        indp=bf.findpos2(psiprob[10],0,359,trand2);
      }

      if(indp>360 || indp<-360) 
      {
        summcpsi[0]++;
        return false;
      }

      tmstr[trandpos].tor[0]=indp+0.5;
      if(tmstr[trandpos].tor[0]<0) tmstr[trandpos].tor[0]+=360;
      else if(tmstr[trandpos].tor[0]>=360) tmstr[trandpos].tor[0]-=360;
      if(trandpos==0) flagtor=geometry.tor2str(tmstr,numseq,3);
      else flagtor=geometry.tor2strp(tmstr,numseq,trandpos);
      if(!flagtor) printf("tor2str wrong in psi %d\n",trandpos);
      flagclash=false;
    }
    else if(attemptNum>=MAX_ITERATION)
    {
      summcpsi[2]++;
      return false;
    }
  }while(flagclash);

  if(tbeta<0) // skip energy calculation
  {
    summcpsi[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }
  *newenergy=energyFunction.calcTotalEnergy(tmstr,numseq);
  if(*newenergy<oldenergy)
  {
    summcpsi[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }
  else
  {
    double trand;
    trand=Random();
    if(trand<exp(-tbeta*(*newenergy-oldenergy)))
    {
      summcpsi[1]++;
      memcpy(decstr,tmstr,numseq*sizeof(point3f));
      return true;
    }
    else
    {
      summcpsi[2]++;
      return false;
    }
  }
}

/* M4: substitute residue geometry by that in topdh, which was generated
 * from 10 residue fragments */
bool FoldDesignMovement::mcfragsweeptopdh(
  point3f *decstr,
  int numseq,
  double tbeta,
  double oldenergy,
  double *newenergy
)
{
  int i,j,k;
  int attemptNum=0;
  int tn;
  float tf,tr;
  int trandpos;
  char tmpname[6];tmpname[5]='\0';
  bool flagclash=true;
  bool flagtor;
  double threshclash=0.01;
  do
  {
    flagclash=false;
    memcpy(tmstr,decstr,numseq*sizeof(point3f));
    trandpos=int((numseq)*Random());
    tn=numtopdh[trandpos];
    if(!flagclash) 
    {
      tr=Random();
      for(j=0;j<tn;j++)
      {
        tf=topdhlist[j][trandpos].prob;
        if(tf>=tr) break;
      }
      if(strcmp(tmstr[trandpos].name,topdhlist[j][trandpos].name)==0 && 
         tmstr[trandpos].resind==topdhlist[j][trandpos].resind) continue;
      memcpy(&tmstr[trandpos],&topdhlist[j][trandpos],sizeof(point3f));
      if(trandpos==0) flagtor=geometry.tor2str(tmstr,numseq,3);
      else flagtor=geometry.tor2strp(tmstr,numseq,trandpos);

      if(!flagtor && flagVerbose) printf("tor2str wrong in topdh %d\n",trandpos);
      flagclash=false;
    }
    else if(attemptNum>=MAX_ITERATION)
    {
      summctop[2]++;
      return false;
    }
  } while(flagclash);

  if (tbeta<0) //Skip energy calculation
  {
    summctop[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }

  *newenergy=energyFunction.calcTotalEnergy(tmstr,numseq);
  if(*newenergy<oldenergy)
  {
    summctop[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }
  else
  {
    double trand;
    trand=Random();
    if(trand<exp(-tbeta*(*newenergy-oldenergy)))
    {
      summctop[1]++;
      memcpy(decstr,tmstr,numseq*sizeof(point3f));
      return true;
    }
    else
    {
      summctop[2]++;
      return false;
    }
  }
}

/* Movement M7 rotates the atoms of a randomly selected segment around
 * the axis connecting the two ending CA atoms, similar to loop modeling
 * movement in DEMO. */
bool FoldDesignMovement::mcfragsweepmid(
  point3f *decstr,
  int numseq,
  double tbeta,
  double oldenergy,
  double *newenergy
)
{ 
  int i,j;
  int attemptNum=0;
  int istart,iend;
  int trandpos,trandpos2;
  double trandang;
  double mati[9];
  point3d  v1;
  bool flagclash=true;
  bool flagtor;
  double threshclash=0.01;
    
  do
  {
    flagclash=false;
    memcpy(tmstr,decstr,numseq*sizeof(point3f));
    do
    {
      trandpos2=int(12*Random())-1;
      trandpos=1+int((numseq-trandpos2-2)*Random());
      trandpos2+=trandpos;
    } while(trandpos<1 || trandpos2>numseq-2);

    if(!flagclash) 
    {
      trandang=180.0*Random()-90.0;
      trandang=trandang*raddeg; //rotangle
      v1.x=tmstr[trandpos2+1].x-tmstr[trandpos-1].x;
      v1.y=tmstr[trandpos2+1].y-tmstr[trandpos-1].y;
      v1.z=tmstr[trandpos2+1].z-tmstr[trandpos-1].z;
      v1=bf.unit(v1); //direction
    
      //method2
      double vcos=cos(trandang);
      double vsin=sin(trandang);
      double vxx=v1.x*v1.x;
      double vxy=v1.x*v1.y;
      double vxz=v1.x*v1.z;
      double vyy=v1.y*v1.y;
      double vyz=v1.y*v1.z;
      double vzz=v1.z*v1.z;
      mati[0]=vxx+(1-vxx)*vcos;
      mati[1]=vxy*(1-vcos)-v1.z*vsin;
      mati[2]=vxz*(1-vcos)+v1.y*vsin;
      mati[3]=vxy*(1-vcos)+v1.z*vsin;
      mati[4]=vyy+(1-vyy)*vcos;
      mati[5]=vyz*(1-vcos)-v1.x*vsin;
      mati[6]=vxz*(1-vcos)-v1.y*vsin;
      mati[7]=vyz*(1-vcos)+v1.x*vsin;
      mati[8]=vzz+(1-vzz)*vcos;

      for(i=trandpos;i<=trandpos2;i++)
      {
        v1.x=(tmstr[i].x-tmstr[trandpos-1].x)*mati[0]
            +(tmstr[i].y-tmstr[trandpos-1].y)*mati[1]
            +(tmstr[i].z-tmstr[trandpos-1].z)*mati[2]+tmstr[trandpos-1].x;
        v1.y=(tmstr[i].x-tmstr[trandpos-1].x)*mati[3]
            +(tmstr[i].y-tmstr[trandpos-1].y)*mati[4]
            +(tmstr[i].z-tmstr[trandpos-1].z)*mati[5]+tmstr[trandpos-1].y;
        v1.z=(tmstr[i].x-tmstr[trandpos-1].x)*mati[6]
            +(tmstr[i].y-tmstr[trandpos-1].y)*mati[7]
            +(tmstr[i].z-tmstr[trandpos-1].z)*mati[8]+tmstr[trandpos-1].z;
        tmstr[i].x=v1.x;
        tmstr[i].y=v1.y;
        tmstr[i].z=v1.z;
        v1.x=(tmstr[i].ptn.x-tmstr[trandpos-1].x)*mati[0]
            +(tmstr[i].ptn.y-tmstr[trandpos-1].y)*mati[1]
            +(tmstr[i].ptn.z-tmstr[trandpos-1].z)*mati[2]+tmstr[trandpos-1].x;
        v1.y=(tmstr[i].ptn.x-tmstr[trandpos-1].x)*mati[3]
            +(tmstr[i].ptn.y-tmstr[trandpos-1].y)*mati[4]
            +(tmstr[i].ptn.z-tmstr[trandpos-1].z)*mati[5]+tmstr[trandpos-1].y;
        v1.z=(tmstr[i].ptn.x-tmstr[trandpos-1].x)*mati[6]
            +(tmstr[i].ptn.y-tmstr[trandpos-1].y)*mati[7]
            +(tmstr[i].ptn.z-tmstr[trandpos-1].z)*mati[8]+tmstr[trandpos-1].z;
        tmstr[i].ptn.x=v1.x;
        tmstr[i].ptn.y=v1.y;
        tmstr[i].ptn.z=v1.z;
        v1.x=(tmstr[i].ptc.x-tmstr[trandpos-1].x)*mati[0]
            +(tmstr[i].ptc.y-tmstr[trandpos-1].y)*mati[1]
            +(tmstr[i].ptc.z-tmstr[trandpos-1].z)*mati[2]+tmstr[trandpos-1].x;
        v1.y=(tmstr[i].ptc.x-tmstr[trandpos-1].x)*mati[3]
            +(tmstr[i].ptc.y-tmstr[trandpos-1].y)*mati[4]
            +(tmstr[i].ptc.z-tmstr[trandpos-1].z)*mati[5]+tmstr[trandpos-1].y;
        v1.z=(tmstr[i].ptc.x-tmstr[trandpos-1].x)*mati[6]
            +(tmstr[i].ptc.y-tmstr[trandpos-1].y)*mati[7]
            +(tmstr[i].ptc.z-tmstr[trandpos-1].z)*mati[8]+tmstr[trandpos-1].z;
        tmstr[i].ptc.x=v1.x;
        tmstr[i].ptc.y=v1.y;
        tmstr[i].ptc.z=v1.z;
      }
      //ptc ptn
      i=trandpos-1;
      v1.x=(tmstr[i].ptc.x-tmstr[trandpos-1].x)*mati[0]
          +(tmstr[i].ptc.y-tmstr[trandpos-1].y)*mati[1]
          +(tmstr[i].ptc.z-tmstr[trandpos-1].z)*mati[2]+tmstr[trandpos-1].x;
      v1.y=(tmstr[i].ptc.x-tmstr[trandpos-1].x)*mati[3]
          +(tmstr[i].ptc.y-tmstr[trandpos-1].y)*mati[4]
          +(tmstr[i].ptc.z-tmstr[trandpos-1].z)*mati[5]+tmstr[trandpos-1].y;
      v1.z=(tmstr[i].ptc.x-tmstr[trandpos-1].x)*mati[6]
          +(tmstr[i].ptc.y-tmstr[trandpos-1].y)*mati[7]
          +(tmstr[i].ptc.z-tmstr[trandpos-1].z)*mati[8]+tmstr[trandpos-1].z;
      tmstr[i].ptc.x=v1.x;
      tmstr[i].ptc.y=v1.y;
      tmstr[i].ptc.z=v1.z;
      i=trandpos2+1;
      v1.x=(tmstr[i].ptn.x-tmstr[trandpos-1].x)*mati[0]
          +(tmstr[i].ptn.y-tmstr[trandpos-1].y)*mati[1]
          +(tmstr[i].ptn.z-tmstr[trandpos-1].z)*mati[2]+tmstr[trandpos-1].x;
      v1.y=(tmstr[i].ptn.x-tmstr[trandpos-1].x)*mati[3]
          +(tmstr[i].ptn.y-tmstr[trandpos-1].y)*mati[4]
          +(tmstr[i].ptn.z-tmstr[trandpos-1].z)*mati[5]+tmstr[trandpos-1].y;
      v1.z=(tmstr[i].ptn.x-tmstr[trandpos-1].x)*mati[6]
          +(tmstr[i].ptn.y-tmstr[trandpos-1].y)*mati[7]
          +(tmstr[i].ptn.z-tmstr[trandpos-1].z)*mati[8]+tmstr[trandpos-1].z;
      tmstr[i].ptn.x=v1.x;
      tmstr[i].ptn.y=v1.y;
      tmstr[i].ptn.z=v1.z;
      istart=trandpos-1;
      iend=trandpos2+2;
      if(iend>=numseq) iend=numseq-1;
      geometry.str2torp(tmstr,numseq,istart,iend);
      if(fabs(tmstr[trandpos-1].ang[2]-angncac)>8*2.818 || 
         fabs(tmstr[trandpos2+1].ang[2]-angncac)>8*2.818)
      {
        summcmid[2]++;
        return false;
      }
      if(istart==0 || istart==1) flagtor=geometry.tor2str(tmstr,numseq,3);
      else flagtor=geometry.tor2strp(tmstr,numseq,istart);
      if(!flagtor && flagVerbose) printf("tor2str wrong in mid %d\n",istart);
      flagclash=false;
    }
    else if(attemptNum>=MAX_ITERATION)
    {
      summcmid[2]++;
      return false;
    }
  } while(flagclash);

  if(tbeta<0) //Skip energy calculation
  {
    summcmid[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }

  *newenergy=energyFunction.calcTotalEnergy(tmstr,numseq);
  if(flagcaca && enelist[14]<enelistbk[14] && 
     (*newenergy<oldenergy || (*newenergy-oldenergy)<-0.01*oldenergy))
  {
    summcmid[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }
  if(*newenergy<oldenergy)
  {
    summcmid[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }
  else
  {
    double trand;
    trand=Random();
    if(trand<exp(-tbeta*(*newenergy-oldenergy)))
    {
      summcmid[1]++;
      memcpy(decstr,tmstr,numseq*sizeof(point3f));
      return true;
    }
    else
    {
      summcmid[2]++;
      return false;
    }
  }
}

/* Move residue k-1 respective to residue k in decoy structure 'tmpstr'
 * mtype - the atom for movement
 *         0: move C[-1]  towards/from N
 *         1: move C[-1]  towards/from CA
 *         2: move CA[-1] towards/from N
 *         3: move CA[-1] towards/from CA
 *         4: move CA[-1] towards/from C[-1]
 *         5: move N[-1]  towards/from CA[-1]
 *         6: move N[-1]  towards/from C[-1]  */
void FoldDesignMovement::singlemoveLMPb(
  point3f *tmstr,
  int k,
  int mtype
)
{
  point3d rp;
  double rnorm;
  switch(mtype)
  {
    case 0: //1 c-1 n
      rp=bf.setv(tmstr[k-1].ptc.x-tmstr[k].ptn.x,
                 tmstr[k-1].ptc.y-tmstr[k].ptn.y,
                 tmstr[k-1].ptc.z-tmstr[k].ptn.z);
      rnorm=bf.norm(rp);
      if(rnorm<epsilon) //Two in the same place
      {
        tmstr[k-1].ptc.x-=0.0001f;
        tmstr[k-1].ptc.y-=0.0001f;
        tmstr[k-1].ptc.z-=0.0001f;
        rp=bf.setv(tmstr[k-1].ptc.x-tmstr[k].ptn.x,
                   tmstr[k-1].ptc.y-tmstr[k].ptn.y,
                   tmstr[k-1].ptc.z-tmstr[k].ptn.z);
        rnorm=bf.norm(rp);  
      }
      if(rnorm<lencn-delcn || rnorm>lencn+delcn)
      {
        rp=bf.scal(rp,lencn/rnorm);
        tmstr[k-1].ptc.x=tmstr[k].ptn.x+rp.x;
        tmstr[k-1].ptc.y=tmstr[k].ptn.y+rp.y;
        tmstr[k-1].ptc.z=tmstr[k].ptn.z+rp.z;
      }
      break;
    case 1: //4 c-1 ca
      rp=bf.setv(tmstr[k-1].ptc.x-tmstr[k].x,
                 tmstr[k-1].ptc.y-tmstr[k].y,
                 tmstr[k-1].ptc.z-tmstr[k].z);
      rnorm=bf.norm(rp);
      if(rnorm<epsilon) //Two in the same place
      {
        tmstr[k-1].ptc.x-=0.0001f;
        tmstr[k-1].ptc.y-=0.0001f;
        tmstr[k-1].ptc.z-=0.0001f;
        rp=bf.setv(tmstr[k-1].ptc.x-tmstr[k].x,
                   tmstr[k-1].ptc.y-tmstr[k].y,
                   tmstr[k-1].ptc.z-tmstr[k].z);
        rnorm=bf.norm(rp);  
      }
      if(rnorm<lencca1-delcca1 || rnorm>lencca1+delcca1)
      {
        rp=bf.scal(rp,lencca1/rnorm);
        tmstr[k-1].ptc.x=tmstr[k].x+rp.x;
        tmstr[k-1].ptc.y=tmstr[k].y+rp.y;
        tmstr[k-1].ptc.z=tmstr[k].z+rp.z;
      }
      break;
    case 2: //2 ca-1 n
      rp=bf.setv(tmstr[k-1].x-tmstr[k].ptn.x,
                 tmstr[k-1].y-tmstr[k].ptn.y,
                 tmstr[k-1].z-tmstr[k].ptn.z);
      rnorm=bf.norm(rp);
      if(rnorm<epsilon) //Two in the same place
      {
        tmstr[k-1].x-=0.0001f;
        tmstr[k-1].y-=0.0001f;
        tmstr[k-1].z-=0.0001f;
        rp=bf.setv(tmstr[k-1].x-tmstr[k].ptn.x,
                   tmstr[k-1].y-tmstr[k].ptn.y,
                   tmstr[k-1].z-tmstr[k].ptn.z);
        rnorm=bf.norm(rp);  
      }
      if(rnorm<lencan1-delcan1 || rnorm>lencan1+delcan1)
      {
        rp=bf.scal(rp,lencan1/rnorm);
        tmstr[k-1].x=tmstr[k].ptn.x+rp.x;
        tmstr[k-1].y=tmstr[k].ptn.y+rp.y;
        tmstr[k-1].z=tmstr[k].ptn.z+rp.z;
      }
      break;
    case 3: //5 ca-1 ca
      rp=bf.setv(tmstr[k-1].x-tmstr[k].x,
                 tmstr[k-1].y-tmstr[k].y,
                 tmstr[k-1].z-tmstr[k].z);
      rnorm=bf.norm(rp);
      if(rnorm<epsilon) //Two in the same place
      {
        tmstr[k-1].x-=0.0001f;
        tmstr[k-1].y-=0.0001f;
        tmstr[k-1].z-=0.0001f;
        rp=bf.setv(tmstr[k-1].x-tmstr[k].x,
                   tmstr[k-1].y-tmstr[k].y,
                   tmstr[k-1].z-tmstr[k].z);
        rnorm=bf.norm(rp);  
      }
      if(rnorm<lencaca-delcaca || rnorm>lencaca+delcaca)
      {
        rp=bf.scal(rp,lencaca/rnorm);
        tmstr[k-1].x=tmstr[k].x+rp.x;
        tmstr[k-1].y=tmstr[k].y+rp.y;
        tmstr[k-1].z=tmstr[k].z+rp.z;
      }
      break;
    case 4: //6 ca c
      rp=bf.setv(tmstr[k-1].x-tmstr[k-1].ptc.x,
                 tmstr[k-1].y-tmstr[k-1].ptc.y,
                 tmstr[k-1].z-tmstr[k-1].ptc.z);
      rnorm=bf.norm(rp);
      if(rnorm<epsilon) //Two in the same place
      {
        tmstr[k-1].x-=0.0001f;
        tmstr[k-1].y-=0.0001f;tmstr[k-1].z-=0.0001f;
        rp=bf.setv(tmstr[k-1].x-tmstr[k-1].ptc.x,
                   tmstr[k-1].y-tmstr[k-1].ptc.y,
                   tmstr[k-1].z-tmstr[k-1].ptc.z);
        rnorm=bf.norm(rp);  
      }
      if(rnorm<lencac-delcac || rnorm>lencac+delcac)
      {
        rp=bf.scal(rp,lencac/rnorm);
        tmstr[k-1].x=tmstr[k-1].ptc.x+rp.x;
        tmstr[k-1].y=tmstr[k-1].ptc.y+rp.y;
        tmstr[k-1].z=tmstr[k-1].ptc.z+rp.z;
      }
      break;
    case 5:  //3 n ca
      rp=bf.setv(tmstr[k-1].ptn.x-tmstr[k-1].x,
                 tmstr[k-1].ptn.y-tmstr[k-1].y,
                 tmstr[k-1].ptn.z-tmstr[k-1].z);
      rnorm=bf.norm(rp);
      if(rnorm<epsilon) //Two in the same place
      {
        tmstr[k-1].ptn.x-=0.0001f;
        tmstr[k-1].ptn.y-=0.0001f;
        tmstr[k-1].ptn.z-=0.0001f;
        rp=bf.setv(tmstr[k-1].ptn.x-tmstr[k-1].x,
                   tmstr[k-1].ptn.y-tmstr[k-1].y,
                   tmstr[k-1].ptn.z-tmstr[k-1].z);
        rnorm=bf.norm(rp);  
      }
      if(rnorm<lennca-delnca || rnorm>lennca+delnca)
      {
        rp=bf.scal(rp,lennca/rnorm);
        tmstr[k-1].ptn.x=tmstr[k-1].x+rp.x;
        tmstr[k-1].ptn.y=tmstr[k-1].y+rp.y;
        tmstr[k-1].ptn.z=tmstr[k-1].z+rp.z;
      }
      break;
    case 6: //7 n c
      rp=bf.setv(tmstr[k-1].ptn.x-tmstr[k-1].ptc.x,
                 tmstr[k-1].ptn.y-tmstr[k-1].ptc.y,
                 tmstr[k-1].ptn.z-tmstr[k-1].ptc.z);
      rnorm=bf.norm(rp);
      if(rnorm<epsilon) //Two in the same place
      {
        tmstr[k-1].ptn.x-=0.0001f;
        tmstr[k-1].ptn.y-=0.0001f;
        tmstr[k-1].ptn.z-=0.0001f;
        rp=bf.setv(tmstr[k-1].ptn.x-tmstr[k-1].ptc.x,
                   tmstr[k-1].ptn.y-tmstr[k-1].ptc.y,
                   tmstr[k-1].ptn.z-tmstr[k-1].ptc.z);
        rnorm=bf.norm(rp);  
      }
      if(rnorm<lennc-delnc || rnorm>lennc+delnc)
      {
        rp=bf.scal(rp,lennc/rnorm);
        tmstr[k-1].ptn.x=tmstr[k-1].ptc.x+rp.x;
        tmstr[k-1].ptn.y=tmstr[k-1].ptc.y+rp.y;
        tmstr[k-1].ptn.z=tmstr[k-1].ptc.z+rp.z;
      }
      break;
    default:        
      break; 
  }
}

/* Move residue k respective to residue k-1 in decoy structure 'tmpstr'
 * mtype - the atom for movement
 *         0: move N  towards/from C[-1]
 *         1: move N  towards/from CA[-1]
 *         2: move CA towards/from C[-1]
 *         3: move CA towards/from CA[-1]
 *         4: move CA towards/from N
 *         5: move C  towards/from CA
 *         6: move C  towards/from N     */
void FoldDesignMovement::singlemoveLMPf(
  point3f *tmstr,
  int k,
  int mtype
)
{
  point3d rp;
  double rnorm;
  switch(mtype)
  {
    case 0: //1 c-1 n
      rp=bf.setv(tmstr[k].ptn.x-tmstr[k-1].ptc.x,
                 tmstr[k].ptn.y-tmstr[k-1].ptc.y,
                 tmstr[k].ptn.z-tmstr[k-1].ptc.z);
      rnorm=bf.norm(rp);
      if(rnorm<epsilon) //Two in the same place
      {
        tmstr[k].ptn.x+=0.0001f;
        tmstr[k].ptn.y+=0.0001f;
        tmstr[k].ptn.z+=0.0001f;
        rp=bf.setv(tmstr[k].ptn.x-tmstr[k-1].ptc.x,
                   tmstr[k].ptn.y-tmstr[k-1].ptc.y,
                   tmstr[k].ptn.z-tmstr[k-1].ptc.z);
        rnorm=bf.norm(rp);  
      }
      if(rnorm<lencn-delcn || rnorm>lencn+delcn)
      {
        rp=bf.scal(rp,lencn/rnorm);
        tmstr[k].ptn.x=tmstr[k-1].ptc.x+rp.x;
        tmstr[k].ptn.y=tmstr[k-1].ptc.y+rp.y;
        tmstr[k].ptn.z=tmstr[k-1].ptc.z+rp.z;
      } 
      break;
    case 1: //2 ca-1 n
      rp=bf.setv(tmstr[k].ptn.x-tmstr[k-1].x,
                 tmstr[k].ptn.y-tmstr[k-1].y,
                 tmstr[k].ptn.z-tmstr[k-1].z);
      rnorm=bf.norm(rp);
      if(rnorm<epsilon) //Two in the same place
      {
        tmstr[k].ptn.x+=0.0001f;
        tmstr[k].ptn.y+=0.0001f;
        tmstr[k].ptn.z+=0.0001f;
        rp=bf.setv(tmstr[k].ptn.x-tmstr[k-1].x,
                   tmstr[k].ptn.y-tmstr[k-1].y,
                   tmstr[k].ptn.z-tmstr[k-1].z);
        rnorm=bf.norm(rp);  
      }
      if(rnorm<lencan1-delcan1 || rnorm>lencan1+delcan1)
      {
        rp=bf.scal(rp,lencan1/rnorm);
        tmstr[k].ptn.x=tmstr[k-1].x+rp.x;
        tmstr[k].ptn.y=tmstr[k-1].y+rp.y;
        tmstr[k].ptn.z=tmstr[k-1].z+rp.z;
      }
      break;
    case 2: //4 c-1 ca
      rp=bf.setv(tmstr[k].x-tmstr[k-1].ptc.x,
                 tmstr[k].y-tmstr[k-1].ptc.y,
                 tmstr[k].z-tmstr[k-1].ptc.z);
      rnorm=bf.norm(rp);
      if(rnorm<epsilon) //Two in the same place
      {
        tmstr[k].x+=0.0001f;
        tmstr[k].y+=0.0001f;
        tmstr[k].z+=0.0001f;
        rp=bf.setv(tmstr[k].x-tmstr[k-1].ptc.x,
                   tmstr[k].y-tmstr[k-1].ptc.y,
                   tmstr[k].z-tmstr[k-1].ptc.z);
        rnorm=bf.norm(rp);  
      }
      if(rnorm<lencca1-delcca1 || rnorm>lencca1+delcca1)
      {
        rp=bf.scal(rp,lencca1/rnorm);
        tmstr[k].x=tmstr[k-1].ptc.x+rp.x;
        tmstr[k].y=tmstr[k-1].ptc.y+rp.y;
        tmstr[k].z=tmstr[k-1].ptc.z+rp.z;
      }
      break;
    case 3: //5 ca-1 ca
      rp=bf.setv(tmstr[k].x-tmstr[k-1].x,
                 tmstr[k].y-tmstr[k-1].y,
                 tmstr[k].z-tmstr[k-1].z);
      rnorm=bf.norm(rp);
      if(rnorm<epsilon) //Two in the same place
      {
        tmstr[k].x+=0.0001f;
        tmstr[k].y+=0.0001f;
        tmstr[k].z+=0.0001f;
        rp=bf.setv(tmstr[k].x-tmstr[k-1].x,
                   tmstr[k].y-tmstr[k-1].y,
                   tmstr[k].z-tmstr[k-1].z);
        rnorm=bf.norm(rp);  
      }
      if(rnorm<lencaca-delcaca || rnorm>lencaca+delcaca)
      {
        rp=bf.scal(rp,lencaca/rnorm);
        tmstr[k].x=tmstr[k-1].x+rp.x;
        tmstr[k].y=tmstr[k-1].y+rp.y;
        tmstr[k].z=tmstr[k-1].z+rp.z;
      }
      break;
    case 4: //3 n ca
      rp=bf.setv(tmstr[k].x-tmstr[k].ptn.x,
                 tmstr[k].y-tmstr[k].ptn.y,
                 tmstr[k].z-tmstr[k].ptn.z);
      rnorm=bf.norm(rp);
      if(rnorm<epsilon) //Two in the same place
      {
        tmstr[k].x+=0.0001f;
        tmstr[k].y+=0.0001f;
        tmstr[k].z+=0.0001f;
        rp=bf.setv(tmstr[k].x-tmstr[k].ptn.x,
                   tmstr[k].y-tmstr[k].ptn.y,
                   tmstr[k].z-tmstr[k].ptn.z);
        rnorm=bf.norm(rp);  
      }
      if(rnorm<lennca-delnca || rnorm>lennca+delnca)
      {
        rp=bf.scal(rp,lennca/rnorm);
        tmstr[k].x=tmstr[k].ptn.x+rp.x;
        tmstr[k].y=tmstr[k].ptn.y+rp.y;
        tmstr[k].z=tmstr[k].ptn.z+rp.z;
      }
      break;
    case 5: //6 ca c
      rp=bf.setv(tmstr[k].ptc.x-tmstr[k].x,
                 tmstr[k].ptc.y-tmstr[k].y,
                 tmstr[k].ptc.z-tmstr[k].z);
      rnorm=bf.norm(rp);
      if(rnorm<epsilon) //Two in the same place
      {
        tmstr[k].ptc.x+=0.0001f;
        tmstr[k].ptc.y+=0.0001f;
        tmstr[k].ptc.z+=0.0001f;
        rp=bf.setv(tmstr[k].ptc.x-tmstr[k].x,
                   tmstr[k].ptc.y-tmstr[k].y,
                   tmstr[k].ptc.z-tmstr[k].z);
        rnorm=bf.norm(rp);  
      }
      if(rnorm<lencac-delcac || rnorm>lencac+delcac)
      {
        rp=bf.scal(rp,lencac/rnorm);
        tmstr[k].ptc.x=tmstr[k].x+rp.x;
        tmstr[k].ptc.y=tmstr[k].y+rp.y;
        tmstr[k].ptc.z=tmstr[k].z+rp.z;
      }
      break;
    case 6: //7 n c
      rp=bf.setv(tmstr[k].ptc.x-tmstr[k].ptn.x,
                 tmstr[k].ptc.y-tmstr[k].ptn.y,
                 tmstr[k].ptc.z-tmstr[k].ptn.z);
      rnorm=bf.norm(rp);
      if(rnorm<epsilon) //Two in the same place
      {
        tmstr[k].ptc.x+=0.0001f;
        tmstr[k].ptc.y+=0.0001f;
        tmstr[k].ptc.z+=0.0001f;
        rp=bf.setv(tmstr[k].ptc.x-tmstr[k].ptn.x,
                   tmstr[k].ptc.y-tmstr[k].ptn.y,
                   tmstr[k].ptc.z-tmstr[k].ptn.z);
        rnorm=bf.norm(rp);
      }
      if(rnorm<lennc-delnc || rnorm>lennc+delnc)
      {
        rp=bf.scal(rp,lennc/rnorm);
        tmstr[k].ptc.x=tmstr[k].ptn.x+rp.x;
        tmstr[k].ptc.y=tmstr[k].ptn.y+rp.y;
        tmstr[k].ptc.z=tmstr[k].ptn.z+rp.z;
      }
      break;
    default: 
      break; 
  }
}

/* Movement M6 is an LMProt perturbation, which first randomly changes the
 * positions of backbone atoms in a selected segment and then tries to
 * restrict all the length of bonds and pseudo-bonds of backbone atoms (and
 * consequentially their bond angles) within the physically allowed region.
 */
bool FoldDesignMovement::mcfragsweepLMP(
  point3f *decstr,
  int numseq,
  double tbeta,
  double oldenergy,
  double *newenergy
)
{
  int i,j,k;
  int attemptNum=0;
  int trandpos,trand2;
  bool flagclash=true;
  double threshclash=0.01;
  point3d rp;
  double rnorm;
  double rmax;
  bool flagdone;
  int numiter;
  int outiter=0;
  int numLMP=0;
  do
  {
    flagclash=false;
    memcpy(tmstr,decstr,numseq*sizeof(point3f));
    trand2=2+int(6*Random());
    trandpos=1+int((numseq-trand2-1)*Random());
    int threshiter=trand2*200;
    if(threshiter>1000) threshiter=1000;
    double treject=Random();
    double ssrat=0;
     
    if(!flagclash) 
    {
      for(i=trandpos;i<trandpos+trand2;i++)
      {
        if(tmstr[i].ssm!='C') ssrat+=1;
      }

      ssrat/=double(trand2);
      if(ssrat>epsilon && treject>0.5*(1-ssrat)+0.05)
      {
        summcLMP[2]++;
        return false;
      }

      //perturbation
      for(i=trandpos;i<trandpos+trand2;i++)
      {
        rmax=0.25+1.25*Random();
        rp=bf.ranv_ergodic(rmax); //rp=bf.ranv(rmax);
        tmstr[i].x+=rp.x;
        tmstr[i].y+=rp.y;
        tmstr[i].z+=rp.z;
        rmax=0.25+1.25*Random();
        rp=bf.ranv_ergodic(rmax); //rp=bf.ranv(rmax);
        tmstr[i].ptn.x+=rp.x;
        tmstr[i].ptn.y+=rp.y;
        tmstr[i].ptn.z+=rp.z;
        rmax=0.25+1.25*Random();
        rp=bf.ranv_ergodic(rmax); //rp=bf.ranv(rmax);
        tmstr[i].ptc.x+=rp.x;
        tmstr[i].ptc.y+=rp.y;
        tmstr[i].ptc.z+=rp.z;
      }
      numiter=0;
      do
      {
        // movement from N to C terminal
        for(k=trandpos;k<trandpos+trand2;k++)
        {
          singlemoveLMPf(tmstr,k,int(Random()*7.0));
          singlemoveLMPf(tmstr,k,int(Random()*7.0));
          singlemoveLMPf(tmstr,k,int(Random()*7.0));
          singlemoveLMPf(tmstr,k,int(Random()*7.0));
          singlemoveLMPf(tmstr,k,int(Random()*7.0));
          singlemoveLMPf(tmstr,k,int(Random()*7.0));
          singlemoveLMPf(tmstr,k,int(Random()*7.0));
        }
        // movement from C to N terminal
        for(k=trandpos+trand2;k>trandpos;k--)
        {
          singlemoveLMPb(tmstr,k,int(Random()*7.0));
          singlemoveLMPb(tmstr,k,int(Random()*7.0));
          singlemoveLMPb(tmstr,k,int(Random()*7.0));
          singlemoveLMPb(tmstr,k,int(Random()*7.0));
          singlemoveLMPb(tmstr,k,int(Random()*7.0));
          singlemoveLMPb(tmstr,k,int(Random()*7.0));
          singlemoveLMPb(tmstr,k,int(Random()*7.0));
        }

        // check whether (pseudo-)bond length is physical
        flagdone=true;
        for(k=trandpos;k<=trandpos+trand2;k++)
        {
          // CA[-1] - CA
          rp=bf.setv(tmstr[k].x-tmstr[k-1].x,
                     tmstr[k].y-tmstr[k-1].y,
                     tmstr[k].z-tmstr[k-1].z);
          rnorm=bf.norm(rp);
          if(rnorm<lencaca-delcaca || rnorm>lencaca+delcaca)
          {
            flagdone=false;
            break;
          }
          // CA[-1] - N
          rp=bf.setv(tmstr[k].ptn.x-tmstr[k-1].x,
                     tmstr[k].ptn.y-tmstr[k-1].y,
                     tmstr[k].ptn.z-tmstr[k-1].z);
          rnorm=bf.norm(rp);
          if(rnorm<lencan1-delcan1 || rnorm>lencan1+delcan1)
          {
            flagdone=false;
            break;
          }
          // C[-1] - N
          rp=bf.setv(tmstr[k].ptn.x-tmstr[k-1].ptc.x,
                     tmstr[k].ptn.y-tmstr[k-1].ptc.y,
                     tmstr[k].ptn.z-tmstr[k-1].ptc.z);
          rnorm=bf.norm(rp);
          if(rnorm<lencn-delcn || rnorm>lencn+delcn)
          {
            flagdone=false;
            break;
          }
          // C[-1] - CA
          rp=bf.setv(tmstr[k].x-tmstr[k-1].ptc.x,
                     tmstr[k].y-tmstr[k-1].ptc.y,
                     tmstr[k].z-tmstr[k-1].ptc.z);
          rnorm=bf.norm(rp);
          if(rnorm<lencca1-delcca1 || rnorm>lencca1+delcca1)
          {
            flagdone=false;
            break;
          }
          if(k==trandpos+trand2) break;
          // N - CA
          rp=bf.setv(tmstr[k].x-tmstr[k].ptn.x,
                     tmstr[k].y-tmstr[k].ptn.y,
                     tmstr[k].z-tmstr[k].ptn.z);
          rnorm=bf.norm(rp);
          if(rnorm<lennca-delnca || rnorm>lennca+delnca)
          {
            flagdone=false;
            break;
          }
          // N - C
          rp=bf.setv(tmstr[k].ptc.x-tmstr[k].ptn.x,
                     tmstr[k].ptc.y-tmstr[k].ptn.y,
                     tmstr[k].ptc.z-tmstr[k].ptn.z);
          rnorm=bf.norm(rp);
          if(rnorm<lennc-delnc || rnorm>lennc+delnc)
          {
            flagdone=false;
            break;
          }
          // CA - C
          rp=bf.setv(tmstr[k].ptc.x-tmstr[k].x,
                     tmstr[k].ptc.y-tmstr[k].y,
                     tmstr[k].ptc.z-tmstr[k].z);
          rnorm=bf.norm(rp);
          if(rnorm<lencac-delcac || rnorm>lencac+delcac)
          {
            flagdone=false;
            break;
          }
        }
        numiter++;
      } while(!flagdone && numiter<threshiter);
      numLMP++;
      if(!flagdone)
      {
        outiter++;
        summcLMP[2]++;
        return false;
      }
      else
      {
        if(trandpos+trand2+2>=numseq) geometry.str2torp(tmstr,numseq,trandpos,numseq-1);
        else geometry.str2torp(tmstr,numseq,trandpos,trandpos+trand2+2);
      }
      flagclash=false;
    }    
    else if(attemptNum>=MAX_ITERATION)
    {
      summcLMP[2]++;
      return false;
    }
  } while(flagclash && outiter<50);
  if(outiter==50)
  {
    summcLMP[2]++;
    return false;
  }

  if(tbeta<0) // skip energy calculation
  {
    summcLMP[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }

  *newenergy=energyFunction.calcTotalEnergy(tmstr,numseq);
  if(flagcaca && enelist[14]<enelistbk[14] && 
     (*newenergy<oldenergy || (*newenergy-oldenergy)<-0.01*oldenergy))
  {
    summcLMP[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }
  if(*newenergy<oldenergy)
  {
    summcLMP[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }
  else
  {
    double trand;
    trand=Random();
    if(trand<exp(-tbeta*(*newenergy-oldenergy)))
    {
      summcLMP[1]++;
      memcpy(decstr,tmstr,numseq*sizeof(point3f));
      return true;
    }
    else
    {
      summcLMP[2]++;
      return false;
    }
  }
}

bool FoldDesignMovement::mcfragsweepCCD3(
  point3f *decstr,
  int numseq,
  int lps, 
  int lpe,
  int lpt,
  point3d pt1,
  point3d pt2,
  point3d pt3,
  double oldenergy,
  double *newenergy
)
{
  int i;
  int trandpos;
  bool flagphi,flagtor;
  int flagpt;
  point3d tp[3];//dest for lpt-1 lpt lpt+1
  tp[0].x=pt1.x;tp[0].y=pt1.y;tp[0].z=pt1.z;
  tp[1].x=pt2.x;tp[1].y=pt2.y;tp[1].z=pt2.z;
  tp[2].x=pt3.x;tp[2].y=pt3.y;tp[2].z=pt3.z;
  point3d pcur,p12,p13;
  int numiter=0;
  double tdist[3],ttheta,tphi,tpsi,tinner;
  memcpy(tmstr2,decstr,numseq*sizeof(point3f));
  flagpt=0;
  int threshiter=25*(lpe-lps+1);
  if(threshiter>600) threshiter=600;
  do
  {
    tdist[0]=10000;
    tdist[1]=10000;
    tdist[2]=10000;
    trandpos=lps+int((lpe-lps+1)*Random());
    if(Random()<0.5) flagphi=true;//change phi n ca
    else flagphi=false;//change psi ca c in pos+1
    if(tmstr2[trandpos].ss3=='H' && tmstr2[trandpos+1].ss3=='H')
    {
      numiter++;
      continue;
    }
    else if(tmstr2[trandpos].ss3=='H')
      flagphi=false;
    else if(tmstr2[trandpos].ss3=='E' && tmstr2[trandpos+1].ss3=='E')
    {
      numiter++;
      continue;
    }
    else if(tmstr2[trandpos].ss3=='E')
    {
      flagphi=false;
    }
    flagpt=(flagpt+1)%3;
    if(flagphi)
    {
      p12=bf.setv(tmstr2[trandpos].x-tmstr2[trandpos].ptn.x,
                  tmstr2[trandpos].y-tmstr2[trandpos].ptn.y,
                  tmstr2[trandpos].z-tmstr2[trandpos].ptn.z);
      p13=bf.setv(tmstr2[lpt+flagpt-1].x-tmstr2[trandpos].ptn.x,
                  tmstr2[lpt+flagpt-1].y-tmstr2[trandpos].ptn.y,
                  tmstr2[lpt+flagpt-1].z-tmstr2[trandpos].ptn.z);
      tinner=bf.angv(p12,p13);

      //In one line no affect when rotate
      if(tinner<epsilon || PI-tinner<epsilon || bf.norm(p12)<epsilon)
      {
        numiter++;
        continue;
      }
      pcur.x=tmstr2[lpt+flagpt-1].x;
      pcur.y=tmstr2[lpt+flagpt-1].y;
      pcur.z=tmstr2[lpt+flagpt-1].z;
      p13.x=tmstr2[trandpos].x;
      p13.y=tmstr2[trandpos].y;
      p13.z=tmstr2[trandpos].z;
      p12.x=tmstr2[trandpos].ptn.x;
      p12.y=tmstr2[trandpos].ptn.y;
      p12.z=tmstr2[trandpos].ptn.z;
      ttheta=bf.phi(pcur.x,pcur.y,pcur.z,
                    p12.x,p12.y,p12.z,
                    p13.x,p13.y,p13.z,
                    tp[flagpt].x,tp[flagpt].y,tp[flagpt].z);
      tphi=tmstr2[trandpos].tor[2]+ttheta;
      if(tphi>360.0) tphi-=360.0;
      tmstr2[trandpos].tor[2]=tphi;
    }
    else
    {
      p12=bf.setv(tmstr2[trandpos].ptc.x-tmstr2[trandpos].x,
                  tmstr2[trandpos].ptc.y-tmstr2[trandpos].y,
                  tmstr2[trandpos].ptc.z-tmstr2[trandpos].z);
      p13=bf.setv(tmstr2[lpt+flagpt-1].x-tmstr2[trandpos].x,
                  tmstr2[lpt+flagpt-1].y-tmstr2[trandpos].y,
                  tmstr2[lpt+flagpt-1].z-tmstr2[trandpos].z);
      tinner=bf.angv(p12,p13);

      //In one line no affect when rotate
      if(tinner<epsilon || PI-tinner<epsilon || bf.norm(p12)<epsilon) 
      {
        numiter++;
        continue;
      }
      pcur.x=tmstr2[lpt+flagpt-1].x;
      pcur.y=tmstr2[lpt+flagpt-1].y;
      pcur.z=tmstr2[lpt+flagpt-1].z;
      p12.x=tmstr2[trandpos].x;
      p12.y=tmstr2[trandpos].y;
      p12.z=tmstr2[trandpos].z;
      p13.x=tmstr2[trandpos].ptc.x;
      p13.y=tmstr2[trandpos].ptc.y;
      p13.z=tmstr2[trandpos].ptc.z;
      ttheta=bf.phi(pcur.x,pcur.y,pcur.z,
                    p12.x,p12.y,p12.z,
                    p13.x,p13.y,p13.z,
                    tp[flagpt].x,tp[flagpt].y,tp[flagpt].z);
      tpsi=tmstr2[trandpos+1].tor[0]+ttheta;
      if(tpsi>360.0) tpsi-=360.0;
      tmstr2[trandpos+1].tor[0]=tpsi;
    }
    numiter++;
    if(trandpos==0) flagtor=geometry.tor2str(tmstr2,numseq,3);
    else flagtor=geometry.tor2strp(tmstr2,numseq,trandpos);
    if(!flagtor) printf("tor2str wrong in CCD3 %d\n",trandpos);
    for(i=0;i<3;i++)
    {
      pcur.x=tmstr2[lpt+i-1].x-tp[i].x;
      pcur.y=tmstr2[lpt+i-1].y-tp[i].y;
      pcur.z=tmstr2[lpt+i-1].z-tp[i].z;
      tdist[i]=pcur.x*pcur.x+pcur.y*pcur.y+pcur.z*pcur.z;
    }
  } while(numiter<threshiter && (tdist[0]>2.50 || tdist[1]>2.50 || tdist[2]>2.50));
  if(numiter<threshiter)
  {
    memcpy(decstr,tmstr2,numseq*sizeof(point3f));
    return true;
  }
  else
  {
    memcpy(decstr,tmstr2,numseq*sizeof(point3f));
    return false;
  }
}

void FoldDesignMovement::getfrombbptable(
  double tval,
  segment *sgval,
  paircont* tbbptb,
  int tnumbbptb
)
{
  //(i i+1]
  int i;
  if(tnumbbptb<2)
  {   
    if(tval>tbbptb[0].dist) printf("wrong %f %f\n",tval,tbbptb[0].dist);
    sgval->init=tbbptb[0].init;
    sgval->term=tbbptb[0].term;
    return;
  }
  i=tnumbbptb/2;
  if(tval>tbbptb[i-1].dist) 
    getfrombbptable(tval,sgval,tbbptb+i,tnumbbptb-i);
  else getfrombbptable(tval,sgval,tbbptb,i);
}

/* In movement M10 (beta pairing), one helix is moved close to another one.
 * The probability of their distance and torsion-angle distribution is the
 * same as that in the beta-pairing energy term. The linkage region between
 * the two beta strands will be rebuilt by mcfragsweepCCD3. As a beta-strand
 * is likely to pair with another one which has similar number of residues
 * to form a beta-sheet, we precalculated the probability for every pair of
 * residues which may form a beta-pair based on their SS types and positions
 * in the SS elements. A pair of residues whose predicted SS types are
 * strands have a higher probability than those with SS types equal to coils
 * and helices. This means even all alpha protein can accept M10 movements,
 * albeit at a much lower acceptance rate than that in a beta protein.
 * During M10, we select the residue pair based on the precalculated
 * probabilities, as implemented by getfrombbptable. The possibilities of
 * forming a beta-pair in antiparallel and parallel sheets are 75% and 25%,
 * respectively. */
bool FoldDesignMovement::mcfragsweepbbp(
  point3f *decstr,
  int numseq,
  double tbeta,
  double oldenergy,
  double *newenergy
)
{
  int i,j,k;
  int attemptNum=0;
  int trand,trand2;
  int typebb;
  bool flagat,flaglr;
  point3d tp[15],ap[15],fp[3];
  double angle1,angle2,angle3;
  double rmat[9],irmat[9],rmat2[9];
  double toldenergy=oldenergy;
  double tnewenergy;
  bool flagclash=true;
  bool flagtor;
  double threshclash=0.01;
  segment sgval;
  if(numbbptb<3)
  {
    summcbbp[2]++;
    return false;
  }

  do
  {
    flagclash=false;
    memcpy(tmstr,decstr,numseq*sizeof(point3f));
    getfrombbptable(Random(),&sgval,bbptb,numbbptb);
    trand=sgval.init;
    trand2=sgval.term;
    if(abs(decstr[trand].indr-trand2)<3 || abs(decstr[trand].indl-trand2)<3)
    { //already paired
      summcbbp[2]++;
      return false;
    }
 
    if(!flagclash) 
    {
      if(Random()<0.75) flagat=false;//anti
      else flagat=true;//para
      if(Random()<0.5) flaglr=false;//left
      else flaglr=true;//right
      if     (flagat  && !flaglr) typebb=0;//paraleft
      else if(flagat  && flaglr)  typebb=1;//pararight
      else if(!flagat && !flaglr) typebb=2;//antileft
      else if(!flagat && flaglr)  typebb=3;//antiright

      if(typebb==0 && (decstr[trand].indl!=-1 || decstr[trand2].indr!=-1))
      {
        flaglr=true;
        typebb=1;
      }
      else if(typebb==1 && (decstr[trand].indr!=-1 || decstr[trand2].indl!=-1))
      {
        flaglr=false;
        typebb=0;
      }
      else if(typebb==2 && (decstr[trand].indl!=-1 || decstr[trand2].indr!=-1))
      {
        flaglr=true;
        typebb=3;
      }
      else if(typebb==3 && (decstr[trand].indr!=-1 || decstr[trand2].indl!=-1))
      {
        flaglr=false;
        typebb=2;
      }

      //flagat
      if(typebb==0 && (decstr[trand].indl!=-1 || decstr[trand2].indr!=-1))
      {
        flagat=false;
        typebb=2;
      }
      else if(typebb==1 && (decstr[trand].indr!=-1 || decstr[trand2].indl!=-1))
      {
        flagat=false;
        typebb=3;
      }
      else if(typebb==2 && (decstr[trand].indl!=-1 || decstr[trand2].indr!=-1))
      {
        flagat=true;
        typebb=0;
      }
      else if(typebb==3 && (decstr[trand].indr!=-1 || decstr[trand2].indl!=-1))
      {
        flagat=true;
        typebb=1;
      }
      if(typebb==0 && (decstr[trand].indr!=-1 && decstr[trand].tpr==3))
      {
        if(Random()<0.90)
        {
          flagat=false;
          typebb=2;
        }
      }
      else if(typebb==1 && (decstr[trand].indl!=-1 && decstr[trand].tpl==2))
      {
        if(Random()<0.90)
        {
          flagat=false;
          typebb=3;
        }
      }
      else if(typebb==2 && (decstr[trand].indr!=-1 && decstr[trand].tpr==1))
      {
        if(Random()<0.90)
        {
          flagat=true;
          typebb=0;
        }
      }
      else if(typebb==3 && (decstr[trand].indl!=-1 && decstr[trand].tpl==0))
      {
        if(Random()<0.90)
        {
          flagat=true;
          typebb=1;
        }
      }
    
      tp[0]=bf.setv(tmstr[trand-1].x,tmstr[trand-1].y,tmstr[trand-1].z);
      tp[1]=bf.setv(tmstr[trand].x,tmstr[trand].y,tmstr[trand].z);
      tp[2]=bf.setv(tmstr[trand+1].x,tmstr[trand+1].y,tmstr[trand+1].z);
      tp[3]=bf.minu(tp[0],tp[1]);
      tp[4]=bf.minu(tp[2],tp[1]);
      tp[5]=bf.prod(tp[3],tp[4]);
      tp[5]=bf.unit(tp[5]);
      if(bf.norm(tp[5])<1e-10)
      {
        summcbbp[2]++;
        return false;
      }

      tp[6]=bf.minu(tp[0],tp[2]);
      tp[6]=bf.unit(tp[6]);
      bf.v2rot(tp[5],rmat);
      bf.copymat(rmat,3,3,irmat);
      bf.rinv(rmat,3);//tp5 to 001
      tp[10]=bf.mmat(rmat,tp[6]);
      angle1=atan2(tp[10].y,tp[10].x);

      /////////////////////////////////////////////////////////////////////
      ap[0]=bf.setv(tmstr[trand2-1].x,tmstr[trand2-1].y,tmstr[trand2-1].z);
      ap[1]=bf.setv(tmstr[trand2].x,tmstr[trand2].y,tmstr[trand2].z);
      ap[2]=bf.setv(tmstr[trand2+1].x,tmstr[trand2+1].y,tmstr[trand2+1].z);
      ap[3]=bf.minu(ap[0],ap[1]);
      ap[4]=bf.minu(ap[2],ap[1]);
      ap[5]=bf.prod(ap[3],ap[4]);
      ap[5]=bf.unit(ap[5]);
      if(bf.norm(ap[5])<1e-10)
      {
        summcbbp[2]++;
        return false;
      }

      ap[6]=bf.minu(ap[0],ap[2]);
      ap[6]=bf.unit(ap[6]);
      if(!flagat) ap[5]=bf.scal(ap[5],-1);
      bf.v2rot(ap[5],rmat);
      bf.rinv(rmat,3);//ap5 to 001
      ap[10]=bf.mmat(rmat,ap[6]); 
      ap[11]=bf.mmat(rmat,ap[3]);
      ap[12]=bf.mmat(rmat,ap[4]);
      angle2=atan2(ap[10].y,ap[10].x);
      if(!flagat) angle3=angle1-angle2+PI+Random()*0.20-0.10;
      else angle3=angle1-angle2+Random()*0.20-0.10;
      bf.a2rot(angle3,rmat2);
      ap[10]=bf.mmat(rmat2,ap[10]);
      angle2=atan2(ap[10].y,ap[10].x);
      ap[13]=bf.mmat(rmat2,ap[11]);
      ap[14]=bf.mmat(rmat2,ap[12]);
      tp[13]=bf.mmat(irmat,ap[13]);
      tp[14]=bf.mmat(irmat,ap[14]);
      tp[11]=tp[13];
      tp[12]=tp[14];
      int type;
      if(flagat && !flaglr)
      { //paraleft norm same type0
        type=1;
        tp[9]=tp[5];
        fp[1]=bf.minu(tp[1],bf.scal(tp[9],4.81+0.1*Random()));
      }
      else if(flagat && flaglr)
      { //pararight norm same type1
        type=2;
        tp[9]=tp[5];
        fp[1]=bf.addv(tp[1],bf.scal(tp[9],4.86+0.1*Random()));
      }
      else if(!flagat && !flaglr)
      { //antileft norm diff type2
        type=3;
        tp[9]=tp[5];
        fp[1]=bf.minu(tp[1],bf.scal(tp[9],4.49+0.1*Random()));
      }
      else if(!flagat && flaglr)
      { //antiright norm diff type3
        type=4;
        tp[9]=tp[5];
        fp[1]=bf.addv(tp[1],bf.scal(tp[9],5.24+0.1*Random()));
      }
      fp[0]=bf.addv(fp[1],tp[11]);
      fp[2]=bf.addv(fp[1],tp[12]);
      mcfragsweepCCD3(tmstr,numseq,trand,trand2-1,trand2,fp[0],fp[1],fp[2],
                      toldenergy,&tnewenergy);
      flagtor=geometry.tor2str(tmstr,numseq,3);
      if(!flagtor) printf("tor2str wrong in bbp\n");
      flagclash=false;
    }
    else if(attemptNum>=MAX_ITERATION)
    {
      summcbbp[2]++;
      return false;
    }
  }while(flagclash);

  if(tbeta<0) // skip energy calculation
  {
    summcbbp[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }
  
  *newenergy=energyFunction.calcTotalEnergy(tmstr,numseq);
  if(*newenergy<oldenergy)
  {
    summcbbp[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }
  else
  {
    double trand;
    trand=Random();
    if(trand<exp(-tbeta*(*newenergy-oldenergy)))
    {
      summcbbp[1]++;
      memcpy(decstr,tmstr,numseq*sizeof(point3f));
      return true;
    }
    else
    {
      summcbbp[2]++;
      return false;
    }
  }
}

/* M3: change omega (CA[-1] - C[-1] - N - CA) backbone torsion angle. */
bool FoldDesignMovement::mcfragsweepome(
  point3f *decstr,
  int numseq,
  double tbeta,
  double oldenergy,
  double *newenergy
)
{
  int i;
  int attemptNum=0;
  int trandpos;
  double trand2;
  bool flagclash=true;
  bool flagtor;
  double threshclash=0.01;
  do
  {
    flagclash=false;
    memcpy(tmstr,decstr,numseq*sizeof(point3f));
    trandpos=int((numseq)*Random());
    if(tmstr[trandpos].ssm!='C' && Random()>0.1)
    {
      summcome[2]++;
      return false;
    }

    if(!flagclash)
    {
      trand2=double(16*Random())-8;
      tmstr[trandpos].tor[1]+=trand2;
      if(tmstr[trandpos].tor[1]<0) tmstr[trandpos].tor[1]+=360;
      else if(tmstr[trandpos].tor[1]>=360) tmstr[trandpos].tor[1]-=360;
      if(tmstr[trandpos].tor[1]>190 || tmstr[trandpos].tor[1]<170)
      {
        summcome[2]++;
        return false;
      }
      if(trandpos==0) flagtor=geometry.tor2str(tmstr,numseq,3);
      else flagtor=geometry.tor2strp(tmstr,numseq,trandpos);

      if(!flagtor) printf("tor2str wrong in ome %d\n",trandpos);
      flagclash=false;
    }
    else if(attemptNum>=MAX_ITERATION)
    {
      summcome[2]++;
      return false;
    }
  } while(flagclash);

  if(tbeta<0) // skip energy calculation
  {
    summcome[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }

  *newenergy=energyFunction.calcTotalEnergy(tmstr,numseq);
  if(*newenergy<oldenergy)
  {
    summcome[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }
  else
  {
    double trand;
    trand=Random();
    if(trand<exp(-tbeta*(*newenergy-oldenergy)))
    {
      summcome[1]++;
      memcpy(decstr,tmstr,numseq*sizeof(point3f));
      return true;
    }
    else
    {
      summcome[2]++;
      return false;
    }
  }
}

/* Movement M1 randomly change one bond length of a randomly selected
 * residue by +-0.24 A */
bool FoldDesignMovement::mcfragsweeplen(
  point3f *decstr,
  int numseq,
  double tbeta,
  double oldenergy,
  double *newenergy
)
{
  int i;
  int attemptNum=0;
  int trandpos;
  double trand2;
  bool flagclash=true;
  bool flagtor;
  double threshclash=0.01;
  do
  {
    flagclash=false;
    memcpy(tmstr,decstr,numseq*sizeof(point3f));
    trandpos=int((numseq)*Random());
 
    if(!flagclash) 
    {
      if(tmstr[trandpos].ssm!='C' && Random()>0.1)
      {
        summclen[2]++;
        return false;
      }
      int lentype=int(Random()*3.0);
      trand2=double(0.48*Random())-0.24;
      tmstr[trandpos].len[lentype]+=trand2;
      if((lentype==0 && fabs(tmstr[trandpos].len[lentype]-lencn)>3*delcn) ||
         (lentype==1 && fabs(tmstr[trandpos].len[lentype]-lennca)>3*delnca) ||
         (lentype==2 && fabs(tmstr[trandpos].len[lentype]-lencac)>3*delcac))
      {
        summclen[2]++;
        return false;
      }
      if(trandpos==0) flagtor=geometry.tor2str(tmstr,numseq,3);
      else flagtor=geometry.tor2strp(tmstr,numseq,trandpos);
      if(!flagtor) printf("tor2str wrong in len %d\n",trandpos);
      flagclash=false;
    }
    else if(attemptNum>=MAX_ITERATION)
    {
      summclen[2]++;
      return false;
    }
  } while(flagclash);

  if(tbeta<0) // skip energy calculation
  {
    summclen[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }
    
  *newenergy=energyFunction.calcTotalEnergy(tmstr,numseq);
  if(*newenergy<oldenergy)
  {
    summclen[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }
  else
  {
    double trand;
    trand=Random();
    if(trand<exp(-tbeta*(*newenergy-oldenergy)))
    {
      summclen[1]++;
      memcpy(decstr,tmstr,numseq*sizeof(point3f));
      return true;
    }
    else
    {
      summclen[2]++;
      return false;
    }
  }
}

/* M2 randomly change one bond angle of a randomly selected residue */
bool FoldDesignMovement::mcfragsweepang(
  point3f *decstr,
  int numseq,
  double tbeta,
  double oldenergy,
  double *newenergy
)
{
  int i;
  int attemptNum=0;
  int trandpos;
  double trand2;
  bool flagclash=true;
  bool flagtor;
  double threshclash=0.01;
  do
  {
    flagclash=false;
    memcpy(tmstr,decstr,numseq*sizeof(point3f));
    trandpos=int((numseq)*Random());
    if(tmstr[trandpos].ssm!='C' && Random()>0.1)
    {
      summcang[2]++;
      return false;
    }

    if(!flagclash) 
    {
      int angtype=int(Random()*3.0);
      trand2=double(20*Random())-10;
      tmstr[trandpos].ang[angtype]+=trand2;
      if((angtype==0 && fabs(tmstr[trandpos].ang[angtype]-angcacn)>2.5*2.009) ||
         (angtype==1 && fabs(tmstr[trandpos].ang[angtype]-angcnca)>2.5*2.227) ||
         (angtype==2 && fabs(tmstr[trandpos].ang[angtype]-angncac)>2*2.818))
      {
        summcang[2]++;
        return false;
      }
      if(trandpos==0) flagtor=geometry.tor2str(tmstr,numseq,3);
      else flagtor=geometry.tor2strp(tmstr,numseq,trandpos);
      if(!flagtor) printf("tor2str wrong in ang %d\n",trandpos);
      flagclash=false;
    }
    else if(attemptNum>=MAX_ITERATION)
    {
      summcang[2]++;
      return false;
    }
  } while(flagclash);

  if(tbeta<0) // skip energy calculation
  {
    summcang[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }

  *newenergy=energyFunction.calcTotalEnergy(tmstr,numseq);
  if(*newenergy<oldenergy)
  {
    summcang[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }
  else
  {
    double trand;
    trand=Random();
    if(trand<exp(-tbeta*(*newenergy-oldenergy)))
    {
      summcang[1]++;
      memcpy(decstr,tmstr,numseq*sizeof(point3f));
      return true;
    }
    else
    {
      summcang[2]++;
      return false;
    }
  }
}

bool FoldDesignMovement::mcfragsweepCCD4(
  point3f *decstr,
  int numseq,
  int lps,
  int lpe,
  point3d pt1,
  int ind1,
  point3d pt2,
  int ind2,
  double oldenergy,
  double *newenergy
)
{
  int trandpos;
  bool flagphi,flagtor;
  int flagpt,indp[2];
  indp[0]=ind1;indp[1]=ind2;
  point3d tp[2]; 
  tp[0].x=pt1.x;tp[0].y=pt1.y;tp[0].z=pt1.z;
  tp[1].x=pt2.x;tp[1].y=pt2.y;tp[1].z=pt2.z;
  point3d pcur,p12,p13;
  int numiter=0;
  double tdist[2],ttheta,tphi,tpsi,tinner;
  memcpy(tmstr2,decstr,numseq*sizeof(point3f));
  flagpt=0;
  int threshiter=25*(lpe-lps+1);
  if(threshiter>600) threshiter=600;
  do
  {
    tdist[0]=10000;
    tdist[1]=10000;
    trandpos=lps+int((lpe-lps+1)*Random());
    if(Random()<0.5) flagphi=true;//change phi n ca
    else flagphi=false;//change psi ca c in pos+1
    if(tmstr2[trandpos].ss3=='H' && tmstr2[trandpos+1].ss3=='H')
    {
      numiter++;
      continue;
    }
    else if(tmstr2[trandpos].ss3=='H')
    {
      flagphi=false;
    }
    else if(tmstr2[trandpos].ss3=='E' && tmstr2[trandpos+1].ss3=='E')
    {
      numiter++;
      continue;
    }
    else if(tmstr2[trandpos].ss3=='E')
    {
      flagphi=false;
    }
    flagpt=(flagpt+1)%2;
    if(flagphi)
    {
      p12=bf.setv(tmstr2[trandpos].x-tmstr2[trandpos].ptn.x,
                  tmstr2[trandpos].y-tmstr2[trandpos].ptn.y,
                  tmstr2[trandpos].z-tmstr2[trandpos].ptn.z);
      p13=bf.setv(tmstr2[indp[flagpt]].x-tmstr2[trandpos].ptn.x,
                  tmstr2[indp[flagpt]].y-tmstr2[trandpos].ptn.y,
                  tmstr2[indp[flagpt]].z-tmstr2[trandpos].ptn.z);
      tinner=bf.angv(p12,p13);

      //In one line no affect when rotate
      if(tinner<epsilon || PI-tinner<epsilon || bf.norm(p12)<epsilon)
      {
        numiter++;
        continue;
      }
      pcur.x=tmstr2[indp[flagpt]].x;
      pcur.y=tmstr2[indp[flagpt]].y;
      pcur.z=tmstr2[indp[flagpt]].z;
      p13.x=tmstr2[trandpos].x;
      p13.y=tmstr2[trandpos].y;
      p13.z=tmstr2[trandpos].z;
      p12.x=tmstr2[trandpos].ptn.x;
      p12.y=tmstr2[trandpos].ptn.y;
      p12.z=tmstr2[trandpos].ptn.z;
      ttheta=bf.phi(pcur.x,pcur.y,pcur.z,
                    p12.x,p12.y,p12.z,
                    p13.x,p13.y,p13.z,
                    tp[flagpt].x,tp[flagpt].y,tp[flagpt].z);
      tphi=tmstr2[trandpos].tor[2]+ttheta;
      if(tphi>360.0) tphi-=360.0;
      tmstr2[trandpos].tor[2]=tphi;   
    }
    else
    {
      p12=bf.setv(tmstr2[trandpos].ptc.x-tmstr2[trandpos].x,
                  tmstr2[trandpos].ptc.y-tmstr2[trandpos].y,
                  tmstr2[trandpos].ptc.z-tmstr2[trandpos].z);
      p13=bf.setv(tmstr2[indp[flagpt]].x-tmstr2[trandpos].x,
                  tmstr2[indp[flagpt]].y-tmstr2[trandpos].y,
                  tmstr2[indp[flagpt]].z-tmstr2[trandpos].z);
      tinner=bf.angv(p12,p13);

      //In one line no affect when rotate
      if(tinner<epsilon || PI-tinner<epsilon || bf.norm(p12)<epsilon)
      {
        numiter++;
        continue;
      }
      pcur.x=tmstr2[indp[flagpt]].x;
      pcur.y=tmstr2[indp[flagpt]].y;
      pcur.z=tmstr2[indp[flagpt]].z;
      p12.x=tmstr2[trandpos].x;
      p12.y=tmstr2[trandpos].y;
      p12.z=tmstr2[trandpos].z;
      p13.x=tmstr2[trandpos].ptc.x;
      p13.y=tmstr2[trandpos].ptc.y;
      p13.z=tmstr2[trandpos].ptc.z;
      ttheta=bf.phi(pcur.x,pcur.y,pcur.z,
                    p12.x,p12.y,p12.z,
                    p13.x,p13.y,p13.z,
                    tp[flagpt].x,tp[flagpt].y,tp[flagpt].z);
      tpsi=tmstr2[trandpos+1].tor[0]+ttheta;
      if(tpsi>360.0) tpsi-=360.0;
      tmstr2[trandpos+1].tor[0]=tpsi;
    }       
    numiter++;
    if(trandpos==0) flagtor=geometry.tor2str(tmstr2,numseq,3);
    else flagtor=geometry.tor2strp(tmstr2,numseq,trandpos);
    if(!flagtor)
    {
      printf("tor2str wrong in CCD4 %d\n",trandpos);
    }
    pcur.x=tmstr2[ind1].x-tp[0].x;
    pcur.y=tmstr2[ind1].y-tp[0].y;
    pcur.z=tmstr2[ind1].z-tp[0].z;
    tdist[0]=pcur.x*pcur.x+pcur.y*pcur.y+pcur.z*pcur.z;
    pcur.x=tmstr2[ind2].x-tp[1].x;
    pcur.y=tmstr2[ind2].y-tp[1].y;
    pcur.z=tmstr2[ind2].z-tp[1].z;
    tdist[1]=pcur.x*pcur.x+pcur.y*pcur.y+pcur.z*pcur.z;
  } while(numiter<threshiter && (tdist[0]>6.50 || tdist[1]>2.50));
  if(numiter<threshiter)
  {
    memcpy(decstr,tmstr2,numseq*sizeof(point3f));
    return true;
  }
  else
  {
    memcpy(decstr,tmstr2,numseq*sizeof(point3f));
    return false;
  }
}

/* In movement M9 (helix packing), one helix is moved close to another one.
 * The probability of their distance and torsion-angle distribution is the
 * same as that in the helix-packing energy term. The linkage region
 * between the two helices will be rebuilt by mcfragsweepCCD4 */
bool FoldDesignMovement::mcfragsweepaaa(
  point3f *decstr,
  int numseq,
  double tbeta,
  double oldenergy,
  double *newenergy
)
{
  int i,j;
  int attemptNum=0;
  int trandpos;
  double trand2,tlength;
  int leng1,leng2,leng3,indbin;
  double rmat[9];
  point3d tp[12];
  double tnewenergy;
  bool flagclash=true;
  bool flagtor;
  double threshclash=0.01;
  sssegment *sse=inputInfo.getSSE();

  if(numsalpha==0)
  {
    summcaaa[2]++;
    return false;
  }

  do
  {
    flagclash=false;
    memcpy(tmstr,decstr,numseq*sizeof(point3f));
    trandpos=int((numsalpha)*Random());
        
    if(!flagclash) 
    {
      leng1=sse[alphasind[trandpos]].term-sse[alphasind[trandpos]].init+1;
      leng2=sse[alphasind[trandpos]+1].term-sse[alphasind[trandpos]+1].init+1;
      leng3=sse[alphasind[trandpos]+2].term-sse[alphasind[trandpos]+2].init+1;
      if(leng1>4 && leng3>4 && leng3<20)
      {
        indbin=0;
        trand2=Random();
        for(i=0;i<30;i++)
        {
          if(trand2<=angleaa[i])
          {
            indbin=i;
            break;
          }
        }
        trand2=(indbin+Random())/30.0*PI;
      }
      else trand2=Random()*PI;
      tp[0]=bf.setv(tmstr[sse[alphasind[trandpos]].init].x,
                    tmstr[sse[alphasind[trandpos]].init].y,
                    tmstr[sse[alphasind[trandpos]].init].z);
      tp[1]=bf.setv(tmstr[sse[alphasind[trandpos]].term].x,
                    tmstr[sse[alphasind[trandpos]].term].y,
                    tmstr[sse[alphasind[trandpos]].term].z);
      tp[2]=bf.setv(tmstr[sse[alphasind[trandpos]+2].init].x,
                    tmstr[sse[alphasind[trandpos]+2].init].y,
                    tmstr[sse[alphasind[trandpos]+2].init].z);//for new begin more flexible
      tp[3]=bf.setv(tmstr[sse[alphasind[trandpos]+2].term].x,
                    tmstr[sse[alphasind[trandpos]+2].term].y,
                    tmstr[sse[alphasind[trandpos]+2].term].z);
      tp[4]=bf.minu(tp[3],tp[2]);
      tlength=bf.norm(tp[4]);
      tp[5]=bf.rana(trand2);
      tp[6]=bf.setv(0,0,1);
      tp[7]=bf.minu(tp[0],tp[1]);
      tp[8]=bf.unit(tp[7]);
      bf.v2rot(tp[8],rmat);
      tp[9]=bf.rotv(tp[5],rmat);
      tp[10]=bf.scal(tp[9],tlength);
      tp[11]=bf.addv(tp[2],tp[10]);//for new end
      mcfragsweepCCD4(tmstr,numseq,sse[alphasind[trandpos]+1].init,
                      sse[alphasind[trandpos]+1].term, tp[2],
                      sse[alphasind[trandpos]+2].init,tp[11],
                      sse[alphasind[trandpos]+2].term, oldenergy,&tnewenergy);
      flagtor=geometry.tor2str(tmstr,numseq,3);
      if(!flagtor) printf("tor2str wrong in aaa\n");
      flagclash=false;
    }
    else if(attemptNum>=MAX_ITERATION)
    {
      summcaaa[2]++;
      return false;
    }
  } while(flagclash);

  if(tbeta<0) // skip energy calculation
  {
    summcaaa[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }
 
  *newenergy=energyFunction.calcTotalEnergy(tmstr,numseq);
  if(*newenergy<oldenergy)
  {
    summcaaa[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }
  else
  {
    double trand;
    trand=Random();
    if(trand<exp(-tbeta*(*newenergy-oldenergy)))
    {
      summcaaa[1]++;
      memcpy(decstr,tmstr,numseq*sizeof(point3f));
      return true;
    }
    else
    {
      summcaaa[2]++;
      return false;
    }
  }
}

/* used by mcfragsweepbtn */
bool FoldDesignMovement::mcfragsweepCCD5(
  point3f *decstr,
  int numseq,
  int lps,
  int lpe,
  point3d pt1,
  int ind1,
  point3d pt2,
  int ind2,
  double oldenergy,
  double *newenergy
)
{
  int trandpos;
  bool flagphi,flagtor;
  int flagpt,indp[2];
  indp[0]=ind1;indp[1]=ind2;
  point3d tp[2]; 
  tp[0].x=pt1.x;tp[0].y=pt1.y;tp[0].z=pt1.z;
  tp[1].x=pt2.x;tp[1].y=pt2.y;tp[1].z=pt2.z;
  point3d pcur,p12,p13;
  int numiter=0;
  double tdist[2],ttheta,tphi,tpsi,tinner;
  memcpy(tmstr2,decstr,numseq*sizeof(point3f));
  flagpt=0;
  int threshiter=100;
  do
  {
    tdist[0]=10000;
    tdist[1]=10000;
    trandpos=lps+int((lpe-lps+1)*Random());
    if(Random()<0.5) flagphi=true;//change phi n ca
    else flagphi=false;//change psi ca c in pos+1
    if(tmstr2[trandpos].ss3=='H' && tmstr2[trandpos+1].ss3=='H')
    {
      numiter++;
      continue;
    }
    else if(tmstr2[trandpos].ss3=='H') flagphi=false;
    else if(tmstr2[trandpos].ss3=='E' && tmstr2[trandpos+1].ss3=='E')
    {
      numiter++;
      continue;
    }
    else if(tmstr2[trandpos].ss3=='E') flagphi=false;
    flagpt=(flagpt+1)%2;
    if(flagphi)
    {
      p12=bf.setv(tmstr2[trandpos].x-tmstr2[trandpos].ptn.x,
                  tmstr2[trandpos].y-tmstr2[trandpos].ptn.y,
                  tmstr2[trandpos].z-tmstr2[trandpos].ptn.z);
      p13=bf.setv(tmstr2[indp[flagpt]].x-tmstr2[trandpos].ptn.x,
                  tmstr2[indp[flagpt]].y-tmstr2[trandpos].ptn.y,
                  tmstr2[indp[flagpt]].z-tmstr2[trandpos].ptn.z);
      tinner=bf.angv(p12,p13);

      //In one line no affect when rotate
      if(tinner<epsilon || PI-tinner<epsilon || bf.norm(p12)<epsilon)
      {
        numiter++;
        continue;
      }
      pcur.x=tmstr2[indp[flagpt]].x;
      pcur.y=tmstr2[indp[flagpt]].y;
      pcur.z=tmstr2[indp[flagpt]].z;
      p13.x=tmstr2[trandpos].x;
      p13.y=tmstr2[trandpos].y;
      p13.z=tmstr2[trandpos].z;
      p12.x=tmstr2[trandpos].ptn.x;
      p12.y=tmstr2[trandpos].ptn.y;
      p12.z=tmstr2[trandpos].ptn.z;
      ttheta=bf.phi(pcur.x,pcur.y,pcur.z,
                    p12.x,p12.y,p12.z,
                    p13.x,p13.y,p13.z,
                    tp[flagpt].x,tp[flagpt].y,tp[flagpt].z);
      tphi=tmstr2[trandpos].tor[2]+ttheta;
      if(tphi>360.0) tphi-=360.0;
      tmstr2[trandpos].tor[2]=tphi;
    }
    else
    {
      p12=bf.setv(tmstr2[trandpos].ptc.x-tmstr2[trandpos].x,
                  tmstr2[trandpos].ptc.y-tmstr2[trandpos].y,
                  tmstr2[trandpos].ptc.z-tmstr2[trandpos].z);
      p13=bf.setv(tmstr2[indp[flagpt]].x-tmstr2[trandpos].x,
                  tmstr2[indp[flagpt]].y-tmstr2[trandpos].y,
                  tmstr2[indp[flagpt]].z-tmstr2[trandpos].z);
      tinner=bf.angv(p12,p13);

      //In one line no affect when rotate
      if(tinner<epsilon || PI-tinner<epsilon || bf.norm(p12)<epsilon)
      {
        numiter++;
        continue;
      }
      pcur.x=tmstr2[indp[flagpt]].x;
      pcur.y=tmstr2[indp[flagpt]].y;
      pcur.z=tmstr2[indp[flagpt]].z;
      p12.x=tmstr2[trandpos].x;
      p12.y=tmstr2[trandpos].y;
      p12.z=tmstr2[trandpos].z;
      p13.x=tmstr2[trandpos].ptc.x;
      p13.y=tmstr2[trandpos].ptc.y;
      p13.z=tmstr2[trandpos].ptc.z;
      ttheta=bf.phi(pcur.x,pcur.y,pcur.z,
                    p12.x,p12.y,p12.z,
                    p13.x,p13.y,p13.z,
                    tp[flagpt].x,tp[flagpt].y,tp[flagpt].z);
      tpsi=tmstr2[trandpos+1].tor[0]+ttheta;
      if(tpsi>360.0) tpsi-=360.0;
      tmstr2[trandpos+1].tor[0]=tpsi;
    }
    numiter++;
    if(trandpos==0) flagtor=geometry.tor2str(tmstr2,indcut,3);
    else flagtor=geometry.tor2strp(tmstr2,indcut,trandpos);
    if(!flagtor) printf("tor2str wrong in CCD5 %d\n",trandpos);
    pcur.x=tmstr2[ind1].x-tp[0].x;
    pcur.y=tmstr2[ind1].y-tp[0].y;
    pcur.z=tmstr2[ind1].z-tp[0].z;
    tdist[0]=pcur.x*pcur.x+pcur.y*pcur.y+pcur.z*pcur.z;
    pcur.x=tmstr2[ind2].x-tp[1].x;
    pcur.y=tmstr2[ind2].y-tp[1].y;
    pcur.z=tmstr2[ind2].z-tp[1].z;
    tdist[1]=pcur.x*pcur.x+pcur.y*pcur.y+pcur.z*pcur.z;
  } while(numiter<threshiter && (tdist[0]>0.50 || tdist[1]>0.50));
  if(numiter<threshiter)
  {
    memcpy(decstr,tmstr2,numseq*sizeof(point3f));
    return true;
  }
  else
  {
    memcpy(decstr,tmstr2,numseq*sizeof(point3f));
    return false;
  }
}

void FoldDesignMovement::getfrombtntable(
  double tval,
  int *tlisteturn,
  double *tprobeturn,
  int *tindet,
  int *tindpos,
  int tnumeturn
)
{
  int i;
  int tpind;
  if(tnumeturn<2)
  {
    if(tval>tprobeturn[0]) printf("wrong %f %f\n",tval,tprobeturn[0]);
    *tindet=0;
    *tindpos=tlisteturn[0];
    return;
  }
  i=tnumeturn/2;
  if(tval>tprobeturn[i-1])
  {
    getfrombtntable(tval,tlisteturn+i,tprobeturn+i,&tpind,tindpos,tnumeturn-i);
    *tindet=tpind+i;
  }
  else getfrombtntable(tval,tlisteturn,tprobeturn,tindet,tindpos,i);
}

/* Movement M11 (beta-turn formation) tries to form a beta-turn motif for
 * every 4-mer segment along the query sequence. The number of M11
 * attempts at each position is proportional to the predicted beta-turn
 * probability from as implemented by getfrombtntable. */
bool FoldDesignMovement::mcfragsweepbtn(
  point3f *decstr,
  int numseq,
  double tbeta,
  double oldenergy,
  double *newenergy
)//beta turn
{
  double trand2;
  int i,j;
  int attemptNum=0;
  int indpos,indi;
  double rmat[9];
  float tleng1,tleng2,tangle1;
  point3d tp[18];
  point3s fp[4];
  double tnewenergy;
  bool flagclash=true;
  bool flagtor;
  double threshclash=0.01;
  if(numeturn==0)
  {
    summcbtn[2]++;
    return false;
  }
  do
  {
    flagclash=false;
    memcpy(tmstr,decstr,numseq*sizeof(point3f));
    trand2=Random();
    getfrombtntable(trand2,listeturn,probeturn,&indi,&indpos,numeturn);

    if(!flagclash) 
    {
      tp[0]=bf.setv(tmstr[indpos].x,tmstr[indpos].y,tmstr[indpos].z);
      tp[1]=bf.setv(tmstr[indpos+1].x,tmstr[indpos+1].y,tmstr[indpos+1].z);
      tp[2]=bf.setv(tmstr[indpos+2].x,tmstr[indpos+2].y,tmstr[indpos+2].z);
      tp[3]=bf.setv(tmstr[indpos+3].x,tmstr[indpos+3].y,tmstr[indpos+3].z);
      tp[4]=bf.minu(tp[0],tp[1]);
      tp[5]=bf.unit(tp[4]);
      bf.v2rot(tp[5],rmat);
      tp[6]=bf.rana(1.58825+0.175*Random()-0.03); // [0.44*pi,0.55*pi)
      tp[7]=bf.rotv(tp[6],rmat);
      tp[8]=bf.minu(tp[2],tp[1]);
      tleng1=bf.norm(tp[8]);
      tp[9]=bf.scal(tp[7],tleng1);
      tp[10]=bf.addv(tp[1],tp[9]);//new third
      tp[11]=bf.minu(tp[3],tp[2]);
      tleng2=bf.norm(tp[11]);
      tangle1=0.994838+0.436*Random()-0.35;
      fp[0].x=tp[0].x;fp[0].y=tp[0].y;fp[0].z=tp[0].z;
      fp[1].x=tp[1].x;fp[1].y=tp[1].y;fp[1].z=tp[1].z;
      fp[2].x=tp[10].x;fp[2].y=tp[10].y;fp[2].z=tp[10].z;
      bf.tor2pos22(fp[0].x,fp[0].y,fp[0].z,
                   fp[1].x,fp[1].y,fp[1].z,
                   tp[10].x,tp[10].y,tp[10].z,
                   tangle1,tleng2,1.65806f+0.12*Random()-0.02f,
                   &fp[3].x,&fp[3].y,&fp[3].z);
      tp[12]=bf.setv(fp[3].x,fp[3].y,fp[3].z);
      mcfragsweepCCD5(tmstr,numseq,indpos,indpos+3,
                      tp[10],indpos+2,tp[12],indpos+3,
                      oldenergy,&tnewenergy);
      flagtor=geometry.tor2str(tmstr,numseq,3);
      if(!flagtor) printf("tor2str wrong in btn\n");
      flagclash=false;
    }
    else if(attemptNum>=MAX_ITERATION)
    {
      summcbtn[2]++;
      return false;
    }
  } while(flagclash);

  if(tbeta<0) // skip energy calculation
  {
    summcbtn[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }
    
  *newenergy=energyFunction.calcTotalEnergy(tmstr,numseq);
  if(*newenergy<oldenergy)
  {
    summcbtn[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }
  else
  {
    double trand;
    trand=Random();
    if(trand<exp(-tbeta*(*newenergy-oldenergy)))
    {
      summcbtn[1]++;
      memcpy(decstr,tmstr,numseq*sizeof(point3f));
      return true;
    }
    else
    {
      summcbtn[2]++;
      return false;
    }
  }
}

/* used by mcfragsweepsft to enclose termini residues after segment
 * shifting. */
bool FoldDesignMovement::mcfragsweepLMP2(
  point3f *tmstr,
  int numseq,
  int poss,
  int pose
)
{   
  //[trandpos, pose)
  int trandpos=poss;
  int trand2=pose-poss;
  int numiter;
  int k;
  point3d rp;
  double rnorm;
  bool flagdone;
  numiter=0;
  do
  {
    for(k=trandpos;k<trandpos+trand2;k++)
    {
      singlemoveLMPf(tmstr,k,int(Random()*7.0));
      singlemoveLMPf(tmstr,k,int(Random()*7.0));
      singlemoveLMPf(tmstr,k,int(Random()*7.0));
      singlemoveLMPf(tmstr,k,int(Random()*7.0));
      singlemoveLMPf(tmstr,k,int(Random()*7.0));
      singlemoveLMPf(tmstr,k,int(Random()*7.0));
      singlemoveLMPf(tmstr,k,int(Random()*7.0));
    }
    ////inverse
    for(k=trandpos+trand2;k>trandpos;k--)
    {
      singlemoveLMPb(tmstr,k,int(Random()*7.0));
      singlemoveLMPb(tmstr,k,int(Random()*7.0));
      singlemoveLMPb(tmstr,k,int(Random()*7.0));
      singlemoveLMPb(tmstr,k,int(Random()*7.0));
      singlemoveLMPb(tmstr,k,int(Random()*7.0));
      singlemoveLMPb(tmstr,k,int(Random()*7.0));
      singlemoveLMPb(tmstr,k,int(Random()*7.0));
    }
        
    ////check
    flagdone=true;
    for(k=trandpos;k<=trandpos+trand2;k++)
    {
      //ca-1 ca
      rp=bf.setv(tmstr[k].x-tmstr[k-1].x,tmstr[k].y-tmstr[k-1].y,tmstr[k].z-tmstr[k-1].z);
      rnorm=bf.norm(rp);
      if(rnorm<lencaca-delcaca || rnorm>lencaca+delcaca)
      {
        flagdone=false;
        break;
      }
      //ca-1 n
      rp=bf.setv(tmstr[k].ptn.x-tmstr[k-1].x,
                 tmstr[k].ptn.y-tmstr[k-1].y,
                 tmstr[k].ptn.z-tmstr[k-1].z);
      rnorm=bf.norm(rp);
      if(rnorm<lencan1-delcan1 || rnorm>lencan1+delcan1)
      {
        flagdone=false;
        break;
      }
      //c-1 n
      rp=bf.setv(tmstr[k].ptn.x-tmstr[k-1].ptc.x,
                 tmstr[k].ptn.y-tmstr[k-1].ptc.y,
                 tmstr[k].ptn.z-tmstr[k-1].ptc.z);
      rnorm=bf.norm(rp);
      if(rnorm<lencn-delcn || rnorm>lencn+delcn)
      {
        flagdone=false;
        break;
      }
      //c-1 ca
      rp=bf.setv(tmstr[k].x-tmstr[k-1].ptc.x,
                 tmstr[k].y-tmstr[k-1].ptc.y,
                 tmstr[k].z-tmstr[k-1].ptc.z);
      rnorm=bf.norm(rp);
      if(rnorm<lencca1-delcca1 || rnorm>lencca1+delcca1)
      {
        flagdone=false;
        break;
      }
      if(k==trandpos+trand2) break;
      //n ca
      rp=bf.setv(tmstr[k].x-tmstr[k].ptn.x,
                 tmstr[k].y-tmstr[k].ptn.y,
                 tmstr[k].z-tmstr[k].ptn.z);
      rnorm=bf.norm(rp);
      if(rnorm<lennca-delnca || rnorm>lennca+delnca)
      {
        flagdone=false;
        break;
      }
      //n c
      rp=bf.setv(tmstr[k].ptc.x-tmstr[k].ptn.x,
                 tmstr[k].ptc.y-tmstr[k].ptn.y,
                 tmstr[k].ptc.z-tmstr[k].ptn.z);
      rnorm=bf.norm(rp);
      if(rnorm<lennc-delnc || rnorm>lennc+delnc)
      {
        flagdone=false;
        break;
      }
      //ca c
      rp=bf.setv(tmstr[k].ptc.x-tmstr[k].x,
                 tmstr[k].ptc.y-tmstr[k].y,
                 tmstr[k].ptc.z-tmstr[k].z);
      rnorm=bf.norm(rp);
      if(rnorm<lencac-delcac || rnorm>lencac+delcac)
      {
        flagdone=false;
        break;
      }
    }
    numiter++;
  } while(!flagdone && numiter<200);
    
  if(!flagdone) return false;
  else return true;
}

/* Movement M8 shifts the residue numbers in a segment forward or backward
 * by one residue, which means that the coordinates  of each residue are
 * copied from their preceded or followed residue in the segment. We then
 * need to delete the unused coordinates of one residue in one terminal and
 * insert new coordinates of another residue in the other terminal, as
 * implemented by mcfragsweepLMP2. This movement can easily adjust the
 * beta-pairing in two well-aligned beta-strands. */
bool FoldDesignMovement::mcfragsweepsft(
  point3f *decstr,
  int numseq,
  double tbeta,
  double oldenergy,
  double *newenergy
)
{
  int i,j;
  int attemptNum=0;
  int istart,iend;
  int trandpos,trandpos2;
  bool flagdirect;
  bool flagclash=true;
  bool flagtor;
  bool flagok1,flagok2;
  bool flagcsc;
  double threshclash=0.01;
  int numse;
  bool exchange_torsion=true; // exchange torsion-angle space coordinate
                              // in addition to cartesian coordinate
  int tor_dim;// index of torsion-angle space dimension
  do
  {
    flagclash=false;
    memcpy(tmstr,decstr,numseq*sizeof(point3f));
    //[trandpos, trandpos2] [3 numseq-4]
    numse=0;
    do
    {
      trandpos2=3+int((numseq-10)*Random());
      trandpos=3+int((numseq-trandpos2-6)*Random());
      trandpos2+=trandpos;
      numse++;
    } while((trandpos<3 || trandpos2>numseq-4 || 
            tmstr[trandpos].ssm!='C' || tmstr[trandpos2].ssm!='C') && numse<300);

    if(tmstr[trandpos].ssm!='C' || tmstr[trandpos2].ssm!='C')
    {
      summcsft[2]++;
      return false;
    }
            
    double treject=Random();
    if(((tmstr[trandpos-1].ssm!='C' || tmstr[trandpos+1].ssm!='C') && treject>0.1) ||
       ((tmstr[trandpos2-1].ssm!='C'|| tmstr[trandpos2+1].ssm!='C')&& treject>0.1))
    {
      summcsft[2]++;
      return false;
    }
    //all coil
    flagcsc=false;
    for(i=trandpos+1;i<trandpos2;i++)
    {
      if(tmstr[i].ssm!='C')
      {
        flagcsc=true;
        break;
      }
    }

    if(!flagclash)
    {
      if(Random()<0.5) flagdirect=0;//back
      else flagdirect=1;//forward

      if(flagdirect)
      {
        for(i=trandpos;i<=trandpos2;i++)
        {
          tmstr[i].x=tmstr[i+1].x;
          tmstr[i].y=tmstr[i+1].y;
          tmstr[i].z=tmstr[i+1].z;
          tmstr[i].ptn.x=tmstr[i+1].ptn.x;
          tmstr[i].ptn.y=tmstr[i+1].ptn.y;
          tmstr[i].ptn.z=tmstr[i+1].ptn.z;
          tmstr[i].ptc.x=tmstr[i+1].ptc.x;
          tmstr[i].ptc.y=tmstr[i+1].ptc.y;
          tmstr[i].ptc.z=tmstr[i+1].ptc.z;
          if(exchange_torsion)
          {
            for(tor_dim=0;tor_dim<3;tor_dim++)
            {
              tmstr[i].len[tor_dim]=tmstr[i+1].len[tor_dim];
              tmstr[i].ang[tor_dim]=tmstr[i+1].ang[tor_dim];
              tmstr[i].tor[tor_dim]=tmstr[i+1].tor[tor_dim];
            }
          }
        }
      }//gap in [trandpos-1,trandpos] dup in [trandpos2,trandpos2+1]
      else
      {
        for(i=trandpos2;i>=trandpos;i--)//for(i=trandpos2;i>=trandpos2;i--)
        {
          tmstr[i].x=tmstr[i-1].x;
          tmstr[i].y=tmstr[i-1].y;
          tmstr[i].z=tmstr[i-1].z;
          tmstr[i].ptn.x=tmstr[i-1].ptn.x;
          tmstr[i].ptn.y=tmstr[i-1].ptn.y;
          tmstr[i].ptn.z=tmstr[i-1].ptn.z;
          tmstr[i].ptc.x=tmstr[i-1].ptc.x;
          tmstr[i].ptc.y=tmstr[i-1].ptc.y;
          tmstr[i].ptc.z=tmstr[i-1].ptc.z;
          if(exchange_torsion)
          {
            for(tor_dim=0;tor_dim<3;tor_dim++)
            {  
              tmstr[i].len[tor_dim]=tmstr[i-1].len[tor_dim];
              tmstr[i].ang[tor_dim]=tmstr[i-1].ang[tor_dim];
              tmstr[i].tor[tor_dim]=tmstr[i-1].tor[tor_dim];
            }
          }
        }
      }//dup in [trandpos-1,trandpos] gap in [trandpos2,trandpos2+1]
      flagok1=mcfragsweepLMP2(tmstr,numseq,trandpos-2,trandpos+2);
      if(!flagok1)
      { 
        summcsft[2]++;
        return false;
      }
      flagok2=mcfragsweepLMP2(tmstr,numseq,trandpos2-1,trandpos2+3);
      if(!flagok2)
      {
        summcsft[2]++;
        return false;
      }
      istart=trandpos-2;
      iend=trandpos+2;
      if(iend>=numseq) iend=numseq-1;
      geometry.str2torp(tmstr,numseq,istart,iend);

      istart=trandpos2-1;
      iend=trandpos2+3;
      if(iend>=numseq) iend=numseq-1;
      geometry.str2torp(tmstr,numseq,istart,iend);

      istart=trandpos-2;
      if(istart==0 || istart==1) flagtor=geometry.tor2str(tmstr,numseq,3);
      else flagtor=geometry.tor2strp(tmstr,numseq,istart);

      if(!flagtor) printf("tor2str wrong in sft %d\n",istart);
      flagclash=false;
    }
    else if(attemptNum>=MAX_ITERATION)
    {
      summcsft[2]++;
      return false;
    }
  } while(flagclash);

  if(tbeta<0) // skip energy calculation
  {
    summcsft[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }

  *newenergy=energyFunction.calcTotalEnergy(tmstr,numseq);
  if(flagcaca && enelist[14]<enelistbk[14] && 
     (*newenergy<oldenergy || (*newenergy-oldenergy)<-0.01*oldenergy))
  {
    summcsft[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }
  if(*newenergy<oldenergy)
  {
    summcsft[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }
  else
  {
    double trand;
    trand=Random();
    if(trand<exp(-tbeta*(*newenergy-oldenergy)))
    {
      summcsft[1]++;
      memcpy(decstr,tmstr,numseq*sizeof(point3f));
      return true;
    }
    else
    {
      summcsft[2]++;
      return false;
    }
  }
}

/* Movement M12 randomly translates a segment,
 * and closes terminal gap by mcfragsweepLMP2 */
bool FoldDesignMovement::mcfragsweeptra(
  point3f *decstr,
  int numseq,
  double tbeta,
  double oldenergy,
  double *newenergy
)
{
  int i,j;
  int attemptNum=0;
  int istart,iend;
  int trandpos,trandpos2;
  bool flagclash=true;
  bool flagtor,flagok1,flagok2;
  double threshclash=0.01;
  point3d rp;
  double rmax;
  do
  {
    flagclash=false;
    memcpy(tmstr,decstr,numseq*sizeof(point3f));
    //[trandpos, trandpos2] [3 numseq-4]
    do
    {
      trandpos2=1+int(19*Random());//1-19
      trandpos=3+int((numseq-trandpos2-6)*Random());//3-l-del-4
      trandpos2+=trandpos;
    } while(trandpos<3 || trandpos2>numseq-4);

    if(!flagclash)
    {
      rmax=0.05+4.95*Random();
      rp=bf.ranv_ergodic(rmax); //rp=bf.ranv(rmax);
      for(i=trandpos;i<=trandpos2;i++)
      {
        tmstr[i].x+=rp.x;
        tmstr[i].y+=rp.y;
        tmstr[i].z+=rp.z;
        tmstr[i].ptn.x+=rp.x;
        tmstr[i].ptn.y+=rp.y;
        tmstr[i].ptn.z+=rp.z;
        tmstr[i].ptc.x+=rp.x;
        tmstr[i].ptc.y+=rp.y;
        tmstr[i].ptc.z+=rp.z;
      }
      flagok1=mcfragsweepLMP2(tmstr,numseq,trandpos-2,trandpos+2);
      if(!flagok1)
      {
        summctra[2]++;
        return false;
      }
      flagok2=mcfragsweepLMP2(tmstr,numseq,trandpos2-1,trandpos2+3);
      if(!flagok2)
      {
        summctra[2]++;
        return false;
      }
      istart=trandpos-2;
      iend=trandpos+2;
      if(iend>=numseq) iend=numseq-1;
      geometry.str2torp(tmstr,numseq,istart,iend);
        
      istart=trandpos2-1;
      iend=trandpos2+3;
      if(iend>=numseq) iend=numseq-1;
      geometry.str2torp(tmstr,numseq,istart,iend);
      istart=trandpos-2;
      if(istart==0 || istart==1) flagtor=geometry.tor2str(tmstr,numseq,3);
      else flagtor=geometry.tor2strp(tmstr,numseq,istart);
      if(!flagtor) printf("tor2str wrong in tra %d\n",istart);
      flagclash=false;
    }
    else if(attemptNum>=MAX_ITERATION)
    {
      summctra[2]++;
      return false;
    }
  } while(flagclash);

  if(tbeta<0) // skip energy calculation
  {
    summctra[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }

  *newenergy=energyFunction.calcTotalEnergy(tmstr,numseq);
  if(flagcaca && enelist[14]<enelistbk[14] && 
     (*newenergy<oldenergy || (*newenergy-oldenergy)<-0.01*oldenergy))
  {
    summctra[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }
  if(*newenergy<oldenergy)
  {
    summctra[0]++;
    memcpy(decstr,tmstr,numseq*sizeof(point3f));
    return true;
  }
  else
  {
    double trand;
    trand=Random();
    if(trand<exp(-tbeta*(*newenergy-oldenergy)))
    {
      summctra[1]++;
      memcpy(decstr,tmstr,numseq*sizeof(point3f));
      return true;
    }
    else
    {
      summctra[2]++;
      return false;
    }
  }
}

/* detemine the movement type according to the given Cummulative
 * Density Function tmov[]. tmov[i] is the probability of selecting
 * movement type < i */
int FoldDesignMovement::getmovetype(
  double tmov[],
  int totmov,
  double trandnum
)
{
  //[0      a      b      c      d     1]
  //substi rotrem tranrem rotmid  extang
  //double tmov[6]={0.0,0.5,0.8,0.85,0.90,1.0001};

  for(int i=0;i<totmov;i++)
    if(trandnum>=tmov[i] && trandnum<tmov[i+1])
      return i;
  return 0;
}

/* entry function for movement and metropolis */
bool FoldDesignMovement::attemptConformationalMove(
  point3f *decstr,
  int numseq,
  double tbeta,
  double oldenergy,
  double *newenergy,
  int *movtype
)
{
  const int ALPHA_PROT=1; 
  const int NOT_ALPHA_PROT=2;
  int totmov=14;
  double *tmov;

  //i_move:            0   1    2    3    4    5    6    7    8    9    10   11   12   13   14 
  //                  sub phi  psi   dh  rot  lmp  bbp  ome  len  ang  aaa  btn  sft  tra
  //speed              1   1    1    1    3    5    4    1    1    1    5    5    3    2
  double tmovab[15] ={0.0,0.34,0.42,0.50,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.92,0.95,1.01,1.0001}; //QE
  double tmovtpb[15]={0.0,0.15,0.17,0.20,0.25,0.35,0.50,0.50,0.55,0.70,0.80,0.80,0.80,0.85,1.0001}; //QT for beta
  double tmovtpa[15]={0.0,0.35,0.37,0.40,0.45,0.55,0.65,0.65,0.70,0.80,0.85,0.85,0.85,0.90,1.0001}; //QT for alpha
  //double tmovrf[14]={0.0,0.10,0.25,0.40,0.50,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.0001}; //combo/model/template
  
  double movetype=Random();
  if(flagLocalMove && inputInfo.getProteinType()==NOT_ALPHA_PROT)
  {
    tmov=tmovtpb; //Local movements for beta proteins
  }    
  else if(flagLocalMove && inputInfo.getProteinType()==ALPHA_PROT)
  {
    tmov=tmovtpa; //Local movements for alpha proteins
  }  
  else
  {
    tmov=tmovab; //Standard move set
  }

  int ntype=getmovetype(tmov,totmov,movetype);
  bool flagAccept; //Whether movement was accepted
  switch(ntype)
  {
    case 0:   // M5: fragment substitution
      flagAccept=mcfragsweepfragsub(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0.34, 0.25, 0.35
    case 1:   // M3: change backbone torsion angle phi
      flagAccept=mcfragsweepphi(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0.08, 0.02, 0.02
    case 2:   // M3: change backbone torsion angle psi
      flagAccept=mcfragsweeppsi(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0.08, 0.03, 0.03
    case 3:   // M4: substitute residue geometry by that in topdh
      flagAccept=mcfragsweeptopdh(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0.10, 0.05, 0.05
    case 4:   // M7: DEMO loop modeling (rotate around loop termini CA)
      flagAccept=mcfragsweepmid(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0.05, 0.10, 0.10
    case 5:   // M6: LMProt perturbation
      flagAccept=mcfragsweepLMP(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0.05, 0.15, 0.10
    case 6:   // M10: beta pairing
      flagAccept=mcfragsweepbbp(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0.05, 0, 0
    case 7:   // M3: change backbone torsion angle omega
      flagAccept=mcfragsweepome(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0.05, 0.05, 0.05 
    case 8:   // M1: change bond length
      flagAccept=mcfragsweeplen(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0.05, 0.15, 0.10
    case 9:   // M2: change bond angle
      flagAccept=mcfragsweepang(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0.05, 0.10, 0.05
    case 10:  // M9: helix packing
      flagAccept=mcfragsweepaaa(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0.02, 0, 0
    case 11:  // M11: beta turn formation 
      flagAccept=mcfragsweepbtn(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0.03, 0, 0
    case 12:  // M8: shift segment by one residue
      flagAccept=mcfragsweepsft(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0.05, 0.05, 0.05
    case 13:  // M12: translate a segment
      flagAccept=mcfragsweeptra(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0, 0.15, 0.1
    default:
      printf("wrong move type %d\n",ntype);
      break; 
  }

  *movtype=ntype;
  summctot[ntype]++;
  i_move_type=ntype;
  
  int i;
  if(flagAccept) // flagAccept=true: accepted move; false: rejected
  {         
    for(i=0;i<20;i++)
      enelistbk[i]=enelist[i];
  }
 
  return flagAccept;
}

/* movement without energy: always accept any movetype */
bool FoldDesignMovement::mcfragsweepcom_noene(
  point3f *decstr,
  int numseq,
  int ntype
)
{
  /* dummy variables */
  double tbeta=-1; // skip energy calculation if temperature is negative
  if(ntype==-1) // mcfragsweepsub double anchored version
  {
    tbeta=-2;
    ntype=0;
  }
  double oldenergy=0;
  double *newenergy;

  bool flagres; // whether movement is accepted
  switch (ntype)
  {
    case 0:   // M5: fragment substitution
      flagres=mcfragsweepfragsub(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0.34, 0.25, 0.35
    case 1:   // M3: change backbone torsion angle phi
      flagres=mcfragsweepphi(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0.08, 0.02, 0.02
    case 2:   // M3: change backbone torsion angle psi
      flagres=mcfragsweeppsi(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0.08, 0.03, 0.03
    case 3:   // M4: substitute residue geometry by that in topdh
      flagres=mcfragsweeptopdh(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0.10, 0.05, 0.05
    case 4:   // M7: DEMO loop modeling (rotate around loop termini CA)
      flagres=mcfragsweepmid(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0.05, 0.10, 0.10
    case 5:   // M6: LMProt perturbation
      flagres=mcfragsweepLMP(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0.05, 0.15, 0.10
    case 6:   // M10: beta pairing
      flagres=mcfragsweepbbp(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0.05, 0, 0
    case 7:   // M3: change backbone torsion angle omega
      flagres=mcfragsweepome(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0.05, 0.05, 0.05 
    case 8:   // M1: change bond length
      flagres=mcfragsweeplen(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0.05, 0.15, 0.10
    case 9:   // M2: change bond angle
      flagres=mcfragsweepang(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0.05, 0.10, 0.05
    case 10:  // M9: helix packing
      flagres=mcfragsweepaaa(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0.02, 0, 0
    case 11:  // M11: beta turn formation 
      flagres=mcfragsweepbtn(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0.03, 0, 0
    case 12:  // M8: shift segment by one residue
      flagres=mcfragsweepsft(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0.05, 0.05, 0.05
    case 13:  // M12: translate a segment
      flagres=mcfragsweeptra(decstr,numseq,tbeta,oldenergy,newenergy);
      break;// tmov = 0, 0.15, 0.1
    default:
      printf("wrong move type %d\n",ntype);
      break; 
  }
  return flagres;
}

void FoldDesignMovement::setpara8( //set temperature for QE/QN/QP/QT
  int numseq
) 
{
  int i,j;
  int N_rep;
  double t1,t2,tt1,tt2,deltai,deltat,tmp;
  double sft=5.0,Lch100=100;

  //printf("n_rep=%d\n",n_rep);
  
  N_rep=n_rep;
  t1=T1; //lowest temperature
  t2=T2; //highest temperature
  
  // T1=10/15;
  
  // -------- vector settings --------------->
  if(temparray)
  {
    delete[]temparray;
    temparray=NULL;
  }
  temparray=new double[N_rep];
  if(ttemparray)
  {
    delete[]ttemparray;
    ttemparray=NULL;
  }
  ttemparray=new double[N_rep];
  if(tttemparray)
  {
    delete[]tttemparray;
    tttemparray=NULL;
  }
  tttemparray=new double[N_rep];
  
  //------- re-scale temperature by length ------------>
  tt1=1+(numseq-Lch100)*S1; //length factor for t1
  tt2=1+(numseq-Lch100)*S2; //length factor for t2
  if(tt1<0.5) tt1=0.5;
  if(tt2>2.5) tt2=2.5;
  t1*=tt1;
  t2*=tt2;
  
  //------- set-up temperature for all replicas ------------>
  deltai=(1.0/t1-1.0/t2)/double(N_rep-1)*2.0/(1.0+sft);
  deltat=(sft-1.0)*deltai/double(N_rep-2);
  printf("t1,t2,deltai,deltat= %f, %f, %f, %f\n",t1,t2,deltai,deltat);

  temparray[N_rep-1]=t2;
  for(i=0;i<N_rep-1;i++)
  {
    temparray[N_rep-2-i]=1.0/(deltai+1.0/temparray[N_rep-1-i]+deltat*i);
    //printf("%5d %5d %8.3f %8.3f\n",i,N_rep-2-i,temparray[N_rep-2-i],temparray[N_rep-1-i]);
    //temparray[i]=t1*pow(t2/t1,float(i)/(N_rep-1));
    //aT_rep(i)=atemp1*(atemp2/atemp1)**(float(i-1)/(N_rep-1))  I-TASSER
  }
  
  printf("4, T_min=%8.3f(%8.3f), T_max=%8.3f(%8.3f)\n",t1,T1,t2,T2);
  for(i=0;i<N_rep;i++) printf("i_rep=%5d, T_i=%8.5f\n",i,temparray[i]);
}

void FoldDesignMovement::settmstr(
  int seqLength
)
{
  if(tmstr)
  {
    delete[]tmstr;
    tmstr=NULL;
  }
  tmstr=new point3f[seqLength];
  if(tmstr2)
  {
    delete[]tmstr2;
    tmstr2=NULL;
  }
  tmstr2=new point3f[seqLength];
  if(fixed_decstr)
  {
    delete[]fixed_decstr;
    fixed_decstr=NULL;
  }
  fixed_decstr=new point3f[seqLength];
}

void FoldDesignMovement::calcabind()
{
  int i;
  sssegment *sse=inputInfo.getSSE();
  int numSSE=inputInfo.getNumSSE();

  if(alphaind)
  {
    delete[]alphaind;
    alphaind=NULL;
  }
  if(betaind)
  {
    delete[]betaind;
    betaind=NULL;
  }
  if(alphasind)
  {
    delete[]alphasind;
    alphasind=NULL;
  }
  numsalpha=0;
  numalpha=0;
  numbeta=0;
  alphaind=new int[numSSE];
  alphasind=new int[numSSE];
  betaind=new int[numSSE];
  for(i=0;i<numSSE;i++)
  {
    if(sse[i].ss=='H') alphaind[numalpha++]=i;
    else if(sse[i].ss=='E') betaind[numbeta++]=i;
  }
  for(i=0;i<numalpha-1;i++)
  {
    if(alphaind[i+1]==alphaind[i]+2) alphasind[numsalpha++]=alphaind[i];
  }
}

void FoldDesignMovement::calcbbptable(
  int numseq
)
{
  int i,j,ii,jj;
  sssegment *sse=inputInfo.getSSE();
  int numSSE=inputInfo.getNumSSE();

  numbbptb=0;
  if(bbptb)
  {
    delete[]bbptb;
    bbptb=NULL;
  }
  numbbptb=numseq*numseq;
  bbptb=new paircont[numbbptb];
  for(i=0;i<numbbptb;i++)
  {
    bbptb[i].init=-1;
    bbptb[i].term=-1;
    bbptb[i].dist=0;
  }
  //depth for beta
  double *betadep=new double[numseq];
  for(i=0;i<numseq;i++) betadep[i]=0;
  for(i=0;i<numSSE;i++)
  {
    if(sse[i].ss=='C')
      for(j=sse[i].init;j<=sse[i].term;j++) betadep[j]=0.6;
    else if(sse[i].ss=='H')
      for(j=sse[i].init;j<=sse[i].term;j++) betadep[j]=0.0;//0.2
  }
  for(i=0;i<numbeta;i++)
  {
    double tlen=(sse[betaind[i]].term+1-sse[betaind[i]].init)/2.0;
    double tmid=(sse[betaind[i]].term+sse[betaind[i]].init)/2.0;
    for(j=sse[betaind[i]].init;j<=sse[betaind[i]].term;j++)
      betadep[j]=tlen-fabs(tmid-j)+6.5;
  }
  //start
  int totnum=0;
  for(i=0;i<numSSE;i++)
  {
    if(sse[i].ss=='H') continue;
    for(ii=sse[i].init;ii<=sse[i].term;ii++)
    {
      for(j=0;j<numSSE;j++)
      {
        if(sse[j].ss=='H') continue;
        for(jj=sse[j].init;jj<=sse[j].term;jj++)
        {
          bbptb[totnum].init=ii;
          bbptb[totnum].term=jj;
          if(fabs(ii-jj)<3)//many ==3 are turns
            bbptb[totnum].dist=0;
          else if(i==j && sse[i].ss=='H')
            bbptb[totnum].dist=0.0;//0.1
          else if(i==j && sse[i].ss=='E')
            bbptb[totnum].dist=0.7;
          else if(i==j && sse[i].ss=='C')
            bbptb[totnum].dist=0.4;
          else if(abs(i-j)==1 && sse[j].ss=='E' && sse[i].ss=='C')
            bbptb[totnum].dist=0.3;
          else if(abs(i-j)==1 && sse[j].ss=='E' && sse[i].ss=='H')
            bbptb[totnum].dist=0.2;
          else if(sse[i].ss=='E' && abs(i-j)==1 && sse[j].ss=='C')
            bbptb[totnum].dist=0.3;
          else if(sse[i].ss=='E' && abs(i-j)==1 && sse[j].ss=='H')
            bbptb[totnum].dist=0.2;
          else if(sse[i].ss=='C' && sse[j].ss=='C')
            bbptb[totnum].dist=0.6;
          else if(sse[i].ss=='C' && sse[j].ss=='H')
            bbptb[totnum].dist=0.2;
          else if(sse[i].ss=='H' && sse[j].ss=='C')
            bbptb[totnum].dist=0.2;
          else if(sse[i].ss=='H' && sse[j].ss=='H')
            bbptb[totnum].dist=0.1;
          else if(sse[i].ss=='H' && sse[j].ss=='E')
            bbptb[totnum].dist=0.3;
          else if(sse[i].ss=='E' && sse[j].ss=='H')
            bbptb[totnum].dist=0.3;
          else if(sse[i].ss=='C' && sse[j].ss=='E')
            bbptb[totnum].dist=(betadep[ii]+betadep[jj])/2.0;
          else if(sse[i].ss=='E' && sse[j].ss=='C')
            bbptb[totnum].dist=(betadep[ii]+betadep[jj])/2.0;
          else if(sse[i].ss=='E' && sse[j].ss=='E')
          {
            bbptb[totnum].dist=betadep[ii]+betadep[jj]+4.0*exp(-fabs(betadep[ii]-betadep[jj]));
            if(abs(i-j)==2) bbptb[totnum].dist+=5.0;
          }
          else
            bbptb[totnum].dist=betadep[ii]+betadep[jj];
          bbptb[totnum].dist*=bborder[abs(bbptb[totnum].init-bbptb[totnum].term)];
          totnum++;
        }
      }
    }
  }
  delete[]betadep;
  bool *flaguse=new bool[numbbptb];
  for(i=0;i<numbbptb;i++)
  {
    flaguse[i]=true;
    if(bbptb[i].dist<epsilon || bbptb[i].init>bbptb[i].term || 
       bbptb[i].init==0 || bbptb[i].term==numseq-1)
    {
      flaguse[i]=false;
    }
  }
  j=0;
  for(i=0;i<numbbptb;i++)
  {
    if(flaguse[i])
    {
      bbptb[j].init=bbptb[i].init;
      bbptb[j].term=bbptb[i].term;
      bbptb[j].dist=bbptb[i].dist;
      j++;
    }
  }
  numbbptb=j;
  delete[]flaguse;
  double totsum=0;
  for(i=0;i<numbbptb;i++) totsum+=bbptb[i].dist;
  for(i=0;i<numbbptb;i++)
  {
    bbptb[i].dist/=totsum;
    bbptb[i].prob=bbptb[i].dist;
  }
  for(i=1;i<numbbptb;i++) bbptb[i].dist+=bbptb[i-1].dist;
}

void FoldDesignMovement::setInitDecoy(
  point3f *decstr,
  int numseq
)//not that flexible
{
  int i,k,j;
  int trand2;
  int istart,iend;
  bool flagtor;
  ParsePDB pp;
  PrintFunc pf;
  int seglength=9;//nosegunit 12 for easy 9 for hard
  double *pin=new double[2*seglength];
  double *pout=new double[2*seglength];

  //first fragment
  trand2=int((topno)*Random());
  for(k=0;k<seglength;k++)
    decstr[k]=fragcont[seglength-1][trand2*seglength+k][0];
  iend=seglength-1;
  int tpend,addlen;
  double minidist,tdist;
  int miniind;
  int ntop=ntopfrag;//1 is better than 5 for abinitio, 5 is better using init.dat
  if(ntop<1) ntop=1;
  else if(ntop>50) ntop=50;
  double mindist[topno];
  int minind[topno];
  while(iend<numseq-1)
  {
    //add length [1 nosegunit-3]
    tpend=iend;
    addlen=1+int((seglength-3)*Random());//[istart,iend] [istart tpend]
    iend=iend+addlen;
    if(iend>=numseq)
    {
      addlen=numseq-1-iend;
      iend=numseq-1;
    }
    istart=iend-seglength+1;
    for(k=0;k<tpend-istart+1;k++)
    {
      pin[k]=decstr[istart+k].tor[0];
      pin[k+tpend-istart+1]=decstr[istart+k].tor[2];
    }
    minidist=1000000;
    for(k=0;k<ntop;k++) mindist[k]=1000000000;
    for(i=0;i<topno;i++)
    {
      for(k=0;k<tpend-istart+1;k++)
      {
        pout[k]=fragcont[seglength-1][i*seglength+k][istart].tor[0];
        pout[k+tpend-istart+1]=fragcont[seglength-1][i*seglength+k][istart].tor[2];
      }
      tdist=pp.lsfrmsdp(pin,pout,tpend-istart+1,360);
      if(tdist<minidist)
      {
        minidist=tdist;
        miniind=i;
      }
      for(k=0;k<ntop;k++)
      {
        if(tdist<mindist[k]) break;
      }
      if(k<ntop)//insert at pos k
      {
        for(j=ntop-2;j>=k;j--)
        {
          mindist[j+1]=mindist[j];
          minind[j+1]=minind[j];
        }
          mindist[k]=tdist;
          minind[k]=i;
      }
    }//i
    miniind=minind[int(Random()*ntop)];//5 for hard
    for(k=0;k<iend-istart+1;k++)
      decstr[istart+k]=fragcont[seglength-1][seglength*miniind+k][istart];
  }
  flagtor=geometry.tor2str(decstr,numseq,3);
  if(!flagtor && flagVerbose) printf("tor2str wrong in init decoy\n");
  delete[]pin;
  delete[]pout;
  pin=NULL;
  pout=NULL;
}

void FoldDesignMovement::setInitDecoyFixed(
  point3f *decstr,
  int numseq
)
{
  int i,j;
  int start_ind,end_ind;
  bool flagtor;
  vector<int> start,end;
  start.push_back(0);
  end.push_back(fixed_pos[0]);
  start_ind=fixed_pos[0];
  for(i=0;i<fixed_pos.size();i++){
    if(fixed_pos[i]==start_ind){
      start_ind++;
    }
    else{
      if(i==fixed_pos.size()-1){
        start.push_back(start_ind);
        end.push_back(numseq);
      }
      else{
        start.push_back(start_ind);
        end.push_back(fixed_pos[i]);
        start_ind=fixed_pos[i+1];
      }
    }
  }

  start_fixed.push_back(fixed_pos[0]);
  start_ind=fixed_pos[0];


  for(i=0;i<fixed_pos.size();i++){
    if(fixed_pos[i]==start_ind){
      start_ind++;
      if(i>=fixed_pos.size()-1){
        end_fixed.push_back(fixed_pos[i]);
      }

    }
    else{
      if(i>=fixed_pos.size()-1){
        end_fixed.push_back(fixed_pos[i]);
      }
      else{
        start_ind--;
        end_fixed.push_back(start_ind);
        start_fixed.push_back(fixed_pos[i]);
        start_ind=fixed_pos[i+1];
      }
    }
  }

  start.push_back(fixed_pos[fixed_pos.size()-1]+1);
  end.push_back(numseq);

  point3d rp;
  double rnorm;
  bool last=false;

  memcpy(fixed_decstr,decstr,numseq*sizeof(point3f));
  for(i=0;i<start.size();i++){
    fixed_decstr[start[i]].len[0]=float(lencn);
    if(end[i]<numseq-1){
      fixed_decstr[end[i]].len[0]=float(lencn);
    }
    if(start[i]==0) flagtor=geometry.tor2str(fixed_decstr,numseq,3);
    else flagtor=geometry.tor2strp(fixed_decstr,numseq,start[i]);
  }
  memcpy(decstr,fixed_decstr,numseq*sizeof(point3f));
}

