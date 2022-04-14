///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "FoldDesignEnergyFunction.h"

FoldDesignEnergyFunction::FoldDesignEnergyFunction()
{
  //---Ramachandran energy tables--->
  H_general_rama=NULL;
  E_general_rama=NULL;
  C_general_rama=NULL;
  H_specific_rama=NULL;
  E_specific_rama=NULL;
  C_specific_rama=NULL;
  G_specific_rama=NULL;
  I_specific_rama=NULL;
  B_specific_rama=NULL;
  S_specific_rama=NULL;
  T_specific_rama=NULL;

  //-----SSE packing energy tables---->
  energyHHPackAngle2=NULL;
  energyHHPackAngle3=NULL;
  energyHHPackAngle4=NULL;
  energyHHPackDist2=NULL;
  energyHHPackDist3=NULL;
  energySSPackAngle2=NULL;
  energySSPackAngle3=NULL;
  energySSPackDist1=NULL;
  energySSPackDist2=NULL;
  energySSPackDist3=NULL;
  energySSPackDist4=NULL;
  energyHSPackAngle2=NULL;
  energyHSPackAngle3=NULL;
  energyHSPackAngle4=NULL;
  energyHSPackDist2=NULL;
  energyHSPackDist3=NULL;
  energyHSPackDist4=NULL;

  kbp=NULL;
  solweight=NULL;
  contactConstr=NULL;
  distanceConstr=NULL;
  distanceWeight=NULL;
  distRestrType=NULL;  
  hbondRamp=false;
  paa=NULL;
}

FoldDesignEnergyFunction::~FoldDesignEnergyFunction()
{
}

void FoldDesignEnergyFunction::setEnergyFunctionParameters(
  const ParseFoldDesignEnergyFiles &energyParameters,
  GeometryCalc inputGeometry,
  InputData input
)
{

  int i;
  geometry=inputGeometry;
  inputInfo=input;

  int seqLength=inputInfo.getSeqLength();

  H_general_rama=energyParameters.H_general_rama;
  E_general_rama=energyParameters.E_general_rama;
  C_general_rama=energyParameters.C_general_rama;
  H_specific_rama=energyParameters.H_specific_rama;
  E_specific_rama=energyParameters.E_specific_rama;
  C_specific_rama=energyParameters.C_specific_rama;
  G_specific_rama=energyParameters.G_specific_rama;
  I_specific_rama=energyParameters.I_specific_rama;
  B_specific_rama=energyParameters.B_specific_rama;
  S_specific_rama=energyParameters.S_specific_rama;
  T_specific_rama=energyParameters.T_specific_rama;

  energyHHPackAngle2=energyParameters.energyHHPackAngle2;
  energyHHPackAngle3=energyParameters.energyHHPackAngle3;
  energyHHPackAngle4=energyParameters.energyHHPackAngle4;
  energyHHPackDist2=energyParameters.energyHHPackDist2;
  energyHHPackDist3=energyParameters.energyHHPackDist3;
  energySSPackAngle2=energyParameters.energySSPackAngle2;
  energySSPackAngle3=energyParameters.energySSPackAngle3;
  energySSPackDist1=energyParameters.energySSPackDist1;
  energySSPackDist2=energyParameters.energySSPackDist2;
  energySSPackDist3=energyParameters.energySSPackDist3;
  energySSPackDist4=energyParameters.energySSPackDist4;
  energyHSPackAngle2=energyParameters.energyHSPackAngle2;
  energyHSPackAngle3=energyParameters.energyHSPackAngle3;
  energyHSPackAngle4=energyParameters.energyHSPackAngle4;
  energyHSPackDist2=energyParameters.energyHSPackDist2;
  energyHSPackDist3=energyParameters.energyHSPackDist3;
  energyHSPackDist4=energyParameters.energyHSPackDist4;

  kbp=energyParameters.kbp;
  solweight=energyParameters.solweight;
  contactConstr=energyParameters.contactConstr;
  distanceConstr=energyParameters.distanceConstr;
  distanceWeight=energyParameters.distanceWeight;
  distRestrType=energyParameters.distRestrType;
  paa=energyParameters.paa;
  for(i=0;i<MAX_ENERGY_TERM_WEIGHTS;i++)
  {
    weights[i]=1.0;//energyParameters.weights[i];
  }

  setContactParameters(seqLength);
  estimateContactNum(seqLength);
  setLongestH();
}

void FoldDesignEnergyFunction::setLongestH()
{
  int numLong=0;
  int ti;
  int i;
  int numSSE=inputInfo.getNumSSE();
  sssegment *sse=inputInfo.getSSE();

  for(i=0;i<numSSE;i++)
  {
    if(sse[i].ss=='H')
    {
      ti=sse[i].term-sse[i].init+1;
      if(ti>numLong)
      {
        numLong=ti;
      }
    }
  }

  longestdist=2.30431+1.42339*numLong;
}

void FoldDesignEnergyFunction::estimateContactNum(
  int seqLength
)
{
  int n1=6;
  int n2=11;
  int Na=0;
  int Nb=0;
  char protClass;
  for(int i=0;i<seqLength;i++)
  {
    if(inputInfo.get3StateSSatPos(i)=='H')
    {
      Na++;
    }
    else if(inputInfo.get3StateSSatPos(i)=='E')
    {
      Nb++;
    }
  }

  if(Na<=n1+4 && Nb>=n2)
  {
    protClass='b';
  }
  else if(Na>=n2 && Nb<=n1-1)
  {
    protClass='a';
  }
  else if(Na>=n2 && Nb>=n2+1)
  {
    protClass='c';
  }
  else
  {
    protClass='l';
  }
  
  double stddevShort,stddevMed,stddevLong;
  if(seqLength<100)
  {
    if(protClass=='b')
    {
      expectedShort=int(0.258580841917*seqLength+14.9753137724);
      expectedMed=int(0.512377321092*seqLength+14.7453150544);
      expectedLong=int(2.21481698903*seqLength-96.9385091355);
      stddevShort=23.7655931423;
      stddevMed=41.3493136278;
      stddevLong=60.612299874;
    }
    else if(protClass=='a')
    {
      expectedShort=int(0.236560401413*seqLength-2.26079303056);
      expectedMed=int(0.23523141911*seqLength-1.05504499874);
      expectedLong=int(1.03708054813*seqLength-42.1633575447);
      stddevShort=16.7518088268;
      stddevMed=20.2354692626;
      stddevLong=39.8141230281;
    }
    else if(protClass=='c')
    {
      expectedShort=int(0.249217831889*seqLength+8.34788441399);
      expectedMed=int(0.392639947404*seqLength+7.52097264883);
      expectedLong=int(1.50449621411*seqLength-47.4875211876);
      stddevShort=17.8238513054;
      stddevMed=32.5036409826;
      stddevLong=52.8871710857;
    }
    else if(protClass=='l')
    {
      expectedShort=int(0.156887830616*seqLength+13.5362211947);
      expectedMed=int(0.168437496527*seqLength+16.8118191804);
      expectedLong=int(1.07204106063*seqLength-34.8330790041);
      stddevShort=19.0297514291;
      stddevMed=33.478349542;
      stddevLong=48.8153461264;
    }
  }
  else
  {
    if(protClass=='b')
    {
      expectedShort=int(0.437493981503*seqLength-6.49025152281);
      expectedMed=int(0.76103436954*seqLength-13.0626314295);
      expectedLong=int(1.74772186844*seqLength-44.9435731352);
      stddevShort=44.4661039753;
      stddevMed=101.11405066;
      stddevLong=124.145576789;
    }
    if(protClass=='a')
    {
      expectedShort=int(0.17124418185*seqLength+3.54989956225);
      expectedMed=int(0.208597610888*seqLength+1.17067837419);
      expectedLong=int(1.25943311791*seqLength-65.2919704946);
      stddevShort=25.2012746685;
      stddevMed=35.210788381;
      stddevLong=89.3872824254;
    }
    if(protClass=='c')
    {
      expectedShort=int(0.267713145805*seqLength+6.32635741825);
      expectedMed=int(0.344876755056*seqLength+12.8738023488);
      expectedLong=int(1.83336461092*seqLength-83.6660696009);
      stddevShort=31.4555693605;
      stddevMed=67.3124557517;
      stddevLong=113.277829653;
    }
    if(protClass=='l')
    {
      expectedShort=int(0.168881494424*seqLength+14.6050424764);
      expectedMed=int(0.154981969163*seqLength+18.1349399225);
      expectedLong=int(1.26918477813*seqLength-57.4085271214);
      stddevShort=26.6455705939;
      stddevMed=35.3825850796;
      stddevLong=95.5240211155;
    }
  }
}

void FoldDesignEnergyFunction::setContactParameters(
  int seqLength
)
{
  int dwell;
 
  if(seqLength<100)
  {
    dwell=6;
  }
  else if(seqLength<120)
  {
    dwell=8;
  }
  else if(seqLength<200)
  {
    dwell=10;
  }
  else
  {
    dwell=12;
  }

  d8=8;
  d10=d8+dwell;
  da=(d8+d10)/2; //middle of d8 and d10
  db=dwell;      //width of first well
  dc=(d10+80)/2; //middle of d10 and 80
  dd=80-d10;     //width of second well
}

void FoldDesignEnergyFunction::calcAllAtomEnergy(
  point3f *decstr,
  int numseq,
  double *enelist
)
{
  const int NUM_ATOMS=8;
  const int CA_INDEX=0,CB_INDEX=4,H_INDEX=5,G_INDEX=6;
  const int GLY_INDEX=5;
  int protType=inputInfo.getProteinType();
  int i,j,k,ii,jj,index;
  int tbin,ind1,ind2,indm,ind5[5],ind6[5];
  int atomid[]={1,17,0,26,2};
  int sgmap[]={2,1,3,4,99,99,99,0};
  point3d pin[8],pjn[8],tp;
  double tdist,tdist2,tdist3,r1,r2;
  double distcut=9.2;
  double sqdistcut=9.2*9.2;
  double scalf=0.00707;
  double cutca=4.0;
  double sqdist,tdists;
  double err;
  float distance,weightD;
  int restrType;
  double *soldat;
  soldat=new double[numseq];
  for(i=0;i<numseq;i++) soldat[i]=0;
  double weight=0,gg=0,gg1=0,weight_s,weight_m,weight_l;
  double b;
  int mk=0;

  int sdist=17;
  int maxdist;
  point3d tcen;
  tcen.x=0;tcen.y=0;tcen.z=0;
  enelist[0]=0;
  enelist[1]=0;  //Excluded volume, ww1
  enelist[4]=0;  //Contact restraints
  enelist[5]=0;  //not used generic side-chain contact potential, 'sgpolarity5.txt', ww5
  enelist[6]=0;  //not used sequence-specific solvation potential, 'sol.txt', ww6
  enelist[7]=0;  //radius of gyration, ww7
  enelist[8]=0;  //distant profile from fragments
  enelist[17]=0; //CA bond-length
  enelist[18]=0;
  enelist[19]=0;
    
  int numShort=0,numMed=0,numLong=0;
  for(i=0;i<numseq;i++)
  {
    ind1=decstr[i].iaa; // identity of 20 amino acids [1-20]
    pin[0]=bf.setv(decstr[i].x,decstr[i].y,decstr[i].z);             //CA
    pin[1]=bf.setv(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z); //N
    pin[2]=bf.setv(decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z); //C
    pin[3]=bf.setv(decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z); //O
    pin[4]=bf.setv(decstr[i].ptb.x,decstr[i].ptb.y,decstr[i].ptb.z); //CB
    pin[5]=bf.setv(decstr[i].pth.x,decstr[i].pth.y,decstr[i].pth.z); //H
    pin[6]=bf.setv(decstr[i].ptg.x,decstr[i].ptg.y,decstr[i].ptg.z); //residue center 
                                                                     //for solvation
    pin[7]=bf.setv(decstr[i].ptsg.x,decstr[i].ptsg.y,decstr[i].ptsg.z); //SG
    ind5[0]=Tab5[ind1][atomid[0]];
    ind5[1]=Tab5[ind1][atomid[1]];
    ind5[2]=Tab5[ind1][atomid[2]];
    ind5[3]=Tab5[ind1][atomid[3]];
    ind5[4]=Tab5[ind1][atomid[4]];
    tcen.x+=pin[0].x;tcen.y+=pin[0].y;tcen.z+=pin[0].z;
    for(j=i+1;j<numseq;j++)
    {
      ind2=decstr[j].iaa;
      pjn[0]=bf.setv(decstr[j].x,decstr[j].y,decstr[j].z); //CA
      pjn[1]=bf.setv(decstr[j].ptn.x,decstr[j].ptn.y,decstr[j].ptn.z); //N
      pjn[2]=bf.setv(decstr[j].ptc.x,decstr[j].ptc.y,decstr[j].ptc.z); //C
      pjn[3]=bf.setv(decstr[j].pto.x,decstr[j].pto.y,decstr[j].pto.z); //O
      pjn[4]=bf.setv(decstr[j].ptb.x,decstr[j].ptb.y,decstr[j].ptb.z); //CB
      pjn[5]=bf.setv(decstr[j].pth.x,decstr[j].pth.y,decstr[j].pth.z); //H
      pjn[6]=bf.setv(decstr[j].ptg.x,decstr[j].ptg.y,decstr[j].ptg.z);
      pjn[7]=bf.setv(decstr[j].ptsg.x,decstr[j].ptsg.y,decstr[j].ptsg.z); //SG
      ind6[0]=Tab5[ind2][atomid[0]];
      ind6[1]=Tab5[ind2][atomid[1]];
      ind6[2]=Tab5[ind2][atomid[2]];
      ind6[3]=Tab5[ind2][atomid[3]];
      ind6[4]=Tab5[ind2][atomid[4]];
      for(ii=0;ii<NUM_ATOMS;ii++)
      {
        if(ii==H_INDEX  || ii==G_INDEX) continue;
        for(jj=0;jj<NUM_ATOMS;jj++)
        {
          if(jj==H_INDEX  || jj==G_INDEX) continue;
          tp=bf.minu(pin[ii],pjn[jj]);
          tdist=tp.x*tp.x+tp.y*tp.y+tp.z*tp.z; //r*r
          tdist2=sqrt(tdist); //r
          int seqsep=fabs(j-i);
          if(ii==CB_INDEX && jj==CB_INDEX)
          {
            if(tdist2<=8.0)
            {
              if(seqsep>=6 && seqsep<12)
              {
                numShort++;
              }
              else if(seqsep<24)
              {
                numMed++;
              }
              else if(seqsep>=24)
              {
                numLong++;
              }
            }
          }

          //-----implement contacts between any two atom types---->     
          if(i<=j)
          {
            index=i*numseq-(i-1)*i/2+j-i;
          }
          else
          {
            index=j*numseq-(j-1)*j/2+i-j;
          }
          if(contactConstr[ii][jj][index]>0.0) //Contacts are specified
          {
            r1=tdist2; //distance between i and j
            weight=contactConstr[ii][jj][index];
            
            if(r1<=d8)
            {
              enelist[19]+=0.0;
            }
            else
            {
              enelist[19]+=weight*(r1-d8)*(r1-d8);
            }

/*
            if(r1<=d8)
            {
              enelist[4]+=-weight;
            }
            else if(r1<d10)
            {
              enelist[4]+=-weight*(1-sin((r1-da)/db*PI))/2;
            }
            else if(r1<80)
            {
              enelist[4]+=weight*(1+sin((r1-dc)/dd*PI))/2;
            }
            else
            {
              enelist[4]+=weight;
            }
*/
          }

          if(distanceConstr[ii][jj][index]>0.0  &&
             distanceWeight[ii][jj][index]>0.0
          )
          {
            distance=distanceConstr[ii][jj][index];
            weightD=distanceWeight[ii][jj][index];
            restrType=distRestrType[ii][jj][index]; 

            if(restrType==0)
	    {
              if(fabs(i-j)>=24)
              {
                enelist[5]+=0.1*weightD*((tdist2-distance)*(tdist2-distance));
              }
              else if(fabs(i-j)>=12)
              {
                enelist[5]+=4*0.1*weightD*((tdist2-distance)*(tdist2-distance));
              }
              else
              {
                enelist[5]+=4*0.1*weightD*((tdist2-distance)*(tdist2-distance));
              }
            }
            else
            {
              if(fabs(i-j)>=24)
              {
                enelist[5]+=-1*weightD*(1.0/(4.0+0.5*(tdist2-distance)*(tdist2-distance)));
              }
              else if(fabs(i-j)>=12)
              {
                enelist[5]+=-1*weightD*(1.5/(2.0+1.0*(tdist2-distance)*(tdist2-distance)));
              }
              else
              {
                enelist[5]+=-1*weightD*(1.0/(1.0+1.0*(tdist2-distance)*(tdist2-distance)));
              }
	    }
          }

          // ^^^^^^^^^^^ contact/distance energy calculation completed ^^^^^^^^^^^^^^

          //---------excluded volume ----->
          tdist-=vdwds[ii][jj];
          if(tdist<0) enelist[1]-=tdist;

          if(ii<5 && jj<5 && tdist2<15.0 && !(ii==4 && ind1==5) && !(jj==4 && ind2==5))
          {
            tbin=int(tdist2/DELTA);
            enelist[0]+=kbp[ind5[ii]][ind6[jj]][tbin]; //RW potential, from 'data.dat'
          }

          //-------- Bond-length ----->
          //if(ii==CA_INDEX && jj==CA_INDEX && j==i+1 && tdist2>cutca) //Ca(i,i+1)
          //{
          //  enelist[17]+=(tdist2-cutca)*(tdist2-cutca); //cutca=4
          //}

          if(ii==CA_INDEX && jj==CA_INDEX)
          {
            indm=i*numseq+j;
            if(paa[indm].ispair)
            {
              if(abs(tdist2-paa[indm].dist)>paa[indm].dstd)
              {
                enelist[8]++;
              }
            }
          }
          if(ii<H_INDEX && jj<H_INDEX && tdist2<15.0 && !(ii==CB_INDEX && ind1==GLY_INDEX) &&
             !(jj==CB_INDEX && ind2==GLY_INDEX))
          { //backbone atoms, ii=4(CB), ind1=5(GLY)
            //------- distance-profile/contact.txt ----->
            //if(ii==CA_INDEX && jj==CA_INDEX)
            //{ //CA distance profile
            //  maxdist=sdist;
            //  if(tbin>maxdist) tbin=maxdist;//as reference state
            //  //helix sheet contact from fragment
            //  indm=i*numseq+j;
            //  if(paa[indm].ispair) enelist[8]-=paa[indm].ratio[tbin]*0.1;
            //}
          }
        } //end of jj
      } //end of ii
      //g

      //-------- Solvation energy ----------->
      tp.x=fabs(pin[6].x-pjn[6].x);
      if(tp.x<=distcut)
      {
        tp.y=fabs(pin[6].y-pjn[6].y);
        if(tp.y<=distcut)
        {
          tp.z=fabs(pin[6].z-pjn[6].z);
          if(tp.z<=distcut)
          {
            sqdist=tp.x*tp.x+tp.y*tp.y+tp.z*tp.z;
            if(sqdist<=sqdistcut)
            {
              if(sqdist<1.0) sqdist=1.0;
              tdists=scalf/sqdist;
              soldat[i]+=tdists*asarea2[ind2];
              soldat[j]+=tdists*asarea2[ind1];
            }
          }
        }
      } //if
    } //j
    soldat[i]+=scalf*asarea2[ind1]/16.00;
  } //i

  //------ re-weight some energy terms ---------->
  double wtdp=0.60;
  if(protType==2) //non-alpha proteins
  {
    wtdp=3.00;
  }
  else if(protType==1) //alpha proteins
  {
    wtdp=0.80;
  }
  //wtdp*=wdpE; //reweight the weight of distProf

  enelist[1]*=0.03;      //ev 0.03
  enelist[4]*=1.0;       //pr1.00
  //enelist[8]*=wtdp*0.01; //distProf, cont0.60  beta 3.00 frag5,>1L   alpha 0.60, frag7,>4L

  //------- radius gyration --------------->
  double trad=0.0;
  tcen.x/=double(numseq);tcen.y/=double(numseq);tcen.z/=double(numseq);
  for(i=0;i<numseq;i++)
  {
    tp.x=decstr[i].x-tcen.x;
    tp.y=decstr[i].y-tcen.y;
    tp.z=decstr[i].z-tcen.z;
    tdist=tp.x*tp.x+tp.y*tp.y+tp.z*tp.z;
    trad+=tdist;
  }
  trad/=double(numseq);
  trad=sqrt(trad);
  double minradius=exp(0.840)*pow(numseq,0.358)-0.5;
  double maxradius=minradius+8.0;
  double tmaxval=longestdist*0.5*sqrt(0.6);
  if(tmaxval>maxradius) maxradius=tmaxval;
  if(trad>maxradius)
  {
    enelist[7]=(trad-maxradius)*(trad-maxradius);
  }
  else if(trad<minradius)
  {
    enelist[7]=(minradius-trad)*(minradius-trad);
  }
  enelist[7]*=0.75;

  //------------ contact number -------------->
  if(numseq>=50){
    enelist[4]+=fabs(expectedShort-numShort);
    enelist[4]+=fabs(expectedMed-numMed);
    enelist[4]+=fabs(expectedLong-numLong);
    enelist[4]*=0.4;
  }

  //------- solvation potential -------------->
  for(i=0;i<numseq;i++)
  {
    ind1=decstr[i].iaa;
    soldat[i]*=scalaa[ind1];
    enelist[18]+=soldat[i]*solweight[i];
  }
  enelist[18]*=0.40;//4.00
  delete[]soldat;
}

double FoldDesignEnergyFunction::calcHbondBackboneEnergy(
  point3f *decstr,
  int numseq
)
{
  int i,j,k,l,m;
  int protType=inputInfo.getProteinType();
  point3d tp[20],ap[20],kp[20];
  point3d pd,pd2,pd3;
  double lamda;
  double diss[20];
  double totenergy=0;
  double totenergy2=0;

  double threshval[5]={19.113828,19.113828,20.723266,20.723266,16.118096};
  for(i=0;i<numseq;i++)
  {
    decstr[i].vpos=0;
    decstr[i].tpos=0;
    decstr[i].ssm='C';
    decstr[i].indl=-1;
    decstr[i].indr=-1;
    decstr[i].tpr=-1;
    decstr[i].tpl=-1;
  }
  //alpha
  for(i=1;i<numseq-4;i++)
  {
    j=i+3;
    if((decstr[i].ss3=='E' &&  decstr[i].stype!='H') ||
       (decstr[j].ss3=='E' &&  decstr[j].stype!='H'))
    {
      continue;
    }
    l=0;
    for(k=0;k<4;k++) if(decstr[i+k].ss3=='H') l++;
    m=0;
    for(k=0;k<4;k++) if(decstr[i+k].stype=='H') m++;
    if(l==3 && m<2 && (decstr[i+1].ss3!='H' || decstr[i+2].ss3!='H'))
    {
      continue;
    }
    else if(l==2 && m<2 &&
            ((decstr[i+1].ss3=='H' && decstr[i+2].ss3=='H') ||
            (decstr[i+1].ss3!='H' && decstr[i+2].ss3!='H')))
    {
      continue;
    }
    if(l==0 && m<3) continue;
    else if(l==1 && m<3) continue;
    else if(l==2 && m<3) lamda=0.05;
    else if(l==3) lamda=0.4;
    else if(l==4) lamda=1.5;
    else if(m==3) lamda=0.6;
    else if(m==4) lamda=1.0;
    else lamda=0.3;

    tp[1]=bf.setv(decstr[i].x,decstr[i].y,decstr[i].z);
    ap[1]=bf.setv(decstr[j].x,decstr[j].y,decstr[j].z);
    tp[7]=bf.minu(tp[1],ap[1]);
    diss[1]=bf.norm(tp[7]);
    if(diss[1]>6.6) continue;

    tp[0]=bf.setv(decstr[i-1].x,decstr[i-1].y,decstr[i-1].z);
    tp[2]=bf.setv(decstr[i+1].x,decstr[i+1].y,decstr[i+1].z);
    ap[0]=bf.setv(decstr[j-1].x,decstr[j-1].y,decstr[j-1].z);
    ap[2]=bf.setv(decstr[j+1].x,decstr[j+1].y,decstr[j+1].z);
    kp[0]=bf.minu(tp[2],ap[0]);
    diss[4]=bf.norm(kp[0]);
    if(diss[4]<3.4 || diss[4]>4.2) continue;
    diss[0]=bf.norm(bf.minu(tp[0],ap[0]));

    //nhoc
    pd.x=decstr[j].pth.x-decstr[i].pto.x;
    pd.y=decstr[j].pth.y-decstr[i].pto.y;
    pd.z=decstr[j].pth.z-decstr[i].pto.z;
    diss[5]=bf.norm(pd);
    if(diss[5]>=5.0 || diss[5]<1.6) continue;

    pd2.x=decstr[i].ptc.x-decstr[i].pto.x;
    pd2.y=decstr[i].ptc.y-decstr[i].pto.y;
    pd2.z=decstr[i].ptc.z-decstr[i].pto.z;
    diss[6]=bf.angv(pd,pd2)*degrad;
    if(diss[6]<70.0 || diss[6]>140.0) continue;

    pd3.x=decstr[j].pth.x-decstr[j].ptn.x;
    pd3.y=decstr[j].pth.y-decstr[j].ptn.y;
    pd3.z=decstr[j].pth.z-decstr[j].ptn.z;
    diss[7]=bf.angv(pd,pd3)*degrad;

    diss[8]=bf.phi(decstr[j].ptn.x,decstr[j].ptn.y,decstr[j].ptn.z,
                   decstr[j].pth.x,decstr[j].pth.y,decstr[j].pth.z,
                   decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,
                   decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z);
    if(diss[8]<120.0 || diss[8]>280.0) continue;

    //second
    j=i+4;
    pd.x=decstr[j].pth.x-decstr[i].pto.x;
    pd.y=decstr[j].pth.y-decstr[i].pto.y;
    pd.z=decstr[j].pth.z-decstr[i].pto.z;
    diss[9]=bf.norm(pd);
    if(diss[9]>=5.0 || diss[9]<1.6) continue;

    pd2.x=decstr[i].ptc.x-decstr[i].pto.x;
    pd2.y=decstr[i].ptc.y-decstr[i].pto.y;
    pd2.z=decstr[i].ptc.z-decstr[i].pto.z;
    diss[10]=bf.angv(pd,pd2)*degrad;
    if(diss[10]<100.0 || diss[10]>170.0) continue;

    pd3.x=decstr[j].pth.x-decstr[j].ptn.x;
    pd3.y=decstr[j].pth.y-decstr[j].ptn.y;
    pd3.z=decstr[j].pth.z-decstr[j].ptn.z;
    diss[11]=bf.angv(pd,pd3)*degrad;
    if(diss[11]<100.0) continue;

    diss[12]=bf.phi(decstr[j].ptn.x,decstr[j].ptn.y,decstr[j].ptn.z,
                    decstr[j].pth.x,decstr[j].pth.y,decstr[j].pth.z,
                    decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,
                    decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z);
    if(diss[12]<100.0 || diss[12]>230.0) continue;

    double tval=0;
    for(k=0;k<1;k++)
    {
      tval+=bf.squgaussian(diss[k],nhochhsigma[k],nhochhmean[k]);
    }
    for(k=0;k<4;k++)
    {
      tval+=bf.squgaussian(diss[5+k],nhocstdval[0][k],nhocmeanval[0][k]);
    }
    for(k=0;k<4;k++)
    {
      tval+=bf.squgaussian(diss[9+k],nhocstdval[1][k],nhocmeanval[1][k]);
    }
    tval=exp(tval);
    if(tval>1e-7)
    {
      if(tval>1e-3) tval=1e-3;
      tval=lamda*(threshval[4]+log(tval));
      totenergy+=tval;
      for(k=0;k<=3;k++) decstr[i+k].ssm='H';
    }
  }

  ////beta
  int inda,indb,indc,indd,inde,indf;
  int delseg;

  for(j=1;j<numseq-4;j++)
  {
    for(k=j+3;k<numseq-1;k++)
    {
      if((decstr[j].ss3=='H' && decstr[j].stype!='E') ||
         (decstr[k].ss3=='H' && decstr[k].stype!='E'))//diff
      {
        continue;
      }
      lamda=0.2;
      if((decstr[j].ss3=='E' && decstr[j].stype=='H') &&
         (decstr[k].ss3=='E' && decstr[k].stype=='H'))
      {
        lamda+=0.1;
      }
      else if((decstr[j].ss3=='E' && decstr[j].stype=='H') ||
              (decstr[k].ss3=='E' && decstr[k].stype=='H'))
      {
        lamda+=0.2;
      }
      else if(decstr[j].ss3=='E' && decstr[k].ss3=='E')
      {
        lamda+=1.6;
      }
      else if(decstr[k].ss3=='E' || decstr[j].ss3=='E')
      {
        lamda+=0.6;
      }
      tp[1]=bf.setv(decstr[j].x,decstr[j].y,decstr[j].z);
      ap[1]=bf.setv(decstr[k].x,decstr[k].y,decstr[k].z);
      tp[7]=bf.minu(tp[1],ap[1]);
      if(bf.norm(tp[7])>8.0) continue;

      inda=j;indb=k;
      tp[4]=bf.setv(decstr[inda].ptn.x,decstr[inda].ptn.y,decstr[inda].ptn.z);
      tp[5]=bf.setv(decstr[inda+1].ptn.x,decstr[inda+1].ptn.y,decstr[inda+1].ptn.z);
      tp[9]=bf.setv(decstr[inda-1].pto.x,decstr[inda-1].pto.y,decstr[inda-1].pto.z);
      tp[10]=bf.setv(decstr[inda].pto.x,decstr[inda].pto.y,decstr[inda].pto.z);
      ap[4]=bf.setv(decstr[indb].ptn.x,decstr[indb].ptn.y,decstr[indb].ptn.z);
      ap[5]=bf.setv(decstr[indb+1].ptn.x,decstr[indb+1].ptn.y,decstr[indb+1].ptn.z);
      ap[9]=bf.setv(decstr[indb-1].pto.x,decstr[indb-1].pto.y,decstr[indb-1].pto.z);
      ap[10]=bf.setv(decstr[indb].pto.x,decstr[indb].pto.y,decstr[indb].pto.z);
      kp[0]=bf.minu(tp[9],ap[4]);  //oi-1 nj
      kp[1]=bf.minu(tp[5],ap[10]); //ni+1 oj
      kp[2]=bf.minu(tp[4],ap[9]);  //ni oj-1
      kp[3]=bf.minu(tp[10],ap[5]); //oi nj+1
      kp[4]=bf.minu(tp[9],ap[5]);  //oi-1 nj+1
      kp[5]=bf.minu(tp[5],ap[9]);  //ni+1 oj-1
      kp[6]=bf.minu(tp[4],ap[10]); //ni oj
      kp[7]=bf.minu(tp[10],ap[4]); //oi nj
      for(m=0;m<8;m++) diss[m]=bf.norm(kp[m]);
      delseg=0;
      diss[9]=diss[0]+diss[1];
      diss[8]=diss[2]+diss[3];
      if(diss[8]<diss[9])
      {
        delseg=1;
        diss[9]=diss[8];
      }
      diss[8]=diss[4]+diss[5];
      if(diss[8]<diss[9])
      {
        delseg=2;
        diss[9]=diss[8];
      }
      diss[8]=diss[6]+diss[7];
      if(diss[8]<diss[9])
      {
        delseg=3;
        diss[9]=diss[8];
      }
      if(diss[2*delseg]>6.8) continue;
      if(diss[2*delseg+1]>6.8) continue;
      if(delseg==0)
      {
        indc=k;indd=j-1;inde=j+1;indf=k;
      }
      else if(delseg==1)
      {
        indc=j;indd=k-1;inde=k+1;indf=j;
      }
      else if(delseg==3)
      {
        indc=k;indd=j;inde=j;indf=k;
      }
      else if(delseg==2)
      {
        indc=k+1;indd=j-1;inde=j+1;indf=k-1;
      }
      pd.x=decstr[indc].pth.x-decstr[indd].pto.x;
      pd.y=decstr[indc].pth.y-decstr[indd].pto.y;
      pd.z=decstr[indc].pth.z-decstr[indd].pto.z;
      diss[0]=bf.norm(pd);
      pd2.x=decstr[indd].ptc.x-decstr[indd].pto.x;
      pd2.y=decstr[indd].ptc.y-decstr[indd].pto.y;
      pd2.z=decstr[indd].ptc.z-decstr[indd].pto.z;
      diss[1]=bf.angv(pd,pd2)*degrad;
      pd3.x=decstr[indc].pth.x-decstr[indc].ptn.x;
      pd3.y=decstr[indc].pth.y-decstr[indc].ptn.y;
      pd3.z=decstr[indc].pth.z-decstr[indc].ptn.z;
      diss[2]=bf.angv(pd,pd3)*degrad;
      ////////
      pd.x=decstr[inde].pth.x-decstr[indf].pto.x;
      pd.y=decstr[inde].pth.y-decstr[indf].pto.y;
      pd.z=decstr[inde].pth.z-decstr[indf].pto.z;
      diss[4]=bf.norm(pd);
      pd2.x=decstr[indf].ptc.x-decstr[indf].pto.x;
      pd2.y=decstr[indf].ptc.y-decstr[indf].pto.y;
      pd2.z=decstr[indf].ptc.z-decstr[indf].pto.z;
      diss[5]=bf.angv(pd,pd2)*degrad;
      pd3.x=decstr[inde].pth.x-decstr[inde].ptn.x;
      pd3.y=decstr[inde].pth.y-decstr[inde].ptn.y;
      pd3.z=decstr[inde].pth.z-decstr[inde].ptn.z;
      diss[6]=bf.angv(pd,pd3)*degrad;

      double tval=0;
      for(m=0;m<3;m++)
      {
        tval+=bf.squgaussian(diss[m],nhocstdval[2+delseg][m],nhocmeanval[2+delseg][m]);
      }
      for(m=0;m<3;m++)
      {
        tval+=bf.squgaussian(diss[4+m],nhocstdval[2+delseg][m],nhocmeanval[2+delseg][m]);
      }
      tval=exp(tval);

      if((delseg==0 && tval>5e-9) || (delseg==1 && tval>5e-9) || 
         (delseg==2 && tval>1e-9) || (delseg==3 && tval>1e-9))
      {
        if(tval>1e-3) tval=1e-3;
        tval=lamda*(threshval[delseg]+log(tval))/threshval[delseg]*threshval[4];
        if(delseg==0 || delseg==2) //Left of j
        {
          if(tval>decstr[j].vpos)
          {
            decstr[j].vpos=tval;
            if(decstr[j].indl>0)
            {
              if(decstr[decstr[j].indl].indr==j) //Set to zero
              {
                decstr[decstr[j].indl].tpos=0;
                decstr[decstr[j].indl].indr=-1;
                if(decstr[decstr[j].indl].indl==-1) decstr[decstr[j].indl].ssm='C';
              }
            }
            decstr[j].indl=k;
            decstr[j].tpl=delseg;
            decstr[k].ssm='E';
            decstr[j].ssm='E';
          }
        }
        else if(delseg==1 || delseg==3) //Right of j
        {
          if(tval>decstr[j].tpos)
          {
            decstr[j].tpos=tval;
            if(decstr[j].indr>0)
            {
              if(decstr[decstr[j].indr].indl==j) //Set to zero
              {
                decstr[decstr[j].indr].vpos=0;
                decstr[decstr[j].indr].indl=-1;
                if(decstr[decstr[j].indr].indr==-1) decstr[decstr[j].indr].ssm='C';
              }
            }
            decstr[j].indr=k;
            decstr[j].tpr=delseg;
            decstr[k].ssm='E';
            decstr[j].ssm='E';
          }
        }
        if(delseg==1 || delseg==2) //Left of k
        {
          if(tval>decstr[k].vpos)
          {
            decstr[k].vpos=tval;
            if(decstr[k].indl>0)
            {
              if(decstr[decstr[k].indl].indr==k) //Set to zero
              {
                decstr[decstr[k].indl].tpos=0;
                decstr[decstr[k].indl].indr=-1;
                if(decstr[decstr[k].indl].indl==-1) decstr[decstr[k].indl].ssm='C';
              }
            }
            decstr[k].indl=j;
            decstr[k].tpl=delseg;
            decstr[j].ssm='E';
            decstr[k].ssm='E';
          }
        }
        else if(delseg==0 || delseg==3) //Right of k
        {
          if(tval>decstr[k].tpos)
          {
            decstr[k].tpos=tval;
            if(decstr[k].indr>0)
            {
              if(decstr[decstr[k].indr].indl==k) //Set to zero
              {
                decstr[decstr[k].indr].vpos=0;
                decstr[decstr[k].indr].indl=-1;
                if(decstr[decstr[k].indr].indr==-1) decstr[decstr[k].indr].ssm='C';
              }
            }
            decstr[k].indr=j;
            decstr[k].tpr=delseg;
            decstr[j].ssm='E';
            decstr[k].ssm='E';
          }
        }
      }
    }
  }

  double wthb=2.0;
  double wttot=3.0;
  if(protType==2)
  {
    wthb=3.0;
    wttot=2.0;
  }
  for(i=0;i<numseq;i++)
  {
    totenergy2+=decstr[i].vpos;
    totenergy2+=decstr[i].tpos;
    if(decstr[i].indl!=-1 && decstr[i].indr!=-1)
    {
      totenergy2+=wthb; //beta 3.0 alpha 2.0
      if((decstr[i].tpl==0 || decstr[i].tpl==1) && (decstr[i].tpr==0 || decstr[i].tpr==1))
        totenergy2+=1.0; //+decstr[i].vpos+decstr[i].tpos;
      else if((decstr[i].tpl==2 || decstr[i].tpl==3) && (decstr[i].tpr==2 || decstr[i].tpr==3))
        totenergy2+=1.5; //+decstr[i].vpos+decstr[i].tpos;
      else if((decstr[i].tpl==0 || decstr[i].tpl==1) && (decstr[i].tpr==2 || decstr[i].tpr==3))
        totenergy2-=0.0;
      else if((decstr[i].tpl==2 || decstr[i].tpl==3) && (decstr[i].tpr==0 || decstr[i].tpr==1))
        totenergy2-=0.0;
    }
    else if(decstr[i].indl!=-1 || decstr[i].indr!=-1)
      totenergy2+=1.5;
  }

  return -(1.0*totenergy+wttot*totenergy2); //beta 1,2  alpha 1,3
}

double FoldDesignEnergyFunction::calcSSConstrEnergy(
  point3f *decstr,
  int numseq
)
{
  const int MATCH_SCORE_COIL=1;
  const int MATCH_SCORE=2;
  const int MISMATCH_SCORE=2;
  const int MISMATCH_SCORE_ONE_COIL=1;
  int i;
  double ssMatch=0.0;

  for(i=0;i<numseq;i++)
  {
    if(decstr[i].ss3=='X')
    {
      continue;
    }
    else if(decstr[i].ssm==decstr[i].ss3 &&
            (decstr[i].ss3=='H' || decstr[i].ss3=='E'))
    {
      ssMatch+=MATCH_SCORE;
    }
    else if((decstr[i].ssm=='H' && decstr[i].ss3=='E') ||
            (decstr[i].ssm=='E' && decstr[i].ss3=='H'))
    {
      ssMatch-=MISMATCH_SCORE;
    }
    else if(decstr[i].ssm==decstr[i].ss3)
    {
      ssMatch+=MATCH_SCORE_COIL;
    }
    else
    {
      ssMatch-=MISMATCH_SCORE_ONE_COIL;
    }
  }

  return -1.0*double(ssMatch/numseq);
}

double FoldDesignEnergyFunction::calcRamaEnergy(
  point3f *decstr,
  int numseq
)
{
  int i;
  double phi,psi;
  double totene=0.0;

  for(i=1;i<numseq-1;i++)
  {
    phi=decstr[i].tor[2];
    psi=decstr[i].tor[0];
    if(phi==-360 || psi==-360) continue;
    if(phi>=360.0)
    {
      phi-=360.0;
    }
    if(psi>=360)
    {
      psi-=360.0;
    }
    //printf("i %d phi %lf psi %lf\n",i,phi,psi);

    if(decstr[i].ss8=='H')
    {
      totene+=H_general_rama[int(phi)][int(psi)];
    }
    else if(decstr[i].ss8=='E')
    {
      totene+=E_general_rama[int(phi)][int(psi)];
    }
    else if(decstr[i].ss8=='C')
    {
      totene+=C_general_rama[int(phi)][int(psi)];
    }
    else if(decstr[i].ss8=='h')
    {
      totene+=H_specific_rama[int(phi)][int(psi)];
    }
    else if (decstr[i].ss8=='e')
    {
      totene+=E_specific_rama[int(phi)][int(psi)];
    }
    else if (decstr[i].ss8=='c')
    {
      totene+=C_specific_rama[int(phi)][int(psi)];
    }
    else if(decstr[i].ss8=='G')
    {
      totene+=G_specific_rama[int(phi)][int(psi)];
    }
    else if (decstr[i].ss8=='I')
    {
      totene+=I_specific_rama[int(phi)][int(psi)];
    }
    else if (decstr[i].ss8=='B')
    {
      totene+=B_specific_rama[int(phi)][int(psi)];
    }
    else if(decstr[i].ss8=='S')
    {
      totene+=S_specific_rama[int(phi)][int(psi)];
    }
    else if (decstr[i].ss8=='T')
    {
      totene+=T_specific_rama[int(phi)][int(psi)];
    }
  }

  return totene;
}

double FoldDesignEnergyFunction::calcBABMotifPenalty(
  point3f *decstr,
  int numseq
)
{
  double totene=0;
  double wt=0;
  int i,j;
  point3d tp[30];
  double tdist[5];
  
  int numBAB=inputInfo.getNumBAB();
  triplebab *babInfo=inputInfo.getBABInfo();
  abssinfo *abss=inputInfo.getABSS();

  for(i=0;i<numBAB;i++)
  {
    // a, b, c means first strand, helix, and second strand, respecitvely
    // initial CA atom of first strand
    tp[0]=bf.setv(decstr[abss[babInfo[i].a].seg.init].x,
                  decstr[abss[babInfo[i].a].seg.init].y,
                  decstr[abss[babInfo[i].a].seg.init].z);
    // Terminal CA atom of first strand
    tp[1]=bf.setv(decstr[abss[babInfo[i].a].seg.term].x,
                  decstr[abss[babInfo[i].a].seg.term].y,
                  decstr[abss[babInfo[i].a].seg.term].z);
    // Initial CA atom of helix
    tp[2]=bf.setv(decstr[abss[babInfo[i].b].seg.init].x,
                  decstr[abss[babInfo[i].b].seg.init].y,
                  decstr[abss[babInfo[i].b].seg.init].z);
    // Terminal CA atom of helix
    tp[3]=bf.setv(decstr[abss[babInfo[i].b].seg.term].x,
                  decstr[abss[babInfo[i].b].seg.term].y,
                  decstr[abss[babInfo[i].b].seg.term].z);
    // Initial CA atom of second strand
    tp[4]=bf.setv(decstr[abss[babInfo[i].c].seg.init].x,
                  decstr[abss[babInfo[i].c].seg.init].y,
                  decstr[abss[babInfo[i].c].seg.init].z);
    // Terminal CA atom of second strand
    tp[5]=bf.setv(decstr[abss[babInfo[i].c].seg.term].x,
                  decstr[abss[babInfo[i].c].seg.term].y,
                  decstr[abss[babInfo[i].c].seg.term].z);
    tp[6]=bf.minu(tp[1],tp[0]);    // vector of first strand
    tp[8]=bf.minu(tp[5],tp[4]);    // vector of second strand
    tdist[0]=bf.angv(tp[6],tp[8]); // angle between two strands
    if(tdist[0]>PI/3.0) continue;
    tdist[4]=PI/3.0-tdist[0];
    tp[9]=bf.addv(tp[1],tp[0]);
    tp[11]=bf.scal(tp[9],0.5);     // midpoint of first strand
    tp[7]=bf.addv(tp[3],tp[2]);
    tp[13]=bf.scal(tp[7],0.5);     // midpoint of helix
    tp[10]=bf.addv(tp[5],tp[4]);
    tp[12]=bf.scal(tp[10],0.5);    // midpoint of second strand

    tp[14]=bf.addv(tp[11],tp[12]);
    tp[15]=bf.scal(tp[14],0.5);    // midpoint between two strands
    tp[16]=bf.minu(tp[13],tp[15]); // from midpoint of strands to helix
    tp[17]=bf.minu(tp[12],tp[11]); // from first strand midpoint to
                                   // second strand midpoint
    tp[18]=bf.prod(tp[17],tp[6]);  // perpendicular to both first strand
                                   // and line connecting strand midpoints
    tp[19]=bf.prod(tp[17],tp[8]);  // perpendicular to both second strand
                                   // and line connecting strand midpoints
    tdist[1]=bf.angv(tp[18],tp[16]);
    tdist[2]=bf.angv(tp[19],tp[16]);
    tdist[3]=(tdist[1]>PI/2.0)+(tdist[2]>PI/2.0);
    if(tdist[3]<2) continue;

    for(j=abss[babInfo[i].a].seg.init;j<=abss[babInfo[i].a].seg.term;j++)
    {
      if(decstr[j].indl>=abss[babInfo[i].c].seg.init &&
        decstr[j].indl<=abss[babInfo[i].c].seg.term)
      {
        totene+=1.0;
      }
      else if(decstr[j].indr>=abss[babInfo[i].c].seg.init &&
              decstr[j].indr<=abss[babInfo[i].c].seg.term)
      {
        totene+=1.0;
      }
    }

    for(j=abss[babInfo[i].c].seg.init;j<=abss[babInfo[i].c].seg.term;j++)
    {
      if (decstr[j].indl>=abss[babInfo[i].a].seg.init &&
          decstr[j].indl<=abss[babInfo[i].a].seg.term)
      {
        totene+=1.0;
      }
      else if(decstr[j].indr>=abss[babInfo[i].a].seg.init &&
              decstr[j].indr<=abss[babInfo[i].a].seg.term)
      {
        totene+=1.0;
      }
    }
  }

  return 3.0*totene;
}


vector<double> FoldDesignEnergyFunction::calcStrandStrandEnergy(
  point3f *decstr,
  int numseq
)
{
  const int DIST_SCALE=10;
  const double OC_SCALE=1.0/1.231015;
  int j,k;
  int seqSep;
  double hwt=2.879385;
  double phi,psi,theta,distanceSSE;
  double lenStrandX,lenStrandY,lenStrandZ,lenStrand,normLenStrandZ,cte_dot;
  point3d cenStrand1,cenStrand2,cenStrand1Temp,cenStrand2Temp,deltaSS;
  point3d sse1Start,sse1End,sse2Start,sse2End,ap[10];
  point3d axisX,axisY,axisZ,vecStrand2,unitSS;
  point3d cte;
  point3d cv[4],ov[4],dissOC[8];
  double sdotOC1,sdotOC2,tempdot[4];
  int sign1ss,sign2ss;
  double energyPhiPsiTheta=0;
  double energyDistTheta=0;
  vector<double> totene(2,0);
  sssegment *sse=inputInfo.getSSE();
  int numsse=inputInfo.getNumSSE();

  for(j=0;j<numsse-1;j++)
  {
    if(sse[j].ss=='E')
    {
      if(sse[j].term-sse[j].init<1) continue;
      ap[0]=bf.setv(decstr[sse[j].init].ptn.x,
                    decstr[sse[j].init].ptn.y,
                    decstr[sse[j].init].ptn.z);
      ap[1]=bf.setv(decstr[sse[j].term].ptc.x,
                    decstr[sse[j].term].ptc.y,
                    decstr[sse[j].term].ptc.z);
      sse1Start.x=ap[0].x;
      sse1Start.y=ap[0].y;
      sse1Start.z=ap[0].z;
      sse1End.x=ap[1].x;
      sse1End.y=ap[1].y;
      sse1End.z=ap[1].z;

      for(k=j+1;k<numsse;k++) if(sse[k].ss=='E')
      {
        if(sse[k].term-sse[k].init<1) continue;
        seqSep=sse[k].init-sse[j].term;
        if(seqSep<2) continue;

        ap[2]=bf.setv(decstr[sse[k].init].ptn.x,
                      decstr[sse[k].init].ptn.y,
                      decstr[sse[k].init].ptn.z);
        ap[3]=bf.setv(decstr[sse[k].term].ptc.x,
                      decstr[sse[k].term].ptc.y,
                      decstr[sse[k].term].ptc.z);
        sse2Start.x=ap[2].x;
        sse2Start.y=ap[2].y;
        sse2Start.z=ap[2].z;
        sse2End.x=ap[3].x;
        sse2End.y=ap[3].y;
        sse2End.z=ap[3].z;

        cenStrand1Temp=bf.addv(sse1Start,sse1End);
        cenStrand1=bf.scal(cenStrand1Temp,0.5);
        cenStrand2Temp=bf.addv(sse2Start,sse2End);
        cenStrand2=bf.scal(cenStrand2Temp,0.5);
        deltaSS=bf.minu(cenStrand2,cenStrand1);
        distanceSSE=bf.norm(deltaSS);
        if(distanceSSE>20) continue;

        axisZ=bf.unit(bf.minu(sse1End,sse1Start));
        axisY=bf.unit(bf.prod(axisZ,deltaSS));
        axisX=bf.unit(bf.prod(axisY,axisZ));
        vecStrand2=bf.minu(sse2End,sse2Start);
        lenStrandX=bf.dotv(vecStrand2,axisX);
        lenStrandY=bf.dotv(vecStrand2,axisY);
        lenStrandZ=bf.dotv(vecStrand2,axisZ);
        unitSS=bf.unit(deltaSS);
        lenStrand=bf.norm(vecStrand2);

        if(lenStrandX != 0 && lenStrandY != 0) 
        {
          phi=atan2(lenStrandY,lenStrandX)*(double)180.0/PI;
          if(phi<0) phi+=360;
        }
        else
        {
          continue;
        }

        normLenStrandZ=lenStrandZ/lenStrand;
        if(lenStrand!=0 && abs(normLenStrandZ)<=1)
        {
          psi=acos(normLenStrandZ)*(double)180.0/PI;
          if(psi<0) psi+=360;
        }
        else
        {
          continue;
        }

        cte=bf.unit(bf.minu(sse1End,cenStrand1));
        cte_dot=bf.dotv(cte,unitSS);
        theta=acos(cte_dot)*(double)180.0/PI;
        if(theta<0) theta+=180;

        if(seqSep<=11 && seqSep>1)
        {
          energyPhiPsiTheta+=energySSPackAngle2[int(phi)][int(psi)][int(theta)];
        }
        else if(seqSep>=12)
        {
          energyPhiPsiTheta+=energySSPackAngle3[int(phi)][int(psi)][int(theta)];
        }

        cv[0]=bf.setv(decstr[sse[j].init].ptc.x,
                      decstr[sse[j].init].ptc.y,
                      decstr[sse[j].init].ptc.z);
        ov[0]=bf.setv(decstr[sse[j].init].pto.x,
                      decstr[sse[j].init].pto.y,
                      decstr[sse[j].init].pto.z);
        cv[1]=bf.setv(decstr[sse[j].init+1].ptc.x,
                      decstr[sse[j].init+1].ptc.y,
                      decstr[sse[j].init+1].ptc.z);
        ov[1]=bf.setv(decstr[sse[j].init+1].pto.x,
                      decstr[sse[j].init+1].pto.y,
                      decstr[sse[j].init+1].pto.z);
        dissOC[0]=bf.minu(cv[0],ov[0]);
        dissOC[1]=bf.minu(ov[1],cv[1]);
        dissOC[2]=bf.scal(dissOC[0],OC_SCALE);
        dissOC[3]=bf.scal(dissOC[1],OC_SCALE);
        tempdot[0]=bf.dotv(dissOC[2],unitSS);
        tempdot[1]=bf.dotv(dissOC[3],unitSS);
        sdotOC1=tempdot[0]+tempdot[1];

        cv[2]=bf.setv(decstr[sse[k].init].ptc.x,
                      decstr[sse[k].init].ptc.y,
                      decstr[sse[k].init].ptc.z);
        ov[2]=bf.setv(decstr[sse[k].init].pto.x,
                      decstr[sse[k].init].pto.y,
                      decstr[sse[k].init].pto.z);
        cv[3]=bf.setv(decstr[sse[k].init+1].ptc.x,
                      decstr[sse[k].init+1].ptc.y,
                      decstr[sse[k].init+1].ptc.z);
        ov[3]=bf.setv(decstr[sse[k].init+1].pto.x,
                      decstr[sse[k].init+1].pto.y,
                      decstr[sse[k].init+1].pto.z);
        dissOC[4]=bf.minu(cv[2],ov[2]);
        dissOC[5]=bf.minu(ov[3],cv[3]);
        dissOC[6]=bf.scal(dissOC[4],OC_SCALE);
        dissOC[7]=bf.scal(dissOC[5],OC_SCALE);
        tempdot[2]=bf.dotv(dissOC[6],unitSS);
        tempdot[3]=bf.dotv(dissOC[7],unitSS);
        sdotOC2=tempdot[2]+tempdot[3];
        sign1ss=(sdotOC1>0.0 ? 2:1);
        sign2ss=(sdotOC2>0.0 ? 2:1);

        if(sign1ss==1 && sign2ss==1)
        {
          energyDistTheta+=energySSPackDist1[int(distanceSSE*DIST_SCALE)][int(theta)];
        }
        else if(sign1ss==1 && sign2ss==2)
        {
          energyDistTheta+=energySSPackDist2[int(distanceSSE*DIST_SCALE)][int(theta)];
        }
        else if(sign1ss==2 && sign2ss==1)
        {
          energyDistTheta+=energySSPackDist3[int(distanceSSE*DIST_SCALE)][int(theta)];
        }
        else
        {
          energyDistTheta+=energySSPackDist4[int(distanceSSE*DIST_SCALE)][int(theta)];
        }
      }
    }
  }
  totene[0]=energyPhiPsiTheta;
  totene[1]=energyDistTheta;

  return totene;
}

vector<double> FoldDesignEnergyFunction::calcHelixHelixEnergy(
  point3f *decstr,
  int numseq
)
{
  const int DIST_SCALE=10;
  int j,k;
  int seqSep;
  double hwt=2.879385;
  double phi,psi,theta,distanceSSE;
  double lenHelixX,lenHelixY,lenHelixZ,lenHelix,normLenHelixZ,cte_dot;
  point3d cenHelix1,cenHelix2,cenHelix1Temp,cenHelix2Temp,deltaHH;
  point3d sse1Start,sse1End,sse2Start,sse2End,ap[10];
  point3d axisX,axisY,axisZ,vecHelix2,unitHH;
  point3d cte;
  double energyPhiPsiTheta=0;
  double energyDistTheta=0;
  vector<double> totene(2,0);
  sssegment *sse=inputInfo.getSSE();
  int numsse=inputInfo.getNumSSE();

  for(j=0;j<numsse-1;j++)
  {
    if(sse[j].ss=='H')
    {
      if(sse[j].term-sse[j].init<2) continue;
      ap[0]=bf.setv(decstr[sse[j].init].x,decstr[sse[j].init].y,decstr[sse[j].init].z);
      ap[1]=bf.setv(decstr[sse[j].init+1].x,decstr[sse[j].init+1].y,decstr[sse[j].init+1].z);
      ap[2]=bf.setv(decstr[sse[j].init+2].x,decstr[sse[j].init+2].y,decstr[sse[j].init+2].z);
      ap[3]=bf.addv(ap[0],ap[2]);
      ap[4]=bf.scal(ap[3],hwt);
      ap[0]=bf.addv(ap[4],ap[1]);
      ap[2]=bf.scal(ap[0],1.0/(1.0+2.0*hwt));
      ap[5]=bf.setv(decstr[sse[j].term].x,decstr[sse[j].term].y,decstr[sse[j].term].z);
      ap[6]=bf.setv(decstr[sse[j].term-1].x,decstr[sse[j].term-1].y,decstr[sse[j].term-1].z);
      ap[7]=bf.setv(decstr[sse[j].term-2].x,decstr[sse[j].term-2].y,decstr[sse[j].term-2].z);
      ap[8]=bf.addv(ap[5],ap[7]);
      ap[9]=bf.scal(ap[8],hwt);
      ap[5]=bf.addv(ap[9],ap[6]);
      ap[7]=bf.scal(ap[5],1.0/(1.0+2.0*hwt));
      sse1Start.x=ap[2].x;
      sse1Start.y=ap[2].y;
      sse1Start.z=ap[2].z;
      sse1End.x=ap[7].x;
      sse1End.y=ap[7].y;
      sse1End.z=ap[7].z;

      for(k=j+1;k<numsse;k++) if(sse[k].ss=='H')
      {
        if(sse[k].term-sse[k].init<2) continue;
        seqSep=sse[k].init-sse[j].term;
        if(seqSep<2) continue;

        ap[0]=bf.setv(decstr[sse[k].init].x,decstr[sse[k].init].y,decstr[sse[k].init].z);
        ap[1]=bf.setv(decstr[sse[k].init+1].x,decstr[sse[k].init+1].y,decstr[sse[k].init+1].z);
        ap[2]=bf.setv(decstr[sse[k].init+2].x,decstr[sse[k].init+2].y,decstr[sse[k].init+2].z);
        ap[3]=bf.addv(ap[0],ap[2]);
        ap[4]=bf.scal(ap[3],hwt);
        ap[0]=bf.addv(ap[4],ap[1]);
        ap[2]=bf.scal(ap[0],1.0/(1.0+2.0*hwt));
        ap[5]=bf.setv(decstr[sse[k].term].x,decstr[sse[k].term].y,decstr[sse[k].term].z);
        ap[6]=bf.setv(decstr[sse[k].term-1].x,decstr[sse[k].term-1].y,decstr[sse[k].term-1].z);
        ap[7]=bf.setv(decstr[sse[k].term-2].x,decstr[sse[k].term-2].y,decstr[sse[k].term-2].z);
        ap[8]=bf.addv(ap[5],ap[7]);
        ap[9]=bf.scal(ap[8],hwt);
        ap[5]=bf.addv(ap[9],ap[6]);
        ap[7]=bf.scal(ap[5],1.0/(1.0+2.0*hwt));
        sse2Start.x=ap[2].x;
        sse2Start.y=ap[2].y;
        sse2Start.z=ap[2].z;
        sse2End.x=ap[7].x;
        sse2End.y=ap[7].y;
        sse2End.z=ap[7].z;

        cenHelix1Temp=bf.addv(sse1Start,sse1End);
        cenHelix1=bf.scal(cenHelix1Temp,0.5);
        cenHelix2Temp=bf.addv(sse2Start,sse2End);
        cenHelix2=bf.scal(cenHelix2Temp,0.5);
        deltaHH=bf.minu(cenHelix2,cenHelix1);
        distanceSSE=bf.norm(deltaHH);
        if(distanceSSE>20) continue;

        axisZ=bf.unit(bf.minu(sse1End,sse1Start));
        axisY=bf.unit(bf.prod(axisZ,deltaHH));
        axisX=bf.unit(bf.prod(axisY,axisZ));
        vecHelix2=bf.minu(sse2End,sse2Start);
        lenHelixX=bf.dotv(vecHelix2,axisX);
        lenHelixY=bf.dotv(vecHelix2,axisY);
        lenHelixZ=bf.dotv(vecHelix2,axisZ);
        unitHH=bf.unit(deltaHH);
        lenHelix=bf.norm(vecHelix2);

        if(lenHelixX != 0 && lenHelixY != 0) 
        {
          phi=atan2(lenHelixY,lenHelixX)*(double)180.0/PI;
          if(phi<0) phi+=360;
        }
        else
        {
          continue;
        }

        normLenHelixZ=lenHelixZ/lenHelix;
        if(lenHelix!=0 && abs(normLenHelixZ)<=1)
        {
          psi=acos(normLenHelixZ)*(double)180.0/PI;
          if(psi<0) psi+=360;
        }
        else
        {
          continue;
        }

        cte=bf.unit(bf.minu(sse1End,cenHelix1));
        cte_dot=bf.dotv(cte,unitHH);
        theta=acos(cte_dot)*(double)180.0/PI;
        if(theta<0) theta+=180;

        if(seqSep<=5 && seqSep>1)
        {
          energyPhiPsiTheta+=energyHHPackAngle2[int(phi)][int(psi)][int(theta)];
        }
        else if(seqSep<=11)
        {
          energyPhiPsiTheta+=energyHHPackAngle3[int(phi)][int(psi)][int(theta)];
        }
        else
        {
          energyPhiPsiTheta+=energyHHPackAngle4[int(phi)][int(psi)][int(theta)];
        }

        if(seqSep<=11 && seqSep>1)
        {
          energyDistTheta+=energyHHPackDist2[int(distanceSSE*DIST_SCALE)][int(theta)];
        }
        else if(seqSep>=12)
        {
          energyDistTheta+=energyHHPackDist3[int(distanceSSE*DIST_SCALE)][int(theta)];
        }
      }
    }
  }
  totene[0]=energyPhiPsiTheta;
  totene[1]=energyDistTheta;

  return totene;
}

vector<double> FoldDesignEnergyFunction::calcHelixStrandEnergy(
  point3f *decstr,
  int numseq
)
{
  const int DIST_SCALE=10;
  int j,k;
  int seqSep;
  double hwt=2.879385;
  double phi,psi,theta,distanceSSE;
  double lenHelixX,lenHelixY,lenHelixZ,lenHelix,normLenHelixZ,cte_dot;
  point3d cenHelix,cenStrand,cenHelixTemp,cenStrandTemp,deltaHS;
  point3d sse1Start,sse1End,sse2Start,sse2End,ap[10];
  point3d axisX,axisY,axisZ,vecHelix,unitHS;
  point3d cte;
  double energyPhiPsiTheta=0;
  double energyDistTheta=0;
  vector<double> totene(2,0);
  sssegment *sse=inputInfo.getSSE();
  int numsse=inputInfo.getNumSSE();

  for(j=0;j<numsse-1;j++)
  {
    if(sse[j].ss=='H')
    {
      if(sse[j].term-sse[j].init<2) continue;
      ap[0]=bf.setv(decstr[sse[j].init].x,decstr[sse[j].init].y,decstr[sse[j].init].z);
      ap[1]=bf.setv(decstr[sse[j].init+1].x,decstr[sse[j].init+1].y,decstr[sse[j].init+1].z);
      ap[2]=bf.setv(decstr[sse[j].init+2].x,decstr[sse[j].init+2].y,decstr[sse[j].init+2].z);
      ap[3]=bf.addv(ap[0],ap[2]);
      ap[4]=bf.scal(ap[3],hwt);
      ap[0]=bf.addv(ap[4],ap[1]);
      ap[2]=bf.scal(ap[0],1.0/(1.0+2.0*hwt));
      ap[5]=bf.setv(decstr[sse[j].term].x,decstr[sse[j].term].y,decstr[sse[j].term].z);
      ap[6]=bf.setv(decstr[sse[j].term-1].x,decstr[sse[j].term-1].y,decstr[sse[j].term-1].z);
      ap[7]=bf.setv(decstr[sse[j].term-2].x,decstr[sse[j].term-2].y,decstr[sse[j].term-2].z);
      ap[8]=bf.addv(ap[5],ap[7]);
      ap[9]=bf.scal(ap[8],hwt);
      ap[5]=bf.addv(ap[9],ap[6]);
      ap[7]=bf.scal(ap[5],1.0/(1.0+2.0*hwt));
      sse1Start.x=ap[2].x;
      sse1Start.y=ap[2].y;
      sse1Start.z=ap[2].z;
      sse1End.x=ap[7].x;
      sse1End.y=ap[7].y;
      sse1End.z=ap[7].z;
      for(k=j+1;k<numsse;k++) if(sse[k].ss=='E')
      {
        if(sse[k].term-sse[k].init<1) continue;
        seqSep=sse[k].init-sse[j].term;
        if(seqSep<2) continue;

        ap[0]=bf.setv(decstr[sse[k].init].ptn.x,
                      decstr[sse[k].init].ptn.y,
                      decstr[sse[k].init].ptn.z);
        ap[1]=bf.setv(decstr[sse[k].term].ptc.x,
                      decstr[sse[k].term].ptc.y,
                      decstr[sse[k].term].ptc.z);
        sse2Start.x=ap[0].x;
        sse2Start.y=ap[0].y;
        sse2Start.z=ap[0].z;
        sse2End.x=ap[1].x;
        sse2End.y=ap[1].y;
        sse2End.z=ap[1].z;

        cenHelixTemp=bf.addv(sse1Start,sse1End);
        cenHelix=bf.scal(cenHelixTemp,0.5);
        cenStrandTemp=bf.addv(sse2Start,sse2End);
        cenStrand=bf.scal(cenStrandTemp,0.5);
        deltaHS=bf.minu(cenStrand,cenHelix);
        distanceSSE=bf.norm(deltaHS);
        if(distanceSSE>20) continue;

        axisZ=bf.unit(bf.minu(sse1End,sse1Start));
        axisY=bf.unit(bf.prod(axisZ,deltaHS));
        axisX=bf.unit(bf.prod(axisY,axisZ));
        vecHelix=bf.minu(sse2End,sse2Start);
        lenHelixX=bf.dotv(vecHelix,axisX);
        lenHelixY=bf.dotv(vecHelix,axisY);
        lenHelixZ=bf.dotv(vecHelix,axisZ);
        unitHS=bf.unit(deltaHS);
        lenHelix=bf.norm(vecHelix);

        if(lenHelixX != 0 && lenHelixY != 0)
        {
          phi=atan2(lenHelixY,lenHelixX)*(double)180.0/PI;
          if(phi<0) phi+=360;
        }
        else
        {
          continue;
        }

        normLenHelixZ=lenHelixZ/lenHelix;
        if(lenHelix!=0 && abs(normLenHelixZ)<=1)
        {
          psi=acos(normLenHelixZ)*(double)180.0/PI;
          if(psi<0) psi+=360;
        }
        else
        {
          continue;
        }
        cte=bf.unit(bf.minu(sse2End,cenStrand));
        cte_dot=bf.dotv(cte,unitHS);
        theta=acos(cte_dot)*(double)180.0/PI;
        if(theta<0) theta+=180;

        if(seqSep<=5 && seqSep>1)
        {
          energyPhiPsiTheta+=energyHSPackAngle2[int(phi)][int(psi)][int(theta)];
        }
        else if(seqSep<=11)
        {
          energyPhiPsiTheta+=energyHSPackAngle3[int(phi)][int(psi)][int(theta)];
        }
        else
        {
          energyPhiPsiTheta+=energyHSPackAngle4[int(phi)][int(psi)][int(theta)];
        }

        if(seqSep<=5 && seqSep>1)
        {
          energyDistTheta+=energyHSPackDist2[int(distanceSSE*DIST_SCALE)][int(theta)];
        }
        else if(seqSep<=11)
        {
          energyDistTheta+=energyHSPackDist3[int(distanceSSE*DIST_SCALE)][int(theta)];
        }
        else
        {
          energyDistTheta+=energyHSPackDist4[int(distanceSSE*DIST_SCALE)][int(theta)];
        }
      }
    }
    else if(sse[j].ss=='E')
    {
      if(sse[j].term-sse[j].init<1) continue;
      ap[0]=bf.setv(decstr[sse[j].init].ptn.x,
                    decstr[sse[j].init].ptn.y,
                    decstr[sse[j].init].ptn.z);
      ap[1]=bf.setv(decstr[sse[j].term].ptc.x,
                    decstr[sse[j].term].ptc.y,
                    decstr[sse[j].term].ptc.z);
      sse2Start.x=ap[0].x;
      sse2Start.y=ap[0].y;
      sse2Start.z=ap[0].z;
      sse2End.x=ap[1].x;
      sse2End.y=ap[1].y;
      sse2End.z=ap[1].z;

      for(k=j+1;k<numsse;k++) if(sse[k].ss=='H')
      {
        if(sse[k].term-sse[k].init<2) continue;
        seqSep=sse[k].init-sse[j].term;
        if(seqSep<2) continue;

        ap[0]=bf.setv(decstr[sse[k].init].x,decstr[sse[k].init].y,decstr[sse[k].init].z);
        ap[1]=bf.setv(decstr[sse[k].init+1].x,decstr[sse[k].init+1].y,decstr[sse[k].init+1].z);
        ap[2]=bf.setv(decstr[sse[k].init+2].x,decstr[sse[k].init+2].y,decstr[sse[k].init+2].z);
        ap[3]=bf.addv(ap[0],ap[2]);
        ap[4]=bf.scal(ap[3],hwt);
        ap[0]=bf.addv(ap[4],ap[1]);
        ap[2]=bf.scal(ap[0],1.0/(1.0+2.0*hwt));
        ap[5]=bf.setv(decstr[sse[k].term].x,decstr[sse[k].term].y,decstr[sse[k].term].z);
        ap[6]=bf.setv(decstr[sse[k].term-1].x,decstr[sse[k].term-1].y,decstr[sse[k].term-1].z);
        ap[7]=bf.setv(decstr[sse[k].term-2].x,decstr[sse[k].term-2].y,decstr[sse[k].term-2].z);
        ap[8]=bf.addv(ap[5],ap[7]);
        ap[9]=bf.scal(ap[8],hwt);
        ap[5]=bf.addv(ap[9],ap[6]);
        ap[7]=bf.scal(ap[5],1.0/(1.0+2.0*hwt));
        sse1Start.x=ap[2].x;
        sse1Start.y=ap[2].y;
        sse1Start.z=ap[2].z;
        sse1End.x=ap[7].x;
        sse1End.y=ap[7].y;
        sse1End.z=ap[7].z;

        cenHelixTemp=bf.addv(sse1Start,sse1End);
        cenHelix=bf.scal(cenHelixTemp,0.5);
        cenStrandTemp=bf.addv(sse2Start,sse2End);
        cenStrand=bf.scal(cenStrandTemp,0.5);
        deltaHS=bf.minu(cenHelix,cenStrand);
        distanceSSE=bf.norm(deltaHS);
        if(distanceSSE>20) continue;

        axisZ=bf.unit(bf.minu(sse1End,sse1Start));
        axisY=bf.unit(bf.prod(axisZ,deltaHS));
        axisX=bf.unit(bf.prod(axisY,axisZ));
        vecHelix=bf.minu(sse2End,sse2Start);
        lenHelixX=bf.dotv(vecHelix,axisX);
        lenHelixY=bf.dotv(vecHelix,axisY);
        lenHelixZ=bf.dotv(vecHelix,axisZ);
        unitHS=bf.unit(deltaHS);
        lenHelix=bf.norm(vecHelix);

        if(lenHelixX!=0 && lenHelixY!=0)
        {
          phi=atan2(lenHelixY,lenHelixX)*(double)180.0/PI;
          if(phi<0) phi+=360;
        }
        else
        {
          continue;
        }

        normLenHelixZ=lenHelixZ/lenHelix;
        if(lenHelix!=0 && abs(normLenHelixZ)<=1)
        {
          psi=acos(normLenHelixZ)*(double)180.0/PI;
          if(psi<0) psi+=360;
        }
        else
        {
          continue;
        }

        cte=bf.unit(bf.minu(sse2End,cenStrand));
        cte_dot=bf.dotv(cte,unitHS);
        theta=acos(cte_dot)*(double)180.0/PI;
        if(theta<0) theta+=180;

        if(seqSep<=5 && seqSep>1)
        {
          energyPhiPsiTheta+=energyHSPackAngle2[int(phi)][int(psi)][int(theta)];
        }
        else if(seqSep<=11)
        {
          energyPhiPsiTheta+=energyHSPackAngle3[int(phi)][int(psi)][int(theta)];
        }
        else
        {
          energyPhiPsiTheta+=energyHSPackAngle4[int(phi)][int(psi)][int(theta)];
        }

        if(seqSep<=5 && seqSep>1)
        {
          energyDistTheta+=energyHSPackDist2[int(distanceSSE*DIST_SCALE)][int(theta)];
        }
        else if(seqSep<=11)
        {
          energyDistTheta+=energyHSPackDist3[int(distanceSSE*DIST_SCALE)][int(theta)];
        }
        else
        {
          energyDistTheta+=energyHSPackDist4[int(distanceSSE*DIST_SCALE)][int(theta)];
        }
      }
    }
  }
  totene[0]=energyPhiPsiTheta;
  totene[1]=energyDistTheta;

  return totene;
}


void FoldDesignEnergyFunction::setHbondRamp(
  int curCycle,
  int numCycles
)
{
  const double CYCLE_SCALE=0.9;

  if(curCycle>int(numCycles*CYCLE_SCALE))
  {
    hbondRamp=true;
  }
  else
  {
    hbondRamp=false;
  }
}

double FoldDesignEnergyFunction::calcTotalEnergy(
  point3f *decstr,
  int numseq
)
{
  const int ALPHA_PROTEIN=1;
  double totalEnergy=0; // total energy E_total = sum of enelist
  int i;
  int zz;
  int protType=inputInfo.getProteinType();
  for(i=0;i<MAX_ENERGY_TERMS;i++) enelist[i]=0;
  geometry.tor2strsg2(decstr,numseq);

  //for(i=0;i<=20;i++) weights[i]=1.0;

  calcAllAtomEnergy(decstr,numseq,enelist);
  enelist[2]=0.03*calcHbondBackboneEnergy(decstr,numseq);
  enelist[6]=0.1*calcRamaEnergy(decstr,numseq);
  enelist[9]=1.0*calcBABMotifPenalty(decstr,numseq);
  enelist[16]=10.0*calcSSConstrEnergy(decstr,numseq);
  

  //-------Calculate secondary structure packing energy--------------->
  vector<double> energyHelixHelix=calcHelixHelixEnergy(decstr,numseq);
  vector<double> energyHelixStrand=calcHelixStrandEnergy(decstr,numseq);
  vector<double> energyStrandStrand=calcStrandStrandEnergy(decstr,numseq);
  enelist[10]=energyHelixHelix[0];
  enelist[11]=energyHelixHelix[1];
  enelist[12]=energyHelixStrand[0];
  enelist[13]=energyHelixStrand[1];
  enelist[14]=energyStrandStrand[0];
  enelist[15]=energyStrandStrand[1];

  //-------Reweight energy terms------------->
  totalEnergy+=0.4*enelist[0]*weights[0];
  totalEnergy+=9.0*enelist[1]*weights[1];  //Excluded volume

  //---Reweight backbone H-bond energy, depending on ramp parameter--->
  if(!hbondRamp) //Weight
  {
    totalEnergy+=8.0*enelist[2]*weights[2];  //Backbone Hbonds
  }
  else
  {
    totalEnergy+=9.0*enelist[2]*weights[3];  //Backbone Hbonds
  }
  totalEnergy+=enelist[4]*weights[4];  //Contact number
  totalEnergy+=10.0*enelist[5]*weights[5];  //Distance constraints
  totalEnergy+=0.65*enelist[6]*weights[6];  //Ramachandran energy
  totalEnergy+=enelist[7]*weights[7];  //Radius of gyration energy
  totalEnergy+=4.0*enelist[8]*weights[8];  //Distprof from fragments
  totalEnergy+=enelist[9]*weights[9];  //Beta-alpha-beta motif packing penalty
  if(protType!=ALPHA_PROTEIN)
  {
    totalEnergy+=0.75*enelist[10]*weights[10]; //Helix-helix packing angle energy
    totalEnergy+=0.5*enelist[11]*weights[11]; //Helix-helix packing distance energy
  }
  else
  {
    totalEnergy+=1.75*enelist[10]*weights[12]; //Helix-helix packing angle energy for alpha prot
    totalEnergy+=2.5*enelist[11]*weights[13]; //Helix-helix packing distance energy for alpha prot
  }
  totalEnergy+=1.0*enelist[12]*weights[14]; //Helix-strand packing angle energy
  totalEnergy+=0.5*enelist[13]*weights[15]; //Helix-strand packing distance energy
  totalEnergy+=0.1*enelist[14]*weights[16]; //Strand-strand packing angle energy
  totalEnergy+=2.0*enelist[15]*weights[17]; //Strand-strand packing distance energy
  totalEnergy+=4.0*enelist[16]*weights[18]; //Secondary structure restraint    
  totalEnergy+=0.0*enelist[17]*weights[19]; //Ca-Ca bond break
  totalEnergy+=2.0*enelist[18]*weights[20]; //Solvation energy
  totalEnergy+=1.0*enelist[19]; //User contact energy

/*
  cout << "-------------- Energy ------------" << endl;
  cout << " enelist 0 "<< enelist[0] << " weights "<< weights[0] << endl;
  cout << " enelist 1 "<< enelist[1] << " weights "<< weights[1] << endl;
  cout << " enelist 2 "<< enelist[2] << " weights "<< weights[2] << endl;
  cout << " enelist 3 "<< enelist[3] << " weights "<< weights[3] << endl;
  cout << " enelist 4 "<< enelist[4] << " weights "<< weights[4] << endl;
  cout << " enelist 5 "<< enelist[5] << " weights "<< weights[5] << endl;
  cout << " enelist 6 "<< enelist[6] << " weights "<< weights[6] << endl;
  cout << " enelist 7 "<< enelist[7] << " weights "<< weights[7] << endl;
  cout << " enelist 8 "<< enelist[8] << " weights "<< weights[8] << endl;
  cout << " enelist 9 "<< enelist[9] << " weights "<< weights[9] << endl;
  cout << " enelist 10 "<< enelist[10] << " weights "<< weights[10] << endl;
  cout << " enelist 11 "<< enelist[11] << " weights "<< weights[11] << endl;
  cout << " enelist 12 "<< enelist[12] << " weights "<< weights[12] << endl;
  cout << " enelist 13 "<< enelist[13] << " weights "<< weights[13] << endl;
  cout << " enelist 14 "<< enelist[14] << " weights "<< weights[14] << endl;
  cout << " enelist 15 "<< enelist[15] << " weights "<< weights[15] << endl;
  cout << " enelist 16 "<< enelist[16] << " weights "<< weights[16] << endl;
  cout << " enelist 17 "<< enelist[17] << " weights "<< weights[17] << endl;
  cout << " enelist 18 "<< enelist[18] << " weights "<< weights[18] << endl;
  cout << " enelist 19 "<< enelist[19] << " weights "<< weights[19] << endl;
  cout << "---------------------------------" << endl;
*/
  
  return totalEnergy;
}
   

