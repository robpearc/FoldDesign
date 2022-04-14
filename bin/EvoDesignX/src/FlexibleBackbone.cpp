///////////////////////////////////////////////////////////////////////////////////////
//Copyright (C) Yang Zhang Lab, University of Michigan - All Rights Reserved
//Unauthorized copying of this file, via any medium is strictly prohibited
//Written by Robin Pearce <robpearc@umich.edu> and Xiaoqiang Huang <xiaoqiah@umich.edu>
////////////////////////////////////////////////////////////////////////////////////////

#include "FlexibleBackbone.h"
#include "RotamerBuilder.h"
#include "RotamerOptimizer.h"
#include "Residue.h"
#include "Utility.h"
#include "kabsch.h"
#include <string.h>
#include <iostream>

extern BOOL FLAG_EXPAND_HYDROXYL_ROT;
extern BOOL FLAG_PPI;
extern BOOL FLAG_PROT_LIG;
extern BOOL FLAG_ENZYME;
extern char DES_CHAINS[10];
static long seed[STREAMS] = {DEFAULT};  /* current state of each stream   */
static int  stream        = 0;          /* stream index, 0 is the default */

double Random(void)
{
  const long Q = MODULUS / MULTIPLIER;
  const long R = MODULUS % MULTIPLIER;
  long t;

  t = MULTIPLIER * (seed[stream] % Q) - R * (seed[stream] / Q);
  if (t > 0) seed[stream] = t;
  else seed[stream] = t + MODULUS;
  return ((double) seed[stream] / MODULUS);
}

/* Sphere surface semi-random sampling for given theta (angle between
 *  *  *  * z-axis and the sampled vector) and random Phi */
XYZ rana(double theta){
  XYZ tp;
  double dphi;
  dphi=(2.0*Random())*PI;
  tp.X=sin(theta)*cos(dphi);
  tp.Y=sin(theta)*sin(dphi);
  tp.Z=cos(theta);
  return tp;
}

XYZ ranv_ergodic(double fac){
  XYZ tp;
  double dr=fabs(fac);               // sample on surface
  if (fac>0) dr=0.00001+fac*Random();// sample inside sphere 
  double v=(2.*Random()-1.);         // [-1, 1)
  double dtheta=acos(v);             // [0, PI)
  double dphi=(2.*Random()-1.)*PI;   // [-PI,PI)
  tp.X=dr*sin(dtheta)*cos(dphi);
  tp.Y=dr*sin(dtheta)*sin(dphi);
  tp.Z=dr*cos(dtheta);
  return tp;
}

int LoadRamaOutlier(double*** rlogduke,int*** ramaduke,char* ramafile){
  FILE* fin=fopen(ramafile,"rb");
  if(fin==NULL){
    return IOError;
  }
  for(int i=0;i<4;i++){
    for(int j=0;j<360;j++){
      fread(ramaduke[i][j],sizeof(int),360,fin);
    }
  }
  for(int i=0;i<4;i++){
    for(int j=0;j<360;j++){
      fread(rlogduke[i][j],sizeof(double),360,fin);
    }
  }
  fclose(fin);
  return Success;
}

int LoadTorDistribution(double** tordis,double* accdis,char* torDisFile){
  FILE* fin=fopen(torDisFile,"rt");
  if(fin==NULL){
    return IOError;
  }
  
  double lineprob[360];
  for(int i=0;i<360;i++){
    lineprob[i]=0;
  }

  char oneline[255];
  float tval;
  int k=0;
  for(int i=0;i<360;i++){
    for(int j=0;j<360;j++){
      fgets(oneline,255,fin);
      sscanf(oneline,"%f",&tval);
      tordis[i][j]=double(tval);
      lineprob[i]+=tordis[i][j];
    }
  }
  fclose(fin);
 
  for(int i=0;i<360;i++){
    for(int j=0;j<360;j++){
      if(k==0) accdis[k]=tordis[i][j];
      else accdis[k]=accdis[k-1]+tordis[i][j];
      k++;
    }
  }
  accdis[k-1]=1.01;
  return Success;
}

int LoadPhiPsiProb(double ***phipsiprob,char *torFileName){       
  FILE* fin=fopen(torFileName,"rt");
  if(fin==NULL){
    return IOError;
  }
   
  char oneline[200];
  for(int i=0;i<8;i++){
    for(int j=0;j<20;j++){
      for(int k=0;k<360;k++){
        fgets(oneline,200,fin);
        sscanf(oneline,"%lf",&phipsiprob[i][j][k]);
      }
    }
  }
  fclose(fin);

  for(int i=0;i<8;i++){       
    for(int j=0;j<20;j++){       
      for(int k=1;k<360;k++){       
        phipsiprob[i][j][k]+=phipsiprob[i][j][k-1];
      }
    }
  }
  return Success;
}

int Loadca2ncbins(float ****cancbins,char *filename){       
  FILE* fin=fopen(filename,"rt");
  if(fin==NULL){
    return IOError;
  }

  int i,j,k,l;
  int ii,jj,kk;
  char oneline[300];
  for(j=0;j<105;j++){               
    for(k=0;k<105;k++){       
      for(l=0;l<160;l++){       
        fgets(oneline,300,fin);
        sscanf(oneline,"%d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f",
                        &ii,&jj,&kk,&cancbins[0][j][k][l],&cancbins[1][j][k][l],&cancbins[2][j][k][l],&cancbins[3][j][k][l],
                        &cancbins[4][j][k][l],&cancbins[5][j][k][l],&cancbins[6][j][k][l],&cancbins[7][j][k][l],&cancbins[8][j][k][l],
                        &cancbins[9][j][k][l],&cancbins[10][j][k][l],&cancbins[11][j][k][l],&cancbins[12][j][k][l]);
      }
    }
  }
  fclose(fin);
  return Success;
}

int Loadbbord(double *bborder,char *filename){
  FILE* fin=fopen(filename,"rt");
  if(fin==NULL){
    return IOError;
  }

  int i,j;
  float tval;
  char oneline[100];
  fgets(oneline,100,fin);
  for(i=0;i<1000;i++){
    fgets(oneline,100,fin);
    sscanf(oneline,"%d %f",&j,&tval);
    bborder[i]=tval*30.0;
  }
  fclose(fin);
  return Success;
}

int LoadSecondaryStructureDesign(ChainFlex *pChainFlex,sssegment *sse,int &numsse,char *ssFile){
  FILE* fin=fopen(ssFile,"rt");
  if(fin==NULL){
    return IOError;
  }
  int i,j,k;
  char tmpstring[255];
  char tmp1;
  const int numseq=ChainFlexGetResidueCount(pChainFlex);
  float a0[numseq],a1[numseq],a2[numseq];
  fgets(tmpstring, 255, fin);
  for(i=0;i<ChainFlexGetResidueCount(pChainFlex);i++){
    ResidueFlex *pResidueFlex = ChainFlexGetResidue(pChainFlex,i);
    fgets(tmpstring, 255, fin);
    sscanf(tmpstring+5,"%c %c %f %f %f",&tmp1,&pResidueFlex->ss,&a0[i],&a1[i],&a2[i]);
  }
  fclose(fin);

  for(i=1;i<numseq-1;i++){
    ResidueFlex *pResidueFlexPrev = ChainFlexGetResidue(pChainFlex,i-1);
    ResidueFlex *pResidueFlex = ChainFlexGetResidue(pChainFlex,i);
    ResidueFlex *pResidueFlexNext = ChainFlexGetResidue(pChainFlex,i+1);

    if(pResidueFlex->ss=='H' && pResidueFlexPrev->ss=='E' && pResidueFlexNext->ss=='E' && a1[i]!=0.999 && a0[i]<a2[i]){
      pResidueFlex->ss='E';
    }
    else if(pResidueFlex->ss=='H' && pResidueFlexPrev->ss=='E' && pResidueFlexNext->ss=='E' && a1[i]!=0.999){
      pResidueFlex->ss='C';
    }
  }
  for(i=0;i<numseq;i++){
    ResidueFlex *pResidueFlexi = ChainFlexGetResidue(pChainFlex,i);
    j=i;
    ResidueFlex *pResidueFlexj = ChainFlexGetResidue(pChainFlex,j);
    while(j<numseq && pResidueFlexj->ss==pResidueFlexi->ss){
      j++;
      pResidueFlexj=ChainFlexGetResidue(pChainFlex,j);
    }
    if((j-i==2 || j-i==1) && pResidueFlexi->ss=='H'){
      for(k=i;k<j;k++){
        ResidueFlex *pResidueFlexk = ChainFlexGetResidue(pChainFlex,k);
        pResidueFlexk->ss='C';
      }
    }
    i=j-1;
  }

  for(i=0;i<numseq;i++){
    ResidueFlex *pResidueFlex1 = ChainFlexGetResidue(pChainFlex,i);
    ResidueFlex *pResidueFlex2 = ChainFlexGetResidue(pChainFlex,i);
    char s1=pResidueFlex1->ss;
    char s2=pResidueFlex2->ss;
    j=i;
    while(j<numseq && s1==s2){
      s1=pResidueFlex2->ss;
      j++;
      pResidueFlex2 = ChainFlexGetResidue(pChainFlex,j);
    }

    sse[numsse].init=i;
    if(j>=numseq){
      sse[numsse].term=j-1;
    }
    else{
      sse[numsse].term=j-2;
    }
    sse[numsse].ss=s2;
    numsse++;
    if(j>=numseq){
      break;
    }
    else{
      i=j-2;
    }
  }
  if(numsse!=0) sse=(sssegment *)realloc(sse,numsse*sizeof(sssegment));

  return Success;
}

void calcbbptable(ChainFlex *pChainFlex,sssegment *sse,int &numbbptb,int numsse,double *bborder,paircont *bbptb,int *alphasind,int &numsalpha){
  int i,j,ii,jj;
  int numbeta=0;
  int numalpha=0;
  numsalpha=0;

  int *betaind=new int[numsse];
  int *alphaind=new int[numsse];
  //alphasind=new int[numsse];
  for(i=0;i<numsse;i++){
    if(sse[i].ss=='H') alphaind[numalpha++]=i;
    else if(sse[i].ss=='E') betaind[numbeta++]=i;
  }

  for(i=0;i<numalpha-1;i++){
    if(alphaind[i+1]==alphaind[i]+2) alphasind[numsalpha++]=alphaind[i];
  }

  for(i=0;i<numbbptb;i++){
    bbptb[i].init=-1;
    bbptb[i].term=-1;
    bbptb[i].dist=0;
  }
  //depth for beta
  double *betadep=new double[ChainFlexGetResidueCount(pChainFlex)];
  for(i=0;i<ChainFlexGetResidueCount(pChainFlex);i++) betadep[i]=0;
  for(i=0;i<numsse;i++){
    if(sse[i].ss=='C')
      for(j=sse[i].init;j<=sse[i].term;j++) betadep[j]=0.6;
    else if(sse[i].ss=='H')
      for(j=sse[i].init;j<=sse[i].term;j++) betadep[j]=0.0;//0.2
  }
  for(i=0;i<numbeta;i++){
    double tlen=(sse[betaind[i]].term+1-sse[betaind[i]].init)/2.0;
    double tmid=(sse[betaind[i]].term+sse[betaind[i]].init)/2.0;
    for(j=sse[betaind[i]].init;j<=sse[betaind[i]].term;j++) betadep[j]=tlen-fabs(tmid-j)+6.5;
  }
  //start
  int totnum=0;
  for(i=0;i<numsse;i++){
    if(sse[i].ss=='H') continue;
    for(ii=sse[i].init;ii<=sse[i].term;ii++){
      for(j=0;j<numsse;j++){
        if(sse[j].ss=='H') continue;
        for(jj=sse[j].init;jj<=sse[j].term;jj++){
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
          else if(sse[i].ss=='E' && sse[j].ss=='E'){
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
  for(i=0;i<numbbptb;i++){
    flaguse[i]=true;
    if(bbptb[i].dist<err || bbptb[i].init>bbptb[i].term || bbptb[i].init==0 || bbptb[i].term==ChainFlexGetResidueCount(pChainFlex)-1)
      flaguse[i]=false;
  }
  j=0;
  for(i=0;i<numbbptb;i++){
    if(flaguse[i]){
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
  for(i=0;i<numbbptb;i++){
    bbptb[i].dist/=totsum;
    bbptb[i].prob=bbptb[i].dist;
  }
  for(i=1;i<numbbptb;i++) bbptb[i].dist+=bbptb[i-1].dist;
}

int loadContact(double **contactCA,double **contactCB,int numseq,char *filename){
  int i,j;
  for(i=0;i<numseq;i++){
    for(j=0;j<numseq;j++){
      contactCA[i][j]=0;   
      contactCB[i][j]=0;   
    }
  }   

  FILE* fin=fopen(filename,"rt");
  if(fin==NULL){
    return IOError;
  }
  
  int istart,iend,atom_type;
  double ialarm;
  char oneline[255];

  fgets(oneline,255,fin); // read file head
  while(!feof(fin)){   
    fgets(oneline,255,fin);
    sscanf(oneline,"%d %d %lf %d",&istart,&iend,&ialarm,&atom_type);
    if(atom_type==1) contactCA[istart-1][iend-1]=ialarm;  //C-alpha
    else contactCB[istart-1][iend-1]=ialarm; //C-beta
  }
  fclose(fin);
  return Success;
}

int loadDistance(double **distanceCB,double **distanceWeight,int numseq,char *filename){
  int i,j;
  for(i=0;i<numseq;i++){
    for(j=0;j<numseq;j++){
      distanceCB[i][j]=0;   
      distanceWeight[i][j]=0;   
    }
  }   

  FILE* fin=fopen(filename,"rt");
  if(fin==NULL){
    return IOError;
  }
  
  int istart,iend,atom_type;
  double weight,dist;
  char oneline[255];

  fgets(oneline,255,fin); // read file head
  while(!feof(fin)){   
    fgets(oneline,255,fin);
    sscanf(oneline,"%d %d %lf %d",&istart,&iend,&weight,&dist);
    distanceCB[istart-1][iend-1]=dist; //C-beta
    distanceWeight[istart-1][iend-1]=weight;
  }
  fclose(fin);
  return Success;
}

double calcDistanceEnergy(Chain *pChain,double **distanceCB,double **distanceWeight){
  int numseq = ChainGetResidueCount(pChain);  
  int i,j;
  double cbDistEnergy=0;

  for(i=0;i<numseq;i++){
    for(j=0;j<numseq;j++){
      if(distanceWeight[i][j]>0.0){
        Residue *pResidue1 = ChainGetResidue(pChain,i);
        Residue *pResidue2 = ChainGetResidue(pChain,j);
        Atom *pAtomCB1 = ResidueGetAtomByName(pResidue1, "CB");
        Atom *pAtomCB2 = ResidueGetAtomByName(pResidue2, "CB");
        XYZ distVec = XYZDifference(&pAtomCB1->xyz,&pAtomCB2->xyz);
        double distance = XYZNormalization(&distVec);
        cbDistEnergy += -1*distanceWeight[i][j]*(1/(1+(distance-distanceCB[i][j])*(distance-distanceCB[i][j])));
      }   
    }
  }
  return cbDistEnergy;
}

void calcContactEnergy(Chain *pChain,double **contactCA,double **contactCB,double &caEnergy,double &cbEnergy){
  int i,j;
  XYZ distVec;
  double distance,weight;
  double d8,dwell,d10,da,db,dc,dd;
  caEnergy=0,cbEnergy=0;
  int numseq = ChainGetResidueCount(pChain);  

  if(numseq<100){
    dwell=6;
  }
  else if(numseq<120){
    dwell=8;
  }
  else if(numseq<200){
    dwell=10;
  }
  else{
    dwell=12;
  }

  d8=8.0;
  d10=d8+dwell;
  da=(d8+d10)/2; // middle of d8 and d10
  db=dwell; //width of first well
  dc=(d10+80)/2; // middle of d10 and 80
  dd=80-d10; // width of second well

  for(i=0;i<numseq;i++){
    for(j=0;j<numseq;j++){
      if(contactCA[i][j]>0.0){
        Residue *pResidue1 = ChainGetResidue(pChain,i);
        Residue *pResidue2 = ChainGetResidue(pChain,j);
        Atom *pAtomCA1 = ResidueGetAtomByName(pResidue1, "CA");
        Atom *pAtomCA2 = ResidueGetAtomByName(pResidue2, "CA");
        distVec = XYZDifference(&pAtomCA1->xyz,&pAtomCA2->xyz);
        distance = XYZNormalization(&distVec);
        weight=contactCA[i][j];    //weight for well depth
        if(distance<=d8) //r<8A
          caEnergy+=-weight;
        else if(distance<d10) // 8<r<10
          caEnergy+=-weight*(1-sin((distance-da)/db*PI))/2;
        else if(distance<80) //10<r<80
          caEnergy+=weight*(1+sin((distance-dc)/dd*PI))/2;
        else //r>80
          caEnergy+=weight;
      }   

      if(contactCB[i][j]>0.0){
        Residue *pResidue1 = ChainGetResidue(pChain,i);
        Residue *pResidue2 = ChainGetResidue(pChain,j);
        Atom *pAtomCB1 = ResidueGetAtomByName(pResidue1, "CB");
        Atom *pAtomCB2 = ResidueGetAtomByName(pResidue2, "CB");
        if(pAtomCB1!=0 && pAtomCB2!=0){
          distVec = XYZDifference(&pAtomCB1->xyz,&pAtomCB2->xyz);
          distance = XYZNormalization(&distVec);
          weight=contactCB[i][j];    //weight for well depth
          if(distance<=d8) //r<8A
            cbEnergy+=-weight;
          else if(distance<d10) // 8<r<10
            cbEnergy+=-weight*(1-sin((distance-da)/db*PI))/2;
          else if(distance<80) //10<r<80
            cbEnergy+=weight*(1+sin((distance-dc)/dd*PI))/2;
          else //r>80
            cbEnergy+=weight;
        }
      }   
    }
  }
}

void swap(int *a, int *b){ 
  int temp = *a; 
  *a = *b; 
  *b = temp; 
} 

void randomizeArray(int arr[], int numDesignSites){
  srand ( time(NULL) ); 
 
  for(int i=numDesignSites-1; i>0; i--){
    int j = rand() % (i+1);
    swap(&arr[i], &arr[j]); 
  } 
}

int getResidueIndex(Residue *newResi){
  int residueIndex;
 
  if(strcmp(ResidueGetName(newResi),"ALA")==0){
    residueIndex = 0;
  }
  else if(strcmp(ResidueGetName(newResi),"CYS")==0){
    residueIndex = 1;
  }
  else if(strcmp(ResidueGetName(newResi),"ASP")==0){
    residueIndex = 2;
  }
  else if(strcmp(ResidueGetName(newResi),"GLU")==0){
    residueIndex = 3;
  }
  else if(strcmp(ResidueGetName(newResi),"PHE")==0){
    residueIndex = 4;
  }
  else if(strcmp(ResidueGetName(newResi),"GLY")==0){
    residueIndex = 5;
  }
  else if(strcmp(ResidueGetName(newResi),"HIS")==0 || 
          strcmp(ResidueGetName(newResi),"HSE")==0 || 
          strcmp(ResidueGetName(newResi),"HSD")==0){
    residueIndex = 6;
  }
  else if(strcmp(ResidueGetName(newResi),"ILE")==0){
    residueIndex = 7;
  }
  else if(strcmp(ResidueGetName(newResi),"LYS")==0){
    residueIndex = 8;
  }
  else if(strcmp(ResidueGetName(newResi),"LEU")==0){
    residueIndex = 9;
  }
  else if(strcmp(ResidueGetName(newResi),"MET")==0){
    residueIndex = 10;
  }
  else if(strcmp(ResidueGetName(newResi),"ASN")==0){
    residueIndex = 11;
  }
  else if(strcmp(ResidueGetName(newResi),"PRO")==0){
    residueIndex = 12;
  }
  else if(strcmp(ResidueGetName(newResi),"GLN")==0){
    residueIndex = 13;
  }
  else if(strcmp(ResidueGetName(newResi),"ARG")==0){
    residueIndex = 14;
  }
  else if(strcmp(ResidueGetName(newResi),"SER")==0){
    residueIndex = 15;
  }
  else if(strcmp(ResidueGetName(newResi),"THR")==0){
    residueIndex = 16;
  }
  else if(strcmp(ResidueGetName(newResi),"VAL")==0){
    residueIndex = 17;
  }
  else if(strcmp(ResidueGetName(newResi),"TRP")==0){
    residueIndex = 18;
  }
  else if(strcmp(ResidueGetName(newResi),"TYR")==0){
    residueIndex = 19;
  }
  return residueIndex;
}

float calcTorsion(float xi,float yi,float zi,
                  float xj,float yj,float zj,
                  float xk,float yk,float zk,
                  float xl,float yl,float zl)
{
  double xij,yij,zij,
  xkj,ykj,zkj,
  xkl,ykl,zkl,
  dxi,dyi,dzi,
  gxi,gyi,gzi,
  bi,bk,ct,
  boi2,boj2,
  z1,z2,s,
  bioj,bjoi;
  float ap;

  /* Calculate the vectors C,B,C                                       */
  xij = xi - xj;
  yij = yi - yj;
  zij = zi - zj;
  xkj = xk - xj;
  ykj = yk - yj;
  zkj = zk - zj;
  xkl = xk - xl;
  ykl = yk - yl;
  zkl = zk - zl;

  /* Calculate the normals to the two planes n1 and n2
     this is given as the cross products:
      AB x BC
     --------- = n1
     |AB x BC|
                   
      BC x CD
     --------- = n2
     |BC x CD|
  */
  dxi = yij * zkj - zij * ykj;     /* Normal to plane 1                */
  dyi = zij * xkj - xij * zkj;
  dzi = xij * ykj - yij * xkj;
  gxi = zkj * ykl - ykj * zkl;     /* Mormal to plane 2                */
  gyi = xkj * zkl - zkj * xkl;
  gzi = ykj * xkl - xkj * ykl;

  /* Calculate the length of the two normals                           */
  bi = dxi * dxi + dyi * dyi + dzi * dzi;
  bk = gxi * gxi + gyi * gyi + gzi * gzi;
  ct = dxi * gxi + dyi * gyi + dzi * gzi;

  boi2 = 1./bi;
  boj2 = 1./bk;
  bi   = (double)sqrt((double)bi);
  bk   = (double)sqrt((double)bk);
  if(bi<err*0.01 || bk<err*0.01) return 180;
  z1   = 1./bi;
  z2   = 1./bk;
  bioj = bi * z2;
  bjoi = bk * z1;
  ct   = ct * z1 * z2;
  if (ct >  1.0)   ct = 1.0;
  if (ct < (-1.0)) ct = -1.0;
  ap   = acos(ct);

  s = xkj * (dzi * gyi - dyi * gzi)
          + ykj * (dxi * gzi - dzi * gxi)
          + zkj * (dyi * gxi - dxi * gyi);

  if (s < 0.0) ap = -ap;

  ap = (ap > 0.0) ? PI-ap : -(PI+ap);

  // angle
  ap *= (float)180.0/PI;
  if(ap<0) ap+=360;
  return ap;
}

void rotate_trans_chain(Chain* pChainNew, Chain* pChainOrig,int chain1_len,int chain2_len,int start, int end, int type)
{
  gsl_matrix *X;
  gsl_matrix *Y;

  //TODO change allocation here
  gsl_matrix *chain1_X = gsl_matrix_alloc(chain1_len*3,3);
  gsl_matrix *chain1_Y = gsl_matrix_alloc(chain1_len*3,3);
  gsl_matrix *chain2_X = gsl_matrix_alloc(chain2_len*3,3);
  gsl_matrix *chain2_Y = gsl_matrix_alloc(chain2_len*3,3);
  gsl_vector *translation = gsl_vector_alloc(3);
  gsl_matrix *U = gsl_matrix_alloc(3,3);

  //type 3 = chain1
  if(type==3){
    X = chain1_X;
    Y = chain1_Y;
  }
  //type 32 = chain2
  else if(type ==32){
    X = chain2_X;
    Y = chain2_Y;
  }
  //read data into X
  int i;
  int j=0;
  double x_center = 0;
  double y_center = 0;
  double z_center = 0;
  for(i=start; i<end; i++){
    //set matrix X.  This one has random movement.
    Residue* newResi = ChainGetResidue(pChainNew,i);
    Atom* pAtomN = ResidueGetAtomByName(newResi, "N");
    Atom* pAtomCA = ResidueGetAtomByName(newResi, "CA");
    Atom* pAtomC = ResidueGetAtomByName(newResi, "C");

    Residue* origResi = ChainGetResidue(pChainOrig,i);
    Atom* pAtomNorig = ResidueGetAtomByName(origResi, "N");
    Atom* pAtomCAorig = ResidueGetAtomByName(origResi, "CA");
    Atom* pAtomCorig = ResidueGetAtomByName(origResi, "C");


    gsl_matrix_set(X,j,0,pAtomN->xyz.X);
    gsl_matrix_set(X,j,1,pAtomN->xyz.Y);
    gsl_matrix_set(X,j,2,pAtomN->xyz.Z);

    gsl_matrix_set(X,j+1,0,pAtomCA->xyz.X);
    gsl_matrix_set(X,j+1,1,pAtomCA->xyz.Y);
    gsl_matrix_set(X,j+1,2,pAtomCA->xyz.Z);

    gsl_matrix_set(X,j+2,0,pAtomC->xyz.X);
    gsl_matrix_set(X,j+2,1,pAtomC->xyz.Y);
    gsl_matrix_set(X,j+2,2,pAtomC->xyz.Z);

    //calculate center of mass
    x_center = pAtomN->xyz.X + pAtomCA->xyz.X + pAtomC->xyz.X;
    y_center = pAtomN->xyz.Y + pAtomCA->xyz.Y + pAtomC->xyz.Y;
    z_center = pAtomN->xyz.Z + pAtomCA->xyz.Z + pAtomC->xyz.Z;

    //set matrix Y.  This one is before movement.
    gsl_matrix_set(Y,j,0,pAtomNorig->xyz.X);
    gsl_matrix_set(Y,j,1,pAtomNorig->xyz.Y);
    gsl_matrix_set(Y,j,2,pAtomNorig->xyz.Z);

    gsl_matrix_set(Y,j+1,0,pAtomCAorig->xyz.X);
    gsl_matrix_set(Y,j+1,1,pAtomCAorig->xyz.Y);
    gsl_matrix_set(Y,j+1,2,pAtomCAorig->xyz.Z);

    gsl_matrix_set(Y,j+2,0,pAtomCorig->xyz.X);
    gsl_matrix_set(Y,j+2,1,pAtomCorig->xyz.Y);
    gsl_matrix_set(Y,j+2,2,pAtomCorig->xyz.Z);

    j=j+3;
  }
  //kabsch transform messes up X and Y vector.  Rotation and translation seeem right though.
  //need to start from original x, remove center of mass and use U matrix and t translation 
  //to put chain in proper place.
  x_center = x_center/j;
  y_center = y_center/j;
  z_center = z_center/j;
  double *s=NULL;
  gsl_matrix_set_zero(U);
  int out = kabsch(end-start, X,Y,U, translation,s);

  j=0;
  for(i=start; i<end; i++){
    Residue* newResi = ChainGetResidue(pChainNew,i);
    Atom* pAtomN = ResidueGetAtomByName(newResi, "N");
    Atom* pAtomCA = ResidueGetAtomByName(newResi, "CA");
    Atom* pAtomC = ResidueGetAtomByName(newResi, "C");

    gsl_matrix_set(X,j,0,pAtomN->xyz.X - x_center);
    gsl_matrix_set(X,j,1,pAtomN->xyz.Y - y_center);
    gsl_matrix_set(X,j,2,pAtomN->xyz.Z - z_center);

    gsl_matrix_set(X,j+1,0,pAtomCA->xyz.X - x_center);
    gsl_matrix_set(X,j+1,1,pAtomCA->xyz.Y - y_center);
    gsl_matrix_set(X,j+1,2,pAtomCA->xyz.Z - z_center);

    gsl_matrix_set(X,j+2,0,pAtomC->xyz.X - x_center);
    gsl_matrix_set(X,j+2,1,pAtomC->xyz.Y - y_center);
    gsl_matrix_set(X,j+2,2,pAtomC->xyz.Z - z_center);

    j=j+3;
  }

  double x,y,z,tx,ty,tz;
  tx = round_nplaces(gsl_vector_get(translation,0),3);
  ty = round_nplaces(gsl_vector_get(translation,1),3);
  tz = round_nplaces(gsl_vector_get(translation,2),3);

  double rx1,rx2,rx3,ry1,ry2,ry3,rz1,rz2,rz3;

  //row 1
  rx1 = round_nplaces(gsl_matrix_get (U, 0, 0),3);
  ry1 = round_nplaces(gsl_matrix_get (U, 0, 1),3);
  rz1 = round_nplaces(gsl_matrix_get (U, 0, 2),3);
  //row 2
  rx2 = round_nplaces(gsl_matrix_get (U, 1, 0),3);
  ry2 = round_nplaces(gsl_matrix_get (U, 1, 1),3);
  rz2 = round_nplaces(gsl_matrix_get (U, 1, 2),3);
  //row 3
  rx3 = round_nplaces(gsl_matrix_get (U, 2, 0),3);
  ry3 = round_nplaces(gsl_matrix_get (U, 2, 1),3);
  rz3 = round_nplaces(gsl_matrix_get (U, 2, 2),3);
  //rotation translation

  double ptn[3];
  double ptca[3];
  double ptc[3];
  int count = 0;

  for(i=start; i<end; i++){
    Residue* newResi = ChainGetResidue(pChainNew,i);
    Atom* pAtomN = ResidueGetAtomByName(newResi, "N");
    Atom* pAtomCA = ResidueGetAtomByName(newResi, "CA");
    Atom* pAtomC = ResidueGetAtomByName(newResi, "C");

    for(j=0;j<3;j++){
      ptn[j] = gsl_matrix_get(X, count, j);
      ptca[j] = gsl_matrix_get(X, count+1, j);
      ptc[j] = gsl_matrix_get(X, count+2, j);
    }
    count = count + 3;

    x = ptca[0]*rx1 + ptca[1]*ry1 + ptca[2]*rz1 +tx;
    y = ptca[0]*rx2 + ptca[1]*ry2 + ptca[2]*rz2 +ty;
    z = ptca[0]*rx3 + ptca[1]*ry3 + ptca[2]*rz3 +tz;
    pAtomCA->xyz.X=x;
    pAtomCA->xyz.Y=y;
    pAtomCA->xyz.Z=z;

    x = ptc[0]*rx1 + ptc[1]*ry1 + ptc[2]*rz1 +tx;
    y = ptc[0]*rx2 + ptc[1]*ry2 + ptc[2]*rz2 +ty;
    z = ptc[0]*rx3 + ptc[1]*ry3 + ptc[2]*rz3 +tz;
    pAtomC->xyz.X=x;
    pAtomC->xyz.Y=y;
    pAtomC->xyz.Z=z;

    x = ptn[0]*rx1 + ptn[1]*ry1 + ptn[2]*rz1 +tx;
    y = ptn[0]*rx2 + ptn[1]*ry2 + ptn[2]*rz2 +ty;
    z = ptn[0]*rx3 + ptn[1]*ry3 + ptn[2]*rz3 +tz;
    pAtomN->xyz.X=x;
    pAtomN->xyz.Y=y;
    pAtomN->xyz.Z=z;
  }
  return;
}

int str2tor(Chain *pChain,ChainFlex *pChainFlex){
  int i;
  XYZ p12,p23; 
  double lencn = 1.338,lennca = 1.460,lencac = 1.525,angncac = 111.008;
  double raddeg = PI/180.0,angcacn = 116.617,angcnca = 121.614,lennc = 2.460;

  Residue *pResi0 = ChainGetResidue(pChain,0);
  ResidueFlex *pResiFlex0 = ChainFlexGetResidue(pChainFlex,0);
  Atom* pAtomN0 = ResidueGetAtomByName(pResi0, "N");
  Atom* pAtomCA0 = ResidueGetAtomByName(pResi0, "CA");
  Atom* pAtomC0 = ResidueGetAtomByName(pResi0, "C");

  if(pResiFlex0->dih[0]<0 || pResiFlex0->len[0]<0){
    pResiFlex0->dih[0]=180.0;
    pResiFlex0->dih[1]=180.0;
    pResiFlex0->dih[2]=180.0;
    pResiFlex0->len[0]=lennc;
    pResiFlex0->ang[0]=angcacn;
    pResiFlex0->ang[1]=angcnca;
  }

  p12 = XYZDifference(&pAtomCA0->xyz,&pAtomN0->xyz);
  p23 = XYZDifference(&pAtomCA0->xyz,&pAtomC0->xyz);

  pResiFlex0->len[1]=XYZNormalization(&p12);
  pResiFlex0->len[2]=XYZNormalization(&p23);
  pResiFlex0->ang[2] = RadToDeg(XYZAngle(&p12, &p23));

  for(i=1;i<ChainGetResidueCount(pChain);i++){
    Residue *pResiPrev = ChainGetResidue(pChain,i-1);
    Residue *pResi = ChainGetResidue(pChain,i);
    ResidueFlex *pResiFlex = ChainFlexGetResidue(pChainFlex,i);
    Atom* pAtomNprev = ResidueGetAtomByName(pResiPrev, "N");
    Atom* pAtomCAprev = ResidueGetAtomByName(pResiPrev, "CA");
    Atom* pAtomCprev = ResidueGetAtomByName(pResiPrev, "C");
    Atom* pAtomN = ResidueGetAtomByName(pResi, "N");
    Atom* pAtomCA = ResidueGetAtomByName(pResi, "CA");
    Atom* pAtomC = ResidueGetAtomByName(pResi, "C");

    p12 = XYZDifference(&pAtomCprev->xyz,&pAtomCAprev->xyz); 
    p23 = XYZDifference(&pAtomCprev->xyz,&pAtomN->xyz);
    pResiFlex->len[0] = XYZNormalization(&p23);
    pResiFlex->ang[0] = RadToDeg(XYZAngle(&p12, &p23));
    pResiFlex->dih[0] = calcTorsion(pAtomNprev->xyz.X,pAtomNprev->xyz.Y,pAtomNprev->xyz.Z,
                                    pAtomCAprev->xyz.X,pAtomCAprev->xyz.Y,pAtomCAprev->xyz.Z,
                                    pAtomCprev->xyz.X,pAtomCprev->xyz.Y,pAtomCprev->xyz.Z,
                                    pAtomN->xyz.X,pAtomN->xyz.Y,pAtomN->xyz.Z);

    p12 = XYZDifference(&pAtomN->xyz,&pAtomCprev->xyz);    
    p23 = XYZDifference(&pAtomN->xyz,&pAtomCA->xyz);    
    pResiFlex->len[1] = XYZNormalization(&p23);
    pResiFlex->ang[1] = RadToDeg(XYZAngle(&p12,&p23));
    pResiFlex->dih[1] = calcTorsion(pAtomCAprev->xyz.X,pAtomCAprev->xyz.Y,pAtomCAprev->xyz.Z,
                                    pAtomCprev->xyz.X,pAtomCprev->xyz.Y,pAtomCprev->xyz.Z,
                                    pAtomN->xyz.X,pAtomN->xyz.Y,pAtomN->xyz.Z,
                                    pAtomCA->xyz.X,pAtomCA->xyz.Y,pAtomCA->xyz.Z);

    p12 = XYZDifference(&pAtomCA->xyz,&pAtomN->xyz);    
    p23 = XYZDifference(&pAtomCA->xyz,&pAtomC->xyz);    
    pResiFlex->len[2] = XYZNormalization(&p23);
    pResiFlex->ang[2] = RadToDeg(XYZAngle(&p12,&p23));
    pResiFlex->dih[2] = calcTorsion(pAtomCprev->xyz.X,pAtomCprev->xyz.Y,pAtomCprev->xyz.Z,
                                    pAtomN->xyz.X,pAtomN->xyz.Y,pAtomN->xyz.Z,
                                    pAtomCA->xyz.X,pAtomCA->xyz.Y,pAtomCA->xyz.Z,
                                    pAtomC->xyz.X,pAtomC->xyz.Y,pAtomC->xyz.Z);
  }
  return Success;
}








/* convert coordinate from torsion angle space to cartesian space
 *  type - atom types for coordinate update
 *  1 and 6: CA only
 *  3: N, CA, C */
/*
BOOL tor2str(Chain *pChain,ChainFlex *pChainFlex){
  int i;
  XYZ pt,pn,pc;
  BOOL flagts;
  BOOL flagwhole=true;

  ResidueFlex *pResidueFlex0 = ChainFlexGetResidue(pChainFlex,0);
  Residue *pResidue0 = ChainGetResidue(pChain,0);
  Atom *pAtomN0 = ResidueGetAtomByName(pResidue0, "N");
  Atom *pAtomCA0 = ResidueGetAtomByName(pResidue0, "CA");
  Atom *pAtomC0 = ResidueGetAtomByName(pResidue0, "C");

  if(pResidueFlex0->len[1]<0) pResidueFlex0->len[1]=float(lennca);
  if(pResidueFlex0->len[2]<0) pResidueFlex0->len[2]=float(lencac);
  if(pResidueFlex0->ang[2]<0) pResidueFlex0->ang[2]=float(angncac);
  pAtomN0->xyz.X=0;
  pAtomN0->xyz.Y=0;
  pAtomN0->xyz.Z=0;
  pAtomCA0->xyz.X=pResidueFlex0->len[1];
  pAtomCA0->xyz.Y=0;
  pAtomCA0->xyz.Z=0;  
  pAtomC0->xyz.X=pResidueFlex0->len[1]-pResidueFlex0->len[2]*cos(pResidueFlex0->ang[2]*raddeg);
  pAtomC0->xyz.Y=pResidueFlex0->len[2]*sin(pResidueFlex0->ang[2]*raddeg);
  pAtomC0->xyz.Z=0;

  if(pResidueFlex0->dih[0]>=0 && pResidueFlex0->dih[2]>=0){
    flagts=convertTorToPos(0,0,0,
                           lennca,0,0,
                           lennca+lencac*sin(angncac*raddeg),-lencac*cos(angncac*raddeg),0,
                           pResidueFlex0->dih[0]*raddeg,float(lencn),
                           float(angcacn)*raddeg,&pn.X,&pn.Y,&pn.Z);
    if(!flagts){
      flagwhole=FALSE;
      //////printf("wrong front coordinates n %d\n",0);
    }
   
    pAtomN0->xyz.X=pn.X;
    pAtomN0->xyz.Y=pn.Y;
    pAtomN0->xyz.Z=pn.Z;
   
    if(pResidueFlex0->dih[1]<0) pResidueFlex0->dih[1]=180;
    flagts=convertTorToPos(lennca,0,0,
                           lennca+lencac*sin(angncac*raddeg),-lencac*cos(angncac*raddeg),0,
                           pAtomN0->xyz.X,pAtomN0->xyz.Y,pAtomN0->xyz.Z,
                           pResidueFlex0->dih[1]*raddeg,float(lennca),
                           angcnca*raddeg,&pt.X,&pt.Y,&pt.Z);
    if(!flagts){
      flagwhole=FALSE;
      //printf("wrong front coordinates ca %d\n",0);
    }
   
    pAtomCA0->xyz.X=pt.X;
    pAtomCA0->xyz.Y=pt.Y;
    pAtomCA0->xyz.Z=pt.Z;
    flagts=convertTorToPos(lennca+lencac*sin(angncac*raddeg),-lencac*cos(angncac*raddeg),0,
                           pAtomN0->xyz.X,pAtomN0->xyz.Y,pAtomN0->xyz.Z, 
                           pAtomCA0->xyz.X,pAtomCA0->xyz.Y,pAtomCA0->xyz.Z,
                           pResidueFlex0->dih[2]*raddeg,float(lencac),
                           angncac*raddeg,&pc.X,&pc.Y,&pc.Z);
    if(!flagts){
      flagwhole=FALSE;
      //printf("wrong front coordinates c %d\n",0);
    }

    pAtomC0->xyz.X=pc.X;
    pAtomC0->xyz.Y=pc.Y;
    pAtomC0->xyz.Z=pc.Z;
  }
  for(i=1;i<ChainGetResidueCount(pChain);i++){
    ResidueFlex *pResidueFlex = ChainFlexGetResidue(pChainFlex,i);
    Residue *pResiduePrev = ChainGetResidue(pChain,i-1);
    Residue *pResidue = ChainGetResidue(pChain,i);

    Atom *pAtomNprev = ResidueGetAtomByName(pResiduePrev, "N");
    Atom *pAtomCAprev = ResidueGetAtomByName(pResiduePrev, "CA");
    Atom *pAtomCprev = ResidueGetAtomByName(pResiduePrev, "C");
    Atom *pAtomN = ResidueGetAtomByName(pResidue, "N");
    Atom *pAtomCA = ResidueGetAtomByName(pResidue, "CA");
    Atom *pAtomC = ResidueGetAtomByName(pResidue, "C");


    if(pResidueFlex->dih[0]<0){
      if(pResidueFlex->ss=='E') pResidueFlex->dih[0]=120;
      else if(pResidueFlex->ss=='H') pResidueFlex->dih[0]=300;
      else pResidueFlex->dih[0]=120;
    }
    if(pResidueFlex->dih[1]<0) pResidueFlex->dih[1]=180;
    if(pResidueFlex->tor[2]<0) pResidueFlex->tor[2]=290;
    if(pResidueFlex->len[0]<0) pResidueFlex->len[0]=float(lencn);
    if(pResidueFlex->len[1]<0) pResidueFlex->len[1]=float(lennca);
    if(pResidueFlex->len[2]<0) pResidueFlex->len[2]=float(lencac);
    if(pResidueFlex->ang[0]<0) pResidueFlex->ang[0]=float(angcacn);
    if(pResidueFlex->ang[1]<0) pResidueFlex->ang[1]=float(angcnca);
    if(pResidueFlex->ang[2]<0) pResidueFlex->ang[2]=float(angncac);
    //0 1 2 n ca c
    //original phi i-1 psi i-1 omega i      
    flagts=convertTorToPos(pAtomNprev->xyz.X,pAtomNprev->xyz.Y,pAtomNprev->xyz.Z,
                           pAtomCAprev->xyz.X,pAtomCAprev->xyz.Y,pAtomCAprev->xyz.Z,
                           pAtomCprev->xyz.X,pAtomCprev->xyz.Y,pAtomCprev->xyz.Z,
                           pResidueFlex->dih[0]*raddeg,pResidueFlex->len[0],
                           pResidueFlex->ang[0]*raddeg,&pn.X,&pn.Y,&pn.Z);
    if(!flagts){
      flagwhole=FALSE;
      //printf("wrong front coordinates n %d %f %f %f\n",
      //i,decstr[i].tor[0],decstr[i].len[0],decstr[i].ang[0]);
    }

    pAtomN->xyz.X=pn.X;
    pAtomN->xyz.Y=pn.Y;
    pAtomN->xyz.Z=pn.Z;

    flagts=convertTorToPos(pAtomCAprev->xyz.X,pAtomCAprev->xyz.Y,pAtomCAprev->xyz.Z, 
                           pAtomCprev->xyz.X,pAtomCprev->xyz.Y,pAtomCprev->xyz.Z,
                           pAtomN->xyz.X,pAtomN->xyz.Y,pAtomN->xyz.Z,
                           pResidueFlex->dih[1]*raddeg,pResidueFlex->len[1],
                           pResidueFlex->ang[1]*raddeg,&pt.X,&pt.Y,&pt.Z);
    if(!flagts){
      flagwhole=FALSE;
      //printf("wrong front coordinates ca %d %f %f %f\n",
      //i,decstr[i].tor[1],decstr[i].len[1],decstr[i].ang[1]);
    }
   
    pAtomCA->xyz.X=pt.X;
    pAtomCA->xyz.Y=pt.Y;
    pAtomCA->xyz.Z=pt.Z;

    flagts=convertTorToPos(pAtomCprev->xyz.X,pAtomCprev->xyz.Y,pAtomCprev->xyz.Z,
                           pAtomN->xyz.X,pAtomN->xyz.Y,pAtomN->xyz.Z,
                           pAtomCA->xyz.X,pAtomCA->xyz.Y,pAtomCA->xyz.Z,
                           pResidueFlex->dih[2]*raddeg,pResidueFlex->len[2],
                           pResidueFlex->ang[2]*raddeg,&pc.X,&pc.Y,&pc.Z);
    if(!flagts){
      flagwhole=FALSE;
      //printf("wrong front coordinates c %d %f %f %f\n",
      //i,decstr[i].tor[2],decstr[i].len[2],decstr[i].ang[2]);
    }

    pAtomC->xyz.X=pc.X;
    pAtomC->xyz.Y=pc.Y;
    pAtomC->xyz.Z=pc.Z;
  }
  return flagwhole;
}
*/

BOOL tor2strp2p2(Chain *pChain,ChainFlex *pChainFlex,int istart,int iend){
  int i;
  double raddeg = PI/180.0;
  XYZ pn;
  BOOL flagts;
  BOOL flagwhole=TRUE;
  int realstart=istart;
  static double cbsta[][3]={
    1.52369,1.92448,122.35124, 1.52962,1.91875,122.28332, 1.53149,1.92096,122.53073, 1.53132,1.92149,122.55859,
    1.53219,1.91936,122.36077, 1.51371,1.90607,121.58025, 1.53172,1.92135,122.58755, 1.54507,1.92240,122.99168,
    1.53180,1.92265,122.48313, 1.53108,1.92040,122.28572, 1.53078,1.91922,122.34940, 1.53148,1.92241,122.84907,
    1.52996,1.94084,115.54662, 1.53057,1.92128,122.53531, 1.53085,1.92046,122.42439, 1.52991,1.91734,122.39611,
    1.54070,1.91400,122.79225, 1.54500,1.92132,123.02119, 1.53172,1.92002,122.56818, 1.53251,1.91842,122.36359
  };//n c ca cb  len ang tor

  if(realstart==0) realstart=1;
  for(i=realstart;i<=iend;i++){
    Residue *pResiPrev = ChainGetResidue(pChain,i-1);
    Residue *pResi = ChainGetResidue(pChain,i);
    Atom* pAtomNprev = ResidueGetAtomByName(pResiPrev, "N");
    Atom* pAtomCAprev = ResidueGetAtomByName(pResiPrev, "CA");
    Atom* pAtomCprev = ResidueGetAtomByName(pResiPrev, "C");
    Atom* pAtomOprev = ResidueGetAtomByName(pResiPrev, "O");
    Atom* pAtomHprev = ResidueGetAtomByName(pResiPrev, "H");
    Atom* pAtomN = ResidueGetAtomByName(pResi, "N");
    Atom* pAtomCA = ResidueGetAtomByName(pResi, "CA");
    Atom* pAtomC = ResidueGetAtomByName(pResi, "C");
    Atom* pAtomO = ResidueGetAtomByName(pResi, "O");
    Atom* pAtomH = ResidueGetAtomByName(pResi, "H");

    /////////////o
    flagts=convertTorToPos(pAtomN->xyz.X,pAtomN->xyz.Y,pAtomN->xyz.Z,
                           pAtomCAprev->xyz.X,pAtomCAprev->xyz.Y,pAtomCAprev->xyz.Z,
                           pAtomCprev->xyz.X,pAtomCprev->xyz.Y,pAtomCprev->xyz.Z,
                           179.6715f*raddeg,1.229f,2.0961f,&pn.X,&pn.Y,&pn.Z);//molprobity
    if(!flagts){       
      flagwhole=FALSE;
      printf("wrong front coordinates d %d\n",i);
    }
    pAtomOprev->xyz.X=pn.X;pAtomOprev->xyz.Y=pn.Y;pAtomOprev->xyz.Z=pn.Z;

    if(pAtomH!=NULL){
      flagts=convertTorToPos(pAtomCprev->xyz.X,pAtomCprev->xyz.Y,pAtomCprev->xyz.Z,
                             pAtomCA->xyz.X,pAtomCA->xyz.Y,pAtomCA->xyz.Z,
                             pAtomN->xyz.X,pAtomN->xyz.Y,pAtomN->xyz.Z,
                             179.8174f*raddeg,0.987f,2.0814f,&pn.X,&pn.Y,&pn.Z);//0.9919f,2.0574f
      if(!flagts){
        flagwhole=FALSE;
        printf("wrong front coordinates d %d\n",i);
      }
      pAtomH->xyz.X=pn.X;pAtomH->xyz.Y=pn.Y;pAtomH->xyz.Z=pn.Z;
    }
  }
  if(iend==ChainGetResidueCount(pChain)-1){
    i=ChainGetResidueCount(pChain);
    Residue *pResiPrev = ChainGetResidue(pChain,i-1);
    Atom* pAtomNprev = ResidueGetAtomByName(pResiPrev, "N");
    Atom* pAtomCAprev = ResidueGetAtomByName(pResiPrev, "CA");
    Atom* pAtomCprev = ResidueGetAtomByName(pResiPrev, "C");
    Atom* pAtomOprev = ResidueGetAtomByName(pResiPrev, "O");
    Atom* pAtomHprev = ResidueGetAtomByName(pResiPrev, "H");

    flagts=convertTorToPos(pAtomNprev->xyz.X,pAtomNprev->xyz.Y,pAtomNprev->xyz.Z,
                           pAtomCAprev->xyz.X,pAtomCAprev->xyz.Y,pAtomCAprev->xyz.Z,
                           pAtomCprev->xyz.X,pAtomCprev->xyz.Y,pAtomCprev->xyz.Z,
                           0.0f,1.2439f, 2.0855f,&pn.X,&pn.Y,&pn.Z);
    if(!flagts){
      flagwhole=FALSE;
      printf("wrong front coordinates o %d\n",i);
    }
    pAtomOprev->xyz.X=pn.X;pAtomOprev->xyz.Y=pn.Y;pAtomOprev->xyz.Z=pn.Z;
  }

  if(istart==0){
    i=0;
    Residue *pResi = ChainGetResidue(pChain,i);
    Atom* pAtomN = ResidueGetAtomByName(pResi, "N");
    Atom* pAtomCA = ResidueGetAtomByName(pResi, "CA");
    Atom* pAtomC = ResidueGetAtomByName(pResi, "C");
    Atom* pAtomO = ResidueGetAtomByName(pResi, "O");
    Atom* pAtomH = ResidueGetAtomByName(pResi, "H");
    
    if(pAtomH!=NULL){
      flagts=convertTorToPos(pAtomC->xyz.X,pAtomC->xyz.Y,pAtomC->xyz.Z,
                             pAtomCA->xyz.X,pAtomCA->xyz.Y,pAtomCA->xyz.Z,
                             pAtomN->xyz.X,pAtomN->xyz.Y,pAtomN->xyz.Z,
                             60*raddeg,0.987f,2.0306f,&pn.X,&pn.Y,&pn.Z);//0.9972f
      if(!flagts){
        flagwhole=FALSE;
        printf("wrong front coordinates %d\n",i);
      }
      pAtomH->xyz.X=pn.X;pAtomH->xyz.Y=pn.Y;pAtomH->xyz.Z=pn.Z;
    }
  }

  for(i=istart;i<=iend;i++){
    //cb
    Residue *pResi = ChainGetResidue(pChain,i);   
    Atom* pAtomN = ResidueGetAtomByName(pResi, "N");
    Atom* pAtomCA = ResidueGetAtomByName(pResi, "CA");
    Atom* pAtomC = ResidueGetAtomByName(pResi, "C");
    Atom* pAtomO = ResidueGetAtomByName(pResi, "O");
    Atom* pAtomH = ResidueGetAtomByName(pResi, "H");
    Atom* pAtomCB = ResidueGetAtomByName(pResi, "CB");

    if(strcmp(ResidueGetName(pResi), "GLY") != 0){
      flagts=convertTorToPos(pAtomN->xyz.X,pAtomN->xyz.Y,pAtomN->xyz.Z,
                             pAtomC->xyz.X,pAtomC->xyz.Y,pAtomC->xyz.Z,
                             pAtomCA->xyz.X,pAtomCA->xyz.Y,pAtomCA->xyz.Z,
                             cbsta[getResidueIndex(pResi)][2]*raddeg,cbsta[getResidueIndex(pResi)][0],cbsta[getResidueIndex(pResi)][1],
                             &pn.X,&pn.Y,&pn.Z);
      if(!flagts){
        flagwhole=FALSE;
        printf("wrong front coordinates cb2 %d\n",i);
      }
      pAtomCB->xyz.X=pn.X;pAtomCB->xyz.Y=pn.Y;pAtomCB->xyz.Z=pn.Z;
    }
  }
  return flagwhole;
}

int str2torp(Chain *pChain,ChainFlex *pChainFlex,int istart,int iend){
  int i;
  XYZ p12,p23;
  int realstart;
  double lencn = 1.338,lennca = 1.460,lencac = 1.525,angncac = 111.008;
  double raddeg = PI/180.0,angcacn = 116.617,angcnca = 121.614,lennc = 2.460;
  
  if(istart==0) realstart=1;
  else realstart=istart;
        
  if(realstart==1){
    Residue *pResi0 = ChainGetResidue(pChain,0);
    ResidueFlex *pResiFlex0 = ChainFlexGetResidue(pChainFlex,0);
    Atom* pAtomN0 = ResidueGetAtomByName(pResi0, "N");
    Atom* pAtomCA0 = ResidueGetAtomByName(pResi0, "CA");
    Atom* pAtomC0 = ResidueGetAtomByName(pResi0, "C");

    if(pResiFlex0->dih[0]<0 || pResiFlex0->len[0]<0){
      pResiFlex0->dih[0]=180.0;
      pResiFlex0->dih[1]=180.0;
      pResiFlex0->dih[2]=180.0;
      pResiFlex0->len[0]=lennc;
      pResiFlex0->ang[0]=angcacn;
      pResiFlex0->ang[1]=angcnca;
    }

    p12 = XYZDifference(&pAtomCA0->xyz,&pAtomN0->xyz);
    p23 = XYZDifference(&pAtomCA0->xyz,&pAtomC0->xyz);
    pResiFlex0->len[1]=XYZNormalization(&p12);
    pResiFlex0->len[2]=XYZNormalization(&p23);
    pResiFlex0->ang[2] = RadToDeg(XYZAngle(&p12, &p23));
  }
  for(i=realstart;i<=iend;i++){
    Residue *pResiPrev = ChainGetResidue(pChain,i-1);
    Residue *pResi = ChainGetResidue(pChain,i);
    ResidueFlex *pResiFlex = ChainFlexGetResidue(pChainFlex,i);
    Atom* pAtomNprev = ResidueGetAtomByName(pResiPrev, "N");
    Atom* pAtomCAprev = ResidueGetAtomByName(pResiPrev, "CA");
    Atom* pAtomCprev = ResidueGetAtomByName(pResiPrev, "C");
    Atom* pAtomN = ResidueGetAtomByName(pResi, "N");
    Atom* pAtomCA = ResidueGetAtomByName(pResi, "CA");
    Atom* pAtomC = ResidueGetAtomByName(pResi, "C");

    p12 = XYZDifference(&pAtomCprev->xyz,&pAtomCAprev->xyz);
    p23 = XYZDifference(&pAtomCprev->xyz,&pAtomN->xyz);
    pResiFlex->len[0] = XYZNormalization(&p23);
    pResiFlex->ang[0] = RadToDeg(XYZAngle(&p12, &p23));
    pResiFlex->dih[0] = calcTorsion(pAtomNprev->xyz.X,pAtomNprev->xyz.Y,pAtomNprev->xyz.Z,
                                    pAtomCAprev->xyz.X,pAtomCAprev->xyz.Y,pAtomCAprev->xyz.Z,
                                    pAtomCprev->xyz.X,pAtomCprev->xyz.Y,pAtomCprev->xyz.Z,
                                    pAtomN->xyz.X,pAtomN->xyz.Y,pAtomN->xyz.Z);

    p12 = XYZDifference(&pAtomN->xyz,&pAtomCprev->xyz);
    p23 = XYZDifference(&pAtomN->xyz,&pAtomCA->xyz);
    pResiFlex->len[1] = XYZNormalization(&p23);
    pResiFlex->ang[1] = RadToDeg(XYZAngle(&p12,&p23));
    pResiFlex->dih[1] = calcTorsion(pAtomCAprev->xyz.X,pAtomCAprev->xyz.Y,pAtomCAprev->xyz.Z,
                                    pAtomCprev->xyz.X,pAtomCprev->xyz.Y,pAtomCprev->xyz.Z,
                                    pAtomN->xyz.X,pAtomN->xyz.Y,pAtomN->xyz.Z,
                                    pAtomCA->xyz.X,pAtomCA->xyz.Y,pAtomCA->xyz.Z);

    p12 = XYZDifference(&pAtomCA->xyz,&pAtomN->xyz);
    p23 = XYZDifference(&pAtomCA->xyz,&pAtomC->xyz);
    pResiFlex->len[2] = XYZNormalization(&p23);
    pResiFlex->ang[2] = RadToDeg(XYZAngle(&p12,&p23));
    pResiFlex->dih[2] = calcTorsion(pAtomCprev->xyz.X,pAtomCprev->xyz.Y,pAtomCprev->xyz.Z,
                                    pAtomN->xyz.X,pAtomN->xyz.Y,pAtomN->xyz.Z,
                                    pAtomCA->xyz.X,pAtomCA->xyz.Y,pAtomCA->xyz.Z,
                                    pAtomC->xyz.X,pAtomC->xyz.Y,pAtomC->xyz.Z);
  }
  return Success;
}

BOOL NewMainChainPos(Chain* pChain, ChainFlex* pChainFlex, int istart){
  int i;
  XYZ pt,pn,pc;
  BOOL flagts;
  BOOL flagWhole=TRUE;
  double lencn = 1.338,lennca = 1.460,lencac = 1.525,angncac = 111.008;
  double raddeg = PI/180.0,angcacn = 116.617,angcnca = 121.614;

  static double cbsta[][3]={
    1.52369,1.92448,122.35124, 1.52962,1.91875,122.28332, 1.53149,1.92096,122.53073, 1.53132,1.92149,122.55859,
    1.53219,1.91936,122.36077, 1.51371,1.90607,121.58025, 1.53172,1.92135,122.58755, 1.54507,1.92240,122.99168,
    1.53180,1.92265,122.48313, 1.53108,1.92040,122.28572, 1.53078,1.91922,122.34940, 1.53148,1.92241,122.84907,
    1.52996,1.94084,115.54662, 1.53057,1.92128,122.53531, 1.53085,1.92046,122.42439, 1.52991,1.91734,122.39611,
    1.54070,1.91400,122.79225, 1.54500,1.92132,123.02119, 1.53172,1.92002,122.56818, 1.53251,1.91842,122.36359
  };//n c ca cb  len ang tor

  if(istart == 0){
    Residue* newResi = ChainGetResidue(pChain,istart);
    ResidueFlex* newResiFlex = ChainFlexGetResidue(pChainFlex,istart);
    Atom* pAtomN = ResidueGetAtomByName(newResi, "N");
    Atom* pAtomCA = ResidueGetAtomByName(newResi, "CA");
    Atom* pAtomC = ResidueGetAtomByName(newResi, "C");

    newResiFlex->len[1] = double(lennca);
    newResiFlex->len[2] = double(lencac);
    newResiFlex->ang[2] = double(angncac);
    pAtomN->xyz.X = 0;
    pAtomN->xyz.Y = 0;
    pAtomN->xyz.Z = 0;
    pAtomCA->xyz.X = double(lennca);
    pAtomCA->xyz.Y = 0;
    pAtomCA->xyz.Z = 0;
    pAtomC->xyz.X = double(lennca) - double(lencac) * cos(double(angncac) * raddeg);
    pAtomC->xyz.Y = double(lencac) * sin(double(angncac) * raddeg);
    pAtomC->xyz.Z = 0;
    if(newResiFlex->dih[0] >= 0 && newResiFlex->dih[2] >= 0){
      flagts=convertTorToPos(0,0,0,
                             lennca,0,0,
                             lennca + lencac * sin(angncac * raddeg),-lencac * cos(angncac * raddeg),0,
                             newResiFlex->dih[0] * raddeg,double(lencn),double(angcacn)*raddeg,
                             &pn.X,&pn.Y,&pn.Z);
      if(!flagts){
        flagWhole=FALSE;
        printf("wrong front coordinates n %d\n",0);
      }
      pAtomN->xyz.X = pn.X; pAtomN->xyz.Y = pn.Y; pAtomN->xyz.Z = pn.Z;
      if(newResiFlex->dih[1]<0) newResiFlex->dih[1]=180;
      flagts=convertTorToPos(lennca,0,0,
                             lennca+lencac*sin(angncac*raddeg),-lencac*cos(angncac*raddeg),0,
                             pAtomN->xyz.X,pAtomN->xyz.Y,pAtomN->xyz.Z,
                             newResiFlex->dih[1]*raddeg,double(lennca),angcnca*raddeg,
                             &pt.X,&pt.Y,&pt.Z);
      if(!flagts){
        flagWhole=FALSE;
        printf("wrong front coordinates ca %d\n",0);
      }
      pAtomCA->xyz.X=pt.X; pAtomCA->xyz.Y=pt.Y; pAtomCA->xyz.Z=pt.Z;
      flagts=convertTorToPos(lennca+lencac*sin(angncac*raddeg),-lencac*cos(angncac*raddeg),0,
                             pAtomN->xyz.X,pAtomN->xyz.Y,pAtomN->xyz.Z,
                             pAtomCA->xyz.X,pAtomCA->xyz.Y,pAtomCA->xyz.Z,
                             newResiFlex->dih[2]*raddeg,double(lencac),angncac*raddeg,
                             &pc.X,&pc.Y,&pc.Z);
       if(!flagts){
         flagWhole = FALSE;
         printf("wrong front coordinates c %d\n",0);
       }
       pAtomC->xyz.X = pc.X; pAtomC->xyz.Y = pc.Y; pAtomC->xyz.Z = pc.Z;
     }
  }
  int realstart = istart;
  if(realstart == 0) realstart = 1;
  for(i = realstart; i < ChainGetResidueCount(pChain); i++){
    Residue* newResiPrev = ChainGetResidue(pChain,i-1);
    Residue* newResi = ChainGetResidue(pChain,i);
    ResidueFlex* newResiFlex = ChainFlexGetResidue(pChainFlex,i);
    Atom* pAtomNPrev = ResidueGetAtomByName(newResiPrev, "N");
    Atom* pAtomCAPrev = ResidueGetAtomByName(newResiPrev, "CA");
    Atom* pAtomCPrev = ResidueGetAtomByName(newResiPrev, "C");
    Atom* pAtomOPrev = ResidueGetAtomByName(newResiPrev, "O");
    Atom* pAtomHPrev = ResidueGetAtomByName(newResiPrev, "H");
    Atom* pAtomN = ResidueGetAtomByName(newResi, "N");
    Atom* pAtomCA = ResidueGetAtomByName(newResi, "CA");
    Atom* pAtomC = ResidueGetAtomByName(newResi, "C");
    Atom* pAtomO = ResidueGetAtomByName(newResi, "O");
    Atom* pAtomH = ResidueGetAtomByName(newResi, "H");

    if(newResiFlex->dih[0]<0){
      if(newResiFlex->ss == 'E'){
        newResiFlex->dih[0] = 120;
      }
      else if(newResiFlex->ss == 'H'){
        newResiFlex->dih[0] = 300;
      }
      else{
        newResiFlex->dih[0] = 120;
      }
    }
    if(newResiFlex->dih[1] < 0){
      newResiFlex->dih[1] = 180;
    }
    if(newResiFlex->dih[2] < 0){
      newResiFlex->dih[2] = 290;
    }
    if(newResiFlex->len[0] < 0){
      newResiFlex->len[0] = double(lencn);
    }
    if(newResiFlex->len[1] < 0){
      newResiFlex->len[1] = double(lennca);
    }
    if(newResiFlex->len[2] < 0){
      newResiFlex->len[2] = double(lencac);
    }
    if(newResiFlex->ang[0] < 0){
      newResiFlex->ang[0] = double(angcacn);
    }
    if(newResiFlex->ang[1] < 0){
      newResiFlex->ang[1] = double(angcnca);
    }
    if(newResiFlex->ang[2] < 0){
      newResiFlex->ang[2] = double(angncac);
    }

    flagts=convertTorToPos(pAtomNPrev->xyz.X,pAtomNPrev->xyz.Y,pAtomNPrev->xyz.Z,
                           pAtomCAPrev->xyz.X,pAtomCAPrev->xyz.Y,pAtomCAPrev->xyz.Z,
                           pAtomCPrev->xyz.X,pAtomCPrev->xyz.Y,pAtomCPrev->xyz.Z,
                           newResiFlex->dih[0]*raddeg,newResiFlex->len[0],newResiFlex->ang[0]*raddeg,
                           &pn.X,&pn.Y,&pn.Z);
    if(!flagts){
      flagWhole = FALSE;
      //printf("wrong front coordinatesp n %d %f %f %f\n",i,decstr[i].dih[0],decstr[i].len[0],decstr[i].ang[0]);
    }
    pAtomN->xyz.X=pn.X;pAtomN->xyz.Y=pn.Y;pAtomN->xyz.Z=pn.Z;

    flagts=convertTorToPos(pAtomCAPrev->xyz.X,pAtomCAPrev->xyz.Y,pAtomCAPrev->xyz.Z,
                           pAtomCPrev->xyz.X,pAtomCPrev->xyz.Y,pAtomCPrev->xyz.Z,
                           pAtomN->xyz.X,pAtomN->xyz.Y,pAtomN->xyz.Z,
                           newResiFlex->dih[1]*raddeg,newResiFlex->len[1],newResiFlex->ang[1]*raddeg,
                           &pt.X,&pt.Y,&pt.Z);
    if(!flagts){
      flagWhole = FALSE;
       //printf("wrong front coordinatesp ca %d %f %f %f\n",i,decstr[i].dih[1],decstr[i].len[1],decstr[i].ang[1]);
    }
    pAtomCA->xyz.X=pt.X;pAtomCA->xyz.Y=pt.Y;pAtomCA->xyz.Z=pt.Z;

    flagts=convertTorToPos(pAtomCPrev->xyz.X,pAtomCPrev->xyz.Y,pAtomCPrev->xyz.Z,
                           pAtomN->xyz.X,pAtomN->xyz.Y,pAtomN->xyz.Z,
                           pAtomCA->xyz.X,pAtomCA->xyz.Y,pAtomCA->xyz.Z,
                           newResiFlex->dih[2]*raddeg,newResiFlex->len[2],newResiFlex->ang[2]*raddeg,
                           &pc.X,&pc.Y,&pc.Z);
    if(!flagts){
      flagWhole = FALSE;
      //printf("wrong front coordinatesp c %d %f %f %f\n",i,decstr[i].dih[2],decstr[i].len[2],decstr[i].ang[2]);
    }
    pAtomC->xyz.X=pc.X;pAtomC->xyz.Y=pc.Y;pAtomC->xyz.Z=pc.Z;
  }
  for(i = realstart; i < ChainGetResidueCount(pChain); i++){
    Residue* newResiPrev = ChainGetResidue(pChain,i-1);
    Residue* newResi = ChainGetResidue(pChain,i);
    Atom* pAtomNPrev = ResidueGetAtomByName(newResiPrev, "N");
    Atom* pAtomCAPrev = ResidueGetAtomByName(newResiPrev, "CA");
    Atom* pAtomCPrev = ResidueGetAtomByName(newResiPrev, "C");
    Atom* pAtomOPrev = ResidueGetAtomByName(newResiPrev, "O");
    Atom* pAtomHPrev = ResidueGetAtomByName(newResiPrev, "H");
    Atom* pAtomN = ResidueGetAtomByName(newResi, "N");
    Atom* pAtomCA = ResidueGetAtomByName(newResi, "CA");
    Atom* pAtomC = ResidueGetAtomByName(newResi, "C");
    Atom* pAtomO = ResidueGetAtomByName(newResi, "O");
    Atom* pAtomH = ResidueGetAtomByName(newResi, "H");
    flagts=convertTorToPos(pAtomN->xyz.X,pAtomN->xyz.Y,pAtomN->xyz.Z,
                           pAtomCAPrev->xyz.X,pAtomCAPrev->xyz.Y,
                           pAtomCAPrev->xyz.Z,pAtomCPrev->xyz.X,pAtomCPrev->xyz.Y,pAtomCPrev->xyz.Z,
                           179.6715f*raddeg,1.229f,2.0961f,
                           &pn.X,&pn.Y,&pn.Z);//molprobity
    if(!flagts){
      flagWhole = FALSE;
      printf("wrong front coordinates d %d\n",i);
    }
    pAtomOPrev->xyz.X=pn.X;pAtomOPrev->xyz.Y=pn.Y;pAtomOPrev->xyz.Z=pn.Z;

    if(pAtomH!=NULL){
      flagts=convertTorToPos(pAtomCPrev->xyz.X,pAtomCPrev->xyz.Y,pAtomCPrev->xyz.Z,
                             pAtomCA->xyz.X,pAtomCA->xyz.Y,pAtomCA->xyz.Z,
                             pAtomN->xyz.X,pAtomN->xyz.Y,pAtomN->xyz.Z,
                             179.8174f*raddeg,0.987f,2.0814f,
                             &pn.X,&pn.Y,&pn.Z);//0.9919f,2.0574f
      if(!flagts){
        flagWhole = FALSE;
        printf("wrong front coordinates d %d\n",i);
      }
      pAtomH->xyz.X = pn.X; pAtomH->xyz.Y = pn.Y; pAtomH->xyz.Z = pn.Z;
    }
  }
  i = ChainGetResidueCount(pChain)-1;
  Residue* newResi = ChainGetResidue(pChain,i);
  Atom* pAtomN = ResidueGetAtomByName(newResi, "N");
  Atom* pAtomCA = ResidueGetAtomByName(newResi, "CA");
  Atom* pAtomC = ResidueGetAtomByName(newResi, "C");
  Atom* pAtomO = ResidueGetAtomByName(newResi, "O");
  Atom* pAtomH = ResidueGetAtomByName(newResi, "H");

  flagts=convertTorToPos(pAtomN->xyz.X,pAtomN->xyz.Y,pAtomN->xyz.Z,
                         pAtomCA->xyz.X,pAtomCA->xyz.Y,pAtomCA->xyz.Z,
                         pAtomC->xyz.X,pAtomC->xyz.Y,pAtomC->xyz.Z,
                         0.0f,1.2439f, 2.0855f,
                         &pn.X,&pn.Y,&pn.Z);
  if(!flagts){
    flagWhole=FALSE;
    printf("wrong front coordinates o %d\n",i);
  }
  pAtomO->xyz.X=pn.X;pAtomO->xyz.Y=pn.Y;pAtomO->xyz.Z=pn.Z;

  if(istart==0){
    i=0;
    Residue* newResi = ChainGetResidue(pChain,i);
    Atom* pAtomN = ResidueGetAtomByName(newResi, "N");
    Atom* pAtomCA = ResidueGetAtomByName(newResi, "CA");
    Atom* pAtomC = ResidueGetAtomByName(newResi, "C");
    Atom* pAtomO = ResidueGetAtomByName(newResi, "O");
    Atom* pAtomH = ResidueGetAtomByName(newResi, "H");

    if(pAtomH!=NULL){
      flagts=convertTorToPos(pAtomC->xyz.X,pAtomC->xyz.Y,pAtomC->xyz.Z,
                             pAtomCA->xyz.X,pAtomCA->xyz.Y,pAtomCA->xyz.Z,
                             pAtomN->xyz.X,pAtomN->xyz.Y,pAtomN->xyz.Z,
                             60*raddeg,0.987f,2.0306f,
                             &pn.X,&pn.Y,&pn.Z);//0.9972f
      if(!flagts){
        flagWhole=FALSE;
        //printf("wrong front coordinates %d\n",i);
      }
      pAtomH->xyz.X = pn.X; pAtomH->xyz.Y = pn.Y; pAtomH->xyz.Z = pn.Z;
    }
  }
  for(i = istart; i < ChainGetResidueCount(pChain); i++){
    Residue* newResi = ChainGetResidue(pChain,i);
    Atom* pAtomN = ResidueGetAtomByName(newResi, "N");
    Atom* pAtomCA = ResidueGetAtomByName(newResi, "CA");
    Atom* pAtomC = ResidueGetAtomByName(newResi, "C");
    Atom* pAtomO = ResidueGetAtomByName(newResi, "O");
    Atom* pAtomH = ResidueGetAtomByName(newResi, "H");

    if(strcmp(ResidueGetName(newResi),"GLY")!=0){
      int residueIndex=getResidueIndex(newResi);
      Atom* pAtomCB = ResidueGetAtomByName(newResi, "CB");

      if(pAtomCB!=NULL){
        flagts=convertTorToPos(pAtomN->xyz.X,pAtomN->xyz.Y,pAtomN->xyz.Z,
                               pAtomC->xyz.X,pAtomC->xyz.Y,pAtomC->xyz.Z,
                               pAtomCA->xyz.X,pAtomCA->xyz.Y,pAtomCA->xyz.Z,
                               cbsta[residueIndex][2]*raddeg,cbsta[residueIndex][0],cbsta[residueIndex][1],
                               &pn.X,&pn.Y,&pn.Z);
        if(!flagts){
          flagWhole=FALSE;
          //printf("wrong front coordinates cb2 %d\n",i);
        }
        pAtomCB->xyz.X = pn.X; pAtomCB->xyz.Y = pn.Y; pAtomCB->xyz.Z = pn.Z;
      }
    }
  }
  return flagWhole;
}

BOOL convertTorToPos(double xi,double yi,double zi,
                     double xj,double yj,double zj,
                     double xk,double yk,double zk,
                     double tang,double tleng,double tinner,
                     double *xl,double *yl,double *zl){
  XYZ p12,p23,e1,e2,e3;
  XYZ p1,p2,p3;
  double tpangle,tcos,tsin,rmat[9];
  double q[4];
  p12.X=xi-xj;p12.Y=yi-yj;p12.Z=zi-zj;
  p23.X=xk-xj;p23.Y=yk-yj;p23.Z=zk-zj;
  if(XYZNormalization(&p12)<err*0.00001 || XYZNormalization(&p23)<err*0.00001 || XYZAngle(&p12,&p23)<err*0.001 || (PI-XYZAngle(&p12,&p23))<err*0.001){
    int imax=XYZMaxNormal(&p23);
    if(imax==0){
      yj+=0.01f;
    }
    else if(imax==1){
      zj+=0.01f;
    }
    else{
      xj+=0.01f;
    }
    p12.X=xi-xj;p12.Y=yi-yj;p12.Z=zi-zj;
    p23.X=xk-xj;p23.Y=yk-yj;p23.Z=zk-zj;
    printf("make adjustment tor2pos22\n");
  }
  e2=XYZUnit(&p23);
  e3=XYZUnit(&p12);
  e1=XYZCrossProduct(&e2,&e3);
  e1=XYZUnit(&e1);
  if(XYZNormalization(&e1)<err || XYZNormalization(&e2)<err){
    printf("wrong in tor2pos22 [%f %f %f] [%f %f %f] [%f %f %f]\n",xi,yi,zi,xj,yj,zj,xk,yk,zk);
    *xl=xk;*yl=yk;*zl=zk;
    return FALSE;
  }
  p1=XYZScale2(&e2,tleng);
  tpangle=(PI-tinner)/2.0;
  tcos=cos(tpangle);
  tsin=sin(tpangle);
  q[0]=tcos;q[1]=tsin*e1.X;q[2]=tsin*e1.Y;q[3]=tsin*e1.Z;
  XYZQToRot(q,rmat);
  p2=XYZMmat(rmat,&p1);

  tpangle=tang/2.0;
  tcos=cos(tpangle);
  tsin=sin(tpangle);
  q[0]=tcos;q[1]=tsin*e2.X;q[2]=tsin*e2.Y;q[3]=tsin*e2.Z;
  XYZQToRot(q,rmat);
  p3=XYZMmat(rmat,&p2);

  *xl=p3.X+xk;*yl=p3.Y+yk;*zl=p3.Z+zk;
  return TRUE;
}

void extractCaDistances(Chain *pChain,pairaa *natDistRestr){
  int i,j,index; 
  XYZ distVec;
  double distance,stdDistance;
  int numseq = ChainGetResidueCount(pChain);

  for(i=0;i<numseq-2;i++){
    for(j=i+2;j<numseq;j++){
      Residue *pResidue1 = ChainGetResidue(pChain,i);
      Residue *pResidue2 = ChainGetResidue(pChain,j);
      Atom *pAtomCA1 = ResidueGetAtomByName(pResidue1, "CA");
      Atom *pAtomCA2 = ResidueGetAtomByName(pResidue2, "CA");
      index=i*numseq-(i-1)*i/2+j-i;
      distVec = XYZDifference(&pAtomCA1->xyz,&pAtomCA2->xyz);
      distance = XYZNormalization(&distVec);
      if(distance<3.9) stdDistance = 3.9;
      else stdDistance = distance;
      natDistRestr[index].dist = distance;
      natDistRestr[index].dstd = stdDistance;
    }
  }
}


double calcConstraintEnergy(Chain *pChain,pairaa *natDistRestr){
  int i,j,index; 
  XYZ distVec;
  double distance,stdDistance;
  double energy = 0;
  int numConstraints = 0;
  int numseq = ChainGetResidueCount(pChain);

  for(i=0;i<numseq-2;i++){
    for(j=i+2;j<numseq;j++){
      Residue *pResidue1 = ChainGetResidue(pChain,i);
      Residue *pResidue2 = ChainGetResidue(pChain,j);
      Atom *pAtomCA1 = ResidueGetAtomByName(pResidue1, "CA");
      Atom *pAtomCA2 = ResidueGetAtomByName(pResidue2, "CA");
      index=i*numseq-(i-1)*i/2+j-i;
      distVec = XYZDifference(&pAtomCA1->xyz,&pAtomCA2->xyz);
      distance = XYZNormalization(&distVec);
      energy += (fabs(distance-natDistRestr[index].dist)/natDistRestr[index].dstd);
      numConstraints++;
    }
  }
  return energy/double(numConstraints);
}


double calcConstraintEnergyWeighted(Chain *pChain,pairaa *natDistRestr,double globalWeight){
  int i,j,index; 
  XYZ distVec;
  double distance,stdDistance;
  double energy = 0;
  int numConstraints = 0;
  int numseq = ChainGetResidueCount(pChain);

  for(i=0;i<numseq-2;i++){
    for(j=i+2;j<numseq;j++){
      Residue *pResidue1 = ChainGetResidue(pChain,i);
      Residue *pResidue2 = ChainGetResidue(pChain,j);
      Atom *pAtomCA1 = ResidueGetAtomByName(pResidue1, "CA");
      Atom *pAtomCA2 = ResidueGetAtomByName(pResidue2, "CA");
      index=i*numseq-(i-1)*i/2+j-i;
      distVec = XYZDifference(&pAtomCA1->xyz,&pAtomCA2->xyz);
      distance = XYZNormalization(&distVec);
      energy += (fabs(distance-natDistRestr[index].dist)/natDistRestr[index].dstd)*(natDistRestr[index].weight);
      numConstraints++;
    }
  }
  return globalWeight*energy/double(numConstraints);
}

int OptimizeSideChains(Structure* pStructure,BBdepRotamerLib* pBBdepRotLib,int numCycles, int deschn){
  for(int i=0; i<numCycles; i++){
    const int numDesignSites=StructureGetDesignSiteCount(pStructure);
    int designSiteArray[numDesignSites];
    for(int j=0; j<numDesignSites; j++){
      designSiteArray[j]=j;
    }

    
    randomizeArray(designSiteArray, numDesignSites); 
   
    for(int j=0; j<numDesignSites; j++){
      ProteinSiteOptimizeRotamerWithBBdepRotLib(pStructure,deschn,designSiteArray[j],pBBdepRotLib);
    }
  }
  return Success;
}

int StructureCalcEnergyFlex(Structure* pStructure,BBdepRotamerLib* pBBdepRotLib,double energyTerms[]){

/*
  for(int i = 0; i < StructureGetChainCount(pStructure); i++){
    Chain *pChainI = StructureGetChain(pStructure,i);
    for(int ir = 0; ir < ChainGetResidueCount(pChainI); ir++){
      Residue *pResIR = ChainGetResidue(pChainI,ir);
      if(pResIR->designSiteType!=Type_ResidueDesignType_Fixed){
        for(int atom1=0; atom1<ResidueGetAtomCount(pResIR); atom1++){
          Atom* pAtom1 = ResidueGetAtom(pResIR, atom1);
          if(pAtom1->isBBAtom==FALSE && pAtom1->isXyzValid) pAtom1->isXyzValid=FALSE;
        }
      }
    }
  }
*/

  for(int i=0;i<StructureGetChainCount(pStructure);i++){
    Chain *pChainI=StructureGetChain(pStructure,i);
    for(int ir=0;ir<ChainGetResidueCount(pChainI);ir++){
      Residue *pResIR=ChainGetResidue(pChainI,ir);
      //if(pResIR->designSiteType==Type_ResidueDesignType_Fixed){
        energyTerms[91]+=pResIR->aapp;
        energyTerms[92]+=pResIR->ramachandran;
        energyTerms[93]+=pResIR->dunbrack;
        EVOEF_AminoAcidReferenceEnergy(pResIR->name,energyTerms);
        EVOEF_EnergyResidueIntraEnergy(pResIR,energyTerms);
      //}
      for(int is=ir+1;is<ChainGetResidueCount(pChainI);is++){
        Residue *pResIS=ChainGetResidue(pChainI,is);
        if(ResidueGetPosInChain(pResIR)+1==ResidueGetPosInChain(pResIS)) EVOEF_EnergyResidueAndNextResidue(pResIR,pResIS,energyTerms);
        else EVOEF_EnergyResidueAndOtherResidueSameChain(pResIR,pResIS,energyTerms);
        //EVOEF_EnergyResidueAndResidueSameChain(pResIR,pResIS,energyTerms);
      }
      for(int k=i+1; k<StructureGetChainCount(pStructure); k++){
        Chain *pChainK = StructureGetChain(pStructure,k);
        for(int ks=0; ks<ChainGetResidueCount(pChainK); ks++){
          Residue *pResKS = ChainGetResidue(pChainK,ks);
          if(ChainGetType(pChainI)==Type_Chain_SmallMol){
            EVOEF_EnergyResidueAndLigandResidue(pResKS,pResIR,energyTerms);
//            if(FLAG_PROT_LIG || FLAG_ENZYME){
//              EVOEF_EnergyResidueAndLigandResidue(pResKS,pResIR,energyTermsBind);
//            }
          }
          else if(ChainGetType(pChainK)==Type_Chain_SmallMol){
            EVOEF_EnergyResidueAndLigandResidue(pResIR,pResKS,energyTerms);
//            if(FLAG_PROT_LIG || FLAG_ENZYME){
//              EVOEF_EnergyResidueAndLigandResidue(pResIR,pResKS,energyTermsBind);
//            }
          }
          else{
            EVOEF_EnergyResidueAndOtherResidueDifferentChain(pResIR,pResKS,energyTerms);
            //EVOEF_EnergyResidueAndOtherResidueDifferentChain(pResIR,pResKS,energyTermsBind);
            //note that this code can easily cause a bug!!! (2/2/2020)
            if((strstr(DES_CHAINS,ChainGetName(pChainI))!=NULL && strstr(DES_CHAINS,ChainGetName(pChainK))==NULL) ||
              (strstr(DES_CHAINS,ChainGetName(pChainI))==NULL && strstr(DES_CHAINS,ChainGetName(pChainK))!=NULL)){
//                if(FLAG_PPI){
//                  EVOEF_EnergyResidueAndOtherResidueDifferentChain(pResIR,pResKS,energyTermsBind);
//                }
            }
          }
        }
      }
    }
  }

  //restore the coordinates
/*  for(int i=0; i<StructureGetChainCount(pStructure); i++){
    Chain *pChainI=StructureGetChain(pStructure,i);
    for(int ir=0; ir<ChainGetResidueCount(pChainI); ir++){
      Residue *pResIR = ChainGetResidue(pChainI,ir);
      if(pResIR->designSiteType!=Type_ResidueDesignType_Fixed){
        for(int atom1=0; atom1<ResidueGetAtomCount(pResIR); atom1++){
          Atom* pAtom1=ResidueGetAtom(pResIR,atom1);
          if(pAtom1->isBBAtom==FALSE) pAtom1->isXyzValid=TRUE;
        }
      }
    }
  }
*/

  return Success;
}

int getAngleOutliers(Chain *pChain,ChainFlex *pChainFlex){
  const int NUM_ANGLES=3;
  int i;
  int istart;
  double stddelta=3.9;
  double stdangle[][4]={
    {111.0,2.7*stddelta,111.0-2.7*stddelta,111.0+2.7*stddelta},//ncac
    {117.2,2.2*stddelta,117.2-2.2*stddelta,117.2+2.2*stddelta},//cacn
    {121.7,2.5*stddelta,121.7-2.5*stddelta,121.7+2.5*stddelta},//cnca

    {113.1,2.5*stddelta,113.1-2.5*stddelta,113.1+2.5*stddelta},//ncac gly
    {116.2,2.0*stddelta,116.2-2.0*stddelta,116.2+2.0*stddelta},//cacn gly 
    {122.3,2.1*stddelta,122.3-2.1*stddelta,122.3+2.1*stddelta},//cnca gly  

    {112.1,2.6*stddelta,112.1-2.6*stddelta,112.1+2.6*stddelta},//ncac pro  
    {117.1,2.8*stddelta,117.1-2.8*stddelta,117.1+2.8*stddelta},//cacn pro 
    {119.3,1.5*stddelta,119.3-1.5*stddelta,119.3+1.5*stddelta},//cnca pro

    {120.1,2.1*stddelta,120.1-2.1*stddelta,120.1+2.1*stddelta},//caco  
    {120.6,1.8*stddelta,120.6-1.8*stddelta,120.6+1.8*stddelta},//caco gly  
    {120.2,2.4*stddelta,120.2-2.4*stddelta,120.2+2.4*stddelta},//caco pro  

    {122.7,1.6*stddelta,122.7-1.6*stddelta,122.7+1.6*stddelta},//ocn  
    {123.2,1.7*stddelta,123.2-1.7*stddelta,123.2+1.7*stddelta},//ocn gly  
    {121.1,1.9*stddelta,121.1-1.9*stddelta,121.1+1.9*stddelta},//ocn pro  
  };

  for(i=0;i<NUM_ANGLES;i++){
    pChainFlex->numAngleOutliers[i]=0;   
  }
  for(int i=1; i<ChainFlexGetResidueCount(pChainFlex); i++){
    ResidueFlex* pResidueFlex = ChainFlexGetResidue(pChainFlex,i);
    Residue* pResidue = ChainGetResidue(pChain,i); 

    for(int j=0;j<NUM_ANGLES;j++){
      pResidueFlex->isAngleOutlierAbove[j] = FALSE;
      pResidueFlex->isAngleOutlierBelow[j] = FALSE;
    } 
 
    if(strcmp(ResidueGetName(pResidue), "GLY") == 0) istart=3;
    else if(strcmp(ResidueGetName(pResidue), "PRO") == 0) istart=6;
    else istart=0;

    if(pResidueFlex->ang[0]<stdangle[istart+1][2]){
      pChainFlex->numAngleOutliers[0]++;
      pResidueFlex->isAngleOutlierBelow[0] = TRUE;
    }
    else if(pResidueFlex->ang[0]>stdangle[istart+1][3]){
      pChainFlex->numAngleOutliers[0]++;
      pResidueFlex->isAngleOutlierAbove[0] = TRUE;
    }

    if(pResidueFlex->ang[1]<stdangle[istart+2][2]){
      pChainFlex->numAngleOutliers[1]++;
      pResidueFlex->isAngleOutlierBelow[1] = TRUE;
    }
    else if(pResidueFlex->ang[1]>stdangle[istart+2][3]){
      pChainFlex->numAngleOutliers[1]++;
      pResidueFlex->isAngleOutlierAbove[1] = TRUE;
    }

    if(pResidueFlex->ang[2]<stdangle[istart+0][2]){
      pChainFlex->numAngleOutliers[2]++;
      pResidueFlex->isAngleOutlierBelow[2] = TRUE;
    }
    else if(pResidueFlex->ang[2]>stdangle[istart+0][3]){
      pChainFlex->numAngleOutliers[2]++;
      pResidueFlex->isAngleOutlierAbove[2] = TRUE;
    }
  }
  return Success;
}

int getBondLengthOutliers(Chain *pChain,ChainFlex *pChainFlex){
  const int NUM_ANGLES=3,NUM_BONDS=4;
  int i;
  int istart;
  double stddelta=3.9;
  double stdlength[][4]={
    {1.459,0.020*stddelta,1.459-0.020*stddelta,1.459+0.020*stddelta},//nca  
    {1.525,0.026*stddelta,1.525-0.026*stddelta,1.525+0.026*stddelta},//cac  
    {1.336,0.023*stddelta,1.336-0.023*stddelta,1.336+0.023*stddelta},//cn  
    {3.813,0.080*stddelta,3.813-0.080*stddelta,3.813+0.080*stddelta},//caca

    {1.456,0.015*stddelta,1.456-0.015*stddelta,1.456+0.015*stddelta},//nca gly  
    {1.514,0.016*stddelta,1.514-0.016*stddelta,1.514+0.016*stddelta},//cac gly  
    {1.326,0.018*stddelta,1.326-0.018*stddelta,1.326+0.018*stddelta},//cn gly  

    {1.468,0.017*stddelta,1.468-0.017*stddelta,1.468+0.017*stddelta},//nca pro  
    {1.524,0.020*stddelta,1.524-0.020*stddelta,1.524+0.020*stddelta},//cac pro  
    {1.338,0.019*stddelta,1.338-0.019*stddelta,1.338+0.019*stddelta},//cn pro  

    {1.229,0.019*stddelta,1.229-0.019*stddelta,1.229+0.019*stddelta},//co 
    {1.232,0.016*stddelta,1.232-0.016*stddelta,1.232+0.016*stddelta},//co gly  
    {1.228,0.020*stddelta,1.228-0.020*stddelta,1.228+0.020*stddelta},//co pro  
  };

  for(i=0;i<NUM_BONDS;i++){
    pChainFlex->numBondLengthOutliers[i]=0;   
  }
  for(int i=1; i<ChainFlexGetResidueCount(pChainFlex); i++){
    ResidueFlex* pResidueFlex = ChainFlexGetResidue(pChainFlex,i);
    Residue* pResidue = ChainGetResidue(pChain,i); 

    for(int j=0;j<NUM_ANGLES;j++){
      pResidueFlex->isBondLengthOutlierAbove[j] = FALSE;
      pResidueFlex->isBondLengthOutlierBelow[j] = FALSE;
    } 
 
    if(strcmp(ResidueGetName(pResidue), "GLY") == 0) istart=4;
    else if(strcmp(ResidueGetName(pResidue), "PRO") == 0) istart=7;
    else istart=0;

    if(pResidueFlex->len[0]<stdlength[istart+2][2]){
      pChainFlex->numBondLengthOutliers[0]++;
      pResidueFlex->isBondLengthOutlierBelow[0] = TRUE;
    }
    else if(pResidueFlex->len[0]>stdlength[istart+2][3]){
      pChainFlex->numBondLengthOutliers[0]++;
      pResidueFlex->isBondLengthOutlierAbove[0] = TRUE;
    }

    if(pResidueFlex->len[1]<stdlength[istart+0][2]){
      pChainFlex->numBondLengthOutliers[1]++;
      pResidueFlex->isBondLengthOutlierBelow[1] = TRUE;
    }
    else if(pResidueFlex->len[1]>stdlength[istart+0][3]){
      pChainFlex->numBondLengthOutliers[1]++;
      pResidueFlex->isBondLengthOutlierAbove[1] = TRUE;
    }

    if(pResidueFlex->len[2]<stdlength[istart+1][2]){
      pChainFlex->numBondLengthOutliers[2]++;
      pResidueFlex->isBondLengthOutlierBelow[2] = TRUE;
    }
    else if(pResidueFlex->len[2]>stdlength[istart+1][3]){
      pChainFlex->numBondLengthOutliers[2]++;
      pResidueFlex->isBondLengthOutlierAbove[2] = TRUE;
    }
  
    //CA i i-1 distance
    if(pResidueFlex->len[3]<stdlength[3][2]){
      pChainFlex->numBondLengthOutliers[3]++;
      pResidueFlex->isBondLengthOutlierBelow[3] = TRUE;
    }
    else if(pResidueFlex->len[3]>stdlength[3][3]){
      pChainFlex->numBondLengthOutliers[3]++;
      pResidueFlex->isBondLengthOutlierAbove[3] = TRUE;
    }
  }
  return Success;
}

int getRamaOutliers(Chain *pChain,ChainFlex *pChainFlex,int ***ramaduke){
  int i;
  int ti,tj;
  int aatype;

  pChainFlex->numRamaOutliers=0;   
  for(int i=1; i<ChainFlexGetResidueCount(pChainFlex)-1; i++){
    ResidueFlex* pResidueFlex = ChainFlexGetResidue(pChainFlex,i);
    Residue* pResidue = ChainGetResidue(pChain,i); 
    Residue* pResidueNext = ChainGetResidue(pChain,i+1); 
    pResidueFlex->isRamaOutlier = FALSE;

    if(strcmp(ResidueGetName(pResidue), "GLY") == 0) aatype=1;
    else if(strcmp(ResidueGetName(pResidue), "PRO") == 0) aatype=2;
    else if(strcmp(ResidueGetName(pResidueNext), "PRO") == 0) aatype=3;
    else aatype=0;
  
    //printf("phi %lf psi %lf\n",pResidue->phipsi[0],pResidue->phipsi[1]);

    ti=int(pResidue->phipsi[0]+180);
    if(ti<0) ti=0;
    else if(ti>359) ti=359;
    tj=int(pResidueNext->phipsi[1]+180);
    if(tj<0) tj=0;
    else if(tj>359) tj=359;
   
    //printf("ti %d tj %d ramaduke %d\n",ti,tj,ramaduke[aatype][ti][tj]);

    if(ramaduke[aatype][ti][tj]>0){      
      pChainFlex->numRamaOutliers++;
      pResidueFlex->isRamaOutlier = TRUE; 
    }
  }
  return Success;
}


int getClashes(Chain *pChain,ChainFlex *pChainFlex){
  int i,j,ii,jj;
  int totclash=0;
  int numseq = ChainGetResidueCount(pChain);
  XYZ disti[4],distj[4];
  XYZ distVec;
  double distance;
  bool *flagpos=new bool[numseq];
  for(i=0;i<numseq;i++){
    flagpos[i]=false;
  }

  double vdwds[8][8]={
    12.96000, 5.29000,13.69000, 8.41000,12.25000,10.24000,12.25000, 1.00000,
    13.69000, 6.25000,12.25000, 6.25000,12.25000, 4.41000,12.96000, 1.00000,
     5.29000, 1.44000, 7.29000, 7.29000, 5.29000, 3.24000, 4.84000, 1.00000,
     6.25000, 4.41000, 5.76000, 5.29000, 6.76000, 2.89000, 4.84000, 1.00000,
    12.25000, 5.29000,12.25000, 7.84000,10.89000, 5.29000,11.56000, 1.00000,
    12.25000, 4.84000, 9.00000, 3.24000,12.25000, 4.84000,12.25000, 1.00000,
    12.25000, 4.84000,10.89000, 7.84000,10.89000, 4.41000,12.25000, 1.00000,
     1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000};

  pChainFlex->numClashes = 0;
  for(i=0;i<numseq;i++){
    Residue *pResidue1 = ChainGetResidue(pChain,i);
    Atom *pAtomCA1 = ResidueGetAtomByName(pResidue1, "CA");
    Atom *pAtomC1 = ResidueGetAtomByName(pResidue1, "C");
    Atom *pAtomN1 = ResidueGetAtomByName(pResidue1, "N");
    Atom *pAtomO1 = ResidueGetAtomByName(pResidue1, "O");

    disti[0]=pAtomCA1->xyz;
    disti[1]=pAtomN1->xyz;
    disti[2]=pAtomC1->xyz;
    disti[3]=pAtomO1->xyz;

    for(j=i+1;j<numseq;j++){
      Residue *pResidue2 = ChainGetResidue(pChain,j);
      Atom *pAtomCA2 = ResidueGetAtomByName(pResidue2, "CA");
      Atom *pAtomC2 = ResidueGetAtomByName(pResidue2, "C");
      Atom *pAtomN2 = ResidueGetAtomByName(pResidue2, "N");
      Atom *pAtomO2 = ResidueGetAtomByName(pResidue2, "O");

      distj[0]=pAtomCA2->xyz;
      distj[1]=pAtomN2->xyz;
      distj[2]=pAtomC2->xyz;
      distj[3]=pAtomO2->xyz;

      for(ii=0;ii<4;ii++){
        for(jj=0;jj<4;jj++){
          if(ii!=0 || jj!=0) continue;
          distVec=XYZDifference(&disti[ii],&distj[jj]);
          distance=distVec.X*distVec.X+distVec.Y*distVec.Y+distVec.Z*distVec.Z;
          distance-=vdwds[ii][jj];
          if(distance<-0.01){    
            flagpos[i]=true;
            flagpos[j]=true;
            totclash++;
          }
        }
      } 
    }
  }

  for(i=0;i<numseq;i++){
    if(flagpos[i]){
      pChainFlex->numClashes++;
      ResidueFlex *pResidueFlex = ChainFlexGetResidue(pChainFlex,i);
      pResidueFlex->isClashed=TRUE; 
    }
  }
  delete[]flagpos;
  return Success;
}

int getMoveType(double tmov[],int totmov,double trandnum){
  int i;
  for(i=0;i<totmov;i++){
    if(trandnum>=tmov[i] && trandnum<tmov[i+1]){
      return i;
    }
  }
  return 0;
}

BOOL makeConformationalMove(Structure *decstr,StructureFlex *pStructureFlex,int chains[],int numchains,double &oldenergy,double &newenergy,
                            BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,
                            double *accdis,double ***phipsiprob,float ****cancbins,sssegment *sse,int numsse,paircont *bbptb,
                            int numbbptb,int *alphasind,int numsalpha,pairaa *natDistRestr,double **contactCA,double **contactCB){
  //int totmov=13;
  int totmov=12;
  //double tmov[15]={0.00,0.10,0.25,0.33,0.38,0.43,0.47,0.57,0.67,0.77,0.85,0.90,0.93,1.00,1.0000};
  //double tmov[15]={0.00,0.10,0.25,0.33,0.38,0.43,0.47,0.57,0.67,0.77,0.85,0.90,0.93,1.00,1.0001};
  double tmov[14]={0.00,0.15,0.30,0.38,0.43,0.48,0.52,0.62,0.72,0.82,0.90,0.95,1.00,1.0001};
  double movetype=Random();
  int ntype=getMoveType(tmov,totmov,movetype);
  bool flagacc=true;
  bool flagres;
  //ntype=10;
  //if(ntype==11)
  
 printf("flagppi %d\n",FLAG_PPI);

  ntype=0;
  switch(ntype)
  {
    case 0:
      flagres=moveAngle(decstr,pStructureFlex,chains,numchains,oldenergy,newenergy,flagacc,bbrotlib,atomParams,resiTopos,ramaduke,natDistRestr,contactCA,contactCB);
      //flagres=moveBackrub(decstr,pStructureFlex,chains,numchains,oldenergy,newenergy,flagacc,bbrotlib,atomParams,resiTopos,ramaduke,natDistRestr,contactCA,contactCB);
      break;
    case 1:
      flagres=moveLength(decstr,pStructureFlex,chains,numchains,oldenergy,newenergy,flagacc,bbrotlib,atomParams,resiTopos,ramaduke,natDistRestr,contactCA,contactCB);
      break;
    case 2:
      flagres=moveOmega(decstr,pStructureFlex,chains,numchains,oldenergy,newenergy,flagacc,bbrotlib,atomParams,resiTopos,ramaduke,natDistRestr,contactCA,contactCB);
      break;
    case 3:
      flagres=moveTor(decstr,pStructureFlex,chains,numchains,oldenergy,newenergy,flagacc,bbrotlib,atomParams,resiTopos,ramaduke,accdis,natDistRestr,contactCA,contactCB);
      break;
    case 4:
      flagres=movePhi(decstr,pStructureFlex,chains,numchains,oldenergy,newenergy,flagacc,bbrotlib,atomParams,resiTopos,ramaduke,accdis,phipsiprob,natDistRestr,contactCA,contactCB); //can remove accdis
      break;
    case 5:
      flagres=movePsi(decstr,pStructureFlex,chains,numchains,oldenergy,newenergy,flagacc,bbrotlib,atomParams,resiTopos,ramaduke,accdis,phipsiprob,natDistRestr,contactCA,contactCB); //can remove accdis
      break;
    case 6:
      flagres=moveBackrub(decstr,pStructureFlex,chains,numchains,oldenergy,newenergy,flagacc,bbrotlib,atomParams,resiTopos,ramaduke,natDistRestr,contactCA,contactCB);
      //flagres=moveAnc(decstr,pStructureFlex,chains,numchains,oldenergy,newenergy,flagacc,bbrotlib,atomParams,resiTopos,ramaduke,accdis,phipsiprob,cancbins,natDistRestr,contactCA,contactCB);
      //flagres=moveAngle(decstr,pStructureFlex,chains,numchains,oldenergy,newenergy,flagacc,bbrotlib,atomParams,resiTopos,ramaduke,natDistRestr,contactCA,contactCB);
      break;
    case 7:
      flagres=moveRotation(decstr,pStructureFlex,chains,numchains,oldenergy,newenergy,flagacc,bbrotlib,atomParams,resiTopos,ramaduke,natDistRestr,contactCA,contactCB);
      break;
    case 8:
      flagres=moveLMP(decstr,pStructureFlex,chains,numchains,oldenergy,newenergy,flagacc,bbrotlib,atomParams,resiTopos,ramaduke,natDistRestr,contactCA,contactCB);
      break;
    case 9:
      flagres=moveShift(decstr,pStructureFlex,chains,numchains,oldenergy,newenergy,flagacc,bbrotlib,atomParams,resiTopos,ramaduke,natDistRestr,contactCA,contactCB);
      break;
    case 10:
      flagres=movePos(decstr,pStructureFlex,chains,numchains,oldenergy,newenergy,flagacc,bbrotlib,atomParams,resiTopos,ramaduke,natDistRestr,contactCA,contactCB);
      break;
    case 11:
      flagres=moveHelix(decstr,pStructureFlex,chains,numchains,oldenergy,newenergy,flagacc,bbrotlib,atomParams,resiTopos,ramaduke,natDistRestr,contactCA,contactCB);
      break;
//    case 12:
//      flagres=moveBetaPack(decstr,pStructureFlex,chains,numchains,oldenergy,newenergy,flagacc,bbrotlib,atomParams,resiTopos,ramaduke,bbptb,numbbptb,natDistRestr,contactCA,contactCB);
      //flagres=moveTor(decstr,pStructureFlex,chains,numchains,oldenergy,newenergy,flagacc,bbrotlib,atomParams,resiTopos,ramaduke,accdis,natDistRestr,contactCA,contactCB);
//      break;
    case 12:
   //   flagres=moveHelixPack(decstr,pStructureFlex,chains,numchains,oldenergy,newenergy,flagacc,bbrotlib,atomParams,resiTopos,ramaduke,sse,alphasind,numsalpha,natDistRestr,contactCA,contactCB);
      break;
    default :
      printf("wrong move type %d\n",ntype);
      break;
  }
  //  *movtype=ntype;
  //int i;
  //  if(flagres){
  //    for(i=0;i<30;i++){
  //      enelistbk[i]=enelist[i];
  //    }
  //  }
  
  printf("old energy: %lf new energy: %lf\n",oldenergy,newenergy);
  return flagres;
}

BOOL moveAngle(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
               bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,pairaa *natDistRestr,
               double **contactCA,double **contactCB){
  const int CA_INDEX=3;
  int i,j;
  int trandpos;
  int numSideChainOptCycles=1;
  BOOL isCAOutlier,flagTor; 
  double meanstda[][4]={
    {117.2,2.2*4.0,117.2-2.2*4.0,117.2+2.2*4.0},//cacn
    {121.7,2.5*4.0,121.7-2.5*4.0,121.7+2.5*4.0},//cnca
    {111.0,2.7*4.0,111.0-2.7*4.0,111.0+2.7*4.0},//ncac
  };
  double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
//  double energyTermsBind[MAX_EVOEF_ENERGY_TERM_NUM]={0};

  int angtype=int(Random()*3.0);
  int wpos;
  int selectChain = deschn[int(Random()*numchains)];

  ChainFlex* pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
  Chain* pChain = StructureGetChain(decstr,selectChain);

  int totwrong = pChainFlex->numAngleOutliers[angtype] + pChainFlex->numBondLengthOutliers[CA_INDEX];
  if(totwrong>0 && Random()<0.1+totwrong/double(ChainFlexGetResidueCount(pChainFlex))){
    wpos=int(Random()*totwrong);
    i=-1;
    isCAOutlier=FALSE;
    for(trandpos=0; trandpos<ChainFlexGetResidueCount(pChainFlex); trandpos++){
      ResidueFlex *pResidueFlex = ChainFlexGetResidue(pChainFlex,trandpos);

      if(pResidueFlex->isRamaOutlier) i++;
      if(i==wpos) break;
      if(pResidueFlex->isAngleOutlierAbove[angtype] || pResidueFlex->isAngleOutlierBelow[angtype]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isBondLengthOutlierAbove[CA_INDEX] || pResidueFlex->isBondLengthOutlierBelow[CA_INDEX]){
        i++;
        isCAOutlier=TRUE;
      }
      if(i==wpos) break;
    }
    if(trandpos<0 || trandpos>=ChainFlexGetResidueCount(pChainFlex)) printf("fatal error %d ang\n",trandpos);
    if(isCAOutlier && trandpos>0 && angtype==2) trandpos--;
  }
  else{
    trandpos=int(ChainFlexGetResidueCount(pChainFlex)*Random());
    ResidueFlex* pResidueFlex = ChainFlexGetResidue(pChainFlex,trandpos);
    if(pResidueFlex->ssm!='C' && Random()>0.1){  
      //nummov[angtype][1]++;
      return FALSE;
    }
  }
  //StructureFlex structureFlexTemp;
  //for(i=0;i<StructureFlexGetChainCount(pStructureFlex);i++){
    //ChainFlex *pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlex chainFlexTemp;
    ChainFlexCreate(&chainFlexTemp);
    ChainFlexCopy(&chainFlexTemp,pChainFlex);
    //StructureFlexAddChain(structureFlexTemp, &chainFlex);
  //}

  Structure structureTemp;
  StructureCreate(&structureTemp);
  for(i=0;i<StructureGetChainCount(decstr);i++){
    Chain *pChainSelect = StructureGetChain(decstr,i);
    Chain chainTemp;
    ChainCreate(&chainTemp);
    ChainCopy(&chainTemp,pChainSelect);
    StructureAddChain(&structureTemp, &chainTemp);
    ChainDestroy(&chainTemp);
  }

  
  if(trandpos==0) {trandpos=1; printf("trandpos is 0\n");}


  Chain *pChainTemp = StructureGetChain(&structureTemp,selectChain);
  ResidueFlex *selectResidueFlex = ChainFlexGetResidue(&chainFlexTemp,trandpos);
  selectResidueFlex->ang[angtype]=meanstda[angtype][0]+meanstda[angtype][1]*(2*Random()-1);
  flagTor=NewMainChainPos(StructureGetChain(&structureTemp,selectChain),&chainFlexTemp,trandpos);


  if(trandpos==0){   //should rebuild all chains in structure to fix problem with NewMainChainPos
    for(int j=trandpos; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  }
  else{
    for(int j=trandpos-1; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  }


  //TODO: optimize interface sidechains
  OptimizeSideChains(&structureTemp,bbrotlib,numSideChainOptCycles,selectChain); 
  StructureCalcEnergyFlex(&structureTemp,bbrotlib,energyTerms);
  EnergyTermWeighting(energyTerms);

  newenergy=energyTerms[0];
  //double constraintEnergy = calcConstraintEnergy(pChainTemp,natDistRestr)*200;
//  newenergy+=calcConstraintEnergy(pChainTemp,natDistRestr)*75;
  //double caEnergy,cbEnergy;
//  newenergy+=calcDistanceEnergy(pChainTemp,distanceCB,distanceWeight);
  //calcContactEnergy(pChainTemp,contactCA,contactCB,caEnergy,cbEnergy);
  //newenergy+=caEnergy;
  //newenergy+=cbEnergy;
  //std::cout << "New Energy: " << newenergy << " Old Energy: " << oldenergy << " Constraint Energy: " << constraintEnergy << std::endl;
  if(newenergy<oldenergy){
    oldenergy=newenergy;
    //nummov[angtype][0]++;
    ChainFlex *origChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlexCopy(origChainFlex,&chainFlexTemp);
    Chain *origChain = StructureGetChain(decstr,selectChain);
    ChainCopy(origChain,StructureGetChain(&structureTemp,selectChain));
    getRamaOutliers(origChain,origChainFlex,ramaduke);
    getAngleOutliers(origChain,origChainFlex);
    getBondLengthOutliers(origChain,origChainFlex);
    getClashes(origChain,origChainFlex);
    StructureDestroy(&structureTemp);
    return TRUE;
  }
  else{
    StructureDestroy(&structureTemp);
    return FALSE;
  }
}

BOOL moveLength(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
               bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,pairaa *natDistRestr,
               double **contactCA,double **contactCB){
  const int CA_INDEX=3;
  int i,j;
  int trandpos;
  int numSideChainOptCycles=1;
  BOOL isCAOutlier,flagTor; 
  double meanstd[][4]={
    {1.336,0.023*4,1.336-0.023*4,1.336+0.023*4},//cn  
    {1.459,0.020*4,1.459-0.020*4,1.459+0.020*4},//nca  
    {1.525,0.026*4,1.525-0.026*4,1.525+0.026*4},//cac  
  };
  double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
  
  int lentype=int(Random()*3.0);
  int wpos;
  int selectChain = deschn[int(Random()*numchains)];

  ChainFlex* pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
  Chain* pChain = StructureGetChain(decstr,selectChain);
  int totwrong = pChainFlex->numBondLengthOutliers[lentype] + pChainFlex->numBondLengthOutliers[CA_INDEX];
  if(totwrong>0 && Random()<0.1+totwrong/double(ChainFlexGetResidueCount(pChainFlex))){
    wpos=int(Random()*totwrong);
    i=-1;
    isCAOutlier=FALSE;
    for(trandpos=0; trandpos<ChainFlexGetResidueCount(pChainFlex); trandpos++){
      ResidueFlex *pResidueFlex = ChainFlexGetResidue(pChainFlex,trandpos);
      if(pResidueFlex->isBondLengthOutlierAbove[lentype] || pResidueFlex->isBondLengthOutlierBelow[lentype]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isBondLengthOutlierAbove[CA_INDEX] || pResidueFlex->isBondLengthOutlierBelow[CA_INDEX]){
        i++;
        isCAOutlier=TRUE;
      }
      if(i==wpos) break;
    }
    if(trandpos<0 || trandpos>=ChainFlexGetResidueCount(pChainFlex)) printf("fatal error %d length\n",trandpos);
    if(isCAOutlier && trandpos>0 && lentype==2) trandpos--;
  }
  else{
    trandpos=int(ChainFlexGetResidueCount(pChainFlex)*Random());
    ResidueFlex* pResidueFlex = ChainFlexGetResidue(pChainFlex,trandpos);
    if(pResidueFlex->ssm!='C' && Random()>0.1){  
      //nummov[angtype][1]++;
      return FALSE;
    }
  }
  //StructureFlex structureFlexTemp;
  //for(i=0;i<StructureFlexGetChainCount(pStructureFlex);i++){
    //ChainFlex *pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlex chainFlexTemp;
    ChainFlexCreate(&chainFlexTemp);
    ChainFlexCopy(&chainFlexTemp,pChainFlex);
    //StructureFlexAddChain(structureFlexTemp, &chainFlex);
  //}

  Structure structureTemp;
  StructureCreate(&structureTemp);
  for(i=0;i<StructureGetChainCount(decstr);i++){
    Chain *pChainSelect = StructureGetChain(decstr,i);
    Chain chainTemp;
    ChainCreate(&chainTemp);
    ChainCopy(&chainTemp,pChainSelect);
    StructureAddChain(&structureTemp, &chainTemp);
    ChainDestroy(&chainTemp);
  }
  Chain *pChainTemp = StructureGetChain(&structureTemp,selectChain);

  ResidueFlex *selectResidueFlex = ChainFlexGetResidue(&chainFlexTemp,trandpos);
  selectResidueFlex->len[lentype]=meanstd[lentype][0]+meanstd[lentype][1]*(2*Random()-1);
  flagTor=NewMainChainPos(StructureGetChain(&structureTemp,selectChain),&chainFlexTemp,trandpos);
  if(trandpos==0){   //should rebuild all chains in structure to fix problem with NewMainChainPos
    for(int j=trandpos; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  }
  else{
    for(int j=trandpos-1; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  }

  //TODO: optimize interface sidechains
  OptimizeSideChains(&structureTemp,bbrotlib,numSideChainOptCycles,selectChain); 
  EnergyTermInitialize(energyTerms);
  StructureCalcEnergyFlex(&structureTemp,bbrotlib,energyTerms);
  EnergyTermWeighting(energyTerms);
  newenergy=energyTerms[0];
//  newenergy+=calcConstraintEnergy(pChainTemp,natDistRestr)*75;
//  newenergy+=calcDistanceEnergy(pChainTemp,distanceCB,distanceWeight);
  if(newenergy<oldenergy){
    oldenergy=newenergy;
    //nummov[angtype][0]++;
    ChainFlex *origChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlexCopy(origChainFlex,&chainFlexTemp);
    Chain *origChain = StructureGetChain(decstr,selectChain);
    ChainCopy(origChain,StructureGetChain(&structureTemp,selectChain));
    getRamaOutliers(origChain,origChainFlex,ramaduke);
    getAngleOutliers(origChain,origChainFlex);
    getBondLengthOutliers(origChain,origChainFlex);
    getClashes(origChain,origChainFlex);
    StructureDestroy(&structureTemp);
    return TRUE;
  }
  else{
    StructureDestroy(&structureTemp);
    return FALSE;
  }
}

BOOL moveOmega(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
               bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,pairaa *natDistRestr,
               double **contactCA,double **contactCB){
  const int CA_INDEX=3;
  int i,j;
  int trandpos;
  int numSideChainOptCycles=1;
  BOOL isCAOutlier,flagTor; 
  double trand2;
  double meanstd[][4]={
    {1.336,0.023*4,1.336-0.023*4,1.336+0.023*4},//cn  
    {1.459,0.020*4,1.459-0.020*4,1.459+0.020*4},//nca  
    {1.525,0.026*4,1.525-0.026*4,1.525+0.026*4},//cac  
  };
  double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
  
  int wpos;
  int selectChain = deschn[int(Random()*numchains)];

  ChainFlex* pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
  Chain* pChain = StructureGetChain(decstr,selectChain);
  int totwrong = pChainFlex->numBondLengthOutliers[CA_INDEX];
  if(totwrong>0 && Random()<0.1+totwrong/double(ChainFlexGetResidueCount(pChainFlex))){
    wpos=int(Random()*totwrong);
    i=-1;
    isCAOutlier=FALSE;
    for(trandpos=0; trandpos<ChainFlexGetResidueCount(pChainFlex); trandpos++){
      ResidueFlex *pResidueFlex = ChainFlexGetResidue(pChainFlex,trandpos);
      if(pResidueFlex->isBondLengthOutlierAbove[CA_INDEX] || pResidueFlex->isBondLengthOutlierBelow[CA_INDEX]){
        i++;
        isCAOutlier=TRUE;
      }
      if(i==wpos) break;
    }
    if(trandpos<0 || trandpos>=ChainFlexGetResidueCount(pChainFlex)) printf("fatal error %d omega\n",trandpos);
  }
  else{
    trandpos=int(ChainFlexGetResidueCount(pChainFlex)*Random());
    ResidueFlex* pResidueFlex = ChainFlexGetResidue(pChainFlex,trandpos);
    if(pResidueFlex->ssm!='C' && Random()>0.1){  
      //nummov[angtype][1]++;
      return FALSE;
    }
  }
  //StructureFlex structureFlexTemp;
  //for(i=0;i<StructureFlexGetChainCount(pStructureFlex);i++){
    //ChainFlex *pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlex chainFlexTemp;
    ChainFlexCreate(&chainFlexTemp);
    ChainFlexCopy(&chainFlexTemp,pChainFlex);
    //StructureFlexAddChain(structureFlexTemp, &chainFlex);
  //}

  Structure structureTemp;
  StructureCreate(&structureTemp);
  for(i=0;i<StructureGetChainCount(decstr);i++){
    Chain *pChainSelect = StructureGetChain(decstr,i);
    Chain chainTemp;
    ChainCreate(&chainTemp);
    ChainCopy(&chainTemp,pChainSelect);
    StructureAddChain(&structureTemp, &chainTemp);
    ChainDestroy(&chainTemp);
  }
  Chain *pChainTemp = StructureGetChain(&structureTemp,selectChain);

  ResidueFlex *selectResidueFlex = ChainFlexGetResidue(&chainFlexTemp,trandpos);

  trand2=double(16*Random())-8;
  selectResidueFlex->dih[1]+=trand2;
  if(selectResidueFlex->dih[1]<0){
    selectResidueFlex->dih[1]+=360;
  }
  else if(selectResidueFlex->dih[1]>=360){
    selectResidueFlex->dih[1]-=360;
  }
  if(selectResidueFlex->dih[1]>190 || selectResidueFlex->dih[1]<170){
    //nummov[6][1]++;
    return FALSE;
  }

  flagTor=NewMainChainPos(StructureGetChain(&structureTemp,selectChain),&chainFlexTemp,trandpos);
  if(trandpos==0){   //should rebuild all chains in structure to fix problem with NewMainChainPos
    for(int j=trandpos; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  }
  else{
    for(int j=trandpos-1; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  }

  //TODO: optimize interface sidechains
  OptimizeSideChains(&structureTemp,bbrotlib,numSideChainOptCycles,selectChain); 
  EnergyTermInitialize(energyTerms);
  StructureCalcEnergyFlex(&structureTemp,bbrotlib,energyTerms);
  EnergyTermWeighting(energyTerms);
  newenergy=energyTerms[0];
//  newenergy+=calcConstraintEnergy(pChainTemp,natDistRestr)*75;
//  newenergy+=calcDistanceEnergy(pChainTemp,distanceCB,distanceWeight);
  if(newenergy<oldenergy){
    oldenergy=newenergy;
    //nummov[angtype][0]++;
    ChainFlex *origChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlexCopy(origChainFlex,&chainFlexTemp);
    Chain *origChain = StructureGetChain(decstr,selectChain);
    ChainCopy(origChain,StructureGetChain(&structureTemp,selectChain));
    getRamaOutliers(origChain,origChainFlex,ramaduke);
    getAngleOutliers(origChain,origChainFlex);
    getBondLengthOutliers(origChain,origChainFlex);
    getClashes(origChain,origChainFlex); 
    StructureDestroy(&structureTemp);
    return TRUE;
  }
  else{
    StructureDestroy(&structureTemp);
    return FALSE;
  }
}

BOOL moveTor(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
             bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,
             double *accdis,pairaa *natDistRestr,double **contactCA,double **contactCB){
  const int CA_INDEX=3;
  int i,j;
  int trandpos;
  int numSideChainOptCycles=1;
  BOOL flagTor; 
  double trand2;
  double meanstd[][4]={
    {1.336,0.023*4,1.336-0.023*4,1.336+0.023*4},//cn  
    {1.459,0.020*4,1.459-0.020*4,1.459+0.020*4},//nca  
    {1.525,0.026*4,1.525-0.026*4,1.525+0.026*4},//cac  
  };
  double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
  
  int wpos;
  int selectChain = deschn[int(Random()*numchains)];

  ChainFlex* pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
  Chain* pChain = StructureGetChain(decstr,selectChain);

  int totwrong = pChainFlex->numRamaOutliers;
  if(totwrong>0 && Random()<0.1+totwrong/double(ChainFlexGetResidueCount(pChainFlex))){
    wpos=int(Random()*totwrong);
    i=-1;
    for(trandpos=0; trandpos<ChainFlexGetResidueCount(pChainFlex); trandpos++){
      ResidueFlex *pResidueFlex = ChainFlexGetResidue(pChainFlex,trandpos);
      if(pResidueFlex->isRamaOutlier) i++;
      if(i==wpos) break;
    }
    if(trandpos<0 || trandpos>=ChainFlexGetResidueCount(pChainFlex)) printf("fatal error %d torsion\n",trandpos);
  }
  else{
    trandpos=int(ChainFlexGetResidueCount(pChainFlex)*Random());
    ResidueFlex* pResidueFlex = ChainFlexGetResidue(pChainFlex,trandpos);
    if(pResidueFlex->ssm!='C' && Random()>0.1){  
      //nummov[angtype][1]++;
      return FALSE;
    }
  }

  if(trandpos==0 || trandpos==ChainFlexGetResidueCount(pChainFlex)-1){
    //nummov[7][1]++;
    return FALSE;
  }

  //StructureFlex structureFlexTemp;
  //for(i=0;i<StructureFlexGetChainCount(pStructureFlex);i++){
    //ChainFlex *pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlex chainFlexTemp;
    ChainFlexCreate(&chainFlexTemp);
    ChainFlexCopy(&chainFlexTemp,pChainFlex);
    //StructureFlexAddChain(structureFlexTemp, &chainFlex);
  //}

  Structure structureTemp;
  StructureCreate(&structureTemp);
  for(i=0;i<StructureGetChainCount(decstr);i++){
    Chain *pChainSelect = StructureGetChain(decstr,i);
    Chain chainTemp;
    ChainCreate(&chainTemp);
    ChainCopy(&chainTemp,pChainSelect);
    StructureAddChain(&structureTemp, &chainTemp);
    ChainDestroy(&chainTemp);
  }
  Chain *pChainTemp = StructureGetChain(&structureTemp,selectChain);

  ResidueFlex *selectResidueFlex = ChainFlexGetResidue(&chainFlexTemp,trandpos);
  ResidueFlex *selectResidueFlexNext = ChainFlexGetResidue(&chainFlexTemp,trandpos+1);
  
  int ii,jj;
  ii=posinarray(accdis,129600,Random());
  if(ii<0) ii=0;
  else if(ii>=129600) ii=129599;
  jj=ii/360;
  ii-=jj*360;
  selectResidueFlex->dih[2]=jj+0.5;
  selectResidueFlexNext->dih[0]=ii+0.5;

  flagTor=NewMainChainPos(StructureGetChain(&structureTemp,selectChain),&chainFlexTemp,trandpos);
  if(trandpos==0){   //should rebuild all chains in structure to fix problem with NewMainChainPos
    for(int j=trandpos; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  }
  else{
    for(int j=trandpos-1; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  }

  //TODO: optimize interface sidechains
  OptimizeSideChains(&structureTemp,bbrotlib,numSideChainOptCycles,selectChain); 
  EnergyTermInitialize(energyTerms);
  StructureCalcEnergyFlex(&structureTemp,bbrotlib,energyTerms);
  EnergyTermWeighting(energyTerms);
  newenergy=energyTerms[0];
//  newenergy+=calcConstraintEnergy(pChainTemp,natDistRestr)*75;
//  newenergy+=calcDistanceEnergy(pChainTemp,distanceCB,distanceWeight);

  if(newenergy<oldenergy){
    oldenergy=newenergy;
    //nummov[angtype][0]++;
    ChainFlex *origChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlexCopy(origChainFlex,&chainFlexTemp);
    Chain *origChain = StructureGetChain(decstr,selectChain);
    ChainCopy(origChain,StructureGetChain(&structureTemp,selectChain));
    getRamaOutliers(origChain,origChainFlex,ramaduke);
    getAngleOutliers(origChain,origChainFlex);
    getBondLengthOutliers(origChain,origChainFlex);
    getClashes(origChain,origChainFlex);
    StructureDestroy(&structureTemp);
    return TRUE;
  }
  else{
    StructureDestroy(&structureTemp);
    return FALSE;
  }
}

BOOL movePhi(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
             bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,
             double *accdis,double ***phipsiprob,pairaa *natDistRestr,double **contactCA,double **contactCB){
  const int CA_INDEX=3;
  int i,j;
  int trandpos;
  int numSideChainOptCycles=1;
  BOOL flagTor; 
  double trand2;
  double meanstd[][4]={
    {1.336,0.023*4,1.336-0.023*4,1.336+0.023*4},//cn  
    {1.459,0.020*4,1.459-0.020*4,1.459+0.020*4},//nca  
    {1.525,0.026*4,1.525-0.026*4,1.525+0.026*4},//cac  
  };
  double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
  
  int wpos;
  int selectChain = deschn[int(Random()*numchains)];

  ChainFlex* pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
  Chain* pChain = StructureGetChain(decstr,selectChain);
  int totwrong = pChainFlex->numRamaOutliers;
  if(totwrong>0 && Random()<0.1+totwrong/double(ChainFlexGetResidueCount(pChainFlex))){
    wpos=int(Random()*totwrong);
    i=-1;
    for(trandpos=0; trandpos<ChainFlexGetResidueCount(pChainFlex); trandpos++){
      ResidueFlex *pResidueFlex = ChainFlexGetResidue(pChainFlex,trandpos);
      if(pResidueFlex->isRamaOutlier) i++;
      if(i==wpos) break;
    }
    if(trandpos<0 || trandpos>=ChainFlexGetResidueCount(pChainFlex)) printf("fatal error %d phi\n",trandpos);
  }
  else{
    trandpos=int(ChainFlexGetResidueCount(pChainFlex)*Random());
    ResidueFlex* pResidueFlex = ChainFlexGetResidue(pChainFlex,trandpos);
    if(pResidueFlex->ssm!='C' && Random()>0.1){  
      //nummov[angtype][1]++;
      return FALSE;
    }
  }
  if(trandpos==0 || trandpos==ChainFlexGetResidueCount(pChainFlex)-1){
    //nummov[7][1]++;
    return FALSE;
  }


  //StructureFlex structureFlexTemp;
  //for(i=0;i<StructureFlexGetChainCount(pStructureFlex);i++){
    //ChainFlex *pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlex chainFlexTemp;
    ChainFlexCreate(&chainFlexTemp);
    ChainFlexCopy(&chainFlexTemp,pChainFlex);
    //StructureFlexAddChain(structureFlexTemp, &chainFlex);
  //}

  Structure structureTemp;
  StructureCreate(&structureTemp);
  for(i=0;i<StructureGetChainCount(decstr);i++){
    Chain *pChainSelect = StructureGetChain(decstr,i);
    Chain chainTemp;
    ChainCreate(&chainTemp);
    ChainCopy(&chainTemp,pChainSelect);
    StructureAddChain(&structureTemp, &chainTemp);
    ChainDestroy(&chainTemp);
  }
  Chain *pChainTemp = StructureGetChain(&structureTemp,selectChain);
  Residue *pResidue = ChainGetResidue(pChainTemp,trandpos);
  ResidueFlex *selectResidueFlex = ChainFlexGetResidue(&chainFlexTemp,trandpos);
  ResidueFlex *selectResidueFlexNext = ChainFlexGetResidue(&chainFlexTemp,trandpos+1);

  int inda=getResidueIndex(pResidue);
  int indp=findpos2(phipsiprob[3][inda],0,359,Random());
  selectResidueFlex->dih[2]=indp+0.5;
  if(selectResidueFlex->dih[2]<0){
    selectResidueFlex->dih[2]+=360;
  }
  else if(selectResidueFlex->dih[2]>=360){
    selectResidueFlex->dih[2]-=360;
  }

  flagTor=NewMainChainPos(StructureGetChain(&structureTemp,selectChain),&chainFlexTemp,trandpos);
  if(trandpos==0){   //should rebuild all chains in structure to fix problem with NewMainChainPos
    for(int j=trandpos; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  }
  else{
    for(int j=trandpos-1; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  }

  //TODO: optimize interface sidechains
  OptimizeSideChains(&structureTemp,bbrotlib,numSideChainOptCycles,selectChain); 
  EnergyTermInitialize(energyTerms);
  StructureCalcEnergyFlex(&structureTemp,bbrotlib,energyTerms);
  EnergyTermWeighting(energyTerms);
  newenergy=energyTerms[0];
//  newenergy+=calcConstraintEnergy(pChainTemp,natDistRestr)*75;
//  newenergy+=calcDistanceEnergy(pChainTemp,distanceCB,distanceWeight);
  if(newenergy<oldenergy){
    oldenergy=newenergy;
    //nummov[angtype][0]++;
    ChainFlex *origChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlexCopy(origChainFlex,&chainFlexTemp);
    Chain *origChain = StructureGetChain(decstr,selectChain);
    ChainCopy(origChain,StructureGetChain(&structureTemp,selectChain));
    getRamaOutliers(origChain,origChainFlex,ramaduke);
    getAngleOutliers(origChain,origChainFlex);
    getBondLengthOutliers(origChain,origChainFlex);
    getClashes(origChain,origChainFlex);
    StructureDestroy(&structureTemp);
    return TRUE;
  }
  else{
    StructureDestroy(&structureTemp);
    return FALSE;
  }
}

BOOL movePsi(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
             bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,
             double *accdis,double ***phipsiprob,pairaa *natDistRestr,double **contactCA,double **contactCB){
  const int CA_INDEX=3;
  int i,j;
  int trandpos;
  int numSideChainOptCycles=1;
  BOOL flagTor; 
  double trand2;
  double meanstd[][4]={
    {1.336,0.023*4,1.336-0.023*4,1.336+0.023*4},//cn  
    {1.459,0.020*4,1.459-0.020*4,1.459+0.020*4},//nca  
    {1.525,0.026*4,1.525-0.026*4,1.525+0.026*4},//cac  
  };
  double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
  
  int wpos;
  int selectChain = deschn[int(Random()*numchains)];

  ChainFlex* pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
  Chain* pChain = StructureGetChain(decstr,selectChain);
  int totwrong = pChainFlex->numRamaOutliers;
  if(totwrong>0 && Random()<0.1+totwrong/double(ChainFlexGetResidueCount(pChainFlex))){
    wpos=int(Random()*totwrong);
    i=-1;
    for(trandpos=0; trandpos<ChainFlexGetResidueCount(pChainFlex); trandpos++){
      ResidueFlex *pResidueFlex = ChainFlexGetResidue(pChainFlex,trandpos);
      if(pResidueFlex->isRamaOutlier) i++;
      if(i==wpos) break;
    }
    if(trandpos<0 || trandpos>=ChainFlexGetResidueCount(pChainFlex)) printf("fatal error %d phi\n",trandpos);
  }
  else{
    trandpos=int(ChainFlexGetResidueCount(pChainFlex)*Random());
    ResidueFlex* pResidueFlex = ChainFlexGetResidue(pChainFlex,trandpos);
    if(pResidueFlex->ssm!='C' && Random()>0.1){  
      //nummov[angtype][1]++;
      return FALSE;
    }
  }
  if(trandpos==0 || trandpos==ChainFlexGetResidueCount(pChainFlex)-1){
    //nummov[7][1]++;
    return FALSE;
  }


  //StructureFlex structureFlexTemp;
  //for(i=0;i<StructureFlexGetChainCount(pStructureFlex);i++){
    //ChainFlex *pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlex chainFlexTemp;
    ChainFlexCreate(&chainFlexTemp);
    ChainFlexCopy(&chainFlexTemp,pChainFlex);
    //StructureFlexAddChain(structureFlexTemp, &chainFlex);
  //}

  Structure structureTemp;
  StructureCreate(&structureTemp);
  for(i=0;i<StructureGetChainCount(decstr);i++){
    Chain *pChainSelect = StructureGetChain(decstr,i);
    Chain chainTemp;
    ChainCreate(&chainTemp);
    ChainCopy(&chainTemp,pChainSelect);
    StructureAddChain(&structureTemp, &chainTemp);
    ChainDestroy(&chainTemp);
  }
  Chain *pChainTemp = StructureGetChain(&structureTemp,selectChain);
  Residue *pResidue = ChainGetResidue(pChainTemp,trandpos);
  ResidueFlex *selectResidueFlex = ChainFlexGetResidue(&chainFlexTemp,trandpos);
  ResidueFlex *selectResidueFlexNext = ChainFlexGetResidue(&chainFlexTemp,trandpos+1);

  int inda=getResidueIndex(pResidue);
  int indp=findpos2(phipsiprob[7][inda],0,359,Random());
  selectResidueFlex->dih[0]=indp+0.5;
  if(selectResidueFlex->dih[0]<0){
    selectResidueFlex->dih[0]+=360;
  }
  else if(selectResidueFlex->dih[0]>=360){
    selectResidueFlex->dih[0]-=360;
  }

  flagTor=NewMainChainPos(StructureGetChain(&structureTemp,selectChain),&chainFlexTemp,trandpos);
  if(trandpos==0){   //should rebuild all chains in structure to fix problem with NewMainChainPos
    for(int j=trandpos; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  }
  else{
    for(int j=trandpos-1; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  }

  //TODO: optimize interface sidechains
  OptimizeSideChains(&structureTemp,bbrotlib,numSideChainOptCycles,selectChain); 
  EnergyTermInitialize(energyTerms);
  StructureCalcEnergyFlex(&structureTemp,bbrotlib,energyTerms);
  EnergyTermWeighting(energyTerms);
  newenergy=energyTerms[0];
//  newenergy+=calcConstraintEnergy(pChainTemp,natDistRestr)*75;
//  newenergy+=calcDistanceEnergy(pChainTemp,distanceCB,distanceWeight);
  if(newenergy<oldenergy){
    oldenergy=newenergy;
    //nummov[angtype][0]++;
    ChainFlex *origChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlexCopy(origChainFlex,&chainFlexTemp);
    Chain *origChain = StructureGetChain(decstr,selectChain);
    ChainCopy(origChain,StructureGetChain(&structureTemp,selectChain));
    getRamaOutliers(origChain,origChainFlex,ramaduke);
    getAngleOutliers(origChain,origChainFlex);
    getBondLengthOutliers(origChain,origChainFlex);
    getClashes(origChain,origChainFlex);
    StructureDestroy(&structureTemp);
    return TRUE;
  }
  else{
    StructureDestroy(&structureTemp);
    return FALSE;
  }
}


BOOL moveAnc(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
             bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,
             double *accdis,double ***phipsiprob,float ****cancbins,pairaa *natDistRestr,double **contactCA,double **contactCB){
  const int CA_INDEX=3;
  int i,j;
  int trandpos;
  int numSideChainOptCycles=1;
  BOOL flagTor; 
  double trand2;
  double meanstd[][4]={
    {1.336,0.023*4,1.336-0.023*4,1.336+0.023*4},//cn  
    {1.459,0.020*4,1.459-0.020*4,1.459+0.020*4},//nca  
    {1.525,0.026*4,1.525-0.026*4,1.525+0.026*4},//cac  
  };
  double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
  double raddeg = PI/180.0;
  
  int wpos;
  int selectChain = deschn[int(Random()*numchains)];

  ChainFlex* pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
  Chain* pChain = StructureGetChain(decstr,selectChain);
  int totwrong = pChainFlex->numRamaOutliers + pChainFlex->numAngleOutliers[0] + pChainFlex->numAngleOutliers[1] 
               + pChainFlex->numAngleOutliers[2] + pChainFlex->numBondLengthOutliers[0] + pChainFlex->numBondLengthOutliers[1]
               + pChainFlex->numBondLengthOutliers[2] + pChainFlex->numRamaOutliers;
  if(totwrong>0 && Random()<0.1+totwrong/double(ChainFlexGetResidueCount(pChainFlex))){
    wpos=int(Random()*totwrong);
    i=-1;
    for(trandpos=0; trandpos<ChainFlexGetResidueCount(pChainFlex); trandpos++){
      ResidueFlex *pResidueFlex = ChainFlexGetResidue(pChainFlex,trandpos);
      if(pResidueFlex->isAngleOutlierAbove[0] || pResidueFlex->isAngleOutlierBelow[0]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isAngleOutlierAbove[1] || pResidueFlex->isAngleOutlierBelow[1]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isAngleOutlierAbove[2] || pResidueFlex->isAngleOutlierBelow[2]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isBondLengthOutlierAbove[0] || pResidueFlex->isBondLengthOutlierBelow[0]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isBondLengthOutlierAbove[1] || pResidueFlex->isBondLengthOutlierBelow[1]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isBondLengthOutlierAbove[2] || pResidueFlex->isBondLengthOutlierBelow[2]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isRamaOutlier) i++;
      if(i==wpos) break;
    }
    trandpos-=1+int(2*Random());
    if(trandpos>ChainGetResidueCount(pChain)-4){
      trandpos=ChainGetResidueCount(pChain)-4;
    }
    else if(trandpos<0){
      trandpos=0;
    }
  }
  else{
    trandpos=int((ChainFlexGetResidueCount(pChainFlex)-3)*Random());
    ResidueFlex* pResidueFlex = ChainFlexGetResidue(pChainFlex,trandpos);
    if(pResidueFlex->ssm!='C' && Random()>0.1){  
      //nummov[angtype][1]++;
      return FALSE;
    }
  }
  if(trandpos==0 || trandpos==ChainFlexGetResidueCount(pChainFlex)-1){
    //nummov[7][1]++;
    return FALSE;
  }

  //StructureFlex structureFlexTemp;
  //for(i=0;i<StructureFlexGetChainCount(pStructureFlex);i++){
    //ChainFlex *pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlex chainFlexTemp;
    ChainFlexCreate(&chainFlexTemp);
    ChainFlexCopy(&chainFlexTemp,pChainFlex);
    //StructureFlexAddChain(structureFlexTemp, &chainFlex);
  //}

  Structure structureTemp;
  StructureCreate(&structureTemp);
  for(i=0;i<StructureGetChainCount(decstr);i++){
    Chain *pChainSelect = StructureGetChain(decstr,i);
    Chain chainTemp;
    ChainCreate(&chainTemp);
    ChainCopy(&chainTemp,pChainSelect);
    StructureAddChain(&structureTemp, &chainTemp);
    ChainDestroy(&chainTemp);
  }
  Chain *pChainTemp = StructureGetChain(&structureTemp,selectChain);
  Residue *pResidue1 = ChainGetResidue(pChainTemp,trandpos);
  Residue *pResidue2 = ChainGetResidue(pChainTemp,trandpos+1);
  Residue *pResidue3 = ChainGetResidue(pChainTemp,trandpos+2);
  Residue *pResidue4 = ChainGetResidue(pChainTemp,trandpos+3);
  Atom *pAtomCA1 = ResidueGetAtomByName(pResidue1, "CA");
  Atom *pAtomCA2 = ResidueGetAtomByName(pResidue2, "CA");  
  Atom *pAtomCA3 = ResidueGetAtomByName(pResidue3, "CA");
  Atom *pAtomCA4 = ResidueGetAtomByName(pResidue4, "CA");
  Atom *pAtomC2 = ResidueGetAtomByName(pResidue2, "C");
  Atom *pAtomN3 = ResidueGetAtomByName(pResidue3, "N");

  XYZ pout;
  double tdist;
  int ii,jj,kk;
  XYZ vector1 = XYZDifference(&pAtomCA2->xyz,&pAtomCA1->xyz);
  XYZ vector2 = XYZDifference(&pAtomCA2->xyz,&pAtomCA3->xyz);
  XYZ vector3 = XYZDifference(&pAtomCA3->xyz,&pAtomCA4->xyz); 
  XYZ vector4 = XYZDifference(&pAtomCA3->xyz,&pAtomCA2->xyz);
  tdist=XYZAngle(&vector1,&vector2);

  ii=int(tdist*50.0)-53;
  if(ii<0) ii=0;
  else if(ii>102) ii=102;
  tdist=XYZAngle(&vector3,&vector4);
  jj=int(tdist*50.0)-53;
  if(jj<0) jj=0;
  else if(jj>102) jj=102;
  tdist=calcTorsion(pAtomCA1->xyz.X,pAtomCA1->xyz.Y,pAtomCA1->xyz.Z,
                    pAtomCA2->xyz.X,pAtomCA2->xyz.Y,pAtomCA2->xyz.Z,
                    pAtomCA3->xyz.X,pAtomCA3->xyz.Y,pAtomCA3->xyz.Z,
                    pAtomCA4->xyz.X,pAtomCA4->xyz.Y,pAtomCA4->xyz.Z)*raddeg;
  kk=int(tdist*25.0);
  convertTorToPos(pAtomCA1->xyz.X,pAtomCA1->xyz.Y,pAtomCA1->xyz.Z,
                  pAtomCA3->xyz.X,pAtomCA3->xyz.Y,pAtomCA3->xyz.Z,
                  pAtomCA2->xyz.X,pAtomCA2->xyz.Y,pAtomCA2->xyz.Z,
                  cancbins[5][ii][jj][kk]*raddeg,cancbins[1][ii][jj][kk],cancbins[3][ii][jj][kk],
                  &pout.X,&pout.Y,&pout.Z);
  pAtomC2->xyz.X=pout.X;
  pAtomC2->xyz.Y=pout.Y;
  pAtomC2->xyz.Z=pout.Z;
  convertTorToPos(pAtomCA4->xyz.X,pAtomCA4->xyz.Y,pAtomCA4->xyz.Z,
                  pAtomCA2->xyz.X,pAtomCA2->xyz.Y,pAtomCA2->xyz.Z,
                  pAtomCA3->xyz.X,pAtomCA3->xyz.Y,pAtomCA3->xyz.Z,
                  cancbins[11][ii][jj][kk]*raddeg,cancbins[7][ii][jj][kk],cancbins[9][ii][jj][kk],
                  &pout.X,&pout.Y,&pout.Z);
  pAtomN3->xyz.X=pout.X;
  pAtomN3->xyz.Y=pout.Y;
  pAtomN3->xyz.Z=pout.Z;

  str2torp(pChainTemp,&chainFlexTemp,trandpos+1,trandpos+3);
  tor2strp2p2(pChainTemp,&chainFlexTemp,trandpos+1,trandpos+2);
  //flagTor=NewMainChainPos(StructureGetChain(&structureTemp,selectChain),&chainFlexTemp,trandpos);

  if(trandpos==0){   //should rebuild all chains in structure to fix problem with NewMainChainPos
    for(int j=trandpos; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  }
  else{
    for(int j=trandpos-1; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  }

  //TODO: optimize interface sidechains
  OptimizeSideChains(&structureTemp,bbrotlib,numSideChainOptCycles,selectChain); 
  EnergyTermInitialize(energyTerms);
  StructureCalcEnergyFlex(&structureTemp,bbrotlib,energyTerms);
  EnergyTermWeighting(energyTerms);
  newenergy=energyTerms[0];
//  newenergy+=calcConstraintEnergy(pChainTemp,natDistRestr)*75;
//  newenergy+=calcDistanceEnergy(pChainTemp,distanceCB,distanceWeight);
  if(newenergy<oldenergy){
    oldenergy=newenergy;
    //nummov[angtype][0]++;
    ChainFlex *origChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlexCopy(origChainFlex,&chainFlexTemp);
    Chain *origChain = StructureGetChain(decstr,selectChain);
    ChainCopy(origChain,StructureGetChain(&structureTemp,selectChain));
    getRamaOutliers(origChain,origChainFlex,ramaduke);
    getAngleOutliers(origChain,origChainFlex);
    getBondLengthOutliers(origChain,origChainFlex);
    getClashes(origChain,origChainFlex);
    StructureDestroy(&structureTemp);
    return TRUE;
  }
  else{
    StructureDestroy(&structureTemp);
    return FALSE;
  }
}

BOOL moveBackrub(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
                 bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,pairaa *natDistRestr,
                 double **contactCA,double **contactCB){
  const int CA_INDEX=3;
  int i,j;
  int trandpos,trandpos2;
  int numSideChainOptCycles=1;
  XYZ v1,v2,b1,b2;
  BOOL flagtor; 
  double raddeg = PI/180.0, angncac = 111.008;
  double trandang;
  double mati[9];
  double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
  
  int wpos;
  int selectChain = deschn[int(Random()*numchains)];

  ChainFlex* pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
  Chain* pChain = StructureGetChain(decstr,selectChain);
  int totwrong = pChainFlex->numRamaOutliers + pChainFlex->numAngleOutliers[2];
  if(totwrong>0 && Random()<0.1+totwrong/double(ChainFlexGetResidueCount(pChainFlex))){
    wpos=int(Random()*totwrong);
    i=-1;
    for(trandpos=0; trandpos<ChainFlexGetResidueCount(pChainFlex); trandpos++){
      ResidueFlex *pResidueFlex = ChainFlexGetResidue(pChainFlex,trandpos);
      if(pResidueFlex->isAngleOutlierAbove[2] || pResidueFlex->isAngleOutlierBelow[2]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isRamaOutlier) i++;
      if(i==wpos) break;
    }
    if(trandpos<1) trandpos=1;
    if(trandpos>ChainFlexGetResidueCount(pChainFlex)-2) trandpos=ChainFlexGetResidueCount(pChainFlex)-2;
  }
  else{
    do{     
      trandpos=1+int((ChainGetResidueCount(pChain)-2)*Random());
    }while(trandpos<1 || trandpos>ChainGetResidueCount(pChain)-2);
    double treject=Random();
    ResidueFlex *pResidueFlex1 = ChainFlexGetResidue(pChainFlex,trandpos-1);
    ResidueFlex *pResidueFlex2 = ChainFlexGetResidue(pChainFlex,trandpos);
    ResidueFlex *pResidueFlex3 = ChainFlexGetResidue(pChainFlex,trandpos+1);

    if((pResidueFlex2->ssm!='C' && pResidueFlex2->ssm==pResidueFlex1->ssm && 
        pResidueFlex2->ssm==pResidueFlex3->ssm && treject>0.1)){       
      //nummov[14][1]++;
      return FALSE;
    }
  }

  //StructureFlex structureFlexTemp;
  //for(i=0;i<StructureFlexGetChainCount(pStructureFlex);i++){
    //ChainFlex *pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlex chainFlexTemp;
    ChainFlexCreate(&chainFlexTemp);
    ChainFlexCopy(&chainFlexTemp,pChainFlex);
    //StructureFlexAddChain(structureFlexTemp, &chainFlex);
  //}

  Structure structureTemp;
  StructureCreate(&structureTemp);
  for(i=0;i<StructureGetChainCount(decstr);i++){
    Chain *pChainSelect = StructureGetChain(decstr,i);
    Chain chainTemp;
    ChainCreate(&chainTemp);
    ChainCopy(&chainTemp,pChainSelect);
    StructureAddChain(&structureTemp, &chainTemp);
    ChainDestroy(&chainTemp);
  }
  Chain *pChainTemp = StructureGetChain(&structureTemp,selectChain);
  Residue *pResidueA = ChainGetResidue(pChainTemp,trandpos-1);
  Residue *pResidueB = ChainGetResidue(pChainTemp,trandpos+1);
  Atom* pAtomNfirst = ResidueGetAtomByName(pResidueA, "N");
  Atom* pAtomCAfirst = ResidueGetAtomByName(pResidueA, "CA");
  Atom* pAtomCfirst = ResidueGetAtomByName(pResidueA, "C");
  Atom* pAtomNlast = ResidueGetAtomByName(pResidueB, "N");
  Atom* pAtomCAlast = ResidueGetAtomByName(pResidueB, "CA");
  Atom* pAtomClast = ResidueGetAtomByName(pResidueB, "C");

  trandang=180.0*Random()-90.0;
  trandang=trandang*raddeg;//rotangle
  v1.X=pAtomCAlast->xyz.X-pAtomCAfirst->xyz.X;
  v1.Y=pAtomCAlast->xyz.Y-pAtomCAfirst->xyz.Y;
  v1.Z=pAtomCAlast->xyz.Z-pAtomCAfirst->xyz.Z;
  v1=XYZUnit(&v1);//direction
  //method2
  double vcos=cos(trandang);
  double vsin=sin(trandang);
  double vxx=v1.X*v1.X;
  double vxy=v1.X*v1.Y;
  double vxz=v1.X*v1.Z;
  double vyy=v1.Y*v1.Y;
  double vyz=v1.Y*v1.Z;
  double vzz=v1.Z*v1.Z;
  mati[0]=vxx+(1-vxx)*vcos;
  mati[1]=vxy*(1-vcos)-v1.Z*vsin;
  mati[2]=vxz*(1-vcos)+v1.Y*vsin;
  mati[3]=vxy*(1-vcos)+v1.Z*vsin;
  mati[4]=vyy+(1-vyy)*vcos;
  mati[5]=vyz*(1-vcos)-v1.X*vsin;
  mati[6]=vxz*(1-vcos)-v1.Y*vsin;
  mati[7]=vyz*(1-vcos)+v1.X*vsin;
  mati[8]=vzz+(1-vzz)*vcos;

  for(i=trandpos;i<=trandpos;i++){
    Residue *pResiduei = ChainGetResidue(pChainTemp,i);
    Atom* pAtomNi = ResidueGetAtomByName(pResiduei, "N");
    Atom* pAtomCAi = ResidueGetAtomByName(pResiduei, "CA");
    Atom* pAtomCi = ResidueGetAtomByName(pResiduei, "C");
 
    v1.X=(pAtomCAi->xyz.X-pAtomCAfirst->xyz.X)*mati[0]+(pAtomCAi->xyz.Y-pAtomCAfirst->xyz.Y)*mati[1]
        +(pAtomCAi->xyz.Z-pAtomCAfirst->xyz.Z)*mati[2]+pAtomCAfirst->xyz.X;
    v1.Y=(pAtomCAi->xyz.X-pAtomCAfirst->xyz.X)*mati[3]+(pAtomCAi->xyz.Y-pAtomCAfirst->xyz.Y)*mati[4]
        +(pAtomCAi->xyz.Z-pAtomCAfirst->xyz.Z)*mati[5]+pAtomCAfirst->xyz.Y;
    v1.Z=(pAtomCAi->xyz.X-pAtomCAfirst->xyz.X)*mati[6]+(pAtomCAi->xyz.Y-pAtomCAfirst->xyz.Y)*mati[7]
        +(pAtomCAi->xyz.Z-pAtomCAfirst->xyz.Z)*mati[8]+pAtomCAfirst->xyz.Z;  
    pAtomCAi->xyz.X=v1.X;pAtomCAi->xyz.Y=v1.Y;pAtomCAi->xyz.Z=v1.Z;

    v1.X=(pAtomNi->xyz.X-pAtomCAfirst->xyz.X)*mati[0]+(pAtomNi->xyz.Y-pAtomCAfirst->xyz.Y)*mati[1]
        +(pAtomNi->xyz.Z-pAtomCAfirst->xyz.Z)*mati[2]+pAtomCAfirst->xyz.X;
    v1.Y=(pAtomNi->xyz.X-pAtomCAfirst->xyz.X)*mati[3]+(pAtomNi->xyz.Y-pAtomCAfirst->xyz.Y)*mati[4]
        +(pAtomNi->xyz.Z-pAtomCAfirst->xyz.Z)*mati[5]+pAtomCAfirst->xyz.Y;
    v1.Z=(pAtomNi->xyz.X-pAtomCAfirst->xyz.X)*mati[6]+(pAtomNi->xyz.Y-pAtomCAfirst->xyz.Y)*mati[7]
        +(pAtomNi->xyz.Z-pAtomCAfirst->xyz.Z)*mati[8]+pAtomCAfirst->xyz.Z;
    pAtomNi->xyz.X=v1.X;pAtomNi->xyz.Y=v1.Y;pAtomNi->xyz.Z=v1.Z;

    v1.X=(pAtomCi->xyz.X-pAtomCAfirst->xyz.X)*mati[0]+(pAtomCi->xyz.Y-pAtomCAfirst->xyz.Y)*mati[1]
        +(pAtomCi->xyz.Z-pAtomCAfirst->xyz.Z)*mati[2]+pAtomCAfirst->xyz.X;
    v1.Y=(pAtomCi->xyz.X-pAtomCAfirst->xyz.X)*mati[3]+(pAtomCi->xyz.Y-pAtomCAfirst->xyz.Y)*mati[4]
        +(pAtomCi->xyz.Z-pAtomCAfirst->xyz.Z)*mati[5]+pAtomCAfirst->xyz.Y;
    v1.Z=(pAtomCi->xyz.X-pAtomCAfirst->xyz.X)*mati[6]+(pAtomCi->xyz.Y-pAtomCAfirst->xyz.Y)*mati[7]
        +(pAtomCi->xyz.Z-pAtomCAfirst->xyz.Z)*mati[8]+pAtomCAfirst->xyz.Z;
    pAtomCi->xyz.X=v1.X;pAtomCi->xyz.Y=v1.Y;pAtomCi->xyz.Z=v1.Z;
  }

  v1.X=(pAtomCfirst->xyz.X-pAtomCAfirst->xyz.X)*mati[0]+(pAtomCfirst->xyz.Y-pAtomCAfirst->xyz.Y)*mati[1]
      +(pAtomCfirst->xyz.Z-pAtomCAfirst->xyz.Z)*mati[2]+pAtomCAfirst->xyz.X;
  v1.Y=(pAtomCfirst->xyz.X-pAtomCAfirst->xyz.X)*mati[3]+(pAtomCfirst->xyz.Y-pAtomCAfirst->xyz.Y)*mati[4]
      +(pAtomCfirst->xyz.Z-pAtomCAfirst->xyz.Z)*mati[5]+pAtomCAfirst->xyz.Y;
  v1.Z=(pAtomCfirst->xyz.X-pAtomCAfirst->xyz.X)*mati[6]+(pAtomCfirst->xyz.Y-pAtomCAfirst->xyz.Y)*mati[7]
      +(pAtomCfirst->xyz.Z-pAtomCAfirst->xyz.Z)*mati[8]+pAtomCAfirst->xyz.Z;
  pAtomCfirst->xyz.X=v1.X;pAtomCfirst->xyz.Y=v1.Y;pAtomCfirst->xyz.Z=v1.Z;

  v1.X=(pAtomNlast->xyz.X-pAtomCAfirst->xyz.X)*mati[0]+(pAtomNlast->xyz.Y-pAtomCAfirst->xyz.Y)*mati[1]
      +(pAtomNlast->xyz.Z-pAtomCAfirst->xyz.Z)*mati[2]+pAtomCAfirst->xyz.X;
  v1.Y=(pAtomNlast->xyz.X-pAtomCAfirst->xyz.X)*mati[3]+(pAtomNlast->xyz.Y-pAtomCAfirst->xyz.Y)*mati[4]
      +(pAtomNlast->xyz.Z-pAtomCAfirst->xyz.Z)*mati[5]+pAtomCAfirst->xyz.Y;
  v1.Z=(pAtomNlast->xyz.X-pAtomCAfirst->xyz.X)*mati[6]+(pAtomNlast->xyz.Y-pAtomCAfirst->xyz.Y)*mati[7]
      +(pAtomNlast->xyz.Z-pAtomCAfirst->xyz.Z)*mati[8]+pAtomCAfirst->xyz.Z;
  pAtomNlast->xyz.X=v1.X;pAtomNlast->xyz.Y=v1.Y;pAtomNlast->xyz.Z=v1.Z;
  
  //middle 2 rotate
  int nbin=72;
  double delta=360/double(nbin);
  int indbin;
  double mindist;
  double tdist;
  for(j=0;j<2;j++){
    mindist=100000000;
    for(i=0;i<nbin;i++){
      trandang=i*delta*raddeg;//rotangle
      Residue *pResidue1 = ChainGetResidue(pChainTemp,trandpos+j-1);
      Residue *pResidue2 = ChainGetResidue(pChainTemp,trandpos+j);
      Atom *pAtomCA1 =  ResidueGetAtomByName(pResidue1, "CA");
      Atom *pAtomC1 =  ResidueGetAtomByName(pResidue1, "C");
      Atom *pAtomN1 =  ResidueGetAtomByName(pResidue1, "N");
      Atom *pAtomCA2 =  ResidueGetAtomByName(pResidue2, "CA"); 
      Atom *pAtomC2 =  ResidueGetAtomByName(pResidue2, "C");
      Atom *pAtomN2 =  ResidueGetAtomByName(pResidue2, "N");

      v1=XYZDifference(&pAtomCA1->xyz,&pAtomCA2->xyz);
      v1=XYZUnit(&v1);//direction
      vcos=cos(trandang);
      vsin=sin(trandang);
      vxx=v1.X*v1.X;
      vxy=v1.X*v1.Y;
      vxz=v1.X*v1.Z;
      vyy=v1.Y*v1.Y;
      vyz=v1.Y*v1.Z;
      vzz=v1.Z*v1.Z;
      mati[0]=vxx+(1-vxx)*vcos;
      mati[1]=vxy*(1-vcos)-v1.Z*vsin;
      mati[2]=vxz*(1-vcos)+v1.Y*vsin;
      mati[3]=vxy*(1-vcos)+v1.Z*vsin;
      mati[4]=vyy+(1-vyy)*vcos;
      mati[5]=vyz*(1-vcos)-v1.X*vsin;
      mati[6]=vxz*(1-vcos)-v1.Y*vsin;
      mati[7]=vyz*(1-vcos)+v1.X*vsin;
      mati[8]=vzz+(1-vzz)*vcos;

      v1.X=(pAtomC1->xyz.X-pAtomCA1->xyz.X)*mati[0]+(pAtomC1->xyz.Y-pAtomCA1->xyz.Y)*mati[1]
          +(pAtomC1->xyz.Z-pAtomCA1->xyz.Z)*mati[2]+pAtomCA1->xyz.X;
      v1.Y=(pAtomC1->xyz.X-pAtomCA1->xyz.X)*mati[3]+(pAtomC1->xyz.Y-pAtomCA1->xyz.Y)*mati[4]
          +(pAtomC1->xyz.Z-pAtomCA1->xyz.Z)*mati[5]+pAtomCA1->xyz.Y;
      v1.Z=(pAtomC1->xyz.X-pAtomCA1->xyz.X)*mati[6]+(pAtomC1->xyz.Y-pAtomCA1->xyz.Y)*mati[7]
          +(pAtomC1->xyz.Z-pAtomCA1->xyz.Z)*mati[8]+pAtomCA1->xyz.Z;
      v2.X=(pAtomN2->xyz.X-pAtomCA1->xyz.X)*mati[0]+(pAtomN2->xyz.Y-pAtomCA1->xyz.Y)*mati[1]
          +(pAtomN2->xyz.Z-pAtomCA1->xyz.Z)*mati[2]+pAtomCA1->xyz.X;
      v2.Y=(pAtomN2->xyz.X-pAtomCA1->xyz.X)*mati[3]+(pAtomN2->xyz.Y-pAtomCA1->xyz.Y)*mati[4]
          +(pAtomN2->xyz.Z-pAtomCA1->xyz.Z)*mati[5]+pAtomCA1->xyz.Y;
      v2.Z=(pAtomN2->xyz.X-pAtomCA1->xyz.X)*mati[6]+(pAtomN2->xyz.Y-pAtomCA1->xyz.Y)*mati[7]
          +(pAtomN2->xyz.Z-pAtomCA1->xyz.Z)*mati[8]+pAtomCA1->xyz.Z;

      tdist=(v1.X-pAtomC1->xyz.X)*(v1.X-pAtomC1->xyz.X)+
            (v1.Y-pAtomC1->xyz.Y)*(v1.Y-pAtomC1->xyz.Y)+
            (v1.Z-pAtomC1->xyz.Z)*(v1.Z-pAtomC1->xyz.Z)+
            (v2.X-pAtomN2->xyz.X)*(v2.X-pAtomN2->xyz.X)+
            (v2.Y-pAtomN2->xyz.Y)*(v2.Y-pAtomN2->xyz.Y)+
            (v2.Z-pAtomN2->xyz.Z)*(v2.Z-pAtomN2->xyz.Z);

      if(tdist<mindist){
        b1.X=v1.X;b1.Y=v1.Y;b1.Z=v1.Z;
        b2.X=v2.X;b2.Y=v2.Y;b2.Z=v2.Z;
        mindist=tdist;
        indbin=i;
      }
    }
    Residue *pResidue1 = ChainGetResidue(pChainTemp,trandpos+j-1);
    Residue *pResidue2 = ChainGetResidue(pChainTemp,trandpos+j);
    Atom *pAtomC1 =  ResidueGetAtomByName(pResidue1, "C");
    Atom *pAtomN2 =  ResidueGetAtomByName(pResidue2, "N");
    pAtomC1->xyz=b1;
    pAtomN2->xyz=b2;
  }
  int istart=trandpos-1;
  int iend=trandpos+2;
  if(iend>=ChainGetResidueCount(pChainTemp)) iend=ChainGetResidueCount(pChainTemp)-1;
  str2torp(pChainTemp,&chainFlexTemp,istart,iend);

  ResidueFlex *pResidueFlex1 = ChainFlexGetResidue(&chainFlexTemp,trandpos-1);
  ResidueFlex *pResidueFlex2 = ChainFlexGetResidue(&chainFlexTemp,trandpos);
  ResidueFlex *pResidueFlex3 = ChainFlexGetResidue(&chainFlexTemp,trandpos+1);
 
  if(fabs(pResidueFlex1->ang[2]-angncac)>10.0 || fabs(pResidueFlex2->ang[2]-angncac)>10.0 || 
     fabs(pResidueFlex3->ang[2]-angncac)>10.0){
    //nummov[14][1]++;
    StructureDestroy(&structureTemp);
    return FALSE;
  }
  if(istart==0 || istart==1) flagtor=NewMainChainPos(pChainTemp,&chainFlexTemp,0);
  else flagtor=NewMainChainPos(pChainTemp,&chainFlexTemp,istart);
  if(!flagtor){
    printf("tor2str wrong in rub %d\n",istart);
  }

  if(trandpos==0){   //should rebuild all chains in structure to fix problem with NewMainChainPos
    for(int j=trandpos; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  }
  else{
    for(int j=trandpos-1; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  }

  //TODO: optimize interface sidechains
  OptimizeSideChains(&structureTemp,bbrotlib,numSideChainOptCycles,selectChain); 
  EnergyTermInitialize(energyTerms);
  StructureCalcEnergyFlex(&structureTemp,bbrotlib,energyTerms);
  EnergyTermWeighting(energyTerms);
  newenergy=energyTerms[0];
  //std::cout << "Old energy: " << oldenergy << " New energy: " << newenergy << std::endl;
//  newenergy+=calcConstraintEnergy(pChainTemp,natDistRestr)*75;
//  newenergy+=calcDistanceEnergy(pChainTemp,distanceCB,distanceWeight);
  if(newenergy<oldenergy){
    oldenergy=newenergy;
    //nummov[angtype][0]++;
    ChainFlex *origChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlexCopy(origChainFlex,&chainFlexTemp);
    Chain *origChain = StructureGetChain(decstr,selectChain);
    ChainCopy(origChain,StructureGetChain(&structureTemp,selectChain));
    getRamaOutliers(origChain,origChainFlex,ramaduke);
    getAngleOutliers(origChain,origChainFlex);
    getBondLengthOutliers(origChain,origChainFlex);
    getClashes(origChain,origChainFlex);
    StructureDestroy(&structureTemp);
    return TRUE;
  }
  else{
    StructureDestroy(&structureTemp);
    return FALSE;
  }
}

BOOL moveRotation(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
                  bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,pairaa *natDistRestr,
                  double **contactCA,double **contactCB){
  const int CA_INDEX=3;
  int i,j;
  int trandpos,trandpos2;
  int numSideChainOptCycles=1;
  XYZ v1;
  BOOL flagTor; 
  double raddeg = PI/180.0;
  double trandang;
  double mati[9];
  double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
  
  int wpos;
  int selectChain = deschn[int(Random()*numchains)];

  ChainFlex* pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
  Chain* pChain = StructureGetChain(decstr,selectChain);
  int totwrong = pChainFlex->numRamaOutliers + pChainFlex->numAngleOutliers[2];
  if(totwrong>0 && Random()<0.1+totwrong/double(ChainFlexGetResidueCount(pChainFlex))){
    wpos=int(Random()*totwrong);
    i=-1;
    for(trandpos=0; trandpos<ChainFlexGetResidueCount(pChainFlex); trandpos++){
      ResidueFlex *pResidueFlex = ChainFlexGetResidue(pChainFlex,trandpos);
      if(pResidueFlex->isAngleOutlierAbove[2] || pResidueFlex->isAngleOutlierBelow[2]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isRamaOutlier) i++;
      if(i==wpos) break;
    }
    if(Random()<0.5)//leftend
    {       
      trandpos2=int(20*Random())-1;
      trandpos2+=trandpos;
    }
    else//rightend
    {       
      trandpos2=trandpos;
      trandpos=int(20*Random())-1;
      trandpos=trandpos2-trandpos;
    }
    if(trandpos<2) trandpos=2;
    if(trandpos2>(ChainGetResidueCount(pChain)-2)) trandpos2=ChainGetResidueCount(pChain)-2;
  }
  //  else if(lp1.nn[12]>0 && Random()<0.7)//clash
  //  {
  //    wpos=int(Random()*lp1.nn[12]);
  //    i=-1;
  //    for(trandpos=0;trandpos<numseq;trandpos++){
  //      if(lp1.indn[12][trandpos]==1) i++;
  //      if(i==wpos) break;
  //    }
  //    trandpos2=int(20*Random())-1;
  //    trandpos-=trandpos2/2;
  //    trandpos2+=trandpos;
  //    if(trandpos<2) trandpos=2;
  //    if(trandpos2>numseq-2) trandpos2=numseq-2;
  //  }
  else{
    //[trandpos, trandpos2] [2 numseq-2] //1 will affect 0 c
    do{
      trandpos2=int(20*Random())-1;
      trandpos=2+int((ChainGetResidueCount(pChain)-trandpos2-3)*Random());
      trandpos2+=trandpos;
    }while(trandpos<2 || trandpos2>ChainGetResidueCount(pChain)-2);
  }
  if(trandpos2<trandpos){
    //nummov[10][1]++;
    return FALSE;
  }

  //StructureFlex structureFlexTemp;
  //for(i=0;i<StructureFlexGetChainCount(pStructureFlex);i++){
    //ChainFlex *pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlex chainFlexTemp;
    ChainFlexCreate(&chainFlexTemp);
    ChainFlexCopy(&chainFlexTemp,pChainFlex);
    //StructureFlexAddChain(structureFlexTemp, &chainFlex);
  //}

  Structure structureTemp;
  StructureCreate(&structureTemp);
  for(i=0;i<StructureGetChainCount(decstr);i++){
    Chain *pChainSelect = StructureGetChain(decstr,i);
    Chain chainTemp;
    ChainCreate(&chainTemp);
    ChainCopy(&chainTemp,pChainSelect);
    StructureAddChain(&structureTemp, &chainTemp);
    ChainDestroy(&chainTemp);
  }
  Chain *pChainTemp = StructureGetChain(&structureTemp,selectChain);
  Residue *pResidueA = ChainGetResidue(pChainTemp,trandpos-1);
  Residue *pResidueB = ChainGetResidue(pChainTemp,trandpos2+1);
  Atom* pAtomNfirst = ResidueGetAtomByName(pResidueA, "N");
  Atom* pAtomCAfirst = ResidueGetAtomByName(pResidueA, "CA");
  Atom* pAtomCfirst = ResidueGetAtomByName(pResidueA, "C");
  Atom* pAtomNlast = ResidueGetAtomByName(pResidueB, "N");
  Atom* pAtomCAlast = ResidueGetAtomByName(pResidueB, "CA");
  Atom* pAtomClast = ResidueGetAtomByName(pResidueB, "C");

  trandang=180.0*Random()-90.0;
  trandang=trandang*raddeg;//rotangle
  v1.X=pAtomCAlast->xyz.X-pAtomCAfirst->xyz.X;
  v1.Y=pAtomCAlast->xyz.Y-pAtomCAfirst->xyz.Y;
  v1.Z=pAtomCAlast->xyz.Z-pAtomCAfirst->xyz.Z;
  v1=XYZUnit(&v1);//direction
  //method2
  double vcos=cos(trandang);
  double vsin=sin(trandang);
  double vxx=v1.X*v1.X;
  double vxy=v1.X*v1.Y;
  double vxz=v1.X*v1.Z;
  double vyy=v1.Y*v1.Y;
  double vyz=v1.Y*v1.Z;
  double vzz=v1.Z*v1.Z;
  mati[0]=vxx+(1-vxx)*vcos;
  mati[1]=vxy*(1-vcos)-v1.Z*vsin;
  mati[2]=vxz*(1-vcos)+v1.Y*vsin;
  mati[3]=vxy*(1-vcos)+v1.Z*vsin;
  mati[4]=vyy+(1-vyy)*vcos;
  mati[5]=vyz*(1-vcos)-v1.X*vsin;
  mati[6]=vxz*(1-vcos)-v1.Y*vsin;
  mati[7]=vyz*(1-vcos)+v1.X*vsin;
  mati[8]=vzz+(1-vzz)*vcos;

  for(i=trandpos;i<=trandpos2;i++){
    Residue *pResiduei = ChainGetResidue(pChainTemp,i);
    Atom* pAtomNi = ResidueGetAtomByName(pResiduei, "N");
    Atom* pAtomCAi = ResidueGetAtomByName(pResiduei, "CA");
    Atom* pAtomCi = ResidueGetAtomByName(pResiduei, "C");
 
    v1.X=(pAtomCAi->xyz.X-pAtomCAfirst->xyz.X)*mati[0]+(pAtomCAi->xyz.Y-pAtomCAfirst->xyz.Y)*mati[1]
        +(pAtomCAi->xyz.Z-pAtomCAfirst->xyz.Z)*mati[2]+pAtomCAfirst->xyz.X;
    v1.Y=(pAtomCAi->xyz.X-pAtomCAfirst->xyz.X)*mati[3]+(pAtomCAi->xyz.Y-pAtomCAfirst->xyz.Y)*mati[4]
        +(pAtomCAi->xyz.Z-pAtomCAfirst->xyz.Z)*mati[5]+pAtomCAfirst->xyz.Y;
    v1.Z=(pAtomCAi->xyz.X-pAtomCAfirst->xyz.X)*mati[6]+(pAtomCAi->xyz.Y-pAtomCAfirst->xyz.Y)*mati[7]
        +(pAtomCAi->xyz.Z-pAtomCAfirst->xyz.Z)*mati[8]+pAtomCAfirst->xyz.Z;  
    pAtomCAi->xyz.X=v1.X;pAtomCAi->xyz.Y=v1.Y;pAtomCAi->xyz.Z=v1.Z;

    v1.X=(pAtomNi->xyz.X-pAtomCAfirst->xyz.X)*mati[0]+(pAtomNi->xyz.Y-pAtomCAfirst->xyz.Y)*mati[1]
        +(pAtomNi->xyz.Z-pAtomCAfirst->xyz.Z)*mati[2]+pAtomCAfirst->xyz.X;
    v1.Y=(pAtomNi->xyz.X-pAtomCAfirst->xyz.X)*mati[3]+(pAtomNi->xyz.Y-pAtomCAfirst->xyz.Y)*mati[4]
        +(pAtomNi->xyz.Z-pAtomCAfirst->xyz.Z)*mati[5]+pAtomCAfirst->xyz.Y;
    v1.Z=(pAtomNi->xyz.X-pAtomCAfirst->xyz.X)*mati[6]+(pAtomNi->xyz.Y-pAtomCAfirst->xyz.Y)*mati[7]
        +(pAtomNi->xyz.Z-pAtomCAfirst->xyz.Z)*mati[8]+pAtomCAfirst->xyz.Z;
    pAtomNi->xyz.X=v1.X;pAtomNi->xyz.Y=v1.Y;pAtomNi->xyz.Z=v1.Z;

    v1.X=(pAtomCi->xyz.X-pAtomCAfirst->xyz.X)*mati[0]+(pAtomCi->xyz.Y-pAtomCAfirst->xyz.Y)*mati[1]
        +(pAtomCi->xyz.Z-pAtomCAfirst->xyz.Z)*mati[2]+pAtomCAfirst->xyz.X;
    v1.Y=(pAtomCi->xyz.X-pAtomCAfirst->xyz.X)*mati[3]+(pAtomCi->xyz.Y-pAtomCAfirst->xyz.Y)*mati[4]
        +(pAtomCi->xyz.Z-pAtomCAfirst->xyz.Z)*mati[5]+pAtomCAfirst->xyz.Y;
    v1.Z=(pAtomCi->xyz.X-pAtomCAfirst->xyz.X)*mati[6]+(pAtomCi->xyz.Y-pAtomCAfirst->xyz.Y)*mati[7]
        +(pAtomCi->xyz.Z-pAtomCAfirst->xyz.Z)*mati[8]+pAtomCAfirst->xyz.Z;
    pAtomCi->xyz.X=v1.X;pAtomCi->xyz.Y=v1.Y;pAtomCi->xyz.Z=v1.Z;
  }

  v1.X=(pAtomCfirst->xyz.X-pAtomCAfirst->xyz.X)*mati[0]+(pAtomCfirst->xyz.Y-pAtomCAfirst->xyz.Y)*mati[1]
      +(pAtomCfirst->xyz.Z-pAtomCAfirst->xyz.Z)*mati[2]+pAtomCAfirst->xyz.X;
  v1.Y=(pAtomCfirst->xyz.X-pAtomCAfirst->xyz.X)*mati[3]+(pAtomCfirst->xyz.Y-pAtomCAfirst->xyz.Y)*mati[4]
      +(pAtomCfirst->xyz.Z-pAtomCAfirst->xyz.Z)*mati[5]+pAtomCAfirst->xyz.Y;
  v1.Z=(pAtomCfirst->xyz.X-pAtomCAfirst->xyz.X)*mati[6]+(pAtomCfirst->xyz.Y-pAtomCAfirst->xyz.Y)*mati[7]
      +(pAtomCfirst->xyz.Z-pAtomCAfirst->xyz.Z)*mati[8]+pAtomCAfirst->xyz.Z;
  pAtomCfirst->xyz.X=v1.X;pAtomCfirst->xyz.Y=v1.Y;pAtomCfirst->xyz.Z=v1.Z;

  v1.X=(pAtomNlast->xyz.X-pAtomCAfirst->xyz.X)*mati[0]+(pAtomNlast->xyz.Y-pAtomCAfirst->xyz.Y)*mati[1]
      +(pAtomNlast->xyz.Z-pAtomCAfirst->xyz.Z)*mati[2]+pAtomCAfirst->xyz.X;
  v1.Y=(pAtomNlast->xyz.X-pAtomCAfirst->xyz.X)*mati[3]+(pAtomNlast->xyz.Y-pAtomCAfirst->xyz.Y)*mati[4]
      +(pAtomNlast->xyz.Z-pAtomCAfirst->xyz.Z)*mati[5]+pAtomCAfirst->xyz.Y;
  v1.Z=(pAtomNlast->xyz.X-pAtomCAfirst->xyz.X)*mati[6]+(pAtomNlast->xyz.Y-pAtomCAfirst->xyz.Y)*mati[7]
      +(pAtomNlast->xyz.Z-pAtomCAfirst->xyz.Z)*mati[8]+pAtomCAfirst->xyz.Z;
  pAtomNlast->xyz.X=v1.X;pAtomNlast->xyz.Y=v1.Y;pAtomNlast->xyz.Z=v1.Z;

  int istart=trandpos-1;
  int iend=trandpos2+2;
  if(iend>=ChainGetResidueCount(pChainTemp)) iend=ChainGetResidueCount(pChainTemp)-1;
  str2torp(pChainTemp,&chainFlexTemp,istart,iend);
  flagTor=tor2strp2p2(pChainTemp,&chainFlexTemp,istart,iend);

  if(trandpos==0){   //should rebuild all chains in structure to fix problem with NewMainChainPos
    for(int j=trandpos; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  }
  else{
    for(int j=trandpos-1; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  }

  //TODO: optimize interface sidechains
  OptimizeSideChains(&structureTemp,bbrotlib,numSideChainOptCycles,selectChain); 
  EnergyTermInitialize(energyTerms);
  StructureCalcEnergyFlex(&structureTemp,bbrotlib,energyTerms);
  EnergyTermWeighting(energyTerms);
  newenergy=energyTerms[0];
  //std::cout << "Old energy: " << oldenergy << " New energy: " << newenergy << std::endl;
//  newenergy+=calcConstraintEnergy(pChainTemp,natDistRestr)*75;
//  newenergy+=calcDistanceEnergy(pChainTemp,distanceCB,distanceWeight);
  if(newenergy<oldenergy){
    oldenergy=newenergy;
    //nummov[angtype][0]++;
    ChainFlex *origChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlexCopy(origChainFlex,&chainFlexTemp);
    Chain *origChain = StructureGetChain(decstr,selectChain);
    ChainCopy(origChain,StructureGetChain(&structureTemp,selectChain));
    getRamaOutliers(origChain,origChainFlex,ramaduke);
    getAngleOutliers(origChain,origChainFlex);
    getBondLengthOutliers(origChain,origChainFlex);
    getClashes(origChain,origChainFlex);
    StructureDestroy(&structureTemp);
    return TRUE;
  }
  else{
    StructureDestroy(&structureTemp);
    return FALSE;
  }
}


void singlemoveLMP(Chain *pChain,int k,int mtype,int flagdir){
  XYZ rp;
  double rnorm;
  double lencn = 1.338f,delcn = 0.005f*15,lencan1 = 2.441f,delcan1 = 0.036f*4;
  double lennca = 1.460f,delnca = 0.004f*20,lencca1 = 2.446f,delcca1 = 0.036f*4;
  double lencaca = 3.813f,delcaca = 0.019f*10,lencac = 1.525f,delcac = 0.004f*20;
  double lennc = 2.460f,delnc = 0.012f*10;

  Residue *pResiPrev = ChainGetResidue(pChain,k-1);
  Residue *pResi = ChainGetResidue(pChain,k);
  Atom *pAtomNprev = ResidueGetAtomByName(pResiPrev, "N"); 
  Atom *pAtomCAprev = ResidueGetAtomByName(pResiPrev, "CA"); 
  Atom *pAtomCprev = ResidueGetAtomByName(pResiPrev, "C"); 
  Atom *pAtomN = ResidueGetAtomByName(pResi, "N"); 
  Atom *pAtomCA = ResidueGetAtomByName(pResi, "CA"); 
  Atom *pAtomC = ResidueGetAtomByName(pResi, "C"); 

  if(mtype==0 && flagdir==0){
    //1 c-1 n
    rp=XYZDifference(&pAtomCprev->xyz,&pAtomN->xyz);
    rnorm=XYZNormalization(&rp);
    if(rnorm<err)//two in the same place
    {
      pAtomN->xyz.X+=0.0001f;pAtomN->xyz.Y+=0.0001f;pAtomN->xyz.Z+=0.0001f;
      rp=XYZDifference(&pAtomCprev->xyz,&pAtomN->xyz);
      rnorm=XYZNormalization(&rp);
    }
    if(rnorm<lencn-delcn || rnorm>lencn+delcn){
      rp=XYZScale2(&rp,lencn/rnorm);
      pAtomN->xyz.X=pAtomCprev->xyz.X+rp.X;
      pAtomN->xyz.Y=pAtomCprev->xyz.Y+rp.Y;
      pAtomN->xyz.Z=pAtomCprev->xyz.Z+rp.Z;
    }
  }
  else if(mtype==1 && flagdir==0){
    //2 ca-1 n
    rp=XYZDifference(&pAtomCAprev->xyz,&pAtomN->xyz);
    rnorm=XYZNormalization(&rp);
    if(rnorm<err)//two in the same place
    {
      pAtomN->xyz.X+=0.0001f;pAtomN->xyz.Y+=0.0001f;pAtomN->xyz.Z+=0.0001f;
      rp=XYZDifference(&pAtomCAprev->xyz,&pAtomN->xyz);
      rnorm=XYZNormalization(&rp);
    }
    if(rnorm<lencan1-delcan1 || rnorm>lencan1+delcan1){
      rp=XYZScale2(&rp,lencan1/rnorm);
      pAtomN->xyz.X=pAtomCAprev->xyz.X+rp.X;
      pAtomN->xyz.Y=pAtomCAprev->xyz.Y+rp.Y;
      pAtomN->xyz.Z=pAtomCAprev->xyz.Z+rp.Z;
    }
  }
  else if(mtype==2 && flagdir==0){
    //3 n ca
    rp=XYZDifference(&pAtomN->xyz,&pAtomCA->xyz);
    rnorm=XYZNormalization(&rp);
    if(rnorm<err)//two in the same place
    {
      pAtomCA->xyz.X+=0.0001f;pAtomCA->xyz.Y+=0.0001f;pAtomCA->xyz.Z+=0.0001f;
      rp=XYZDifference(&pAtomN->xyz,&pAtomCA->xyz);
      rnorm=XYZNormalization(&rp);
    }
    if(rnorm<lennca-delnca || rnorm>lennca+delnca){
      rp=XYZScale2(&rp,lennca/rnorm);
      pAtomCA->xyz.X=pAtomN->xyz.X+rp.X;
      pAtomCA->xyz.Y=pAtomN->xyz.Y+rp.Y;
      pAtomCA->xyz.Z=pAtomN->xyz.Z+rp.Z;
    }
  }
  else if(mtype==3 && flagdir==0){
    //4 c-1 ca
    rp=XYZDifference(&pAtomCprev->xyz,&pAtomCA->xyz);
    rnorm=XYZNormalization(&rp);
    if(rnorm<err)//two in the same place
    {
      pAtomCA->xyz.X+=0.0001f;pAtomCA->xyz.Y+=0.0001f;pAtomCA->xyz.Z+=0.0001f;
      rp=XYZDifference(&pAtomCprev->xyz,&pAtomCA->xyz);
      rnorm=XYZNormalization(&rp);
    }
    if(rnorm<lencca1-delcca1 || rnorm>lencca1+delcca1){
      rp=XYZScale2(&rp,lencca1/rnorm);
      pAtomCA->xyz.X=pAtomCprev->xyz.X+rp.X;
      pAtomCA->xyz.Y=pAtomCprev->xyz.Y+rp.Y;
      pAtomCA->xyz.Z=pAtomCprev->xyz.Z+rp.Z;
    }
  }
  else if(mtype==4 && flagdir==0){
    //5 ca-1 ca                               
    rp=XYZDifference(&pAtomCAprev->xyz,&pAtomCA->xyz);
    rnorm=XYZNormalization(&rp);
    if(rnorm<err)//two in the same place
    {
      pAtomCA->xyz.X+=0.0001f;pAtomCA->xyz.Y+=0.0001f;pAtomCA->xyz.Z+=0.0001f;
      rp=XYZDifference(&pAtomCAprev->xyz,&pAtomCA->xyz);
      rnorm=XYZNormalization(&rp);      
    }
    if(rnorm<lencaca-delcaca || rnorm>lencaca+delcaca){
      rp=XYZScale2(&rp,lencaca/rnorm);
      pAtomCA->xyz.X=pAtomCAprev->xyz.X+rp.X;
      pAtomCA->xyz.Y=pAtomCAprev->xyz.Y+rp.Y;
      pAtomCA->xyz.Z=pAtomCAprev->xyz.Z+rp.Z;   
    }
  }
  else if(mtype==5 && flagdir==0){
    //6 ca c
    rp=XYZDifference(&pAtomCA->xyz,&pAtomC->xyz);
    rnorm=XYZNormalization(&rp);
    if(rnorm<err)//two in the same place
    {
      pAtomC->xyz.X+=0.0001f;pAtomC->xyz.Y+=0.0001f;pAtomC->xyz.Z+=0.0001f;
      rp=XYZDifference(&pAtomCA->xyz,&pAtomC->xyz);
      rnorm=XYZNormalization(&rp);
    }
    if(rnorm<lencac-delcac || rnorm>lencac+delcac){
      rp=XYZScale2(&rp,lencac/rnorm);
      pAtomC->xyz.X=pAtomCA->xyz.X+rp.X;
      pAtomC->xyz.Y=pAtomCA->xyz.Y+rp.Y;
      pAtomC->xyz.Z=pAtomCA->xyz.Z+rp.Z;
    }
  }
  else if(mtype==6 && flagdir==0){
    //7 n c
    rp=XYZDifference(&pAtomN->xyz,&pAtomC->xyz);
    rnorm=XYZNormalization(&rp);
    if(rnorm<err)//two in the same place
    {
      pAtomC->xyz.X+=0.0001f;pAtomC->xyz.Y+=0.0001f;pAtomC->xyz.Z+=0.0001f;
      rp=XYZDifference(&pAtomN->xyz,&pAtomC->xyz);
      rnorm=XYZNormalization(&rp);      
    }
    if(rnorm<lennc-delnc || rnorm>lennc+delnc){
      rp=XYZScale2(&rp,lennc/rnorm);    
      pAtomC->xyz.X=pAtomN->xyz.X+rp.X;
      pAtomC->xyz.Y=pAtomN->xyz.Y+rp.Y;
      pAtomC->xyz.Z=pAtomN->xyz.Z+rp.Z;       
    }
  }
  else if(mtype==0 && flagdir==1){
    //1 c-1 n
    rp=XYZDifference(&pAtomN->xyz,&pAtomCprev->xyz);
    rnorm=XYZNormalization(&rp);
    if(rnorm<err)//two in the same place
    {
      pAtomCprev->xyz.X-=0.0001f;pAtomCprev->xyz.Y-=0.0001f;pAtomCprev->xyz.Z-=0.0001f;
      rp=XYZDifference(&pAtomN->xyz,&pAtomCprev->xyz);
      rnorm=XYZNormalization(&rp);      
    }
    if(rnorm<lencn-delcn || rnorm>lencn+delcn){
      rp=XYZScale2(&rp,lencn/rnorm);
      pAtomCprev->xyz.X=pAtomN->xyz.X+rp.X;
      pAtomCprev->xyz.Y=pAtomN->xyz.Y+rp.Y; 
      pAtomCprev->xyz.Z=pAtomN->xyz.Z+rp.Z;        
    }
  }
  else if(mtype==1 && flagdir==1){
    //4 c-1 ca
    rp=XYZDifference(&pAtomCA->xyz,&pAtomCprev->xyz);
    rnorm=XYZNormalization(&rp);
    if(rnorm<err)//two in the same place
    {
      pAtomCprev->xyz.X-=0.0001f;pAtomCprev->xyz.Y-=0.0001f;pAtomCprev->xyz.Z-=0.0001f;
      rp=XYZDifference(&pAtomCA->xyz,&pAtomCprev->xyz);
      rnorm=XYZNormalization(&rp);      
    }
    if(rnorm<lencca1-delcca1 || rnorm>lencca1+delcca1){
      rp=XYZScale2(&rp,lencca1/rnorm);
      pAtomCprev->xyz.X=pAtomCA->xyz.X+rp.X;
      pAtomCprev->xyz.Y=pAtomCA->xyz.Y+rp.Y;
      pAtomCprev->xyz.Z=pAtomCA->xyz.Z+rp.Z;        
    }
  }
  else if(mtype==2 && flagdir==1){
    //6 ca-1 c-1
    rp=XYZDifference(&pAtomCprev->xyz,&pAtomCAprev->xyz);
    rnorm=XYZNormalization(&rp);
    if(rnorm<err)//two in the same place
    {
      pAtomCAprev->xyz.X-=0.0001f;pAtomCAprev->xyz.Y-=0.0001f;pAtomCAprev->xyz.Z-=0.0001f;
      rp=XYZDifference(&pAtomCprev->xyz,&pAtomCAprev->xyz);
      rnorm=XYZNormalization(&rp);      
    }
    if(rnorm<lencac-delcac || rnorm>lencac+delcac){
      rp=XYZScale2(&rp,lencac/rnorm);
      pAtomCAprev->xyz.X=pAtomCprev->xyz.X+rp.X;
      pAtomCAprev->xyz.Y=pAtomCprev->xyz.Y+rp.Y;
      pAtomCAprev->xyz.Z=pAtomCprev->xyz.Z+rp.Z;
    }
  }
  else if(mtype==3 && flagdir==1){
    //2 ca-1 n
    rp=XYZDifference(&pAtomN->xyz,&pAtomCAprev->xyz);
    rnorm=XYZNormalization(&rp);
    if(rnorm<err)//two in the same place
    {
      pAtomCAprev->xyz.X-=0.0001f;pAtomCAprev->xyz.Y-=0.0001f;pAtomCAprev->xyz.Z-=0.0001f;
      rp=XYZDifference(&pAtomN->xyz,&pAtomCAprev->xyz);
      rnorm=XYZNormalization(&rp);      
    }
    if(rnorm<lencan1-delcan1 || rnorm>lencan1+delcan1){
      rp=XYZScale2(&rp,lencan1/rnorm);
      pAtomCAprev->xyz.X=pAtomN->xyz.X+rp.X;
      pAtomCAprev->xyz.Y=pAtomN->xyz.Y+rp.Y;
      pAtomCAprev->xyz.Z=pAtomN->xyz.Z+rp.Z;
    }
  }
  else if(mtype==4 && flagdir==1){
    //5 ca-1 ca
    rp=XYZDifference(&pAtomCA->xyz,&pAtomCAprev->xyz);
    rnorm=XYZNormalization(&rp);
    if(rnorm<err)//two in the same place
    {
      pAtomCAprev->xyz.X-=0.0001f;pAtomCAprev->xyz.Y-=0.0001f;pAtomCAprev->xyz.Z-=0.0001f;
      rp=XYZDifference(&pAtomCA->xyz,&pAtomCAprev->xyz);
      rnorm=XYZNormalization(&rp);      
    }
    if(rnorm<lencaca-delcaca || rnorm>lencaca+delcaca){
      rp=XYZScale2(&rp,lencaca/rnorm);
      pAtomCAprev->xyz.X=pAtomCA->xyz.X+rp.X;
      pAtomCAprev->xyz.Y=pAtomCA->xyz.Y+rp.Y;
      pAtomCAprev->xyz.Z=pAtomCA->xyz.Z+rp.Z;   
    }
  }
  else if(mtype==5 && flagdir==1){
    //3 n-1 ca-1
    rp=XYZDifference(&pAtomCAprev->xyz,&pAtomNprev->xyz);
    rnorm=XYZNormalization(&rp);
    if(rnorm<err)//two in the same place
    {
      pAtomNprev->xyz.X-=0.0001f;pAtomNprev->xyz.Y-=0.0001f;pAtomNprev->xyz.Z-=0.0001f;
      rp=XYZDifference(&pAtomCAprev->xyz,&pAtomNprev->xyz);
      rnorm=XYZNormalization(&rp);      
    }
    if(rnorm<lennca-delnca || rnorm>lennca+delnca){
      rp=XYZScale2(&rp,lennca/rnorm);
      pAtomNprev->xyz.X=pAtomCAprev->xyz.X+rp.X;
      pAtomNprev->xyz.Y=pAtomCAprev->xyz.Y+rp.Y;
      pAtomNprev->xyz.Z=pAtomCAprev->xyz.Z+rp.Z;
    }
  }
  else if(mtype==6 && flagdir==1){
    //7 n c
    rp=XYZDifference(&pAtomCprev->xyz,&pAtomNprev->xyz);
    rnorm=XYZNormalization(&rp);
    if(rnorm<err)//two in the same place
    {
      pAtomNprev->xyz.X-=0.0001f;pAtomNprev->xyz.Y-=0.0001f;pAtomNprev->xyz.Z-=0.0001f;
      rp=XYZDifference(&pAtomCprev->xyz,&pAtomNprev->xyz);
      rnorm=XYZNormalization(&rp);      
    }
    if(rnorm<lennc-delnc || rnorm>lennc+delnc){
      rp=XYZScale2(&rp,lennc/rnorm);
      pAtomNprev->xyz.X=pAtomCprev->xyz.X+rp.X;
      pAtomNprev->xyz.Y=pAtomCprev->xyz.Y+rp.Y;
      pAtomNprev->xyz.Z=pAtomCprev->xyz.Z+rp.Z;
    }
  }
}

BOOL moveLMP(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
             bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,pairaa *natDistRestr,
             double **contactCA,double **contactCB){
  const int CA_INDEX=3;
  XYZ rp;
  int i,j;
  int trandpos;
  int numSideChainOptCycles=1;
  BOOL flagTor,flagdone; 
  double trand2;
  double rnorm,rmax;
  double lencn = 1.338f,delcn = 0.005f*15,lencan1 = 2.441f,delcan1 = 0.036f*4;
  double lennca = 1.460f,delnca = 0.004f*20,lencca1 = 2.446f,delcca1 = 0.036f*4;
  double lencaca = 3.813f,delcaca = 0.019f*10,lencac = 1.525f,delcac = 0.004f*20;
  double lennc = 2.460f,delnc = 0.012f*10;
  double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
  
  int wpos;
  int selectChain = deschn[int(Random()*numchains)];

  ChainFlex* pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
  Chain* pChain = StructureGetChain(decstr,selectChain);
  int totwrong = pChainFlex->numRamaOutliers + pChainFlex->numAngleOutliers[0] + pChainFlex->numAngleOutliers[1]
               + pChainFlex->numAngleOutliers[2] + pChainFlex->numBondLengthOutliers[0] + pChainFlex->numBondLengthOutliers[1]
               + pChainFlex->numBondLengthOutliers[2] + pChainFlex->numBondLengthOutliers[3] + pChainFlex->numRamaOutliers;

  if(totwrong>0 && Random()<0.1+totwrong/double(ChainFlexGetResidueCount(pChainFlex))){
    wpos=int(Random()*totwrong);
    i=-1;
    for(trandpos=0; trandpos<ChainFlexGetResidueCount(pChainFlex); trandpos++){
      ResidueFlex *pResidueFlex = ChainFlexGetResidue(pChainFlex,trandpos);

      if(pResidueFlex->isAngleOutlierAbove[0] || pResidueFlex->isAngleOutlierBelow[0]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isAngleOutlierAbove[1] || pResidueFlex->isAngleOutlierBelow[1]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isAngleOutlierAbove[2] || pResidueFlex->isAngleOutlierBelow[2]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isBondLengthOutlierAbove[0] || pResidueFlex->isBondLengthOutlierBelow[0]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isBondLengthOutlierAbove[1] || pResidueFlex->isBondLengthOutlierBelow[1]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isBondLengthOutlierAbove[2] || pResidueFlex->isBondLengthOutlierBelow[2]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isBondLengthOutlierAbove[3] || pResidueFlex->isBondLengthOutlierBelow[3]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isRamaOutlier) i++;
      if(i==wpos) break;
    }
    trand2=2+int(6*Random());
    trandpos-=trand2/2;
    trand2+=trandpos;
    if(trandpos<1) trandpos=1;
    if(trand2>ChainGetResidueCount(pChain)-2) trand2=ChainGetResidueCount(pChain)-2;
    trand2-=trandpos;
    if(trand2<2){
      //nummov[11][1]++;
      return FALSE;
    }
  }
  else{
    trand2=2+int(6*Random());
    trandpos=1+int((ChainGetResidueCount(pChain)-trand2-1)*Random());
  }

  //StructureFlex structureFlexTemp;
  //for(i=0;i<StructureFlexGetChainCount(pStructureFlex);i++){
    //ChainFlex *pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlex chainFlexTemp;
    ChainFlexCreate(&chainFlexTemp);
    ChainFlexCopy(&chainFlexTemp,pChainFlex);
    //StructureFlexAddChain(structureFlexTemp, &chainFlex);
  //}

  Structure structureTemp;
  StructureCreate(&structureTemp);
  for(i=0;i<StructureGetChainCount(decstr);i++){
    Chain *pChainSelect = StructureGetChain(decstr,i);
    Chain chainTemp;
    ChainCreate(&chainTemp);
    ChainCopy(&chainTemp,pChainSelect);
    StructureAddChain(&structureTemp, &chainTemp);
    ChainDestroy(&chainTemp);
  }
  Chain *pChainTemp = StructureGetChain(&structureTemp,selectChain);

  int threshiter=trand2*200;
  if(threshiter>1000) threshiter=1000;
  //perturbation [trandpos,trandpos+trand2)
  for(i=trandpos;i<trandpos+trand2;i++){   
    Residue *pResiduePrev = ChainGetResidue(pChainTemp,i-1);
    Residue *pResidue = ChainGetResidue(pChainTemp,i);
    Atom *pAtomNprev = ResidueGetAtomByName(pResiduePrev, "N");
    Atom *pAtomCAprev = ResidueGetAtomByName(pResiduePrev, "CA");
    Atom *pAtomCprev = ResidueGetAtomByName(pResiduePrev, "C");
    Atom *pAtomN = ResidueGetAtomByName(pResidue, "N");
    Atom *pAtomCA = ResidueGetAtomByName(pResidue, "CA");
    Atom *pAtomC = ResidueGetAtomByName(pResidue, "C");

    rmax=0.25+0.25*Random();
    rp=ranv_ergodic(rmax);
    pAtomCA->xyz.X+=rp.X;
    pAtomCA->xyz.Y+=rp.Y;
    pAtomCA->xyz.Z+=rp.Z;
    rmax=0.25+0.25*Random();
    rp=ranv_ergodic(rmax);
    pAtomN->xyz.X+=rp.X;
    pAtomN->xyz.Y+=rp.Y;
    pAtomN->xyz.Z+=rp.Z;
    rmax=0.25+0.25*Random();
    rp=ranv_ergodic(rmax);
    pAtomC->xyz.X+=rp.X;
    pAtomC->xyz.Y+=rp.Y;
    pAtomC->xyz.Z+=rp.Z;
  }
  int numiter=0;
  do{
    for(int k=trandpos;k<trandpos+trand2;k++){
      singlemoveLMP(pChainTemp,k,int(Random()*7.0),0);
      singlemoveLMP(pChainTemp,k,int(Random()*7.0),0);
      singlemoveLMP(pChainTemp,k,int(Random()*7.0),0);
      singlemoveLMP(pChainTemp,k,int(Random()*7.0),0);
      singlemoveLMP(pChainTemp,k,int(Random()*7.0),0);
      singlemoveLMP(pChainTemp,k,int(Random()*7.0),0);
      singlemoveLMP(pChainTemp,k,int(Random()*7.0),0);
    }
    ////////////////////////////////////////////////////////////////////////////////////////inverse         
    for(int k=trandpos+trand2;k>trandpos;k--){
      singlemoveLMP(pChainTemp,k,int(Random()*7.0),1);
      singlemoveLMP(pChainTemp,k,int(Random()*7.0),1);
      singlemoveLMP(pChainTemp,k,int(Random()*7.0),1);
      singlemoveLMP(pChainTemp,k,int(Random()*7.0),1);
      singlemoveLMP(pChainTemp,k,int(Random()*7.0),1);
      singlemoveLMP(pChainTemp,k,int(Random()*7.0),1);
      singlemoveLMP(pChainTemp,k,int(Random()*7.0),1);
    }
    ///////////////////////////////////////////////////////////////////////////////////////check
    flagdone=TRUE;
    for(int k=trandpos;k<=trandpos+trand2;k++){
      Residue *pResiduePrev = ChainGetResidue(pChainTemp,k-1);
      Residue *pResidue = ChainGetResidue(pChainTemp,k);
      Atom *pAtomNprev = ResidueGetAtomByName(pResiduePrev, "N");
      Atom *pAtomCAprev = ResidueGetAtomByName(pResiduePrev, "CA");
      Atom *pAtomCprev = ResidueGetAtomByName(pResiduePrev, "C");
      Atom *pAtomN = ResidueGetAtomByName(pResidue, "N");
      Atom *pAtomCA = ResidueGetAtomByName(pResidue, "CA");
      Atom *pAtomC = ResidueGetAtomByName(pResidue, "C");

      //ca-1 ca
      rp=XYZDifference(&pAtomCAprev->xyz,&pAtomCA->xyz);
      rnorm=XYZNormalization(&rp);
      if(rnorm<lencaca-delcaca || rnorm>lencaca+delcaca){
        flagdone=FALSE;
        break;
      }
      //ca-1 n
      rp=XYZDifference(&pAtomCAprev->xyz,&pAtomN->xyz);
      rnorm=XYZNormalization(&rp);
      if(rnorm<lencan1-delcan1 || rnorm>lencan1+delcan1){
        flagdone=FALSE;
        break;
      }
      //c-1 n
      rp=XYZDifference(&pAtomCprev->xyz,&pAtomN->xyz);
      rnorm=XYZNormalization(&rp);
      if(rnorm<lencn-delcn || rnorm>lencn+delcn){
        flagdone=FALSE;
        break;
      }
      //c-1 ca
      rp=XYZDifference(&pAtomCprev->xyz,&pAtomCA->xyz);
      rnorm=XYZNormalization(&rp);
      if(rnorm<lencca1-delcca1 || rnorm>lencca1+delcca1){
        flagdone=FALSE;
        break;
      }
      if(k==trandpos+trand2){
        break;
      }
      //n ca
      rp=XYZDifference(&pAtomN->xyz,&pAtomCA->xyz);
      rnorm=XYZNormalization(&rp);
      if(rnorm<lennca-delnca || rnorm>lennca+delnca){
        flagdone=FALSE;
        break;
      }
      //n c
      rp=XYZDifference(&pAtomN->xyz,&pAtomC->xyz);
      rnorm=XYZNormalization(&rp);
      if(rnorm<lennc-delnc || rnorm>lennc+delnc){
        flagdone=FALSE;
        break;
      }
      //ca c
      rp=XYZDifference(&pAtomCA->xyz,&pAtomC->xyz);
      rnorm=XYZNormalization(&rp);
      if(rnorm<lencac-delcac || rnorm>lencac+delcac){
        flagdone=FALSE;
        break;
      }
    }
    numiter++;
  }while(!flagdone && numiter<threshiter);
  if(!flagdone){
    //nummov[11][1]++;
    return FALSE;
  }
  
  str2torp(pChainTemp,&chainFlexTemp,trandpos,trandpos+trand2);
  tor2strp2p2(pChainTemp,&chainFlexTemp,trandpos,trandpos+trand2);

  if(trandpos==0){   //should rebuild all chains in structure to fix problem with NewMainChainPos
    for(int j=trandpos; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  }
  else{
    for(int j=trandpos-1; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  }

  //TODO: optimize interface sidechains
  OptimizeSideChains(&structureTemp,bbrotlib,numSideChainOptCycles,selectChain); 
  EnergyTermInitialize(energyTerms);
  StructureCalcEnergyFlex(&structureTemp,bbrotlib,energyTerms);
  EnergyTermWeighting(energyTerms);
  newenergy=energyTerms[0];
//  newenergy+=calcConstraintEnergy(pChainTemp,natDistRestr)*75;
//  newenergy+=calcDistanceEnergy(pChainTemp,distanceCB,distanceWeight);
  if(newenergy<oldenergy){
    oldenergy=newenergy;
    //nummov[angtype][0]++;
    ChainFlex *origChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlexCopy(origChainFlex,&chainFlexTemp);
    Chain *origChain = StructureGetChain(decstr,selectChain);
    ChainCopy(origChain,StructureGetChain(&structureTemp,selectChain));
    getRamaOutliers(origChain,origChainFlex,ramaduke);
    getAngleOutliers(origChain,origChainFlex);
    getBondLengthOutliers(origChain,origChainFlex);
    getClashes(origChain,origChainFlex);
    StructureDestroy(&structureTemp);
    return TRUE;
  }
  else{
    StructureDestroy(&structureTemp);
    return FALSE;
  }
}

BOOL mcfragsweepLMP2(Chain *pChain,int poss,int pose){
  int trandpos=poss;
  int trand2=pose-poss;
  int numiter=0;
  int k;
  XYZ rp;
  double rnorm;
  BOOL flagdone;
  double lencn = 1.338f,delcn = 0.005f*15,lencan1 = 2.441f,delcan1 = 0.036f*4;
  double lennca = 1.460f,delnca = 0.004f*20,lencca1 = 2.446f,delcca1 = 0.036f*4;
  double lencaca = 3.813f,delcaca = 0.019f*10,lencac = 1.525f,delcac = 0.004f*20;
  double lennc = 2.460f,delnc = 0.012f*10;

  do{
    for(k=trandpos;k<trandpos+trand2;k++){
      singlemoveLMP(pChain,k,int(Random()*7.0),0);
      singlemoveLMP(pChain,k,int(Random()*7.0),0);
      singlemoveLMP(pChain,k,int(Random()*7.0),0);
      singlemoveLMP(pChain,k,int(Random()*7.0),0);
      singlemoveLMP(pChain,k,int(Random()*7.0),0);
      singlemoveLMP(pChain,k,int(Random()*7.0),0);
      singlemoveLMP(pChain,k,int(Random()*7.0),0);
    }
    ////////////////////////////////////////////////////////////////////////////////////////inverse     
    for(k=trandpos+trand2;k>trandpos;k--){
      singlemoveLMP(pChain,k,int(Random()*7.0),1);
      singlemoveLMP(pChain,k,int(Random()*7.0),1);
      singlemoveLMP(pChain,k,int(Random()*7.0),1);
      singlemoveLMP(pChain,k,int(Random()*7.0),1);
      singlemoveLMP(pChain,k,int(Random()*7.0),1);
      singlemoveLMP(pChain,k,int(Random()*7.0),1);
      singlemoveLMP(pChain,k,int(Random()*7.0),1);
    }
    ///////////////////////////////////////////////////////////////////////////////////////check
    flagdone=TRUE;
    for(k=trandpos;k<=trandpos+trand2;k++){
      Residue *pResiduePrev = ChainGetResidue(pChain,k-1);
      Residue *pResidue = ChainGetResidue(pChain,k);
      Atom *pAtomNprev = ResidueGetAtomByName(pResiduePrev, "N");
      Atom *pAtomCAprev = ResidueGetAtomByName(pResiduePrev, "CA");
      Atom *pAtomCprev = ResidueGetAtomByName(pResiduePrev, "C");
      Atom *pAtomN = ResidueGetAtomByName(pResidue, "N");
      Atom *pAtomCA = ResidueGetAtomByName(pResidue, "CA");
      Atom *pAtomC = ResidueGetAtomByName(pResidue, "C");

      //ca-1 ca
      rp=XYZDifference(&pAtomCAprev->xyz,&pAtomCA->xyz);
      rnorm=XYZNormalization(&rp);
      if(rnorm<lencaca-delcaca || rnorm>lencaca+delcaca){
        flagdone=FALSE;
        break;
      }
      //ca-1 n
      rp=XYZDifference(&pAtomCAprev->xyz,&pAtomN->xyz);
      rnorm=XYZNormalization(&rp);
      if(rnorm<lencan1-delcan1 || rnorm>lencan1+delcan1){
        flagdone=FALSE;
        break;
      }
      //c-1 n
      rp=XYZDifference(&pAtomCprev->xyz,&pAtomN->xyz);
      rnorm=XYZNormalization(&rp);
      if(rnorm<lencn-delcn || rnorm>lencn+delcn){
        flagdone=FALSE;
        break;
      }
      //c-1 ca
      rp=XYZDifference(&pAtomCprev->xyz,&pAtomCA->xyz);
      rnorm=XYZNormalization(&rp);
      if(rnorm<lencca1-delcca1 || rnorm>lencca1+delcca1){
        flagdone=FALSE;
        break;
      }
      if(k==trandpos+trand2){
        break;
      }
      //n ca
      rp=XYZDifference(&pAtomN->xyz,&pAtomCA->xyz);
      rnorm=XYZNormalization(&rp);
      if(rnorm<lennca-delnca || rnorm>lennca+delnca){
        flagdone=FALSE;
        break;
      }
      //n c
      rp=XYZDifference(&pAtomN->xyz,&pAtomC->xyz);
      rnorm=XYZNormalization(&rp);
      if(rnorm<lennc-delnc || rnorm>lennc+delnc){
        flagdone=FALSE;
        break;
      }
      //ca c
      rp=XYZDifference(&pAtomCA->xyz,&pAtomC->xyz);
      rnorm=XYZNormalization(&rp);
      if(rnorm<lencac-delcac || rnorm>lencac+delcac){
        flagdone=FALSE;
        break;
      }
    }
    numiter++;
  }while(!flagdone && numiter<2000);
  if(!flagdone){
    return FALSE;
  }
  else{
    return TRUE;
  }
}

BOOL moveShift(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
             bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,pairaa *natDistRestr,
             double **contactCA,double **contactCB){
  const int CA_INDEX=3;
  int i,j;
  int flagdirect;
  int trandpos,trandpos2,iend,istart,numse;
  int numSideChainOptCycles=1;
  BOOL flagTor,flagok1,flagok2; 
  double trand2;
  double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
  
  int wpos;
  int selectChain = deschn[int(Random()*numchains)];

  ChainFlex* pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
  Chain* pChain = StructureGetChain(decstr,selectChain);

  int totwrong = pChainFlex->numRamaOutliers + pChainFlex->numAngleOutliers[0] + pChainFlex->numAngleOutliers[1]
               + pChainFlex->numAngleOutliers[2] + pChainFlex->numBondLengthOutliers[0] + pChainFlex->numBondLengthOutliers[1]
               + pChainFlex->numBondLengthOutliers[2] + pChainFlex->numBondLengthOutliers[3] + pChainFlex->numRamaOutliers;

  if(totwrong>0 && Random()<0.1+totwrong/double(ChainFlexGetResidueCount(pChainFlex))){
    wpos=int(Random()*totwrong);
    i=-1;
    for(trandpos=0; trandpos<ChainFlexGetResidueCount(pChainFlex); trandpos++){
      ResidueFlex *pResidueFlex = ChainFlexGetResidue(pChainFlex,trandpos);

      if(pResidueFlex->isAngleOutlierAbove[0] || pResidueFlex->isAngleOutlierBelow[0]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isAngleOutlierAbove[1] || pResidueFlex->isAngleOutlierBelow[1]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isAngleOutlierAbove[2] || pResidueFlex->isAngleOutlierBelow[2]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isBondLengthOutlierAbove[0] || pResidueFlex->isBondLengthOutlierBelow[0]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isBondLengthOutlierAbove[1] || pResidueFlex->isBondLengthOutlierBelow[1]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isBondLengthOutlierAbove[2] || pResidueFlex->isBondLengthOutlierBelow[2]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isBondLengthOutlierAbove[3] || pResidueFlex->isBondLengthOutlierBelow[3]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isRamaOutlier) i++;
      if(i==wpos) break;
    }
    if(Random()<0.5)//leftend
    {
      trandpos2=3+int(20*Random());
      trandpos2+=trandpos;
    }
    else//rightend
    {
      trandpos2=trandpos;
      trandpos=3+int(20*Random());
      trandpos=trandpos2-trandpos;
    }
    if(trandpos<3) trandpos=3;
    if(trandpos2>ChainGetResidueCount(pChain)-4) trandpos2=ChainGetResidueCount(pChain)-4;
    if(trandpos2-trandpos<3){
      //nummov[12][1]++;
      return FALSE;
    }
  }
  else{
    //[trandpos, trandpos2] [3 numseq-4]
    numse=0;
    ResidueFlex *pResidueFlex1,*pResidueFlex2,*pResidueFlex3,*pResidueFlex4,*pResidueFlex5,*pResidueFlex6; 
    do{
      //trandpos2=3+int((numseq-10)*Random());
      trandpos2=3+int(20*Random());//problem when short protein
      trandpos=3+int((ChainGetResidueCount(pChain)-trandpos2-6)*Random());
      trandpos2+=trandpos;
      pResidueFlex1 = ChainFlexGetResidue(pChainFlex,trandpos-1);
      pResidueFlex2 = ChainFlexGetResidue(pChainFlex,trandpos);
      pResidueFlex3 = ChainFlexGetResidue(pChainFlex,trandpos+1);
      pResidueFlex4 = ChainFlexGetResidue(pChainFlex,trandpos2-1);
      pResidueFlex5 = ChainFlexGetResidue(pChainFlex,trandpos2);
      pResidueFlex6 = ChainFlexGetResidue(pChainFlex,trandpos2+1);
      numse++;
    }while((trandpos<3 || trandpos2>ChainGetResidueCount(pChain)-4 || pResidueFlex2->ssm!='C' || pResidueFlex5->ssm!='C') && numse<300);
    if(pResidueFlex2->ssm!='C' || pResidueFlex5->ssm!='C'){
      //nummov[12][1]++;
      return FALSE;
    }
    double treject=Random();
    if(((pResidueFlex1->ssm!='C' || pResidueFlex3->ssm!='C') && treject>0.1) || 
       ((pResidueFlex4->ssm!='C' || pResidueFlex6->ssm!='C') && treject>0.1) ||
       trandpos<3 || trandpos2>ChainGetResidueCount(pChain)-4)
    {
      //nummov[12][1]++;
      return FALSE;
    }
  }
  if(trandpos2<trandpos){
    //nummov[12][1]++;
    return FALSE;
  }

  //StructureFlex structureFlexTemp;
  //for(i=0;i<StructureFlexGetChainCount(pStructureFlex);i++){
    //ChainFlex *pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlex chainFlexTemp;
    ChainFlexCreate(&chainFlexTemp);
    ChainFlexCopy(&chainFlexTemp,pChainFlex);
    //StructureFlexAddChain(structureFlexTemp, &chainFlex);
  //}

  Structure structureTemp;
  StructureCreate(&structureTemp);
  for(i=0;i<StructureGetChainCount(decstr);i++){
    Chain *pChainSelect = StructureGetChain(decstr,i);
    Chain chainTemp;
    ChainCreate(&chainTemp);
    ChainCopy(&chainTemp,pChainSelect);
    StructureAddChain(&structureTemp, &chainTemp);
    ChainDestroy(&chainTemp);
  }
  Chain *pChainTemp = StructureGetChain(&structureTemp,selectChain);

  if(Random()<0.5) flagdirect=0;//back
  else flagdirect=1;//forward

  if(flagdirect){
    for(i=trandpos;i<=trandpos2;i++){
      Residue *pResidue = ChainGetResidue(pChainTemp,i);
      Residue *pResidueNext = ChainGetResidue(pChainTemp,i+1);
      Atom *pAtomCA = ResidueGetAtomByName(pResidue, "CA");
      Atom *pAtomN = ResidueGetAtomByName(pResidue, "N");
      Atom *pAtomC = ResidueGetAtomByName(pResidue, "C");
      Atom *pAtomCAnext = ResidueGetAtomByName(pResidueNext, "CA");
      Atom *pAtomNnext = ResidueGetAtomByName(pResidueNext, "N");
      Atom *pAtomCnext = ResidueGetAtomByName(pResidueNext, "C");

      pAtomCA->xyz.X = pAtomCAnext->xyz.X;
      pAtomCA->xyz.Y = pAtomCAnext->xyz.Y;
      pAtomCA->xyz.Z = pAtomCAnext->xyz.Z;
      pAtomN->xyz.X = pAtomNnext->xyz.X;
      pAtomN->xyz.Y = pAtomNnext->xyz.Y;
      pAtomN->xyz.Z = pAtomNnext->xyz.Z;
      pAtomC->xyz.X = pAtomCnext->xyz.X;
      pAtomC->xyz.Y = pAtomCnext->xyz.Y;
      pAtomC->xyz.Z = pAtomCnext->xyz.Z;
    }
  }//gap in [trandpos-1,trandpos] dup in [trandpos2,trandpos2+1]
  else{
    for(i=trandpos2;i>=trandpos;i--){
      Residue *pResidue = ChainGetResidue(pChainTemp,i);
      Residue *pResiduePrev = ChainGetResidue(pChainTemp,i-1);
      Atom *pAtomCA = ResidueGetAtomByName(pResidue, "CA");
      Atom *pAtomN = ResidueGetAtomByName(pResidue, "N");
      Atom *pAtomC = ResidueGetAtomByName(pResidue, "C");
      Atom *pAtomCAprev = ResidueGetAtomByName(pResiduePrev, "CA");
      Atom *pAtomNprev = ResidueGetAtomByName(pResiduePrev, "N");
      Atom *pAtomCprev = ResidueGetAtomByName(pResiduePrev, "C");

      pAtomCA->xyz.X = pAtomCAprev->xyz.X;
      pAtomCA->xyz.Y = pAtomCAprev->xyz.Y;
      pAtomCA->xyz.Z = pAtomCAprev->xyz.Z;
      pAtomN->xyz.X = pAtomNprev->xyz.X;
      pAtomN->xyz.Y = pAtomNprev->xyz.Y;
      pAtomN->xyz.Z = pAtomNprev->xyz.Z;
      pAtomC->xyz.X = pAtomCprev->xyz.X;
      pAtomC->xyz.Y = pAtomCprev->xyz.Y;
      pAtomC->xyz.Z = pAtomCprev->xyz.Z;
    }
  }//dup in [trandpos-1,trandpos] gap in [trandpos2,trandpos2+1]
  flagok1=mcfragsweepLMP2(pChain,trandpos-2,trandpos+2);
  if(!flagok1){
    //nummov[12][1]++;
    return FALSE;
  }
  flagok2=mcfragsweepLMP2(pChain,trandpos2-1,trandpos2+3);
  if(!flagok2){
    //nummov[12][1]++;
    return FALSE;
  }
  istart=trandpos-2;
  iend=trandpos2+3;
  if(iend>=ChainGetResidueCount(pChain)) iend=ChainGetResidueCount(pChain)-1;

  str2torp(pChainTemp,&chainFlexTemp,istart,iend);
  tor2strp2p2(pChainTemp,&chainFlexTemp,istart,iend);

  //flagTor=NewMainChainPos(StructureGetChain(&structureTemp,selectChain),&chainFlexTemp,trandpos);
  if(trandpos==0){   //should rebuild all chains in structure to fix problem with NewMainChainPos
    for(int j=trandpos; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  }
  else{
    for(int j=trandpos-1; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  }

  //TODO: optimize interface sidechains
  OptimizeSideChains(&structureTemp,bbrotlib,numSideChainOptCycles,selectChain); 
  EnergyTermInitialize(energyTerms);
  StructureCalcEnergyFlex(&structureTemp,bbrotlib,energyTerms);
  EnergyTermWeighting(energyTerms);
  newenergy=energyTerms[0];
//  newenergy+=calcConstraintEnergy(pChainTemp,natDistRestr)*75;
//  newenergy+=calcDistanceEnergy(pChainTemp,distanceCB,distanceWeight);
  if(newenergy<oldenergy){
    oldenergy=newenergy;
    //nummov[angtype][0]++;
    ChainFlex *origChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlexCopy(origChainFlex,&chainFlexTemp);
    Chain *origChain = StructureGetChain(decstr,selectChain);
    ChainCopy(origChain,StructureGetChain(&structureTemp,selectChain));
    getRamaOutliers(origChain,origChainFlex,ramaduke);
    getAngleOutliers(origChain,origChainFlex);
    getBondLengthOutliers(origChain,origChainFlex);
    getClashes(origChain,origChainFlex);
    StructureDestroy(&structureTemp);
    return TRUE;
  }
  else{
    StructureDestroy(&structureTemp);
    return FALSE;
  }
}

BOOL movePos(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
             bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,pairaa *natDistRestr,
             double **contactCA,double **contactCB){
  const int CA_INDEX=3;
  int i,j;
  int flagdirect;
  int trandpos,trandpos2,numse;
  int numSideChainOptCycles=1;
  BOOL flagTor,flagok1,flagok2; 
  double trand2;
  double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
  
  int wpos;
  int selectChain = deschn[int(Random()*numchains)];

  ChainFlex* pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
  Chain* pChain = StructureGetChain(decstr,selectChain);

  int totwrong = pChainFlex->numRamaOutliers + pChainFlex->numAngleOutliers[0] + pChainFlex->numAngleOutliers[1]
               + pChainFlex->numAngleOutliers[2] + pChainFlex->numBondLengthOutliers[0] + pChainFlex->numBondLengthOutliers[1]
               + pChainFlex->numBondLengthOutliers[2] + pChainFlex->numBondLengthOutliers[3] + pChainFlex->numRamaOutliers
               + pChainFlex->numClashes; 
  if(totwrong>0 && Random()<0.1+totwrong/double(ChainFlexGetResidueCount(pChainFlex))){
    wpos=int(Random()*totwrong);
    i=-1;
    for(trandpos=0; trandpos<ChainFlexGetResidueCount(pChainFlex); trandpos++){
      ResidueFlex *pResidueFlex = ChainFlexGetResidue(pChainFlex,trandpos);

      if(pResidueFlex->isClashed) i++;
      if(i==wpos) break;
      if(pResidueFlex->isAngleOutlierAbove[0] || pResidueFlex->isAngleOutlierBelow[0]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isAngleOutlierAbove[1] || pResidueFlex->isAngleOutlierBelow[1]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isAngleOutlierAbove[2] || pResidueFlex->isAngleOutlierBelow[2]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isBondLengthOutlierAbove[0] || pResidueFlex->isBondLengthOutlierBelow[0]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isBondLengthOutlierAbove[1] || pResidueFlex->isBondLengthOutlierBelow[1]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isBondLengthOutlierAbove[2] || pResidueFlex->isBondLengthOutlierBelow[2]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isBondLengthOutlierAbove[3] || pResidueFlex->isBondLengthOutlierBelow[3]) i++;
      if(i==wpos) break;
      if(pResidueFlex->isRamaOutlier) i++;
      if(i==wpos) break;
    }
    if(trandpos>ChainGetResidueCount(pChain)-1){
      trandpos=ChainGetResidueCount(pChain)-1;
    }
    else if(trandpos<1){
      trandpos=1;
    }
  }
  else// trandpos [1,L-1]
  {
    trandpos=1+int((ChainGetResidueCount(pChain)-1)*Random());
    ResidueFlex *pResidueFlex = ChainFlexGetResidue(pChainFlex,trandpos);
    if(pResidueFlex->ssm!='C' && Random()>0.1){
      //nummov[18][1]++;
      return FALSE;
    }
  }

  //StructureFlex structureFlexTemp;
  //for(i=0;i<StructureFlexGetChainCount(pStructureFlex);i++){
    //ChainFlex *pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlex chainFlexTemp;
    ChainFlexCreate(&chainFlexTemp);
    ChainFlexCopy(&chainFlexTemp,pChainFlex);
    //StructureFlexAddChain(structureFlexTemp, &chainFlex);
  //}

  Structure structureTemp;
  StructureCreate(&structureTemp);
  for(i=0;i<StructureGetChainCount(decstr);i++){
    Chain *pChainSelect = StructureGetChain(decstr,i);
    Chain chainTemp;
    ChainCreate(&chainTemp);
    ChainCopy(&chainTemp,pChainSelect);
    StructureAddChain(&structureTemp, &chainTemp);
    ChainDestroy(&chainTemp);
  }
  Chain *pChainTemp = StructureGetChain(&structureTemp,selectChain);

  wpos=int(3*Random());
  XYZ tp;
  tp=ranv_ergodic(0.2);
  Residue *pResidue = ChainGetResidue(pChainTemp,trandpos);
  Atom *selectAtom;
  if(wpos==0) selectAtom = ResidueGetAtomByName(pResidue, "N");
  else if(wpos==1) selectAtom = ResidueGetAtomByName(pResidue, "CA");
  else if(wpos==2) selectAtom = ResidueGetAtomByName(pResidue, "C");
  selectAtom->xyz.X+=tp.X;
  selectAtom->xyz.Y+=tp.Y;
  selectAtom->xyz.Z+=tp.Z;
  int istart=trandpos-1;
  if(istart<0) istart=0;
  int iend=trandpos+1;
  if(iend>ChainGetResidueCount(pChainTemp)-1) iend=ChainGetResidueCount(pChainTemp)-1;
  str2torp(pChainTemp,&chainFlexTemp,istart,iend);
  tor2strp2p2(pChainTemp,&chainFlexTemp,istart,iend);
  //flagTor=NewMainChainPos(StructureGetChain(&structureTemp,selectChain),&chainFlexTemp,trandpos);

  if(trandpos==0){   //should rebuild all chains in structure to fix problem with NewMainChainPos
    for(int j=trandpos; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  }
  else{
    for(int j=trandpos-1; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  }

  //TODO: optimize interface sidechains
  OptimizeSideChains(&structureTemp,bbrotlib,numSideChainOptCycles,selectChain); 
  EnergyTermInitialize(energyTerms);
  StructureCalcEnergyFlex(&structureTemp,bbrotlib,energyTerms);
  EnergyTermWeighting(energyTerms);
  newenergy=energyTerms[0];
//  newenergy+=calcConstraintEnergy(pChainTemp,natDistRestr)*75;
//  newenergy+=calcDistanceEnergy(pChainTemp,distanceCB,distanceWeight);
  if(newenergy<oldenergy){
    oldenergy=newenergy;
    //nummov[angtype][0]++;
    ChainFlex *origChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlexCopy(origChainFlex,&chainFlexTemp);
    Chain *origChain = StructureGetChain(decstr,selectChain);
    ChainCopy(origChain,StructureGetChain(&structureTemp,selectChain));
    getRamaOutliers(origChain,origChainFlex,ramaduke);
    getAngleOutliers(origChain,origChainFlex);
    getBondLengthOutliers(origChain,origChainFlex);
    getClashes(origChain,origChainFlex);
    StructureDestroy(&structureTemp);
    return TRUE;
  }
  else{
    StructureDestroy(&structureTemp);
    return FALSE;
  }
}

BOOL mcfragsweepCCD3(Chain *pChain,ChainFlex *pChainFlex,int lps, int lpe,int lpt,XYZ pt1,
                     XYZ pt2,XYZ pt3,double oldenergy,double *newenergy){
  int i;
  int trandpos;
  BOOL flagphi,flagtor;
  int flagpt;
  XYZ tp[3];//dest for lpt-1 lpt lpt+1
  tp[0].X=pt1.X;tp[0].Y=pt1.Y;tp[0].Z=pt1.Z;
  tp[1].X=pt2.X;tp[1].Y=pt2.Y;tp[1].Z=pt2.Z;
  tp[2].X=pt3.X;tp[2].Y=pt3.Y;tp[2].Z=pt3.Z;
  XYZ pcur,p12,p13;
  int numiter=0;
  double tdist[3],ttheta,tphi,tpsi,tinner;

  flagpt=0;
  int threshiter=25*(lpe-lps+1);
  if(threshiter>600) threshiter=600;
  do{
    tdist[0]=10000;
    tdist[1]=10000;
    tdist[2]=10000;
    trandpos=lps+int((lpe-lps+1)*Random());
    ResidueFlex *pResidueFlex = ChainFlexGetResidue(pChainFlex,trandpos);
    ResidueFlex *pResidueFlexNext = ChainFlexGetResidue(pChainFlex,trandpos+1);

    if(Random()<0.5) flagphi=TRUE;//change phi n ca
    else flagphi=FALSE;//change psi ca c in pos+1
    if(pResidueFlex->ss=='H' && pResidueFlexNext->ss=='H'){
      numiter++;
      continue;
    }
    else if(pResidueFlex->ss=='H'){
      flagphi=FALSE;
    }
    else if(pResidueFlex->ss=='E' && pResidueFlexNext->ss=='E'){
      numiter++;
      continue;
    }
    else if(pResidueFlex->ss=='E'){
      flagphi=FALSE;
    }
    flagpt=(flagpt+1)%3;
    Residue *pResidue1 = ChainGetResidue(pChain,trandpos);
    Residue *pResidue2 = ChainGetResidue(pChain,lpt+flagpt-1);
    Atom *pAtomCA1 = ResidueGetAtomByName(pResidue1, "CA");
    Atom *pAtomN1 = ResidueGetAtomByName(pResidue1, "N");
    Atom *pAtomC1 = ResidueGetAtomByName(pResidue1, "C");
    Atom *pAtomCA2 = ResidueGetAtomByName(pResidue2, "CA");
    Atom *pAtomN2 = ResidueGetAtomByName(pResidue2, "N");
    Atom *pAtomC2 = ResidueGetAtomByName(pResidue2, "C");
    if(flagphi){

      p12=XYZDifference(&pAtomN1->xyz,&pAtomCA1->xyz);  
      p13=XYZDifference(&pAtomN1->xyz,&pAtomCA2->xyz);  
      tinner=XYZAngle(&p12,&p13);
      if(tinner<err || PI-tinner<err || XYZNormalization(&p12)<err)//in one line no affect when rotate
      {
        numiter++;
        continue;
      }
      pcur.X=pAtomCA2->xyz.X;pcur.Y=pAtomCA2->xyz.Y;pcur.Z=pAtomCA2->xyz.Z;
      p13.X=pAtomCA1->xyz.X;p13.Y=pAtomCA1->xyz.Y;p13.Z=pAtomCA1->xyz.Z;
      p12.X=pAtomN1->xyz.X;p12.Y=pAtomN1->xyz.Y;p12.Z=pAtomN1->xyz.Z;
      ttheta=calcTorsion(pcur.X,pcur.Y,pcur.Z,p12.X,p12.Y,p12.Z,p13.X,p13.Y,p13.Z,tp[flagpt].X,tp[flagpt].Y,tp[flagpt].Z);
      tphi=pResidueFlex->dih[2]+ttheta;
      if(tphi>360.0) tphi-=360.0;
      pResidueFlex->dih[2]=tphi;  
    }
    else{
      p12=XYZDifference(&pAtomCA1->xyz,&pAtomC1->xyz);
      p13=XYZDifference(&pAtomCA1->xyz,&pAtomCA2->xyz);
      tinner=XYZAngle(&p12,&p13);
      if(tinner<err || PI-tinner<err || XYZNormalization(&p12)<err)//in one line no affect when rotate
      {
        numiter++;
        continue;
      }
      pcur.X=pAtomCA2->xyz.X;pcur.Y=pAtomCA2->xyz.Y;pcur.Z=pAtomCA2->xyz.Z;
      p13.X=pAtomC1->xyz.X;p13.Y=pAtomC1->xyz.Y;p13.Z=pAtomC1->xyz.Z;
      p12.X=pAtomCA1->xyz.X;p12.Y=pAtomCA1->xyz.Y;p12.Z=pAtomCA1->xyz.Z;
      ttheta=calcTorsion(pcur.X,pcur.Y,pcur.Z,p12.X,p12.Y,p12.Z,p13.X,p13.Y,p13.Z,tp[flagpt].X,tp[flagpt].Y,tp[flagpt].Z);
      tpsi=pResidueFlexNext->dih[0]+ttheta;
      if(tpsi>360.0) tpsi-=360.0;
      pResidueFlexNext->dih[0]=tpsi;
    }
    numiter++;

    /////here here here 
    flagtor=NewMainChainPos(pChain,pChainFlex,trandpos);
    //if(trandpos==0) flagtor=tor2str(pChain,pChainFlex);
    //else flagtor=ps.tor2strp(tmstr2,numseq,trandpos);
    if(!flagtor) printf("tor2str wrong in CCD3 %d\n",trandpos);
    for(i=0;i<3;i++){
      Residue *pResidue = ChainGetResidue(pChain,lpt+flagpt-1);
      Atom *pAtomCA = ResidueGetAtomByName(pResidue, "CA");
      pcur.X=pAtomCA->xyz.X-tp[i].X;
      pcur.Y=pAtomCA->xyz.Y-tp[i].Y;
      pcur.Z=pAtomCA->xyz.Z-tp[i].Z;
      tdist[i]=pcur.X*pcur.X+pcur.Y*pcur.Y+pcur.Z*pcur.Z;
    }
  }while(numiter<threshiter && (tdist[0]>2.50 || tdist[1]>2.50 || tdist[2]>2.50));

  if(numiter<threshiter){
    return TRUE;
  }
  else{
    return FALSE;
  }
}

BOOL mcfragsweepCCD4(Chain *pChain,ChainFlex *pChainFlex,int lps,int lpe,XYZ pt1,int ind1,XYZ pt2,int ind2,double oldenergy,double *newenergy){
  int trandpos;
  BOOL flagphi,flagtor;
  int flagpt,indp[2];
  indp[0]=ind1;indp[1]=ind2;
  XYZ tp[2];
  tp[0].X=pt1.X;tp[0].Y=pt1.Y;tp[0].Z=pt1.Z;
  tp[1].X=pt2.X;tp[1].Y=pt2.Y;tp[1].Z=pt2.Z;
  XYZ pcur,p12,p13;
  int numiter=0;
  double tdist[2],ttheta,tphi,tpsi,tinner;
  flagpt=0;
  int threshiter=25*(lpe-lps+1);
  if(threshiter>600) threshiter=600;
  do{
    tdist[0]=10000;
    tdist[1]=10000;
    trandpos=lps+int((lpe-lps+1)*Random());
    if(Random()<0.5) flagphi=TRUE;//change phi n ca
    else flagphi=FALSE;//change psi ca c in pos+1
    ResidueFlex *pResidueFlex1 = ChainFlexGetResidue(pChainFlex,trandpos);
    ResidueFlex *pResidueFlex2 = ChainFlexGetResidue(pChainFlex,trandpos+1);

    if(pResidueFlex1->ss=='H' && pResidueFlex2->ss=='H'){   
      numiter++;
      continue;
    }
    else if(pResidueFlex1->ss=='H'){   
      flagphi=FALSE;
    }
    else if(pResidueFlex1->ss=='E' && pResidueFlex2->ss=='E'){   
      numiter++;
      continue;
    }
    else if(pResidueFlex1->ss=='E'){
      flagphi=FALSE;
    }
    flagpt=(flagpt+1)%2;
    if(flagphi){
      Residue *pResidue1 = ChainGetResidue(pChain,trandpos);
      Residue *pResidue2 = ChainGetResidue(pChain,indp[flagpt]);
      Atom *pAtomCA1 = ResidueGetAtomByName(pResidue1, "CA");      
      Atom *pAtomN1 = ResidueGetAtomByName(pResidue1, "N");     
      Atom *pAtomCA2 = ResidueGetAtomByName(pResidue2, "CA");      
      Atom *pAtomN2 = ResidueGetAtomByName(pResidue2, "N");     
       
      p12=XYZDifference(&pAtomN1->xyz,&pAtomCA1->xyz);
      p13=XYZDifference(&pAtomN1->xyz,&pAtomCA2->xyz);
      tinner=XYZAngle(&p12,&p13);
      if(tinner<err || PI-tinner<err || XYZNormalization(&p12)<err){//in one line no affect when rotate
        numiter++;
        continue;
      }
      pcur.X=pAtomCA2->xyz.X;pcur.Y=pAtomCA2->xyz.Y;pcur.Z=pAtomCA2->xyz.Z;
      p13.X=pAtomCA1->xyz.X;p13.Y=pAtomCA1->xyz.Y;p13.Z=pAtomCA1->xyz.Z;
      p12.X=pAtomN1->xyz.X;p12.Y=pAtomN1->xyz.Y;p12.Z=pAtomN1->xyz.Z;
      ttheta=calcTorsion(pcur.X,pcur.Y,pcur.Z,
                         p12.X,p12.Y,p12.Z,
                         p13.X,p13.Y,p13.Z,
                         tp[flagpt].X,tp[flagpt].Y,tp[flagpt].Z);
      tphi=pResidueFlex1->dih[2]+ttheta;
      if(tphi>360.0) tphi-=360.0;
      pResidueFlex1->dih[2]=tphi;
    }
    else{
      Residue *pResidue1 = ChainGetResidue(pChain,trandpos);
      Residue *pResidue2 = ChainGetResidue(pChain,indp[flagpt]);
      Atom *pAtomCA1 = ResidueGetAtomByName(pResidue1, "CA");      
      Atom *pAtomC1 = ResidueGetAtomByName(pResidue1, "C");     
      Atom *pAtomCA2 = ResidueGetAtomByName(pResidue2, "CA");      
      Atom *pAtomN2 = ResidueGetAtomByName(pResidue2, "N");     

      p12=XYZDifference(&pAtomCA1->xyz,&pAtomC1->xyz);
      p13=XYZDifference(&pAtomCA1->xyz,&pAtomCA2->xyz);
      tinner=XYZAngle(&p12,&p13);
      if(tinner<err || PI-tinner<err || XYZNormalization(&p12)<err){//in one line no affect when rotate
        numiter++;
        continue;
      }
      pcur.X=pAtomCA2->xyz.X;pcur.Y=pAtomCA2->xyz.Y;pcur.Z=pAtomCA2->xyz.Z;
      p12.X=pAtomCA1->xyz.X;p12.Y=pAtomCA1->xyz.Y;p12.Z=pAtomCA1->xyz.Z;
      p13.X=pAtomC1->xyz.X;p13.Y=pAtomC1->xyz.Y;p13.Z=pAtomC1->xyz.Z;
      ttheta=calcTorsion(pcur.X,pcur.Y,pcur.Z,
                         p12.X,p12.Y,p12.Z,
                         p13.X,p13.Y,p13.Z,
                         tp[flagpt].X,tp[flagpt].Y,tp[flagpt].Z);
      tpsi=pResidueFlex2->dih[0]+ttheta;
      if(tpsi>360.0) tpsi-=360.0;
      pResidueFlex2->dih[0]=tpsi;
    }
    numiter++;
    flagtor=NewMainChainPos(pChain,pChainFlex,trandpos);
    if(!flagtor){
      printf("tor2str wrong in CCD4 %d\n",trandpos);
    }
    Residue *pResidue0 = ChainGetResidue(pChain,0);
    Residue *pResidue1 = ChainGetResidue(pChain,ind1);
    Residue *pResidue2 = ChainGetResidue(pChain,1);
    Residue *pResidue3 = ChainGetResidue(pChain,ind2);
    Atom *pAtomCA0 = ResidueGetAtomByName(pResidue0, "CA");      
    Atom *pAtomCA1 = ResidueGetAtomByName(pResidue1, "CA");      
    Atom *pAtomCA2 = ResidueGetAtomByName(pResidue2, "CA");      
    Atom *pAtomCA3 = ResidueGetAtomByName(pResidue3, "CA");      

    pcur.X=pAtomCA1->xyz.X-pAtomCA0->xyz.X;
    pcur.Y=pAtomCA1->xyz.Y-pAtomCA0->xyz.Y;
    pcur.Z=pAtomCA1->xyz.Z-pAtomCA0->xyz.Z;
    tdist[0]=pcur.X*pcur.X+pcur.Y*pcur.Y+pcur.Z*pcur.Z;

    pcur.X=pAtomCA3->xyz.X-pAtomCA2->xyz.X;
    pcur.Y=pAtomCA3->xyz.Y-pAtomCA2->xyz.Y;
    pcur.Z=pAtomCA3->xyz.Z-pAtomCA2->xyz.Z;
    tdist[1]=pcur.X*pcur.X+pcur.Y*pcur.Y+pcur.Z*pcur.Z;
    // printf("%d %d %f %f %f\n",numiter,threshiter,tdist[0],tdist[1],tdist[2]);
  }while(numiter<threshiter && (tdist[0]>6.50 || tdist[1]>2.50));
  // printf("%d [%d %d %d] pos %d %d dist %f theta %f\n",numiter,lps,lpe,lpt,trandpos, flagphi, tdist,ttheta);
  if(numiter<threshiter){
    return TRUE;
  }
  else{
    return FALSE;
  }
}

void getfrombbptable(double tval,segment *sgval,paircont* tbbptb,int tnumbbptb){
  //(i i+1]
  int i;
  if(tnumbbptb<2){
    if(tval>tbbptb[0].dist) printf("wrong %f %f\n",tval,tbbptb[0].dist);
    sgval->init=tbbptb[0].init;
    sgval->term=tbbptb[0].term;
    return;
  }
  i=tnumbbptb/2;
  if(tval>tbbptb[i-1].dist) getfrombbptable(tval,sgval,tbbptb+i,tnumbbptb-i);
  else getfrombbptable(tval,sgval,tbbptb,i);
}

BOOL moveBetaPack(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
                  bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,
                  paircont *bbptb,int numbbptb,pairaa *natDistRestr,double **contactCA,double **contactCB){
  const int CA_INDEX=3;
  int i,j;
  int flagdirect;
  int trandpos,trandpos2,numse;
  int numSideChainOptCycles=1;
  BOOL flagTor,flagok1,flagok2; 
  double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};
  
  int trand,trand2;
  int typebb;
  bool flagat,flaglr;
  XYZ tp[15],ap[15],fp[3];
  double angle1,angle2,angle3;
  double rmat[9],irmat[9],rmat2[9];
  double toldenergy=oldenergy;
  double tnewenergy;
  bool flagclash=false;
  bool flagtor;
  double threshclash=0.01;
  segment sgval;
  if(numbbptb<3){   
    //summcbbp[2]++;
    return FALSE;
  }

  int selectChain = deschn[int(Random()*numchains)]; 
  ChainFlex* pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);

  //StructureFlex structureFlexTemp;
  //for(i=0;i<StructureFlexGetChainCount(pStructureFlex);i++){
    //ChainFlex *pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlex chainFlexTemp;
    ChainFlexCreate(&chainFlexTemp);
    ChainFlexCopy(&chainFlexTemp,pChainFlex);
    //StructureFlexAddChain(structureFlexTemp, &chainFlex);
  //}

  Structure structureTemp;
  StructureCreate(&structureTemp);
  for(i=0;i<StructureGetChainCount(decstr);i++){
    Chain *pChainSelect = StructureGetChain(decstr,i);
    Chain chainTemp;
    ChainCreate(&chainTemp);
    ChainCopy(&chainTemp,pChainSelect);
    StructureAddChain(&structureTemp, &chainTemp);
    ChainDestroy(&chainTemp);
  }
  Chain *pChainTemp = StructureGetChain(&structureTemp,selectChain);

  do{
    getfrombbptable(Random(),&sgval,bbptb,numbbptb);
    trand=sgval.init;
    trand2=sgval.term;

    ResidueFlex *pResidueFlex1 = ChainFlexGetResidue(&chainFlexTemp,trand);
    ResidueFlex *pResidueFlex2 = ChainFlexGetResidue(&chainFlexTemp,trand2);

    if(abs(pResidueFlex1->indr-trand2)<3 || abs(pResidueFlex1->indl-trand2)<3){ //already paired
      //summcbbp[2]++;
      return FALSE;
    }
    if(Random()<0.75) flagat=FALSE;//anti
    else flagat=TRUE;//para
    if(Random()<0.5)
    flaglr=FALSE;//left
    else flaglr=TRUE;//right
    if     (flagat  && !flaglr) typebb=0;//paraleft
    else if(flagat  && flaglr)  typebb=1;//pararight
    else if(!flagat && !flaglr) typebb=2;//antileft
    else if(!flagat && flaglr)  typebb=3;//antiright

    if(typebb==0 && (pResidueFlex1->indl!=-1 || pResidueFlex2->indr!=-1)){
      flaglr=TRUE;
      typebb=1;
    }
    else if(typebb==1 && (pResidueFlex1->indr!=-1 || pResidueFlex2->indl!=-1)){
      flaglr=FALSE;
      typebb=0;
    }
    else if(typebb==2 && (pResidueFlex1->indl!=-1 || pResidueFlex2->indr!=-1)){
      flaglr=TRUE;
      typebb=3;
    }
    else if(typebb==3 && (pResidueFlex1->indr!=-1 || pResidueFlex2->indl!=-1)){
      flaglr=FALSE;
      typebb=2;
    }

    //flagat
    if(typebb==0 && (pResidueFlex1->indl!=-1 || pResidueFlex2->indr!=-1)){
      flagat=FALSE;
      typebb=2;
    }
    else if(typebb==1 && (pResidueFlex1->indr!=-1 || pResidueFlex2->indl!=-1)){
      flagat=FALSE;
      typebb=3;
    }
    else if(typebb==2 && (pResidueFlex1->indl!=-1 || pResidueFlex2->indr!=-1)){
      flagat=TRUE;
      typebb=0;
    }
    else if(typebb==3 && (pResidueFlex1->indr!=-1 || pResidueFlex2->indl!=-1)){
      flagat=TRUE;
      typebb=1;
    }
    if(typebb==0 && (pResidueFlex1->indr!=-1 && pResidueFlex1->tpr==3)){
      if(Random()<0.90){
        flagat=FALSE;
        typebb=2;
      }
    }
    else if(typebb==1 && (pResidueFlex1->indl!=-1 && pResidueFlex1->tpl==2)){
      if(Random()<0.90){
        flagat=FALSE;
        typebb=3;
      }
    }
    else if(typebb==2 && (pResidueFlex1->indr!=-1 && pResidueFlex1->tpr==1)){
      if(Random()<0.90){
        flagat=TRUE;
        typebb=0;
      }
    }
    else if(typebb==3 && (pResidueFlex1->indl!=-1 && pResidueFlex1->tpl==0)){
      if(Random()<0.90){
        flagat=TRUE;
        typebb=1;
      }
    }

    Residue *pResidue1 = ChainGetResidue(pChainTemp,trand-1);
    Residue *pResidue2 = ChainGetResidue(pChainTemp,trand);
    Residue *pResidue3 = ChainGetResidue(pChainTemp,trand+1);

    Atom *pAtomCA1 = ResidueGetAtomByName(pResidue1, "CA");
    Atom *pAtomCA2 = ResidueGetAtomByName(pResidue2, "CA");
    Atom *pAtomCA3 = ResidueGetAtomByName(pResidue3, "CA");

    tp[0].X=pAtomCA1->xyz.X;tp[0].Y=pAtomCA1->xyz.Y;tp[0].Z=pAtomCA1->xyz.Z;
    tp[1].X=pAtomCA2->xyz.X;tp[1].Y=pAtomCA2->xyz.Y;tp[1].Z=pAtomCA2->xyz.Z;
    tp[2].X=pAtomCA3->xyz.X;tp[2].Y=pAtomCA3->xyz.Y;tp[2].Z=pAtomCA3->xyz.Z;
    tp[3]=XYZDifference(&tp[1],&tp[0]);
    tp[4]=XYZDifference(&tp[1],&tp[2]);
    tp[5]=XYZDifference(&tp[4],&tp[3]);
    tp[5]=XYZUnit(&tp[5]);
    if(XYZNormalization(&tp[5])<1e-10){
      //summcbbp[2]++;
      return FALSE;
    }

    tp[6]=XYZDifference(&tp[2],&tp[0]);
    tp[6]=XYZUnit(&tp[6]);
    v2rot(tp[5],rmat);
    copyMat(rmat,3,3,irmat);
    rinv(rmat,3);//tp5 to 001
    tp[10]=XYZMmat(rmat,&tp[6]);
    angle1=atan2(tp[10].Y,tp[10].X);

    /////////////////////////////////////////////////////////////////////
    Residue *pResidue4 = ChainGetResidue(pChainTemp,trand2-1);
    Residue *pResidue5 = ChainGetResidue(pChainTemp,trand2);
    Residue *pResidue6 = ChainGetResidue(pChainTemp,trand2+1);

    Atom *pAtomCA4 = ResidueGetAtomByName(pResidue4, "CA");
    Atom *pAtomCA5 = ResidueGetAtomByName(pResidue5, "CA");
    Atom *pAtomCA6 = ResidueGetAtomByName(pResidue6, "CA");
 
    ap[0].X=pAtomCA4->xyz.X;ap[0].Y=pAtomCA4->xyz.Y;ap[0].Z=pAtomCA4->xyz.Z;
    ap[1].X=pAtomCA5->xyz.X;ap[1].Y=pAtomCA5->xyz.Y;ap[1].Z=pAtomCA5->xyz.Z;
    ap[2].X=pAtomCA6->xyz.X;ap[2].Y=pAtomCA6->xyz.Y;ap[2].Z=pAtomCA6->xyz.Z;

    ap[3]=XYZDifference(&ap[1],&ap[0]);
    ap[4]=XYZDifference(&ap[1],&ap[2]);
    ap[5]=XYZCrossProduct(&ap[3],&ap[4]);
    ap[5]=XYZUnit(&ap[5]);
    if(XYZNormalization(&ap[5])<1e-10){
      //summcbbp[2]++;
      return FALSE;
    }

    ap[6]=XYZDifference(&ap[2],&ap[0]);
    ap[6]=XYZUnit(&ap[6]);
    if(!flagat) ap[5]=XYZScale2(&ap[5],-1);
    v2rot(ap[5],rmat);
    rinv(rmat,3);//ap5 to 001
    ap[10]=XYZMmat(rmat,&ap[6]);
    ap[11]=XYZMmat(rmat,&ap[3]);
    ap[12]=XYZMmat(rmat,&ap[4]);
    angle2=atan2(ap[10].Y,ap[10].X);
    if(!flagat) angle3=angle1-angle2+PI+Random()*0.20-0.10;
    else angle3=angle1-angle2+Random()*0.20-0.10;
    a2rot(angle3,rmat2);
    ap[10]=XYZMmat(rmat2,&ap[10]);
    angle2=atan2(ap[10].Y,ap[10].X);
    ap[13]=XYZMmat(rmat2,&ap[11]);
    ap[14]=XYZMmat(rmat2,&ap[12]);
    tp[13]=XYZMmat(irmat,&ap[13]);
    tp[14]=XYZMmat(irmat,&ap[14]);
    tp[11]=tp[13];
    tp[12]=tp[14];
    int type;
    if(flagat && !flaglr){ //paraleft norm same type0
      type=1;
      tp[9]=tp[5];
      tp[5]=XYZScale2(&tp[9],4.81+0.1*Random());
      fp[1]=XYZDifference(&tp[5],&tp[1]);
    }
    else if(flagat && flaglr){ //pararight norm same type1
      type=2;
      tp[9]=tp[5];
      tp[5]=XYZScale2(&tp[9],4.86+0.1*Random());
      fp[1]=XYZSum(&tp[1],&tp[5]);
    }
    else if(!flagat && !flaglr){ //antileft norm diff type2
      type=3;
      tp[9]=tp[5];
      tp[5]=XYZScale2(&tp[9],4.49+0.1*Random());
      fp[1]=XYZDifference(&tp[5],&tp[1]);
    }
    else if(!flagat && flaglr){ //antiright norm diff type3
      type=4;
      tp[9]=tp[5];
      tp[5]=XYZScale2(&tp[9],5.24+0.1*Random());
      fp[1]=XYZSum(&tp[1],&tp[5]);
    }
    fp[0]=XYZSum(&fp[1],&tp[11]);
    fp[2]=XYZSum(&fp[1],&tp[12]);
    mcfragsweepCCD3(pChainTemp,&chainFlexTemp,trand,trand2-1,trand2,fp[0],fp[1],fp[2],toldenergy,&tnewenergy);
    flagtor=NewMainChainPos(pChainTemp,&chainFlexTemp,0);
    if(!flagtor) printf("tor2str wrong in bbp\n");
    flagclash=TRUE;
  }while(!flagclash);

  //flagTor=NewMainChainPos(StructureGetChain(&structureTemp,selectChain),&chainFlexTemp,trandpos);
  //if(trandpos==0){   //should rebuild all chains in structure to fix problem with NewMainChainPos
    for(int j=0; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  //}
  //else{
  //  for(int j=trandpos-1; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
  //    ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
  //    ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
  //    if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
  //  }
  //}

  //TODO: optimize interface sidechains
  OptimizeSideChains(&structureTemp,bbrotlib,numSideChainOptCycles,selectChain); 
  EnergyTermInitialize(energyTerms);
  StructureCalcEnergyFlex(&structureTemp,bbrotlib,energyTerms);
  EnergyTermWeighting(energyTerms);
  newenergy=energyTerms[0];
//  newenergy+=calcConstraintEnergy(pChainTemp,natDistRestr)*75;
//  newenergy+=calcDistanceEnergy(pChainTemp,distanceCB,distanceWeight);
  if(newenergy<oldenergy){
    oldenergy=newenergy;
    //nummov[angtype][0]++;
    ChainFlex *origChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlexCopy(origChainFlex,&chainFlexTemp);
    Chain *origChain = StructureGetChain(decstr,selectChain);
    ChainCopy(origChain,StructureGetChain(&structureTemp,selectChain));
    getRamaOutliers(origChain,origChainFlex,ramaduke);
    getAngleOutliers(origChain,origChainFlex);
    getBondLengthOutliers(origChain,origChainFlex);
    getClashes(origChain,origChainFlex);
    StructureDestroy(&structureTemp);
    return TRUE;
  }
  else{
    StructureDestroy(&structureTemp);
    return FALSE;
  }
}

BOOL moveHelixPack(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
                   bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,
                   sssegment *sse,int *alphasind,int numsalpha,pairaa *natDistRestr,double **contactCA,double **contactCB){
  int numSideChainOptCycles=1;
  int i;
  int trandpos;
  double trand2,tlength;
  int leng1,leng2,leng3,indbin;
  double rmat[9];
  XYZ tp[12];
  double tnewenergy;
  BOOL flagclash=FALSE;
  BOOL flagtor;
  const double angleaa[30]={0.01371,0.05788,0.115,0.17441,0.25209,0.3115,0.38842,0.45163,0.52094,0.5773,
                            0.62909,0.68469,0.73343,0.77684,0.80807,0.8431,0.86443,0.89337,0.9086,0.92764,
                            0.93906,0.95582,0.9642,0.97258,0.97791,0.98857,0.99238,0.99466,0.99999,1.0};
  double threshclash=0.01;
  double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};

  if(numsalpha==0){   
    //summcaaa[2]++;
    return FALSE;
  }

  int selectChain = deschn[int(Random()*numchains)];
  ChainFlex* pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);

  //StructureFlex structureFlexTemp;
  //for(i=0;i<StructureFlexGetChainCount(pStructureFlex);i++){
    //ChainFlex *pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlex chainFlexTemp;
    ChainFlexCreate(&chainFlexTemp);
    ChainFlexCopy(&chainFlexTemp,pChainFlex);
    //StructureFlexAddChain(structureFlexTemp, &chainFlex);
  //}

  Structure structureTemp;
  StructureCreate(&structureTemp);
  for(i=0;i<StructureGetChainCount(decstr);i++){
    Chain *pChainSelect = StructureGetChain(decstr,i);
    Chain chainTemp;
    ChainCreate(&chainTemp);
    ChainCopy(&chainTemp,pChainSelect);
    StructureAddChain(&structureTemp, &chainTemp);
    ChainDestroy(&chainTemp);
  }
  Chain *pChainTemp = StructureGetChain(&structureTemp,selectChain);

  do{
    trandpos=int((numsalpha)*Random());
    leng1=sse[alphasind[trandpos]].term-sse[alphasind[trandpos]].init+1;
    leng2=sse[alphasind[trandpos]+1].term-sse[alphasind[trandpos]+1].init+1;
    leng3=sse[alphasind[trandpos]+2].term-sse[alphasind[trandpos]+2].init+1;
    if(leng1>4 && leng3>4 && leng3<20){
      indbin=0;
      trand2=Random();
      for(i=0;i<30;i++){
        if(trand2<=angleaa[i]){
          indbin=i;
          break;
        }
      }
      trand2=(indbin+Random())/30.0*PI;
    }
    else trand2=Random()*PI;
    Residue *pResidue1 = ChainGetResidue(pChainTemp,sse[alphasind[trandpos]].init);
    Residue *pResidue2 = ChainGetResidue(pChainTemp,sse[alphasind[trandpos]].term);
    Residue *pResidue3 = ChainGetResidue(pChainTemp,sse[alphasind[trandpos]+2].init);
    Residue *pResidue4 = ChainGetResidue(pChainTemp,sse[alphasind[trandpos]+2].term);
    Atom *pAtomCA1 = ResidueGetAtomByName(pResidue1, "CA");
    Atom *pAtomCA2 = ResidueGetAtomByName(pResidue2, "CA");
    Atom *pAtomCA3 = ResidueGetAtomByName(pResidue3, "CA");
    Atom *pAtomCA4 = ResidueGetAtomByName(pResidue4, "CA");

    tp[0].X=pAtomCA1->xyz.X;tp[0].Y=pAtomCA1->xyz.Y;tp[0].Z=pAtomCA1->xyz.Z;
    tp[1].X=pAtomCA2->xyz.X;tp[1].Y=pAtomCA2->xyz.Y;tp[1].Z=pAtomCA2->xyz.Z;
    tp[2].X=pAtomCA3->xyz.X;tp[2].Y=pAtomCA3->xyz.Y;tp[2].Z=pAtomCA3->xyz.Z;
    tp[3].X=pAtomCA4->xyz.X;tp[3].Y=pAtomCA4->xyz.Y;tp[3].Z=pAtomCA4->xyz.Z;
    tp[4]=XYZDifference(&tp[2],&tp[3]);
    tlength=XYZNormalization(&tp[4]);
    tp[5]=rana(trand2);
    tp[6].X=0;tp[6].Y=0;tp[6].Z=1;
    tp[7]=XYZDifference(&tp[1],&tp[0]);
    tp[8]=XYZUnit(&tp[7]);
    v2rot(tp[8],rmat);
    tp[9]=rotv(&tp[5],rmat);
    tp[10]=XYZScale2(&tp[9],tlength);
    tp[11]=XYZSum(&tp[2],&tp[10]);//for new end
    mcfragsweepCCD4(pChainTemp,&chainFlexTemp,sse[alphasind[trandpos]+1].init,
                    sse[alphasind[trandpos]+1].term, tp[2],
                    sse[alphasind[trandpos]+2].init,tp[11],
                    sse[alphasind[trandpos]+2].term, oldenergy,&tnewenergy);
    flagtor=NewMainChainPos(pChainTemp,&chainFlexTemp,0);
    if(!flagtor) printf("tor2str wrong in aaa\n");
    flagclash=TRUE;
  }while(!flagclash);

  //flagTor=NewMainChainPos(StructureGetChain(&structureTemp,selectChain),&chainFlexTemp,trandpos);
  //if(ChainGetType(pChainTemp)!=Type_Chain_Protein) continue;
  //if(trandpos==0){   //should rebuild all chains in structure to fix problem with NewMainChainPos
    for(int j=0; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  //}
  //else{
  //  for(int j=trandpos-1; j<ChainGetResidueCount(pChainTemp); j++){
  //    ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
  //    ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
  //    ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
  //    if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
  //  }
  //}

  //TODO: optimize interface sidechains
  OptimizeSideChains(&structureTemp,bbrotlib,numSideChainOptCycles,selectChain); 
  EnergyTermInitialize(energyTerms);
  StructureCalcEnergyFlex(&structureTemp,bbrotlib,energyTerms);
  EnergyTermWeighting(energyTerms);
  newenergy=energyTerms[0];
//  newenergy+=calcConstraintEnergy(pChainTemp,natDistRestr)*75;
//  newenergy+=calcDistanceEnergy(pChainTemp,distanceCB,distanceWeight);
  //std::cout << "New Energy: " << newenergy << " Old Energy: " << oldenergy << std::endl;
  if(newenergy<oldenergy){
    oldenergy=newenergy;
    //nummov[angtype][0]++;
    ChainFlex *origChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlexCopy(origChainFlex,&chainFlexTemp);
    Chain *origChain = StructureGetChain(decstr,selectChain);
    ChainCopy(origChain,StructureGetChain(&structureTemp,selectChain));
    getRamaOutliers(origChain,origChainFlex,ramaduke);
    getAngleOutliers(origChain,origChainFlex);
    getBondLengthOutliers(origChain,origChainFlex);
    getClashes(origChain,origChainFlex);
    StructureDestroy(&structureTemp);
    return TRUE;
  }
  else{
    StructureDestroy(&structureTemp);
    return FALSE;
  }
}

BOOL genesse(ChainFlex *pChainFlex,sssegment *sse,int &numsse){
  int i,j;
  numsse=0;
  for(i=0;i<ChainFlexGetResidueCount(pChainFlex);i++){
    j=i;
    ResidueFlex *pResidueFlexi = ChainFlexGetResidue(pChainFlex,i);
    ResidueFlex *pResidueFlexj = ChainFlexGetResidue(pChainFlex,j);

    while(j<ChainFlexGetResidueCount(pChainFlex) && pResidueFlexj->ssm==pResidueFlexi->ssm){
      j++;
      pResidueFlexj = ChainFlexGetResidue(pChainFlex,j);
    }
    sse[numsse].init=i;
    sse[numsse].term=j-1;
    sse[numsse].ss=pResidueFlexi->ssm;
    numsse++;
    i=j-1;
  }
  if(numsse!=0) sse=(sssegment *)realloc(sse,numsse*sizeof(sssegment));
  return TRUE;
}

BOOL moveHelix(Structure *decstr,StructureFlex *pStructureFlex,int deschn[],int numchains,double &oldenergy,double &newenergy,
               bool flagacc,BBdepRotamerLib *bbrotlib,AtomParamsSet *atomParams,ResiTopoSet *resiTopos,int ***ramaduke,pairaa *natDistRestr,
               double **contactCA,double **contactCB){
  int numSideChainOptCycles=1;
  int i,j;
  static int torideal[][2]={
    148,159, 148,158, 148,159, 148,159, 149,156,
    148,159, 148,158, 148,157, 149,158, 148,158,
    148,158, 148,159, 150,161, 148,159, 148,159,
    148,158, 148,157, 148,157, 149,157, 149,156,

    104, 75, 119, 65, 135, 53, 113, 71, 117, 71,
     90, 90, 121, 66, 119, 64, 113, 71, 126, 62,
    115, 70, 125, 64, 149, 70, 116, 69, 118, 67,
    110, 76, 119, 64, 119, 64, 120, 68, 114, 76,

    148, 72, 148, 67, 136,  0, 150,164, 148, 67,
     41,  0, 148, 68, 145, 64, 147,170, 143, 74,
    148, 70,  25, 21, 149, 73, 135,  0, 149, 68,
    144, 78, 142, 81, 146, 66, 144, 71, 148, 65
  };//h20 e20 c20
  double raddeg = PI/180.0;
  double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM]={0};

  int selectChain = deschn[int(Random()*numchains)];
  ChainFlex* pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
  Chain* pChain = StructureGetChain(decstr,selectChain);

  //StructureFlex structureFlexTemp;
  //for(i=0;i<StructureFlexGetChainCount(pStructureFlex);i++){
    //ChainFlex *pChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlex chainFlexTemp;
    ChainFlexCreate(&chainFlexTemp);
    ChainFlexCopy(&chainFlexTemp,pChainFlex);
    //StructureFlexAddChain(structureFlexTemp, &chainFlex);
  //}

  Structure structureTemp;
  StructureCreate(&structureTemp);
  for(i=0;i<StructureGetChainCount(decstr);i++){
    Chain *pChainSelect = StructureGetChain(decstr,i);
    Chain chainTemp;
    ChainCreate(&chainTemp);
    ChainCopy(&chainTemp,pChainSelect);
    StructureAddChain(&structureTemp, &chainTemp);
    ChainDestroy(&chainTemp);
  }
  Chain *pChainTemp = StructureGetChain(&structureTemp,selectChain);

  sssegment *sse=new sssegment[ChainGetResidueCount(pChain)];  //TODO: change to calc
  int numsse;
  genesse(pChainFlex,sse,numsse); //change calc ssm using hbond energy
  int del=3;
  bool flagdir;
  bool flagtor;
  int tothelix=0;
  for(i=0;i<numsse;i++){
    if(sse[i].ss=='H') tothelix++;
  }
printf("before\n");

  if(tothelix==0){
    //nummov[13][1]++;
    return FALSE;
  }
printf("after\n");

  j=int(tothelix*Random());
  tothelix=0;
  for(i=0;i<numsse;i++){
    if(sse[i].ss=='H' && tothelix==j) break;
    else if(sse[i].ss=='H') tothelix++;
  }
  if(Random()<0.5) flagdir=FALSE;//head
  else flagdir=TRUE;
  if((!flagdir && sse[i].init==0) || (flagdir && sse[i].term==ChainGetResidueCount(pChain)-1)){
    //nummov[13][1]++;
    return FALSE;
  }

  XYZ pt,pn,pc;
  j=-1;
  if(!flagdir){
    tothelix=sse[i].init-1;
    if(tothelix<0){
      //nummov[13][1]++;
      return FALSE;
    }
    if(i>0 && sse[i-1].ss=='C') j=sse[i-1].init;
    ResidueFlex *pResidueFlex1 = ChainFlexGetResidue(&chainFlexTemp,tothelix);
    ResidueFlex *pResidueFlex2 = ChainFlexGetResidue(&chainFlexTemp,tothelix+1);
    Residue *pResidue1 = ChainGetResidue(pChainTemp,tothelix);
    Residue *pResidue2 = ChainGetResidue(pChainTemp,tothelix+1);
    int aaInd1 = getResidueIndex(pResidue1);
    int aaInd2 = getResidueIndex(pResidue2);

    pResidueFlex2->dih[2]=2*torideal[aaInd2][0]+1;  
    pResidueFlex2->dih[1]=180.0;
    pResidueFlex2->dih[0]=2*torideal[aaInd1][1]+1;                

    Atom *pAtomN1 = ResidueGetAtomByName(pResidue1, "N");
    Atom *pAtomCA1 = ResidueGetAtomByName(pResidue1, "CA");
    Atom *pAtomC1 = ResidueGetAtomByName(pResidue1, "C");
    Atom *pAtomN2 = ResidueGetAtomByName(pResidue2, "N");
    Atom *pAtomCA2 = ResidueGetAtomByName(pResidue2, "CA");
    Atom *pAtomC2 = ResidueGetAtomByName(pResidue2, "C");
    flagtor=convertTorToPos(pAtomC2->xyz.X,pAtomC2->xyz.Y,pAtomC2->xyz.Z,
                            pAtomCA2->xyz.X,pAtomCA2->xyz.Y,pAtomCA2->xyz.Z, 
                            pAtomN2->xyz.X,pAtomN2->xyz.Y,pAtomN2->xyz.Z,
                            pResidueFlex2->dih[2]*raddeg,pResidueFlex2->len[0],pResidueFlex2->ang[1]*raddeg,
                            &pn.X,&pn.Y,&pn.Z);
    if(!flagtor){
      printf("wrong %d back c\n",tothelix);
    }
   
    pAtomC1->xyz.X=pn.X;pAtomC1->xyz.Y=pn.Y;pAtomC1->xyz.Z=pn.Z;

    flagtor=convertTorToPos(pAtomCA2->xyz.X,pAtomCA2->xyz.Y,pAtomCA2->xyz.Z,
                            pAtomN2->xyz.X,pAtomN2->xyz.Y,pAtomN2->xyz.Z,
                            pAtomC1->xyz.X,pAtomC1->xyz.Y,pAtomC1->xyz.Z,
                            pResidueFlex2->dih[1]*raddeg,pResidueFlex1->len[2],pResidueFlex2->ang[0]*raddeg,
                            &pt.X,&pt.Y,&pt.Z);
    if(!flagtor){
      printf("wrong %d back ca\n",tothelix);
    }
    pAtomCA1->xyz.X=pt.X;pAtomCA1->xyz.Y=pt.Y;pAtomCA1->xyz.Z=pt.Z;

    flagtor=convertTorToPos(pAtomN2->xyz.X,pAtomN2->xyz.Y,pAtomN2->xyz.Z,
                            pAtomC1->xyz.X,pAtomC1->xyz.Y,pAtomC1->xyz.Z,
                            pAtomCA1->xyz.X,pAtomCA1->xyz.Y,pAtomCA1->xyz.Z,
                            pResidueFlex2->dih[0]*raddeg,pResidueFlex1->len[1],pResidueFlex1->ang[2]*raddeg,
                            &pc.X,&pc.Y,&pc.Z);
    if(!flagtor){
      printf("wrong %d back n\n",tothelix);
    }
   
    pAtomN1->xyz.X=pc.X;pAtomN1->xyz.Y=pc.Y;pAtomN1->xyz.Z=pc.Z;

    //[j,tothelix-1]
    if(tothelix==0){
      str2torp(pChainTemp,&chainFlexTemp,0,1);
      //printf("head0 %d %d\n",0,1);
    }
    else if(tothelix==1){
      str2torp(pChainTemp,&chainFlexTemp,0,2);
      //printf("head1 %d %d\n",0,2);
    }
    else{
      if(j==-1){
        j=tothelix-del;
        if(j<1) j=1;
      }
      else{
        j=tothelix-j;
        if(j>=del) j=int(j*Random());
        j=tothelix-j;
      }
      mcfragsweepLMP2(pChainTemp,j,tothelix);
      //printf("head2 %d %d\n",j-1,tothelix);
      str2torp(pChainTemp,&chainFlexTemp,j-1,tothelix);
    }
  }
  else{
    tothelix=sse[i].term+1;
    if(tothelix>ChainGetResidueCount(pChain)-1){
      //nummov[13][1]++;
      return FALSE;
    }
    if(i<ChainGetResidueCount(pChain)-1 && sse[i+1].ss=='C') j=sse[i+1].term;
    if(j==ChainGetResidueCount(pChain)-1) j--;
    ResidueFlex *pResidueFlex1 = ChainFlexGetResidue(&chainFlexTemp,tothelix-1);
    ResidueFlex *pResidueFlex2 = ChainFlexGetResidue(&chainFlexTemp,tothelix);
    Residue *pResidue1 = ChainGetResidue(pChainTemp,tothelix-1);
    Residue *pResidue2 = ChainGetResidue(pChainTemp,tothelix);
    int aaInd1 = getResidueIndex(pResidue1);
    int aaInd2 = getResidueIndex(pResidue2);
    pResidueFlex2->dih[2]=2*torideal[aaInd2][0]+1;
    pResidueFlex2->dih[1]=180.0;
    pResidueFlex2->dih[0]=2*torideal[aaInd1][1]+1;
    Atom *pAtomN1 = ResidueGetAtomByName(pResidue1, "N");
    Atom *pAtomCA1 = ResidueGetAtomByName(pResidue1, "CA");
    Atom *pAtomC1 = ResidueGetAtomByName(pResidue1, "C");
    Atom *pAtomN2 = ResidueGetAtomByName(pResidue2, "N");
    Atom *pAtomCA2 = ResidueGetAtomByName(pResidue2, "CA");
    Atom *pAtomC2 = ResidueGetAtomByName(pResidue2, "C");

    flagtor=convertTorToPos(pAtomN1->xyz.X,pAtomN1->xyz.Y,pAtomN1->xyz.Z,
                            pAtomCA1->xyz.X,pAtomCA1->xyz.Y,pAtomCA1->xyz.Z, 
                            pAtomC1->xyz.X,pAtomC1->xyz.Y,pAtomC1->xyz.Z, 
                            pResidueFlex2->dih[0]*raddeg,pResidueFlex2->len[0],pResidueFlex2->ang[0]*raddeg,
                            &pn.X,&pn.Y,&pn.Z);
    if(!flagtor){
      printf("wrong %d for n\n",tothelix);
    }
    
    pAtomN2->xyz.X=pn.X;pAtomN2->xyz.Y=pn.Y;pAtomN2->xyz.Z=pn.Z;

    flagtor=convertTorToPos(pAtomCA1->xyz.X,pAtomCA1->xyz.Y,pAtomCA1->xyz.Z, 
                            pAtomC1->xyz.X,pAtomC1->xyz.Y,pAtomC1->xyz.Z, 
                            pAtomN2->xyz.X,pAtomN2->xyz.Y,pAtomN2->xyz.Z, 
                            pResidueFlex2->dih[1]*raddeg,pResidueFlex2->len[1],pResidueFlex2->ang[1]*raddeg,
                            &pt.X,&pt.Y,&pt.Z);
    if(!flagtor){
      printf("wrong %d for ca\n",tothelix);
    }

    pAtomCA2->xyz.X=pt.X;pAtomCA2->xyz.Y=pt.Y;pAtomCA2->xyz.Z=pt.Z;

    flagtor=convertTorToPos(pAtomC1->xyz.X,pAtomC1->xyz.Y,pAtomC1->xyz.Z,
                            pAtomN2->xyz.X,pAtomN2->xyz.Y,pAtomN2->xyz.Z,
                            pAtomCA2->xyz.X,pAtomCA2->xyz.Y,pAtomCA2->xyz.Z,                          
                            pResidueFlex2->dih[2]*raddeg,pResidueFlex2->len[2],pResidueFlex2->ang[2]*raddeg,
                            &pc.X,&pc.Y,&pc.Z);
    if(!flagtor){
      printf("wrong %d for c\n",tothelix);
    }
    
    pAtomC2->xyz.X=pc.X;pAtomC2->xyz.Y=pc.Y;pAtomC2->xyz.Z=pc.Z;

    //[tothelix+1,j]
    if(tothelix==ChainGetResidueCount(pChain)-1){
      str2torp(pChainTemp,&chainFlexTemp,ChainGetResidueCount(pChain)-2,ChainGetResidueCount(pChain)-1);
      //printf("tail0 %d %d\n",numseq-2,numseq-1);
    }
    else if(tothelix==ChainGetResidueCount(pChain)-2){
      str2torp(pChainTemp,&chainFlexTemp,ChainGetResidueCount(pChain)-3,ChainGetResidueCount(pChain)-1);
      //printf("tail1 %d %d\n",numseq-3,numseq-1);
    }
    else{
      if(j==-1){
        j=tothelix+del;
        if(j>ChainGetResidueCount(pChain)-2) j=ChainGetResidueCount(pChain)-2;
      }
      else{
        j=j-tothelix;
        if(j>=del) j=int(j*Random());
        j=tothelix+j;
      }
      mcfragsweepLMP2(pChainTemp,tothelix+1,j+1);
      //printf("tail3 %d %d\n",tothelix,j+1);
      str2torp(pChainTemp,&chainFlexTemp,tothelix,j+1);
    }
  }

  //TODO: Fix when start
  //flagTor=NewMainChainPos(StructureGetChain(&structureTemp,selectChain),&chainFlexTemp,trandpos);
  //if(ChainGetType(pChainTemp)!=Type_Chain_Protein) continue;
  //if(trandpos==0){   //should rebuild all chains in structure to fix problem with NewMainChainPos
    for(j=0; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
      ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
      //ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
      if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
    }
  //}
  //else{
  //  for(int j=trandpos-1; j<ChainGetResidueCount(pChainTemp); j++){
      //ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,deschn,j,pBBdepRotLib,atomParams,resiTopos); //for design positions
  //    ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&structureTemp,selectChain,j,bbrotlib,atomParams,resiTopos); //for repack residues
  //    ProteinSiteAddCrystalRotamer(&structureTemp,selectChain,j,resiTopos);
  //    if(FLAG_EXPAND_HYDROXYL_ROT) ProteinSiteExpandHydroxylRotamers(&structureTemp,selectChain,j,resiTopos);
  //  }
  //}

  //TODO: optimize interface sidechains
  OptimizeSideChains(&structureTemp,bbrotlib,numSideChainOptCycles,selectChain); 
  EnergyTermInitialize(energyTerms);
  StructureCalcEnergyFlex(&structureTemp,bbrotlib,energyTerms);
  EnergyTermWeighting(energyTerms);
  newenergy=energyTerms[0];
//  newenergy+=calcConstraintEnergy(pChainTemp,natDistRestr)*75;
//  newenergy+=calcDistanceEnergy(pChainTemp,distanceCB,distanceWeight);
  //EnergyHbondnhoc(pChainTemp,&chainFlexTemp,1);
  delete[]sse;
  sse=NULL; 
  if(newenergy<oldenergy){
    oldenergy=newenergy;
    //nummov[angtype][0]++;
    ChainFlex *origChainFlex = StructureFlexGetChain(pStructureFlex,selectChain);
    ChainFlexCopy(origChainFlex,&chainFlexTemp);
    Chain *origChain = StructureGetChain(decstr,selectChain);
    ChainCopy(origChain,StructureGetChain(&structureTemp,selectChain));
    getRamaOutliers(origChain,origChainFlex,ramaduke);
    getAngleOutliers(origChain,origChainFlex);
    getBondLengthOutliers(origChain,origChainFlex);
    getClashes(origChain,origChainFlex);
    StructureDestroy(&structureTemp);
    return TRUE;
  }
  else{
    StructureDestroy(&structureTemp);
    return FALSE;
  }
}


