///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "CommonPara.h"
#include "BasicFunc.h"
#include "ParsePDB.h"
#include "SegCluster.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SegCluster::SegCluster()
{
  nosample=0;
  fv=NULL;
  ck=0;
  cfv=NULL;
  rmin=0;
  rmax=1000000;
}

SegCluster::~SegCluster()
{
  int i;
  if(fv)
  {
    delete[]fv;
    fv=NULL;
  }
  if(cfv)
  {
    for(i=0;i<ck;i++)
    {
      delete[]cfv[i].ind2vec;
    }
    delete[]cfv;
    cfv=NULL;
  }
}

double SegCluster::calckmeanmodels(
  int num,
  int nk,
  double **inmat,
  modelinfo *minf,
  bool **flagmodels
)
{
  int i,j;
  double totvar;
  int numiter=0;
  int totchange=0;
  int totcenterchange=0;
  int *indcenter=new int[nk];
  int deli=num/nk;
  int stai=0;
  if(deli==0)
  {
    printf("Not enough samples %d %d\n",num,nk);
    return 0;
  }
  for(i=0;i<nk;i++)
  {
    for(j=0;j<num;j++)
    {
      flagmodels[i][j]=false;
    }
    minf[i].ind=stai;
    stai+=deli;
  }
  calckmeanclass(num,nk,inmat,minf,flagmodels);
  for(i=0;i<nk;i++)
  {
    indcenter[i]=calckmeancenter(num,inmat,flagmodels[i]);
    minf[i].ind=indcenter[i];
  }
  calckmeanvar(num,nk,inmat,minf,flagmodels);
  do
  {
    totchange=calckmeanclass(num,nk,inmat,minf,flagmodels);
    totcenterchange=0;
    for(i=0;i<nk;i++)
    {
      indcenter[i]=calckmeancenter(num,inmat,flagmodels[i]);
      if(indcenter[i]!=minf[i].ind)
      {
        totcenterchange++;
      }
      minf[i].ind=indcenter[i];	
    }
    totvar=calckmeanvar(num,nk,inmat,minf,flagmodels);
    numiter++;
  } while(numiter<200 && totchange>num*0.01);
  delete[]indcenter;
  /////////////////////////////////////////////largest clu in front
  modelinfo tm;
  bool *tmflag=new bool[num];
  for(i=0;i<nk-1;i++)
  {
    for(j=i+1;j<nk;j++)
    {
      if(minf[i].totnei<minf[j].totnei)
      {
        tm=minf[i];
        minf[i]=minf[j];
        minf[j]=tm;
        memcpy(tmflag,flagmodels[i],num*sizeof(bool));
        memcpy(flagmodels[i],flagmodels[j],num*sizeof(bool));
        memcpy(flagmodels[j],tmflag,num*sizeof(bool));
      }
    }
  }
  delete[]tmflag;
	
  return totvar;
}

double SegCluster::calckmeanmodels2(
  int num,
  int nk,
  double **inmat,
  modelinfo *minf,
  bool **flagmodels
)//no initial center 
{
  int i,j;
  double totvar;
  int numiter=0;
  int totchange=0;
  int totcenterchange=0;
  int *indcenter=new int[nk];
  do
  {
    totchange=calckmeanclass(num,nk,inmat,minf,flagmodels);
    totcenterchange=0;
    for(i=0;i<nk;i++)
    {
      indcenter[i]=calckmeancenter(num,inmat,flagmodels[i]);
      if(indcenter[i]!=minf[i].ind)
      {
        totcenterchange++;
      }
      minf[i].ind=indcenter[i];	
    }
    totvar=calckmeanvar(num,nk,inmat,minf,flagmodels);
    numiter++;
  } while(numiter<200 && totchange>num*0.01);
  delete[]indcenter;
  /////////////////////////////////////////////largest clu in front
  modelinfo tm;
  bool *tmflag=new bool[num];
  for(i=0;i<nk-1;i++)
  {
    for(j=i+1;j<nk;j++)
    {
      if(minf[i].totnei<minf[j].totnei)
      {
        tm=minf[i];
        minf[i]=minf[j];
        minf[j]=tm;
        memcpy(tmflag,flagmodels[i],num*sizeof(bool));
        memcpy(flagmodels[i],flagmodels[j],num*sizeof(bool));
        memcpy(flagmodels[j],tmflag,num*sizeof(bool));
      }  
    }
  }
  delete[]tmflag;
	
  return totvar;
}

void SegCluster::calckmeanmodelsiter(
  int num,
  int nk,
  double **inmat,
  modelinfo *minf,
  bool **flagmodels
)
{
  int i,j;
  int ii,jj;
  int numiter=0;	
  int minii,minij,maxiv;
  double totvar;
  int indvar=-1;
  double minvar;
  bool *flagclusterii=new bool[num];
  bool *flagclusterjj=new bool[num];
  minvar=calckmeanmodels(num,nk,inmat,minf,flagmodels);
  calckmeansplit(num,nk,inmat,minf,flagmodels,&minii,&minij,&maxiv);
  if(maxiv>=0)
  {
    calckmeanfarest(num,inmat,flagmodels[maxiv],&ii,&jj,flagclusterii,flagclusterjj);
  }
  while(numiter<20 && maxiv>=0 && inmat[minf[minii].ind][minf[minij].ind]<1.0*(minf[minii].cut+minf[minij].cut))
  {
    for(i=0;i<num;i++)
    {
      if(flagmodels[minij][i])
      {
        flagmodels[minii][i]=true;
      }
    }
    minf[minii].totnei+=minf[minij].totnei;
    j=calckmeancenter(num,inmat,flagmodels[minii]);
    minf[minii].ind=j;

    //minf[minij].ind=ii;
    memcpy(flagmodels[minij],flagclusterii,num*sizeof(bool));
    j=calckmeancenter(num,inmat,flagmodels[minij]);
    minf[minij].ind=j;

    //minf[maxiv].ind=jj;
    memcpy(flagmodels[maxiv],flagclusterjj,num*sizeof(bool));
    j=calckmeancenter(num,inmat,flagmodels[maxiv]);
    minf[maxiv].ind=j;

    totvar=calckmeanmodels2(num,nk,inmat,minf,flagmodels);
    calckmeansplit(num,nk,inmat,minf,flagmodels,&minii,&minij,&maxiv);
    if(maxiv>=0)
    {
      calckmeanfarest(num,inmat,flagmodels[maxiv],&ii,&jj,flagclusterii,flagclusterjj);
    }
    numiter++;  
    if(totvar<minvar)
    {
      minvar=totvar;
      indvar=numiter;
    }
  }
	
  //redo
  calckmeanmodels(num,nk,inmat,minf,flagmodels);
  calckmeansplit(num,nk,inmat,minf,flagmodels,&minii,&minij,&maxiv);
  if(maxiv>=0)
  {
    calckmeanfarest(num,inmat,flagmodels[maxiv],&ii,&jj,flagclusterii,flagclusterjj);
  }
  for(int k=0;k<indvar;k++)
  {
    for(i=0;i<num;i++)
    {
      if(flagmodels[minij][i])
      {
        flagmodels[minii][i]=true;
      }
    }
    minf[minii].totnei+=minf[minij].totnei;
    j=calckmeancenter(num,inmat,flagmodels[minii]);
    minf[minii].ind=j;

    //minf[minij].ind=ii;
    memcpy(flagmodels[minij],flagclusterii,num*sizeof(bool));
    j=calckmeancenter(num,inmat,flagmodels[minij]);
    minf[minij].ind=j;

    //minf[maxiv].ind=jj;
    memcpy(flagmodels[maxiv],flagclusterjj,num*sizeof(bool));
    j=calckmeancenter(num,inmat,flagmodels[maxiv]);
    minf[maxiv].ind=j;

    totvar=calckmeanmodels2(num,nk,inmat,minf,flagmodels);
    calckmeansplit(num,nk,inmat,minf,flagmodels,&minii,&minij,&maxiv);
    if(maxiv>=0)
    {
      calckmeanfarest(num,inmat,flagmodels[maxiv],&ii,&jj,flagclusterii,flagclusterjj);
    }
    if(totvar<minvar)
    {
      minvar=totvar;
      indvar=numiter;
    }
  }
	
  delete[]flagclusterii;
  delete[]flagclusterjj;
}

void SegCluster::calckmeansplit(
  int num,
  int nk,
  double **inmat,
  modelinfo *minf,
  bool **flagmodels,
  int *minii,
  int *minij,
  int *maxiv
)
{
  int i,j;
  double minidist,maxivar;
  double **matcen;
  matcen=new double*[nk];
  for(i=0;i<nk;i++)
  {
    matcen[i]=new double[nk];
  }
  for(i=0;i<nk;i++)
  {
    for(j=0;j<nk;j++)
    {
      matcen[i][j]=inmat[minf[i].ind][minf[j].ind];   
    } 
  }
  minidist=1000000000;
  for(i=0;i<nk;i++)
  {
    for(j=i+1;j<nk;j++)
    {
      if(matcen[i][j]<minidist)
      {
        minidist=matcen[i][j];
        *minii=i;
        *minij=j;
      }
    }
  }
  maxivar=-1000000;
  for(i=0;i<nk;i++)
  {
    if(i!=*minii && i!=*minij)
    {
      if(minf[i].cut>maxivar && minf[i].totnei>1 && minf[i].totnei>0.05*num)
      {
        maxivar=minf[i].cut;
        *maxiv=i;
      }
    }
  }
  if(maxivar<0)
  {
    for(i=0;i<nk;i++)
    {
      if(i!=*minii && i!=*minij)
      {
        if(minf[i].cut>maxivar && minf[i].totnei>1)
        {
          maxivar=minf[i].cut;
          *maxiv=i;
        }
      }
    }
  }
  if(maxivar<0 && minf[0].totnei>1)
  {
    *maxiv=0;
    maxivar=minf[0].cut;
    minidist=1000000000;
    for(i=1;i<nk;i++)
    {
      for(j=i+1;j<nk;j++)
      {
        if(matcen[i][j]<minidist)
        {
          minidist=matcen[i][j];
          *minii=i;
          *minij=j;
        }
      }
    }
  }
  if(maxivar<0 || minidist>10000000)
  {
    printf("No cluster for split %d %d\n",num,nk);
    *maxiv=-1;
  }
  for(i=0;i<nk;i++)
  {
    delete[]matcen[i];
  }
  delete[]matcen;
}

void SegCluster::calckmeanfarest(
  int num,
  double **inmat,
  bool *flagcluster,
  int *ii,
  int *jj,
  bool *flagclusterii,
  bool *flagclusterjj
)
{
  int i,j;
  double tdist=-1000000;
  for(i=0;i<num;i++)
  {	
    if(flagcluster[i])
    {
      for(j=i+1;j<num;j++)
      {
        if(flagcluster[j])
        {
          if(tdist<inmat[i][j])
          {
            tdist=inmat[i][j];
            *ii=i;
            *jj=j;
          }
        }
      }
    }
  }
  for(i=0;i<num;i++)
  {
    flagclusterii[i]=false;
    flagclusterjj[i]=false;
    if(flagcluster[i])
    {
      if(inmat[*ii][i]<inmat[*jj][i])
      {
        flagclusterii[i]=true;
      }
      else
      {
        flagclusterjj[i]=true;
      }
      if(i==*ii)
      {
        flagclusterii[i]=true;
        flagclusterjj[i]=false;
      }
      else if(i==*jj)
      {
        flagclusterii[i]=false;
        flagclusterjj[i]=true;
      }
    }
  }
}

double SegCluster::calckmeanvar(
  int num,
  int nk,
  double **inmat,
  modelinfo *minf,
  bool **flagmodels
)
{
  int i,j;
  double totvar=0;
  double tvar;
  int tnum;
  for(i=0;i<nk;i++)
  {
    tvar=0;
    tnum=0;
    if(minf[i].totnei==1)
    {
      minf[i].cut=0;
      continue;
    }
    for(j=0;j<num;j++)
    {
      if(flagmodels[i][j])
      {
        tnum++;
        tvar+=inmat[minf[i].ind][j]*inmat[minf[i].ind][j];
      }
    }
    totvar+=tvar;
    tvar/=double(tnum-1);
    minf[i].cut=sqrt(tvar);
  }
  totvar/=double(num-nk);
  totvar=sqrt(totvar);
	
  return totvar;
}

int SegCluster::calckmeanclass(
  int num,
  int nk,
  double **inmat,
  modelinfo *minf,
  bool **flagmodels
)//know center calc neighbors totchange
{
  int i,j;
  int totchange=0;
  double tmin;
  int indmin;
  for(i=0;i<nk;i++)
  {
    minf[i].totnei=0;
  }
  for(i=0;i<num;i++)
  {
    tmin=1000000000;
    for(j=0;j<nk;j++)
    {
      if(inmat[i][minf[j].ind]<tmin)
      {
        tmin=inmat[i][minf[j].ind];
        indmin=j;
      }
      if(i==minf[j].ind)
      {
        indmin=j;
        break;
      }
    }
    if(!flagmodels[indmin][i])
    {
      for(j=0;j<nk;j++)
      {
        flagmodels[j][i]=false;
      }
      flagmodels[indmin][i]=true;
      totchange++;
    }
    minf[indmin].totnei++;
  }
	
  return totchange;
}

int SegCluster::calckmeancenter(
  int num,
  double **inmat,
  bool *flagcluster
)//each cluster calc new center index
{
  int i,j,indmin;
  double ttot;
  double tmin=1000000000000;
  for(i=0;i<num;i++)
  {	
    if(flagcluster[i])
    {
      ttot=0;
      for(j=0;j<num;j++)
      {
        if(flagcluster[j])
        {
          ttot+=inmat[i][j]*inmat[i][j]*inmat[i][j];
        }
      }
      if(ttot<tmin)
      {
        indmin=i;
        tmin=ttot;
      }
    }
  }
	
  return indmin;
}

int SegCluster::calcgaumodels(
  int num,
  double **inmat,
  double sigma,
  int *vind,
  double probmin,
  double probmax,
  double deltadist,
  int nummodel,
  modelinfo *minf,
  bool **flagmodels
)
{
  int i,j;
  bool *flagexclude=new bool[num];
  bool *flagcluster=new bool[num];
  for(i=0;i<num;i++)
  {
    flagexclude[i]=true;
  }
  int totremain=num;
  int totneigh,indmod;
  double tmpdist=sigma;
  for(i=0;i<nummodel;i++)
  {
    tmpdist=sigma;
    if(totremain<5)
    {
      break;
    }
    for(j=0;j<num;j++)
    {
      if(flagexclude[vind[j]])
      {
        indmod=vind[j];
        break;
      }
    }
    int iternum=0;
    do
    {
      totneigh=calcnumcluster(num,indmod,inmat,tmpdist,flagexclude,flagcluster);
      if(totneigh/double(totremain)>probmax)
      {
        tmpdist-=deltadist;
      }
      else if(totneigh/double(totremain)<probmin)
      {
        tmpdist+=deltadist;
      }
      iternum++;
      if(iternum>1000)
      {
        break;
      }
    } while(totneigh/double(totremain)<probmin || totneigh/double(totremain)>probmax);
    if(iternum>1000)
    {
      break;
    }
    minf[i].ind=indmod;
    minf[i].totnei=totneigh;
    minf[i].cut=tmpdist;
    memcpy(flagmodels[i],flagcluster,num*sizeof(bool));
    totremain-=totneigh;
    for(j=0;j<num;j++)
    {
      if(flagcluster[j])
      {
        flagexclude[j]=false;
      }
    }
  }
  delete[]flagexclude;
  delete[]flagcluster;
	
  return i;
}

int SegCluster::calcmodels(
  int num, 
  double **inmat,
  double probmin,
  double probmax,
  double initdist,
  double deltadist,
  int nummodel,
  modelinfo *minf,
  bool **flagmodels
)
{
  int i,j;
  double ddeltadist=deltadist;
  bool *flagexclude=new bool[num];
  bool *flagcluster=new bool[num];
  for(i=0;i<num;i++)
  {
    flagexclude[i]=true;
  }
  double tmpdist=initdist;
  bool flagsl[2];
  int totremain=num;
  int totneigh,indmod;
  for(i=0;i<nummodel;i++)
  {
    tmpdist=initdist;
    if(totremain<5)
    {
      break;
    }
    int iternum=0;
    do
    {
      flagsl[1]=flagsl[0];
      indmod=calcmaxneighbor(num,inmat,tmpdist,flagexclude,&totneigh,flagcluster);
      if(totneigh/double(totremain)>probmax && tmpdist>rmin)
      {
        flagsl[0]=false;
        if(iternum>1000)
        {
          if(flagsl[0]!=flagsl[1])//deltadist is big
          {
            ddeltadist/=1.5;
          }
          else//deltadist is too small, drop too slow
          {
            ddeltadist*=1.2;
          }
          iternum=0;
        }
        tmpdist-=ddeltadist;
        if(tmpdist<rmin) tmpdist=rmin;
      }
      else if(totneigh/double(totremain)<probmin && tmpdist<rmax)
      {
        flagsl[0]=true;
        if(iternum>1000)
        {
          if(flagsl[0]!=flagsl[1])//deltadist is big
          {
            ddeltadist/=1.5;
          }
          else//deltadist is too small, rise too slow
          {
            ddeltadist*=1.2;
          }
          iternum=0;
        }
        tmpdist+=ddeltadist;
        if(tmpdist>rmax) tmpdist=rmax;
      }
      else
      {
        break;
      }
      iternum++;
      if(iternum>1000)
      {
        break;
      }
      if(ddeltadist<epsilon || ddeltadist>10000000)
      {
        break;
      }
    } while(totneigh/double(totremain)<probmin || totneigh/double(totremain)>probmax);
    minf[i].ind=indmod;
    minf[i].totnei=totneigh;
    minf[i].cut=tmpdist;
    fprintf(stderr,"cluster %2d size %7d cutoff %10.5f ratio %.3f niter %5d remain %7d\n",
            i, minf[i].totnei, minf[i].cut,totneigh/double(totremain),iternum,totremain);
    memcpy(flagmodels[i],flagcluster,num*sizeof(bool));
    totremain-=totneigh;
    for(j=0;j<num;j++)
    {
      if(flagcluster[j])
      {
        flagexclude[j]=false;
      }
    }
  }
  delete[]flagexclude;
  delete[]flagcluster;
	
  return i;
}

int SegCluster::calcmaxneighbor(
  int num, 
  double **inmat,
  double distcut,
  bool *flagexclude,
  int *totneigh,
  bool *flagcluster
)
{
  int i;
  int maxneigbor=-1;
  int indmax;
  int tmpnum;
  for(i=0;i<num;i++) if(flagexclude[i])
  {
    tmpnum=calcnumcluster(num,i,inmat,distcut,flagexclude,flagcluster);
    if(tmpnum>maxneigbor)
    {
      indmax=i;
      maxneigbor=tmpnum;
      *totneigh=maxneigbor;
    }
  }
  calcnumcluster(num,indmax,inmat,distcut,flagexclude,flagcluster);
	
  return indmax;
}

//return which one
int SegCluster::calcmaxneighbor(
  int num,
  char *distname,
  double distcut,
  bool *flagexclude,
  int *totneigh,
  bool *flagcluster
)
{
  int i;
  int maxneigbor=-1;
  int indmax;
  int tmpnum;
  for(i=0;i<num;i++) if(flagexclude[i])
  {
    tmpnum=calcnumcluster(num,i,distname,distcut,flagexclude,flagcluster);
    if(tmpnum>maxneigbor)
    {
      indmax=i;
      maxneigbor=tmpnum;
      *totneigh=maxneigbor;
    }
  }
  calcnumcluster(num,indmax,distname,distcut,flagexclude,flagcluster);

  return indmax;
}

int SegCluster::calcnumcluster(
  int num,
  int indi,
  double **inmat,
  double distcut,
  bool *flagexclude,
  bool *flagcluster
)
{
  int tot=0;
  int i;
  for(i=0;i<num;i++)
  {
    flagcluster[i]=false;
    if(flagexclude[i] && inmat[indi][i]<=distcut)
    {
      tot++;
      flagcluster[i]=true;
    }
  }

  return tot;
}

int SegCluster::calcnumcluster(
  int num,
  int indi,
  char *distname,
  double distcut,
  bool *flagexclude,
  bool *flagcluster
)
{
  int tot=0;
  int i;
  FILE *file;
  float tdist;
  char tmpname[150];
  sprintf(tmpname,"%s%d.txt",distname,indi);
  file=fopen(tmpname,"rt");
  for(i=0;i<num;i++)
  {
    flagcluster[i]=false;
    fgets(tmpname,150,file);
    sscanf(tmpname,"%f",&tdist);
    if(tdist<=distcut && flagexclude[i])
    {
      tot++;
      flagcluster[i]=true;
    }
  }
  fclose(file);

  return tot;
}

void SegCluster::calcgaussian(
  int num,
  double **inmat,
  double sigma,
  double *outvect
)
{
  int i,j;
  double tdist;
  BasicFunc bf;
  for(i=0;i<num;i++)
  {
    outvect[i]=1;
    for(j=i+1;j<num;j++)
    {
      tdist=bf.fungaussian(inmat[i][j],1,sigma,0);
      outvect[i]+=tdist;
      outvect[j]+=tdist;
    }
  }
}

void SegCluster::calcfragdistmat(
  int indi,
  char *fragname,
  double **outmat
)
{
  int i,j,k;
  FILE *file,*file2;
  point3f ptmp;
  double tdist,pmat[9],ptrans[3];
  ParsePDB pp;
  char namefrag[200];
  double pin[3*nosegunit],pout[3*nosegunit];
  sprintf(namefrag,"%s%d.topse",fragname,indi);
  file=fopen(namefrag,"rt");
  if(file==NULL)
  {
    printf("not found frag file %s\n",namefrag);
    return;
  }
  for(i=0;i<topno;i++)
  {
    fgets(namefrag,200,file);
    for(j=0;j<nosegunit;j++)
    {
      fgets(namefrag,200,file);
      sscanf(namefrag,"%c %f %f %f %c %f %f %f %f %f %f %f %f %f %f %f %f",
             &ptmp.residueid,&ptmp.x,&ptmp.y,&ptmp.z,&ptmp.stype,&ptmp.phi,
             &ptmp.leng,&ptmp.angl,&ptmp.tor[0],&ptmp.len[0],&ptmp.ang[0],
             &ptmp.tor[1],&ptmp.len[1],&ptmp.ang[1],&ptmp.tor[2],&ptmp.len[2],
             &ptmp.ang[2]);
      pin[j]=ptmp.x;
      pin[nosegunit+j]=ptmp.y;
      pin[2*nosegunit+j]=ptmp.z;
    }
    sprintf(namefrag,"%s%d.topse",fragname,indi);
    file2=fopen(namefrag,"rt");
    for(j=0;j<topno;j++)
    {		
      fgets(namefrag,200,file2);
      for(k=0;k<nosegunit;k++)
      {
        fgets(namefrag,200,file2);
        sscanf(namefrag,"%c %f %f %f %c %f %f %f %f %f %f %f %f %f %f %f %f",
               &ptmp.residueid,&ptmp.x,&ptmp.y,&ptmp.z,&ptmp.stype,&ptmp.phi,
               &ptmp.leng,&ptmp.angl,&ptmp.tor[0],&ptmp.len[0],&ptmp.ang[0],
               &ptmp.tor[1],&ptmp.len[1],&ptmp.ang[1],&ptmp.tor[2],&ptmp.len[2],
               &ptmp.ang[2]);
        pout[k]=ptmp.x;
        pout[nosegunit+k]=ptmp.y;
        pout[2*nosegunit+k]=ptmp.z;
      }	
      tdist=pp.lsfrmsd(pin,pout,nosegunit,pmat,ptrans);
      outmat[i][j]=tdist;
      outmat[j][i]=tdist;   
    }
    fclose(file2);
  }
  fclose(file);
}

void SegCluster::calcdistmat(
  int num,
  dihedral *intor,
  double **outmat
)
{
  int i,j;
  double tdist,delx,dely,tx,ty;
  for(i=0;i<num;i++)
  {
    outmat[i][i]=0;
    for(j=i+1;j<num;j++)
    {
      delx=fabs(intor[i].phi-intor[j].phi);
      tx=fabs(intor[i].phi-intor[j].phi+360);
      if(delx>tx)
      {
        delx=tx;
      }
      tx=fabs(intor[i].phi-intor[j].phi-360);
      if(delx>tx)
      {
        delx=tx;
      }
      dely=fabs(intor[i].psi-intor[j].psi);
      ty=fabs(intor[i].psi-intor[j].psi+360);
      if(dely>ty)
      {
        dely=ty;
      }
      ty=fabs(intor[i].psi-intor[j].psi-360);
      if(dely>ty)
      {
        dely=ty;
      }
      tdist=sqrt(delx*delx+dely*dely);
      outmat[i][j]=tdist;
      outmat[j][i]=tdist;
    }
  }
}

bool SegCluster::extractone(
  char *filename,
  int indnum,
  char *outname
)
{
  int j;
  FILE *file,*fileout;
  char oneline[300];
  file=fopen(filename,"rt");
  if(!file)
  {
    printf("not found remc %s\n",filename);
    return false;
  }
  do
  {
    fgets(oneline,300,file);
    j=-10000;
    if(!(oneline[0]=='M' && oneline[1]=='O' && oneline[2]=='D' && 
         oneline[3]=='E' && oneline[4]=='L'))
    {
      continue;
    }
    sscanf(oneline+9,"%d",&j);
  } while(j!=indnum && !feof(file));
  fileout=fopen(outname,"wt");
  fprintf(fileout,"%s",oneline);
  while(!feof(file) && !(oneline[0]=='E' && oneline[1]=='N' && 
        oneline[2]=='D' && oneline[3]=='M'))
  {
    fgets(oneline,300,file);
    if(strlen(oneline)>14 && oneline[13]=='S' && oneline[14]=='G') continue;
    fprintf(fileout,"%s",oneline);
  }
  fclose(fileout);
  fclose(file);

  return true;
}

int SegCluster::extractremc(
  int istart,
  int iend,
  char *prename,
  int nummc,
  char *outname
)
{
  int i,k;
  FILE *file,*fileout;
  char pdbname[200];
  int indmc=0;
  for(i=istart;i<=iend;i++)
  {
    sprintf(pdbname,"%s%d.pdb",prename,i);
    file=fopen(pdbname,"rt");
    if(!file)
    {
      printf("not found remc %s\n",pdbname);
      return -1;
    }
    for(k=0;k<nummc;k++)
    {
      sprintf(pdbname,"%s%d.pdb",outname,indmc);
      printf("%s\n",pdbname);
      fileout=fopen(pdbname,"wt");
      do
      {
        fgets(pdbname,200,file);
        fprintf(fileout,"%s",pdbname);
      } while(!(pdbname[0]=='E' && pdbname[1]=='N' && pdbname[2]=='D' && 
              pdbname[3]=='M' && pdbname[4]=='D') && !feof(file));
      if(feof(file) && !(pdbname[0]=='E' && pdbname[1]=='N' && pdbname[2]=='D' && 
         pdbname[3]=='M' && pdbname[4]=='D')) 
      {
        fclose(fileout);
        break;
      }
      fclose(fileout);
      indmc++;
    }
    fclose(file);
  }

  return indmc;
}

int SegCluster::extractremc(
  char *filename,
  int indstart,
  char *outname
)
{
  FILE *file,*fileout;
  int indmc=0;
  char oneline[300];
  file=fopen(filename,"rt");
  if(!file)
  {
    printf("not found remc %s\n",filename);
    return indmc;
  }
  while(!feof(file))
  {
    sprintf(oneline,"%s%d.pdb",outname,indstart+indmc);
    fileout=fopen(oneline,"wt");
    do
    {
      fgets(oneline,300,file);
      fprintf(fileout,"%s",oneline);
    } while(!(oneline[0]=='E' && oneline[1]=='N' && oneline[2]=='D' && oneline[3]=='M' && 
            oneline[4]=='D') && !feof(file));
    if(feof(file) && !(oneline[0]=='E' && oneline[1]=='N' && oneline[2]=='D' && 
       oneline[3]=='M' && oneline[4]=='D')) 
    {
      fclose(fileout);
      break;
    }
    fclose(fileout);
    indmc++;
  }
  fclose(file);

  return indmc;
}

int SegCluster::extractremc(
  int istart,
  int iend,
  int numseq,
  char *prename,
  int nummc,
  char *outname
)
{
  int i,k;
  FILE *file,*fileout;
  char pdbname[200];
  int indmc=0;
  for(i=istart;i<=iend;i++)
  {
    sprintf(pdbname,"%s/%d.pdb",prename,i);
    file=fopen(pdbname,"rt");
    if(!file)  
    {
      printf("not found remc %s\n",pdbname);
      return -1;
    }
    for(k=0;k<nummc;k++)
    {
      sprintf(pdbname,"%s/%d.pdb",outname,indmc);
      fileout=fopen(pdbname,"wt");
      do
      {
        fgets(pdbname,200,file);
        fprintf(fileout,"%s",pdbname);
      } while(!(pdbname[0]=='E' && pdbname[1]=='N' && pdbname[2]=='D' && pdbname[3]=='M' && 
              pdbname[4]=='D')&& !feof(file));
      if(feof(file) && !(pdbname[0]=='E' && pdbname[1]=='N' && pdbname[2]=='D' && 
         pdbname[3]=='M' && pdbname[4]=='D')) 
      {
        fclose(fileout);
        break;
      }
      fclose(fileout);
      indmc++; 
    }
    fclose(file);
  }

  return indmc;
}

void SegCluster::outputdist(
  int num,
  int indi,
  char *proname,
  double sigma,
  char *outname,
  double *gauval
)
{
  ParsePDB pp1,pp2;
  int i,j;
  double tdist;
  FILE *file;
  char tmpname[150];
  double totval=0;
  BasicFunc bf;
  sprintf(tmpname,"%s.%d",proname,indi+1);
  pp1.loadpdb(tmpname);
  pp1.extractbb(0,-1,1);
  double pmat[9],ptrans[3];
  double *pout=new double[pp1.numbb*3];
  double *pin=new double[pp1.numbb*3];
  for(i=0;i<pp1.numbb;i++)
  {
    if(pp1.bb[i].indca==-1) pp1.bb[i].indca=0;
    pin[i]=pp1.proseq[pp1.bb[i].indca].x;
    pin[i+pp1.numbb]=pp1.proseq[pp1.bb[i].indca].y;
    pin[i+2*pp1.numbb]=pp1.proseq[pp1.bb[i].indca].z;
  }
  sprintf(tmpname,"%s%d.txt",outname,indi);
  file=fopen(tmpname,"wt");
  for(j=0;j<num;j++)
  {
    if(j==indi)
    {
      fprintf(file,"%f\n",0.0);
      continue;
    }
    sprintf(tmpname,"%s.%d",proname,j+1);
    pp2.loadpdb(tmpname);
    pp2.extractbb(0,-1,1);	
    for(i=0;i<pp1.numbb;i++)
    {
      if(pp2.bb[i].indca==-1) pp2.bb[i].indca=0;
      pout[i]=pp2.proseq[pp2.bb[i].indca].x;
      pout[i+pp1.numbb]=pp2.proseq[pp2.bb[i].indca].y;
      pout[i+2*pp1.numbb]=pp2.proseq[pp2.bb[i].indca].z;
    }
    tdist=pp1.lsfrmsd(pin,pout,pp1.numbb,pmat,ptrans);
    tdist+=fabs(dscore[indi]-dscore[j]);
    fprintf(file,"%f\n",tdist);
    totval+=bf.fungaussian(tdist,1,sigma,0);
  }
  fclose(file);
  delete[]pin;
  delete[]pout;
  *gauval=totval;
}

double SegCluster::eucdist(
  int p, 
  int c
)
{  
  double dist;                         
  int i;                                 
  dist=0;
  for(i=0;i<vecleng;i++)
  {
    dist += (cfv[c].f[i]-fv[p].f[i])*(cfv[c].f[i]-fv[p].f[i]);
  } 

  return dist;
}

int SegCluster::closestcluster(
  int pat
)
{
  int i,clustid;
  double mindist,d;
  mindist=100000000;
  clustid=-1;
  for(i=0;i<ck;i++) 
  {
    d=eucdist(pat,i);
    if(d<mindist)
    {
      mindist=d;
      clustid=i;
    } 
  } 

  return clustid;
}

void SegCluster::distrisamples()
{
  int i,pat,clustid,memberindex;
  for(i=0;i<ck;i++)
  {
    cfv[i].num=0;
  }
  for(pat=0;pat<nosample;pat++) 
  {
    clustid=closestcluster(pat);
    memberindex=cfv[clustid].num;
    cfv[clustid].ind2vec[memberindex]=pat;
    cfv[clustid].num++;
  } 
}

bool SegCluster::calcnewcenters()
{
  int vectid,i,j,k;
  bool convflag;
  double tmp[vecleng];
  convflag=true;
  for(i=0;i<ck;i++)
  {           
    for(j=0;j<vecleng;j++) 
    {           
      tmp[j]=0.0;
    } 
    for(j=0;j<cfv[i].num;j++) 
    { 
      vectid=cfv[i].ind2vec[j];
      for(k=0;k<vecleng;k++)
      {        
        tmp[k]+=fv[vectid].f[k];    
      }	 
    } 
    for(k=0;k<vecleng;k++)
    { 
      tmp[k]=tmp[k]/cfv[i].num;
      if(tmp[k]!=cfv[i].f[k])
      {
        convflag=false;
      }
      cfv[i].f[k]=tmp[k];
    } 	  
  } 

  return convflag;
}

void SegCluster::runkmeans()
{
  bool converged;
  numiter=0;
  converged=false;
  while(converged==false) 
  {
    distrisamples();
    converged=calcnewcenters();
    numiter++;
  }
}
