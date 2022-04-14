///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "ParsePDB.h"
#include "BasicFunc.h"
#include "SegCluster.h"
#include "Fragments.h"

Fragments::Fragments()
{
  int i;
  namepdb=NULL;
  homolist=NULL;
  numhomo=0;
  numsel=0;
  dhcennum=NULL;
  for(i=0;i<20;i++)
  {
     fragcont[i]=NULL;
  }
  for(i=0;i<nosegdh;i++)
  {
    fragdh[i]=NULL;
  }
  oldseglength=0;
}


Fragments::~Fragments()
{
  int i,j;
  if(namepdb)
  {
    delete[]namepdb;
    namepdb=NULL;
  }
  if(homolist)
  {
    delete[]homolist;
    homolist=NULL;
  }

  if(oldseglength>0)
  {
    if(fragcont[0])
    {
      for(j=0;j<oldseglength*topno;j++)
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
  for(i=0;i<nosegdh;i++)
  {
    if(fragdh[i]) delete[]fragdh[i];
    fragdh[i]=NULL;
  }
  if(dhcennum)
  {
    delete[]dhcennum;
    dhcennum=NULL;
  }

}

bool Fragments::loadPISCES(char *selfile)
{
  FILE *file;
  file= fopen(selfile, "r");
  if(!file)
  {
    printf("Error when loading PISCES file %s\n",selfile);
    return false;
  }
  int allocnumsel=500;
  numsel=0;
  if(namepdb)
  {
    delete[]namepdb;
  }
  namepdb=new namep[allocnumsel];
  char oneline[300];
  fgets(oneline,300,file);
  while(!feof(file))
  {
    fgets(oneline,300,file);
    sscanf(oneline,"%s %d",namepdb[numsel].name,&namepdb[numsel].seqnum);
    sscanf(oneline+27,"%f %f %f",&namepdb[numsel].res,&namepdb[numsel].rfac,
           &namepdb[numsel].rval);
    numsel++;
    if(numsel>=allocnumsel)
    {
      allocnumsel*=2;
      namepdb=(namep *)realloc(namepdb,allocnumsel*sizeof(namep));
    }
  }
  fclose(file);
  numsel--;
  if(numsel!=0)
  namepdb=(namep *)realloc(namepdb,numsel*sizeof(namep));
  return true;
}

void Fragments::calcfragdh(
  char *dataDir,
  int seqlength,
  char *outname,
  int segleng,
  int cuttop
)
{
  int i,j,k;
  char fragFileName[STD_FILE_NAME_LENGTH+1];
  if(!fragcont[0])
  {
    sprintf(fragFileName,"%s/%dseqfra.topse",dataDir,segleng);
    loadfragd(seqlength,segleng,cuttop,fragFileName);
  }
  for(i=0;i<nosegdh;i++)
  {
    if(fragdh[i]) delete[]fragdh[i];
    fragdh[i]=new dihedral[seqlength];
  }
  if(dhcennum)
  {
    delete[]dhcennum;
    dhcennum=new int[seqlength];
  }
  else
  {
    dhcennum=new int[seqlength];
  }

  FILE *file;
  ParsePDB pp;
  SegCluster sc;
  dihedral **dh=new dihedral*[seqlength];
  int *dhnum=new int[seqlength];
  int istart;
  for(i=0;i<seqlength;i++)
  {
    dhnum[i]=0;
    dh[i]=new dihedral[cuttop*segleng];
  }

  for(i=0;i<seqlength-segleng+1;i++)
  {
    for(j=0;j<cuttop;j++)
    {
      istart=j*segleng;
      for(k=0;k<segleng;k++)
      {
        if(fragcont[0][istart+k][i].tor[0]>=0 &&
           fragcont[0][istart+k][i].tor[1]>=0 && fragcont[0][istart+k][i].tor[2]>=0)
        {
          strcpy(dh[i+k][dhnum[i+k]].name,fragcont[0][istart+k][i].name);
          dh[i+k][dhnum[i+k]].ind=fragcont[0][istart+k][i].resind;
          dh[i+k][dhnum[i+k]].ss=fragcont[0][istart+k][i].stype;
          dh[i+k][dhnum[i+k]].aa=fragcont[0][istart+k][i].residueid;
          dh[i+k][dhnum[i+k]].tp[0].x=fragcont[0][istart+k][i].x;
          dh[i+k][dhnum[i+k]].tp[0].y=fragcont[0][istart+k][i].y;
          dh[i+k][dhnum[i+k]].tp[0].z=fragcont[0][istart+k][i].z;
          dh[i+k][dhnum[i+k]].caphi=fragcont[0][istart+k][i].phi;
          dh[i+k][dhnum[i+k]].calen=fragcont[0][istart+k][i].leng;
          dh[i+k][dhnum[i+k]].cainn=fragcont[0][istart+k][i].angl;
          dh[i+k][dhnum[i+k]].psi=fragcont[0][istart+k][i].tor[0];
          dh[i+k][dhnum[i+k]].leng[0]=fragcont[0][istart+k][i].len[0];
          dh[i+k][dhnum[i+k]].inner[0]=fragcont[0][istart+k][i].ang[0];
          dh[i+k][dhnum[i+k]].omega=fragcont[0][istart+k][i].tor[1];
          dh[i+k][dhnum[i+k]].leng[1]=fragcont[0][istart+k][i].len[1];
          dh[i+k][dhnum[i+k]].inner[1]=fragcont[0][istart+k][i].ang[1];
          dh[i+k][dhnum[i+k]].phi=fragcont[0][istart+k][i].tor[2];
          dh[i+k][dhnum[i+k]].leng[2]=fragcont[0][istart+k][i].len[2];
          dh[i+k][dhnum[i+k]].inner[2]=fragcont[0][istart+k][i].ang[2];
          dhnum[i+k]++;
        }
      }//k
    }//j
  }//i

  for(j=0;j<seqlength;j++)
  {
    for(i=0;i<nosegdh;i++)
    {
      strcpy(fragdh[i][j].name,"UNKWO");
      fragdh[i][j].ss='C';
      fragdh[i][j].aa='X';
      fragdh[i][j].ind=0;
      fragdh[i][j].caphi=180;
      fragdh[i][j].calen=lencaca;
      fragdh[i][j].cainn=180;
      fragdh[i][j].phi=180;
      fragdh[i][j].psi=180;
      fragdh[i][j].omega=180;
      fragdh[i][j].leng[0]=lencn;
      fragdh[i][j].leng[1]=lennca;
      fragdh[i][j].leng[2]=lencac;
      fragdh[i][j].inner[0]=angcacn;
      fragdh[i][j].inner[1]=angcnca;
      fragdh[i][j].inner[2]=angncac;
      fragdh[i][j].tp[0].x=0;fragdh[i][j].tp[0].y=0;fragdh[i][j].tp[0].z=0;
      fragdh[i][j].tp[1].x=0;fragdh[i][j].tp[1].y=0;fragdh[i][j].tp[1].z=0;
      fragdh[i][j].tp[2].x=0;fragdh[i][j].tp[2].y=0;fragdh[i][j].tp[2].z=0;
      fragdh[i][j].tp[3].x=0;fragdh[i][j].tp[3].y=0;fragdh[i][j].tp[3].z=0;
      fragdh[i][j].clunum=0;
    }
    if(dhnum[j]<=nosegdh)
    {
      printf("few data for clustering residue %d %d\n",j,dhnum[j]);
      dhcennum[j]=dhnum[j];
      for(i=0;i<dhnum[j];i++)
      {
        fragdh[i][j]=dh[j][i];
        fragdh[i][j].clunum=(i+1.0)/double(dhcennum[j]);
      }
      for(i=0;i<dhcennum[j]-1;i++)
      {
        for(k=i+1;k<dhcennum[j];k++)
        {
          if(strcmp(fragdh[i][j].name,fragdh[k][j].name)==0 && fragdh[i][j].ind==fragdh[k][j].ind)
          {
            if(k<dhcennum[j]-1)
            {
              fragdh[k][j]=fragdh[dhcennum[j]-1][j];
              dhcennum[j]--;
              k--;
            }
            else
            {
              dhcennum[j]--;
            }
          }
        }
      }
      continue;
    }

    double **inmat=new double*[dhnum[j]];
    int nummodel=nosegdh;
    int actnum,actnum2,actnum3;
    bool **flagmodels=new bool*[nummodel];
    modelinfo *minf=new modelinfo[nummodel];
    modelinfo *allminf=new modelinfo[3*nummodel];
    modelinfo tpminf;
    int numcluster=0;
    for(i=0;i<dhnum[j];i++)
    {
      inmat[i]=new double[dhnum[j]];
    }
    sc.calcdistmat(dhnum[j],dh[j],inmat);
    double *vect=new double[dhnum[j]];
    int *vind=new int[dhnum[j]];
    BasicFunc bf;
    for(i=0;i<dhnum[j];i++)
    {
      vind[i]=i;
    }
    for(i=0;i<nummodel;i++)
    {
      flagmodels[i]=new bool[dhnum[j]];
    }
    sc.calcgaussian(dhnum[j],inmat,100,vect);
    bf.partbub3(vect,dhnum[j],vind,dhnum[j]);
    actnum=sc.calcgaumodels(dhnum[j],inmat,60,vind,0.35,0.85,1.0,nummodel,minf,flagmodels);
    for(i=0;i<actnum;i++)
    {
      allminf[numcluster]=minf[i];
      numcluster++;
    }
    delete[]vect;
    delete[]vind;
    actnum2=sc.calcmodels(dhnum[j],inmat,0.35,0.85,5,1.0,nummodel,minf,flagmodels);
    for(i=0;i<actnum2;i++)
    {
      allminf[numcluster]=minf[i];
      numcluster++;
    }
    actnum3=nosegdh/2-3;
    if(dhnum[j]<actnum3) actnum3=dhnum[j];
    sc.calckmeanmodelsiter(dhnum[j],actnum3,inmat,minf,flagmodels);//good for psi
    for(i=0;i<actnum3;i++)
    {
      allminf[numcluster]=minf[i];
      numcluster++;
    }
    for(i=0;i<numcluster-1;i++)
    {
      for(k=i+1;k<numcluster;k++)
      {
        if((allminf[i].ind==allminf[k].ind) ||
           (strcmp(dh[j][allminf[i].ind].name,dh[j][allminf[k].ind].name)==0 && 
            dh[j][allminf[i].ind].ind==dh[j][allminf[k].ind].ind))
        {
          int numax=allminf[i].totnei;
          if(numax<allminf[k].totnei) numax=allminf[k].totnei;
          allminf[i].totnei=numax;
          if(k<numcluster-1)
          {
            allminf[k]=allminf[numcluster-1];
            numcluster--;
            k--;
          }
          else
          {
            numcluster--;
          }
        }
      }
    }
    for(i=0;i<numcluster-1;i++)
    {
      for(k=i+1;k<numcluster;k++)
      {
        if(allminf[i].totnei<allminf[k].totnei)
        {
          tpminf=allminf[i];
          allminf[i]=allminf[k];
          allminf[k]=tpminf;
        }
      }
    }
    if(numcluster>nosegdh) numcluster=nosegdh;
    double totnum=0.0;
    for(i=0;i<numcluster;i++)
    {
      fragdh[i][j]=dh[j][allminf[i].ind];
      fragdh[i][j].clunum=allminf[i].totnei;
      totnum+=allminf[i].totnei;
    }
    i=0;
    fragdh[i][j].clunum/=totnum;
    for(i=1;i<numcluster;i++)
    {
      fragdh[i][j].clunum/=totnum;
      fragdh[i][j].clunum+=fragdh[i-1][j].clunum;
    }
    dhcennum[j]=numcluster;
    for(i=0;i<dhnum[j];i++)
    {
      delete[]inmat[i];
    }
    delete[]inmat;
    delete[]minf;
    delete[]allminf;
    for(i=0;i<nummodel;i++)
    {
      delete[]flagmodels[i];
    }
    delete[]flagmodels;
  }
  for(i=0;i<seqlength;i++)
  {
    delete[]dh[i];
  }
  delete[]dh;
  delete[]dhnum;
  file=fopen(outname,"wt");
  for(int l=0;l<seqlength;l++)
  {
    fprintf(file,"%d %d\n",l,dhcennum[l]);
    for(j=0;j<dhcennum[l];j++)
    {
      fprintf(file,"%c %8.3f %8.3f %8.3f %c %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %4d %s %5.3f\n",
              fragdh[j][l].aa,fragdh[j][l].tp[0].x,fragdh[j][l].tp[0].y,fragdh[j][l].tp[0].z,fragdh[j][l].ss,
              fragdh[j][l].caphi,fragdh[j][l].calen,fragdh[j][l].cainn,fragdh[j][l].psi,fragdh[j][l].leng[0],
              fragdh[j][l].inner[0],fragdh[j][l].omega,fragdh[j][l].leng[1],fragdh[j][l].inner[1],
              fragdh[j][l].phi,fragdh[j][l].leng[2],fragdh[j][l].inner[2],fragdh[j][l].ind,fragdh[j][l].name,fragdh[j][l].clunum);
    }
  }
  fclose(file);
}


bool Fragments::loadfragd(
  int numseq,
  int seglength,
  int ntop,
  char *filename
)
{
  int j,k,l,tmpind,istart;
  FILE *file;
  char tmpname[300],tmpstr[8];

  file=fopen(filename,"rt");
  if(!file)
  {
    printf("cannot load fragments %s\n",filename);
    return false;
  }
        
  if(fragcont[0])
  {
    for(j=0;j<oldseglength*topno;j++)
    {
      delete[]fragcont[0][j];
    }
    delete[]fragcont[0];
  }
  fragcont[0]=new point3f*[seglength*topno];
  oldseglength=seglength;
  for(j=0;j<seglength*topno;j++)
  {
    fragcont[0][j]=new point3f[numseq-seglength+1];
  }
               
  for(j=0;j<numseq-seglength+1;j++)
  {
    for(l=0;l<ntop;l++)
    {
      istart=l*seglength;
      fgets(tmpname,300,file);
      sscanf(tmpname,"%s %d %f",tmpstr,&tmpind,&fragcont[0][istart][j].ptsg.x);
      for(k=0;k<seglength;k++)
      {
        fgets(tmpname,300,file);
        sscanf(tmpname,"%c %f %f %f %c %f %f %f %f %f %f %f %f %f %f %f %f",
               &fragcont[0][istart+k][j].residueid,&fragcont[0][istart+k][j].x,&fragcont[0][istart+k][j].y,&fragcont[0][istart+k][j].z,
               &fragcont[0][istart+k][j].stype,&fragcont[0][istart+k][j].phi,&fragcont[0][istart+k][j].leng,&fragcont[0][istart+k][j].angl,
               &fragcont[0][istart+k][j].tor[0],&fragcont[0][istart+k][j].len[0],&fragcont[0][istart+k][j].ang[0],
               &fragcont[0][istart+k][j].tor[1],&fragcont[0][istart+k][j].len[1],&fragcont[0][istart+k][j].ang[1],
               &fragcont[0][istart+k][j].tor[2],&fragcont[0][istart+k][j].len[2],&fragcont[0][istart+k][j].ang[2]);
        fragcont[0][istart+k][j].resind=tmpind+k;
        strcpy(fragcont[0][istart+k][j].name,tmpstr);
      }
      if(j==numseq-seglength && l==ntop-1)
      {
        printf("loading fragments %2d %d %d %d %f %f %f %s\n",0,j,l,seglength,
               fragcont[0][istart+seglength-1][j].tor[2],fragcont[0][istart+seglength-1][j].len[2],
               fragcont[0][istart+seglength-1][j].ang[2],fragcont[0][istart+seglength-1][j].name);
      }
    }               
  }
  fclose(file);
        
  return true;
}

void Fragments::scoreFragments(
  seginfo **seginfall,
  int seqlength,
  InputData inputInfo,
  char *dbfile,
  int segleng,
  int topleng,
  double *wtterm
)
{
  int i,j,k,l;
  int seqLengthQuery;
  bool flagfile;
  char namestr[150];
  ssef  *qssef=new ssef[seqlength];
  ssef  *ss8=new ssef[seqlength];

  ss8=inputInfo.get8StateSS();

  memcpy(qssef,ss8,seqlength*sizeof(ssef));

  char *tpseqdata=new char[seqlength];
  for(i=0;i<seqlength;i++)
  {
    tpseqdata[i]=qssef[i].res;
  }
  int snum;
  int ii;
  int *tpsegnum=new int[seqlength-segleng+1];
  for(ii=0;ii<seqlength-segleng+1;ii++)
  {
    tpsegnum[ii]=0;
  }
  int blocksize=2008;
  int blocknum=numsel/blocksize;
  if(numsel%blocksize!=0) blocknum++;
  int inds,inde,jj,kk;
  float tmpdist,tdist;
  spdf  *sspdf[blocksize];
  spdf  *sspdf2[blocksize];
  spdf  *sspdf3[blocksize];
  ssef  *sssef[blocksize];
  double *sdatexp[blocksize];
  double *sdatphi[blocksize];
  double *sdatpsi[blocksize];
  spdf  *strpdf[blocksize];

  bool flaghomo;
  inds=0;
  inde=inds+blocksize;
  if(inde>numsel) inde=numsel;
  FILE *filedb=fopen(dbfile,"rb");
  if(!filedb)
  {
    fprintf(stderr,"no dbfile\n");
    return;
  }

  for(jj=0;jj<blocknum;jj++)
  {

    printf("doing block %4d/%4d\n",inde,numsel);
    for(i=inds;i<inde;i++)
    {
      sspdf[i-inds]=new spdf[namepdb[i].seqnum];
      sspdf2[i-inds]=new spdf[namepdb[i].seqnum];
      sspdf3[i-inds]=new spdf[namepdb[i].seqnum];
      sssef[i-inds]=new ssef[namepdb[i].seqnum];
      sdatexp[i-inds]=new double[namepdb[i].seqnum];
      sdatphi[i-inds]=new double[namepdb[i].seqnum];
      sdatpsi[i-inds]=new double[namepdb[i].seqnum];
      strpdf[i-inds]=new spdf[namepdb[i].seqnum];
    }

    for(i=inds;i<inde;i++)
    {
      fread(sspdf[i-inds],sizeof(spdf),namepdb[i].seqnum,filedb);//chk
      fread(sspdf2[i-inds],sizeof(spdf),namepdb[i].seqnum,filedb);//mtx
      fread(sspdf3[i-inds],sizeof(spdf),namepdb[i].seqnum,filedb);//pro
      fread(sssef[i-inds],sizeof(ssef),namepdb[i].seqnum,filedb);
      fread(sdatexp[i-inds],sizeof(double),namepdb[i].seqnum,filedb);
      fread(sdatphi[i-inds],sizeof(double),namepdb[i].seqnum,filedb);
      fread(sdatpsi[i-inds],sizeof(double),namepdb[i].seqnum,filedb);
      fread(strpdf[i-inds],sizeof(spdf),namepdb[i].seqnum,filedb);
    }
        
    for(i=inds;i<inde;i++)
    {
      flaghomo=false;
      for(j=0;j<numhomo;j++)
      {
        if(i==homolist[j])
        {
          flaghomo=true;
          break;
        }
      }
      if(flaghomo) continue;
      snum=namepdb[i].seqnum;
      for(j=0;j<snum-segleng+1;j++)
      {
        for(kk=0;kk<seqlength-segleng+1;kk++)
        {
          tmpdist=0;
          for(k=0;k<segleng;k++)
          {
            if(qssef[kk+k].ss=='X' || qssef[kk+k].ss==sssef[i-inds][j+k].ss) continue;
            if(qssef[kk+k].ss=='c' && sssef[i-inds][j+k].ss=='C') continue;
            if(qssef[kk+k].ss=='e' && sssef[i-inds][j+k].ss=='E') continue;
            if(qssef[kk+k].ss=='h' && sssef[i-inds][j+k].ss=='H') continue;
            if(qssef[kk+k].ss=='C' && (sssef[i-inds][j+k].ss=='T' || sssef[i-inds][j+k].ss=='S' ||
               sssef[i-inds][j+k].ss=='C')) continue;
            if(qssef[kk+k].ss=='E' && (sssef[i-inds][j+k].ss=='E' || sssef[i-inds][j+k].ss=='B')) continue;
            if(qssef[kk+k].ss=='H' && (sssef[i-inds][j+k].ss=='H' || sssef[i-inds][j+k].ss=='G'
               || sssef[i-inds][j+k].ss=='I')) continue;

            if((qssef[kk+k].ss=='G' || qssef[kk+k].ss=='I' || qssef[kk+k].ss=='h' || 
               qssef[kk+k].ss=='H') && (sssef[i-inds][j+k].ss=='E' || sssef[i-inds][j+k].ss=='B'))
            {
              tmpdist+=4.0;
              continue;
            }
            else if((qssef[kk+k].ss=='B' || qssef[kk+k].ss=='e' || qssef[kk+k].ss=='E') && 
                    (sssef[i-inds][j+k].ss=='H' || sssef[i-inds][j+k].ss=='G' || sssef[i-inds][j+k].ss=='I'))
            {
              tmpdist+=4.0;
              continue;
            }
            else if(qssef[kk+k].ss=='G' && (sssef[i-inds][j+k].ss=='H' || sssef[i-inds][j+k].ss=='I'))
            {
              tmpdist+=0.50;
              continue;
            }
            else if(qssef[kk+k].ss=='I' && (sssef[i-inds][j+k].ss=='G' || sssef[i-inds][j+k].ss=='H'))
            {
              tmpdist+=0.50;
              continue;
            }
            else if(qssef[kk+k].ss=='h' && (sssef[i-inds][j+k].ss=='G' || sssef[i-inds][j+k].ss=='I'))
            {
              tmpdist+=0.50;
              continue;
            }
            else if((qssef[kk+k].ss=='B' && sssef[i-inds][j+k].ss=='E') || 
                    (qssef[kk+k].ss=='e' && sssef[i-inds][j+k].ss=='B'))
            {
              tmpdist+=0.25;
              continue;
            }
            else if(qssef[kk+k].ss=='T' && (sssef[i-inds][j+k].ss=='S' ||
                    sssef[i-inds][j+k].ss=='C'))
            {
              tmpdist+=0.25;
              continue;
            }
            else if(qssef[kk+k].ss=='S' && (sssef[i-inds][j+k].ss=='T' ||
                    sssef[i-inds][j+k].ss=='C'))
            {
              tmpdist+=0.25;
              continue;
            }
            else if(qssef[kk+k].ss=='c' && (sssef[i-inds][j+k].ss=='T' ||
                    sssef[i-inds][j+k].ss=='S'))
            {
              tmpdist+=0.25;
              continue;
            }
            else if(qssef[kk+k].ss=='T' && (sssef[i-inds][j+k].ss=='H' || sssef[i-inds][j+k].ss=='G' || 
                    sssef[i-inds][j+k].ss=='I'))
            {
              tmpdist+=1.0;
              continue;
            }
            else
            {
              tmpdist+=2.0;
              continue;
            }
          }
          if(tpsegnum[kk]==0)
          {
            strcpy(seginfall[kk][tpsegnum[kk]].name,namepdb[i].name);
            seginfall[kk][tpsegnum[kk]].fdist=tmpdist;
            seginfall[kk][tpsegnum[kk]].start=j;
            seginfall[kk][tpsegnum[kk]].indtmp=i;
            for(ii=0;ii<segleng;ii++)
            {
              seginfall[kk][tpsegnum[kk]].css[ii]=sssef[i-inds][j+ii].ss;
              seginfall[kk][tpsegnum[kk]].caa[ii]=sssef[i-inds][j+ii].res;
              seginfall[kk][tpsegnum[kk]].solv[ii]=sdatexp[i-inds][j+ii];
            }
            tpsegnum[kk]++;
          }
          else if(tpsegnum[kk]==topleng)
          {
            //compare new one with the topno
            if(tmpdist>=seginfall[kk][topleng-1].fdist)
            {

            }
            else
            {
              k=topleng-1;
              while(k>=0 && tmpdist<seginfall[kk][k].fdist)
              {
                k--;
              }
              k++; //in pos k
              for(l=topleng-2;l>=k;l--)
              {
                seginfall[kk][l+1]=seginfall[kk][l];    
              }
              strcpy(seginfall[kk][k].name,namepdb[i].name);
              seginfall[kk][k].fdist=tmpdist;
              seginfall[kk][k].start=j;
              seginfall[kk][k].indtmp=i;
              for(ii=0;ii<segleng;ii++)
              {
                seginfall[kk][k].css[ii]=sssef[i-inds][j+ii].ss;
                seginfall[kk][k].caa[ii]=sssef[i-inds][j+ii].res;
                seginfall[kk][k].solv[ii]=sdatexp[i-inds][j+ii];
              }
            }
          }
          else if(tpsegnum[kk]<topleng)
          {
            //compare and insert
            k=tpsegnum[kk]-1;
            while(k>=0 && tmpdist<seginfall[kk][k].fdist)
            {
              k--;
            }
            k++; //in pos k
            for(l=tpsegnum[kk]-1;l>=k;l--)
            {
              seginfall[kk][l+1]=seginfall[kk][l];
            }
            strcpy(seginfall[kk][k].name,namepdb[i].name);
            seginfall[kk][k].fdist=tmpdist;
            seginfall[kk][k].start=j;
            seginfall[kk][k].indtmp=i;
            for(ii=0;ii<segleng;ii++)
            {
              seginfall[kk][k].css[ii]=sssef[i-inds][j+ii].ss;
              seginfall[kk][k].caa[ii]=sssef[i-inds][j+ii].res;
              seginfall[kk][k].solv[ii]=sdatexp[i-inds][j+ii];
            }
            tpsegnum[kk]++;
          }
        }//each pos in query
      }//each pos in temp
    }//each template
    for(i=inds;i<inde;i++)
    {
      delete[]sspdf[i-inds];
      delete[]sspdf2[i-inds];
      delete[]sspdf3[i-inds];
      delete[]sssef[i-inds];
      delete[]sdatexp[i-inds];
      delete[]sdatphi[i-inds];
      delete[]sdatpsi[i-inds];
      delete[]strpdf[i-inds];
    }
    inds+=blocksize;
    inde=inds+blocksize;
    if(inde>numsel) inde=numsel;
  }
  fclose(filedb);
  delete[]tpsegnum;
  delete[]qssef;
}

bool Fragments::generateFragments(
  char *dataDir,
  char *listName,
  InputData inputInfo,
  char *mtxDbName,
  char *dbFeatName,
  char *outName, 
  int fragSize
)//mkfragmentsdes
{
  int i,j,k,m;
  char solvname[600];
  FILE *file;
  int qleng;
  int topnum=200;
  int seqLength;

  bool flagps=loadPISCES(listName);
  if(!flagps)
  {
    printf("PISCES can not be opened %s\n",listName);
    return false;
  }

  seqLength=inputInfo.getSeqLength();
  vector<double> fragsolv(seqLength,0.0);
  vector<double> fragsolvavg(seqLength,0.0);
  vector<int> fragcount(seqLength,0);
  double wtterm[]={0.30,2.00,0.50,0.40,0.00};//min best
  seginfo **seginfall=new seginfo*[seqLength];
  for(j=0;j<seqLength;j++)
  {
    seginfall[j]=new seginfo[topno];
  }

  m=fragSize;
  printf("begin segments for %d %d\n",m,seqLength-m+1);
  scoreFragments(seginfall,seqLength,inputInfo,mtxDbName,m,topnum,wtterm);
  printf("read feat db\n");
  dihedral **allfeat=new dihedral*[numsel];
  for(i=0;i<numsel;i++)
  {
    allfeat[i]=new dihedral[namepdb[i].seqnum];
  }
  file=fopen(dbFeatName,"rb");
  if(!file)
  {
    printf("no feat db file %s\n",dbFeatName);
    for(i=0;i<numsel;i++)
    {
      delete[]allfeat[i];
    }
    delete[]allfeat;
    for(j=0;j<seqLength;j++)
    {
      delete[]seginfall[j];
    }
    delete[]seginfall;
    return false;
  }
  for(i=0;i<numsel;i++)
  {
    fread(allfeat[i],sizeof(dihedral),namepdb[i].seqnum,file);
  }
  fclose(file);
  file= fopen(outName,"wt");
  if(!file)
  {
    printf("topse can not be opened %s\n",outName);
    for(i=0;i<numsel;i++)
    {
      delete[]allfeat[i];
    }
    delete[]allfeat;
    for(j=0;j<seqLength;j++)
    {
      delete[]seginfall[j];
    }
    delete[]seginfall;
    return false;
  }
  for(i=0;i<seqLength;i++) 
  {
    int index = m-1;
    for(j=i-m+1;j<i+1;j++) 
    {
      if(j>=0 && j<seqLength-m+1) 
      {
        for(int l=0;l<topnum;l++) 
        {
          fragsolv[j+index]+=seginfall[j][l].solv[index];
          fragcount[j+index]+=1;
        }
      }
      index--;
    }
  }

  sprintf(solvname,"%s/%d_solv.txt",dataDir,m);
  ofstream solvfile;
  solvfile.open(solvname);
  solvfile << seqLength << endl;
  for(i=0;i<seqLength;i++) 
  {
    fragsolvavg[i]=fragsolv[i]/fragcount[i];
    solvfile << i << " " << fragsolvavg[i] << endl;
  }
  solvfile.close();

  int iii;
  for(j=0;j<seqLength-m+1;j++)
  {
    for(i=0;i<topnum;i++)
    {
      fprintf(file,"%c%c%c%c%c %4d %9.6f %4d %3d\n",seginfall[j][i].name[0],
              seginfall[j][i].name[1],seginfall[j][i].name[2],seginfall[j][i].name[3],
              seginfall[j][i].name[4],seginfall[j][i].start,seginfall[j][i].fdist,j,i);

      //get coordinate
      iii=seginfall[j][i].indtmp;
      for(k=seginfall[j][i].start;k<seginfall[j][i].start+m;k++)
      {
        fprintf(file,"%c %8.3f %8.3f %8.3f %c ",seginfall[j][i].caa[k-seginfall[j][i].start],
                allfeat[iii][k].tp[0].x,allfeat[iii][k].tp[0].y,allfeat[iii][k].tp[0].z,
                seginfall[j][i].css[k-seginfall[j][i].start]);
        //phi for ca length inner phi psi omega
        if(k<3)
        {
          fprintf(file,"%8.3f %8.3f %8.3f ",-360.0,-1.0,-360.0);
        }
        else
        {
          fprintf(file,"%8.3f %8.3f %8.3f ",allfeat[iii][k-2].caphi,allfeat[iii][k-2].calen,
                  allfeat[iii][k-2].cainn);
        }
        if(k<1)
        {
          fprintf(file,"%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",-360.0,-1.0,
                  -360.0,-360.0,-1.0,-360.0,allfeat[iii][k].phi,allfeat[iii][k].leng[2],
                  allfeat[iii][k].inner[2]);
        }
        else
        {
          fprintf(file,"%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                  allfeat[iii][k-1].psi,allfeat[iii][k-1].leng[0],allfeat[iii][k-1].inner[0],
                  allfeat[iii][k-1].omega,allfeat[iii][k-1].leng[1],allfeat[iii][k-1].inner[1],
                  allfeat[iii][k].phi,allfeat[iii][k].leng[2],allfeat[iii][k].inner[2]);
        }
      }//k
    }//i
  }//j
  fclose(file);
  
  for(i=0;i<numsel;i++)
  {
    delete[]allfeat[i];
  }
  delete[]allfeat;
  allfeat=NULL;
  for(j=0;j<seqLength;j++)
  {
    delete[]seginfall[j];
  }
  delete[]seginfall; 
  seginfall=NULL;

  return true;
}

