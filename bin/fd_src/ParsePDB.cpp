///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "ParsePDB.h"
#include "BasicFunc.h"
#include "PrintFunc.h"

////////////////////////////////////////////////////////////////////////////////
// Construction/Destruction
////////////////////////////////////////////////////////////////////////////////

ParsePDB::ParsePDB()
{
  headnumchain=-1;
  indexchain=0;
  int i;
  for(i=0;i<13;i++)
  {
    cancbins[i]=NULL;
  }

  for(i=0;i<maxnumchains;i++)
  {	
    aminoinchain[i]=NULL;
  }
  pt.ismtrix=false;
  pt.isorigx=false;
  pt.isscale=false;
  helixl=NULL;
  sheetl=NULL;
  turnl=NULL;
  linkl=NULL;
  numhelix=0;numsheet=0;numturn=0;numlink=0;
  allochelix=30;allocsheet=30;allocturn=30;alloclink=30;
  ssbondl=NULL;sltbrgl=NULL;cispepl=NULL;hydbndl=NULL;sitel=NULL;
  numssbond=0;numsltbrg=0;numcispep=0;numhydbnd=0;numsite=0;
  allocssbond=30;allocsltbrg=30;alloccispep=30;allochydbnd=30;allocsite=30;

  numpromod=0;
  allocpromod=1;
  promod=NULL;
  flagfirstmodel=true;
  flagbeginmodel=true;
  numproseq=0;
  allocproseq=800;
  proseq=NULL;
  flagter=false;
  numbb=0;
  bb=NULL;

  numsse=0;
  sses=NULL;
  numhssp=0;
  hssps=NULL;
  hbond=NULL;
  cacont=NULL;
  bsheet=NULL;
  numbsheet=0;
  ahelix=NULL;
  numahelix=0;
  bstrand=NULL;
  numbstrand=0;
  abss=NULL;
  numabss=0;
  diheang=NULL;
  numtor=0;

  bbid=NULL;
  dsspid=NULL;
  cp=NULL;
  numcon=0;
  edget.nedges=0;
  edget.edges=NULL;
  numbone=0;
  pbone=NULL;
  bsp=NULL;
  btmat=NULL;
  btscore=NULL;
  etraj=NULL;
  numtraj=0;
  numcont[0]=0;
  numcont[1]=0;
  segcont[0]=NULL;
  segcont[1]=NULL;
  resnumber=0;
  resmap=NULL;
  hinfo.dddd[0]=-1;hinfo.dddd[1]=-1;hinfo.dddd[2]=-1;
  hinfo.drev[0]=-1;hinfo.drev[1]=-1;hinfo.drev[2]=-1;
  hinfo.dnew[0]=-1;hinfo.dnew[1]=-1;hinfo.dnew[2]=-1;
  hinfo.fobs=false;
  hinfo.frev=false;
  hinfo.species=-1;
  hinfo.iept=-1;
  strcpy(hinfo.name,"XXXX");
  hinfo.fres=-1;
  cca=NULL;
  afd=NULL;
  distmat=NULL;
}

ParsePDB::~ParsePDB()
{
  int i,j,k;
  for(i=0;i<13;i++)
  {
    if(cancbins[i])
    {
      for(j=0;j<105;j++)
      {
        for(k=0;k<105;k++)
          delete[]cancbins[i][j][k];
        delete[]cancbins[i][j];
      }
      delete[]cancbins[i];
    }
  }

  for(i=0;i<maxnumchains;i++)
  {
    if(aminoinchain[i]!=NULL)
    {	
      free(aminoinchain[i]);
      aminoinchain[i]=NULL;
    }
  }
  if(helixl!=NULL)
  {
    free(helixl);
    helixl=NULL;
  }
  if(sheetl!=NULL)
  {
    free(sheetl);
    sheetl=NULL;
  }
  if(turnl!=NULL)
  {
    free(turnl);
    turnl=NULL;
  }
  if(linkl!=NULL)
  {
    free(linkl);
    linkl=NULL;
  }
	
  if(ssbondl!=NULL)
  {
    free(ssbondl);
    ssbondl=NULL;
  }
  if(sltbrgl!=NULL)
  {
    free(sltbrgl);
    sltbrgl=NULL;
  }
  if(cispepl!=NULL)
  {
    free(cispepl);
    cispepl=NULL;
  }
  if(hydbndl!=NULL)
  {
    free(hydbndl);
    hydbndl=NULL;
  }
  if(sitel!=NULL)
  {
    free(sitel);
    sitel=NULL;
  }
  if(promod!=NULL)
  {
    free(promod);
    promod=NULL;
  }
  if(proseq!=NULL)
  {
    free(proseq);
    proseq=NULL;
  }
  if(bb!=NULL)
  {
    delete[]bb;
    bb=NULL;
  }
  if(sses)
  {
    delete[]sses;
    sses=NULL;
  }
  if(hssps)
  {
    delete[]hssps;
    hssps=NULL;
  }
  if(hbond)
  {
    delete[]hbond;
    hbond=NULL;
  }
  if(cacont)
  {
    delete[]cacont;
    cacont=NULL;
  }
  if(bsheet)
  {
    delete[]bsheet;
    bsheet=NULL;
  }
  if(ahelix)
  {
    delete[]ahelix;
    ahelix=NULL;
  }
  if(bstrand)
  {
    delete[]bstrand;
    bstrand=NULL;
  }
  if(abss)
  {
    delete[]abss;
    abss=NULL;
  }
  if(diheang)
  {
    delete[]diheang;
    diheang=NULL;
  }
  if(bbid)
  {
    delete[]bbid;
    bbid=NULL;
  }
  if(dsspid)
  {
    delete[]dsspid;
    dsspid=NULL;
  }
  if(cp)
  {
    delete[]cp;
    cp=NULL;
  }
  if(edget.edges!=NULL)
  {
    free(edget.edges);
    edget.nedges=0;
    edget.edges=NULL;
  }
  if(pbone)
  {
    delete[]pbone;
    pbone=NULL;
  }
  if(bsp)
  {
    delete[]bsp;
    bsp=NULL;
  }
  if(btmat)
  {
    delete[]btmat;
    btmat=NULL;
  }
  if(btscore)
  {
    delete[]btscore;
    btscore=NULL;
  }
  if(etraj)
  {
    delete[]etraj;
    etraj=NULL;
  }
  numtraj=0;
  if(segcont[0])
  {
    delete[]segcont[0];
    segcont[0]=NULL;
  }
  numcont[0]=0;
  if(segcont[1])
  {
    delete[]segcont[1];
    segcont[1]=NULL;
  }
  numcont[1]=0;
  if(resmap)
  {
    delete[]resmap;
    resmap=NULL;
  }
  resnumber=0;
  if(cca)
  {
    delete[]cca;
    cca=NULL;
  }
  if(afd)
  {
    delete[]afd;
    afd=NULL;
  }
  if(distmat)
  {
    delete[]distmat;
    distmat=NULL;
  }
}


/* read chain content 'pdb_lines' into 'decstr_vec'. Only the most
 * basic coordinate information are parsed. Note that secondary
 * structures (ss2, ssm, stype) are all initialized to 'C'.
 * 'decstr_vec' must be pre-initialized to chain size */
void ParsePDB::loadPDBChain(
  const vector<string> &pdb_lines,
  point3f *decstr_vec
)
{
  string prev_resi="";
  string resi="";
  string atom_name="    ";
  string resn="UNK";
  float x,y,z;
  int idx=-1; // residue index, starting from 0
  int a; // amino acid type

  for (int i=0;i<pdb_lines.size();i++)
  {
    atom_name=pdb_lines[i].substr(12,4);
    resn=pdb_lines[i].substr(17,3);
    x=atof(pdb_lines[i].substr(30,8).c_str());
    y=atof(pdb_lines[i].substr(38,8).c_str());
    z=atof(pdb_lines[i].substr(46,8).c_str());
    resi=pdb_lines[i].substr(22,5);
	
    if (resi!=prev_resi) // new residue
    {
      idx+=1;
      prev_resi=resi;

      decstr_vec[idx].aaa  ='X';
      for (a=0;a<26;a++)
      {
        if (resn==aad3[a])
        {
          decstr_vec[idx].aaa=aad1[a];
          break;
        }
      }
      decstr_vec[idx].stype='H';
      decstr_vec[idx].ss3  ='C';
      decstr_vec[idx].ssm  ='C';
      decstr_vec[idx].resind=atoi(resi.c_str());
    }

    if (atom_name==" CA ")
    {
      decstr_vec[idx].x=x;
      decstr_vec[idx].y=y;
      decstr_vec[idx].z=z;
    }
    else if (atom_name==" N  ")
    {
      decstr_vec[idx].ptn.x=x;
      decstr_vec[idx].ptn.y=y;
      decstr_vec[idx].ptn.z=z;
    }
    else if (atom_name==" C  ")
    {
      decstr_vec[idx].ptc.x=x;
      decstr_vec[idx].ptc.y=y;
      decstr_vec[idx].ptc.z=z;
    }
    else if (atom_name==" H  ")
    {
      decstr_vec[idx].pth.x=x;
      decstr_vec[idx].pth.y=y;
      decstr_vec[idx].pth.z=z;
    }
    else if (atom_name==" O  ")
    {
      decstr_vec[idx].pto.x=x;
      decstr_vec[idx].pto.y=y;
      decstr_vec[idx].pto.z=z;
    }
    else if (atom_name==" CB ")
    {
      decstr_vec[idx].ptb.x=x;
      decstr_vec[idx].ptb.y=y;
      decstr_vec[idx].ptb.z=z;
    }
    else if (atom_name==" SC ")
    {
      decstr_vec[idx].ptsg.x=x;
      decstr_vec[idx].ptsg.y=y;
      decstr_vec[idx].ptsg.z=z;
    }
  }

  /* clean up */
  prev_resi.clear();
  resi.clear();
  atom_name.clear();
}

/* read chain content 'pdb_lines' into 'decstr_vec'. Only the most
 * basic coordinate information are parsed. Note that secondary
 * structures (ss2, ssm, stype) are all initialized to 'C'.
 * 'decstr_vec' must be pre-initialized to chain size */
void ParsePDB::loadPDBChainFull(
  const vector<string>&pdb_lines,
  point3ffull *decstr_vec)
{
  string prev_resi="";
  string resi="";
  string atom_name="    ";
  string resn="UNK";
  float x,y,z;
  int idx=-1; // residue index, starting from 0
  int a; // amino acid type
  int count;
  for (int i=0;i<pdb_lines.size();i++)
  {
    atom_name=pdb_lines[i].substr(12,4);
    resn=pdb_lines[i].substr(17,3);
    x=atof(pdb_lines[i].substr(30,8).c_str());
    y=atof(pdb_lines[i].substr(38,8).c_str());
    z=atof(pdb_lines[i].substr(46,8).c_str());
    resi=pdb_lines[i].substr(22,5);
	
    if (resi!=prev_resi) // new residue
    {
      idx+=1;
      prev_resi=resi;
      count=0;
      decstr_vec[idx].aaa  ='X';
      for (a=0;a<26;a++)
      {
        if (resn==aad3[a])
        {
          decstr_vec[idx].aaa=aad1[a];
          break;
        }
      }
      decstr_vec[idx].stype='H';
      decstr_vec[idx].ss3  ='C';
      decstr_vec[idx].ssm  ='C';
      decstr_vec[idx].resind=atoi(resi.c_str());
    }

    if (atom_name==" CA ")
    {
      decstr_vec[idx].x=x;
      decstr_vec[idx].y=y;
      decstr_vec[idx].z=z;
    }
    else if (atom_name==" N  ")
    {
      decstr_vec[idx].ptn.x=x;
      decstr_vec[idx].ptn.y=y;
      decstr_vec[idx].ptn.z=z;
    }
    else if (atom_name==" C  ")
    {
      decstr_vec[idx].ptc.x=x;
      decstr_vec[idx].ptc.y=y;
      decstr_vec[idx].ptc.z=z;
    }
    else if (atom_name==" H  ")
    {
      decstr_vec[idx].pth.x=x;
      decstr_vec[idx].pth.y=y;
      decstr_vec[idx].pth.z=z;
    }
    else if (atom_name==" O  ")
    {
      decstr_vec[idx].pto.x=x;
      decstr_vec[idx].pto.y=y;
      decstr_vec[idx].pto.z=z;
    }
    else if (atom_name==" CB ")
    {
      decstr_vec[idx].ptb.x=x;
      decstr_vec[idx].ptb.y=y;
      decstr_vec[idx].ptb.z=z;
    }
    else if (atom_name==" SC ")
    {
      decstr_vec[idx].ptsg.x=x;
      decstr_vec[idx].ptsg.y=y;
      decstr_vec[idx].ptsg.z=z;
    } 
    else 
    {
      decstr_vec[idx].ptscx[count]=x;
      decstr_vec[idx].ptscy[count]=y;
      decstr_vec[idx].ptscz[count]=z;
      count++;
    }
  }

  /* clean up */
  prev_resi.clear();
  resi.clear();
  atom_name.clear();
}

/* read chain content 'pdb_lines' into 'decstr_vec'. Only the most
 * basic coordinate information are parsed. Note that secondary
 * structures (ss2, ssm, stype) are all initialized to 'C'.
 * 'decstr_vec' must be pre-initialized to chain size */
void ParsePDB::loadPDBChain(
  const vector<string> &pdb_lines,
  vector<point3f> &decstr_vec
)
{
  string prev_resi="";
  string resi="";
  string atom_name="    ";
  string resn="UNK";
  float x,y,z;
  int idx=-1; // residue index, starting from 0
  int a; // amino acid type

  for (int i=0;i<pdb_lines.size();i++)
  {
    atom_name=pdb_lines[i].substr(12,4);
    resn=pdb_lines[i].substr(17,3);
    x=atof(pdb_lines[i].substr(30,8).c_str());
    y=atof(pdb_lines[i].substr(38,8).c_str());
    z=atof(pdb_lines[i].substr(46,8).c_str());
    resi=pdb_lines[i].substr(22,5);

    if (resi!=prev_resi) // new residue
    {
      idx+=1;
      prev_resi=resi;

      decstr_vec[idx].aaa  ='X';
      for (a=0;a<26;a++)
      {
        if (resn==aad3[a])
        {
          decstr_vec[idx].aaa=aad1[a];
          break;
        }
      }
      decstr_vec[idx].stype='H';
      decstr_vec[idx].ss3  ='C';
      decstr_vec[idx].ssm  ='C';
      decstr_vec[idx].resind=atoi(resi.c_str());
    }

    if (atom_name==" CA ")
    {
      decstr_vec[idx].x=x;
      decstr_vec[idx].y=y;
      decstr_vec[idx].z=z;
    }
    else if (atom_name==" N  ")
    {
      decstr_vec[idx].ptn.x=x;
      decstr_vec[idx].ptn.y=y;
      decstr_vec[idx].ptn.z=z;
    }
    else if (atom_name==" C  ")
    {
      decstr_vec[idx].ptc.x=x;
      decstr_vec[idx].ptc.y=y;
      decstr_vec[idx].ptc.z=z;
    }
    else if (atom_name==" H  ")
    {
      decstr_vec[idx].pth.x=x;
      decstr_vec[idx].pth.y=y;
      decstr_vec[idx].pth.z=z;
    }
    else if (atom_name==" O  ")
    {
      decstr_vec[idx].pto.x=x;
      decstr_vec[idx].pto.y=y;
      decstr_vec[idx].pto.z=z;
    }
    else if (atom_name==" CB ")
    {
      decstr_vec[idx].ptb.x=x;
      decstr_vec[idx].ptb.y=y;
      decstr_vec[idx].ptb.z=z;
    }
    else if (atom_name==" SC ")
    {
      decstr_vec[idx].ptsg.x=x;
      decstr_vec[idx].ptsg.y=y;
      decstr_vec[idx].ptsg.z=z;
    }
  }

  /* clean up */
  prev_resi.clear();
  resi.clear();
  atom_name.clear();
}

/* directly read PDB file 'pdbfile' into 'decstr_vec'. Only the most
 * basic coordinate information are parsed. Note that secondary
 * structures (ss2, ssm, stype) are all initialized to 'C' */
int ParsePDB::loadPDBChain(
  const char *pdbfile,
  vector<point3f>&decstr_vec
)
{    
  int seqnum=0;

  /* read file */
  vector <string> pdb_lines;
  string line;
  char prev_chainID=0;
  char chainID=0;
  string prev_resi="";
  string resi="";

  ifstream fp(pdbfile,ios::in);
  while(fp.good())
  {
    getline(fp,line);
    if(line.substr(0,3)=="END") break;
    if(line.length()<53 || line.substr(0,6)!="ATOM  " ||
       (line[16]!=' ' && line[16]!='A')) 
    {
      continue;
    }
    chainID=line[22];
    if(!prev_chainID) prev_chainID=chainID;
    if(chainID!=prev_chainID) break;

    resi=line.substr(22,5);
    if(resi!=prev_resi)
    {
      seqnum+=1;
      prev_resi=resi;
    }

    pdb_lines.push_back(line);
  }
  fp.close();

  /* parse file content into decstr_vec */
  point3f tmp_residue;
  decstr_vec.assign(seqnum,tmp_residue);
  loadPDBChain(pdb_lines,decstr_vec);

  /* clean up */
  line.clear();
  pdb_lines.clear();

  return seqnum;
}

bool ParsePDB::loadpdb(
  char *pdbfile
)
{
  FILE *filein;
  char oneline[600];
  bool modelend=false;
  filein = fopen(pdbfile, "rt");
  if(filein==NULL)
  {
    printf( "PDB file can not be opened %s\n", pdbfile);
    return false;
  }
  int i;
  for(i=0;i<5;i++) pdbname[i]=pdbfile[i];

  //init
  headnumchain=-1;
  indexchain=0;
  for(i=0;i<maxnumchains;i++)
  {
    if(aminoinchain[i]!=NULL)
    {
      delete[]aminoinchain[i];
      aminoinchain[i]=NULL;
    }
  }

  numhelix=0;numsheet=0;numturn=0;numlink=0;
  allochelix=30;allocsheet=30;allocturn=30;alloclink=30;
  numssbond=0;numsltbrg=0;numcispep=0;numhydbnd=0;numsite=0;
  allocssbond=30;allocsltbrg=30;alloccispep=30;allochydbnd=30;allocsite=30;
  numpromod=0;allocpromod=1;
  numproseq=0;allocproseq=800;
  pt.ismtrix=false;
  pt.isorigx=false;
  pt.isscale=false;
  if(helixl!=NULL)
  {
    free(helixl);
    helixl=NULL;
  }
  if(sheetl!=NULL)
  {
    free(sheetl);
    sheetl=NULL;
  }
  if(turnl!=NULL)
  {
    free(turnl);
    turnl=NULL;
  }
  if(linkl!=NULL)
  {
    free(linkl);
    linkl=NULL;
  }
  if(ssbondl!=NULL)
  {
    free(ssbondl);
    ssbondl=NULL;
  }
  if(sltbrgl!=NULL)
  {
    free(sltbrgl);
    sltbrgl=NULL;
  }
  if(cispepl!=NULL)
  {
    free(cispepl);
    cispepl=NULL;
  }
  if(hydbndl!=NULL)
  {
    free(hydbndl);
    hydbndl=NULL;
  }
  if(sitel!=NULL)
  {
    free(sitel);
    sitel=NULL;
  }
  if(promod!=NULL)
  {
    free(promod);
    promod=NULL;
  }
  if(proseq!=NULL)
  {
    free(proseq);
    proseq=NULL;
  }
  // allocpromod 1
  promod=new promodel[allocpromod];
  flagfirstmodel=true;
  flagbeginmodel=true;
  for(i=0;i<allocpromod;i++)
  {
    promod[i].nchain=0;
    for(int j=0;j<maxnumchains;j++)
    {
      promod[i].procha[j].chainseg.init=-1;
      promod[i].procha[j].chainseg.term=-1;
    }
  }
  //read
  char oldline[600];oldline[0]='\0';
  while(!feof(filein) && !modelend)
  {
    fgets(oneline,600,filein);
    if(strcmp(oneline,oldline)!=0 || oldline[0]=='\0')
      dealline(oneline);
    strcpy(oldline,oneline);
  }
  fclose(filein);
  
  //realloc
  if(promod[numpromod].procha[promod[numpromod].nchain].chainseg.term==-1)
  {
    promod[numpromod].procha[promod[numpromod].nchain].chainseg.term=numproseq-1;
  }
  if(numhelix!=0)
    helixl=(helixlink*)realloc(helixl, numhelix*sizeof(helixlink)); 
  if(numsheet!=0)
    sheetl=(sheetlink*)realloc(sheetl, numsheet*sizeof(sheetlink));
  if(numturn!=0)
    turnl=(turnlink*)realloc(turnl, numturn*sizeof(turnlink));
  if(numlink!=0)
    linkl=(linklink*)realloc(linkl, numlink*sizeof(linklink));
  if(numssbond!=0)
    ssbondl=(ssbondlink*)realloc(ssbondl, numssbond*sizeof(ssbondlink));
  if(numhydbnd!=0)
    hydbndl=(hydbndlink*)realloc(hydbndl, numhydbnd*sizeof(hydbndlink));
  if(numsltbrg!=0)
    sltbrgl=(sltbrglink*)realloc(sltbrgl, numsltbrg*sizeof(sltbrglink));
  if(numcispep!=0)
    cispepl=(cispeplink*)realloc(cispepl, numcispep*sizeof(cispeplink));
  if(numsite!=0)
    sitel=(sitelink*)realloc(sitel, numsite*sizeof(sitelink));
  if(numproseq!=0)
    proseq=(atom*)realloc(proseq, numproseq*sizeof(atom));
  if(numpromod!=-1)
    promod=(promodel*)realloc(promod, (numpromod+1)*sizeof(promodel));
  
  if(hinfo.dddd[2]!=-1)
  {
    hinfo.dnew[0]=hinfo.dddd[0];
    hinfo.dnew[1]=hinfo.dddd[1];
    hinfo.dnew[2]=hinfo.dddd[2];
  }
  if(hinfo.drev[2]!=-1)
  {
    if(hinfo.dnew[2]<hinfo.drev[2])
    {
      hinfo.dnew[2]=hinfo.drev[2];
      hinfo.dnew[1]=hinfo.drev[1];
      hinfo.dnew[0]=hinfo.drev[0];
    }
    else if(hinfo.dnew[2]==hinfo.drev[2])
    {
      if(hinfo.dnew[1]<hinfo.drev[1])
      {
        hinfo.dnew[1]=hinfo.drev[1];
        hinfo.dnew[0]=hinfo.drev[0];
      }
      else if(hinfo.dnew[1]==hinfo.drev[1])
      {
        if(hinfo.dnew[0]<hinfo.drev[0])
        {
          hinfo.dnew[0]=hinfo.drev[0];
        }
      }
    }
  }

  return true;
}

void ParsePDB::extractbb(
  int model, 
  int chain,
  int type
)
{
  //printf("%d %d %d\n",numproseq,model,numpromod);
  if(numproseq==0 || promod==NULL)
  {
    return;
  }
  if(model<-1 || model>numpromod)
  {
    printf("no this model\n");
    return;
  }
  if(type!=1 && type!=3 && type!=36 && type!=79)
  {
    printf("type should be 1 or 3 or 36 or 79\n");
    return;
  }
  if((model>=0 && chain> promod[model].nchain)  || chain<-1)
  {
    printf("no this chain\n");
    return;
  }
  
  if(bb)
  {
    delete[]bb;
    bb=NULL;
  }
  numbb=0;
  if(numproseq==0) return;
  int seqinit,seqterm;
  if(model==-1)
  {
    seqinit=promod[0].procha[0].chainseg.init;
    seqterm=promod[numpromod].procha[promod[numpromod].nchain].chainseg.term;
  }
  else if(chain==-1)
  {
    seqinit=promod[model].procha[0].chainseg.init;
    seqterm=promod[model].procha[promod[model].nchain].chainseg.term;
  }	
  else
  {
    seqinit=promod[model].procha[chain].chainseg.init;
    seqterm=promod[model].procha[chain].chainseg.term;
  }
  if(seqinit==-1 || seqterm==-1)
  {
    printf("no atoms in the model\n");
    return;
  }
  int allocbb=1+numproseq/4;
  bb=new boneinfo[allocbb];
  int i,j,k,l;
  char tpname[6];
  int indaa;
  BasicFunc bf;
  if(type==1)// ca
  {
    for(i=seqinit;i<=seqterm;i++)
    {
      j=i;
      while(j<=seqterm && proseq[i].resno==proseq[j].resno && proseq[i].ins==proseq[j].ins)
      {
        j++;
      }
      bb[numbb].indca=-1;
      bb[numbb].sst=0;
      if(proseq[i].residueid<0 && proseq[i].residueid>-9)
      {
        bb[numbb].indca=-1;
        bb[numbb].sst=0;
        for(k=i;k<j;k++)
        {
          if(/*proseq[k].simpletype==1 &&*/ proseq[k].detail==6)
          {
            bb[numbb].indca=k;			
            break;
          }
        }
        if(bb[numbb].indca==-1)
        {
          //printf("no ca atom in %d %d\n",proseq[i].seqno, proseq[i].resno);
        }
        bb[numbb].istart=i;
        bb[numbb].iend=j-1;
        bb[numbb].resid='0';
        bb[numbb].resind=proseq[i].resno;
        numbb++;
        if(numbb==allocbb)
        {
          allocbb*=2;
          bb = (boneinfo *)realloc(bb, allocbb*sizeof(boneinfo));
        }
        i=j-1;
        continue;
      }
	  
      for(k=i;k<j;k++)
      {
        if(/*proseq[k].simpletype==1 &&*/ proseq[k].detail==0)
        {
          bb[numbb].indca=k;			
          break;
        }
      }
      if(bb[numbb].indca==-1)
      {
        //printf("no ca atom in res %d %d\n",proseq[i].seqno, proseq[i].resno);
      }
      bb[numbb].istart=i;
      bb[numbb].iend=j-1;
      bb[numbb].resid=aad1[proseq[i].residueid];
      bb[numbb].resind=proseq[i].resno;
      if(bb[numbb].indca==-1 /*&& (proseq[i].residueid>19 || bb[numbb].iend-bb[numbb].istart<2)*/)
      {
        numbb--;
        //printf("eliminate the res %d %d\n",proseq[i].seqno, proseq[i].resno);
      }
      numbb++;
      if(numbb==allocbb)
      {
        allocbb*=2;
        bb=(boneinfo *)realloc(bb, allocbb*sizeof(boneinfo));
      }
      i=j-1;
    }
    if(numbb!=0)
      bb=(boneinfo *)realloc(bb, numbb*sizeof(boneinfo));
  }
  else if(type==3)//n ca c
  {
    for(i=seqinit;i<=seqterm;i++)
    {
      j=i;
      while(j<=seqterm && proseq[i].resno==proseq[j].resno  && proseq[i].ins==proseq[j].ins)
      {
        j++;
      }
	  
      bb[numbb].indn=-1;bb[numbb].indca=-1;bb[numbb].indc=-1;bb[numbb].indo=-1;bb[numbb].indcb=-1;bb[numbb].indf=-1;
      bb[numbb].indh=-1;bb[numbb].indha=-1;
      bb[numbb].sst=0;
      if(proseq[i].residueid<0 && proseq[i].residueid>-9)
      {
        bb[numbb].indca=-1;
        bb[numbb].sst=0;
        for(k=i;k<j;k++)
        {
          if(/*proseq[k].simpletype==1 &&*/ proseq[k].detail==6)
          {
            bb[numbb].indca=k;			
            break;
          }
        }
        if(bb[numbb].indca==-1)
        {
          //printf("no ca atom in nucleic acid %d %d\n",proseq[i].seqno, proseq[i].resno);
        }
        bb[numbb].istart=i;
        bb[numbb].iend=j-1;
        bb[numbb].resid='0';
        bb[numbb].resind=proseq[i].resno;
        numbb++;
        if(numbb==allocbb)
        {
          allocbb*=2;
          bb=(boneinfo *)realloc(bb, allocbb*sizeof(boneinfo));
        }
        i=j-1;      
        continue;
      }
	  
      for(k=i;k<j;k++)
      {
        if(/*proseq[k].simpletype==1 && */proseq[k].detail==2)
        {
          bb[numbb].indn=k;	
        }
        else if(/*proseq[k].simpletype==1 && */proseq[k].detail==0)
        {
          bb[numbb].indca=k;
        }
        //else 
        else if(/*proseq[k].simpletype==1 && */proseq[k].detail==1)
        {
          bb[numbb].indc=k;
        }
        else if(/*proseq[k].simpletype==1 && */proseq[k].detail==3)
        {
          bb[numbb].indo=k;
        }
        else if(/*proseq[k].simpletype==1 && */proseq[k].detail==7 && proseq[k].detailtype[0]=='C' &&
                proseq[k].detailtype[1]=='B')
        {
          bb[numbb].indcb=k;
        }
        else if(proseq[k].detailtype[0]==ffff[proseq[k].residueid][0] &&
                proseq[k].detailtype[1]==ffff[proseq[k].residueid][1] &&
                proseq[k].detailtype[2]==ffff[proseq[k].residueid][2])
        {
          bb[numbb].indf=k;
        }
        else if(/*proseq[k].simpletype==1 && */proseq[k].detail==5 && proseq[k].detailtype[0]=='H' &&
                proseq[k].detailtype[1]==' ')
        {
          bb[numbb].indh=k;
        }
        else if(/*proseq[k].simpletype==1 && */
		proseq[k].residueid==5 && proseq[k].detailtype[0]=='1' &&
		proseq[k].detailtype[1]=='H' && proseq[k].detailtype[2]=='A')
        {
          bb[numbb].indha=k;
        }
        else if(/*proseq[k].simpletype==1 && */ proseq[k].detailtype[0]=='H' && proseq[k].detailtype[1]=='A')
        {
          bb[numbb].indha=k;
        }
      }
      if(bb[numbb].indca==-1)
      {
        //printf("no ca atom in res %d %d\n",proseq[i].seqno, proseq[i].resno);
      }
      bb[numbb].istart=i;
      bb[numbb].iend=j-1;
      bb[numbb].resid=aad1[proseq[i].residueid];
      bb[numbb].resind=proseq[i].resno;
      if(bb[numbb].indca==-1 /*&& (proseq[i].residueid>19 || bb[numbb].iend-bb[numbb].istart<2)*/)
      {
        numbb--;
        //printf("eliminate the res %d %d\n",proseq[i].seqno, proseq[i].resno);
      }
      numbb++;
      if(numbb==allocbb)
      {
        allocbb*=2;
        bb=(boneinfo *)realloc(bb, allocbb*sizeof(boneinfo));
      }
      i=j-1;
      //notice that i will add 1 in the for()
    }
    if(numbb!=0)
      bb=(boneinfo *)realloc(bb, numbb*sizeof(boneinfo));
    effbb=0;
    for(i=0;i<numbb;i++)
    {
      if(bb[i].indca!=-1 && bb[i].indc!=-1 && bb[i].indn!=-1 && bb[i].indo!=-1)
      {
        effbb++; 
      }
    }
  }
  else if(type==36)//all heavy atoms
  {
    for(i=seqinit;i<=seqterm;i++)
    {
      j=i;
      while(j<=seqterm && proseq[i].resno==proseq[j].resno  && proseq[i].ins==proseq[j].ins)
      {
        j++;
      }
      for(k=0;k<36;k++)
      {
        bb[numbb].indall[k]=-1;
      }
      bb[numbb].indn=-1;bb[numbb].indca=-1;bb[numbb].indc=-1;bb[numbb].indo=-1;bb[numbb].indcb=-1;bb[numbb].indf=-1;
      bb[numbb].indh=-1;bb[numbb].indha=-1;
      bb[numbb].sst=0;
      for(k=i;k<j;k++)
      {
        for(l=0;l<proseq[k].detailtype[5];l++)
        {
          tpname[l]=proseq[k].detailtype[l];
        }
        tpname[proseq[k].detailtype[5]]='\0';
        indaa=bf.atomid(tpname);
        if(indaa>=36) continue;
        //printf("%s %d\n",tpname,indaa);
        bb[numbb].indall[indaa]=k;
        if(/*proseq[k].simpletype==1 && */proseq[k].detail==2)
        {
          bb[numbb].indn=k;	
        }
        else if(/*proseq[k].simpletype==1 && */proseq[k].detail==0)
        {
          bb[numbb].indca=k;
        }
        //else 
        else if(/*proseq[k].simpletype==1 && */proseq[k].detail==1)
        {
          bb[numbb].indc=k;
        }
        else if(/*proseq[k].simpletype==1 && */proseq[k].detail==3)
        {
          bb[numbb].indo=k;
        }
        else if(/*proseq[k].simpletype==1 && */proseq[k].detail==7 && proseq[k].detailtype[0]=='C' &&
                proseq[k].detailtype[1]=='B')
        {
          bb[numbb].indcb=k;
        }
        else if(proseq[k].detailtype[0]==ffff[proseq[k].residueid][0] &&
                proseq[k].detailtype[1]==ffff[proseq[k].residueid][1] &&
                proseq[k].detailtype[2]==ffff[proseq[k].residueid][2])
        {
          bb[numbb].indf=k;
        }
        else if(/*proseq[k].simpletype==1 && */proseq[k].detail==5 && proseq[k].detailtype[0]=='H' &&
                proseq[k].detailtype[1]==' ')
        {
          bb[numbb].indh=k;
        }
        else if(/*proseq[k].simpletype==1 && */
                proseq[k].residueid==5 && proseq[k].detailtype[0]=='1' &&
                proseq[k].detailtype[1]=='H' && proseq[k].detailtype[2]=='A')
        {
          bb[numbb].indha=k;
        }
        else if(/*proseq[k].simpletype==1 && */ proseq[k].detailtype[0]=='H' && proseq[k].detailtype[1]=='A')
        {
          bb[numbb].indha=k;
        }
      } 
      if(bb[numbb].indca==-1)
      {
        //printf("no ca atom in res %d %d\n",proseq[i].seqno, proseq[i].resno);
      }
      bb[numbb].istart=i;
      bb[numbb].iend=j-1;
      bb[numbb].resid=aad1[proseq[i].residueid];
      bb[numbb].resind=proseq[i].resno;
      if(bb[numbb].indca==-1 && (proseq[i].residueid>19 || bb[numbb].iend-bb[numbb].istart<2))
      {
        numbb--;
        //printf("eliminate the res %d %d\n",proseq[i].seqno, proseq[i].resno);
      }
      numbb++;
      if(numbb==allocbb)
      {
        allocbb*=2;
        bb=(boneinfo *)realloc(bb, allocbb*sizeof(boneinfo));
      }
      i=j-1;
      //notice that i will add 1 in the for()
    }
    if(numbb!=0)
      bb=(boneinfo *)realloc(bb, numbb*sizeof(boneinfo));
    effbb=0;
    for(i=0;i<numbb;i++)
    {
      if(bb[i].indca!=-1 && bb[i].indc!=-1 && bb[i].indn!=-1 && bb[i].indo!=-1)
      {
        effbb++;
      }
    }
  }
  else if(type==79)//all atoms
  {
    for(i=seqinit;i<=seqterm;i++)
    {
      j=i;
      while(j<=seqterm && proseq[i].resno==proseq[j].resno  && proseq[i].ins==proseq[j].ins)
      {
        j++;
      }
      for(k=0;k<83;k++)
      {
        bb[numbb].indall[k]=-1;
      }
      bb[numbb].indn=-1;bb[numbb].indca=-1;bb[numbb].indc=-1;bb[numbb].indo=-1;bb[numbb].indcb=-1;bb[numbb].indf=-1;
      bb[numbb].indh=-1;bb[numbb].indha=-1;
      bb[numbb].sst=0;
      for(k=i;k<j;k++)
      {
        for(l=0;l<proseq[k].detailtype[5];l++)
        {
          tpname[l]=proseq[k].detailtype[l];
        }
        tpname[proseq[k].detailtype[5]]='\0';
        indaa=bf.atomid(tpname);
        if(indaa>=83) continue;
        //printf("%s %d\n",tpname,indaa);
        bb[numbb].indall[indaa]=k;
        if(/*proseq[k].simpletype==1 && */proseq[k].detail==2)
        {
          bb[numbb].indn=k;	
        }
        else if(/*proseq[k].simpletype==1 && */proseq[k].detail==0)
        {
          bb[numbb].indca=k;
        }
        //else 
        else if(/*proseq[k].simpletype==1 && */proseq[k].detail==1)
        {
          bb[numbb].indc=k;
        }
        else if(/*proseq[k].simpletype==1 && */proseq[k].detail==3)
        {
          bb[numbb].indo=k;
        }
        else if(/*proseq[k].simpletype==1 && */proseq[k].detail==7 && proseq[k].detailtype[0]=='C' &&
		proseq[k].detailtype[1]=='B')
        {
          bb[numbb].indcb=k;
        }
        else if(proseq[k].detailtype[0]==ffff[proseq[k].residueid][0] &&
                proseq[k].detailtype[1]==ffff[proseq[k].residueid][1] &&
                proseq[k].detailtype[2]==ffff[proseq[k].residueid][2])
        {
          bb[numbb].indf=k;
        }
        else if(/*proseq[k].simpletype==1 && */proseq[k].detail==5 && proseq[k].detailtype[0]=='H' &&
                proseq[k].detailtype[1]==' ')
        {
          bb[numbb].indh=k;
        }
        else if(/*proseq[k].simpletype==1 && */
                proseq[k].residueid==5 && proseq[k].detailtype[0]=='1' &&
                proseq[k].detailtype[1]=='H' && proseq[k].detailtype[2]=='A')
        {
          bb[numbb].indha=k;
        }
        else if(/*proseq[k].simpletype==1 && */ proseq[k].detailtype[0]=='H' && proseq[k].detailtype[1]=='A')
        {
          bb[numbb].indha=k;
        }
      }
      if(bb[numbb].indca==-1)
      {
        //printf("no ca atom in res %d %d\n",proseq[i].seqno, proseq[i].resno);
      }
      bb[numbb].istart=i;
      bb[numbb].iend=j-1;
      bb[numbb].resid=aad1[proseq[i].residueid];
      bb[numbb].resind=proseq[i].resno;
      if(bb[numbb].indca==-1 && (proseq[i].residueid>19 || bb[numbb].iend-bb[numbb].istart<2))
      {
        numbb--;
        //printf("eliminate the res %d %d\n",proseq[i].seqno, proseq[i].resno);
      }
      numbb++;
      if(numbb==allocbb)
      {
        allocbb*=2;
        bb=(boneinfo *)realloc(bb, allocbb*sizeof(boneinfo));
      }
      i=j-1;
      //notice that i will add 1 in the for()
    }
    if(numbb!=0)
      bb=(boneinfo *)realloc(bb, numbb*sizeof(boneinfo));
    effbb=0;
    for(i=0;i<numbb;i++)
    {
      if(bb[i].indca!=-1 && bb[i].indc!=-1 && bb[i].indn!=-1 && bb[i].indo!=-1)
      {
        effbb++;
      }
    }
  }
  if(numbb!=(proseq[seqterm].resno-proseq[seqinit].resno+1))
  {
    //printf("Lack some residues %d-%d=%d\n",proseq[seqterm].resno-proseq[seqinit].resno+1,numbb,
    //proseq[seqterm].resno-proseq[seqinit].resno+1-numbb);
  }
  if(resmap)
  {
    delete[]resmap;
  }
  resnumber=0;
  for(i=0;i<numbb;i++)
  {
    if(bb[i].resind>resnumber)
      resnumber=bb[i].resind;
  }
  resmap=new int[1000+resnumber+1];
  for(i=0;i<1000+resnumber+1;i++)
  {
    resmap[i]=-1;
  }
  for(i=0;i<numbb;i++)
  {
    resmap[1000+bb[i].resind]=i;
  }
}

int ParsePDB::dealline(
  char *string
)
{
  int i,j;
  if(string[0]=='H' && string[1]=='E' && string[2]=='A' && string[3]=='D' && string[4]=='E' && string[5]=='R')
  {
    //HEADER    Created by PULCHRA version  0.99900
    char *ptr=NULL;
    ptr=strstr(string+6,"PULCHRA");
    if(!ptr && strlen(string)>60)
    {	      
      char mmm[8];
      sscanf(string+62,"%s",hinfo.name);
      mmm[3]='\0';mmm[0]=string[53];mmm[1]=string[54];mmm[2]=string[55];
      sscanf(string+50,"%d",&hinfo.dddd[0]);
      sscanf(string+57,"%d",&hinfo.dddd[2]);
      if(hinfo.dddd[2]>=50) hinfo.dddd[2]+=1900;
      else hinfo.dddd[2]+=2000;
      hinfo.dddd[1]=-1;
      for(i=0;i<12;i++)
      {
        if(strcmp(mmm,strmonth[i])==0)
        {
          hinfo.dddd[1]=i;
          break;
        }
      }
    }
  }
  else if(string[0]=='O' && string[1]=='B' && string[2]=='S' && string[3]=='L' && string[4]=='T' && string[5]=='E')
  {
    hinfo.fobs=true;	 
  }
  else if(string[0]=='S' && string[1]=='O' && string[2]=='U' && string[3]=='R' && string[4]=='C' && string[5]=='E')
  { 
    char *ptr=NULL;
    ptr=strstr(string+7,"HUMAN");
    if(ptr && hinfo.species==-1)
    {
      hinfo.species=0;
    }	
    ptr=strstr(string+7,"HOMO SAPIENS");
    if(ptr)
    {
      hinfo.species=1;
    }	
  }
  else if(string[0]=='E' && string[1]=='X' && string[2]=='P' && string[3]=='D' && string[4]=='T' && string[5]=='A')
  {
    char *ptr=NULL;
    ptr=strstr(string+9,"DIFFRACTION");
    if(ptr)
    {
      hinfo.iept=0;
    }
    ptr=strstr(string+9,"NMR");
    if(ptr)
    {
      hinfo.iept=1;
    }
    ptr=strstr(string+9,"ELECTRON MICROSCOPY");
    if(ptr)
    {
      hinfo.iept=2;
    }		
    ptr=strstr(string+9,"ELECTRON CRYSTALLOGRAPHY");
    if(ptr)
    {
      hinfo.iept=4;
    }
    ptr=strstr(string+9,"SOLUTION SCATTERING");
    if(ptr)
    {
      hinfo.iept=5;
    }
    ptr=strstr(string+9,"INFRARED SPECTROSCOPY");
    if(ptr)
    {
      hinfo.iept=6;
    }
    ptr=strstr(string+9,"THEORETICAL");//put in the end
    if(ptr)
    {
      hinfo.iept=3;
    }
  }
  else if(!hinfo.frev && string[0]=='R' && string[1]=='E' && string[2]=='V' && string[3]=='D' && string[4]=='A' && string[5]=='T')
  {
    char mmm[8];
    mmm[3]='\0';mmm[0]=string[16];mmm[1]=string[17];mmm[2]=string[18];
    sscanf(string+13,"%d",&hinfo.drev[0]);
    sscanf(string+20,"%d",&hinfo.drev[2]);
    if(hinfo.drev[2]>=50) hinfo.drev[2]+=1900;
    else hinfo.drev[2]+=2000;
    hinfo.drev[1]=-1;
    for(i=0;i<12;i++)
    {
      if(strcmp(mmm,strmonth[i])==0)
      {
        hinfo.drev[1]=i;
        break;
      }   
    }
    hinfo.frev=true;
  }
  else if(string[0]=='R' && string[1]=='E' && string[2]=='M' && string[3]=='A' && string[4]=='R' && string[5]=='K' && string[9]=='2')
  {
    if(string[11]=='R' && string[12]=='E' && string[13]=='S' && string[31]=='A' && string[32]=='N' && string[33]=='G')
    {
      sscanf(string+22,"%f",&hinfo.fres);
    }
    else if(string[11]=='R' && string[12]=='E' && string[13]=='S' && string[28]=='A' && string[29]=='N' && string[30]=='G')
    {
      sscanf(string+22,"%f",&hinfo.fres);
    }
    else if(string[11]=='R' && string[12]=='E' && string[13]=='S' && string[23]=='N' && string[24]=='O' && string[25]=='T')
    {

    }
  }
  //SEQRES  28 A  360  ALA THR ALA VAL LYS ILE THR LEU LEU                  2ACH  74
  //SEQRES   1 B   40  SER ALA LEU VAL GLU THR ARG THR ILE VAL ARG PHE ASN  2ACH  75
  else if (string[0]=='S' && string[1]=='E' && string[2]=='Q' && string[3]=='R' && string[4]=='E' && string[5]=='S')
  {
    char temnum[5];
    temnum[4]='\0';
    char temid;
    int  tempnum;
    bool isanew;
    temid=string[11];
    if(headnumchain==-1)
    {
      isanew=true;
    }
    else if(headidchain[headnumchain]!=temid)
    {
      isanew=true;
    }
    else
    {
      isanew=false;
    }
    for(i=0;i<4;i++)
    {
      temnum[i]=string[13+i];
    }
    tempnum=atoi(temnum);
    if(headnumchain>=maxnumchains-1)
    {
      printf("Too many chains\n");
      return false;
    }
    else if(isanew)//start of a new chain
    {
      if(headnumchain!=-1)
      {
        if(indexchain!=numinchain[headnumchain])
        {		
          //sprintf(tmpstr,"ideal%d,actual%d,ideal%d,chain%d",tempnum,indexchain,numinchain[headnumchain],headnumchain);
          //printf("Not equal as claimed %s\n",tmpstr);
        }
      }
      indexchain=0;
      char temami[3];	
      int aminum;
      aminoinchain[++headnumchain]=new char[tempnum];
      headidchain[headnumchain]=temid;
      numinchain[headnumchain]=tempnum;
      for(i=0;i<13;i++)
      {
        for(j=0;j<3;j++)
        {
          temami[j]=string[19+i*4+j];
        }	
        aminum=aminoid(temami);
        if(aminum>=0 && aminum<24 && indexchain<tempnum)
        {
          aminoinchain[headnumchain][indexchain]=aminum;
          indexchain++;
        }
        //else if(aminum==50 && indexchain<tempnum)
        //{
        //  aminoinchain[headnumchain][indexchain]=aminum;
        //  indexchain++;
        //}
        else if(temami[0]==' ' && temami[1]==' ' && temami[2]==' '  && indexchain!=tempnum)
        {
          printf("Empty but not end1\n");
        }
        else if(aminum>=0 && aminum<24 && !(temami[0]==' ' && temami[1]==' ' && temami[2]==' ')
                && indexchain==tempnum)
        {
          printf("Too many amino1\n");
        }
      }	
    }
    else// push   
    {
      char temami[3];		
      int aminum;
      for(i=0;i<13;i++)
      {
        for(j=0;j<3;j++)
        {
          temami[j]=string[19+i*4+j];
        }
        aminum=aminoid(temami);
        if(aminum>=0 && aminum<24 && indexchain<tempnum)
        {
          aminoinchain[headnumchain][indexchain]=aminum;
          indexchain++;
        }
        //else if(aminum==50 && indexchain<tempnum)
        //{
        //  aminoinchain[headnumchain][indexchain]=aminum;
        //  indexchain++;
        //}
        else if(temami[0]==' ' && temami[1]==' ' && temami[2]==' ' && indexchain!=tempnum)
        {
          printf("Empty but not end2\n");
        }
        else if(aminum>=0 && aminum<24 && !(temami[0]==' ' && temami[1]==' ' && temami[2]==' ') 
                &&  indexchain==tempnum)
        {
          printf("Too many amino2\n");
        }
      }
    }	
    return 1;
  }
  else if (string[0]=='O' && string[1]=='R' && string[2]=='I' && string[3]=='G' && string[4]=='X' ) 
  {
    pt.isorigx=true;
    char cid[2];
    cid[1]='\0';
    cid[0]=string[5];
    int id=atoi(cid)-1;
    double txyz[3],ttran;
    char temxyz[10],temtran[10];   
    for(i=0;i<3;i++)
    {
      for(j=0;j<10;j++)
      {
        temxyz[j]=string[10+10*i+j];
      }
      txyz[i]=atof(temxyz);
      pt.origxn[id][i]=txyz[i];
    }
    for(j=0;j<10;j++)
    {
      temtran[j]=string[45+j];
    }    
    ttran=atof(temtran);
    pt.origxn[id][3]=ttran;
    return 2;
  }
  else if (string[0]=='S' && string[1]=='C' && string[2]=='A' && string[3]=='L' && string[4]=='E' ) 
  {
    pt.isscale=true;
    char cid[2];
    cid[1]='\0';
    cid[0]=string[5];
    int id=atoi(cid)-1;
    double txyz[3],ttran;
    char temxyz[10],temtran[10];
    for(i=0;i<3;i++)
    {
      for(j=0;j<10;j++)
      {
        temxyz[j]=string[10+10*i+j];
      }
      txyz[i]=atof(temxyz);
      pt.scalen[id][i]=txyz[i];
    }
    for(j=0;j<10;j++)
    {
      temtran[j]=string[45+j];
    }
    ttran=atof(temtran);
    pt.scalen[id][3]=ttran;
    return 3;
  }
  else if (string[0]=='M' && string[1]=='T' && string[2]=='R' && string[3]=='I' && string[4]=='X' ) 
  {
    pt.ismtrix=true;
    char cid[2];
    cid[1]='\0';
    cid[0]=string[5];
    int id=atoi(cid)-1;
    double txyz[3],ttran;
    char temxyz[10],temtran[10];
    for(i=0;i<3;i++)
    {
      for(j=0;j<10;j++)
      {
        temxyz[j]=string[10+10*i+j];
      }
      txyz[i]=atof(temxyz);
      pt.mtrixn[id][i]=txyz[i];
    }
    for(j=0;j<10;j++)
    {
      temtran[j]=string[45+j];
    }
    ttran=atof(temtran);
    pt.mtrixn[id][3]=ttran;
    return 4;
  }
  //HELIX    1   1 ILE A  146  ALA A  148  5                                   3  
  //HELIX    1 A1  GLY     15  GLN     25  1                                1CRP  86
  else if (string[0]=='H' && string[1]=='E' && string[2]=='L' && string[3]=='I' && string[4]=='X' ) 
  {
    if(numhelix==0)
    {
      helixl=new helixlink[allochelix];
    }
    else if(numhelix==allochelix-1)
    {
      allochelix*=2;
      helixl=(helixlink*)realloc(helixl, allochelix*sizeof(helixlink)); 		
    }
    char temami[3],temchainid,temamiindex[5],temclass[3];
    temamiindex[4]='\0';
    temclass[2]='\0';
    int aminum,amiindex,helclass;
    //start
    for(i=0;i<3;i++)
    {
      temami[i]=string[15+i];
    }
    temchainid=string[19];
    aminum=aminoid(temami);
    for(i=0;i<4;i++)
    {
      temamiindex[i]=string[21+i];
    }
    amiindex=atoi(temamiindex);
    helixl[numhelix].initchainid=temchainid;
    helixl[numhelix].initid=aminum;
    helixl[numhelix].initindex=amiindex;
    //end
    for(i=0;i<3;i++)
    {
      temami[i]=string[27+i];
    }
    temchainid=string[31];
    aminum=aminoid(temami);
    for(i=0;i<4;i++)
    {
      temamiindex[i]=string[33+i];
    }
    amiindex=atoi(temamiindex);
    helixl[numhelix].endchainid=temchainid;
    helixl[numhelix].endid=aminum;
    helixl[numhelix].endindex=amiindex;
    for(i=0;i<2;i++)
    {
      temclass[i]=string[38+i];
    }
    helclass=atoi(temclass);
    int temleng;
    char temlength[6];
    temlength[5]='\0';
    for(i=0;i<5;i++)
    {
      temlength[i]=string[71+i];
    }
    temleng=atoi(temlength);
    helixl[numhelix].helixlength=temleng;
    helixl[numhelix].helixtype=helclass;
    numhelix++;
    return 5;
  }
  //SHEET    1 S1  5 ASP    38  ILE    46  0                                1CRP  91
  //SHEET    2 S1  5 GLU    49  ASP    57 -1  O  LEU    53   N  LYS    42   1CRP  92
  //SHEET    1   B 2 PHE A 215  THR A 218  0                                        
  //SHEET    2   B 2 THR A 221  PHE A 223 -1  N  PHE A 223   O  PHE A 215           		
  else if (string[0]=='S' && string[1]=='H' && string[2]=='E' && string[3]=='E' && string[4]=='T' ) 
  {
    if(numsheet==0)
    {
      sheetl=new sheetlink[allocsheet];
    }
    else if(numsheet==allocsheet-1)
    {
      allocsheet*=2;
      sheetl=(sheetlink*)realloc(sheetl, allocsheet*sizeof(sheetlink)); 		
    }
    char temami[3],temchainid,temamiindex[5],temclass[3];
    temamiindex[4]='\0';  
    temclass[2]='\0';
    int aminum,amiindex,sheclass;
    //start
    for(i=0;i<3;i++)
    {
      temami[i]=string[17+i];
    }
    temchainid=string[21];
    aminum=aminoid(temami);
    for(i=0;i<4;i++)
    {
      temamiindex[i]=string[22+i];
    }
    amiindex=atoi(temamiindex);
    sheetl[numsheet].initchainid=temchainid;
    sheetl[numsheet].initid=aminum;
    sheetl[numsheet].initindex=amiindex;
    //end
    for(i=0;i<3;i++)
    {
      temami[i]=string[28+i];
    }
    temchainid=string[32]; 
    aminum=aminoid(temami);
    for(i=0;i<4;i++)
    {
      temamiindex[i]=string[33+i];
    }
    amiindex=atoi(temamiindex);
    sheetl[numsheet].endchainid=temchainid;
    sheetl[numsheet].endid=aminum;
    sheetl[numsheet].endindex=amiindex;
    for(i=0;i<2;i++)
    {
      temclass[i]=string[38+i];
    }
    sheclass=atoi(temclass);
    sheetl[numsheet].sheettype=sheclass;

    //curr
    for(i=0;i<4;i++)
    {
      sheetl[numsheet].curratom[i]=string[41+i];
    }
    for(i=0;i<3;i++)
    {
      temami[i]=string[45+i];
    }
    temchainid=string[49];
    aminum=aminoid(temami);
    for(i=0;i<4;i++)
    {
      temamiindex[i]=string[50+i];
    }
    amiindex=atoi(temamiindex);
    sheetl[numsheet].currchainid=temchainid;
    sheetl[numsheet].currid=aminum;
    sheetl[numsheet].currindex=amiindex;
    //prev
    for(i=0;i<4;i++)
    {
      sheetl[numsheet].prevatom[i]=string[56+i];
    }
    for(i=0;i<3;i++)
    {
      temami[i]=string[60+i];
    }
    temchainid=string[64];
    aminum=aminoid(temami);
    for(i=0;i<4;i++)
    {
      temamiindex[i]=string[65+i];
    }
    amiindex=atoi(temamiindex);
    sheetl[numsheet].prevchainid=temchainid;
    sheetl[numsheet].previd=aminum;
    sheetl[numsheet].previndex=amiindex;
    numsheet++;
    return 6;
  }
  //TURN     1 T1A SER A  11  ASP A  14     TYPE III                        8ATC 377
  //TURN     2 T2A ASP A 129  ASN A 132     TYPE I                          8ATC 378
  else if (string[0]=='T' && string[1]=='U' && string[2]=='R' && string[3]=='N' && string[4]==' ' ) 
  {
    if(numturn==0)
    {
      turnl=new turnlink[allocturn];
    }
    else if(numturn==allocturn-1)
    {
      allocturn*=2;
      turnl=(turnlink*)realloc(turnl, allocturn*sizeof(turnlink)); 		
    }
    char temami[3],temchainid,temamiindex[5];
    temamiindex[4]='\0';
    int aminum,amiindex;
    //start
    for(i=0;i<3;i++)
    {
      temami[i]=string[15+i];
    }
    temchainid=string[19];
    aminum=aminoid(temami);
    for(i=0;i<4;i++)
    {
      temamiindex[i]=string[20+i];
    }
    amiindex=atoi(temamiindex);
    turnl[numturn].initchainid=temchainid;
    turnl[numturn].initid=aminum;
    turnl[numturn].initindex=amiindex;
    //end
    for(i=0;i<3;i++)
    {
      temami[i]=string[26+i];
    }
    temchainid=string[30];
    aminum=aminoid(temami);
    for(i=0;i<4;i++)
    {
      temamiindex[i]=string[31+i];
    }
    amiindex=atoi(temamiindex);
    turnl[numturn].endchainid=temchainid;
    turnl[numturn].endid=aminum;
    turnl[numturn].endindex=amiindex;
    numturn++;
    return 7;
  }
  //LINK         C   ACE C 100                 N   PTR C 101  
  else if (string[0]=='L' && string[1]=='I' && string[2]=='N' && string[3]=='K' && string[4]==' ' ) 
  {
    if(numlink==0)
    {
      linkl=new linklink[alloclink];
    }
    else if(numlink==alloclink-1)
    {
      alloclink*=2;
      linkl=(linklink*)realloc(linkl, alloclink*sizeof(linklink)); 		
    }
    char temchainid,temamiindex[5];
    temamiindex[4]='\0';
    int amiindex;
    //start
    for(i=0;i<4;i++)
    {
      linkl[numlink].initatom[i]=string[12+i];
    }
    for(i=0;i<3;i++)
    {
      linkl[numlink].initid[i]=string[17+i];
    }
    temchainid=string[21];
    for(i=0;i<4;i++)
    {
      temamiindex[i]=string[22+i];
    }
    amiindex=atoi(temamiindex);
    linkl[numlink].initchainid=temchainid;
    linkl[numlink].initindex=amiindex;
    //end
    for(i=0;i<4;i++)
    {
      linkl[numlink].endatom[i]=string[42+i];
    }
    for(i=0;i<3;i++)
    {
      linkl[numlink].endid[i]=string[47+i];
    }
    temchainid=string[51];
    for(i=0;i<4;i++)
    {
      temamiindex[i]=string[53+i];
    }
    amiindex=atoi(temamiindex);
    linkl[numlink].endchainid=temchainid;
    linkl[numlink].endindex=amiindex;
    numlink++;
    return 8;
  }
  //SSBOND   1 CYS E   48    CYS E   51                          2555  
  //SSBOND   1 CYS G  119    CYS G  205  
  else if (string[0]=='S' && string[1]=='S' && string[2]=='B' && string[3]=='O' && string[4]=='N' && string[5]=='D') 
  {
    if(numssbond==0)
    {
      ssbondl=new ssbondlink[allocssbond];
    }
    else if(numssbond==allocssbond-1)
    {
      allocssbond*=2;
      ssbondl=(ssbondlink*)realloc(ssbondl, allocssbond*sizeof(ssbondlink)); 		
    }
    //start
    ssbondl[numssbond].initchainid=string[15];
    char temnumc[5];
    temnumc[4]='\0';
    int temnumi;
    for(i=0;i<4;i++)
    {
      temnumc[i]=string[17+i];
    }
    temnumi=atoi(temnumc);
    ssbondl[numssbond].initindex=temnumi;
    //end
    ssbondl[numssbond].endchainid=string[29];
    for(i=0;i<4;i++)
    {
      temnumc[i]=string[31+i];
    }
    temnumi=atoi(temnumc);
    ssbondl[numssbond].endindex=temnumi;
    numssbond++;
    return 9;
  }
  //HYDBND       NH2 ARG    111                 OD1 ASP    149   1555		
  else if (string[0]=='H' && string[1]=='Y' && string[2]=='D' && string[3]=='B' && string[4]=='N' && string[5]=='D') 
  {
    if(numhydbnd==0)
    {
      hydbndl=new hydbndlink[allochydbnd];
    }
    else if(numhydbnd==allochydbnd-1)
    {
      allochydbnd*=2;
      hydbndl=(hydbndlink*)realloc(hydbndl, allochydbnd*sizeof(hydbndlink)); 		
    }
    //start
    for(i=0;i<4;i++)
    {
      hydbndl[numhydbnd].initatom[i]=string[12+i];
    }
    hydbndl[numhydbnd].initchainid=string[21];
    char temami[3],temseq[5];
    temseq[4]='\0';
    int temnum,temseqi;
    for(i=0;i<3;i++)
    {
      temami[i]=string[17+i];
    }
    temnum=aminoid(temami);
    hydbndl[numhydbnd].initid=temnum;
    for(i=0;i<4;i++)
    {
      temseq[i]=string[22+i];
    }
    temseqi=atoi(temseq);
    hydbndl[numhydbnd].initindex=temseqi;
    //end
    for(i=0;i<4;i++)
    {
      hydbndl[numhydbnd].endatom[i]=string[43+i];
    }
    hydbndl[numhydbnd].endchainid=string[52];
    for(i=0;i<3;i++)
    {
      temami[i]=string[48+i];
    }
    temnum=aminoid(temami);
    hydbndl[numhydbnd].endid=temnum;
    for(i=0;i<4;i++)
    {
      temseq[i]=string[53+i];
    }
    temseqi=atoi(temseq);
    hydbndl[numhydbnd].endindex=temseqi;
		
    numhydbnd++;
    return 10;
  }
  //SLTBRG       O   GLU    10                 NZ  LYS  115             3654		
  else if (string[0]=='S' && string[1]=='L' && string[2]=='T' && string[3]=='B' && string[4]=='R' && string[5]=='G') 
  {
    if(numsltbrg==0)
    {
      sltbrgl=new sltbrglink[allocsltbrg];
    }
    else if(numsltbrg==allocsltbrg-1)
    {
      allocsltbrg*=2;
      sltbrgl=(sltbrglink*)realloc(sltbrgl, allocsltbrg*sizeof(sltbrglink)); 		
    }
    //start
    for(i=0;i<4;i++)
    {
      sltbrgl[numsltbrg].initatom[i]=string[12+i];
    }
    sltbrgl[numsltbrg].initchainid=string[21];
    char temami[3],temseq[5];
    temseq[4]='\0';
    int temnum,temseqi;
    for(i=0;i<3;i++)
    {
      temami[i]=string[17+i];
    }
    temnum=aminoid(temami);
    sltbrgl[numsltbrg].initid=temnum;
    for(i=0;i<4;i++)
    {
      temseq[i]=string[22+i];
    }
    temseqi=atoi(temseq); 
    sltbrgl[numsltbrg].initindex=temseqi;
    //end
    for(i=0;i<4;i++)
    {
      sltbrgl[numsltbrg].endatom[i]=string[42+i];
    }
    sltbrgl[numsltbrg].endchainid=string[51];
    for(i=0;i<3;i++)
    {
      temami[i]=string[47+i];
    }
    temnum=aminoid(temami);
    sltbrgl[numsltbrg].endid=temnum;
    for(i=0;i<4;i++)
    {
      temseq[i]=string[52+i];
    }
    temseqi=atoi(temseq);
    sltbrgl[numsltbrg].endindex=temseqi;

    numsltbrg++;
    return 11;
  }
  //CISPEP   2 THR D   92    PRO D   93          0       359.80		
  else if (string[0]=='C' && string[1]=='I' && string[2]=='S' && string[3]=='P' && string[4]=='E' && string[5]=='P') 
  {
    if(numcispep==0)
    {
      cispepl=new cispeplink[alloccispep];
    }
    else if(numcispep==alloccispep-1)
    {
      alloccispep*=2;
      cispepl=(cispeplink*)realloc(cispepl, alloccispep*sizeof(cispeplink)); 		
    }
    //start
    cispepl[numcispep].initchainid=string[15];
    char temnumc[5],temami[3];
    temnumc[4]='\0';
    int temnumi,temamii;
    for(i=0;i<3;i++)
    {
      temami[i]=string[11+i];
    }
    temamii=aminoid(temami);
    cispepl[numcispep].initid=temamii;
    for(i=0;i<4;i++)
    {
      temnumc[i]=string[17+i];
    }
    temnumi=atoi(temnumc);
    cispepl[numcispep].initindex=temnumi;
    //end
    cispepl[numcispep].endchainid=string[29];
    for(i=0;i<3;i++) 
    {
      temami[i]=string[25+i];
    }
    temamii=aminoid(temami);
    cispepl[numcispep].endid=temamii;
    for(i=0;i<4;i++)
    {
      temnumc[i]=string[31+i];
    }
    temnumi=atoi(temnumc);
    cispepl[numcispep].endindex=temnumi;

    numcispep++;
    return 12;
  }
	
  //SITE     3 PAA  9 GLN A 231                                             8ATC 385
  //SITE     1 ZNB  4 CYS B 109  CYS B 114  CYS B 138  CYS B 141            8ATC 386		
  else if (string[0]=='S' && string[1]=='I' && string[2]=='T' && string[3]=='E' && string[4]==' ' ) 
  {
    if(numsite==0)
    {
      sitel=new sitelink[allocsite];
    }
    else if(numsite==allocsite-1)
    {
      allocsite*=2;
      sitel=(sitelink*)realloc(sitel, allocsite*sizeof(sitelink)); 		
    }
    //start 
    for(i=0;i<3;i++)
    {
      sitel[numsite].sitename[i]=string[11+i];
    }
    char temnumc[3];
    temnumc[2]='\0';
    int temnum;
    for(i=0;i<2;i++)
    {
      temnumc[i]=string[15+i];
    }
    temnum=atoi(temnumc);
    sitel[numsite].sitenum=temnum;
    char temami[3],temseq[5];
    temseq[4]='\0';
    int temamii,temseqi;
    for(i=0;i<4;i++)
    {
      for(j=0;j<3;j++)
      {
        temami[j]=string[18+11*i+j];
      }
      temamii=aminoid(temami);
      sitel[numsite].aminoid[i]=temamii;
      sitel[numsite].chainid[i]=string[22+11*i];
      for(j=0;j<4;j++)
      {
        temseq[j]=string[23+11*i+j];
      }
      temseqi=atoi(temseq);
      sitel[numsite].amiindex[i]=temseqi;
    }
			
    numsite++;
    return 13;
  }
  //MODEL        1		
  else if (string[0]=='M' && string[1]=='O' && string[2]=='D' && string[3]=='E' && string[4]=='L' ) 
  {
    //alloc space
    char tempnum[5];
    tempnum[4]='\0';
    int tempnumi;
    for(i=0;i<4;i++)
    { 
      tempnum[i]=string[10+i];
    }
    tempnumi=atoi(tempnum);
		
    if(numpromod==0 && flagfirstmodel)
    {
      flagfirstmodel=false;
      return 42;
    }
    if(numpromod==allocpromod-1)
    {
      allocpromod*=2;
      promod=(promodel*)realloc(promod, allocpromod*sizeof(promodel));
      for(i=numpromod+1;i<allocpromod;i++)
      {
        promod[i].nchain=0;
        for(int j=0;j<maxnumchains;j++)
        {
          promod[i].procha[j].chainseg.init=-1;
          promod[i].procha[j].chainseg.term=-1;
        }
      }
    }
    //end of model
    if(promod[numpromod].procha[promod[numpromod].nchain].chainseg.term==-1)
    {
      promod[numpromod].procha[promod[numpromod].nchain].chainseg.term=numproseq-1;
    }
    flagbeginmodel=true;
    //only here change numpromod
    numpromod++;
    return 14;
  }
  else if (string[0]=='E' && string[1]=='N' && string[2]=='D' && string[3]=='M' && string[4]=='D' && string[5]=='L') 
  {
    //do nothing
    char tempnum[6],tempami[3],tempindex[5],tempid;
    tempnum[5]='\0';
    //tempami[3]='\0';
    tempindex[4]='\0';
    int tempnumi,tempindexi;
    for(i=0;i<5;i++)
    {
      tempnum[i]=string[6+i];
    }
    tempnumi=atoi(tempnum);
    for(i=0;i<3;i++)
    {
      tempami[i]=string[17+i];
    }
    tempid=string[21];
    for(i=0;i<4;i++)
    {
      tempindex[i]=string[22+i];
    }
	
    tempindexi=atoi(tempindex);
    //end of model
    if(promod[numpromod].procha[promod[numpromod].nchain].chainseg.term==-1)
    {
      promod[numpromod].procha[promod[numpromod].nchain].chainseg.term=numproseq-1;
    }
    return 15;
  }
  //HETATM 1415  O2  BLE P   1      13.775  30.147  14.862  1.09 20.95           O
  //TER    1416      BLE P   1  
  else if (string[0]=='T' && string[1]=='E' && string[2]=='R') 
  {
    //end a chain
    flagter=true;
    return 16;
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //    ATOM   1660  CE2 TYR B 205      43.549  -4.115  45.779  1.00 19.60           C 
  //    HETATM 2292 2H   HOH   550      23.244   1.464  47.034  1.00  0.00           H 
  else if((string[0]=='A' && string[1]=='T' && string[2]=='O' && string[3]=='M' && string[4]==' ' ) ||
	  (string[0]=='H' && string[1]=='E' && string[2]=='T' && string[3]=='A' && string[4]=='T' && string[5]=='M'))
  {
    if (string[16]!=' ' && string[16]!='A' && numproseq>0)
    {	
      if(proseq[numproseq-1].detailtype[0]==string[12] &&proseq[numproseq-1].detailtype[1]==string[13] &&
         proseq[numproseq-1].detailtype[2]==string[14] &&proseq[numproseq-1].detailtype[3]==string[15] )
      {
        return 41;		
      }			
    }		
    if(numproseq==0)
    {
      proseq=new atom[allocproseq];
    }
    else if(numproseq==allocproseq-1)
    {
      allocproseq*=2;
      proseq=(atom*)realloc(proseq, allocproseq*sizeof(atom)); 		
    }
    //begin
    if(string[0]=='A' && string[1]=='T' && string[2]=='O' && string[3]=='M')
    {
      proseq[numproseq].simpletype=1;
    }
    else if(string[0]=='H' && string[1]=='E' && string[2]=='T' && string[3]=='A' && string[4]=='T' && string[5]=='M')
    {
      proseq[numproseq].simpletype=2;
    }
    char tempseq[6];
    tempseq[5]='\0';
    //indexnumber
    i=6,j=0;
    for(i=6;i<11;i++)
    {
      tempseq[j++]=string[i];
    }
    proseq[numproseq].seqno=atoi(tempseq);
    //type
    for(j=0;j<4;j++)
    {
      proseq[numproseq].detailtype[j]=' ';
    }
    //start none empty	
    for(i=12;i<16;i++)
    {
      if(string[i]!=' ')
      {
        break;
      }
    }
    proseq[numproseq].detailtype[4]=i-12;
    i=12;j=0;
    while(i<16)
    {
      if(string[i]!=' ')
      {
        proseq[numproseq].detailtype[j++]=string[i];
      }
      i++;
    }	
    proseq[numproseq].detailtype[5]=j;
    proseq[numproseq].detail=-1;
    if(proseq[numproseq].detailtype[0]=='C')
    {  
      if(proseq[numproseq].detailtype[1]=='A' && proseq[numproseq].detailtype[2]==' ')
      {
        proseq[numproseq].detail=0;//ca		
      }       
      else if(proseq[numproseq].detailtype[1]==' ')
      {
        proseq[numproseq].detail=1;//c
      }
      else proseq[numproseq].detail=7;//other cb
    }
    else if(proseq[numproseq].detailtype[0]=='O')
    {
      if(proseq[numproseq].detailtype[1]==' ')
      {
        proseq[numproseq].detail=3;
      }
      else
      {
        proseq[numproseq].detail=11;
      }
    }
    else if(proseq[numproseq].detailtype[0]=='N')
    {
      if(proseq[numproseq].detailtype[1]==' ')
      {
        proseq[numproseq].detail=2;
      }
      else 
      {
        proseq[numproseq].detail=8;//other ne
      }        
    }
    else if(proseq[numproseq].detailtype[0]=='S')
    {
      proseq[numproseq].detail=4;
    }
    else if(proseq[numproseq].detailtype[0]=='P')
    {
      proseq[numproseq].detail=6;
    }
    else if(proseq[numproseq].detailtype[0]=='H' ||
	    proseq[numproseq].detailtype[1]=='H' ||
            string[13]=='H')
    {
      if(proseq[numproseq].detailtype[0]=='H' && proseq[numproseq].detailtype[1]==' ')
        proseq[numproseq].detail=5;
      else proseq[numproseq].detail=12;
    }
    else if(string[13]=='F' && string[14]=='E')
    {
      proseq[numproseq].detail=9;
    }
    else
    {
      proseq[numproseq].detail=10;
    }
    proseq[numproseq].alt=string[16];
    //amino
    i=17,j=0;
    for(i=17;i<20;i++)
    {
      proseq[numproseq].residue[j++]=string[i];
    }
    proseq[numproseq].residueid=aminoid(proseq[numproseq].residue);
    //chainid
    i=21;
    proseq[numproseq].chainid=string[i];	
    //seq num
    char tempresnum[5]; 
    tempresnum[4]='\0';//important added
    i=22,j=0;
    for(i=22;i<26;i++)
    {
      tempresnum[j++]=string[i];
    }
    proseq[numproseq].resno =atoi(tempresnum);
    proseq[numproseq].ins=string[26];
    //sscanf(string+22,"%d",&proseq[numproseq].resno); //both are ok
    char tempx[8],tempy[8],tempz[8];
    //x
    i=30;j=0;
    for(i=30;i<38;i++)
    {
      tempx[j++]=string[i];
    }
    proseq[numproseq].x=float(atof(tempx));   
    //y
    i=38;j=0;
    for(i=38;i<46;i++)
    {
      tempy[j++]=string[i];
    }
    proseq[numproseq].y=float(atof(tempy));  
    //z
    i=46;j=0;
    for(i=46;i<54;i++)
    {
      tempz[j++]=string[i];
    }
    proseq[numproseq].z=float(atof(tempz));  
    //sscanf(string+30,"%f",&proseq[numproseq].x);
    //sscanf(string+38,"%f",&proseq[numproseq].y);
    //sscanf(string+46,"%f",&proseq[numproseq].z);
    //sscanf(string+54,"%f",&proseq[numproseq].occu);
    //sscanf(string+60,"%f",&proseq[numproseq].tempe);

    //deal with promod[numpromod]
    if(promod[numpromod].nchain==0 && flagbeginmodel)
    {
      promod[numpromod].procha[promod[numpromod].nchain].chainid=string[21];
      promod[numpromod].procha[promod[numpromod].nchain].chainseg.init=numproseq;	
      flagbeginmodel=false;
    }
    else if((flagter) || (string[21]!=proseq[numproseq-1].chainid &&
	    !(string[0]=='H' && string[1]=='E' && string[2]=='T' && string[3]=='A' && 
            proseq[numproseq-1].simpletype==2))//3d5a
            || (string[21]!=promod[numpromod].procha[promod[numpromod].nchain].chainid && 
            string[0]=='A' && string[1]=='T' && string[2]=='O' && string[3]=='M')//3fi0
    )
    {
      promod[numpromod].procha[promod[numpromod].nchain].chainseg.term=numproseq-1;
      promod[numpromod].nchain=promod[numpromod].nchain+1;
      promod[numpromod].procha[promod[numpromod].nchain].chainid=string[21];
      promod[numpromod].procha[promod[numpromod].nchain].chainseg.init=numproseq;
    }	
    if(flagter)
    {
      flagter=false;
    } 
    numproseq++;

    return 17;
  }
	
  return 0;
}

bool ParsePDB::loadca2ncbins(
  char *fileName
)
{
  if(cancbins[0]) return true;
  FILE *file;
  file=fopen(fileName,"rb");
  if(!file)
  {
    printf("no ca2nc file %s\n",fileName);
    return false;
  }
  printf("loading ca2nc bins\n");
  int i,j,k;
  for(i=0;i<13;i++)
  {
    if(!cancbins[i])
    {
      cancbins[i]=new float**[105];
      for(j=0;j<105;j++)
      {
        cancbins[i][j]=new float *[105];
        for(k=0;k<105;k++)
        {
          cancbins[i][j][k]=new float[160];
        }
      }
    }
  }
  for(i=0;i<13;i++)
  {
    if(i%2==0) continue;
    for(j=0;j<105;j++)
    {
      for(k=0;k<105;k++)
      {
        fread(cancbins[i][j][k],sizeof(float),160,file);
      }
    }
  }
  fclose(file);

  return true;
}

void ParsePDB::ca2nc(
  point3f *decstr,
  int numseq
)
{
  int j;
  point3d tp[5];
  point3s pout;
  int ii,jj,kk;
  BasicFunc bf;
  double tdist;
        
  ////////////////////////////head
  ii=26;jj=26;kk=21;
  tp[1]=bf.setv(decstr[0].x,decstr[0].y,decstr[0].z);
  tp[2]=bf.setv(decstr[1].x,decstr[1].y,decstr[1].z);
  tp[3]=bf.setv(decstr[2].x,decstr[2].y,decstr[2].z);
  bf.tor2pos22(tp[3].x,tp[3].y,tp[3].z,
               tp[2].x,tp[2].y,tp[2].z,
               tp[1].x,tp[1].y,tp[1].z,
               float(kk)/25.0f,3.813f,float(ii+53.0)/50.0f,
               &pout.x,&pout.y,&pout.z);
  tp[0]=bf.setv(pout.x,pout.y,pout.z);
  bf.tor2pos22(tp[0].x,tp[0].y,tp[0].z,
               tp[2].x,tp[2].y,tp[2].z,
               tp[1].x,tp[1].y,tp[1].z,
               cancbins[5][ii][jj][kk]*raddeg,cancbins[1][ii][jj][kk],cancbins[3][ii][jj][kk],
               &pout.x,&pout.y,&pout.z);
  decstr[0].ptc.x=pout.x;
  decstr[0].ptc.y=pout.y;
  decstr[0].ptc.z=pout.z;

  bf.tor2pos22(tp[3].x,tp[3].y,tp[3].z,
               tp[1].x,tp[1].y,tp[1].z,
               tp[2].x,tp[2].y,tp[2].z,
               cancbins[11][ii][jj][kk]*raddeg,cancbins[7][ii][jj][kk],cancbins[9][ii][jj][kk],
               &pout.x,&pout.y,&pout.z);
  decstr[1].ptn.x=pout.x;
  decstr[1].ptn.y=pout.y;
  decstr[1].ptn.z=pout.z;

  bf.tor2pos22(decstr[1].ptn.x,decstr[1].ptn.y,decstr[1].ptn.z,
               decstr[0].ptc.x,decstr[0].ptc.y,decstr[0].ptc.z,
               tp[1].x,tp[1].y,tp[1].z,
               180.0f*raddeg,1.460f,111.008f*raddeg,
               &pout.x,&pout.y,&pout.z);
  decstr[0].ptn.x=pout.x;
  decstr[0].ptn.y=pout.y;
  decstr[0].ptn.z=pout.z;

  ////////////////////////////////////tail
  tp[0]=bf.setv(decstr[numseq-3].x,decstr[numseq-3].y,decstr[numseq-3].z);
  tp[1]=bf.setv(decstr[numseq-2].x,decstr[numseq-2].y,decstr[numseq-2].z);
  tp[2]=bf.setv(decstr[numseq-1].x,decstr[numseq-1].y,decstr[numseq-1].z);
  bf.tor2pos22(tp[0].x,tp[0].y,tp[0].z,
               tp[1].x,tp[1].y,tp[1].z,
               tp[2].x,tp[2].y,tp[2].z,
               float(kk)/25.0f,3.813f,float(ii+53.0)/50.0f,
               &pout.x,&pout.y,&pout.z);
  tp[3]=bf.setv(pout.x,pout.y,pout.z);
  bf.tor2pos22(tp[0].x,tp[0].y,tp[0].z,
               tp[2].x,tp[2].y,tp[2].z,
               tp[1].x,tp[1].y,tp[1].z,
               cancbins[5][ii][jj][kk]*raddeg,cancbins[1][ii][jj][kk],cancbins[3][ii][jj][kk],
               &pout.x,&pout.y,&pout.z);
  decstr[numseq-2].ptc.x=pout.x;
  decstr[numseq-2].ptc.y=pout.y;
  decstr[numseq-2].ptc.z=pout.z;

  bf.tor2pos22(tp[3].x,tp[3].y,tp[3].z,
               tp[1].x,tp[1].y,tp[1].z,
               tp[2].x,tp[2].y,tp[2].z,
               cancbins[11][ii][jj][kk]*raddeg,cancbins[7][ii][jj][kk],cancbins[9][ii][jj][kk],
               &pout.x,&pout.y,&pout.z);
  decstr[numseq-1].ptn.x=pout.x;
  decstr[numseq-1].ptn.y=pout.y;
  decstr[numseq-1].ptn.z=pout.z;

  bf.tor2pos22(decstr[numseq-2].ptc.x,decstr[numseq-2].ptc.y,decstr[numseq-2].ptc.z,
                decstr[numseq-1].ptn.x,decstr[numseq-1].ptn.y,decstr[numseq-1].ptn.z,
                tp[2].x,tp[2].y,tp[2].z,
                180.0f*raddeg,1.525f,111.008f*raddeg,
                &pout.x,&pout.y,&pout.z);
  decstr[numseq-1].ptc.x=pout.x;
  decstr[numseq-1].ptc.y=pout.y;
  decstr[numseq-1].ptc.z=pout.z;
  for(j=0;j<numseq-3;j++)
  {
    tp[0]=bf.setv(decstr[j].x,decstr[j].y,decstr[j].z);
    tp[1]=bf.setv(decstr[j+1].x,decstr[j+1].y,decstr[j+1].z);
    tp[2]=bf.setv(decstr[j+2].x,decstr[j+2].y,decstr[j+2].z);
    tp[3]=bf.setv(decstr[j+3].x,decstr[j+3].y,decstr[j+3].z);
    tdist=bf.angv(bf.minu(tp[0],tp[1]),bf.minu(tp[2],tp[1]));
    ii=int(tdist*50.0)-53;
    if(ii<0) ii=0;
    else if(ii>102) ii=102;

    tdist=bf.angv(bf.minu(tp[3],tp[2]),bf.minu(tp[1],tp[2]));
    jj=int(tdist*50.0)-53;
    if(jj<0) jj=0;
    else if(jj>102) jj=102;

    tdist=bf.phi(tp[0].x,tp[0].y,tp[0].z,
                 tp[1].x,tp[1].y,tp[1].z,
                 tp[2].x,tp[2].y,tp[2].z,
                 tp[3].x,tp[3].y,tp[3].z)*raddeg;
    kk=int(tdist*25.0);

    bf.tor2pos22(tp[0].x,tp[0].y,tp[0].z,
                 tp[2].x,tp[2].y,tp[2].z,
                 tp[1].x,tp[1].y,tp[1].z,
                 cancbins[5][ii][jj][kk]*raddeg,cancbins[1][ii][jj][kk],cancbins[3][ii][jj][kk],
                 &pout.x,&pout.y,&pout.z);
    decstr[j+1].ptc.x=pout.x;
    decstr[j+1].ptc.y=pout.y;
    decstr[j+1].ptc.z=pout.z;

    bf.tor2pos22(tp[3].x,tp[3].y,tp[3].z,
                 tp[1].x,tp[1].y,tp[1].z,
                 tp[2].x,tp[2].y,tp[2].z,
                 cancbins[11][ii][jj][kk]*raddeg,cancbins[7][ii][jj][kk],cancbins[9][ii][jj][kk],
                 &pout.x,&pout.y,&pout.z);
    bf.tor2pos22(decstr[j+1].ptc.x,decstr[j+1].ptc.y,decstr[j+1].ptc.z,
                 tp[1].x,tp[1].y,tp[1].z,
                 tp[2].x,tp[2].y,tp[2].z,
                 180.0f*raddeg,1.460f,0.275f,
                 &pout.x,&pout.y,&pout.z);
    decstr[j+2].ptn.x=pout.x;
    decstr[j+2].ptn.y=pout.y;
    decstr[j+2].ptn.z=pout.z;
  }
  //get tor[0]
  int i;
  for(i=1;i<numseq;i++)
  {
    decstr[i].tor[0]=bf.phi(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,
                            decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
                            decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                            decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z);
  }
  decstr[0].tor[0]=180;

  //add o
  bool flagts;
  point3s pn;
  for(i=1;i<numseq;i++)
  {
    flagts=bf.tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                        decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
                        decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                        179.6715f*raddeg,1.229f,2.0961f,
                        &pn.x,&pn.y,&pn.z);
    decstr[i-1].pto.x=pn.x;decstr[i-1].pto.y=pn.y;decstr[i-1].pto.z=pn.z;
  }
  i=numseq;
  flagts=bf.tor2pos22(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,
                      decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
                      decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                      0.0f,1.2439f,2.0855f,
                      &pn.x,&pn.y,&pn.z);
  decstr[i-1].pto.x=pn.x;decstr[i-1].pto.y=pn.y;decstr[i-1].pto.z=pn.z;
}

double ParsePDB::lsfrmsdp(
  double pin[], 
  double pout[], 
  int np,
  double delta
)
{
  int i;
  double tsum=0;
  double sf=2*PI/delta;
  double tdist,tdist2;
  for(i=0;i<np;i++)
  {
    tdist=acos(cos((pin[i]-pout[i])*sf))/sf;
    tdist2=acos(cos((pin[i+np]-pout[i+np])*sf))/sf;
    tsum+=tdist*tdist+tdist2*tdist2;
  }
  tsum/=double(np);
  tsum=sqrt(tsum)/180.0;

  return tsum;
}

double ParsePDB::lsfrmsd(
  double pin[], 
  double pout[], 
  int np, 
  double pmat[], 
  double ptrans[]
)
{
  int i;
  transpara(pin, pout, np, pmat, ptrans);
  double *ptmp=new double[3*np];
  double rmsd=0;
  BasicFunc bf;
  bf.trmul(pmat,pout,3,3,np,ptmp);
  for(i=0;i<np;i++)
  {
    ptmp[i]+=ptrans[0];
    ptmp[np+i]+=ptrans[1];
    ptmp[2*np+i]+=ptrans[2];
  }
  for(i=0;i<np;i++)
  {
    rmsd+=(ptmp[i]-pin[i])*(ptmp[i]-pin[i])+(ptmp[np+i]-pin[np+i])*(ptmp[np+i]-pin[np+i])+
          (ptmp[2*np+i]-pin[2*np+i])*(ptmp[2*np+i]-pin[2*np+i]);
  }
  rmsd=sqrt(rmsd/double(np));
  delete[]ptmp;

  return rmsd;
}

int ParsePDB::transpara( 
  double pin[], 
  double pout[], 
  int np, 
  double pmat[], 
  double ptrans[]
)
{
  int i,j;
  int rval;
  BasicFunc bf;
  //center
  double cpin[3],cpout[3];
  for(i=0;i<3;i++)
  {       
    cpin[i]=0;
    cpout[i]=0;
  }
  for(i=0;i<np;i++)
  {       
    for(j=0;j<3;j++)
    {       
      cpin[j]+=pin[j*np+i];
      cpout[j]+=pout[j*np+i];
    }
  }
  for(i=0;i<3;i++)
  {       
    cpin[i]/=double(np);
    cpout[i]/=double(np);
  }
  for(i=0;i<np;i++)
  {       
    for(j=0;j<3;j++)
    {       
      pin[j*np+i]-=cpin[j];
      pout[j*np+i]-=cpout[j];
    }
  }
  ///start
  double *tpin=new double[3*np];
  bf.tranmat(pin,3,np,tpin);
  double hmat[9],umat[9],vmat[9],vt[9],ut[9];
  double eps=0.000001;

  bf.trmul(pout,tpin,3,np,3,hmat);
  rval=bf.muav(hmat,3,3,umat,vmat,eps,4);
  if(rval==-1)
  {       
    delete[]tpin;
    return rval;
  }
  bf.tranmat(umat,3,3,ut);
  bf.tranmat(vmat,3,3,vt);
  double sign=bf.sdet(hmat,3);// Kabsch algorithm
  int bb;
  if(sign>=0)
  {
    bb=1;
  }
  else
  {
    bb=-1;
  }
  double bmat[9]={1,0,0,0,1,0,0,0,bb};
  double tmpmat[9];
  bf.trmul(vt,bmat,3,3,3,tmpmat);
  bf.trmul(tmpmat,ut,3,3,3,pmat);
  bf.trmul(pmat,cpout,3,3,1,ptrans);
  for(i=0;i<3;i++)
  {
    ptrans[i]=cpin[i]-ptrans[i];
  }
  for(i=0;i<np;i++)
  {
    for(j=0;j<3;j++)
    {
      pin[j*np+i]+=cpin[j];
      pout[j*np+i]+=cpout[j];
    }
  }
  delete[]tpin;
  
  return rval;
}

int ParsePDB::aminoid(
  char *aminoname
)
{
  int i;
  //empty
  if(aminoname[0]==' ' && aminoname[1]==' ' && aminoname[2]==' ')
  {
    //printf("empty amino name\n");
    return 23;
  }
  //one in 26
  for(i=0;i<26;i++)
  {
    if(aminoname[0]==aad3[i][0] && aminoname[1]==aad3[i][1] && aminoname[2]==aad3[i][2])
    {
      if(i==24)
      {
        i=23;
      }
      return i;
    }
  }
  for(i=0;i<8;i++)
  {
    if(aminoname[0]==dnarnares[i][0] && aminoname[1]==dnarnares[i][1] && aminoname[2]==dnarnares[i][2])
    {
      return -i-1;
    }
  }
  //unknown
  //printf("unknown amino name %c%c%c\n",aminoname[0],aminoname[1],aminoname[2]);
  if(aminoname[0]=='M' && aminoname[1]=='S' && aminoname[2]=='E') return 10;
  
  return 23;
}

