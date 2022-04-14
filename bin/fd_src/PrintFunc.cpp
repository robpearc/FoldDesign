///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

/* PrintFunc - auxiliary functions for printing purposes */
#include "PrintFunc.h"
#include "BasicFunc.h"

PrintFunc::PrintFunc(){}
PrintFunc::~PrintFunc(){}

/* print secondary structure of a chain 'decstr' with 'resnum' residues */
void PrintFunc::printSS(
  const point3f *decstr,
  const int resnum
)
{
  string txt="ss3\t";
  for (int i=0;i<resnum;i++) txt+=decstr[i].ss3;
  txt+="\nssm\t";
  for (int i=0;i<resnum;i++) txt+=decstr[i].ssm;
  txt+="\nstype\t";
  for (int i=0;i<resnum;i++) txt+=decstr[i].stype;
  cout<<txt<<endl;
  txt.clear();
}

/* write decoy conformations in SPICKER format
 * denergy - total energy of decoy
 * iind  - model index
 * ftype - whether create new file (0) or append to old file (1) */
bool PrintFunc::writetra(
  point3f *decstr, 
  int numseq, 
  char *outname,
  int iind, 
  double denergy, 
  int ftype
)
{
  int i;
  FILE *fp;

  if(ftype==0) fp=fopen(outname,"wt");
  else if(ftype==1) fp=fopen(outname,"a+");

  if(!fp)
  {
    printf("cannot write file %s\n",outname);
    return false;
  }
  
  fprintf(fp,"%d %.1f %d %d\n",numseq,denergy,iind,iind);
  for(i=0;i<numseq;i++)
    fprintf(fp,"%10.3f %10.3f %10.3f\n",decstr[i].x,decstr[i].y,decstr[i].z);
  fclose(fp);

  return true;
}

/* similar to writetra but write PDB file instead 
 * bbres    - list of residue index
 * eneterms - list of energy terms
 * numene   - number of energy terms
 * iind     - model index
 * atype    - atom type to print. e.g. atype=7 for printing N, CA, C
 *            +  1: CA
 *            +  2: N
 *            +  4: C
 *            +  8: O
 *            + 16: H
 *            + 32: CB
 *            + 64: SC
 *            +128: CT
 * ftype    - 0: overwrite existing file; 1: append to existing file
 */
bool PrintFunc::writepdb(
  point3f *decstr,
  int numseq,
  int *bbres,
  const char *outname,
  int iind, 
  double denergy, 
  double *eneterms,
  int numene,
  int atype,
  int ftype)
{
  FILE *fp;
  if(ftype==0) fp=fopen(outname,"wt");
  else if(ftype==1) fp=fopen(outname,"a+");
  if(!fp)
  {
    printf("cannot write file %s\n",outname);
    return false;
  }

  bool parse_CT=(atype>=128); atype%=128;
  bool parse_SC=(atype>= 64); atype%= 64;
  bool parse_CB=(atype>= 32); atype%= 32;
  bool parse_H =(atype>= 16); atype%= 16;
  bool parse_O =(atype>=  8); atype%=  8;
  bool parse_C =(atype>=  4); atype%=  4;
  bool parse_N =(atype>=  2); atype%=  2;
  bool parse_CA=(atype>=  1);

  fprintf(fp,"MODEL     %4d %4d %10.3f",iind,numseq,denergy);
  int i;
  for(i=0;i<numene;i++) fprintf(fp," %10.4f",eneterms[i]);
  fprintf(fp,"\n");

  int indatom=1;
  for (int r=0;r<numseq;r++)
  {
    const char *resn;
    for(i=0;i<26;i++)
    {
      if(decstr[r].aaa==aad1[i]) 
      {
        resn=aad3[i];
        break;
      }
    }
    if(i==26) resn="UNK";
    if(parse_N)
      fprintf(fp, "ATOM  %5d %s %s  %4d    %8.3f%8.3f%8.3f\n",
              indatom++, " N  ", resn, bbres[r],
              decstr[r].ptn.x, decstr[r].ptn.y, decstr[r].ptn.z);
    if(parse_CA)
      fprintf(fp, "ATOM  %5d %s %s  %4d    %8.3f%8.3f%8.3f\n",
              indatom++, " CA ", resn, bbres[r],
              decstr[r].x,     decstr[r].y,     decstr[r].z);
    if(parse_C)
      fprintf(fp, "ATOM  %5d %s %s  %4d    %8.3f%8.3f%8.3f\n",
              indatom++, " C  ", resn, bbres[r],
              decstr[r].ptc.x, decstr[r].ptc.y, decstr[r].ptc.z);
    if(parse_O)
      fprintf(fp, "ATOM  %5d %s %s  %4d    %8.3f%8.3f%8.3f\n",
              indatom++, " O  ", resn, bbres[r],
              decstr[r].pto.x, decstr[r].pto.y, decstr[r].pto.z);
    if(parse_CB)
      fprintf(fp, "ATOM  %5d %s %s  %4d    %8.3f%8.3f%8.3f\n",
              indatom++, " CB ", resn, bbres[r],
              decstr[r].ptb.x, decstr[r].ptb.y, decstr[r].ptb.z);
    if(parse_H)
      fprintf(fp, "ATOM  %5d %s %s  %4d    %8.3f%8.3f%8.3f\n",
              indatom++, " H  ", resn, bbres[r],
              decstr[r].pth.x, decstr[r].pth.y, decstr[r].pth.z);
    if(parse_SC)
      fprintf(fp, "ATOM  %5d %s %s  %4d    %8.3f%8.3f%8.3f\n",
              indatom++, " SC ", resn, bbres[r],
              decstr[r].ptsg.x,decstr[r].ptsg.y,decstr[r].ptsg.z);
    if(parse_CT)
      fprintf(fp, "ATOM  %5d %s %s  %4d    %8.3f%8.3f%8.3f\n",
              indatom++, " CT ", resn, bbres[r],
              decstr[r].ptg.x, decstr[r].ptg.y, decstr[r].ptg.z);
  }
  fprintf(fp,"TER\nENDMDL\n");
  fclose(fp);

  return true;
}

bool PrintFunc::writepdb(
  point3f *decstr,
  int numseq,
  char *outname,
  int iind,
  double denergy,
  double *eneterms,
  int numene,
  int atype,
  int ftype
)
{
  int i;
  FILE *file;
  BasicFunc bf;
  int indaa,indaaf,indatom=1;
  if(ftype==0)
  {
    file=fopen(outname,"wt");
  }
  else if(ftype==1)
  {
    file=fopen(outname,"a+");
  }
  if(!file)
  {
    printf("cannot write file %s\n",outname);
    return false;
  }
  if(atype!=71)
  {
    fprintf(file,"MODEL     %4d %4d %10.3f",iind,numseq,denergy);
    for(i=0;i<numene;i++) fprintf(file," %10.4f",eneterms[i]);
    fprintf(file,"\n");
  }
  if(atype==8)
  {
    for(i=0;i<numseq;i++)
    {
      indaa=bf.aminoid(decstr[i].aaa);
      if(indaa<0 || indaa>19) indaa=5;
      indaaf=bf.aminoid(decstr[i].residueid);
      if(indaaf<0 || indaaf>19) indaaf=5;
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','N',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],decstr[i].ss3,
              decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].x,decstr[i].y,decstr[i].z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','O',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='G')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','C','B',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                i+1,decstr[i].ptb.x,decstr[i].ptb.y,decstr[i].ptb.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='A' && decstr[i].aaa!='G')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','S','G',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                i+1,decstr[i].ptsg.x,decstr[i].ptsg.y,decstr[i].ptsg.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='P')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','H',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                i+1,decstr[i].pth.x,decstr[i].pth.y,decstr[i].pth.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='G')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','H','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                i+1,decstr[i].ptha.x,decstr[i].ptha.y,decstr[i].ptha.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      else
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,'1','H','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                i+1,decstr[i].ptha.x,decstr[i].ptha.y,decstr[i].ptha.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
    }
    fprintf(file,"TER\nENDMDL\n");
  }
  else if(atype==71)
  {
    for(i=0;i<numseq;i++)
    {
      indaa=bf.aminoid(decstr[i].aaa);
      if(indaa<0 || indaa>19) indaa=5;
      indaaf=bf.aminoid(decstr[i].residueid);
      if(indaaf<0 || indaaf>19) indaaf=5;
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f  1.00  0.00      S_00 %3d\n",
              indatom++,' ','N',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,i+1);
      if(decstr[i].aaa!='P')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f  1.00  0.00      S_00 %3d\n",
                indatom++,' ','H',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                i+1,decstr[i].pth.x,decstr[i].pth.y,decstr[i].pth.z,i+1);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f  1.00  0.00      S_00 %3d\n",
              indatom++,' ','C','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].x,decstr[i].y,decstr[i].z,i+1);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f  1.00  0.00      S_00 %3d\n",
              indatom++,' ','H','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].ptha.x,decstr[i].ptha.y,decstr[i].ptha.z,i+1);
      if(decstr[i].aaa!='G')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f  1.00  0.00      S_00 %3d\n",
                indatom++,' ','C','B',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                i+1,decstr[i].ptb.x,decstr[i].ptb.y,decstr[i].ptb.z,i+1);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f  1.00  0.00      S_00 %3d\n",
              indatom++,' ','C',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,i+1);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f  1.00  0.00      S_00 %3d\n",
              indatom++,' ','O',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,i+1);
    }
    fprintf(file,"TER\n");
  }
  else if(atype==7)
  {
    for(i=0;i<numseq;i++)
    {
      indaa=bf.aminoid(decstr[i].aaa);
      if(indaa<0 || indaa>19) indaa=5;
      indaaf=bf.aminoid(decstr[i].residueid);
      if(indaaf<0 || indaaf>19) indaaf=5;
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','N',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].x,decstr[i].y,decstr[i].z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','O',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='G')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','C','B',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                i+1,decstr[i].ptb.x,decstr[i].ptb.y,decstr[i].ptb.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='G')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','S','G',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                i+1,decstr[i].ptsg.x,decstr[i].ptsg.y,decstr[i].ptsg.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='P')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','H',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                i+1,decstr[i].pth.x,decstr[i].pth.y,decstr[i].pth.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
    }
    fprintf(file,"TER\nENDMDL\n");
  }
  else if(atype==61)
  {
    for(i=0;i<numseq;i++)
    {
      indaa=bf.aminoid(decstr[i].aaa);
      if(indaa<0 || indaa>19) indaa=5;
      indaaf=bf.aminoid(decstr[i].residueid);
      if(indaaf<0 || indaaf>19) indaaf=5;
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','N',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].x,decstr[i].y,decstr[i].z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','O',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='G')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','C','B',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                i+1,decstr[i].ptb.x,decstr[i].ptb.y,decstr[i].ptb.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='G')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','H',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                i+1,decstr[i].pth.x,decstr[i].pth.y,decstr[i].pth.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
    }
    fprintf(file,"TER\nENDMDL\n");
  }
  else if(atype==6)
  {
    for(i=0;i<numseq;i++)
    {
      indaa=bf.aminoid(decstr[i].aaa);
      if(indaa<0 || indaa>19) indaa=5;
      indaaf=bf.aminoid(decstr[i].residueid);
      if(indaaf<0 || indaaf>19) indaaf=5;
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','N',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].x,decstr[i].y,decstr[i].z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','O',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='G')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','C','B',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                i+1,decstr[i].ptb.x,decstr[i].ptb.y,decstr[i].ptb.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='G')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','S','G',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                i+1,decstr[i].ptsg.x,decstr[i].ptsg.y,decstr[i].ptsg.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
    }
    fprintf(file,"TER\nENDMDL\n");
  }
  else if(atype==5)
  {
    for(i=0;i<numseq;i++)
    {
      indaa=bf.aminoid(decstr[i].aaa);
      if(indaa<0 || indaa>19) indaa=5;
      indaaf=bf.aminoid(decstr[i].residueid);
      if(indaaf<0 || indaaf>19) indaaf=5;
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','N',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].x,decstr[i].y,decstr[i].z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','O',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='G')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','C','B',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                i+1,decstr[i].ptb.x,decstr[i].ptb.y,decstr[i].ptb.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
    }
    fprintf(file,"TER\nENDMDL\n");
  }
  else if(atype==51)
  {
    for(i=0;i<numseq;i++)
    {
      indaa=bf.aminoid(decstr[i].aaa);
      if(indaa<0 || indaa>19) indaa=5;
      indaaf=bf.aminoid(decstr[i].residueid);
      if(indaaf<0 || indaaf>19) indaaf=5;
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','N',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].x,decstr[i].y,decstr[i].z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','O',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='G')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','H',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                i+1,decstr[i].pth.x,decstr[i].pth.y,decstr[i].pth.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
    }
    fprintf(file,"TER\nENDMDL\n");
  }
  else if(atype==4)
  {
    for(i=0;i<numseq;i++)
    {
      indaa=bf.aminoid(decstr[i].aaa);
      if(indaa<0 || indaa>19) indaa=5;
      indaaf=bf.aminoid(decstr[i].residueid);
      if(indaaf<0 || indaaf>19) indaaf=5;
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','N',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].x,decstr[i].y,decstr[i].z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','O',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='G')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','S','G',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                i+1,decstr[i].ptsg.x,decstr[i].ptsg.y,decstr[i].ptsg.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
    }
    fprintf(file,"TER\nENDMDL\n");
  }
  else if(atype==41)
  {
    for(i=0;i<numseq;i++)
    {
      indaa=bf.aminoid(decstr[i].aaa);
      if(indaa<0 || indaa>19) indaa=5;
      indaaf=bf.aminoid(decstr[i].residueid);
      if(indaaf<0 || indaaf>19) indaaf=5;
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','N',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].x,decstr[i].y,decstr[i].z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','O',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
    }
    fprintf(file,"TER\nENDMDL\n");
  }
  else if(atype==3)
  {
    for(i=0;i<numseq;i++)
    {
      indaa=bf.aminoid(decstr[i].aaa);
      if(indaa<0 || indaa>19) indaa=5;
      indaaf=bf.aminoid(decstr[i].residueid);
      if(indaaf<0 || indaaf>19) indaaf=5;
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','N',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].x,decstr[i].y,decstr[i].z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
    }
    fprintf(file,"TER\nENDMDL\n");
  }
  else if(atype==1)
  {
    for(i=0;i<numseq;i++)
    {
      indaa=bf.aminoid(decstr[i].aaa);
      if(indaa<0 || indaa>19) indaa=5;
      indaaf=bf.aminoid(decstr[i].residueid);
      if(indaaf<0 || indaaf>19) indaaf=5;
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].x,decstr[i].y,decstr[i].z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
    }
    fprintf(file,"TER\nENDMDL\n");
  }
  else if(atype==0)//combo
  {
    for(i=0;i<numseq;i++)
    {
      indaa=decstr[i].iaa;
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f\n",
              i+1,' ','C','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              i+1,decstr[i].x,decstr[i].y,decstr[i].z);
    }
    fprintf(file,"TER\nENDMDL\n");
  }
  else
  {
    printf("no option %d\n",atype);
  }
  fclose(file);
  return true;
}


bool PrintFunc::writepdb(
  point3f *decstr,
  int numseq,
  int *bbres,
  char *outname,
  int iind,
  double denergy,
  double *eneterms,
  int numene,
  int atype,
  int ftype
)
{
  int i;
  FILE *file;
  BasicFunc bf;
  int indaa,indaaf,indatom=1;
  if(ftype==0)
  {
    file=fopen(outname,"wt");
  }
  else if(ftype==1)
  {
    file=fopen(outname,"a+");
  }
  if(!file)
  {
    printf("cannot write file %s\n",outname);
    return false;
  }
  if(atype!=71)
  {
    fprintf(file,"MODEL     %4d %4d %10.3f",iind,numseq,denergy);
    for(i=0;i<numene;i++) fprintf(file," %10.4f",eneterms[i]);
    fprintf(file,"\n");
  }
  if(atype==8)
  {
    for(i=0;i<numseq;i++)
    {
      indaa=bf.aminoid(decstr[i].aaa);
      if(indaa<0 || indaa>19) indaa=5;
      indaaf=bf.aminoid(decstr[i].residueid);
      if(indaaf<0 || indaaf>19) indaaf=5;
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','N',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].x,decstr[i].y,decstr[i].z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','O',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='G')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','C','B',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                bbres[i],decstr[i].ptb.x,decstr[i].ptb.y,decstr[i].ptb.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='A' && decstr[i].aaa!='G')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','S','G',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                bbres[i],decstr[i].ptsg.x,decstr[i].ptsg.y,decstr[i].ptsg.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='P')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','H',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                bbres[i],decstr[i].pth.x,decstr[i].pth.y,decstr[i].pth.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='G')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','H','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                bbres[i],decstr[i].ptha.x,decstr[i].ptha.y,decstr[i].ptha.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      else
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,'1','H','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                bbres[i],decstr[i].ptha.x,decstr[i].ptha.y,decstr[i].ptha.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
    }
    fprintf(file,"TER\nENDMDL\n");
  }
  else if(atype==71)
  {
    for(i=0;i<numseq;i++)
    {
      indaa=bf.aminoid(decstr[i].aaa);
      if(indaa<0 || indaa>19) indaa=5;
      indaaf=bf.aminoid(decstr[i].residueid);
      if(indaaf<0 || indaaf>19) indaaf=5;
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f  1.00  0.00      S_00 %3d\n",
              indatom++,' ','N',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,i+1);
      if(decstr[i].aaa!='P')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f  1.00  0.00      S_00 %3d\n",
                indatom++,' ','H',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                bbres[i],decstr[i].pth.x,decstr[i].pth.y,decstr[i].pth.z,i+1);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f  1.00  0.00      S_00 %3d\n",
              indatom++,' ','C','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].x,decstr[i].y,decstr[i].z,i+1);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f  1.00  0.00      S_00 %3d\n",
              indatom++,' ','H','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].ptha.x,decstr[i].ptha.y,decstr[i].ptha.z,i+1);
      if(decstr[i].aaa!='G')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f  1.00  0.00      S_00 %3d\n",
                indatom++,' ','C','B',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                bbres[i],decstr[i].ptb.x,decstr[i].ptb.y,decstr[i].ptb.z,i+1);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f  1.00  0.00      S_00 %3d\n",
              indatom++,' ','C',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,i+1);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f  1.00  0.00      S_00 %3d\n",
              indatom++,' ','O',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,i+1);
    }
    fprintf(file,"TER\n");
  }
  else if(atype==7)
  {
    for(i=0;i<numseq;i++)
    {
      indaa=bf.aminoid(decstr[i].aaa);
      if(indaa<0 || indaa>19) indaa=5;
      indaaf=bf.aminoid(decstr[i].residueid);
      if(indaaf<0 || indaaf>19) indaaf=5;
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','N',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].x,decstr[i].y,decstr[i].z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','O',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='G')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','C','B',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                bbres[i],decstr[i].ptb.x,decstr[i].ptb.y,decstr[i].ptb.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='G')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','S','G',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                bbres[i],decstr[i].ptsg.x,decstr[i].ptsg.y,decstr[i].ptsg.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='P')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','H',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                bbres[i],decstr[i].pth.x,decstr[i].pth.y,decstr[i].pth.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
    }
    fprintf(file,"TER\nENDMDL\n");
  }
  else if(atype==61)
  {
    for(i=0;i<numseq;i++)
    {
      indaa=bf.aminoid(decstr[i].aaa);
      if(indaa<0 || indaa>19) indaa=5;
      indaaf=bf.aminoid(decstr[i].residueid);
      if(indaaf<0 || indaaf>19) indaaf=5;
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','N',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].x,decstr[i].y,decstr[i].z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','O',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='G')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','C','B',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                bbres[i],decstr[i].ptb.x,decstr[i].ptb.y,decstr[i].ptb.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='G')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','H',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                bbres[i],decstr[i].pth.x,decstr[i].pth.y,decstr[i].pth.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
    }
    fprintf(file,"TER\nENDMDL\n");
  }
  else if(atype==6)
  {
    for(i=0;i<numseq;i++)
    {
      indaa=bf.aminoid(decstr[i].aaa);
      if(indaa<0 || indaa>19) indaa=5;
      indaaf=bf.aminoid(decstr[i].residueid);
      if(indaaf<0 || indaaf>19) indaaf=5;
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','N',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].x,decstr[i].y,decstr[i].z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','O',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='G')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','C','B',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                bbres[i],decstr[i].ptb.x,decstr[i].ptb.y,decstr[i].ptb.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='G')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','S','G',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                bbres[i],decstr[i].ptsg.x,decstr[i].ptsg.y,decstr[i].ptsg.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
    }
    fprintf(file,"TER\nENDMDL\n");
  }
  else if(atype==5)
  {
    for(i=0;i<numseq;i++)
    {
      indaa=bf.aminoid(decstr[i].aaa);
      if(indaa<0 || indaa>19) indaa=5;
      indaaf=bf.aminoid(decstr[i].residueid);
      if(indaaf<0 || indaaf>19) indaaf=5;
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','N',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].x,decstr[i].y,decstr[i].z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','O',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='G')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','C','B',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                bbres[i],decstr[i].ptb.x,decstr[i].ptb.y,decstr[i].ptb.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
    }
    fprintf(file,"TER\nENDMDL\n");
  }
  else if(atype==51)
  {
    for(i=0;i<numseq;i++)
    {
      indaa=bf.aminoid(decstr[i].aaa);
      if(indaa<0 || indaa>19) indaa=5;
      indaaf=bf.aminoid(decstr[i].residueid);
      if(indaaf<0 || indaaf>19) indaaf=5;
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','N',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].x,decstr[i].y,decstr[i].z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','O',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='G')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','H',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                bbres[i],decstr[i].pth.x,decstr[i].pth.y,decstr[i].pth.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
    }
    fprintf(file,"TER\nENDMDL\n");
  }
  else if(atype==4)
  {
    for(i=0;i<numseq;i++)
    {
      indaa=bf.aminoid(decstr[i].aaa);
      if(indaa<0 || indaa>19) indaa=5;
      indaaf=bf.aminoid(decstr[i].residueid);
      if(indaaf<0 || indaaf>19) indaaf=5;
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','N',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].x,decstr[i].y,decstr[i].z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','O',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      if(decstr[i].aaa!='G')
        fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
                indatom++,' ','S','G',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
                bbres[i],decstr[i].ptsg.x,decstr[i].ptsg.y,decstr[i].ptsg.z,
                decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
                decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
                decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
    }
    fprintf(file,"TER\nENDMDL\n");
  }
  else if(atype==41)
  {
    for(i=0;i<numseq;i++)
    {
      indaa=bf.aminoid(decstr[i].aaa);
      if(indaa<0 || indaa>19) indaa=5;
      indaaf=bf.aminoid(decstr[i].residueid);
      if(indaaf<0 || indaaf>19) indaaf=5;
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','N',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].x,decstr[i].y,decstr[i].z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','O',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].pto.x,decstr[i].pto.y,decstr[i].pto.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
    }
    fprintf(file,"TER\nENDMDL\n");
  }
  else if(atype==3)
  {
    for(i=0;i<numseq;i++)
    {
      indaa=bf.aminoid(decstr[i].aaa);
      if(indaa<0 || indaa>19) indaa=5;
      indaaf=bf.aminoid(decstr[i].residueid);
      if(indaaf<0 || indaaf>19) indaaf=5;
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','N',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].x,decstr[i].y,decstr[i].z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C',' ',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
    }
    fprintf(file,"TER\nENDMDL\n");
  }
  else if(atype==1)
  {
    for(i=0;i<numseq;i++)
    {
      indaa=bf.aminoid(decstr[i].aaa);
      if(indaa<0 || indaa>19) indaa=5;
      indaaf=bf.aminoid(decstr[i].residueid);
      if(indaaf<0 || indaaf>19) indaaf=5;
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %4d %c%c%c %c%c%c%c%c %c%c%c\n",
              indatom++,' ','C','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].x,decstr[i].y,decstr[i].z,
              decstr[i].resind,aad3[indaaf][0], aad3[indaaf][1],aad3[indaaf][2],decstr[i].name[0],
              decstr[i].name[1],decstr[i].name[2],decstr[i].name[3],decstr[i].name[4],
              decstr[i].ss3,decstr[i].stype,decstr[i].ssm);
    }
    fprintf(file,"TER\nENDMDL\n");
  }
  else if(atype==0)//combo
  {
    for(i=0;i<numseq;i++)
    {
      indaa=decstr[i].iaa;
      fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f\n",
              i+1,' ','C','A',' ', aad3[indaa][0], aad3[indaa][1],aad3[indaa][2],' ',
              bbres[i],decstr[i].x,decstr[i].y,decstr[i].z);
    }
    fprintf(file,"TER\nENDMDL\n");
  }
  else
  {
    printf("no option %d\n",atype);
  }
  fclose(file);
  return true;
}

/* write decoy conformations in SPICKER format
 * denergy - total energy of decoy
 * iind  - model index
 * ftype - whether create new file (0) or append to old file (1) */
bool PrintFunc::writetradesign(
  point3f *decstr, 
  int numseq, 
  char *outname,
  int iind, 
  double denergy, 
  int ss_score, 
  int ftype
)
{
  int i;
  FILE *fp;

  if(ftype==0) fp=fopen(outname,"wt");
  else if(ftype==1) fp=fopen(outname,"a+");

  if(!fp)
  {
    printf("cannot write file %s\n",outname);
    return false;
  }
  
  fprintf(fp,"%d %.1f %d %d %d\n",numseq,denergy,ss_score,iind,iind);
  for(i=0;i<numseq;i++)
    fprintf(fp,"%10.3f %10.3f %10.3f\n",decstr[i].x,decstr[i].y,decstr[i].z);
  fclose(fp);

  return true;
}

