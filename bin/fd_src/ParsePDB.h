///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#if !defined(AFX_PARSEPDB_H__0A306B68_CEEB_4EEF_8D55_4CCE700E11A3__INCLUDED_)
#define AFX_PARSEPDB_H__0A306B68_CEEB_4EEF_8D55_4CCE700E11A3__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <vector>
#include "CommonPara.h"

class ParsePDB  
{
  public:
    ParsePDB();
    virtual ~ParsePDB();

    char pdbname[6];
    int headnumchain;//number of chains in header file
    char headidchain[maxnumchains];// id of chains
    int indexchain,numinchain[maxnumchains];// num of aminos in each chain
    char *aminoinchain[maxnumchains];// aminos in each chain
    paratrans pt;
    int numhelix,numsheet,numturn,numlink;
    int allochelix,allocsheet,allocturn,alloclink;
    helixlink *helixl;
    sheetlink *sheetl;
    turnlink *turnl;
    linklink *linkl;
    headinfo hinfo;
    int numssbond,numsltbrg,numcispep,numhydbnd,numsite;
    int allocssbond,allocsltbrg,alloccispep,allochydbnd,allocsite;
    ssbondlink *ssbondl;
    sltbrglink *sltbrgl;
    cispeplink *cispepl;
    hydbndlink *hydbndl;
    sitelink *sitel;
    int numpromod;
    int allocpromod;
    promodel *promod;
    bool flagfirstmodel;
    bool flagbeginmodel;
    int numproseq;
    int allocproseq;
    atom *proseq;
    bool flagter;
    int aminoid(char *aminoname);
    int aminoid(char aminoname);
    int dealline(char *string);
    void loadPDBChain(const vector<string>&pdb_lines,vector<point3f>&decstr_vec);
    int loadPDBChain(const char *pdbfile,vector<point3f>&decstr_vec);
    void loadPDBChain(const vector<string>&pdb_lines,point3f *decstr_vec);
    void loadPDBChainFull(const vector<string>&pdb_lines,point3ffull *decstr_vec);
    //void loadpdb_chainhbond(const vector<string>&pdb_lines,point3f *decstr_vec);
    bool loadpdb(char *pdbfile);
    int numbb,effbb;
    boneinfo *bb;
    void extractbb(int model, int chain,int type);
    double lsfrmsdp(double pin[],double pout[],int np,double delta);//phi psi pair
    double lsfrmsd(double pin[],double pout[],int np,double pmat[],double ptrans[]);
    int transpara(double pin[],double pout[],int np,double pmat[],double ptrans[]);
    int numhssp;
    hssp *hssps;
    int numsse;
    sseinfo *sses;//for dssp
    hybond *hbond;
    caneigh *cacont;
    betasheet *bsheet;
    int numbsheet;
    alphahelix *ahelix;
    int numahelix;
    betastrand *bstrand;
    int numbstrand;
    int numabss;
    abssinfo *abss;

    bool loadca2ncbins(char *fileName);
    void ca2nc(point3f *decstr,int numseq);
    float ***cancbins[13];
    int *bbid,*dsspid;
    dihedral *diheang;
    int numtor;
    paircont *cp;
    int numcon;
    edgetable edget; 
    int numbone;
    point3s *pbone;
    int numtraj,lengbb;
    double *etraj;
    bspline *bsp;
    double *btscore;
    bool *btmat;
    int numcont[2];
    segment *segcont[2];
    int resnumber;
    int *resmap;
    closeca *cca;
    double *afd;
    double *distmat;
};

#endif // !defined(AFX_PARSEPDB_H__0A306B68_CEEB_4EEF_8D55_4CCE700E11A3__INCLUDED_)
