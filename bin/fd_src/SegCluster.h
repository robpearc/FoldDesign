///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#if !defined(AFX_SEGCLUSTER_H__5D48C60B_C78A_456A_95B9_B9ED36745C7E__INCLUDED_)
#define AFX_SEGCLUSTER_H__5D48C60B_C78A_456A_95B9_B9ED36745C7E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class SegCluster  
{
  public:
    SegCluster();
    int nosample;
    int ck;
    featvec *fv;
    cfeatvec *cfv;
    int numiter;
    void creatfv(int num);
    void creatcfv(int knum);
    void normalizefv();
    void runkmeans();
    double eucdist(int p, int c);
    int closestcluster(int pat);
    void distrisamples();
    bool  calcnewcenters();
    void recoverfeat();

    void outputdist(int num,int indi,char *proname,double sigma,char *outname,double *gauval);
    int calcnumcluster(int num,int indi,char *distname,double distcut,bool *flagexclude,bool *flagcluster);
    int calcmaxneighbor(int num,char *distname,double distcut,bool *flagexclude,int *totneigh,
                        bool *flagcluster);
    int calcmodels(int num,char *distname,double probmin,double probmax,double initdist,double deltadist,
                   int nummodel,modelinfo *minf,bool **flagmodels);
    void calccuremodels(int num,int nk,double **inmat,modelinfo *minf,bool **flagmodels);
    double calccuredist(int num,int nk,double **inmat,int ci,int cj,int *labmat);//minidist in two cluster
    double calccuredist2(int num,int nk,double **inmat,int ci,int cj,int *labmat);//center dist in two better
    double calccuremerge(int num,int nk,double **inmat,int totclu,int *ci,int *cj,int *labmat);
    double calckmeanmodels(int num,int nk,double **inmat,modelinfo *minf,bool **flagmodels);
    double calckmeanmodels2(int num,int nk,double **inmat,modelinfo *minf,bool **flagmodels);//no initial center
    void calckmeansplit(int num,int nk,double **inmat,modelinfo *minf,bool **flagmodels,int *minii,int *minij,
                        int *maxiv);
    void calckmeanfarest(int num,double **inmat,bool *flagcluster,int *ii,int *jj,
                         bool *flagclusterii,bool *flagclusterjj);
    void calckmeanmodelsiter(int num,int nk,double **inmat,modelinfo *minf,bool **flagmodels);
    double calckmeanvar(int num,int nk,double **inmat,modelinfo *minf,bool **flagmodels);//tot var
    int calckmeanclass(int num,int nk,double **inmat,modelinfo *minf,bool **flagmodels);//totchange
    int calckmeancenter(int num,double **inmat,bool *flagcluster);//center index	
    void calcdistmat(int num,dihedral *intor,double **outmat);
    void calcfragdistmat(int indi,char *fragname,double **outmat);//topse
    void calcgaussian(int num,double **inmat,double sigma,double *outvect);
    int calcgaumodels(int num,double **inmat,double sigma,int *vind,
                      double probmin,double probmax,double deltadist,int nummodel,modelinfo *minf,
                      bool **flagmodels);
    int calcnumcluster(int num,int indi,double **inmat,double distcut,bool *flagexclude,bool *flagcluster);
    int calcmaxneighbor(int num, double **inmat,double distcut,bool *flagexclude,int *totneigh,
                        bool *flagcluster);
    int calcmodels(int num, double **inmat,double probmin,double probmax,double initdist,double deltadist,
                   int nummodel,modelinfo *minf,bool **flagmodels);
    int calcmodels(int num, double **inmat,double probmin,double probmax,double initdist,double deltadist,
                   int nummodel,bool *flagexclude,modelinfo *minf,bool **flagmodels);
    int calcmodels(char *listname,char *outname);//get templates
    int calcmodels(char *listname,char *outname,double idcut);//nonhomologous
    void calcmeanstd(int num, double **inmat,double *tmean,double *tstd,double *tmin,double *tmax);
    void calcmeanstd(int num, double **inmat,bool *flagexclude,double *tmean,double *tstd,double *tmin,
                     double *tmax);
    void calcmeanstd(int num,char *distname,double *tmean,double *tstd,double *tmin,double *tmax);
    double rmin,rmax;
    double *dscore;
    int extractremc(int istart,int iend,char *prename,int nummc,char *outname);
    int extractremc(int istart,int iend,int numseq,char *prename,int nummc,char *outname);//same as above
    int extractremc(char *filename,int indstart,char *outname);
    bool extractone(char *filename,int indnum,char *outname);
    void calcinfca(point3d *pbb,int numbone,int numrec,bool *posinf,int *num1,int *num2);
    virtual ~SegCluster();
};

#endif // !defined(AFX_SEGCLUSTER_H__5D48C60B_C78A_456A_95B9_B9ED36745C7E__INCLUDED_)
