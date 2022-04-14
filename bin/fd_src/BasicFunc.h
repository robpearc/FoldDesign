///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef BASICFUNC_H
#define BASICFUNC_H

#include "CommonPara.h"

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class BasicFunc  
{
  public:
    BasicFunc(); 
    virtual ~BasicFunc();

    double sdet(double a[],int n);
    void copymat(double a[],int m,int n,double b[]);
    void tranmat(double a[],int m,int n,double b[]);
    void a2rot(double ang,double rot[]);
    double phi(double xi,double yi,double zi,double xj,double yj,double zj,double xk,
               double yk,double zk,double xl,double yl,double zl);
    float phi(float xi,float yi,float zi,float xj,float yj,float zj,float xk,
              float yk,float zk,float xl,float yl,float zl);
    bool tor2pos22(float xi,float yi,float zi,float xj,float yj,float zj,float xk,
                   float yk,float zk,float tang,float tleng,float tinner,
                   float *xl,float *yl,float *zl);//quaternion
    double maxinthree(float fx,float fy,float fz);
    double footpoint(point3d ps1,point3d pe1, point3d tp,point3d *tfp,double *tdist);
    int linecross(point3d ps1,point3d pe1,point3d ps2,point3d pe2,
                  point3d *pc1,point3d *pc2,double  *dist);
    point3d prod(point3d p1,point3d p2);
    point3d minu(point3d p1,point3d p2);
    point3d addv(point3d p1,point3d p2);
    point3d scal(point3d p1,double f);
    point3d setv(double tx,double ty,double tz);
    double dotv(point3d p1,point3d p2);
    point3d ranv(double fac);
    point3d ranv_ergodic(double fac);
    point3d rana(double theta);
    double angv(point3d p1,point3d p2);
    double calcCosineAngle(point3d p1,point3d p2,point3d p3,point3d p4);
    double calcKappa(point3d p1,point3d p2, point3d p3, point3d p4);
    int maxnormal(point3d p1);
    double  norm(point3d p1);
    point3d unit(point3d p1);
    point3d rotv(point3d p1,double rot[]);
    point3d mmat(double mat[3][3],point3d p1);
    point3d mmat(double mat[9],point3d p1);
    bool rinv(double a[], int n);
    bool rinv(double **a, int n);
    void v2rot(point3d p1,double rot[]);
    double squgaussian(double x,double sigma,double miu);
    double fungaussian(double x,double fac,double sigma,double miu);
    int posinarray(double *p,int n,double cval);
    int findpos(double *p,int is,int ie,double pt);//from small to large
    int findpos2(double *p,int is,int ie,double pt);//from small to large
    int aminoid(char *aminoname);
    int aminoid(char aminoname);
    int atomid(char *atomname);
    void  q2rot(double *q,double rmax[]);
    void cross(double a[3], double b[3], double c[3]);
    void partbub3(double *p,int np,int *ind,int nind);
    void partbub2(double *p,int np,int *ind,int nind);
    void rbub(double *p,int n,int *ind);
    void rqck(double *p,int n,int *ind);
    void rsplit(char **p,int n,int *m,int *ind);
    void rsplit(double *p,int n,int *m, int *ind);
    void trmul(double a[],double b[],int m,int n,int k,double c[]);
    int muav(double a[],int m,int n,double u[],double v[],double eps,int ka);
    void sss(double fg[2],double cs[2]);
    void ppp(double a[],double e[],double s[],double v[],int m,int n);

 
  private:
    bool flagVerbose;

};

#endif 
