///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "BasicFunc.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

BasicFunc::BasicFunc()
{
  flagVerbose=true;
}

BasicFunc::~BasicFunc()
{
}

double BasicFunc::squgaussian(
  double x,
  double sigma,
  double miu
)
{
  return -(x-miu)*(x-miu)/(2.0*sigma*sigma);
}

double BasicFunc::fungaussian(
  double x,
  double fac,
  double sigma,
  double miu
)
{
  return fac*exp(-(x-miu)*(x-miu)/2.0/sigma/sigma);
}

void BasicFunc::copymat(
  double a[],
  int m,
  int n,
  double b[]
)
{
  int i,j;
  for(i=0;i<m;i++)
  {
    for(j=0;j<n;j++)
    {
      b[i*n+j]=a[i*n+j];
    }
  }
}

void BasicFunc::tranmat(
  double a[],
  int m,
  int n,
  double b[]
)
{
  int i,j;
  for(i=0;i<m;i++)
  {
    for(j=0;j<n;j++)
    {
      b[j*m+i]=a[i*n+j];
    }
  }
}

double BasicFunc::sdet(
  double a[], //content in a[] will be changed
  int n
)
{
  int i,j,k,is,js,l,u,v;
  double f,det,q,d;
  f=1.0; det=1.0;
  for (k=0; k<=n-2; k++)
  { 
    q=0.0;
    for (i=k; i<=n-1; i++)
    {
      for (j=k; j<=n-1; j++)
      {
        l=i*n+j; 
        d=fabs(a[l]);
        if (d>q) 
        {
          q=d; is=i; js=j;
        }
      }
    }
    if (q+1.0==1.0)
    {
      det=0.0; 
      return det;
    }
    if (is!=k)
    {
      f=-f;
      for (j=k; j<=n-1; j++)
      {
        u=k*n+j;
        v=is*n+j;
        d=a[u]; 
        a[u]=a[v];
        a[v]=d;
      }
    }
    if (js!=k)
    {
      f=-f;
      for (i=k; i<=n-1; i++)
      {
        u=i*n+js; v=i*n+k;
        d=a[u]; a[u]=a[v]; a[v]=d;
      }
    }
    l=k*n+k;
    det=det*a[l];
    for (i=k+1; i<=n-1; i++)
    {
      d=a[i*n+k]/a[l];
      for (j=k+1; j<=n-1; j++)
      {
        u=i*n+j;
        a[u]=a[u]-d*a[k*n+j];
      }
    }
  }
  det=f*det*a[n*n-1];
  
  return det;
}

double BasicFunc::dotv(
  point3d p1,
  point3d p2
)
{
  return p1.x*p2.x+p1.y*p2.y+p1.z*p2.z;
}

void  BasicFunc::q2rot(
  double *q,
  double rmax[]
)
{
  int i,j;
  for(i=0;i<3;i++)
  {
    for(j=0;j<3;j++)
    {
      rmax[i*3+j]=0;
    }
  }
  rmax[0]=q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3];
  rmax[1]=2*(q[1]*q[2]-q[0]*q[3]);
  rmax[2]=2*(q[1]*q[3]+q[0]*q[2]);
  rmax[3]=2*(q[1]*q[2]+q[0]*q[3]);
  rmax[4]=q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3];
  rmax[5]=2*(q[2]*q[3]-q[0]*q[1]);
  rmax[6]=2*(q[1]*q[3]-q[0]*q[2]);
  rmax[7]=2*(q[2]*q[3]+q[0]*q[1]);
  rmax[8]=q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3];
}

bool BasicFunc::rinv(
  double a[], 
  int n
)
{
  int *is,*js,i,j,k,l,u,v;
  double d,p;
  is=(int *)malloc(n*sizeof(int));
  js=(int *)malloc(n*sizeof(int));
  for(k=0; k<=n-1; k++)
  {
    d=0.0;
    for(i=k; i<=n-1; i++)
    {
      for(j=k; j<=n-1; j++)
      { 
        l=i*n+j;
        p=fabs(a[l]);
        if(p>d) 
        { 
          d=p; is[k]=i; js[k]=j;
        }
      }
    }
    if(d+1.0==1.0)
    { 
      free(is);
      free(js);
      printf("err**not inv\n");
      return(false);
    }
    if(is[k]!=k)
    {
      for(j=0; j<=n-1; j++)
      { 
        u=k*n+j; v=is[k]*n+j;
        p=a[u]; a[u]=a[v]; a[v]=p;
      }
    }
    if(js[k]!=k)
    {
      for(i=0; i<=n-1; i++)
      { 
        u=i*n+k; v=i*n+js[k];
        p=a[u]; a[u]=a[v]; a[v]=p;
      }
    }
    l=k*n+k;
    a[l]=1.0/a[l];
    for(j=0; j<=n-1; j++)
    {
      if(j!=k)
      { 
        u=k*n+j; a[u]=a[u]*a[l];
      }
    }
    for(i=0; i<=n-1; i++)
    {
      if(i!=k)
      {
        for(j=0; j<=n-1; j++)
        {
          if(j!=k)
          {
            u=i*n+j;
            a[u]=a[u]-a[i*n+k]*a[k*n+j];
          }
        }
      }
    }
    for(i=0; i<=n-1; i++)
    {
      if(i!=k)
      { 
        u=i*n+k; a[u]=-a[u]*a[l];
      }
    }
  }
  for(k=n-1; k>=0; k--)
  { 
    if(js[k]!=k)
    {
      for(j=0; j<=n-1; j++)
      {
        u=k*n+j; v=js[k]*n+j;
        p=a[u]; a[u]=a[v]; a[v]=p;
      }
    }
    if(is[k]!=k)
    {
      for(i=0; i<=n-1; i++)
      {
        u=i*n+k; v=i*n+is[k];
        p=a[u]; a[u]=a[v]; a[v]=p;
      }
    }
  }
  free(is); 
  free(js);

  return true;
}

int BasicFunc::linecross(
  point3d ps1,
  point3d pe1,
  point3d ps2,
  point3d pe2,
  point3d *pc1,
  point3d *pc2,
  double *dist
)
{
  point3d pse1,pse2,pss12,tfp;
  double k1,k2;
  double ang1,tdist;
  double dp1,dp2,dp3,dp4,minidp;
  int inddp;
  pse1=minu(pe1,ps1);
  pse2=minu(pe2,ps2);
  pss12=minu(ps1,ps2);
  ang1=angv(pse1,pse2);
  if(ang1<epsilon || 180-ang1<epsilon) //parallel
  {
    k1=footpoint(ps1,pe1,ps2,&tfp,&tdist);
    if(tdist>=0 && tdist<=norm(pse1))
    {
      *pc1=ps2;
      *pc2=tfp;
      return 9; //good
    }
    k1=footpoint(ps1,pe1,pe2,&tfp,&tdist);
    if(tdist>=0 && tdist<=norm(pse1))
    {
      *pc1=pe2;
      *pc2=tfp;
      return 9; //good
    }
    k1=footpoint(ps2,pe2,ps1,&tfp,&tdist);
    if(tdist>=0 && tdist<=norm(pse2))
    {
      *pc1=ps1;
      *pc2=tfp;
      return 9; //good
    }
    k1=footpoint(ps2,pe2,pe1,&tfp,&tdist);
    if(tdist>=0 && tdist<=norm(pse2))
    {
      *pc1=pe1;
      *pc2=tfp;
      return 9; //good
    }
    dp1=norm(minu(ps1,ps2));
    inddp=1;
    minidp=dp1;
    *pc1=ps1;*pc2=ps2;
    dp2=norm(minu(ps1,pe2));
    dp3=norm(minu(pe1,ps2));
    dp4=norm(minu(pe1,pe2));
    if(dp2<minidp)
    {
      *pc1=ps1;*pc2=pe2;
      inddp=2;
      minidp=dp2;
    }
    if(dp3<minidp)
    {
      *pc1=pe1;*pc2=ps2;
      inddp=3;
      minidp=dp3;
    }
    if(dp4<minidp)
    {
      *pc1=pe1;*pc2=pe2;
      inddp=4;
      minidp=dp4;
    }
    *dist=minidp;
    return 9+inddp;
  }
  ////coplane
  double tmat[9];
  tmat[0]=pe1.x-ps1.x;
  tmat[1]=pe1.y-ps1.y;
  tmat[2]=pe1.z-ps1.z;
  tmat[3]=ps2.x-ps1.x;
  tmat[4]=ps2.y-ps1.y;
  tmat[5]=ps2.z-ps1.z;
  tmat[6]=pe2.x-ps1.x;
  tmat[7]=pe2.y-ps1.y;
  tmat[8]=pe2.z-ps1.z;
  double tdet=sdet(tmat,3);
  point3d pdi=prod(pse1,pse2);
  int indmax=maxnormal(pdi);
  if(fabs(tdet)<epsilon) //coplane
  {
    if(indmax==2)
      k1=(pss12.x*pse2.y-pss12.y*pse2.x)/(pse1.y*pse2.x-pse1.x*pse2.y);
    else if(indmax==1)
      k1=(pss12.x*pse2.z-pss12.z*pse2.x)/(pse1.z*pse2.x-pse1.x*pse2.z);
    else
      k1=(pss12.z*pse2.y-pss12.y*pse2.z)/(pse1.y*pse2.z-pse1.z*pse2.y);
    *pc1=addv(ps1,scal(pse1,k1));
    *pc2=*pc1;
    *dist=0;
    k2=norm(minu(*pc1,ps2))/norm(pse2);
    if(k1>=0 && k1<=1 && k2>=0 && k2<=1)
      return 14; //good
    else if(k1<0 && k2>=0 && k2<=1)
      return 15;
    else if(k1>1 && k2>=0 && k2<=1)
      return 16;
    else if(k1>=0 && k1<=1 && k2<0)
      return 17;
    else if(k1<0 && k2<0)
      return 18;
    else if(k1>1 && k2<0)
      return 19;
    else if(k1>=0 && k1<=1 && k2>1)
      return 20;
    else if(k1<0 && k2>1)
      return 21;
    else if(k1>1 && k2>1)
      return 22;
  }
  ////
  tfp=minu(ps2,ps1);
  *dist=dotv(tfp,pdi)/norm(pdi);
  point3d pfc1=addv(ps1,scal(pdi,*dist/norm(pdi)));
  point3d pfc2=addv(pe1,scal(pdi,*dist/norm(pdi)));
  pse1=minu(pfc2,pfc1);
  pss12=minu(pfc1,ps2);
  if(indmax==2)
    k1=(pss12.x*pse2.y-pss12.y*pse2.x)/(pse1.y*pse2.x-pse1.x*pse2.y);
  else if(indmax==1)
    k1=(pss12.x*pse2.z-pss12.z*pse2.x)/(pse1.z*pse2.x-pse1.x*pse2.z);
  else
    k1=(pss12.z*pse2.y-pss12.y*pse2.z)/(pse1.y*pse2.z-pse1.z*pse2.y);
  *pc2=addv(pfc1,scal(pse1,k1));
  *pc1=minu(*pc2,scal(pdi,*dist/norm(pdi)));
  *dist=fabs(*dist);
  k2=norm(minu(*pc1,ps2))/norm(pse2);
  if(k1>=0 && k1<=1 && k2>=0 && k2<=1)
    return 0;//good
  else if(k1<0 && k2>=0 && k2<=1)
    return 1;
  else if(k1>1 && k2>=0 && k2<=1)
    return 2;
  else if(k1>=0 && k1<=1 && k2<0)
    return 3;
  else if(k1<0 && k2<0)
    return 4;
  else if(k1>1 && k2<0)
    return 5;
  else if(k1>=0 && k1<=1 && k2>1)
    return 6;
  else if(k1<0 && k2>1)
    return 7;
  else if(k1>1 && k2>1)
    return 8;

  return 100;
}

//vertical dist tfp foot tdist-to ps1
double BasicFunc::footpoint(
  point3d ps1,
  point3d pe1, 
  point3d tp,
  point3d *tfp,
  double *tdist
)
{
  point3d pse,pstp;
  double ang1,k;
  pse=minu(pe1,ps1);
  pstp=minu(tp,ps1);
  ang1=angv(pstp,pse);
  if(ang1<epsilon || 180-ang1<epsilon)
  {
    *tfp=tp;
    if(ang1<epsilon) *tdist=1;
    else *tdist=-1;
    *tdist*=norm(pstp);
    return 0.0;
  }
  k=dotv(pse,pstp)/dotv(pse,pse);
  *tfp=addv(ps1,scal(pse,k));
  *tdist=k*norm(pse);
  pstp=minu(tp,*tfp);
    
  return norm(pstp);
}

double BasicFunc::maxinthree(
  float fx,
  float fy,
  float fz
)
{
  double tmax=fx;
  if(tmax<fy) tmax=fy;
  if(tmax<fz) tmax=fz;
    
  return tmax;
}

int BasicFunc::maxnormal(
  point3d p1
)
{
  int i;
  double tem,temp[3];
  temp[0]=fabs(p1.x);
  temp[1]=fabs(p1.y);
  temp[2]=fabs(p1.z);
  i=0;
  tem=temp[0];
  if(temp[1]>tem)
  {
    i=1;
    tem=temp[1];
  }
  if(temp[2]>tem)
  {
    i=2;
    tem=temp[2];
  }
    
  return i;
}

point3d BasicFunc::unit(
  point3d p1
)
{
  point3d temp;
  double t=sqrt(p1.x*p1.x+p1.y*p1.y+p1.z*p1.z);
  if(t<1e-40)
  {
    temp.x=0;
    temp.y=0;
    temp.z=0;
    return temp;
  }
  temp.x=p1.x/t;
  temp.y=p1.y/t;
  temp.z=p1.z/t;
  
  return temp;
}

void BasicFunc::a2rot(
  double ang,
  double rot[]
)
{
  //the rotmatrix for around [0,0,1] angle
  double cosphi=cos(ang);
  double sinphi=sin(ang);
  rot[0]=cosphi;rot[1]=-sinphi;rot[2]=0;
  rot[3]=sinphi;rot[4]=cosphi; rot[5]=0;
  rot[6]=0;     rot[7]=0;      rot[8]=1;
}

void BasicFunc::v2rot(
  point3d p1,
  double rot[]
)
{
  //the rotmatrix for [0,0,1] to p1
  double pr=sqrt(p1.x*p1.x+p1.y*p1.y+p1.z*p1.z);
  double theta=acos(p1.z/pr);
  double phi=atan2(p1.y,p1.x)+PI/2.0;
  double cosphi=cos(phi);
  double sinphi=sin(phi);
  double costheta=cos(theta);
  double sintheta=sin(theta);
  rot[0]=cosphi;rot[1]=-sinphi*costheta;rot[2]=sinphi*sintheta;
  rot[3]=sinphi;rot[4]=cosphi*costheta; rot[5]=-cosphi*sintheta;
  rot[6]=0;     rot[7]=sintheta;        rot[8]=costheta;
}

point3d BasicFunc::rotv(
  point3d p1,
  double rot[]
)
{
  point3d tv;
  tv.x=rot[0]*p1.x+rot[1]*p1.y+rot[2]*p1.z;
  tv.y=rot[3]*p1.x+rot[4]*p1.y+rot[5]*p1.z;
  tv.z=rot[6]*p1.x+rot[7]*p1.y+rot[8]*p1.z;
    
  return tv;
}

point3d BasicFunc::mmat(
  double mat[9],
  point3d p1
)
{
  point3d temp;
  temp.x=p1.x*mat[0]+p1.y*mat[1]+p1.z*mat[2];
  temp.y=p1.x*mat[3]+p1.y*mat[4]+p1.z*mat[5];
  temp.z=p1.x*mat[6]+p1.y*mat[7]+p1.z*mat[8];
    
  return temp;
}

point3d BasicFunc::setv(
  double tx,
  double ty,
  double tz
)
{
  point3d temp;
  temp.x=tx;
  temp.y=ty;
  temp.z=tz;
    
  return temp;
}

point3d BasicFunc::prod(
  point3d p1,
  point3d p2
)//from p1 to p2
{
  point3d temp;
  temp.x=p1.y*p2.z-p1.z*p2.y;
  temp.y=p1.z*p2.x-p1.x*p2.z;
  temp.z=p1.x*p2.y-p1.y*p2.x;
  
  return temp;
}

point3d BasicFunc::minu(
  point3d p1,
  point3d p2
)
{
  point3d temp;
  temp.x=p1.x-p2.x;
  temp.y=p1.y-p2.y;
  temp.z=p1.z-p2.z;
    
  return temp;
}

point3d BasicFunc::addv(
  point3d p1,
  point3d p2
)
{
  point3d temp;
  temp.x=p1.x+p2.x;
  temp.y=p1.y+p2.y;
  temp.z=p1.z+p2.z;
    
  return temp;
}

point3d BasicFunc::scal(
  point3d p1,
  double f
)
{
  point3d temp;
  temp.x=f*p1.x;
  temp.y=f*p1.y;
  temp.z=f*p1.z;
    
  return temp;
}

double BasicFunc::norm(
  point3d p1
)
{
  return sqrt(p1.x*p1.x+p1.y*p1.y+p1.z*p1.z);
}

double BasicFunc::angv(
  point3d p1,
  point3d p2
)
{
  double t1=norm(p1);
  double t2=norm(p2);
  if(t1<epsilon || t2<epsilon) return 0;
  double t3=dotv(p1,p2)/t1/t2;
  if(t3<-1.0) t3=-1.0;
  else if(t3>1.0) t3=1.0;

  return acos(t3);
}

//p1 calpha, p2 prevprev calpha, p3 nextnext calpha, calpha
double BasicFunc::calcKappa(
  point3d p1, //ith Ca atom
  point3d p2, //i-2 Ca atom
  point3d p3, //i+2 Ca atom
  point3d p4  //ith Ca atom
)
{
  double result = 360;
  double ckap = calcCosineAngle(p1, p2, p3, p4);
  double skap = sqrt(1 - ckap * ckap);
  result = atan2(skap, ckap) * 180 / PI;
    
  return result;
}


double BasicFunc::calcCosineAngle(
  point3d p1,
  point3d p2,
  point3d p3,
  point3d p4
)
{
  point3d v12 = minu(p1,p2);
  point3d v34 = minu(p3,p4);

  double result = 0;

  double x = dotv(v12, v12) * dotv(v34, v34);
  if(x > 0)
  {
    result = dotv(v12, v34) / sqrt(x);
  }

  return result;
}

float BasicFunc::phi(
  float xi,float yi,float zi,
  float xj,float yj,float zj,
  float xk,float yk,float zk,
  float xl,float yl,float zl
)
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
    
  /* Calculate the vectors C,B,C */
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
  dxi = yij * zkj - zij * ykj;     /* Normal to plane 1 */
  dyi = zij * xkj - xij * zkj;
  dzi = xij * ykj - yij * xkj;
  gxi = zkj * ykl - ykj * zkl;     /* Mormal to plane 2 */
  gyi = xkj * zkl - zkj * xkl;
  gzi = ykj * xkl - xkj * ykl;

  /* Calculate the length of the two normals */
  bi = dxi * dxi + dyi * dyi + dzi * dzi;
  bk = gxi * gxi + gyi * gyi + gzi * gzi;
  ct = dxi * gxi + dyi * gyi + dzi * gzi;

  boi2 = 1./bi;
  boj2 = 1./bk;
  bi   = (double)sqrt((double)bi);
  bk   = (double)sqrt((double)bk);
  if(bi<epsilon*0.01 || bk<epsilon*0.01) return 180;
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

  return(ap);
}

double BasicFunc::phi(
  double xi,double yi,double zi,
  double xj,double yj,double zj,
  double xk,double yk,double zk,
  double xl,double yl,double zl
)
{
  double xij,yij,zij,
         xkj,ykj,zkj,
         xkl,ykl,zkl,
         dxi,dyi,dzi,
         gxi,gyi,gzi,
         bi,bk,ct,
         boi2,boj2,
         z1,z2,ap,s,
         bioj,bjoi;
    
  /* Calculate the vectors C,B,C */
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
  dxi = yij * zkj - zij * ykj;     /* Normal to plane 1 */
  dyi = zij * xkj - xij * zkj;
  dzi = xij * ykj - yij * xkj;
  gxi = zkj * ykl - ykj * zkl;     /* Normal to plane 2 */
  gyi = xkj * zkl - zkj * xkl;
  gzi = ykj * xkl - xkj * ykl;

  /* Calculate the length of the two normals */
  bi = dxi * dxi + dyi * dyi + dzi * dzi;
  bk = gxi * gxi + gyi * gyi + gzi * gzi;
  ct = dxi * gxi + dyi * gyi + dzi * gzi;

  boi2 = 1./bi;
  boj2 = 1./bk;
  bi   = (double)sqrt((double)bi);
  bk   = (double)sqrt((double)bk);
  if(bi<epsilon*0.01 || bk<epsilon*0.01) return 180;
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
  ap *= (double)180.0/PI;
  if(ap<0) ap+=360;

  return(ap);
}

/* Sphere surface semi-random sampling for given theta (angle between
 * z-axis and the sampled vector) and random Phi */
point3d BasicFunc::rana(
  double theta
)
{
  point3d tp;
  double dphi;
  dphi=(2.0*Random())*PI;
  tp.x=sin(theta)*cos(dphi);
  tp.y=sin(theta)*sin(dphi);
  tp.z=cos(theta);

  return(tp);
}

/* This function is supposed to sample points within a sphere centered at
 * origin and with radius being the absolute value of 'fac'. If 'fac' is
 * positive, points are sampled inside the sphere; if 'fac' is on negative,
 * points are just sampled just on surface. Note that this is NOT ergodic
 * sampling, as regions close to the two poles are much more likely to be
 * sampled compared to region close to equator. Use BasicFunc::ranv_ergodic
 * instead. */
point3d BasicFunc::ranv(
  double fac
)
{
  point3d tp;
  double dr=fabs(fac);               // sample on surface
  if (fac>0) dr=0.00001+fac*Random();// sample inside sphere 
  double dtheta=(1.0*Random())*PI;   // [0, PI)
  double dphi=(2.0*Random()-1.0)*PI; // [-PI,PI)
  tp.x=dr*sin(dtheta)*cos(dphi);
  tp.y=dr*sin(dtheta)*sin(dphi);
  tp.z=dr*cos(dtheta);

  return(tp);
}

/* This function has the same purpose as BasicFunc::ranv, except that the
 * non-ergodicity issue in BasicFunc::ranv is fixed. Absolute value of
 * 'fac' is sampling radius. If 'fac' is negative, points are sampled on
 * surface of sphere; if 'fac' is positive, points are sampled inside the
 * sphere */
point3d BasicFunc::ranv_ergodic(
  double fac
)
{
  point3d tp;
  double dr=fabs(fac);               // sample on surface
  if (fac>0) dr=0.00001+fac*Random();// sample inside sphere 
  double v=(2.*Random()-1.);         // [-1, 1)
  double dtheta=acos(v);             // [0, PI)
  double dphi=(2.*Random()-1.)*PI;   // [-PI,PI)
  tp.x=dr*sin(dtheta)*cos(dphi);
  tp.y=dr*sin(dtheta)*sin(dphi);
  tp.z=dr*cos(dtheta);

  return(tp);
}

bool BasicFunc::tor2pos22(      
  float xi,float yi,float zi,
  float xj,float yj,float zj, 
  float xk,float yk,float zk,
  float tang,float tleng,float tinner, 
  float *xl,float *yl,float *zl
)
{ 
  point3d p12,p23,e1,e2,e3;
  point3d p1,p2,p3;
  double tpangle,tcos,tsin,rmat[9];
  double q[4];
  p12.x=xi-xj;p12.y=yi-yj;p12.z=zi-zj;
  p23.x=xk-xj;p23.y=yk-yj;p23.z=zk-zj;

  if(norm(p12)<epsilon*0.00001 || norm(p23)<epsilon*0.00001 || 
     angv(p12,p23)<epsilon*0.001 || (PI-angv(p12,p23))<epsilon*0.001)
  {
    int imax=maxnormal(p23);
    if(imax==0) yj+=0.01f;
    else if(imax==1) zj+=0.01f;
    else xj+=0.01f;
    p12.x=xi-xj;p12.y=yi-yj;p12.z=zi-zj;
    p23.x=xk-xj;p23.y=yk-yj;p23.z=zk-zj;
    if(flagVerbose)
    {
      printf("make adjustment tor2pos22\n");
    }
  }
  e2=unit(p23);
  e3=unit(p12);
  e1=prod(e2,e3);
  e1=unit(e1);
  if(norm(e1)<epsilon || norm(e2)<epsilon)
  {
    if(flagVerbose)
    {
      printf("wrong in tor2pos22 [%f %f %f] [%f %f %f] [%f %f %f]\n",
             xi,yi,zi,xj,yj,zj,xk,yk,zk);
    }
    *xl=xk;*yl=yk;*zl=zk;
    return false;
  }
  p1=scal(e2,tleng);
  tpangle=(PI-tinner)/2.0;
  tcos=cos(tpangle);
  tsin=sin(tpangle);
  q[0]=tcos;q[1]=tsin*e1.x;q[2]=tsin*e1.y;q[3]=tsin*e1.z;
  q2rot(q,rmat);
  p2=mmat(rmat,p1);
    
  tpangle=tang/2.0;
  tcos=cos(tpangle);
  tsin=sin(tpangle);
  q[0]=tcos;q[1]=tsin*e2.x;q[2]=tsin*e2.y;q[3]=tsin*e2.z;
  q2rot(q,rmat);
  p3=mmat(rmat,p2);
    
  *xl=p3.x+xk;*yl=p3.y+yk;*zl=p3.z+zk;

  return true;
}

int BasicFunc::findpos(
  double *p,
  int is,
  int ie,
  double pt
)//[is,ie]
{       
  if(is==ie)
  {       
    if(pt<p[ie])
    {       
      printf("some wrong order\n");
    }
    return is;
  }
  else if(is==ie-1)
  {       
    if(pt>=p[is])
    {       
      return is;
    }
    else
    {       
      return ie;
    }
  }
  else//is!=im!=ie
  {       
    int im=(is+ie)/2;
    if(pt==p[im])
    {       
      return(im);
    }
    else if(pt>p[im])
    {       
      return(findpos(p,is,im,pt));
    }
    else
    {       
      return(findpos(p,im+1,ie,pt));
    }
  }
}


int BasicFunc::findpos2(
  double *p,
  int is,
  int ie,
  double pt
)
{
  if(is==ie) return is;
  else if(is==ie-1)
  {
    if(pt<=p[is]) return is;
    else return ie;
  }
  else//is!=im!=ie
  {
    int im=(is+ie)/2;
    if(pt==p[im]) return(im);
    else if(pt<p[im]) return(findpos2(p,is,im,pt));
    else return(findpos2(p,im+1,ie,pt));
  }
}

int BasicFunc::posinarray(
  double *p,
  int n,
  double cval
)
{
  //cval in [p[i],p[i+1]) output i
  if(n==0) return 0;
  else if(n==1)
  {
    if(cval<p[0]) return -1;
    else return 0;
  }
  else if(n==2)
  {
    if(cval<p[0]) return -1;
    else if(cval>p[1]) return 1;
    else return 0;
  }
  int m=n/2;
  if(cval<p[m]) return posinarray(p,m+1,cval);
  else return m+posinarray(p+m,n-m,cval);
}

int BasicFunc::atomid(
  char *atomname
)
{
  int i;
  for(i=0;i<83;i++)
  {
    if(strcmp(atomname,sTab4[i])==0) return i;
  }
  
  return 83;
}

int BasicFunc::aminoid(
  char aminoname
)
{
  int i;
  //empty
  if(aminoname==' ') return 23;
  // one in 26
  for(i=0;i<26;i++)
  {
    if(aminoname==aad1[i]) return i;
  }
  
  return 23;
}

void BasicFunc::partbub3(
  double *p, 
  int np, 
  int *ind, 
  int nind
)
{
  int threshnind=2*nind;
  if(threshnind<200) threshnind=200;
  if(np<=nind)//all
  {
    rbub(p,np,ind);
    return;
  }
  else if(np<threshnind)//np is small
  {
    partbub2(p,np,ind,nind);
    return;
  }
  int *i,i0;
  i=&i0;
  rsplit(p,np,i,ind);
  if(i0+1>=nind)
  {
    partbub3(p,i0,ind,nind);
  }
  else
  {
    rqck(p,i0,ind);
    partbub3(p+i0+1, np-i0-1, ind+i0+1, nind-i0-1);
  }
}

void BasicFunc::partbub2(
  double *p, 
  int np, 
  int *ind, 
  int nind
)
{
  int i,j;
  int temind;
  double temp;
  for(i=0;i<nind;i++)
  {
    for(j=i+1;j<nind;j++)
    {
      if(p[j]>p[i])
      {
        temp=p[i];p[i]=p[j];p[j]=temp;
        temind=ind[i];ind[i]=ind[j];ind[j]=temind;
      }
    }
  }

  int ipos;
  for(i=nind;i<np;i++)
  {
    if(p[i]>p[nind-1])
    {
      ipos=findpos(p,0,nind-1,p[i]);
      for(j=nind-2;j>=ipos;j--)
      {
        p[j+1]=p[j];
        ind[j+1]=ind[j];
      }
      p[ipos]=p[i];
      ind[ipos]=ind[i];
    }//if
  }
}

void BasicFunc::rbub(
  double *p,
  int n, 
  int *ind
)
{
  int m,k,j,i;
  double d;
  k=0;
  m=n-1;
  while (k<m)
  {
    j=m-1; m=0;
    for(i=k; i<=j; i++)
    {
      if(p[i]>p[i+1])
      {
        d=p[i]; p[i]=p[i+1]; p[i+1]=d; m=i;
        int ii=ind[i];ind[i]=ind[i+1];ind[i+1]=ii;//
      }
    }
    j=k+1; k=0;
    for(i=m;i>=j;i--)
    {
      if(p[i-1]>p[i])
      {
        d=p[i]; p[i]=p[i-1]; p[i-1]=d; k=i;
        int ii=ind[i];ind[i]=ind[i-1];ind[i-1]=ii;//
      }
    }
  }
  return;
}

void BasicFunc::rqck(
  double *p,
  int n, 
  int *ind
)
{
  int m,i0,*i;
  double *s;
  i=&i0;
  if (n>10)
  {
    rsplit(p,n,i,ind);
    m=i0;
    rqck(p,m,ind);
    s=p+(i0+1);
    m=n-(i0+1);
    rqck(s,m,ind+(i0+1));//
  }
  else rbub(p,n,ind);
    
  return;
}

void BasicFunc::rsplit(
  char **p,
  int n,
  int *m, 
  int *ind
)
{
  int i,j,k,l;
  char t[10];
  i=0;
  j=n-1;
  k=(i+j)/2;
  if(strcmp(p[i],p[j])>=0 && strcmp(p[j],p[k])>=0) l=j;
  else if(strcmp(p[i],p[k])>=0 && strcmp(p[k],p[j])>=0) l=k;
  else l=i;
  strcpy(t,p[l]);
  strcpy(p[l],p[i]);
  int ii=ind[l];//
  ind[l]=ind[i];//

  while(i!=j)
  {
    while((i<j)&& strcmp(p[j],t)>=0) j=j-1;
    if(i<j)
    {
      strcpy(p[i],p[j]);
      ind[i]=ind[j];//
      i=i+1;
      while((i<j)&& strcmp(p[i],t)<=0) i=i+1;
      if(i<j)
      {
        strcpy(p[j],p[i]);
        ind[j]=ind[i];//
        j=j-1;
      }
    }
  }
  strcpy(p[i],t);
  ind[i]=ii;//
  *m=i;
    
  return;
}

void BasicFunc::rsplit(
  double *p,
  int n,
  int *m, 
  int *ind
)
{
  int i,j,k,l;
  double t;
  i=0;
  j=n-1;
  k=(i+j)/2;
  if((p[i]>=p[j])&&(p[j]>=p[k])) l=j;
  else if((p[i]>=p[k])&&(p[k]>=p[j])) l=k;
  else l=i;
  t=p[l];
  p[l]=p[i];
  int ii=ind[l];//
  ind[l]=ind[i];//

  while(i!=j)
  {
    while((i<j)&&(p[j]>=t)) j=j-1;
    if(i<j)
    {
      p[i]=p[j];
      ind[i]=ind[j];//
      i=i+1;
      while((i<j)&&(p[i]<=t)) i=i+1;
      if(i<j)
      {
        p[j]=p[i];
        ind[j]=ind[i];//
        j=j-1;
      }
    }
  }
  p[i]=t;
  ind[i]=ii;//
  *m=i;

  return;
}

void BasicFunc::trmul(
  double a[],
  double b[],
  int m,
  int n,
  int k,
  double c[]
)
{
  int i,j,l,u;
  for(i=0;i<=m-1;i++)
  {
    for(j=0;j<=k-1;j++)
    {
      u=i*k+j; c[u]=0.0;
      for(l=0;l<=n-1;l++)
      {
        c[u]=c[u]+a[i*n+l]*b[l*k+j];
      }
    }
  }

  return;
}

int BasicFunc::muav(
  double a[],
  int m,
  int n,
  double u[],
  double v[],
  double eps,
  int ka
)
{
  int i,j,k,l,it,ll,kk,ix,iy,mm,nn,iz,m1,ks;
  double d,dd,t,sm,sm1,em1,sk,ek,b,c,shh,fg[2],cs[2];
  double *s,*e,*w;
  s=(double *)malloc(ka*sizeof(double));
  e=(double *)malloc(ka*sizeof(double));
  w=(double *)malloc(ka*sizeof(double));
  it=60; k=n;
  if (m-1<n) k=m-1;
  l=m;
  if (n-2<m) l=n-2;
  if (l<0) l=0;
  ll=k;
  if (l>k) ll=l;
  if (ll>=1)
  {
    for (kk=1; kk<=ll; kk++)
    {
      if (kk<=k)
      {                   
        d=0.0;
        for (i=kk; i<=m; i++)
        {                       
          ix=(i-1)*n+kk-1; d=d+a[ix]*a[ix];
        }
        s[kk-1]=sqrt(d);
        if (s[kk-1]!=0.0)
        {       
          ix=(kk-1)*n+kk-1;
          if (a[ix]!=0.0)
          {                           
            s[kk-1]=fabs(s[kk-1]);
            if (a[ix]<0.0) s[kk-1]=-s[kk-1];
          }
          for (i=kk; i<=m; i++)
          {                           
            iy=(i-1)*n+kk-1;
            a[iy]=a[iy]/s[kk-1];
          }
          a[ix]=1.0+a[ix];
        }
        s[kk-1]=-s[kk-1];
      }
      if (n>=kk+1)
      {
        for (j=kk+1; j<=n; j++)
        {
          if ((kk<=k)&&(s[kk-1]!=0.0))
          {
            d=0.0;
            for (i=kk; i<=m; i++)
            {
              ix=(i-1)*n+kk-1;
              iy=(i-1)*n+j-1;
              d=d+a[ix]*a[iy];
            }
            d=-d/a[(kk-1)*n+kk-1];
            for (i=kk; i<=m; i++)
            {
              ix=(i-1)*n+j-1;
              iy=(i-1)*n+kk-1;
              a[ix]=a[ix]+d*a[iy];
            }
          }
          e[j-1]=a[(kk-1)*n+j-1];
        }
      }
      if (kk<=k)
      {
        for (i=kk; i<=m; i++)
        {
          ix=(i-1)*m+kk-1; iy=(i-1)*n+kk-1;
          u[ix]=a[iy];
        }
      }
      if (kk<=l)
      {
        d=0.0;
        for (i=kk+1; i<=n; i++) d=d+e[i-1]*e[i-1];
        e[kk-1]=sqrt(d);
        if (e[kk-1]!=0.0)
        {
          if (e[kk]!=0.0)
          {
            e[kk-1]=fabs(e[kk-1]);
            if (e[kk]<0.0) e[kk-1]=-e[kk-1];
          }
          for (i=kk+1; i<=n; i++)
             e[i-1]=e[i-1]/e[kk-1];
          e[kk]=1.0+e[kk];
        }
        e[kk-1]=-e[kk-1];
        if ((kk+1<=m)&&(e[kk-1]!=0.0))
        {
          for (i=kk+1; i<=m; i++) w[i-1]=0.0;
          for (j=kk+1; j<=n; j++)
            for (i=kk+1; i<=m; i++)
              w[i-1]=w[i-1]+e[j-1]*a[(i-1)*n+j-1];
          for (j=kk+1; j<=n; j++)
            for (i=kk+1; i<=m; i++)
            {
              ix=(i-1)*n+j-1;
              a[ix]=a[ix]-w[i-1]*e[j-1]/e[kk];
            }
        }
        for (i=kk+1; i<=n; i++)
          v[(i-1)*n+kk-1]=e[i-1];
      }
    }
  }
  mm=n;
  if (m+1<n) mm=m+1;
  if (k<n) s[k]=a[k*n+k];
  if (m<mm) s[mm-1]=0.0;
  if (l+1<mm) e[l]=a[l*n+mm-1];
  e[mm-1]=0.0;
  nn=m;
  if (m>n) nn=n;
  if (nn>=k+1)
  { 
    for (j=k+1; j<=nn; j++)
    { 
      for (i=1; i<=m; i++)
        u[(i-1)*m+j-1]=0.0;
      u[(j-1)*m+j-1]=1.0;
    }
  }
  if (k>=1)
  { 
    for (ll=1; ll<=k; ll++)
    { 
      kk=k-ll+1; iz=(kk-1)*m+kk-1;
      if (s[kk-1]!=0.0)
      { 
        if (nn>=kk+1)
          for (j=kk+1; j<=nn; j++)
          { 
            d=0.0;
            for (i=kk; i<=m; i++)
            { 
              ix=(i-1)*m+kk-1;
              iy=(i-1)*m+j-1;
              d=d+u[ix]*u[iy]/u[iz];
            }
            d=-d;
            for (i=kk; i<=m; i++)
            { 
              ix=(i-1)*m+j-1;
              iy=(i-1)*m+kk-1;
              u[ix]=u[ix]+d*u[iy];
            }
          }
        for (i=kk; i<=m; i++)
        { 
          ix=(i-1)*m+kk-1; u[ix]=-u[ix];
        }
        u[iz]=1.0+u[iz];
        if (kk-1>=1)
        for (i=1; i<=kk-1; i++)
          u[(i-1)*m+kk-1]=0.0;
      }
      else
      { 
        for (i=1; i<=m; i++)
          u[(i-1)*m+kk-1]=0.0;
        u[(kk-1)*m+kk-1]=1.0;
      }
    }
  }
  for (ll=1; ll<=n; ll++)
  { 
    kk=n-ll+1; iz=kk*n+kk-1;
    if ((kk<=l)&&(e[kk-1]!=0.0))
    { 
      for (j=kk+1; j<=n; j++)
      { 
        d=0.0;
        for (i=kk+1; i<=n; i++)
        { 
          ix=(i-1)*n+kk-1; iy=(i-1)*n+j-1;
          d=d+v[ix]*v[iy]/v[iz];
        }
        d=-d;
        for (i=kk+1; i<=n; i++)
        { 
          ix=(i-1)*n+j-1; iy=(i-1)*n+kk-1;
          v[ix]=v[ix]+d*v[iy];
        }
      }
    }
    for (i=1; i<=n; i++)
      v[(i-1)*n+kk-1]=0.0;
    v[iz-n]=1.0;
  }
  for (i=1; i<=m; i++)
    for (j=1; j<=n; j++)
      a[(i-1)*n+j-1]=0.0;
  m1=mm; it=60;
  while (1==1)
  {
    if (mm==0)
    {
      ppp(a,e,s,v,m,n);
      free(s); free(e); free(w); return(1);
    }
    if (it==0)
    {
      ppp(a,e,s,v,m,n);
      free(s); free(e); free(w); return(-1);
    }
    kk=mm-1;
    while ((kk!=0)&&(fabs(e[kk-1])!=0.0))
    {
      d=fabs(s[kk-1])+fabs(s[kk]);
      dd=fabs(e[kk-1]);
      if (dd>eps*d) kk=kk-1;
      else e[kk-1]=0.0;
    }
    if (kk==mm-1)
    {
      kk=kk+1;
      if (s[kk-1]<0.0)
      {
        s[kk-1]=-s[kk-1];
        for (i=1; i<=n; i++)
        { 
          ix=(i-1)*n+kk-1; 
          v[ix]=-v[ix];
        }
      }
      while ((kk!=m1)&&(s[kk-1]<s[kk]))
      {
        d=s[kk-1]; s[kk-1]=s[kk]; s[kk]=d;
        if (kk<n)
          for (i=1; i<=n; i++)
          { 
            ix=(i-1)*n+kk-1; iy=(i-1)*n+kk;
            d=v[ix]; v[ix]=v[iy]; v[iy]=d;
          }
        if (kk<m)
          for (i=1; i<=m; i++)
          { 
            ix=(i-1)*m+kk-1; iy=(i-1)*m+kk;
            d=u[ix]; u[ix]=u[iy]; u[iy]=d;
          }
          kk=kk+1;
      }
      it=60;
      mm=mm-1;
    }
    else
    { 
      ks=mm;
      while ((ks>kk)&&(fabs(s[ks-1])!=0.0))
      { 
        d=0.0;
        if (ks!=mm) d=d+fabs(e[ks-1]);
        if (ks!=kk+1) d=d+fabs(e[ks-2]);
        dd=fabs(s[ks-1]);
        if (dd>eps*d) ks=ks-1;
        else s[ks-1]=0.0;
      }
      if (ks==kk)
      { 
        kk=kk+1;
        d=fabs(s[mm-1]);
        t=fabs(s[mm-2]);
        if (t>d) d=t;
        t=fabs(e[mm-2]);
        if (t>d) d=t;
        t=fabs(s[kk-1]);
        if (t>d) d=t;
        t=fabs(e[kk-1]);
        if (t>d) d=t;
        sm=s[mm-1]/d; sm1=s[mm-2]/d;
        em1=e[mm-2]/d;
        sk=s[kk-1]/d; ek=e[kk-1]/d;
        b=((sm1+sm)*(sm1-sm)+em1*em1)/2.0;
        c=sm*em1; c=c*c; shh=0.0;
        if ((b!=0.0)||(c!=0.0))
        { 
          shh=sqrt(b*b+c);
          if (b<0.0) shh=-shh;
          shh=c/(b+shh);
        }
        fg[0]=(sk+sm)*(sk-sm)-shh;
        fg[1]=sk*ek;
        for (i=kk; i<=mm-1; i++)
        { 
          sss(fg,cs);
          if (i!=kk) e[i-2]=fg[0];
          fg[0]=cs[0]*s[i-1]+cs[1]*e[i-1];
          e[i-1]=cs[0]*e[i-1]-cs[1]*s[i-1];
          fg[1]=cs[1]*s[i];
          s[i]=cs[0]*s[i];
          if ((cs[0]!=1.0)||(cs[1]!=0.0))
            for (j=1; j<=n; j++)
            { 
               ix=(j-1)*n+i-1;
               iy=(j-1)*n+i;
               d=cs[0]*v[ix]+cs[1]*v[iy];
               v[iy]=-cs[1]*v[ix]+cs[0]*v[iy];
               v[ix]=d;
            }
          sss(fg,cs);
          s[i-1]=fg[0];
          fg[0]=cs[0]*e[i-1]+cs[1]*s[i];
          s[i]=-cs[1]*e[i-1]+cs[0]*s[i];
          fg[1]=cs[1]*e[i];
          e[i]=cs[0]*e[i];
          if (i<m)
            if ((cs[0]!=1.0)||(cs[1]!=0.0))
              for (j=1; j<=m; j++)
              { 
                ix=(j-1)*m+i-1;
                iy=(j-1)*m+i;
                d=cs[0]*u[ix]+cs[1]*u[iy];
                u[iy]=-cs[1]*u[ix]+cs[0]*u[iy];
                u[ix]=d;
              }
        }
        e[mm-2]=fg[0];
        it=it-1;
      }
      else
      { 
        if (ks==mm)
        { 
          kk=kk+1;
          fg[1]=e[mm-2]; e[mm-2]=0.0;
          for (ll=kk; ll<=mm-1; ll++)
          { 
            i=mm+kk-ll-1;
            fg[0]=s[i-1];
            sss(fg,cs);
            s[i-1]=fg[0];
            if (i!=kk)
            { 
              fg[1]=-cs[1]*e[i-2];
              e[i-2]=cs[0]*e[i-2];
            }
            if ((cs[0]!=1.0)||(cs[1]!=0.0))
              for (j=1; j<=n; j++)
              { 
                ix=(j-1)*n+i-1;
                iy=(j-1)*n+mm-1;
                d=cs[0]*v[ix]+cs[1]*v[iy];
                v[iy]=-cs[1]*v[ix]+cs[0]*v[iy];
                v[ix]=d;
              }
          }
        }
        else
        { 
          kk=ks+1;
          fg[1]=e[kk-2];
          e[kk-2]=0.0;
          for (i=kk; i<=mm; i++)
          { 
            fg[0]=s[i-1];
            sss(fg,cs);
            s[i-1]=fg[0];
            fg[1]=-cs[1]*e[i-1];
            e[i-1]=cs[0]*e[i-1];
            if ((cs[0]!=1.0)||(cs[1]!=0.0))
              for (j=1; j<=m; j++)
              { 
                ix=(j-1)*m+i-1;
                iy=(j-1)*m+kk-2;
                d=cs[0]*u[ix]+cs[1]*u[iy];
                u[iy]=-cs[1]*u[ix]+cs[0]*u[iy];
                u[ix]=d;
              }
          }
        }
      }
    }
  }
  return(1);
}

void BasicFunc::ppp(
  double a[],
  double e[],
  double s[],
  double v[],
  int m,
  int n
)
{
  int i,j,p,q;
  double d;
  if (m>=n) i=n;
  else i=m;
  for (j=1; j<=i-1; j++)
  {           
    a[(j-1)*n+j-1]=s[j-1];
    a[(j-1)*n+j]=e[j-1];
  }
  a[(i-1)*n+i-1]=s[i-1];
  if (m<n) a[(i-1)*n+i]=e[i-1];
  for (i=1; i<=n-1; i++)
  {           
    for (j=i+1; j<=n; j++)
    {       
      p=(i-1)*n+j-1;
      q=(j-1)*n+i-1;
      d=v[p]; 
      v[p]=v[q];
      v[q]=d;
    }
  }

  return;
}

void BasicFunc::sss(
  double fg[2],
  double cs[2]
)
{ 
  double r,d;
  if ((fabs(fg[0])+fabs(fg[1]))==0.0)
  {
    cs[0]=1.0;
    cs[1]=0.0;
    d=0.0;
  }
  else
  {
    d=sqrt(fg[0]*fg[0]+fg[1]*fg[1]);
    if (fabs(fg[0])>fabs(fg[1]))
    {
      d=fabs(d);
      if (fg[0]<0.0)
      {
        d=-d;
      }
    }
    if (fabs(fg[1])>=fabs(fg[0]))
    {
      d=fabs(d);
      if (fg[1]<0.0)
      {
        d=-d;
      }
    }
    cs[0]=fg[0]/d;
    cs[1]=fg[1]/d;
  }
  r=1.0;
  if (fabs(fg[0])>fabs(fg[1]))
  {
    r=cs[1];
  }
  else
  {
    if (cs[0]!=0.0)
    {
      r=1.0/cs[0];
    }
  }
  fg[0]=d; fg[1]=r;

  return;
}

