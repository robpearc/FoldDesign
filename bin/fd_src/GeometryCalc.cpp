///////////////////////////////////////////////////////////////////////////////////////////////
//         Copyright (C) Regents of the University of Michigan - All Rights Reserved         //
//          Unauthorized copying of this file, via any medium is strictly prohibited         //
//                                                                                           //
//                          Author: Robin Pearce <robpearc@umich.edu>                        //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "BasicFunc.h"
#include "GeometryCalc.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

GeometryCalc::GeometryCalc()
{
  sgpos=NULL;
  sgpos2=NULL;
  flagVerbose=true;
}

GeometryCalc::~GeometryCalc()
{
}

bool GeometryCalc::configureGeometry(
  const ParseGeometryFiles &geometryParameters
)
{
  sgpos=geometryParameters.sgpos;
}

/* Convert coordinate from torsion angle space to cartesian space
 * type - atom types for coordinate update
 *        1 and 6: CA only
 *        3: N, CA, C */
bool GeometryCalc::tor2str(
  point3f *decstr,
  int seqnum,
  int type
)
{
  int i;
  point3s pt,pn,pc;
  BasicFunc bf;
  bool flagts;
  bool flagWhole=true;
  if(type==1) //all 0 n 1 ca 2 c
  {
    if(decstr[1].leng<0) decstr[1].leng=float(lennca);
    if(decstr[2].leng<0) decstr[2].leng=float(lencac);
    if(decstr[2].angl<0) decstr[2].angl=float(angncac);
    decstr[0].x=0;
    decstr[0].y=0;
    decstr[0].z=0;
    decstr[1].x=decstr[1].leng;
    decstr[1].y=0;
    decstr[1].z=0;
    decstr[2].x=decstr[1].leng-decstr[2].leng*cos(decstr[2].angl*raddeg);
    decstr[2].y=decstr[2].leng*sin(decstr[2].angl*raddeg);
    decstr[2].z=0;
    for(i=3;i<seqnum;i++)
    {
      flagts=bf.tor2pos22(decstr[i-3].x,decstr[i-3].y,decstr[i-3].z,
                          decstr[i-2].x,decstr[i-2].y,decstr[i-2].z,
                          decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
                          decstr[i].phi*raddeg,decstr[i].leng,
                          decstr[i].angl*raddeg,&pt.x,&pt.y,&pt.z);
      if(!flagts)
      {
        flagWhole=false;
        if(flagVerbose)
        {
          printf("Wrong front coordinates %d\n",i);
        }
      }
      decstr[i].x=pt.x;
      decstr[i].y=pt.y;
      decstr[i].z=pt.z;
    }
  }
  else if(type==6) //Ca 
  {
    if(decstr[1].leng<0) decstr[1].leng=float(lencaca);
    if(decstr[2].leng<0) decstr[2].leng=float(lencaca);
    if(decstr[2].angl<0) decstr[2].angl=106.422f;
    decstr[0].x=0;
    decstr[0].y=0;
    decstr[0].z=0;
    decstr[1].x=decstr[1].leng;
    decstr[1].y=0;
    decstr[1].z=0;
    decstr[2].x=decstr[1].leng-decstr[2].leng*cos(decstr[2].angl*raddeg);
    decstr[2].y=decstr[2].leng*sin(decstr[2].angl*raddeg);
    decstr[2].z=0;
    for(i=3;i<seqnum;i++)
    {
      flagts=bf.tor2pos22(decstr[i-3].x,decstr[i-3].y,decstr[i-3].z,
                          decstr[i-2].x,decstr[i-2].y,decstr[i-2].z,
                          decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
                          decstr[i].phi*raddeg,decstr[i].leng,
                          decstr[i].angl*raddeg,&pt.x,&pt.y,&pt.z);
      if(!flagts)
      {
        flagWhole=false;
        if(flagVerbose)
        {
          printf("Wrong front coordinates %d\n",i);
        }
      }
      decstr[i].x=pt.x;
      decstr[i].y=pt.y;
      decstr[i].z=pt.z;
    }
  }
  else if(type==3)
  {
    if(decstr[0].len[1]<0) decstr[0].len[1]=float(lennca);
    if(decstr[0].len[2]<0) decstr[0].len[2]=float(lencac);
    if(decstr[0].ang[2]<0) decstr[0].ang[2]=float(angncac);
    decstr[0].ptn.x=0;
    decstr[0].ptn.y=0;
    decstr[0].ptn.z=0;
    decstr[0].x=decstr[0].len[1];
    decstr[0].y=0;
    decstr[0].z=0;
    decstr[0].ptc.x=decstr[0].len[1]-decstr[0].len[2]*cos(decstr[0].ang[2]*raddeg);
    decstr[0].ptc.y=decstr[0].len[2]*sin(decstr[0].ang[2]*raddeg);
    decstr[0].ptc.z=0;
    if(decstr[0].tor[0]>=0 && decstr[0].tor[2]>=0)
    {
      flagts=bf.tor2pos22(0,0,0,
                          lennca,0,0,
                          lennca+lencac*sin(angncac*raddeg),-lencac*cos(angncac*raddeg),0,
                          decstr[0].tor[0]*raddeg,float(lencn),
                          float(angcacn)*raddeg,&pn.x,&pn.y,&pn.z);
      if(!flagts)
      {
        flagWhole=false;
        if(flagVerbose)
        {
          printf("Wrong front coordinates n %d\n",0);
        }
      }
      decstr[0].ptn.x=pn.x;
      decstr[0].ptn.y=pn.y;
      decstr[0].ptn.z=pn.z;
      if(decstr[0].tor[1]<0) decstr[0].tor[1]=180;
      flagts=bf.tor2pos22(lennca,0,0,
                          lennca+lencac*sin(angncac*raddeg),-lencac*cos(angncac*raddeg),0,
                          decstr[0].ptn.x,decstr[0].ptn.y,decstr[0].ptn.z,
                          decstr[0].tor[1]*raddeg,float(lennca),
                          angcnca*raddeg,&pt.x,&pt.y,&pt.z);
      if(!flagts)
      {
        flagWhole=false;
        if(flagVerbose)
        {
          printf("Wrong front coordinates ca %d\n",0);
        }
      }
      decstr[0].x=pt.x;
      decstr[0].y=pt.y;
      decstr[0].z=pt.z;
      flagts=bf.tor2pos22(lennca+lencac*sin(angncac*raddeg),-lencac*cos(angncac*raddeg),0,
                          decstr[0].ptn.x,decstr[0].ptn.y,decstr[0].ptn.z,
                          decstr[0].x,decstr[0].y,decstr[0].z,
                          decstr[0].tor[2]*raddeg,float(lencac),
                          angncac*raddeg,&pc.x,&pc.y,&pc.z);
      if(!flagts)
      {
        flagWhole=false;
        if(flagVerbose)
        {
          printf("Wrong front coordinates c %d\n",0);
        }
      }
      decstr[0].ptc.x=pc.x;
      decstr[0].ptc.y=pc.y;
      decstr[0].ptc.z=pc.z;
    }
    for(i=1;i<seqnum;i++)
    {
      if(decstr[i].tor[0]<0) 
      {
        if(decstr[i].stype=='E') decstr[i].tor[0]=120;
        else if(decstr[i].stype=='H') decstr[i].tor[0]=300;
        else decstr[i].tor[0]=120;
      }
      if(decstr[i].tor[1]<0) decstr[i].tor[1]=180;
      if(decstr[i].tor[2]<0) decstr[i].tor[2]=290;
      if(decstr[i].len[0]<0) decstr[i].len[0]=float(lencn);
      if(decstr[i].len[1]<0) decstr[i].len[1]=float(lennca);
      if(decstr[i].len[2]<0) decstr[i].len[2]=float(lencac);
      if(decstr[i].ang[0]<0) decstr[i].ang[0]=float(angcacn);
      if(decstr[i].ang[1]<0) decstr[i].ang[1]=float(angcnca);
      if(decstr[i].ang[2]<0) decstr[i].ang[2]=float(angncac);
      //0 1 2 n ca c
      //original phi i-1 psi i-1 omega i
      flagts=bf.tor2pos22(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,
                          decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
                          decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                          decstr[i].tor[0]*raddeg,decstr[i].len[0],
                          decstr[i].ang[0]*raddeg,&pn.x,&pn.y,&pn.z);
      if(!flagts)
      {
        flagWhole=false;
        if(flagVerbose)
        {
          printf("wrong front coordinates n %d %f %f %f\n",
                 i,decstr[i].tor[0],decstr[i].len[0],decstr[i].ang[0]);
        }
      }
      decstr[i].ptn.x=pn.x;
      decstr[i].ptn.y=pn.y;
      decstr[i].ptn.z=pn.z;
      flagts=bf.tor2pos22(decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
                          decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                          decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                          decstr[i].tor[1]*raddeg,decstr[i].len[1],
                          decstr[i].ang[1]*raddeg,&pt.x,&pt.y,&pt.z);
      if(!flagts)
      {
        flagWhole=false;
        if(flagVerbose)
        {
          printf("wrong front coordinates ca %d %f %f %f\n",
                 i,decstr[i].tor[1],decstr[i].len[1],decstr[i].ang[1]);
        }
      }
      decstr[i].x=pt.x;
      decstr[i].y=pt.y;
      decstr[i].z=pt.z;
      flagts=bf.tor2pos22(decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                          decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                          decstr[i].x,decstr[i].y,decstr[i].z,
                          decstr[i].tor[2]*raddeg,decstr[i].len[2],
                          decstr[i].ang[2]*raddeg,&pc.x,&pc.y,&pc.z);
      if(!flagts)
      {
        flagWhole=false;
        if(flagVerbose)
        {
          printf("wrong front coordinates c %d %f %f %f\n",
                 i,decstr[i].tor[2],decstr[i].len[2],decstr[i].ang[2]);
        }
      }
      decstr[i].ptc.x=pc.x;
      decstr[i].ptc.y=pc.y;
      decstr[i].ptc.z=pc.z;
    }
  }
 
  return flagWhole;
}

/* Convert coordinate from torsion-angle space to cartesian space */
bool GeometryCalc::tor2strp(
  point3f *decstr,
  int seqnum,
  int istart
)
{
  int i;
  point3s pt,pn,pc;
  BasicFunc bf;
  bool flagts;
  bool flagWhole=true;
  for(i=istart;i<seqnum;i++)
  {
    if(decstr[i].tor[0]<0) 
    {
      if(decstr[i].stype=='E')      decstr[i].tor[0]=120;
      else if(decstr[i].stype=='H') decstr[i].tor[0]=300;
      else                          decstr[i].tor[0]=120;
    }
    if(decstr[i].tor[1]<0) decstr[i].tor[1]=180;
    if(decstr[i].tor[2]<0) decstr[i].tor[2]=290;
    if(decstr[i].len[0]<0) decstr[i].len[0]=float(lencn);
    if(decstr[i].len[1]<0) decstr[i].len[1]=float(lennca);
    if(decstr[i].len[2]<0) decstr[i].len[2]=float(lencac);
    if(decstr[i].ang[0]<0) decstr[i].ang[0]=float(angcacn);
    if(decstr[i].ang[1]<0) decstr[i].ang[1]=float(angcnca);
    if(decstr[i].ang[2]<0) decstr[i].ang[2]=float(angncac);
    //0 1 2 n ca c
    //original phi i-1 psi i-1 omega i
    flagts=bf.tor2pos22(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,
                        decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
                        decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                        decstr[i].tor[0]*raddeg,decstr[i].len[0],
                        decstr[i].ang[0]*raddeg,&pn.x,&pn.y,&pn.z);
    if(!flagts)
    {
      flagWhole=false;
      if(flagVerbose)
      {
        printf("Wrong front coordinatesp n %d %f %f %f\n",
               i,decstr[i].tor[0],decstr[i].len[0],decstr[i].ang[0]);
      } 
    }
    decstr[i].ptn.x=pn.x;
    decstr[i].ptn.y=pn.y;
    decstr[i].ptn.z=pn.z;
    //Ca
    flagts=bf.tor2pos22(decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
                        decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                        decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                        decstr[i].tor[1]*raddeg,decstr[i].len[1],
                        decstr[i].ang[1]*raddeg,&pt.x,&pt.y,&pt.z);
    if(!flagts)
    {
      flagWhole=false;
      if(flagVerbose)
      {
        printf("Wrong front coordinatesp ca %d %f %f %f\n",
               i,decstr[i].tor[1],decstr[i].len[1],decstr[i].ang[1]);
      }
    }
    decstr[i].x=pt.x;
    decstr[i].y=pt.y;
    decstr[i].z=pt.z;
    //c
    flagts=bf.tor2pos22(decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                        decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                        decstr[i].x,decstr[i].y,decstr[i].z,
                        decstr[i].tor[2]*raddeg,decstr[i].len[2],
                        decstr[i].ang[2]*raddeg,&pc.x,&pc.y,&pc.z);
    if(!flagts)
    {
      flagWhole=false;
      if(flagVerbose)
      {
        printf("Wrong front coordinatesp c %d %f %f %f\n",
               i,decstr[i].tor[2],decstr[i].len[2],decstr[i].ang[2]);
      }
    }
    decstr[i].ptc.x=pc.x;
    decstr[i].ptc.y=pc.y;
    decstr[i].ptc.z=pc.z;
  }

  return flagWhole;
}

/* Convert coordinate from cartesian space to torsion-angle space */
void GeometryCalc::str2torp(
  point3f *decstr,
  int seqnum, //Not used
  int iStart,
  int iEnd
)
{
  int i;
  point3d p12,p23;
  BasicFunc bf;
  int realStart;
  if(iStart==0) realStart=1;
  else realStart=iStart;
  if(realStart==1)
  {
    if(decstr[0].tor[0]<0 || decstr[0].len[0]<0)
    {
      decstr[0].tor[0]=180.0;
      decstr[0].tor[1]=180.0;
      decstr[0].tor[2]=180.0;
      decstr[0].len[0]=lennc;
      decstr[0].ang[0]=angcacn;
      decstr[0].ang[1]=angcnca;
    }
    p12=bf.setv(decstr[0].ptn.x-decstr[0].x,
                decstr[0].ptn.y-decstr[0].y,
                decstr[0].ptn.z-decstr[0].z);
    p23=bf.setv(decstr[0].ptc.x-decstr[0].x,
                decstr[0].ptc.y-decstr[0].y,
                decstr[0].ptc.z-decstr[0].z);
    decstr[0].len[1]=bf.norm(p12);
    decstr[0].len[2]=bf.norm(p23);
    decstr[0].ang[2]=bf.angv(p12,p23)*degrad;
  }
  for(i=realStart;i<=iEnd;i++)
  {
    p12=bf.setv(decstr[i-1].x-decstr[i-1].ptc.x,
                decstr[i-1].y-decstr[i-1].ptc.y,
                decstr[i-1].z-decstr[i-1].ptc.z);
    p23=bf.setv(decstr[i].ptn.x-decstr[i-1].ptc.x,
                decstr[i].ptn.y-decstr[i-1].ptc.y,
                decstr[i].ptn.z-decstr[i-1].ptc.z);
    decstr[i].len[0]=bf.norm(p23);
    decstr[i].ang[0]=bf.angv(p12,p23)*degrad;
    decstr[i].tor[0]=bf.phi(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,
                            decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
                            decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                            decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z);
    p12=bf.setv(decstr[i-1].ptc.x-decstr[i].ptn.x,
                decstr[i-1].ptc.y-decstr[i].ptn.y,
                decstr[i-1].ptc.z-decstr[i].ptn.z);
    p23=bf.setv(decstr[i].x-decstr[i].ptn.x,
                decstr[i].y-decstr[i].ptn.y,
                decstr[i].z-decstr[i].ptn.z);
    decstr[i].len[1]=bf.norm(p23);
    decstr[i].ang[1]=bf.angv(p12,p23)*degrad;
    decstr[i].tor[1]=bf.phi(decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
                            decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                            decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                            decstr[i].x,decstr[i].y,decstr[i].z);
    p12=bf.setv(decstr[i].ptn.x-decstr[i].x,
                decstr[i].ptn.y-decstr[i].y,
                decstr[i].ptn.z-decstr[i].z);
    p23=bf.setv(decstr[i].ptc.x-decstr[i].x,
                decstr[i].ptc.y-decstr[i].y,
                decstr[i].ptc.z-decstr[i].z);
    decstr[i].len[2]=bf.norm(p23);
    decstr[i].ang[2]=bf.angv(p12,p23)*degrad;
    decstr[i].tor[2]=bf.phi(decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                            decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                            decstr[i].x,decstr[i].y,decstr[i].z,
                            decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z);
  }
}

void GeometryCalc::str2tor(
  point3f *decstr,
  int seqnum,
  int type
)
{
  int i;
  point3d p12,p23;
  BasicFunc bf;
  if(type==1 || type==6)
  {
    if(decstr[0].leng<0 || decstr[0].phi<0)
    {
      if(type==1)
      {
        decstr[0].leng=lennc;
        decstr[0].angl=angcacn;
        decstr[1].angl=angcnca;
      }
      else
      {
        decstr[0].leng=lencaca;
        decstr[0].angl=106.422f;
        decstr[1].angl=106.422f;
      }
      decstr[0].phi=180.0;
      decstr[1].phi=180.0;
      decstr[2].phi=180.0;
    }

    p12=bf.setv(decstr[0].x-decstr[1].x,decstr[0].y-decstr[1].y,decstr[0].z-decstr[1].z);
    p23=bf.setv(decstr[2].x-decstr[1].x,decstr[2].y-decstr[1].y,decstr[2].z-decstr[1].z);
    decstr[1].leng=bf.norm(p12);
    decstr[2].leng=bf.norm(p23);
    decstr[2].angl=bf.angv(p12,p23);
    for(i=3;i<seqnum;i++)
    {
      p12=bf.setv(decstr[i-2].x-decstr[i-1].x,
                  decstr[i-2].y-decstr[i-1].y,
                  decstr[i-2].z-decstr[i-1].z);
      p23=bf.setv(decstr[i].x-decstr[i-1].x,
                  decstr[i].y-decstr[i-1].y,
                  decstr[i].z-decstr[i-1].z);
      decstr[i].leng=bf.norm(p23);
      decstr[i].angl=bf.angv(p12,p23);
      decstr[i].phi=bf.phi(decstr[i-3].x,decstr[i-3].y,decstr[i-3].z,
                           decstr[i-2].x,decstr[i-2].y,decstr[i-2].z,
                           decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
                           decstr[i].x,decstr[i].y,decstr[i].z);
    }
  }
  else if(type==3)
  {
    if(decstr[0].tor[0]<0 || decstr[0].len[0]<0)
    {
      decstr[0].tor[0]=180.0;
      decstr[0].tor[1]=180.0;
      decstr[0].tor[2]=180.0;
      decstr[0].len[0]=lennc;
      decstr[0].ang[0]=angcacn;	
      decstr[0].ang[1]=angcnca;
    }
		
    p12=bf.setv(decstr[0].ptn.x-decstr[0].x,
                decstr[0].ptn.y-decstr[0].y,
                decstr[0].ptn.z-decstr[0].z);
    p23=bf.setv(decstr[0].ptc.x-decstr[0].x,
                decstr[0].ptc.y-decstr[0].y,
                decstr[0].ptc.z-decstr[0].z);
    decstr[0].len[1]=bf.norm(p12);
    decstr[0].len[2]=bf.norm(p23);
    decstr[0].ang[2]=bf.angv(p12,p23)*degrad;
    for(i=1;i<seqnum;i++)
    {
      p12=bf.setv(decstr[i-1].x-decstr[i-1].ptc.x,
                  decstr[i-1].y-decstr[i-1].ptc.y,
                  decstr[i-1].z-decstr[i-1].ptc.z);
      p23=bf.setv(decstr[i].ptn.x-decstr[i-1].ptc.x,
                  decstr[i].ptn.y-decstr[i-1].ptc.y,
                  decstr[i].ptn.z-decstr[i-1].ptc.z);
      decstr[i].len[0]=bf.norm(p23);
      decstr[i].ang[0]=bf.angv(p12,p23)*degrad;
      decstr[i].tor[0]=bf.phi(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,
                              decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
                              decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                              decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z);
			
      p12=bf.setv(decstr[i-1].ptc.x-decstr[i].ptn.x,
                  decstr[i-1].ptc.y-decstr[i].ptn.y,
                  decstr[i-1].ptc.z-decstr[i].ptn.z);
      p23=bf.setv(decstr[i].x-decstr[i].ptn.x,
                  decstr[i].y-decstr[i].ptn.y,
                  decstr[i].z-decstr[i].ptn.z);
      decstr[i].len[1]=bf.norm(p23);
      decstr[i].ang[1]=bf.angv(p12,p23)*degrad;
      decstr[i].tor[1]=bf.phi(decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
                              decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                              decstr[i].ptn.x,decstr[i].ptn.y,
                              decstr[i].ptn.z,decstr[i].x,decstr[i].y,decstr[i].z);
		
      p12=bf.setv(decstr[i].ptn.x-decstr[i].x,decstr[i].ptn.y-decstr[i].y,decstr[i].ptn.z-decstr[i].z);
      p23=bf.setv(decstr[i].ptc.x-decstr[i].x,decstr[i].ptc.y-decstr[i].y,decstr[i].ptc.z-decstr[i].z);
      decstr[i].len[2]=bf.norm(p23);
      decstr[i].ang[2]=bf.angv(p12,p23)*degrad;
      decstr[i].tor[2]=bf.phi(decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                              decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,decstr[i].x,
                              decstr[i].y,decstr[i].z,decstr[i].ptc.x,decstr[i].ptc.y,
                              decstr[i].ptc.z);
    }
  }	
}

bool GeometryCalc::tor2strsg(
  point3f *decstr,
  int seqnum
)
{
  int i;
  int tind;
  point3s pn;
  BasicFunc bf;
  bool flagts;
  bool flagWhole=true;
  for(i=1;i<seqnum;i++)
  {
    //atom o
    flagts=bf.tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                        decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
                        decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                        179.6715f*raddeg,1.229f,2.0961f,&pn.x,&pn.y,&pn.z);
    decstr[i-1].pto.x=pn.x;
    decstr[i-1].pto.y=pn.y;
    decstr[i-1].pto.z=pn.z;
    //atom h
    flagts=bf.tor2pos22(decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                        decstr[i].x,decstr[i].y,decstr[i].z,
                        decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                        179.8174f*raddeg,0.987f,2.0814f,&pn.x,&pn.y,&pn.z);//0.9919f,2.0574f
    if(!flagts)
    {
      flagWhole=false;
      if(flagVerbose)
      {
        printf("Wrong front coordinates d %d\n",i);
      }
    }
    decstr[i].pth.x=pn.x;
    decstr[i].pth.y=pn.y;
    decstr[i].pth.z=pn.z;
  }
  //atom o
  i=seqnum;
  flagts=bf.tor2pos22(decstr[i-1].ptn.x,decstr[i-1].ptn.y,decstr[i-1].ptn.z,
                      decstr[i-1].x,decstr[i-1].y,decstr[i-1].z,
                      decstr[i-1].ptc.x,decstr[i-1].ptc.y,decstr[i-1].ptc.z,
                      0.0f,1.2439f, 2.0855f,&pn.x,&pn.y,&pn.z);
  if(!flagts)
  {
    flagWhole=false;
    if(flagVerbose)
    {
      printf("Wrong front coordinates o %d\n",i);
    }
  }
  decstr[i-1].pto.x=pn.x;
  decstr[i-1].pto.y=pn.y;
  decstr[i-1].pto.z=pn.z;
  //atom h 
  i=0;
  flagts=bf.tor2pos22(decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
                      decstr[i].x,decstr[i].y,decstr[i].z,
                      decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                      60*raddeg,0.987f,2.0306f,&pn.x,&pn.y,&pn.z);//0.9972f
  if(!flagts)
  {
    flagWhole=false;
    if(flagVerbose)
    {
      printf("Wrong front coordinates %d\n",i);
    }
  }
  decstr[i].pth.x=pn.x;
  decstr[i].pth.y=pn.y;
  decstr[i].pth.z=pn.z;
  //atom cb new
  for(i=0;i<seqnum;i++)
  {
    tind=bf.aminoid(decstr[i].aaa);
    if(tind>19) tind=19;
    flagts=bf.tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                        decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
                        decstr[i].x,decstr[i].y,decstr[i].z,
                        cbsta[tind][2]*raddeg,cbsta[tind][0],cbsta[tind][1],
                        &pn.x,&pn.y,&pn.z);
    if(!flagts)
    {
      flagWhole=false;
      if(flagVerbose)
      {
        printf("Wrong front coordinates cb2 %d\n",i);
      }
    }
    decstr[i].ptb.x=pn.x;
    decstr[i].ptb.y=pn.y;
    decstr[i].ptb.z=pn.z;
    if(decstr[i].aaa=='G')
    {
      decstr[i].ptb.x=decstr[i].x;
      decstr[i].ptb.y=decstr[i].y;
      decstr[i].ptb.z=decstr[i].z;
    }
  }

  return flagWhole;
}

bool GeometryCalc::tor2strsg2(
  point3f *decstr,
  int seqnum
)
{
  bool flagWhole=tor2strsg(decstr,seqnum);
  int i;
  int tind;
  point3s pn;
  BasicFunc bf;
  bool flagts;

  //atom sg
  int ti,tj;
  int cutnum=72;
  double delta=5.0;
  for(i=0;i<seqnum;i++)
  {
    if(decstr[i].aaa=='G')
    {
      decstr[i].ptsg.x=decstr[i].x;
      decstr[i].ptsg.y=decstr[i].y;
      decstr[i].ptsg.z=decstr[i].z;
      continue;
    }
    tind=bf.aminoid(decstr[i].aaa);
    if(tind>19) tind=19;
    if(i<seqnum-1)
    {
      ti=int(decstr[i].tor[2]/delta);
      tj=int(decstr[i+1].tor[0]/delta);
      if(ti>=0 && ti<cutnum && tj>=0 && tj<cutnum)
      {
        flagts=bf.tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                            decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
                            decstr[i].x,decstr[i].y,decstr[i].z,
                            sgpos[2][tind][ti][tj],sgpos[0][tind][ti][tj],
                            sgpos[1][tind][ti][tj],&pn.x,&pn.y,&pn.z);
        if(!flagts)
        {
          flagWhole=false;
          //printf("wrong front coordinates full sg %d\n",i);
        }
        decstr[i].ptsg.x=pn.x;
        decstr[i].ptsg.y=pn.y;
        decstr[i].ptsg.z=pn.z;
        continue;
      }
    }
    flagts=bf.tor2pos22(decstr[i].ptn.x,decstr[i].ptn.y,decstr[i].ptn.z,
                        decstr[i].ptc.x,decstr[i].ptc.y,decstr[i].ptc.z,
                        decstr[i].x,decstr[i].y,decstr[i].z,
                        sglatavg[tind][2]*raddeg,sglatavg[tind][0],sglatavg[tind][1],
                        &pn.x,&pn.y,&pn.z);
    if(!flagts)
    {
      flagWhole=false;
      //printf("wrong front coordinates sg %d\n",i);
    }
    decstr[i].ptsg.x=pn.x;
    decstr[i].ptsg.y=pn.y;
    decstr[i].ptsg.z=pn.z;
  }
  //atom ct
  for(i=0;i<seqnum;i++)
  {
    tind=bf.aminoid(decstr[i].aaa);
    pn.x=0;pn.y=0;pn.z=0;
    pn.x+=decstr[i].ptn.x;
    pn.y+=decstr[i].ptn.y;
    pn.z+=decstr[i].ptn.z;
    pn.x+=decstr[i].x;
    pn.y+=decstr[i].y;
    pn.z+=decstr[i].z;
    pn.x+=decstr[i].ptc.x;
    pn.y+=decstr[i].ptc.y;
    pn.z+=decstr[i].ptc.z;
    pn.x+=decstr[i].pto.x;
    pn.y+=decstr[i].pto.y;
    pn.z+=decstr[i].pto.z;
    pn.x+=decstr[i].ptsg.x*sgatomnum[tind];
    pn.y+=decstr[i].ptsg.y*sgatomnum[tind];
    pn.z+=decstr[i].ptsg.z*sgatomnum[tind];
    decstr[i].ptg.x=pn.x/double(sgatomnum[tind]+4.0);
    decstr[i].ptg.y=pn.y/double(sgatomnum[tind]+4.0);
    decstr[i].ptg.z=pn.z/double(sgatomnum[tind]+4.0);
  }

  return flagWhole;
}

