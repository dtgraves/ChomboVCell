#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include <stdio.h>

#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "IndexTM.H"
#include "WrappedGShop.H"
#include "BoundaryAreaRefCrit.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "PolyGeom.H"
#include "EBAMRIO.H"

#include "PlaneIF.H"
#include "SphereIF.H"
#include "UnionIF.H"
#include "IntersectionIF.H"
#include "ComplementIF.H"

#include "MitochondriaIF.H"

#include "EBMenagerieUtils.H"

#include "UsingNamespace.H"

void getDomainOriginDx(Box      & a_domain,
                       RealVect & a_origin,
                       Real     & a_dx)
{
  ParmParse pp;
  Vector<int> n_cell(SpaceDim);
  pp.getarr("n_cell",n_cell,0,SpaceDim);

  CH_assert(n_cell.size() == SpaceDim);

  IntVect lo = IntVect::Zero;
  IntVect hi;

  for (int ivec = 0; ivec < SpaceDim; ivec++)
  {
    if (n_cell[ivec] <= 0)
    {
      pout() << "Bogus number of cells input = " << n_cell[ivec];
      exit(1);
    }

    hi[ivec] = n_cell[ivec] - 1;
  }

  a_domain.setSmall(lo);
  a_domain.setBig(hi);
  Vector<Real> prob_lo(SpaceDim,1.0);
  Real prob_hi;

  pp.getarr("prob_lo",prob_lo,0,SpaceDim);
  pp.get("prob_hi",prob_hi);

  a_dx = (prob_hi-prob_lo[0])/n_cell[0];

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    a_origin[idir] = prob_lo[idir];
  }
}

BaseIF* makeMultiSphereGeometry(MFIndexSpace   & a_mfIndexSpace,
                                const Box      & a_domain,
                                const RealVect & a_origin,
                                const Real     & a_dx)
{
  // parse input file
  ParmParse pp;

  int numSpheres;
  pp.get("num_spheres",numSpheres);

  Vector<RealVect> center(numSpheres);
  Vector<Real>     radius(numSpheres,-1.0);
  Vector<int>      phaseInside(numSpheres,0);
  Vector<int>      used(numSpheres,0);

  for (int i = 0; i < numSpheres; i++)
  {
    char centerStr[100];
    sprintf(centerStr,"sphere_center_%d",i+1);

    // ParmParse doesn't get RealVects, so work-around with Vector<Real>
    Vector<Real> vectorCenter;
    pp.getarr(centerStr,vectorCenter,0,SpaceDim);
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      center[i][idir] = vectorCenter[idir];
    }

    char radiusStr[100];
    sprintf(radiusStr,"sphere_radius_%d",i+1);
    pp.get(radiusStr,radius[i]);

    char phaseInsideStr[100];
    sprintf(phaseInsideStr,"sphere_phase_inside_%d",i+1);
    pp.get(phaseInsideStr,phaseInside[i]);
  }

  Vector<BaseIF*> spheres(numSpheres);

  for (int i = 0; i < numSpheres; i++)
  {
    bool inside = (phaseInside[i] == 0);
    spheres[i] = new SphereIF(radius[i],center[i],inside);
  }

  int numIntersectionLists;
  pp.get("num_intersection_lists",numIntersectionLists);

  Vector<Vector<int> > intersectionLists(numIntersectionLists);
  for (int i = 0; i < numIntersectionLists; i++)
  {
    char intersectionListNumStr[100];
    sprintf(intersectionListNumStr,"intersection_list_num_%d",i+1);

    int intersectionListNum;
    pp.get(intersectionListNumStr,intersectionListNum);

    intersectionLists[i].resize(intersectionListNum);

    char intersectionListStr[100];
    sprintf(intersectionListStr,"intersection_list_%d",i+1);

    pp.getarr(intersectionListStr,intersectionLists[i],0,intersectionListNum);

    for (int j = 0; j < intersectionListNum; j++)
    {
      intersectionLists[i][j] -= 1;
    }
  }

  Vector<BaseIF*> unionList;

  for (int i = 0; i < numIntersectionLists; i++)
  {
    Vector<BaseIF*> curSpheres(intersectionLists[i].size());

    for (int j = 0; j < intersectionLists[i].size(); j++)
    {
      if (used[intersectionLists[i][j]] == 0)
      {
        curSpheres[j] = spheres[intersectionLists[i][j]];
        used[intersectionLists[i][j]] = 1; 
      }
      else
      {
        MayDay::Error("Sphere used on more than one intersection list");
      }
    }

    IntersectionIF* curIntersection = new IntersectionIF(curSpheres);
    unionList.push_back(curIntersection);
  }

  for (int i = 0; i < numSpheres; i++)
  {
    if (used[i] == 0)
    {
      unionList.push_back(spheres[i]);
      used[i] = 1;
    }
  }

  RefCountedPtr<BaseIF> everythingPhase0( new UnionIF(unionList));


  // Complement for MF
  bool complement = true;
  RefCountedPtr<BaseIF> everythingPhase1(new ComplementIF(*everythingPhase0,complement));

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  Vector<GeometryService*> geometries(2,NULL);

  bool useNewGeometry = true;
  pp.query("use_new_geometry_gen", useNewGeometry);
  
  if(useNewGeometry)
    {
      Real threshold = sqrt(2.);
      pp.query("ba_threshold", threshold);
      pout() << "using boundaryarea wrapped gshop refinement criterion with threshold ba frac=" <<  threshold << endl;
      RefCountedPtr<WGSRefinementCriterion> refcrit(new BoundaryAreaRefCrit(threshold));
      pout() << "using new wrapped gshop geometry generation" << endl;
      int minRefine = 0;
      int maxRefine = 3;
      WrappedGShop* workshop0 = new WrappedGShop(everythingPhase0,RealVect::Zero, a_dx, a_domain, minRefine, maxRefine);
      workshop0->setRefinementCriterion(refcrit);
      workshop0->m_phase = 0;
      geometries[0] = workshop0;

      WrappedGShop* workshop1 = new WrappedGShop(everythingPhase1,RealVect::Zero, a_dx, a_domain, minRefine, maxRefine);
      workshop0->setRefinementCriterion(refcrit);
      workshop1->m_phase = 1;
      geometries[1] = workshop1;
    }
  else
    {
      pout() << "using old style geometry generation" << endl;
      GeometryShop* workshop0 = new GeometryShop(*everythingPhase0,0,vectDx);
      workshop0->m_phase = 0;
      geometries[0] = workshop0;

      GeometryShop* workshop1 = new GeometryShop(*everythingPhase1,0,vectDx);
      workshop1->m_phase = 1;
      geometries[1] = workshop1;
    }
  int maxBoxSize;
  pp.get("maxboxsize",maxBoxSize);

  // int maxCoarsenings = 0;

  // This generates the new EBIS
  a_mfIndexSpace.define(a_domain,a_origin,a_dx,geometries,maxBoxSize/* ,maxCoarsenings */);

  for (int i = 0; i < numIntersectionLists; i++)
  {
    delete unionList[i];
  }

  for (int i = 0; i < numSpheres; i++)
  {
    delete spheres[i];
  }

  BaseIF* retval =  everythingPhase0->newImplicitFunction();
  delete geometries[0];
  delete geometries[1];

  return retval;

}

BaseIF* makeCalciumOneGeometry(MFIndexSpace   & a_mfIndexSpace,
                               const Box      & a_domain,
                               const RealVect & a_origin,
                               const Real     & a_dx)
{
  // parse input file
  ParmParse pp;

  Real radius = 3;
  bool inside = true;

  RealVect centerL(D_DECL(-9.40,6.65,0.00));
  SphereIF sphereL(radius,centerL,inside);

  RealVect normal1rect1L(D_DECL(1.0,0.0,0.0));
  RealVect point1rect1L(D_DECL(-10.0,0.0,0.0));
  PlaneIF line1rect1L(normal1rect1L,point1rect1L,inside);

  RealVect normal2rect1L(D_DECL(-1.0,0.0,0.0));
  RealVect point2rect1L(D_DECL(-8.8,0.0,0.0));
  PlaneIF line2rect1L(normal2rect1L,point2rect1L,inside);

  RealVect normal3rect1L(D_DECL(0.0,1.0,0.0));
  RealVect point3rect1L(D_DECL(0.0,9.0,0.0));
  PlaneIF line3rect1L(normal3rect1L,point3rect1L,inside);

  RealVect normal4rect1L(D_DECL(0.0,-1.0,0.0));
  RealVect point4rect1L(D_DECL(0.0,14.0,0.0));
  PlaneIF line4rect1L(normal4rect1L,point4rect1L,inside);

  Vector<BaseIF*> linesrect1L(4);
  linesrect1L[0] = &line1rect1L;
  linesrect1L[1] = &line2rect1L;
  linesrect1L[2] = &line3rect1L;
  linesrect1L[3] = &line4rect1L;

  IntersectionIF rect1L(linesrect1L);

  RealVect normal1rect2L(D_DECL(1.0,0.0,0.0));
  RealVect point1rect2L(D_DECL(-10.0,0.0,0.0));
  PlaneIF line1rect2L(normal1rect2L,point1rect2L,inside);

  RealVect normal2rect2L(D_DECL(-1.0,0.0,0.0));
  RealVect point2rect2L(D_DECL(-5.0,0.0,0.0));
  PlaneIF line2rect2L(normal2rect2L,point2rect2L,inside);

  RealVect normal3rect2L(D_DECL(0.0,1.0,0.0));
  RealVect point3rect2L(D_DECL(0.0,13.4,0.0));
  PlaneIF line3rect2L(normal3rect2L,point3rect2L,inside);

  RealVect normal4rect2L(D_DECL(0.0,-1.0,0.0));
  RealVect point4rect2L(D_DECL(0.0,14.6,0.0));
  PlaneIF line4rect2L(normal4rect2L,point4rect2L,inside);

  Vector<BaseIF*> linesrect2L(4);
  linesrect2L[0] = &line1rect2L;
  linesrect2L[1] = &line2rect2L;
  linesrect2L[2] = &line3rect2L;
  linesrect2L[3] = &line4rect2L;

  IntersectionIF rect2L(linesrect2L);

  RealVect normal1rect3L(D_DECL(1.0,0.0,0.0));
  RealVect point1rect3L(D_DECL(-5.0,0.0,0.0));
  PlaneIF line1rect3L(normal1rect3L,point1rect3L,inside);

  RealVect normal2rect3L(D_DECL(-1.0,0.0,0.0));
  RealVect point2rect3L(D_DECL(-2.5,0.0,0.0));
  PlaneIF line2rect3L(normal2rect3L,point2rect3L,inside);

  RealVect normal3rect3L(D_DECL(0.0,1.0,0.0));
  RealVect point3rect3L(D_DECL(0.0,11.0,0.0));
  PlaneIF line3rect3L(normal3rect3L,point3rect3L,inside);

  RealVect normal4rect3L(D_DECL(0.0,-1.0,0.0));
  RealVect point4rect3L(D_DECL(0.0,17.0,0.0));
  PlaneIF line4rect3L(normal4rect3L,point4rect3L,inside);

  Vector<BaseIF*> linesrect3L(4);
  linesrect3L[0] = &line1rect3L;
  linesrect3L[1] = &line2rect3L;
  linesrect3L[2] = &line3rect3L;
  linesrect3L[3] = &line4rect3L;

  IntersectionIF rect3L(linesrect3L);

  Vector<BaseIF*> SRLparts(4);
  SRLparts[0] = &sphereL;
  SRLparts[1] = &rect1L;
  SRLparts[2] = &rect2L;
  SRLparts[3] = &rect3L;

  UnionIF SRL(SRLparts);

  RealVect centerR(D_DECL( 9.40,6.65,0.00));
  SphereIF sphereR(radius,centerR,inside);

  RealVect normal1rect1R(D_DECL(-1.0,0.0,0.0));
  RealVect point1rect1R(D_DECL(10.0,0.0,0.0));
  PlaneIF line1rect1R(normal1rect1R,point1rect1R,inside);

  RealVect normal2rect1R(D_DECL(1.0,0.0,0.0));
  RealVect point2rect1R(D_DECL(8.8,0.0,0.0));
  PlaneIF line2rect1R(normal2rect1R,point2rect1R,inside);

  RealVect normal3rect1R(D_DECL(0.0,1.0,0.0));
  RealVect point3rect1R(D_DECL(0.0,9.0,0.0));
  PlaneIF line3rect1R(normal3rect1R,point3rect1R,inside);

  RealVect normal4rect1R(D_DECL(0.0,-1.0,0.0));
  RealVect point4rect1R(D_DECL(0.0,14.0,0.0));
  PlaneIF line4rect1R(normal4rect1R,point4rect1R,inside);

  Vector<BaseIF*> linesrect1R(4);
  linesrect1R[0] = &line1rect1R;
  linesrect1R[1] = &line2rect1R;
  linesrect1R[2] = &line3rect1R;
  linesrect1R[3] = &line4rect1R;

  IntersectionIF rect1R(linesrect1R);

  RealVect normal1rect2R(D_DECL(-1.0,0.0,0.0));
  RealVect point1rect2R(D_DECL(10.0,0.0,0.0));
  PlaneIF line1rect2R(normal1rect2R,point1rect2R,inside);

  RealVect normal2rect2R(D_DECL(1.0,0.0,0.0));
  RealVect point2rect2R(D_DECL(5.0,0.0,0.0));
  PlaneIF line2rect2R(normal2rect2R,point2rect2R,inside);

  RealVect normal3rect2R(D_DECL(0.0,1.0,0.0));
  RealVect point3rect2R(D_DECL(0.0,13.4,0.0));
  PlaneIF line3rect2R(normal3rect2R,point3rect2R,inside);

  RealVect normal4rect2R(D_DECL(0.0,-1.0,0.0));
  RealVect point4rect2R(D_DECL(0.0,14.6,0.0));
  PlaneIF line4rect2R(normal4rect2R,point4rect2R,inside);

  Vector<BaseIF*> linesrect2R(4);
  linesrect2R[0] = &line1rect2R;
  linesrect2R[1] = &line2rect2R;
  linesrect2R[2] = &line3rect2R;
  linesrect2R[3] = &line4rect2R;

  IntersectionIF rect2R(linesrect2R);

  RealVect normal1rect3R(D_DECL(-1.0,0.0,0.0));
  RealVect point1rect3R(D_DECL(5.0,0.0,0.0));
  PlaneIF line1rect3R(normal1rect3R,point1rect3R,inside);

  RealVect normal2rect3R(D_DECL(1.0,0.0,0.0));
  RealVect point2rect3R(D_DECL(2.5,0.0,0.0));
  PlaneIF line2rect3R(normal2rect3R,point2rect3R,inside);

  RealVect normal3rect3R(D_DECL(0.0,1.0,0.0));
  RealVect point3rect3R(D_DECL(0.0,11.0,0.0));
  PlaneIF line3rect3R(normal3rect3R,point3rect3R,inside);

  RealVect normal4rect3R(D_DECL(0.0,-1.0,0.0));
  RealVect point4rect3R(D_DECL(0.0,17.0,0.0));
  PlaneIF line4rect3R(normal4rect3R,point4rect3R,inside);

  Vector<BaseIF*> linesrect3R(4);
  linesrect3R[0] = &line1rect3R;
  linesrect3R[1] = &line2rect3R;
  linesrect3R[2] = &line3rect3R;
  linesrect3R[3] = &line4rect3R;

  IntersectionIF rect3R(linesrect3R);

  Vector<BaseIF*> SRRparts(4);
  SRRparts[0] = &sphereR;
  SRRparts[1] = &rect1R;
  SRRparts[2] = &rect2R;
  SRRparts[3] = &rect3R;

  UnionIF SRR(SRRparts);

  UnionIF SR(SRL,SRR);

  RealVect normal1XTRA(D_DECL(-1.0,0.0,0.0));
  RealVect point1XTRA(D_DECL(-21.0,0.0,0.0));
  PlaneIF line1XTRA(normal1XTRA,point1XTRA,inside);

  RealVect normal2XTRA(D_DECL(1.0,0.0,0.0));
  RealVect point2XTRA(D_DECL(21.0,0.0,0.0));
  PlaneIF line2XTRA(normal2XTRA,point2XTRA,inside);

  RealVect normal3XTRA(D_DECL(0.0,-1.0,0.0));
  RealVect point3XTRA(D_DECL(0.0,0.5,0.0));
  PlaneIF line3XTRA(normal3XTRA,point3XTRA,inside);

  RealVect normal4XTRA(D_DECL(0.0,1.0,0.0));
  RealVect point4XTRA(D_DECL(0.0,19.5,0.0));
  PlaneIF line4XTRA(normal4XTRA,point4XTRA,inside);

  RealVect normal1rect1XTRA(D_DECL(1.0,0.0,0.0));
  RealVect point1rect1XTRA(D_DECL(-1.5,0.0,0.0));
  PlaneIF line1rect1XTRA(normal1rect1XTRA,point1rect1XTRA,inside);

  RealVect normal2rect1XTRA(D_DECL(-1.0,0.0,0.0));
  RealVect point2rect1XTRA(D_DECL(1.5,0.0,0.0));
  PlaneIF line2rect1XTRA(normal2rect1XTRA,point2rect1XTRA,inside);

  RealVect normal3rect1XTRA(D_DECL(0.0,1.0,0.0));
  RealVect point3rect1XTRA(D_DECL(0.0,8.0,0.0));
  PlaneIF line3rect1XTRA(normal3rect1XTRA,point3rect1XTRA,inside);

  RealVect normal4rect1XTRA(D_DECL(0.0,-1.0,0.0));
  RealVect point4rect1XTRA(D_DECL(0.0,21.5,0.0));
  PlaneIF line4rect1XTRA(normal4rect1XTRA,point4rect1XTRA,inside);

  Vector<BaseIF*> linesrect1XTRA(4);
  linesrect1XTRA[0] = &line1rect1XTRA;
  linesrect1XTRA[1] = &line2rect1XTRA;
  linesrect1XTRA[2] = &line3rect1XTRA;
  linesrect1XTRA[3] = &line4rect1XTRA;

  IntersectionIF rect1XTRA(linesrect1XTRA);

  Vector<BaseIF*> XTRAparts(5);
  XTRAparts[0] = &line1XTRA;
  XTRAparts[1] = &line2XTRA;
  XTRAparts[2] = &line3XTRA;
  XTRAparts[3] = &line4XTRA;
  XTRAparts[4] = &rect1XTRA;

  UnionIF XTRA(XTRAparts);

  RefCountedPtr<BaseIF> everythingPhase0( new UnionIF(SR,XTRA));

  // Complement for MF
  bool complement = true;
  RefCountedPtr<BaseIF> everythingPhase1(new ComplementIF(*everythingPhase0,complement));

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  Vector<GeometryService*> geometries(2,NULL);
  bool useNewGeometry = true;
  pp.query("use_new_geometry_gen", useNewGeometry);
  if(useNewGeometry)
    {
      Real threshold = sqrt(2.);
      pp.query("ba_threshold", threshold);
      pout() << "using boundaryarea wrapped gshop refinement criterion with threshold ba frac=" <<  threshold << endl;
      RefCountedPtr<WGSRefinementCriterion> refcrit(new BoundaryAreaRefCrit(threshold));
      pout() << "using new wrapped gshop geometry generation" << endl;
      int minRefine = 0;
      int maxRefine = 3;
      WrappedGShop* workshop0 = new WrappedGShop(everythingPhase0,RealVect::Zero, a_dx, a_domain, minRefine, maxRefine);
      workshop0->setRefinementCriterion(refcrit);
      workshop0->m_phase = 0;
      geometries[0] = workshop0;

      WrappedGShop* workshop1 = new WrappedGShop(everythingPhase1,RealVect::Zero, a_dx, a_domain, minRefine, maxRefine);
      workshop0->setRefinementCriterion(refcrit);
      workshop1->m_phase = 1;
      geometries[1] = workshop1;
    }
  else
    {
      pout() << "using old style geometry generation" << endl;
      GeometryShop* workshop0 = new GeometryShop(*everythingPhase0,0,vectDx);
      workshop0->m_phase = 0;
      geometries[0] = workshop0;

      GeometryShop* workshop1 = new GeometryShop(*everythingPhase1,0,vectDx);
      workshop1->m_phase = 1;
      geometries[1] = workshop1;
    }

  int maxBoxSize;
  pp.get("maxboxsize",maxBoxSize);

  // int maxCoarsenings = 0;

  // This generates the new EBIS
  a_mfIndexSpace.define(a_domain,a_origin,a_dx,geometries,maxBoxSize/* ,maxCoarsenings */);

  BaseIF* retval = everythingPhase0->newImplicitFunction();
  delete geometries[0];
  delete geometries[1];
  return retval;
}

BaseIF* makeMitochondriaGeometry(MFIndexSpace   & a_mfIndexSpace,
                                 const Box      & a_domain,
                                 const RealVect & a_origin,
                                 const Real     & a_dx)
{
  // parse input file
  ParmParse pp;

  RefCountedPtr<BaseIF> everythingPhase0( new MitochondriaIF1());

  // Complement for MF
  bool complement = true;
  RefCountedPtr<BaseIF> everythingPhase1(new ComplementIF(*everythingPhase0,complement));

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  Vector<GeometryService*> geometries(2,NULL);

  bool useNewGeometry = true;
  pp.query("use_new_geometry_gen", useNewGeometry);
  
  if(useNewGeometry)
    {
      Real threshold = sqrt(2.);
      pp.query("ba_threshold", threshold);
      pout() << "using boundaryarea wrapped gshop refinement criterion with threshold ba frac=" <<  threshold << endl;
      RefCountedPtr<WGSRefinementCriterion> refcrit(new BoundaryAreaRefCrit(threshold));
      pout() << "using new wrapped gshop geometry generation" << endl;
      int minRefine = 0;
      int maxRefine = 3;
      WrappedGShop* workshop0 = new WrappedGShop(everythingPhase0,RealVect::Zero, a_dx, a_domain, minRefine, maxRefine);
      workshop0->setRefinementCriterion(refcrit);
      workshop0->m_phase = 0;
      geometries[0] = workshop0;

      WrappedGShop* workshop1 = new WrappedGShop(everythingPhase1,RealVect::Zero, a_dx, a_domain, minRefine, maxRefine);
      workshop0->setRefinementCriterion(refcrit);
      workshop1->m_phase = 1;
      geometries[1] = workshop1;
    }
  else
    {
      pout() << "using old style geometry generation" << endl;
      GeometryShop* workshop0 = new GeometryShop(*everythingPhase0,0,vectDx);
      workshop0->m_phase = 0;
      geometries[0] = workshop0;

      GeometryShop* workshop1 = new GeometryShop(*everythingPhase1,0,vectDx);
      workshop1->m_phase = 1;
      geometries[1] = workshop1;
    }
  int maxBoxSize;
  pp.get("maxboxsize",maxBoxSize);

  // int maxCoarsenings = 0;

  // This generates the new EBIS
  a_mfIndexSpace.define(a_domain,a_origin,a_dx,geometries,maxBoxSize/* ,maxCoarsenings */);

  BaseIF* retval =  everythingPhase0->newImplicitFunction();
  delete geometries[0];
  delete geometries[1];
  return retval;
}

void getConnectedComponents(Vector<RefCountedPtr<EBIndexSpace> > & a_allComponents,
                            MFIndexSpace                         & a_mfIndexSpace)
{
  for (int phase = 0; phase < 2; phase++)
  {
    // Select one index-space
    Chombo_EBIS::alias((EBIndexSpace*)a_mfIndexSpace.EBIS(phase));

    EBIndexSpace* ebisPtr = Chombo_EBIS::instance();

    Vector<RefCountedPtr<EBIndexSpace> > components;
    components = ebisPtr->connectedComponents();
    
    pout() << "phase = " << phase << ", ncomp = " << components.size() << endl;

    for (int comp = 0; comp < components.size(); comp++)
    {
      a_allComponents.push_back(components[comp]);
    }
  }
}

void makeLayout(DisjointBoxLayout& a_dbl,
                const Box&         a_domain)
{
  ParmParse pp;

  Vector<Box> vbox(1,a_domain);

  int maxBoxSize;
  pp.get("maxboxsize",maxBoxSize);

  domainSplit(a_domain,vbox,maxBoxSize);

  Vector<int> procAssign;
  LoadBalance(procAssign,vbox);

  a_dbl.define(vbox,procAssign);
}

void makeEBISL(EBISLayout&              a_ebisl,
               const DisjointBoxLayout& a_grids,
               const Box&               a_domain,
               const int&               a_nghost)
{
  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  CH_assert(ebisPtr->isDefined());

  ebisPtr->fillEBISLayout(a_ebisl,a_grids,a_domain,a_nghost);
}

void fillData(LevelData<EBCellFAB>& a_level,
              const RealVect&       a_origin,
              const Real&           a_dx,
              const BaseIF&         a_implicit)
{
  // Set the data in each VoF to the value of the implicit function at the
  // cell center of the cell containing the VoF
  for (DataIterator dit = a_level.dataIterator(); dit.ok(); ++dit)
  {
    EBCellFAB& ebcell = a_level[dit()];
    const EBISBox& ebisbox = ebcell.getEBISBox();
    FArrayBox& data = ebcell.getFArrayBox();
    const Box& box = ebcell.box();

    for (BoxIterator bit(box); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();

      RealVect coord(a_origin);
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        coord[idir] += a_dx * (iv[idir] + 0.5);
      }

      Real value = a_implicit.value(coord);

      Vector<VolIndex> vofs = ebisbox.getVoFs(iv);
      int size = vofs.size();

      if (size == 0)
      {
        data(iv,0) = value;
      }
      else
      {
        for (int i = 0; i < size; i++)
        {
          ebcell(vofs[i],0) = value;
        }
      }
    }
  }
}

void fillData(LevelData<EBCellFAB>& a_level,
              const RealVect&       a_origin,
              const Real&           a_dx,
              const Real&           a_value)
{
  // Set the data in each VoF to the value of "a_value
  for (DataIterator dit = a_level.dataIterator(); dit.ok(); ++dit)
  {
    EBCellFAB& ebcell = a_level[dit()];
    const EBISBox& ebisbox = ebcell.getEBISBox();
    FArrayBox& data = ebcell.getFArrayBox();
    const Box& box = ebcell.box();

    for (BoxIterator bit(box); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();

      RealVect coord(a_origin);
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        coord[idir] += a_dx * (iv[idir] + 0.5);
      }

      Vector<VolIndex> vofs = ebisbox.getVoFs(iv);
      int size = vofs.size();

      if (size == 0)
      {
        data(iv,0) = a_value;
      }
      else
      {
        for (int i = 0; i < size; i++)
        {
          ebcell(vofs[i],0) = a_value;
        }
      }
    }
  }
}

void fillData(LevelData<EBCellFAB>& a_level,
              const RealVect&       a_origin,
              const RealVect&       a_dx,
              const BaseIF&         a_implicit)
{
  // Set the data in each VoF to the value of the implicit function at the
  // cell center of the cell containing the VoF
  for (DataIterator dit = a_level.dataIterator(); dit.ok(); ++dit)
  {
    EBCellFAB& ebcell = a_level[dit()];
    const EBISBox& ebisbox = ebcell.getEBISBox();
    FArrayBox& data = ebcell.getFArrayBox();
    const Box& box = ebcell.box();

    for (BoxIterator bit(box); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();

      RealVect coord(a_origin);
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        coord[idir] += a_dx[idir] * (iv[idir] + 0.5);
      }

      Real value = a_implicit.value(coord);

      Vector<VolIndex> vofs = ebisbox.getVoFs(iv);
      int size = vofs.size();

      if (size == 0)
      {
        data(iv,0) = value;
      }
      else
      {
        for (int i = 0; i < size; i++)
        {
          ebcell(vofs[i],0) = value;
        }
      }
    }
  }
}

void fillData(LevelData<EBCellFAB>& a_data,
              const RealVect&       a_origin,
              const RealVect&       a_dx,
              const BaseIF&         a_implicit,
              const IntVect&        a_ghost)
{
  // Set the data in each VoF to the value of the implicit function at the
  // cell center of the cell containing the VoF
  const DisjointBoxLayout& dbl = a_data.disjointBoxLayout();
  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit)
  {
    const Box& dblBox = dbl.get(dit());
    EBCellFAB& ebcell = a_data[dit()];
    const EBISBox& ebisbox = ebcell.getEBISBox();
    FArrayBox& data = ebcell.getFArrayBox();
    Box ghostedBox = dblBox;
    ghostedBox.grow(a_ghost);

    for (BoxIterator bit(ghostedBox); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();

      RealVect coord(a_origin);
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        coord[idir] += a_dx[idir] * (iv[idir] + 0.5);
      }

      Real value = a_implicit.value(coord);

      if(!dblBox.contains(iv))
      {//do this for domain ghosts (interior ghosts get over-written by an exchange)
        data(iv,0) = value;
      }
      else
      {
        Vector<VolIndex> vofs = ebisbox.getVoFs(iv);
        int size = vofs.size();
        if (size == 0)
        {
          data(iv,0) = value;
        }
        else
        {
          for (int i = 0; i < size; i++)
          {
            ebcell(vofs[i],0) = value;
          }
        }
      }
    }
  }
  a_data.exchange();
}

void makeLayouts(Vector<DisjointBoxLayout>& a_dbl,
                 const Vector<int>&         a_refRatio,
                 const Vector<RealVect>&    a_dx,
                 const RealVect&            a_origin,
                 const Vector<Box>&         a_domain)
{
  pout() << "Begin making layouts " << endl;
  ParmParse pp;

  int maxBoxSize;
  pp.get("maxboxsize",maxBoxSize);

  Real fillRatio;
  pp.get("fill_ratio",fillRatio);

  int blockFactor;
  pp.get("block_factor",blockFactor);

  int bufferSize;
  pp.get("buffer_size",bufferSize);

  int numLevels = a_domain.size();
  a_dbl.resize(numLevels);

//   for(int ilev = 0;ilev<numLevels;ilev++)
//   {
//     Vector<Box> vbox(1,a_domain[ilev]);
//
//     domainSplit(a_domain[ilev],vbox,maxBoxSize);
//
//     Vector<int> procAssign;
//     LoadBalance(procAssign,vbox);
//
//     a_dbl[ilev].define(vbox,procAssign);
//   }

  ///////

  // first, construct tags
  Vector<IntVectSet> tagsVect(numLevels);
  Vector<Vector<Box> > boxes;
  Vector<Vector<Box> > vectBoxes(numLevels);

  for(int ilev = 0; ilev < numLevels; ilev++)
  {
    domainSplit(a_domain[ilev], vectBoxes[ilev],
                maxBoxSize);
  }

  if(numLevels != 1)
  {
    // do tagging
    tagCells(tagsVect,a_refRatio,a_dx,a_origin,a_domain);

    BRMeshRefine meshrefine(a_domain[0], a_refRatio,
                            fillRatio,   blockFactor,
                            bufferSize,  maxBoxSize);

    int top_level = numLevels - 1;

    // for now, just assume that all levels get regridded at the same time
    int lbase = 0;
    int new_finest_level = meshrefine.regrid(boxes,
                                             tagsVect,
                                             lbase,
                                             top_level,
                                             vectBoxes);

    CH_assert(new_finest_level == numLevels-1);

    for(int ilev = 0; ilev < numLevels; ilev++)
    {
      Vector<int> procAssign;
      LoadBalance(procAssign,boxes[ilev]);

      a_dbl[ilev].define(boxes[ilev],procAssign);
    }
  }
  else
  {
    Vector<int> procAssign;
    LoadBalance(procAssign,vectBoxes[0]);

    a_dbl[0].define(vectBoxes[0],procAssign);
  }
}

void tagCells(Vector<IntVectSet>&     a_tags,
              const Vector<int>&      a_refRatio,
              const Vector<RealVect>& a_dx,
              const RealVect&         a_origin,
              const Vector<Box>&      a_domain)
{
  int numLevels = a_tags.size();

  for (int lev=0; lev<numLevels-1; lev++)//only tag on coarser levels
  {
    IntVectSet& levelTags = a_tags[lev];

    tagCellsLevel(levelTags, a_domain[lev],a_dx[lev],a_origin);

    //fix the no-tag case (refine based on coarser)
    const int numPts = levelTags.numPts();
    if(numPts == 0)
    {
      if(lev==0)
      {//tag based on an arbitrary box
        Box tagBox = a_domain[lev];
        tagBox.coarsen(a_refRatio[lev]);
        levelTags |= tagBox;
        MayDay::Warning("no tags for coarsest level::arbitrary refinement specified");
      }
      else
      {//tag based on coarser tags
        MayDay::Warning("tagging level based on coarser tags");
        IntVectSet refTags = a_tags[lev-1];
        refTags.refine(a_refRatio[lev]);
        levelTags |= refTags;
      }
    }
  }
}

void tagCellsLevel(IntVectSet&     a_tags,
                   const Box&      a_domain,
                   const RealVect& a_dx,
                   const RealVect& a_origin)
{
  // parse input file
  ParmParse pp;

  //tag boxes based on low and high points in physical coordinates specified in the input file
  int numBox = 2;
  Vector<const char*> nameLo(numBox);
  Vector<const char*> nameHi(numBox);
  nameLo[0] = "tag_box1_lo";
  nameHi[0] = "tag_box1_hi";
  nameLo[1] = "tag_box2_lo";
  nameHi[1] = "tag_box2_hi";

  for (int iBox = 0; iBox < numBox; iBox++)
  {
    Vector<Real> lo(SpaceDim);
    Vector<Real> hi(SpaceDim);

    pp.getarr(nameLo[iBox],lo,0,SpaceDim);
    pp.getarr(nameHi[iBox],hi,0,SpaceDim);

    //figure out the lo and hi iv's of the tag box
    IntVect ivLo,ivHi;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      Real xivLo = (lo[idir] - a_origin[idir])/a_dx[idir];
      Real xivHi = (hi[idir] - a_origin[idir])/a_dx[idir];
      ivLo[idir] = int(xivLo - 0.5);
      ivHi[idir] = int(xivHi - 0.5);
    }

    //make the box and add it to a_tags
    Box tagBox = Box(ivLo,ivHi);
    tagBox &= a_domain;
    a_tags |= tagBox;
  }
}

void makeEBISL(Vector<EBISLayout>&              a_ebisl,
               const Vector<DisjointBoxLayout>& a_grids,
               const Vector<Box>&               a_domain,
               const int&                       a_nghost)
{
  pout() << "Begin making EBISLs " << endl;

  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  CH_assert(ebisPtr->isDefined());

  int numLevels = a_domain.size();
  a_ebisl.resize(numLevels);
  for(int ilev = 0;ilev<numLevels;ilev++)
  {
    ebisPtr->fillEBISLayout(a_ebisl[ilev],a_grids[ilev],a_domain[ilev],a_nghost);
  }
}


void createEBDistributionFiles()
{
#ifdef CH_USE_HDF5
  pout() << "Begin createEBDistributionFiles " << endl;

  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  CH_assert(ebisPtr->isDefined());

  ProblemDomain domain;
  Vector<LevelData<EBCellFAB>*> allBoxes, IrregBoxes;
  int numLevels = ebisPtr->numLevels();

  domain = ebisPtr->getBox(numLevels-1);
  Real dx = ebisPtr->dx(numLevels-1);
  Vector<int> refRatio(numLevels, 2);
  Vector<DisjointBoxLayout> grids;
  Vector<std::string> names(1, std::string("procID"));

  for(int i=numLevels-1; i>=0; i--)
  {
    DisjointBoxLayout dbl = ebisPtr->levelGrids(i);
    grids.push_back(dbl);
    EBISLayout ebisl;
    ebisPtr->fillEBISLayout(ebisl, dbl, ebisPtr->getBox(i), 1);
    EBCellFactory factory(ebisl);
    allBoxes.push_back(new LevelData<EBCellFAB>(dbl, 1, IntVect::Zero, factory));
    LevelData<EBCellFAB>& ld = *(allBoxes.back());
    for(DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      int p = procID();
      ld[dit].setVal(p);
    }
  }

  writeEBHDF5(std::string("EBIndexSpace.hdf5"),
              grids,
              allBoxes,
              names,
              domain.domainBox(),
              dx,
              0,
              0,
              refRatio,
              numLevels,
              false,
              Vector<Real>(1,0));

  for(int i=0; i<numLevels; ++i)
  {
    delete allBoxes[i];
  }
#endif
}
