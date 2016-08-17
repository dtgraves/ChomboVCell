#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _EBMENAGERIEUTILS_H_
#define _EBMENAGERIEUTILS_H_

#include "Box.H"
#include "RealVect.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "MFIndexSpace.H"
#include "GeometryShop.H"

#include "BaseIF.H"

#include "UsingNamespace.H"

///
BaseIF* makeMultiSphereGeometry(MFIndexSpace   & a_mfIndexSpace,
                                const Box      & a_domain,
                                const RealVect & a_origin,
                                const Real     & a_dx);

BaseIF* makeCalciumOneGeometry(MFIndexSpace   & a_mfIndexSpace,
                               const Box      & a_domain,
                               const RealVect & a_origin,
                               const Real     & a_dx);

BaseIF* makeMitochondriaGeometry(MFIndexSpace   & a_mfIndexSpace,
                                 const Box      & a_domain,
                                 const RealVect & a_origin,
                                 const Real     & a_dx);


void getDomainOriginDx(Box      & a_domain,
                       RealVect & a_origin,
                       Real     & a_dx);

///
/**
   This makes each connected component into its own EBIS (assumes that the
   geometry has already been constructed as a MFIndexSpace.
 */
void getConnectedComponents(Vector<RefCountedPtr<EBIndexSpace> > & a_allComponents,
                            MFIndexSpace                         & a_mfIndexSpace);

/***************/
/// Define grids by splitting up domain
/***************/
void makeLayout(DisjointBoxLayout & a_dbl1,
                const Box         & a_domain);

/***************/
// Make the corresponding layout
/***************/
void makeEBISL(EBISLayout              & a_ebisl,
               const DisjointBoxLayout & a_grids,
               const Box               & a_domain,
               const int               & a_nghost);

/***************/
// Put data in the layout
/***************/
void fillData(LevelData<EBCellFAB> & a_level,
              const RealVect       & a_origin,
              const Real           & a_dx,
              const BaseIF         & a_implicit);

/***************/
// Put constant data in the layout
/***************/
void fillData(LevelData<EBCellFAB> & a_level,
              const RealVect       & a_origin,
              const Real           & a_dx,
              const Real           & a_value);

/***************/
// Put data in the layout
/***************/
void fillData(LevelData<EBCellFAB> & a_level,
              const RealVect       & a_origin,
              const RealVect       & a_dx,
              const BaseIF         & a_implicit);

/***************/
/// Define grids by splitting up domain
/***************/
void makeLayouts(Vector<DisjointBoxLayout> & a_dbl,
                 const Vector<int>         & a_refRatio,
                 const Vector<RealVect>    & a_dx,
                 const RealVect            & a_origin,
                 const Vector<Box>         & a_domain);

/***************/
/// Tag cells for Local Refinement
/***************/
void tagCells(Vector<IntVectSet>     & a_tags,
              const Vector<int>      & a_refRatio,
              const Vector<RealVect> & a_dx,
              const RealVect         & a_origin,
              const Vector<Box>      & a_domain);

void tagCellsLevel(IntVectSet     & a_tags,
                   const Box      & a_domain,
                   const RealVect & a_dx,
                   const RealVect & a_origin);

/***************/
// Make the corresponding layout
/***************/
void makeEBISL(Vector<EBISLayout>              & a_ebisl,
               const Vector<DisjointBoxLayout> & a_grids,
               const Vector<Box>               & a_domain,
               const int                       & a_nghost);

/***************/
// Put data in the layout
/***************/
void fillData(LevelData<EBCellFAB> & a_data,
              const RealVect       & a_origin,
              const RealVect       & a_dx,
              const BaseIF         & a_implicit,
              const IntVect        & a_ghost);

void createEBDistributionFiles();

#endif
