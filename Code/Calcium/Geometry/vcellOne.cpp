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
#include "DebugDump.H"

#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "PolyGeom.H"
#include "GeometryShop.H"
#include "MFIndexSpace.H"

#include "EBFABView.H"

#include "SphereIF.H"
#include "CH_Timer.H"
#include "UnionIF.H"
#include "IntersectionIF.H"
#include "ComplementIF.H"

#include "EBMenagerieUtils.H"

#include "UsingNamespace.H"

int main (int argc, char** argv)
{
#ifdef CH_MPI
  MPI_Init(&argc,&argv);
#endif

  // Begin forever present scoping trick
  {
    char* in_file = (char*)"vcellOne.inputs";

    if (argc > 1)
    {
      in_file = argv[1];
    }

    // Parse input file
    ParmParse pp(0,NULL,NULL,in_file);

    Box domain;
    RealVect origin;
    Real dx;

    BaseIF* implicit;
    MFIndexSpace mfIndexSpace;

    Vector<RefCountedPtr<EBIndexSpace> > allComponents;
    
    getDomainOriginDx( domain, origin, dx);
    implicit = makeCalciumOneGeometry(mfIndexSpace, domain, origin, dx);
    getConnectedComponents(allComponents, mfIndexSpace);

    for (int comp = 0; comp < allComponents.size(); comp++)
    {
      // Select one index-space
      Chombo_EBIS::alias(allComponents[comp]);

      // Make grids
      DisjointBoxLayout grids;
      makeLayout(grids,domain);

      // Create ebislayout
      int nghost = 0;
      EBISLayout ebisl;
      makeEBISL(ebisl,grids,domain,nghost);

      // Make a LevelData
      int nComps = 1;

      IntVect ghost = IntVect::Unit;
      ghost *= nghost;

      RefCountedPtr<DataFactory<EBCellFAB> > rcpFactory(new EBCellFactory(ebisl));
      LevelData<EBCellFAB> level(grids,nComps,ghost,*rcpFactory);

      // Put the component number in as data
      fillData(level,origin,dx,comp);

      // Write the data and the EB out
      const char* basename = "vcellOne";

      char name[1000];
      sprintf(name,"%s_%03d_%dd.hdf5",basename,comp,SpaceDim);
#ifdef CH_USE_HDF5
      writeEBLevelname(&level,name);
#endif
    }

    // Done with this object
    delete implicit;
  } // End scoping trick

  CH_TIMER_REPORT();

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return 0;
}
