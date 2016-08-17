#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>
#include <cstdio>
#include <iostream>

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "LevelData.H"
#include "DebugDump.H"
#include "FABView.H"

#include "PolyGeom.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "EBArith.H"
#include "RedistStencil.H"
#include "VoFIterator.H"
#include "EBLevelRedist.H"
#include "EBAMRDataOps.H"
#include "EBAMRIO.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBAMRPoissonOp.H"
#include "EBBackwardEuler.H"
#include "EBDebugDump.H"
#include "EBFABView.H"

#include "GeometryShop.H"
#include "AllRegularService.H"
#include "TiltedCylinderIF.H"

#include "PoissonUtilities.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

/************/
/************/
void generateData(Vector< LevelData<EBCellFAB>* >& a_datum,
                  Vector< DisjointBoxLayout     >& a_grids,
                  Vector< EBISLayout            >& a_ebisl,
                  const Vector<Vector<Box>      >& a_boxes,
                  const PoissonParameters        & a_params,
                  const int                      & a_refToFinestCalc)
{
  ParmParse pp;

  defineGrids(a_grids, a_ebisl, a_boxes, a_params.refRatio, a_params);

  pout() << "filling data" << endl;
  Vector<LevelData<EBCellFAB>* > phiOld, source;
  Vector<LevelData<EBCellFAB>* >& phiNew = a_datum;

  Real domainValue;
  pp.get("domain_value",domainValue);

  newAndFillSmoothCellData(phiOld, a_grids, a_ebisl, a_params, 1, true, domainValue);
  newAndFillSmoothCellData(phiNew, a_grids, a_ebisl, a_params, 1, true, domainValue);
  newAndFillSmoothCellData(source, a_grids, a_ebisl, a_params, 1, true, 0.0);

  pout() << "defining  backward Euler  solver" << endl;
  RefCountedPtr< AMRMultiGrid<LevelData<EBCellFAB> > >
    solver = RefCountedPtr< AMRMultiGrid<LevelData<EBCellFAB> > > (new AMRMultiGrid<LevelData<EBCellFAB> >);

  Real beta;
  pp.get("diffusion_coeff", beta);

  RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > > factory;
  Real alpha = 1;
  BiCGStabSolver<LevelData<EBCellFAB> > bottomSolver;
  bottomSolver.m_verbosity = 0;
  defineSolver(*solver, a_grids, a_ebisl, bottomSolver, a_params, 0.0,alpha, beta, factory);
  EBBackwardEuler               integrator(solver, *factory, a_params.coarsestDomain, a_params.refRatio, -1, 3);

  int numLevels = a_grids.size();

  Real dx = a_params.coarsestDx[0];

  Real dt;
  pp.get("dt", dt);

  int nstop;
  pp.get("max_steps", nstop);
  nstop /= a_refToFinestCalc;

  bool zeroPhi;
  pp.get("zero_phi",zeroPhi);

  Real time = 0;
  int maxlev = a_datum.size() - 1;

  Vector<string> names;
  names.push_back("stuff");

  bool replaceCovered = false;
  Vector<Real> coveredValues;

  Real maxval = 0;
  Real minval = 0;

  solver->init(phiOld, source, maxlev, 0);
  int istep;
  char filename[1000];
  for (istep = 0; istep < nstop; istep++)
    {
      EBAMRDataOps::getMaxMin(maxval, minval, phiOld, 0, a_params.refRatio);
      pout() << endl;
      pout() << "----- step = " << istep << ", time = "  << time << ", max = " << maxval << ", min = " <<  minval << endl;
      for (DataIterator dit = phiOld[0]->dataIterator(); dit.ok(); ++dit) {
        pout() << "----- value (4,4,3) = " 
               << (*(phiOld[0]))[dit()](VolIndex(IntVect(D_DECL(4,4,3)),0),0)
               << endl;
      }
      pout() << endl;

      sprintf(filename,"diffusion.%04d.%1dd.hdf5",istep,SpaceDim);
      writeEBHDF5(filename,
                  a_grids,
                  phiOld,
                  names,
                  a_params.coarsestDomain,
                  dx,
                  dt,
                  time,
                  a_params.refRatio,
                  numLevels,
                  replaceCovered,
                  coveredValues);

      integrator.oneStep(phiNew, phiOld, source, dt, 0, maxlev, zeroPhi);
      time += dt;

      EBAMRDataOps::assign(phiOld,phiNew);
    }

  EBAMRDataOps::getMaxMin(maxval, minval, phiOld, 0, a_params.refRatio);
  pout() << endl;
  pout() << "----- step = " << istep << ", time = "  << time << ", max = " << maxval << ", min = " <<  minval << endl;
  for (DataIterator dit = phiOld[0]->dataIterator(); dit.ok(); ++dit) {
    pout() << "----- value (4,4,3) = " 
           << (*(phiOld[0]))[dit()](VolIndex(IntVect(D_DECL(4,4,3)),0),0)
           << endl;
  }
  pout() << endl;

  sprintf(filename,"diffusion.%04d.%1dd.hdf5",istep,SpaceDim);
  writeEBHDF5(filename,
              a_grids,
              phiOld,
              names,
              a_params.coarsestDomain,
              dx,
              dt,
              time,
              a_params.refRatio,
              numLevels,
              replaceCovered,
              coveredValues);

  //phinew = datum lives on.
  for (int ilev = 0; ilev <= maxlev; ilev++)
    {
      delete phiOld[ilev];
      delete source[ilev];
    }
}

/*****/
void solutionErrorTest(int testverbosity,
                       int fileout)
{
  if (testverbosity >= 1)
    pout() << "getting parameters" << endl;

  PoissonParameters param;
  getPoissonParameters(param);

  if (testverbosity >= 1)
    pout() << "generating geometry" << endl;

  definePoissonGeometry(param);

  Box domain = param.coarsestDomain.domainBox();

  if (testverbosity >= 1)
    pout() << "getting grids" << endl;

  Vector<DisjointBoxLayout> grids;
  Vector<EBISLayout>        ebisl;
  Vector<Vector<Box> >      boxes;

  getSimpleBoxes(boxes, param.refRatio, domain);

  Vector< LevelData<EBCellFAB> *> solution;

  int ref = 1;

  if (testverbosity >= 1)
    pout() << "generating fine solution" << endl;
  generateData(solution, grids, ebisl, boxes, param, ref);

  for (int ilev = 0; ilev < param.refRatio.size(); ilev++)
    {
      delete solution[ilev];
    }
}

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
  {
    // setChomboMPIErrorHandler();
#endif

    // Check for an input file
    char* inFile = NULL;

    if (a_argc > 1)
      {
        inFile = a_argv[1];
      }
    else
      {
        pout() << "Usage: <executable name> <inputfile>" << endl;
        pout() << "No input file specified" << endl;
        return -1;
      }

    ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);

    int testverbosity;
    pp.get("testverbosity", testverbosity);

    int fileout;
    pp.get("do_error_output", fileout);

    solutionErrorTest(testverbosity, fileout);

#ifdef CH_MPI
    dumpmemoryatexit();
  }
  MPI_Finalize();
#endif
}

