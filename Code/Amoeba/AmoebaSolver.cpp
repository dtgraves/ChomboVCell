#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BRMeshRefine.H"
#include "RefCountedPtr.H"
#include "AMRMultiGrid.H"

#include "EBIndexSpace.H"
#include "MFIndexSpace.H"
#include "EBEllipticLoadBalance.H"
#include "EBLevelDataOps.H"
#include "EBAMRDataOps.H"
#include "EBAMRIO.H"
#include "EBLevelGrid.H"
#include "EBQuadCFInterp.H"
#include "EBAMRPoissonOpFactory.H"
#include "BaseDomainBC.H"
#include "DirichletConductivityDomainBC.H"
#include "DirichletConductivityEBBC.H"
#include "NeumannConductivityDomainBC.H"
#include "NeumannConductivityEBBC.H"
#include "DirichletPoissonDomainBC.H"
#include "DirichletPoissonEBBC.H"
#include "NeumannPoissonDomainBC.H"
#include "NeumannPoissonEBBC.H"
#include "EBMenagerieUtils.H"
#include "AmoebaSolver.H"
#include "AMRBoxesAndRanksIO.H"
#include "ParmParse.H"
#include "AmoebaSolverF_F.H"

AmoebaParams::
AmoebaParams()
{
  ParmParse pp;
  pp.get("const_coeff", m_constCoeff);

  pp.get("source_scaling",m_sourceScaling);
  pp.get("sink_scaling",m_sinkScaling);


  vector<Real> loCorner(SpaceDim);
  pp.getarr("prob_lo",loCorner,0,SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_loCorner[idir] = loCorner[idir];
    }


  pp.get("num_levels",m_numLevels);

  vector<int> numCells(SpaceDim);
  pp.getarr("num_cells",numCells,0,SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_numCells[idir] = numCells[idir];
    }
  

  Real prob_hi;
  pp.get("prob_hi", prob_hi);
  
  Real domainLength =prob_hi-m_loCorner[0];
  if(domainLength <= 0.0)
    {
      MayDay::Error("prob_hi bigger than prob_hi");
    }
  m_dx = domainLength/m_numCells[0];

  // Compute the indices of lo and hi end of the coarsest domain
  IntVect loEnd;
  IntVect hiEnd;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      loEnd[idir] = 0;
      hiEnd[idir] = m_numCells[idir] - 1;
    }

  m_refRatio.resize(m_numLevels);
  if (m_numLevels > 1)
    {
      vector<int> refRatio(m_numLevels-1);
      pp.getarr("ref_ratio",refRatio,0,m_numLevels-1);

      for (int ilev = 0; ilev < m_numLevels-1; ilev++)
        {
          m_refRatio[ilev] = refRatio[ilev];
        }
      m_refRatio[m_numLevels-1] = refRatio[m_numLevels-2];
    }
  else
    {
      m_refRatio[0] = 2;
    }

  // Make a box that size
  Box domainBox(loEnd,hiEnd);

  // Define a non-periodic problem domain
  m_coarsestDomain.define(domainBox);

  pp.get("dt",m_dt);
  pp.get("end_time",m_endTime);

  pp.get("output_interval",m_outputInterval);
  pp.get("output_prefix",m_outputPrefix);

  pp.get("maxboxsize",m_maxBoxSize);
  pp.get("block_factor",m_blockFactor);

  pp.get("fill_ratio",m_fillRatio);
  m_nestingRadius = 2;

  pp.get("mg_num_cycles",m_mgNumCycles);
  pp.get("mg_num_smooths",m_mgNumSmooths);
  pp.get("mg_relax_type",m_mgRelaxType);
  pp.get("mg_lazy_relax",m_mgLazyRelax);
  pp.get("mg_toler",m_mgToler);
  pp.get("mg_hang_toler",m_mgHangToler);
  pp.get("mg_iter_max",m_mgIterMax);
  pp.get("mg_num_precond_iter",m_mgNumPrecondIter);

  m_numGhostEBISLayout = 4;

  m_numGhostSoln   =3*IntVect::Unit;
  m_numGhostSource =  IntVect::Zero;

  pp.get("num_comp",m_ncomp);
}

// Print all the solver parameters
void 
AmoebaParams::
print()
{

  pout() << "Solver Parameters: \n";
  pout() << "source scaling = " << m_sourceScaling << "\n";
  pout() << "sink   scaling = " << m_sinkScaling   << "\n";
  pout() << "\n";
  pout() << "lo corner = " << m_loCorner << "\n";
  pout() << "\n";
  pout() << "dx = " << m_dx       << "\n";
  pout() << "\n";
  pout() << "num levels = " << m_numLevels << "\n";
  pout() << "num cells  = " << m_numCells  << "\n";
  if (m_numLevels > 1)
    {
      pout() << "ref ratio  = " << m_refRatio  << "\n";
    }
  pout() << "\n";
  pout() << "dt       = " << m_dt      << "\n";
  pout() << "end time = " << m_endTime << "\n";
  pout() << "\n";
  pout() << "output interval = " << m_outputInterval << "\n";
  pout() << "output prefix   = " << m_outputPrefix   << "\n";
  pout() << "\n";
  pout() << "max box size = " << m_maxBoxSize  << "\n";
  pout() << "block factor = " << m_blockFactor << "\n";
  pout() << "fill ratio = " << m_fillRatio << "\n";
  pout() << "mg num cycles       = " << m_mgNumCycles      << "\n";
  pout() << "mg num smooths      = " << m_mgNumSmooths     << "\n";
  pout() << "mg relax type       = " << m_mgRelaxType      << "\n";
  pout() << "mg lazy relax       = " << m_mgLazyRelax      << "\n";
  pout() << "mg toler            = " << m_mgToler          << "\n";
  pout() << "mg hang toler       = " << m_mgHangToler      << "\n";
  pout() << "mg iter max         = " << m_mgIterMax        << "\n";
  pout() << "mg num precond iter = " << m_mgNumPrecondIter << "\n";
  pout() << "\n";
}

AmoebaSolver::
AmoebaSolver(const AmoebaParams & a_params)
{
  // Save a copy of the solver parameters
  m_params = a_params;

  m_volumes.resize(0);
}

AmoebaSolver::
~AmoebaSolver()
{
  // Delete data holder for the solution at the old time (on all AMR levels)
  for (int ivol = 0; ivol < m_volumes.size(); ivol++)
    {
      for (int ilev = 0; ilev < m_solnOld[ivol].size(); ilev++)
        {
          delete m_solnOld[ivol][ilev];
          delete m_solnNew[ivol][ilev];
          delete m_scalOld[ivol][ilev];
          delete m_scalNew[ivol][ilev];
          delete m_scalRHS[ivol][ilev];
          delete m_soursin[ivol][ilev];
          delete m_extrStn[ivol][ilev];

          m_solnOld[ivol][ilev] = NULL;
          m_solnNew[ivol][ilev] = NULL;
          m_scalOld[ivol][ilev] = NULL;
          m_scalNew[ivol][ilev] = NULL;
          m_scalRHS[ivol][ilev] = NULL;
          m_soursin[ivol][ilev] = NULL;
          m_extrStn[ivol][ilev] = NULL;
        }
    }
  m_solnOld.resize(0);
  m_solnNew.resize(0);
  m_scalOld.resize(0);
  m_scalNew.resize(0);
  m_scalRHS.resize(0);
  m_soursin.resize(0);
  m_bounVal.resize(0);
  m_scalBou.resize(0);
  m_extrStn.resize(0);
 }

///
void 
AmoebaSolver::
getExtrapStencils(Vector<RefCountedPtr<BaseIndex  > >& a_destVoFs,
                  Vector<RefCountedPtr<BaseStencil> >& a_stencils,
                  const IntVectSet                   & a_cfivs,
                  const DataIndex                    & a_dit, 
                  int ivol, int ilev, Real a_dx)
{
  const EBISBox& ebisBox = m_ebisl[ivol][ilev][a_dit];
  const     Box&    grid = m_grids[ilev][a_dit];
  IntVectSet ivs = ebisBox.getIrregIVS(grid);
  VoFIterator vofit(ivs, ebisBox.getEBGraph());
  const Vector<VolIndex>& vofs = vofit.getVector();
  a_destVoFs.resize(vofs.size());
  a_stencils.resize(vofs.size());
  for(int ivof = 0; ivof < vofs.size(); ivof++)
    {
      VoFStencil  extrapStenc;
      //distance of extrapolation = dx*boundaryCentroid
      RealVect dist = ebisBox.bndryCentroid(vofs[ivof]);
      dist *= a_dx;
      
      Vector<VoFStencil> pointStenc;
      Vector<Real> distance;
      bool dropOrder = false;
      IntVectSet& cfivs = (IntVectSet&)a_cfivs;
      EBArith::johanStencil(dropOrder, pointStenc, distance, vofs[ivof], ebisBox,  a_dx*RealVect::Unit, cfivs);
      if(!dropOrder)
        {
          CH_assert(distance.size() >= 2);
          CH_assert(pointStenc.size() >= 2);
          //value at boundary is extrapolated along a ray.
          //we have data at points (phi) and distance to those points(x)
          //phi_boundary = (x1 phi0 - x0 phi1)/(x1-x0)
          pointStenc[0] *=  distance[1];
          pointStenc[1] *= -distance[0];
          extrapStenc += pointStenc[0];  //now val = x1 phi0
          extrapStenc += pointStenc[1];  //now val = x1 phi0 - x0 phi1
          extrapStenc *= 1.0/(distance[1]-distance[0]);  //now val = (x1 phi0 - x0 phi1)/(x1-x0)
        }
      else
        {
          EBArith::getExtrapolationStencil(extrapStenc, dist, a_dx*RealVect::Unit, vofs[ivof], ebisBox, -1, &cfivs, 0);
        }
      a_destVoFs[ivof] = RefCountedPtr<BaseIndex  >(new   VolIndex(vofs[ivof]));
      a_stencils[ivof] = RefCountedPtr<BaseStencil>(new VoFStencil(extrapStenc));
    }
}
///
void 
AmoebaSolver::
initStencils()
{
  m_extrStn.resize(m_volumes.size());
  for(int ivol = 0; ivol < m_volumes.size(); ivol++)
    {
      m_extrStn[ivol].resize(m_params.m_numLevels);
      Vector<EBLevelGrid>  eblg;
      Vector<RefCountedPtr<EBQuadCFInterp> > quadCFI;
      int inco = m_params.m_ncomp;
      getEBLGAndQuadCFI(eblg, quadCFI, ivol, inco);
      Real dxlev = m_params.m_dx;
      for(int ilev = 0; ilev < m_params.m_numLevels; ilev++)
        {
          m_extrStn[ivol][ilev] = new LayoutData< RefCountedPtr< AggStencil< EBCellFAB, BaseIVFAB<Real> > > >(m_grids[ilev]);
          for(DataIterator dit =m_grids[ilev].dataIterator(); dit.ok(); ++dit)
            {
              Vector< RefCountedPtr<BaseIndex>   > vofs;
              Vector< RefCountedPtr<BaseStencil> > stencils;
              
              getExtrapStencils(vofs, stencils, (*eblg[ilev].getCFIVS())[dit()], dit(), ivol, ilev,  dxlev);
              const       EBCellFAB& srcData = (*m_solnOld[ivol][ilev])[dit()];
              const BaseIVFAB<Real>& dstData = (*m_dataBou[ivol][ilev])[dit()];

              (*m_extrStn[ivol][ilev])[dit()] = RefCountedPtr< AggStencil < EBCellFAB, BaseIVFAB<Real> > >
                (new AggStencil < EBCellFAB, BaseIVFAB<Real> > (vofs, stencils, srcData, dstData));

            }
          if(ilev < m_params.m_numLevels-1)
            dxlev /= m_params.m_refRatio[ilev];
        }
    }
}
///
void 
AmoebaSolver::
defineSolver(RefCountedPtr<EBBackwardEuler>& a_integrator, int a_ivol, int a_ivar)
{
  CH_TIME("AmoebaSolver::defineSolver");

  // This is the multigrid solver used for backward Euler
  RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB > > >
    solver(new AMRMultiGrid<LevelData<EBCellFAB> >);


  // Set the verbosity of the bottom solver for multigrid
  m_bottomSolver.m_verbosity = 0;
  
  RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > > operatorFactory;
  if(!m_params.m_constCoeff)
    {
      RefCountedPtr<EBConductivityOpFactory> opfact;
  
      getVariableCoeffOpFactory(opfact, a_ivol, a_ivar);

      operatorFactory = opfact;
    }
  else
    {
      RefCountedPtr<EBAMRPoissonOpFactory> opfact;
  
      getConstantCoeffOpFactory(opfact, a_ivol, a_ivar);

      operatorFactory = opfact;
    }

  // Define the multigrid solver and set various parameters
  solver->define(m_params.m_coarsestDomain,
                 *operatorFactory,
                 &m_bottomSolver,
                 m_params.m_numLevels);

  Real normThresh = 1.0e-30;

  solver->setSolverParameters(m_params.m_mgNumSmooths,
                              m_params.m_mgNumSmooths,
                              m_params.m_mgNumSmooths,
                              m_params.m_mgNumCycles,
                              m_params.m_mgIterMax,
                              m_params.m_mgToler,
                              m_params.m_mgHangToler,
                              normThresh);

  solver->m_verbosity = 3;
  solver->init(m_scalOld[a_ivol],m_scalRHS[a_ivol],m_params.m_numLevels-1,0);

  // Create the backward Euler solver based on the multigrid solver
  a_integrator = RefCountedPtr<EBBackwardEuler>
    (new EBBackwardEuler(solver,*operatorFactory,
                         m_params.m_coarsestDomain,
                         m_params.m_refRatio,
                         m_params.m_numLevels) );
}
///
void 
AmoebaSolver::
init()
{
  CH_TIME("AmoebaSolver::init");

  // Initialize the geometry
  initGeometry();

  // Initialize the index space
  initIndexSpace();

  // Initialize the data
  initData();

  //initialize stencils for extrapolation
  initStencils();
}
///
void 
AmoebaSolver::
extrapolateDataToBoundary()
{
  CH_TIME("extrapolate_data_to_irreg_boundary");

  //Taylor series extrapolation of data to boundary
  //do this for all variables 
  int isrc = 0; int idst = 0; int inco = m_params.m_ncomp;
  for(int ivol = 0; ivol < m_volumes.size(); ivol++)
    {
      Vector<EBLevelGrid>  eblg;
      Vector<RefCountedPtr<EBQuadCFInterp> > quadCFI;
      getEBLGAndQuadCFI(eblg, quadCFI, ivol, inco);
      for(int ilev = 0; ilev < m_params.m_numLevels; ilev++)
        {
          if(ilev > 0)
            {
              quadCFI[ilev]->interpolate((*m_solnOld[ivol][ilev  ]), 
                                         (*m_solnOld[ivol][ilev-1]), 
                                         Interval(0, inco-1));

            }
          for(DataIterator dit =m_grids[ilev].dataIterator(); dit.ok(); ++dit)
            {
              //the last arguments say to loop over  all the variables with the stencil
              //and the false is not just increment the solution but replace the value there.
              (*m_extrStn[ivol][ilev])[dit()]->apply((*m_dataBou[ivol][ilev])[dit()], 
                                                     (*m_solnOld[ivol][ilev])[dit()], 
                                                     isrc, idst, inco, false);
            }
        }
    }
}
///
void 
AmoebaSolver::
setBoundaryValues()
{
  CH_TIME("set_boundary_values");

  //use data extrapolated from the boundary to set boundary data values.
  for(int ivol = 0; ivol < m_volumes.size(); ivol++)
    {
      Real dxlev = m_params.m_dx;
      for(int ilev = 0; ilev < m_params.m_numLevels; ilev++)
        {
          for(DataIterator dit =m_grids[ilev].dataIterator(); dit.ok(); ++dit)
            {
              const EBGraph& ebgraph =m_ebisl[ivol][ilev][dit()].getEBGraph();
              const IntVectSet& ivs = (*m_bounVal[ivol][ilev])[dit()].getIVS();
              for(VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit)
                {
                  for(int ivar = 0; ivar < m_params.m_ncomp; ivar++)
                    {
                      // if there is another volume to grab from, grab away and average what I 
                      // get.   as an example I am setting the boundary flux to (otherval - myval)/dx
                      // as a proxy to a gradient.   if there is no other value, i will set the value of the flux
                      //to zero but I will throw a warning (because I don't think it ought to happen)
                      Real thisVal = (*m_dataBou[ivol][ilev])[dit()](vofit(), ivar);
                      Real otherVals = 0;   
                      int nvals = 0;
                      bool found = false;
                      for(int jvol = 0; jvol < m_volumes.size(); jvol++)
                        {
                          if(jvol != ivol)
                            {
                              const IntVectSet& ivs = (*m_dataBou[jvol][ilev])[dit()].getIVS();
                              if(ivs.contains(vofit().gridIndex()))
                                {
                                  found = true;
                                  otherVals += (*m_dataBou[jvol][ilev])[dit()](vofit(), ivar);
                                  nvals ++;
                                }
                            }
                        }
                      Real flux = 0;
                      if(found)
                        {
                          flux = ((otherVals/nvals) - thisVal);
                          //flux = ((otherVals/nvals) + thisVal)/2;
                        }
                      else
                        {
                          //no other volumes found.   in practice, usually a regular cell 
                          //that has been labeled irregular
                        }

                      (*m_bounVal[ivol][ilev])[dit()](vofit(), ivar) = flux;
                    }
                }
            }
          dxlev /= m_params.m_refRatio[ilev];
        }
    }
}
///
void 
AmoebaSolver::
advanceOneVariable(int a_ivar)
{
  CH_TIME("advanceOneVariable");
  Interval zeroint(0,0);
  Interval ivarint(a_ivar,a_ivar);
  //solver is for a single variable.  copy solution and rhs to scratch space
  for(int ivol = 0; ivol < m_volumes.size(); ivol++)
    {
      CH_TIME("copy_to_scratch");
      for(int ilev = 0; ilev < m_params.m_numLevels; ilev++)
        {
          m_solnOld[ivol][ilev]->copyTo(ivarint, *m_scalOld[ivol][ilev], zeroint);
          m_solnNew[ivol][ilev]->copyTo(ivarint, *m_scalNew[ivol][ilev], zeroint);
          m_soursin[ivol][ilev]->copyTo(ivarint, *m_scalRHS[ivol][ilev], zeroint);
          m_bounVal[ivol][ilev]->copyTo(ivarint, *m_scalBou[ivol][ilev], zeroint);
        }
    }

  // Advance one time step
  for(int ivol = 0; ivol < m_volumes.size(); ivol++)
    {

      pout() << "advancing volume " << ivol << ", variable  "<< a_ivar <<  " in time " << endl;

      CH_TIME("integrator::onestep");


      m_integrator[ivol][a_ivar]->oneStep(m_scalNew[ivol],
                                          m_scalOld[ivol],
                                          m_scalRHS[ivol],
                                          m_params.m_dt,
                                          0,
                                          m_params.m_numLevels-1,
                                          true);

    }

  //copy stuff back from scratch
  for(int ivol = 0; ivol < m_volumes.size(); ivol++)
    {
      CH_TIME("copy_from_scratch");
      EBAMRDataOps::assign(m_solnNew[ivol], m_scalNew[ivol], ivarint, zeroint);
    }
}
///
void 
AmoebaSolver::
run()
{
  CH_TIME("AmoebaSolver::run");

  // Compute the number of time steps
  int numSteps = m_params.m_endTime / m_params.m_dt;

  // Iterate until the end time is reached
  int step;
  for (step = 0; step < numSteps; step++)
    {
      // Set and print the current time
      Real time = step * m_params.m_dt;
      m_time = time;
      pout() << "time = " << time << "\n";


      // Set the source/sink term (combined here)
      setSource();

      /// Taylor series extrapolate data to the boundary for use in boundary conditions
      extrapolateDataToBoundary();

      //set boundary condtions on embedded boundary
      setBoundaryValues();

      // Write the solution if it's the right time
      if (m_params.m_outputInterval > 0 &&
          step % m_params.m_outputInterval == 0)
        {
          writeOutput(step,time);
        }

      if((step == 0 ) || (!m_params.m_constCoeff))
        {
          m_integrator.resize(m_volumes.size());
          for(int ivol = 0; ivol < m_volumes.size(); ivol++)
            {
              m_integrator[ivol].resize(m_params.m_ncomp);
              for(int ivar = 0; ivar <  m_params.m_ncomp; ivar++)
                {
                  defineSolver(m_integrator[ivol][ivar], ivol, ivar);
                }
            }
        }
      //advance the solution
      for(int ivar = 0; ivar < m_params.m_ncomp; ivar++)
        {
          advanceOneVariable(ivar);
        }

      // Copy the new solution to the old solution
      for(int ivol = 0; ivol < m_volumes.size(); ivol++)
        {
          CH_TIME("AmoebaSolver::copy");
          EBAMRDataOps::assign(m_solnOld[ivol],m_solnNew[ivol]);
        }
    }

  // Write the solution at the end
  if (m_params.m_outputInterval > 0)
    {
      writeOutput(step,m_params.m_endTime);
    }
}
///
void 
AmoebaSolver::
getDiffusionConstants()
{
  ParmParse pp;
  pp.getarr("diffusion_constants", m_diffusionConstants, 0, m_volumes.size());
  if(!m_params.m_constCoeff)
    {
      pp.getarr("diffusion_timeCoefs", m_diffusionTimeCoefs, 0, m_volumes.size());
      m_diffusionGradients.resize(m_volumes.size());
      for(int ivol = 0; ivol < m_diffusionConstants.size(); ivol++)
        {
          char gradname[100];
          sprintf(gradname, "diffusion_gradient%d", ivol);
          Vector<Real> vecgrad;
          pp.getarr(gradname, vecgrad, 0, SpaceDim);
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              m_diffusionGradients[ivol][idir] = vecgrad[idir];
            }
        }
      for(int ivol = 0; ivol < m_diffusionConstants.size(); ivol++)
        pout() << "diffusion gradient[" << ivol << "] = " << m_diffusionGradients[ivol] << "\n";
      for(int ivol = 0; ivol < m_diffusionConstants.size(); ivol++)
        pout() << "diffusion timecoef[" << ivol << "] = " << m_diffusionTimeCoefs[ivol] << "\n";
    }
  for(int ivol = 0; ivol < m_diffusionConstants.size(); ivol++)
    pout() << "diffusion constant[" << ivol << "] = " << m_diffusionConstants[ivol] << "\n";
}
///
void 
AmoebaSolver::
initGeometry()
{
  CH_TIME("AmoebaSolver::initGeometry");

  // Compute the finest domain and cell size
  ProblemDomain fineDomain = m_params.m_coarsestDomain;
  RealVect fineDx = m_params.m_dx * RealVect::Unit;

  for (int ilev = 0; ilev < m_params.m_numLevels-1; ilev++)
    {
      fineDomain.refine(m_params.m_refRatio[ilev]);
      fineDx /= m_params.m_refRatio[ilev];
    }

  // Create the implicit function for one of the known geometries
  BaseIF* geometry;
  MFIndexSpace mfIndexSpace;

  geometry = makeMultiSphereGeometry(mfIndexSpace, fineDomain.domainBox(), m_params.m_loCorner, fineDx[0]);
  getConnectedComponents(m_volumes, mfIndexSpace);
  getDiffusionConstants();

}
///
void 
AmoebaSolver::
initIndexSpace()
{
  CH_TIME("AmoebaSolver::initIndexSpace");

  // Create a set of grids at the coarsest level
  Vector<int> coarsestProcs;
  Vector<Box> coarsestBoxes;
  
  domainSplit(m_params.m_coarsestDomain,
              coarsestBoxes,
              m_params.m_maxBoxSize,
              m_params.m_blockFactor);

  // Order them using a space filling curve
  mortonOrdering(coarsestBoxes);

  LoadBalance(coarsestProcs,coarsestBoxes);

  // Create more of the infrastructure for geometry
  m_grids.resize(m_params.m_numLevels);
  m_ebisl.resize(m_volumes.size());
  m_grids[0] = DisjointBoxLayout(coarsestBoxes,
                                 coarsestProcs,
                                 m_params.m_coarsestDomain);

  for(int ivol = 0; ivol < m_volumes.size(); ivol++)
    {
      m_ebisl[ivol].resize(m_params.m_numLevels);
      
  
      m_volumes[ivol]->fillEBISLayout(m_ebisl[ivol][0],
                                      m_grids[0],
                                      m_params.m_coarsestDomain,
                                      m_params.m_numGhostEBISLayout);
    }

  if (m_params.m_numLevels > 1)
    {
      // If there is more than one level, the finer levels need to created by
      // "BRMeshRefine"
      BRMeshRefine meshRefine;

      meshRefine.define(m_params.m_coarsestDomain,
                        m_params.m_refRatio,
                        m_params.m_fillRatio,
                        m_params.m_blockFactor,
                        m_params.m_nestingRadius,
                        m_params.m_maxBoxSize);

      // Compute the second finest domain
      ProblemDomain secondFinestDomain = m_params.m_coarsestDomain;

      for (int ilev = 0; ilev < m_params.m_numLevels-2; ilev++)
        {
          secondFinestDomain.refine(m_params.m_refRatio[ilev]);
        }

      // Tags for creating the finer levels
      Vector<IntVectSet> tags(m_params.m_numLevels);

      for(int ivol = 0; ivol < m_volumes.size(); ivol++)
        {
          // Get the depth of the second to finest level
          int depth = m_volumes[ivol]->getLevel(secondFinestDomain);
          IntVectSet tagsVol = m_volumes[ivol]->irregCells(depth);
          tagsVol.grow(2);
          tags[m_params.m_numLevels-2] |= tagsVol;
        }


      Vector<Vector<Box> > oldBoxes(m_params.m_numLevels);
      Vector<Vector<Box> > newBoxes;

      // Need the boxes on the coarsest level and the tags on the second to
      // finest level to make all the boxes on all the levels
      oldBoxes[0] = coarsestBoxes;

      // Make all the boxes on all the levels
      meshRefine.regrid(newBoxes,
                        tags,
                        0,
                        m_params.m_numLevels-1,
                        oldBoxes);

      // Go through all the new levels, Morton order the boxes, load balance the
      // result, and create the data structures needed
      ProblemDomain curDomain = m_params.m_coarsestDomain;
      for (int ilev = 1; ilev < m_params.m_numLevels; ilev++)
        {
          Vector<int> curProcs;

          curDomain.refine(m_params.m_refRatio[ilev-1]);

          // Order them using a space filling curve
          mortonOrdering(newBoxes[ilev]);

          // Do load balancing using timing information
          LoadBalance(curProcs, newBoxes[ilev]);


          // Create more of the infrastructure for geometry
          m_grids[ilev] = DisjointBoxLayout(newBoxes[ilev],
                                            curProcs,
                                            curDomain);

          for(int ivol = 0; ivol < m_volumes.size(); ivol++)
            {
              m_volumes[ivol]->fillEBISLayout(m_ebisl[ivol][ilev],
                                              m_grids[ilev],
                                              curDomain,
                                              m_params.m_numGhostEBISLayout);
            }
        }
    }
  ParmParse pp;
  if(pp.contains("abr_file_output_name"))
    {
      string abrFile;
      pp.get("abr_file_output_name",  abrFile);
      Vector< Vector<Box> > boxes(m_params.m_numLevels);
      Vector< Vector<int> > ranks(m_params.m_numLevels);
      for(int ilev = 0; ilev < m_params.m_numLevels; ilev++)
        {
          boxes[ilev] = m_grids[ilev].boxArray();
          ranks[ilev] = m_grids[ilev].procIDs();
        }
      Box baseProbBox = m_params.m_coarsestDomain.domainBox();
      writeABRfile(boxes, ranks, 
                   m_params.m_refRatio, 
                   m_params.m_numLevels,
                   numProc(),
                   baseProbBox,
                   abrFile);
    }
}
/////
void
AmoebaSolver::
getEBLGAndQuadCFI(Vector<EBLevelGrid>                   & a_ebLevelGrids,
                  Vector<RefCountedPtr<EBQuadCFInterp> >& a_quadCFInterp,
                  int a_ivol, int a_ncomp)
{
  a_ebLevelGrids.resize(m_grids.size());
  a_quadCFInterp.resize(m_grids.size());

  // Define the data holders and interpolators
  ProblemDomain levelDomain = m_params.m_coarsestDomain;
  ProblemDomain coarserDomain;
  for (int ilev = 0; ilev < m_grids.size(); ilev++)
    {
      a_ebLevelGrids[ilev].define(m_grids[ilev],m_ebisl[a_ivol][ilev],levelDomain);

      if (ilev > 0)
        {
          int numVariables = a_ncomp;

          a_quadCFInterp[ilev] = RefCountedPtr<EBQuadCFInterp>
            (new EBQuadCFInterp(m_grids[ilev],
                                m_grids[ilev-1],
                                m_ebisl[a_ivol][ilev],
                                m_ebisl[a_ivol][ilev-1],
                                coarserDomain,
                                m_params.m_refRatio[ilev-1],
                                numVariables,
                                *(a_ebLevelGrids[ilev].getCFIVS()),
                                &(*m_volumes[a_ivol])));
        }

      coarserDomain = levelDomain;

      if (ilev < m_grids.size()-1)
        {
          levelDomain.refine(m_params.m_refRatio[ilev]);
        }
    }
}

//////
void 
AmoebaSolver::
setConductivityCoefs(LevelData<EBCellFAB>          &    a_aco,
                     LevelData<EBFluxFAB>          &    a_bco,
                     LevelData<BaseIVFAB<Real> >   &    a_bcoIrreg,
                     const Real                    &    a_dx,
                     int a_ivol, int a_ivar)
{
  CH_TIME("DiffusionSolver::setConductivityCoefs");
  const DisjointBoxLayout& dbl = a_aco.disjointBoxLayout();
  for(DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      a_aco[dit()].setVal(1.0);
      const EBISBox& ebisBox = a_aco[dit()].getEBISBox();
      const Box&         box = dbl.get(dit());
      IntVectSet ivsIrreg = ebisBox.getIrregIVS(box);
      BaseIVFAB<Real>& irregb = a_bcoIrreg[dit()];

      for(VoFIterator vofit(ivsIrreg, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          RealVect bndryCentroid = ebisBox.bndryCentroid(vofit());
          const IntVect& iv = vofit().gridIndex();
          int ihere = 0;
          if((iv[0] ==32)  && iv[1]==15)
            {
              ihere = 1;
            }
          bndryCentroid *= a_dx;
          RealVect bndryLoc = EBArith::getVofLocation(vofit(), a_dx*RealVect::Unit, bndryCentroid);
          Real bcopoint;
          FORT_POINTGETBCOAMOEBA(CHF_REAL(bcopoint),
                                 CHF_REALVECT(bndryLoc),
                                 CHF_REAL(m_time),
                                 CHF_CONST_INT(a_ivol),
                                 CHF_CONST_INT(a_ivar),
                                 CHF_REAL(    m_diffusionConstants[a_ivol]),
                                 CHF_REALVECT(m_diffusionGradients[a_ivol]),
                                 CHF_REAL(    m_diffusionTimeCoefs[a_ivol])
                                 );

          irregb(vofit(), 0) = bcopoint;
        }
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          EBFaceFAB&       bcodir = a_bco[dit()][idir];
          BaseFab<Real>& svbcodir =   bcodir.getSingleValuedFAB();
          //do all as regular faces first
          FORT_GETBCOAMOEBA(CHF_FRA1(svbcodir,0),
                            CHF_CONST_REAL(m_params.m_dx),
                            CHF_CONST_REAL(m_time),
                            CHF_CONST_INT(a_ivol),
                            CHF_CONST_INT(a_ivar),
                            CHF_BOX(svbcodir.box()),
                            CHF_REAL(    m_diffusionConstants[a_ivol]),
                            CHF_REALVECT(m_diffusionGradients[a_ivol]),
                            CHF_REAL(    m_diffusionTimeCoefs[a_ivol]),
                            CHF_CONST_INT(idir)
                            );
          
          for(FaceIterator faceit(ivsIrreg, ebisBox.getEBGraph(), idir, FaceStop::SurroundingWithBoundary); faceit.ok(); ++faceit)
            {
              const FaceIndex& face =  faceit();
              if(box.contains(face.gridIndex(Side::Lo)) || box.contains(face.gridIndex(Side::Hi)))
                {
                  RealVect bndryCentroid = ebisBox.centroid(face);
                  bndryCentroid *= a_dx;
                  RealVect faceLoc = EBArith::getFaceLocation(face, a_dx*RealVect::Unit, bndryCentroid);
                  Real bcopoint;
                  FORT_POINTGETBCOAMOEBA(CHF_REAL(bcopoint),
                                         CHF_REALVECT(faceLoc),
                                         CHF_REAL(m_time),
                                         CHF_CONST_INT(a_ivol),
                                         CHF_CONST_INT(a_ivar),
                                         CHF_REAL(    m_diffusionConstants[a_ivol]),
                                         CHF_REALVECT(m_diffusionGradients[a_ivol]),
                                         CHF_REAL(    m_diffusionTimeCoefs[a_ivol])
                                         );
                  bcodir(face, 0) = bcopoint;
                }
            }
        } //end loop over directions
    }
}
///
void
AmoebaSolver::
defineVariableCoeffs(Vector<RefCountedPtr<LevelData<EBCellFAB> > >&           a_aco,
                     Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&           a_bco,
                     Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&    a_bcoIrreg,
                     int a_ivol, int a_ivar)
{
  CH_TIME("AmoebaSolver::defineConductivityCoef");
  a_aco.resize(        m_params.m_numLevels);
  a_bco.resize(        m_params.m_numLevels);
  a_bcoIrreg.resize(   m_params.m_numLevels);
  Real dxLev = m_params.m_dx;
  ProblemDomain domLev =   m_params.m_coarsestDomain;
  for(int ilev = 0; ilev < m_params.m_numLevels; ilev++)
    {
      LayoutData<IntVectSet> irregSets(m_grids[ilev]);
      int nghost = 1;
      for(DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          //has to correspond to number of ghost cells
          Box grownBox = grow(m_grids[ilev].get(dit()), nghost);
          grownBox &= domLev;
          irregSets[dit()] = m_ebisl[a_ivol][ilev][dit()].getIrregIVS(grownBox);
        }
      EBFluxFactory        ebfluxfact(m_ebisl[a_ivol][ilev]);
      EBCellFactory        ebcellfact(m_ebisl[a_ivol][ilev]);
      BaseIVFactory<Real>  baseivfact(m_ebisl[a_ivol][ilev], irregSets);

      a_aco[ilev]         = RefCountedPtr<LevelData<EBCellFAB       > >(new LevelData<EBCellFAB       >(m_grids[ilev], 1, nghost*IntVect::Unit, ebcellfact));
      a_bco[ilev]         = RefCountedPtr<LevelData<EBFluxFAB       > >(new LevelData<EBFluxFAB       >(m_grids[ilev], 1, nghost*IntVect::Unit, ebfluxfact));
      a_bcoIrreg[ilev]    = RefCountedPtr<LevelData<BaseIVFAB<Real> > >(new LevelData<BaseIVFAB<Real> >(m_grids[ilev], 1, nghost*IntVect::Unit, baseivfact));

      setConductivityCoefs(*a_aco[ilev], *a_bco[ilev], *a_bcoIrreg[ilev], dxLev, a_ivol, a_ivar);
      dxLev /=      m_params.m_refRatio[ilev];
      domLev.refine(m_params.m_refRatio[ilev]);
    }
}
///
void
AmoebaSolver::
getConstantCoeffOpFactory(RefCountedPtr<EBAMRPoissonOpFactory>& a_factory, int a_ivol, int a_ivar)
{
  // Set up the no flux domain and embedded boundary conditions
  RefCountedPtr<NeumannPoissonDomainBCFactory> domBC(new NeumannPoissonDomainBCFactory());
  RefCountedPtr<NeumannPoissonEBBCFactory>      ebBC(new NeumannPoissonEBBCFactory());

  //RefCountedPtr<DirichletPoissonDomainBCFactory> domBC(new DirichletPoissonDomainBCFactory());
  //RefCountedPtr<DirichletPoissonEBBCFactory>      ebBC(new DirichletPoissonEBBCFactory());
  //ebBC->setOrder(1);

  domBC->setValue(0.0);
  ebBC->setValue(0.);

  Vector<EBLevelGrid>  eblg;
  Vector<RefCountedPtr<EBQuadCFInterp> > quadCFI;
  getEBLGAndQuadCFI(eblg, quadCFI, a_ivol);

  //coefficients come in through the =coefficients.
    //  pout() << "using multicolored gauss seidel" << endl;
  int relaxType= m_params.m_mgRelaxType;
  if(relaxType == 1)
    {
      pout() << "using multi-colored gauss seidel relaxation" << endl;
    }
  else if(relaxType == 2)
    {
      pout() << "using gsrb fast relaxation" << endl;
    }
  else
    {
      MayDay::Error("bogus relaxType for variable coefficients");
    }

  int precond = m_params.m_mgNumPrecondIter;
  Real alpha = 1.;
  a_factory = RefCountedPtr<EBAMRPoissonOpFactory>
    (new EBAMRPoissonOpFactory(eblg, m_params.m_refRatio, quadCFI, 
                               m_params.m_dx*RealVect::Unit,  RealVect::Zero,  precond, relaxType, 
                               domBC, ebBC, alpha, m_diffusionConstants[a_ivol], m_time,
                               m_params.m_numGhostSoln, m_params.m_numGhostSource));

  a_factory->setData(m_scalBou[a_ivol]);

}
///
void 
AmoebaSolver::
getVariableCoeffOpFactory(RefCountedPtr<EBConductivityOpFactory>& a_factory, int a_ivol, int a_ivar)
{
  // Set up the no flux domain and embedded boundary conditions
  RefCountedPtr<NeumannConductivityDomainBCFactory> domBC(new NeumannConductivityDomainBCFactory());
  domBC->setValue(0.0);

  //RefCountedPtr<DirichletConductivityDomainBCFactory> domBC(new DirichletConductivityDomainBCFactory());
  //domBC->setValue(0.0);

  //RefCountedPtr<DirichletConductivityEBBCFactory>      ebBC(new DirichletConductivityEBBCFactory());
  //ebBC->setOrder(1);
  //ebBC->setValue(0.0);

  RefCountedPtr<NeumannConductivityEBBCFactory>      ebBC(new NeumannConductivityEBBCFactory());
  ebBC->setValue(0.);

  Vector<EBLevelGrid>  eblg;
  Vector<RefCountedPtr<EBQuadCFInterp> > quadCFI;
  getEBLGAndQuadCFI(eblg, quadCFI, a_ivol);
  Vector<RefCountedPtr<LevelData<EBCellFAB> > >           aco;
  Vector<RefCountedPtr<LevelData<EBFluxFAB> > >           bco;
  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >    bcoIrreg;

  defineVariableCoeffs(aco, bco, bcoIrreg, a_ivol, a_ivar);

  //coefficients come in through the =coefficients.
  Real unity = 1.0;
  //  int relaxType = 0;
  //  pout() << "using multicolored gauss seidel" << endl;
  int relaxType= m_params.m_mgRelaxType;
  if(relaxType == 1)
    {
      pout() << "using multi-colored gauss seidel relaxation" << endl;
    }
  else if(relaxType == 0)
    {
      pout() << "using point jacobi relaxation" << endl;
    }
  else if(relaxType == 2)
    {
      pout() << "using gsrb fast relaxation" << endl;
    }
  else
    {
      MayDay::Error("bogus relaxType for variable coefficients");
    }

  a_factory = RefCountedPtr<EBConductivityOpFactory>
    (new EBConductivityOpFactory(eblg, quadCFI, unity, unity, aco, bco, bcoIrreg,
                                 m_params.m_dx,  m_params.m_refRatio, domBC, ebBC,
                                 m_params.m_numGhostSoln, m_params.m_numGhostSource, relaxType));
  a_factory->setData(m_scalBou[a_ivol]);

}

////////
void 
AmoebaSolver::
initData()
{

  CH_TIME("AmoebaSolver::initData");

  ParmParse pp;
  pp.getarr("initial_value",m_params.m_initialValue, 0, m_volumes.size());

  for(int ivol = 0; ivol < m_volumes.size(); ivol++)
    {
      pout() << "initial value for volume " << ivol << " =  " << m_params.m_initialValue[ivol] << endl;
    }

  // Make pointers for data holders at each AMR level

  m_solnOld.resize(m_volumes.size());
  m_solnNew.resize(m_volumes.size());
  m_soursin.resize(m_volumes.size());
  m_bounVal.resize(m_volumes.size());
  m_dataBou.resize(m_volumes.size());
  m_scalBou.resize(m_volumes.size());
  m_scalOld.resize(m_volumes.size());
  m_scalNew.resize(m_volumes.size());
  m_scalRHS.resize(m_volumes.size());
  m_irrSets.resize(m_volumes.size());

  for(int ivol = 0; ivol < m_volumes.size(); ivol++)
    {
      m_solnOld[ivol].resize(m_params.m_numLevels);
      m_solnNew[ivol].resize(m_params.m_numLevels);
      m_soursin[ivol].resize(m_params.m_numLevels);
      m_bounVal[ivol].resize(m_params.m_numLevels);
      m_dataBou[ivol].resize(m_params.m_numLevels);
      m_scalBou[ivol].resize(m_params.m_numLevels);
      m_scalOld[ivol].resize(m_params.m_numLevels);
      m_scalNew[ivol].resize(m_params.m_numLevels);
      m_scalRHS[ivol].resize(m_params.m_numLevels);
      m_irrSets[ivol].resize(m_params.m_numLevels);

      for (int ilev = 0; ilev < m_params.m_numLevels; ilev++)
        {
          m_irrSets[ivol][ilev] = RefCountedPtr<LayoutData<IntVectSet> >(new LayoutData<IntVectSet>(m_grids[ilev]));
          for(DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
            {
              (*m_irrSets[ivol][ilev])[dit()] = m_ebisl[ivol][ilev][dit()].getIrregIVS(m_grids[ilev][dit()]);
            }

          // At each level, allocate and initialize (as needed) space for the old
          // solution, the new solution, and the source/sink terms
          EBCellFactory        ebCellFactory(m_ebisl[ivol][ilev]);
          BaseIVFactory<Real>  bivfabFactory(m_ebisl[ivol][ilev], *m_irrSets[ivol][ilev]);
          m_dataBou[ivol][ilev] = RefCountedPtr<LevelData< BaseIVFAB<Real> > >(new LevelData< BaseIVFAB<Real> >(m_grids[ilev], m_params.m_ncomp, IntVect::Zero, bivfabFactory));
          m_bounVal[ivol][ilev] = RefCountedPtr<LevelData< BaseIVFAB<Real> > >(new LevelData< BaseIVFAB<Real> >(m_grids[ilev], m_params.m_ncomp, IntVect::Zero, bivfabFactory));
          m_scalBou[ivol][ilev] = RefCountedPtr<LevelData< BaseIVFAB<Real> > >(new LevelData< BaseIVFAB<Real> >(m_grids[ilev],                1, IntVect::Zero, bivfabFactory));
          m_solnOld[ivol][ilev] = new LevelData<EBCellFAB>(m_grids[ilev], m_params.m_ncomp, m_params.m_numGhostSoln,   ebCellFactory);
          m_solnNew[ivol][ilev] = new LevelData<EBCellFAB>(m_grids[ilev], m_params.m_ncomp, m_params.m_numGhostSoln,   ebCellFactory);
          m_soursin[ivol][ilev] = new LevelData<EBCellFAB>(m_grids[ilev], m_params.m_ncomp, m_params.m_numGhostSource, ebCellFactory);
          m_scalOld[ivol][ilev] = new LevelData<EBCellFAB>(m_grids[ilev],                1, m_params.m_numGhostSoln,   ebCellFactory);
          m_scalNew[ivol][ilev] = new LevelData<EBCellFAB>(m_grids[ilev],                1, m_params.m_numGhostSoln,   ebCellFactory);
          m_scalRHS[ivol][ilev] = new LevelData<EBCellFAB>(m_grids[ilev],                1, m_params.m_numGhostSource, ebCellFactory);
          EBLevelDataOps::setVal(*(m_solnOld[ivol][ilev]),m_params.m_initialValue[ivol]);
          EBLevelDataOps::setVal(*(m_solnNew[ivol][ilev]),m_params.m_initialValue[ivol]);

          // At each level, count and report the number of irregular cells and the
          // number of computational cells
          Real totalFortran = 0.0;
          Real totalIrreg   = 0.0;

          DisjointBoxLayout     fineGrids   = m_grids[ilev];
          LevelData<EBCellFAB>* fineDataPtr = m_solnOld[ivol][ilev];

          for (DataIterator dit=fineGrids.dataIterator(); dit.ok(); ++dit)
            {
              const Box&        curBox      = fineGrids[dit()];
              const EBISBox&    curEBISBox  = (*fineDataPtr)[dit()].getEBISBox();
              const IntVectSet& curIrregIVS = curEBISBox.getIrregIVS(curBox);

              totalFortran += curBox.numPts();
              totalIrreg   += curIrregIVS.numPts();
            }

#ifdef CH_MPI
          Real mpiTotal;
          int status = MPI_Allreduce(&totalFortran, &mpiTotal, 1, MPI_CH_REAL,
                                     MPI_SUM, Chombo_MPI::comm);

          if (status != MPI_SUCCESS)
            {
              MayDay::Error("MPI error summing 'totalFortran'");
            }

          totalFortran = mpiTotal;

          status = MPI_Allreduce(&totalIrreg, &mpiTotal, 1, MPI_CH_REAL,
                                 MPI_SUM, Chombo_MPI::comm);

          if (status != MPI_SUCCESS)
            {
              MayDay::Error("MPI error summing 'totalIrreg'");
            }

          totalIrreg = mpiTotal;
#endif

          long long totalFortranInt = totalFortran;
          long long totalIrregInt   = totalIrreg;

          std::ios::fmtflags origFlags = pout().flags();
          int origWidth = pout().width();
          int origPrecision = pout().precision();

          pout() << "volume " << ivol << ":" << "\n";
          pout() << "Level " << ilev << ":" << "\n";
          pout() << setiosflags(ios::right);
          pout() << "  Total computation cells: " << setw(10) << totalFortranInt << "\n";
          pout() << "  Total irregular cells:   " << setw(10) << totalIrregInt << "\n";
          pout() << "\n";

          pout().flags(origFlags);
          pout().width(origWidth);
          pout().precision(origPrecision);
        } //end loop over levels
    } //end loop over volumes
}
///
void 
AmoebaSolver::
setSource()
{
  CH_TIME("AmoebaSolver::setSource");

  Real dx = m_params.m_dx;

  for (int ivol = 0; ivol < m_volumes.size(); ivol++)
    {
      // Set the source/sink term on each level
      for (int ilev = 0; ilev < m_params.m_numLevels; ilev++)
        {
          // Grids on this level
          DisjointBoxLayout& curGrids = m_grids[ilev];

          // Data holders for the source/sink and old solution on this level
          LevelData<EBCellFAB>* sourceLevel = m_soursin[ivol][ilev];
          LevelData<EBCellFAB>* solnLevel   = m_solnOld[ivol][ilev];

          // Set the source/sink term for each grid
          for (DataIterator dit=curGrids.dataIterator(); dit.ok(); ++dit)
            {
              // Current grid
              const Box& curBox = curGrids[dit()];

              // Single grid data holders
              EBCellFAB& sourceEBCellFAB = (*sourceLevel)[dit()];
              EBCellFAB& solnEBCellFAB   = (*solnLevel)[dit()];

              // Single grid data holders for single valued cells  
              BaseFab<double>& sourceFAB = sourceEBCellFAB.getSingleValuedFAB();
              BaseFab<double>& solnFAB   = solnEBCellFAB.getSingleValuedFAB();

              // Check the geometry type to determine the source/sink to use
              // Call Fortran to set the source/sink for all the single valued cells
              for(int ivar = 0; ivar < m_params.m_ncomp; ivar++)
                {
                  FORT_AMOEBASIMPLESOURCESINK(CHF_FRA1(sourceFAB,ivar),
                                              CHF_CONST_FRA1(solnFAB,ivar),
                                              CHF_CONST_REAL(m_params.m_sourceScaling),
                                              CHF_CONST_REAL(m_params.m_sinkScaling),
                                              CHF_CONST_REAL(dx),
                                              CHF_CONST_INT(ivol),
                                              CHF_CONST_INT(ivar),
                                              CHF_CONST_REALVECT(m_params.m_loCorner),
                                              CHF_BOX(curBox));
                }

              // Get all the machinery need to update multi-valued cells
              const EBISBox& curEBISBox = sourceEBCellFAB.getEBISBox();
              const EBGraph& curEBGraph = curEBISBox.getEBGraph();
              IntVectSet multiCells = curEBISBox.getMultiCells(curBox);

              // Set the source/sink for all the multi-value cells
              for (VoFIterator vofit(multiCells,curEBGraph); vofit.ok(); ++vofit)
                {
                  const VolIndex& volIndex = vofit();
                  const IntVect&  iv = volIndex.gridIndex();

                  Real source;
                  
                  for(int ivar = 0; ivar < m_params.m_ncomp; ivar++)
                    {
                      // Call Fortran for a single value and assign the result
                      FORT_AMOEBASIMPLESOURCESINKPOINT(CHF_REAL(source),
                                                       CHF_CONST_REAL(solnEBCellFAB(volIndex,0,1)),
                                                       CHF_CONST_REAL(m_params.m_sourceScaling),
                                                       CHF_CONST_REAL(m_params.m_sinkScaling),
                                                       CHF_CONST_REAL(dx),
                                                       CHF_CONST_INT(ivol),
                                                       CHF_CONST_INT(ivar),
                                                       CHF_CONST_REALVECT(m_params.m_loCorner),
                                                       CHF_CONST_INTVECT(iv));

                      sourceEBCellFAB(volIndex,ivar) = source;
                    }

                }
            }

          dx /= m_params.m_refRatio[ilev];
        }
    }
}
///
void 
AmoebaSolver::
writeOutput(const int  & a_step,
            const Real & a_time)
{
  CH_TIME("AmoebaSolver::writeOutput");

  for(int ivol = 0; ivol < m_volumes.size(); ivol++)
    {

      // Generate the full output file name
      char suffix[128];
      sprintf(suffix,".%06d.vol%d.%dd.hdf5",a_step,ivol,SpaceDim);

      string filename = m_params.m_outputPrefix + suffix;

      // Name the solution variable
      Vector<string> names(m_params.m_ncomp);
      for(int ivar = 0; ivar < m_params.m_ncomp; ivar++)
        {
          char charstr[100];
          sprintf(charstr, "G.var%d",ivar);
          names[ivar] = string(charstr);
        }

      // Leave the values in the covered cells alone
      bool replaceCovered = false;
      Vector<Real> dummy;

      // Write the output file
      writeEBHDF5(filename,
                  m_grids,
                  m_solnOld[ivol],
                  names,
                  m_params.m_coarsestDomain,
                  m_params.m_dx,
                  m_params.m_dt,
                  a_time,
                  m_params.m_refRatio,
                  m_params.m_numLevels,
                  replaceCovered,
                  dummy);
      if(a_step == 0)
        {
          char rhssuffix[128];
          sprintf(rhssuffix,".vol%d.%dd.hdf5",ivol,SpaceDim);
          string rhsname = string("rhs.") +  m_params.m_outputPrefix + rhssuffix;
          // Name the solution variable
          Vector<string> rhsnames(m_params.m_ncomp);
          for(int ivar = 0; ivar < m_params.m_ncomp; ivar++)
            {
              char charstr[100];
              sprintf(charstr, "rhs.var%d",ivar);
              rhsnames[ivar] = string(charstr);
            }
          // Write the output file
          writeEBHDF5(rhsname,
                      m_grids,
                      m_soursin[ivol],
                      names,
                      m_params.m_coarsestDomain,
                      m_params.m_dx,
                      m_params.m_dt,
                      a_time,
                      m_params.m_refRatio,
                      m_params.m_numLevels,
                      replaceCovered,
                      dummy);
        }
    }
}
///
