#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef __AmoebaSolver_H_
#define __AmoebaSolver_H_

#include "Vector.H"
#include "AggStencil.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "BiCGStabSolver.H"
#include "EBQuadCFInterp.H"
#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "EBLevelGrid.H"
#include "EBConductivityOpFactory.H"
#include "EBConductivityOp.H"
#include "EBAMRPoissonOp.H"
#include "EBAMRPoissonOpFactory.H"
#include "EBBackwardEuler.H"

/// A class to hold all the solver parameters
class AmoebaParams
{
public:
  /// Null constructor
  AmoebaParams();

  /// Print the parameters
  void print();


  /// Strengths of the source and sink
  Real m_sourceScaling;
  Real m_sinkScaling;

  /// Solution at t = 0 (different constant for each volume)
  Vector<Real> m_initialValue;

  /// The physical domain corner and cell size (at the coarsest level)
  RealVect m_loCorner;
  Real m_dx;
    
  /// The number of AMR levels
  int m_numLevels;
  /// Number of cells in each direction at the coarsest level
  IntVect m_numCells;
  /// The refinement ratio between each pair of levels
  Vector<int> m_refRatio;

  /// The coarsest problem domain
  ProblemDomain m_coarsestDomain;

  /// Time step and time to end the simulation
  Real m_dt;
  Real m_endTime;

  /// How often to output and a prefix for the output name
  int m_outputInterval;
  string m_outputPrefix;

  /// Used to determine box sizes and alignments
  int m_maxBoxSize;
  int m_blockFactor;

  /// Used to generate AMR grids
  Real   m_fillRatio;
  int    m_nestingRadius;

  /// Multigrid parameters
  int  m_mgNumCycles;
  int  m_mgNumSmooths;
  int  m_mgRelaxType;
  bool m_mgLazyRelax;
  Real m_mgToler;
  Real m_mgHangToler;
  int  m_mgIterMax;
  int  m_mgNumPrecondIter;

  bool m_constCoeff;
  int m_ncomp;
  /// Ghost cell information
  int m_numGhostEBISLayout;
  IntVect m_numGhostSoln;
  IntVect m_numGhostSource;
};

/// A solver for the diffusion equation with source and sink terms
class AmoebaSolver
{
public:
  /// Construct with the solver parameters
  AmoebaSolver(const AmoebaParams & a_params);

  /// Destructor
  ~AmoebaSolver();

  /// Initialize the domain, geometry, and solver
  void init();

  /// Run the solver to completion
  void run();

protected:
  // Initialize the geometry
  void initGeometry();

  // Initialize the index space
  void initIndexSpace();

  // Initialize the data
  void initData();

  //set up extrapolation stencil holders
  void initStencils();

  // Initialize the solver
  void defineSolver(RefCountedPtr<EBBackwardEuler>& a_integrator, int a_ivol, int ivar);

  // Set the current source term
  void setSource();

  /// extrapolate data along rays to the boundary for use in boundary conditions
  void  extrapolateDataToBoundary();

  // Write an output file
  void writeOutput(const int  & a_step,
                   const Real & a_time);

  void getEBLGAndQuadCFI(Vector<EBLevelGrid>                   & a_ebLevelGrids,
                         Vector<RefCountedPtr<EBQuadCFInterp> >& a_quadCFInterp,
                         int ivol, int ncomp =1);

  void setConductivityCoefs(LevelData<EBCellFAB>          &    a_aco,
                            LevelData<EBFluxFAB>          &    a_bco,
                            LevelData<BaseIVFAB<Real> >   &    a_bcoIrreg,
                            const Real                    &    a_dx,
                            int a_ivol, int a_ivar);
  
  void defineVariableCoeffs(Vector<RefCountedPtr<LevelData<EBCellFAB> > >&           a_aco,
                            Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&           a_bco,
                            Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&    a_bcoIrreg,
                            int a_ivol, int a_ivar);

  void getExtrapStencils(Vector<RefCountedPtr<BaseIndex  > >& a_dstVoFs,
                         Vector<RefCountedPtr<BaseStencil> >& a_stencil,
                         const IntVectSet                   & a_cfivs,
                         const DataIndex                    & a_dit,
                         int ivol, int ilev, Real a_dx);

  void getVariableCoeffOpFactory(RefCountedPtr<EBConductivityOpFactory>& a_factory, int a_ivol, int a_ivar);

  void getConstantCoeffOpFactory(RefCountedPtr<EBAMRPoissonOpFactory>& a_factory, int a_ivol, int a_ivar);

  void getDiffusionConstants();
  void initSolverVariableCoeff();
  void advanceOneVariable(int a_ivar);
  void setBoundaryValues();

  /// All the solver parameters
  AmoebaParams m_params;

  /// The bottom solver used by multigrid when it can no longer coarsen
  BiCGStabSolver<LevelData<EBCellFAB> > m_bottomSolver;

  /// All the grids at all AMR levels --- all volumes have the same layout.
  Vector<DisjointBoxLayout> m_grids;

  /// EBISLayouts for all AMR levels for all volumes
  Vector< Vector<EBISLayout> > m_ebisl;
  
  //each connected volume
  Vector< RefCountedPtr<EBIndexSpace> > m_volumes;

  /// The solution at the old and new times
  Vector<  Vector<LevelData<EBCellFAB>* > > m_solnOld;
  Vector<  Vector<LevelData<EBCellFAB>* > > m_solnNew;
  

  /// scalars for temporaries
  Vector<  Vector<LevelData<EBCellFAB>* > > m_scalOld;
  Vector<  Vector<LevelData<EBCellFAB>* > > m_scalNew;
  Vector<  Vector<LevelData<EBCellFAB>* > > m_scalRHS;

  /// The source/sink terms (combined into one "source")
  Vector<  Vector<LevelData<EBCellFAB>* > > m_soursin;

  /// data extrapolated from the field to the irregular face.
  Vector< Vector< RefCountedPtr< LevelData<BaseIVFAB<Real> > > > > m_dataBou;

  /// boundary condition at the cut faces
  Vector< Vector< RefCountedPtr< LevelData<BaseIVFAB<Real> > > > > m_bounVal;
  Vector< Vector< RefCountedPtr< LevelData<BaseIVFAB<Real> > > > > m_scalBou;
  //sets of irregular cells
  Vector< Vector< RefCountedPtr< LayoutData<IntVectSet> > > > m_irrSets;

  
  /// Diffusion coefficient = m dot x + b + K time This is b for each connected volume.
  Vector<Real>     m_diffusionConstants;

  /// Diffusion coefficient = m dot x + b + K time.  This is m for each connected volume.
  Vector<RealVect> m_diffusionGradients;

  /// Diffusion coefficient = m dot x + b + K time.  This is K for each connected volume.
  Vector<Real> m_diffusionTimeCoefs;

  //since coeffs are a function of time, keep the time as member data
  Real m_time;
  
  Vector< Vector<RefCountedPtr<EBBackwardEuler> > > m_integrator;
  
  //this is the stencil that extrapolates data to the irregular boundary
  Vector< Vector< LayoutData< RefCountedPtr< AggStencil< EBCellFAB, BaseIVFAB<Real> > > >* > > m_extrStn;
private:
};

#endif
