#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef __DIFFUSIONPARAMS_H_
#define __DIFFUSIONPARAMS_H_

#include <string>
using std::string;

#include "IntVect.H"
#include "RealVect.H"
#include "ProblemDomain.H"

#include "UsingNamespace.H"

/// A class to hold all the solver parameters
class DiffusionParams
{
public:
  /// Null constructor
  DiffusionParams();

  /// Print the parameters
  void print();

  /// Diffusion constant
  Real m_diffusionConstant;

  /// Strengths of the source and sink
  Real m_sourceScaling;
  Real m_sinkScaling;

  /// Solution at t = 0 (constant)
  Real m_initialValue;

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
  string m_abrFile;

  /// Multigrid parameters
  int  m_mgNumCycles;
  int  m_mgNumSmooths;
  int  m_mgRelaxType;
  bool m_mgLazyRelax;
  Real m_mgToler;
  Real m_mgHangToler;
  int  m_mgIterMax;
  int  m_mgNumPrecondIter;

  /// Ghost cell information
  IntVect m_numGhostSoln;
  IntVect m_numGhostSource;
};

#endif
