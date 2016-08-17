#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _MITOCHONDRIA_H_
#define _MITOCHONDRIA_H_

#include "BaseIF.H"
#include "TransformIF.H"

#include "UsingNamespace.H"

class SphereDimpleIF1: public BaseIF
{
public:
  SphereDimpleIF1();

  virtual ~SphereDimpleIF1();

  virtual Real value(const RealVect & a_point) const;

  virtual BaseIF* newImplicitFunction() const;
};

class CylinderDimpleIF1: public BaseIF
{
public:
  CylinderDimpleIF1();

  virtual ~CylinderDimpleIF1();

  virtual Real value(const RealVect & a_point) const;

  virtual BaseIF* newImplicitFunction() const;
};

class MitochondriaIF1: public BaseIF
{
public:
  MitochondriaIF1();

  virtual ~MitochondriaIF1();

  virtual Real value(const RealVect & a_point) const;

  virtual BaseIF* newImplicitFunction() const;

private:
  BaseIF* m_dimple;
#if 0
  TransformIF* m_dimple2;
  TransformIF* m_dimple3;
#endif
};

#endif
