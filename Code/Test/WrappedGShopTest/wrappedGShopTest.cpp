#include <parstream.H>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;
#include <assert.h>
#include <BRMeshRefine.H>
#include <EBISLayout.H>
#include <EBLevelGrid.H>
#include <GeometryShop.H>
#include <MFIndexSpace.H>
#include <ComplementIF.H>
#include <SphereIF.H>
#include <PlaneIF.H>
#include <IntersectionIF.H>
#include <RefCountedPtr.H>
#include <BoundaryAreaRefCrit.H>
#include <WrappedGShop.H>

class MyIF : public BaseIF
{
	Real value(const RealVect& a_point) const
	{
#if CH_SPACEDIM == 1
		double x = a_point[0];
		double y = 6.6250 + 4 * 5.0/16;
		double z = 6.9375 + 8 * 5.0/16;
#elif CH_SPACEDIM == 2
		double x = a_point[0];
		double y = a_point[1];
		double z = 6.9375 + 8 * 5.0/16;
#elif CH_SPACEDIM == 3
		double x = a_point[0];
		double y = a_point[1];
		double z = a_point[2];
#endif
		return  (z - 8.5 + (3.0 / (x + 5.0))) * (z - 8.5 + (3.0 / (x + 5.0))) + (y - 8.0) * (y - 8.0) - 1.0;
	}

	Real value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
	                 const IndexTM<Real,GLOBALDIM>& a_point) const
	{
		Real retval= LARGEREALVAL;
		int maxDir = a_partialDerivative.maxDir(false);
		int derivativeOrder = a_partialDerivative.sum();

#if CH_SPACEDIM == 1
		double x = a_point[0];
		double y = 0.0;
		double z = 0.0;
#elif CH_SPACEDIM == 2
		double x = a_point[0];
		double y = a_point[1];
		double z = 0.0;
#elif CH_SPACEDIM == 3
		double x = a_point[0];
		double y = a_point[1];
		double z = a_point[2];
#endif

		if (derivativeOrder == 0)
		{
			retval = value(a_point);
		}
		else if (derivativeOrder == 1)
		{
			if (maxDir == 0)
			{
				retval = (z - 8.5 + (3.0 / (x + 5.0)))  *  ((-6)/((x + 5.0) * (x + 5.0)));
			}
			else if (maxDir == 1)
			{
				 retval =  2 * (y - 8);
			}
			else
			{
				 retval = 2*(z - 8.5 + (3.0 / (x + 5.0)));
			}
		}
		else
		{
			retval = 0;
      // pout() << "MyIF:  a_point " << a_point << " " << a_partialDerivative << " "  << retval << endl;
		}
//		pout() << "MyIF:  a_point " << a_point << " " << a_partialDerivative << " "  << retval << endl;
		return retval;
	}

	BaseIF* newImplicitFunction() const
	{
		MyIF* ifPtr = new MyIF();
	  return static_cast<BaseIF*>(ifPtr);
	}
};

int main(int argc, char *argv[])
{
	const int NUM_PHASES = 2;
	int maxBoxSize = 32;
	int blockFactor = 8;
	int numGhostEBISLayout = 4;
	const RealVect domainOrigin(D_DECL(4.375, 6.625, 6.9375));
  const int ncells = 16;
	const double dx = 5.0/ncells;
	const RealVect dxDy(D_DECL(dx, dx, dx));
	const Box domainBox(IntVect(D_DECL(0,0,0)), IntVect(D_DECL(ncells/2 - 1, ncells/2 - 1, ncells - 1)));
	ProblemDomain domain(domainBox);

	MyIF *myIf = new MyIF();
	RefCountedPtr<BaseIF> everythingPhase0(myIf);
	bool complement = true;
	RefCountedPtr<BaseIF> everythingPhase1(new ComplementIF(*everythingPhase0, complement));

	Vector<GeometryService*> geometries(NUM_PHASES, NULL);

	Real threshold = sqrt(2.);
	RefCountedPtr<WGSRefinementCriterion> refcrit(new BoundaryAreaRefCrit(threshold));
	pout() << "using new wrapped gshop geometry generation" << endl;
	int minRefine = 0;
	int maxRefine = 10;
  LocalCoordMoveSwitch::s_turnOffMoveLocalCoords = true;
	WrappedGShop* workshop0 = new WrappedGShop(everythingPhase0, domainOrigin, dx, domain, minRefine, maxRefine);
	workshop0->setRefinementCriterion(refcrit);
	workshop0->m_phase = 0;
	geometries[0] = workshop0;

	WrappedGShop* workshop1 = new WrappedGShop(everythingPhase1, domainOrigin, dx, domain, minRefine, maxRefine);
	workshop1->setRefinementCriterion(refcrit);
	workshop1->m_phase = 1;
	geometries[1] = workshop1;

	MFIndexSpace mfIndexSpace;

	mfIndexSpace.define(domainBox, domainOrigin, dx, geometries, maxBoxSize);

	pout() << "Collecting Connected Components from MFIndexSpace" << endl;

	Vector<int> coarsestProcs;
	Vector<Box> coarsestBoxes;

	domainSplit(domain, coarsestBoxes, maxBoxSize, blockFactor);
	mortonOrdering(coarsestBoxes);
	LoadBalance(coarsestProcs, coarsestBoxes);

	DisjointBoxLayout grid = DisjointBoxLayout(coarsestBoxes, coarsestProcs, domain);
	Vector< Vector<RefCountedPtr<EBIndexSpace> > > phaseVolumes;
	Vector< Vector<EBISLayout> > ebisLayouts;
	Vector<int> numPhaseVolumes;
	phaseVolumes.resize(NUM_PHASES);
	ebisLayouts.resize(NUM_PHASES);
	numPhaseVolumes.resize(NUM_PHASES);

	for (int iphase = 0; iphase < NUM_PHASES; iphase ++)
	{
		// Select one index-space
		Chombo_EBIS::alias((EBIndexSpace*) mfIndexSpace.EBIS(iphase));
		EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
		Vector<RefCountedPtr<EBIndexSpace> > volumes = ebisPtr->connectedComponents();
		int numVolumes = volumes.size();
		numPhaseVolumes[iphase] = numVolumes;

		phaseVolumes[iphase].resize(numVolumes);
		ebisLayouts[iphase].resize(numVolumes);
		pout() << "SmallVolumeFraction_Example_OneFluid_copy2: PHASE " << iphase << ", # volumes=" << numVolumes << endl;
		for (int ivol = 0; ivol < numVolumes; ++ ivol)
		{
			phaseVolumes[iphase][ivol] = volumes[ivol];
			phaseVolumes[iphase][ivol]->fillEBISLayout(ebisLayouts[iphase][ivol],
									grid,
									domain,
									numGhostEBISLayout);
		}
	}

	for (int iphase = 0; iphase < NUM_PHASES; ++ iphase)
	{
		for (int ivol = 0; ivol < numPhaseVolumes[iphase]; ++ ivol)
		{
			for(DataIterator dit = grid.dataIterator(); dit.ok(); ++ dit)
			{
				const Box& currBox = grid[dit()];
				const EBISBox& currEBISBox = ebisLayouts[iphase][ivol][dit()];

				const EBGraph& currEBGraph = currEBISBox.getEBGraph();
				IntVectSet irregCells = currEBISBox.getIrregIVS(currBox);
				for (VoFIterator vofit(irregCells,currEBGraph); vofit.ok(); ++ vofit)
				{
					const VolIndex& vof = vofit();
					const IntVect& gridIndex = vof.gridIndex();
					if ((gridIndex*16)/ncells == IntVect(D_DECL(5,3,6)) ||
              (gridIndex*16)/ncells == IntVect(D_DECL(4,1,1))
             )
					{
						const RealVect& bndryCentroid = currEBISBox.bndryCentroid(vof);
						double volFrac = currEBISBox.volFrac(vof);
						pout() << "(iphase, ivol) = (" << iphase << "," << ivol << ") @" << gridIndex << ", volFrac=" << volFrac << ", bndryCentroid=" << bndryCentroid << endl;
					}
				}
			}
		}
	}

	delete geometries[0];
	delete geometries[1];
}
