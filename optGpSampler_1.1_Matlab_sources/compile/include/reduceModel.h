/**
	HEADER FILE

	@author: WL Megchelenbrink
	@Last update: 28-03-2014
	--------
	@Changelog:
	--------
	28-03-2014:: File created

	--------
*/

#ifndef REDUCE_MODEL_H
#define REDUCE_MODEL_H

#include "CbModel.h"
#include "CbModelCreator.h"

//using namespace arma;
#define FLUX_BOUND_ERR 1e-7

class ReduceModel
{
	protected:
		
		std::string solverName;
		const double tol;

		boost::shared_ptr<CbModel> reductionStep(const boost::shared_ptr<CbModel>& model);
		boost::shared_ptr<CbModel> setFixedFluxes(const boost::shared_ptr<CbModel>& model, const double maxBiomass, std::vector<int>* fixedReactions, std::vector<double>* fixedFluxes);
		double getMaximumRange(const boost::shared_ptr<CbModel>& model);
		
	public:
		
		ReduceModel(std::string solverName, const double tol);
		boost::shared_ptr<CbModel> reduce(const boost::shared_ptr<CbModel>& model, std::vector<int>* fixedReactions, std::vector<double>* fixedFluxes);
		bool areModelsConsistent(const boost::shared_ptr<CbModel>& model, const boost::shared_ptr<CbModel>& reducedModel);

};

#endif // REDUCE_MODEL_H
