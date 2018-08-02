#include "reduceModel.h"
#include <armadillo>

using namespace arma;
using namespace std;

#include "mex.h"


#ifdef __cplusplus 
    extern "C" bool utIsInterruptPending();
#else
    extern bool utIsInterruptPending();
#endif

ReduceModel::ReduceModel(string solverName, const double tol) : solverName(solverName), tol(tol) {}




//double ReduceModel::getMaximumRange(const boost::shared_ptr<CbModel>& model)
//{
//	double maxRange = 0;
//	for(int i=0;  i < model->getNrRxns(); i++)
//	{
//		maxRange = max(maxRange, model->getUb(i) - model->getLb(i));
//		//maxFlux = max(maxFlux, abs(model->getUb(i)));
//	}
//
//	return maxRange;
//}


boost::shared_ptr<CbModel> ReduceModel::reduce(const boost::shared_ptr<CbModel>& model, vector<int>* fixedReactions, vector<double>* fixedFluxes)
{
	
	bool maxReduction = false;

	boost::shared_ptr<CbModel> prevModel = model->clone();
	boost::shared_ptr<CbModel> reducedModel;
	
	//printf(" --------------------- INIT STEP ----------------------------- \n");
	//printf("MAX FBA orginal model=%f \n", model->optimize("max"));
	//printf("MAX FBA orginal prevModel=%f \n", prevModel->optimize("max"));
	//printf(" --------------------- ******* ----------------------------- \n\n");

	int i = 0;
	while(!maxReduction)
	{
		if(utIsInterruptPending())
			return model;

		// Now find dead reactions and metabolites
		reducedModel = reductionStep(prevModel);
	
		// Propagate constraints
		reducedModel->propagate();

//		printf("Reduction iteration %d:: Nr rxns previous model=%d, nRxns current model = %d, (%d) \n", i, prevModel->getNrRxns(), reducedModel->getNrRxns(), model->getNrRxns());
//		printf("   >>>  MAX FBA prevModel=%f, reducedModel=%f \n", prevModel->optimize("max"), reducedModel->optimize("max"));


		// Model reduction converged
		if(reducedModel->getNrRxns() == prevModel->getNrRxns())
			maxReduction = true;

		// Update
		prevModel = reducedModel->clone();
		i++;
	}

	// Set fixed reactions
	reducedModel = setFixedFluxes(reducedModel, reducedModel->optimize("max"), fixedReactions, fixedFluxes);
	reducedModel->propagate();
	
	// @TODO:: sort fixed reactions and fluxes


	return reducedModel;

}

/*
 Reduce the model by removing all fixed reaction and associated metabolites
*/
boost::shared_ptr<CbModel> ReduceModel::reductionStep(const boost::shared_ptr<CbModel>& model)
{
	// Get original stoichiometric matrix
	mat S(*model->getS());

	vector<double> lb;
	vector<double> ub;
	vector<double> c;
	vector<double> b;
	vector<string> mets;
	vector<string> rxns;
	
	vector<u32> keepDir;

//	double maxRange = getMaximumRange(model);
	//double relativeTol = tol * maxRange;

//	printf("ok maximum range = %4.4f, tol=%f, tolerance=%f \n", maxRange, tol, relativeTol);
	
	for(int i=0; i < model->getNrRxns(); i++)
	{
		// Check if the flux range is very small
	//	double relativeRange = (model->getUb(i) - model->getLb(i)) / maxRange;
		
		//if(i > 670 && i < 675)
		//{
		//	printf("i=%d, lb=%f, ub=%f, relRange=%f \n", i, model->getLb(i), model->getUb(i), relativeRange);
		//}

		if(abs(model->getLb(i)) < tol)
			model->setLb(i, 0.0);

		if(abs(model->getUb(i)) < tol)
			model->setUb(i, 0.0);


		
		if(model->getUb(i) - model->getLb(i) < tol && abs(model->getUb(i)) < tol)
		{
			// Near zero flux, make zero
		///	printf("setting bounds to zero for i=%d [%f, %f] \n", i, model->getLb(i), model->getUb(i));
		
			
		//	model->setLb(i, 0.0);
		//	model->setUb(i, 0.0);
		}	
		
		if(abs(model->getLb(i)) >= tol || abs(model->getUb(i)) >= tol)
		{
			// Add non-zero fluxes to the reduced model
			keepDir.push_back(i);
			lb.push_back(model->getLb(i));
			ub.push_back(model->getUb(i));
			c.push_back(model->getObjectiveCoef(i));
			rxns.push_back(model->getRxnName(i));
		}
	}

	
	S = S.cols(uvec(keepDir));

	uvec activeMetabolites = sum(abs(S), 1) != 0;
	for(int i=0; i < S.n_rows; i++) 
	{
		if(activeMetabolites[i]) 
		{
			b.push_back(model->getRhs(i));
			mets.push_back(model->getMetName(i));
		}
	}
	
	
	// Keep only non-orphanized metabolites
	S = S.rows(find(activeMetabolites));

	// Create a new CbModel
	CbModelCreator mc(solverName);
	boost::shared_ptr<CbModel> reducedModel = mc.create(sp_mat(S), &lb, &ub, &c, &b);
	
	reducedModel->setRxnNames(&rxns);
	reducedModel->setMetNames(&mets);
	
	return reducedModel;
}



boost::shared_ptr<CbModel> ReduceModel::setFixedFluxes(const boost::shared_ptr<CbModel>& model, const double maxBiomass, vector<int>* fixedReactions, vector<double>* fixedFluxes)
{
	// Only store fixed reactions that are not zero
	for(int i=0; i < model->getNrRxns(); i++)
	{
		if(model->getUb(i) - model->getLb(i) < tol &&  abs(model->getUb(i)) > tol)
		{
			double oldLb = model->getLb(i);
			model->setLb(i, model->getUb(i));
			
			if(model->optimize("max") - tol > maxBiomass || model->optimize("max") + tol < maxBiomass)
			{
				// The model becomes different if we make this flux fixed; roll back the change
				model->setLb(i, oldLb);
			}
			else
			{
				fixedReactions->push_back(i); 

				// Usually lb and ub will be the same, except for large tolerance, so we take the average
				fixedFluxes->push_back( (model->getUb(i) + model->getLb(i)) / 2 ); 
			}
		}
	}

	return model;
}


bool ReduceModel::areModelsConsistent(const boost::shared_ptr<CbModel>& model, const boost::shared_ptr<CbModel>& reducedModel)
{
	if( abs(model->optimize("max") - reducedModel->optimize("max")) > 1e-4 || abs(model->optimize("min") - reducedModel->optimize("min")) > 1e-4)
	{
		printf("Min/max objective original model ::[%f, %f] \n", model->optimize("min"), model->optimize("max"));
		printf("Min/max objective reduced model :: [%f, %f] \n", reducedModel->optimize("min"), reducedModel->optimize("max"));
		return false;
	}
	else
	{
		printf("---- Models are consistent ----  \n");
		printf("Min/max objective original model ::[%f, %f] \n", model->optimize("min"), model->optimize("max"));
		printf("Min/max objective reduced model :: [%f, %f] \n", reducedModel->optimize("min"), reducedModel->optimize("max"));
	
	}

	return true;
}


