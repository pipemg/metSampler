/**
	HEADER FILE
	Constraint-based metabolic model class using the IBM ILOG Cplex solver. 
	
	@author: WL Megchelenbrink
	@Last update: 28-03-2014
	--------
	@Changelog:
	--------
	28-03-2014:: File created

	--------
*/
#ifndef CBMODELGUROBI_H
#define CBMODELGUROBI_H

#include "gurobi_c++.h"
#include "CbModel.h"

extern "C" CbModel* createGurobiModel(const arma::sp_mat S, std::vector<double>* lb, std::vector<double>* ub, std::vector<double>* c, std::vector<double>* b);

class CbModelGurobi : public CbModel
{
	protected:

		GRBEnv env;
		GRBModel model;
		GRBVar* rxns;
		
		DllExport bool setSolverSpecificModel(vector<double>* b);
		DllExport CbModelGurobi& clone(const CbModelGurobi& source);
		
	public:	
		
		// Default constructor, only used in combination with assignment operator
		DllExport CbModelGurobi();
		DllExport ~CbModelGurobi();
		DllExport boost::shared_ptr<CbModel> clone();
		DllExport void describe();

		// Preferred constructor
		DllExport CbModelGurobi(const arma::sp_mat S, vector<double>* lb, vector<double>* ub, vector<double>* c, vector<double>* b); 
		
		// Copy and assignment constructor
		DllExport CbModelGurobi(const CbModelGurobi& cbModelSource); 
		DllExport CbModelGurobi& operator= (const CbModelGurobi &CbModelSource);
		
		using  CbModel::getLb;
		using  CbModel::getUb;
		using  CbModel::getObjectiveCoefs;
		using  CbModel::getRhs;



		// Flux bound getters and setters
		DllExport void setLb(vector<double>* lb);
		DllExport void setUb(vector<double>* ub);
		DllExport void setLb(int rxnId, double lbValue);
		DllExport void setUb(int rxnId, double ubValue);
		DllExport void setBounds(vector<double>* lb, vector<double>* ub);
		DllExport void setBounds(int rxnId, double lb, double ub);
		DllExport double getLb(int rxnId);
		DllExport double getUb(int rxnId);
		DllExport void updateGRBModel();

		// Objective getters and setters
		DllExport void setObjectiveCoef(int rxnId, double objCoefficient);
		DllExport void setObjectiveCoefs(vector<double>* c);
		DllExport void getObjectiveCoefs(vector<double>* c);
		DllExport double getObjectiveCoef(int rxnId);
		DllExport double getRhs(int rxnId);


		// Optimize the model and return the objective value (FBA)
		DllExport double optimize(string objectiveSense, vector<double>* fluxVector = NULL);

		// Propagate the flux bounds (FVA)
		DllExport bool propagate(bool resetObjectiveCoefs=true);
};

#endif // CBMODELGUROBI_H


	
