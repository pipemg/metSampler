/**
	HEADER FILE
	Constraint-based metabolic model class using the IBM ILOG Cplex solver. 
	
	@author: WL Megchelenbrink
	@Last update: 17-04-2014
	--------
	@Changelog:
	--------
	28-03-2014:: File created

	--------
*/

#ifndef CBMODELCPLEX_H
#define CBMODELCPLEX_H

#define IL_STD

#include <ilcplex/ilocplex.h>
#include <ilcplex/cplex.h>
#include "CbModel.h"


extern "C" CbModel* createCplexModel(const arma::sp_mat S, std::vector<double>* lb, std::vector<double>* ub, std::vector<double>* c, std::vector<double>* b);

ILOSTLBEGIN


class CbModelCplex : public CbModel
{
	protected:
		
		IloEnv env;
		IloModel model;
		IloNumVarArray rxns;					
		IloObjective obj;
		IloCplex cplex;

		std::map<int, std::string> rxnIdsAndCplexNames;
		
		
		DllExport bool setSolverSpecificModel(vector<double>* b);
		DllExport CbModelCplex& clone(const CbModelCplex& source);	// Old way of cloning
		


		DllExport int getRxnIdByCplexName(std::string cplexRxnName) const;
		DllExport void setRxnIdAndCplexName(int rxnId, std::string cplexRxnName);

	public:	

		DllExport boost::shared_ptr<CbModel> clone(); 
		DllExport void describe();

		// Default constructor, only used in combination with assignment operator
		DllExport CbModelCplex();
		DllExport ~CbModelCplex();

		// Preferred constructor
		DllExport CbModelCplex(const arma::sp_mat S, vector<double>* lb, vector<double>* ub, vector<double>* c, vector<double>* b); 
		
		// Copy and assignment constructor
		DllExport CbModelCplex(const CbModelCplex& cbModelSource); 
		DllExport CbModelCplex& operator= (const CbModelCplex &CbModelSource);

		///////////////////////
//		DllExport boost::shared_ptr<std::vector<double>> getLb();
//		DllExport boost::shared_ptr<std::vector<double>> getUb();
//		DllExport boost::shared_ptr<std::vector<double>> getObjectiveCoefs();
//		DllExport boost::shared_ptr<std::vector<double>> getRhs();
//		DllExport boost::shared_ptr<std::vector<std::string>> getMetNames();
//		DllExport boost::shared_ptr<std::vector<std::string>> getRxnNames();
		using  CbModel::getLb;
		using  CbModel::getUb;
		using  CbModel::getObjectiveCoefs;
		using  CbModel::getRhs;
		
		//using  boost::shared_ptr<std::vector<double>> CbModel::getUb();

		// Flux bound getters and setters
		DllExport void setLb(vector<double>* lb);
		DllExport void setUb(vector<double>* ub);
		DllExport void setLb(int rxnId, double lbValue);
		DllExport void setUb(int rxnId, double ubValue);
		DllExport void setBounds(vector<double>* lb, vector<double>* ub);
		DllExport void setBounds(int rxnId, double lb, double ub);
		DllExport double getLb(int rxnId);
		DllExport double getUb(int rxnId);

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


#endif // CBMODELCPLEX_H
