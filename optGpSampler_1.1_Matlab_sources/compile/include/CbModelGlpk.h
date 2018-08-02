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

#ifndef CBMODELGLPK_H
#define CBMODELGLPK_H

#include <glpk.h>
#include "CbModel.h"

extern "C" CbModel* createGlpkModel(const arma::sp_mat S, std::vector<double>* lb, std::vector<double>* ub, std::vector<double>* c, std::vector<double>* b);

class CbModelGlpk : public CbModel
{
	protected:
		
		glp_prob* lp;
	
		/** PROTECTED METHODS **/
		DllExport bool setSolverSpecificModel(std::vector<double>* b);
		DllExport CbModelGlpk& clone(const CbModelGlpk& source);			
		

	public:	

		// Default constructor, only used in combination with assignment operator
		DllExport CbModelGlpk();
		DllExport ~CbModelGlpk();
		DllExport boost::shared_ptr<CbModel> clone();
		DllExport void describe();

		// Preferred constructor
		DllExport CbModelGlpk(const arma::sp_mat S, std::vector<double>* lb, std::vector<double>* ub, std::vector<double>* c, std::vector<double>* b); 
		
		// Copy and assignment constructor
		DllExport CbModelGlpk(const CbModelGlpk& cbModelSource); 
		DllExport CbModelGlpk& operator= (const CbModelGlpk &CbModelSource);

		using  CbModel::getLb;
		using  CbModel::getUb;
		using  CbModel::getObjectiveCoefs;
		using  CbModel::getRhs;

		// Flux bound getters and setters
		DllExport void setLb(std::vector<double>* lb);
		DllExport void setUb(std::vector<double>* ub);
		DllExport void setLb(int rxnId, double lbValue);
		DllExport void setUb(int rxnId, double ubValue);
		DllExport void setBounds(std::vector<double>* lb, std::vector<double>* ub);
		DllExport void setBounds(int rxnId, double lb, double ub);
		DllExport double getLb(int rxnId);
		DllExport double getUb(int rxnId);
		DllExport double getObjectiveCoef(int rxnId);
		DllExport double getRhs(int rxnId);
			
		// Objective getters and setters
		DllExport void setObjectiveCoef(int rxnId, double objCoefficient);
		DllExport void setObjectiveCoefs(std::vector<double>* c);
		DllExport void getObjectiveCoefs(std::vector<double>* c);

		// Optimize the model and return the objective value (FBA)
		DllExport double optimize(std::string objectiveSense, std::vector<double>* fluxVector = NULL);

		// Propagate the flux bounds (FVA)
		DllExport bool propagate(bool resetObjectiveCoefs=true);
};

#endif //CBMODELGLPK_H