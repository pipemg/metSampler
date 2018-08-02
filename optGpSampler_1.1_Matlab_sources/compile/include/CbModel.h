/**
	HEADER FILE
	Abstract constraint-based metabolic model class. 
	This class is extended by the solver specific classes:
	CbModelCplex, CbModelGurobi and CbModelGlpk

	@author: WL Megchelenbrink
	@Last update: 28-03-2014
	--------
	@Changelog:
	--------
	28-03-2014:: File created

	--------
*/

#ifndef CBMODEL_H
#define CBMODEL_H

#if defined __unix__
	#define DllExport
#else   
	#define DllExport __declspec(dllexport)
#endif

//#include <omp.h>
#include <armadillo>

#include <vector>
#include <map>
#include <memory>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

//using namespace arma;
#define FLUX_BOUND_ERR 1e-7


class CbModel
{
	protected:
	
		arma::sp_mat S;
		int nMets, nRxns;
		const std::string solverName;

		std::vector<std::string> metNames;
		std::vector<std::string> rxnNames;
		
		virtual bool setSolverSpecificModel(std::vector<double>* b) = 0;
	
		//virtual CbModel* clone();

	//	virtual CbModel(const CbModel& cbModelSource) = 0; 
		//virtual CbModel* clone(const CbModel& source)  = 0;
		
		//CbModelCplex& CbModelCplex::clone(const CbModelCplex& source);
		
		

	public:
		int getNrRxns();	
		int getNrMets();
	//	void setRxnNames(vector<string>* rxnNames);	//@todo implement
	//	void setMetNames(vector<string>* metNames); //@todo implement

		/* *************** VIRTUAL FUNCTIONS BELOW *************** */
		// Flux bound setters

		virtual boost::shared_ptr<CbModel> clone() = 0;
		virtual void describe() = 0;


		CbModel(const std::string solverName);
		CbModel(const std::string solverName, arma::sp_mat S);
		
		virtual ~CbModel() = 0;

		virtual void setBounds(int rxnId, double lb, double ub) = 0;
		virtual void setBounds(std::vector<double>* lb, std::vector<double>* ub) = 0; 
		virtual void setLb(std::vector<double>* lb) = 0;
		virtual void setUb(std::vector<double>* ub) = 0;
		virtual void setLb(int rxnId, double lb) = 0;
		virtual void setUb(int rxnId, double ub) = 0;

		// Flux bound getters
		virtual double getLb(int rxnId) = 0;
		virtual double getUb(int rxnId) = 0;

		boost::shared_ptr<std::vector<double> > getLb();
		boost::shared_ptr<std::vector<double> > getUb();
		boost::shared_ptr<std::vector<double> > getObjectiveCoefs();
		boost::shared_ptr<std::vector<double> > getRhs();
		boost::shared_ptr<std::vector<std::string> > getMetNames();
		boost::shared_ptr<std::vector<std::string> > getRxnNames();
	
		// Objective coefficients setters
		virtual	void setObjectiveCoef(int rxnId, double objCoefficient) = 0;
		virtual void setObjectiveCoefs(std::vector<double>* c) = 0;
	
		virtual double getObjectiveCoef(int rxnId) = 0;
		virtual double getRhs(int rxnId) = 0;


		virtual double optimize(std::string objectiveSense, std::vector<double>* fluxVector = NULL) = 0;
		virtual bool propagate(bool resetObjectiveCoefs=true) = 0;
		
		
		void setMinimalGrowthRateConstraint(int rnxId, double growthRate, std::string sense);
		void printBounds();
		
		std::string getRxnName(int rxnId);
		std::string getMetName(int metId);

		void setRxnName(int rxnId, std::string rxnName);
		void setMetName(int rxnId, std::string metName);
		void setRxnNames(std::vector<std::string>* rxnNames);
		void setMetNames(std::vector<std::string>* metNames);

		int getRxnIdByName(std::string rxnName);
		int getMetIdByName(std::string metName);
//		
		arma::sp_mat* getS();



/** TODO, these go into a separate class	
		void printModel();
		void printAllModelConstraint();
		
		void printObjCoefficients();
		void printMinMax();
	*/
};

//class CbModelGurobi : public CbModel;
//class CbModelGlkp : public CbModel;

#endif // CBMODEL_H
