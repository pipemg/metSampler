/**
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

#include "CbModel.h"

using namespace std;
using namespace arma;

//#include "mex.h"

CbModel::CbModel(const string solverName) : solverName(solverName)
{

	nMets = 0;
	nRxns = 0;
};

CbModel::CbModel(const string solverName, sp_mat S) : solverName(solverName)
{
	
	this->S = S;
	this->nMets = S.n_rows;
	this->nRxns = S.n_cols;
	metNames.resize(S.n_rows);
	rxnNames.resize(S.n_cols);
}

CbModel::~CbModel()
{
}


arma::sp_mat* CbModel::getS()
{
	return &S;
}


boost::shared_ptr<vector<double> > CbModel::getLb()
{
	boost::shared_ptr<vector<double> > lb(new vector<double>(nRxns));

	for(int i=0; i < nRxns; i++)
	{
		lb->at(i) = getLb(i);
	}

	return lb;
}

boost::shared_ptr<vector<double> > CbModel::getUb()
{
	boost::shared_ptr<vector<double> > ub(new vector<double>(nRxns));

	for(int i=0; i < nRxns; i++)
	{
		ub->at(i) = getUb(i);
	}

	return ub;
}


boost::shared_ptr<vector<double> > CbModel::getObjectiveCoefs()
{
	boost::shared_ptr<vector<double> > c(new vector<double>(nRxns));

	for(int i=0; i < nRxns; i++)
	{
		c->at(i) = getObjectiveCoef(i);
	}

	return c;
}


boost::shared_ptr<vector<double> > CbModel::getRhs()
{
	boost::shared_ptr<vector<double> > b(new vector<double>(nMets));

	for(int i=0; i < nMets; i++)
	{
		b->at(i) = getRhs(i);
	}

	return b;
}


boost::shared_ptr<std::vector<string> > CbModel::getMetNames()
{
	boost::shared_ptr<std::vector<string> > _metNames(new vector<string>(nMets));
	for(int i=0; i < nMets; i++)
	{
		_metNames->at(i) = getMetName(i);
	}

	return _metNames;
}

boost::shared_ptr<std::vector<string> > CbModel::getRxnNames()
{
	boost::shared_ptr<std::vector<string> > _rxnNames(new vector<string>(nRxns));
	for(int i=0; i < nRxns; i++)
	{
		_rxnNames->at(i) = getRxnName(i);
	}

	return _rxnNames;
}




/**
	Returns the number of reactions in the stoichiometric matrix
*/
int CbModel::getNrRxns()
{
	return nRxns;
}

/**
	Returns the number of metabolites in the stoichiometric matrix
*/
int CbModel::getNrMets()
{
	return nMets;
}


void CbModel::printBounds()
{
	//printf(" === MODEL BOUNDS === \n");
	for(int i=0; i < nRxns; i++) 
	{
		//printf("%03d :: %s [%8.2f, %8.2f] \n", i+1, getRxnNameById(i), getLb(i), getUb(i));
	}
}


/**
	Sets the minimal growth rate for an FVA
*/
void CbModel::setMinimalGrowthRateConstraint(int rxnId, double growthPercentage, string objectiveSense)
{
	// Add the minimum required biomass output as a constraint
	setObjectiveCoef(rxnId, 1.0);
	double maxGrowthRate = optimize(objectiveSense);

	if(strcmp(objectiveSense.c_str(), "max") == 0)
	{
		//printf("setting lb[%d] = %f \n", rxnId, (growthPercentage / 100.0) * maxGrowthRate);
		setLb(rxnId,  (growthPercentage / 100.0) * maxGrowthRate);
	}
	else
	{
		setUb(rxnId,  (growthPercentage / 100.0) * maxGrowthRate);
	}
}


string CbModel::getRxnName(int rxnId)
{
	return rxnNames[rxnId];
}

string CbModel::getMetName(int metId)
{
	return metNames[metId];
}

void CbModel::setRxnName(int rxnId, std::string rxnName)
{
	rxnNames[rxnId] = rxnName;
}

void CbModel::setMetName(int metId, std::string metName)
{
	metNames[metId] = metName;
}

void CbModel::setRxnNames(vector<string>* rxnNames)
{
	this->rxnNames = (*rxnNames);
}

void CbModel::setMetNames(vector<string>* metNames)
{
	this->metNames = (*metNames);
}

int CbModel::getRxnIdByName(string rxnName)
{
	vector<string>::iterator it = find(rxnNames.begin(), rxnNames.end(), rxnName);
	if(it == rxnNames.end())
			return -1;
	
	return int(it - rxnNames.begin());
}

int CbModel::getMetIdByName(string metName)
{
	vector<string>::iterator it = find(metNames.begin(), metNames.end(), metName);
	if(it == metNames.end())
			return -1;
	
	return int(it - metNames.begin());
}













/*
int CbModel::setRxnIdAndName(string rxnName)
{
	std::map<string,int>::iterator it=rxnIdsAndNames.find(rxnName);
	return rxnIdsAndNames.;
}


void CbModel::setRxnIdAndName(string rxnName, int rxnId)
{
	 rxnIdsAndNames.insert (std::pair<int,string>(rxnId,rxnName));
}*/



		
/*
void CbModel::printModel()
{
	this->printAllModelConstraint();
	this->printBounds();
	this->printObjCoefficients();
	this->printMinMax();
}
		
// Print all model extractables
void CbModel::printAllModelConstraint()
{
	printf(" === MODEL CONSTRAINTS == \n");
	for(IloModel::Iterator it(model); it.ok(); ++it)
	{
		cout << (*it) << endl;
	}
}


		
void CbModel::printObjCoefficients()
{
	printf(" === MODEL OBJECTIVE COEFFICIENTS === \n");
	for (IloExpr::LinearIterator it = IloExpr(obj.getExpr()).getLinearIterator(); it.ok();++it)
	{
		printf("Env=%s, Id=%d,  coef=%4.2f, bounds=[%8.2f, %8.2f] \n",it.getVar().getEnv().getName(), it.getVar().getId(), it.getCoef(), it.getVar().getLb(), it.getVar().getUb()); 
	}
		
	printf("Environment for Cplex:: %s \n", cplex.getEnv().getName());
}

void CbModel::printMinMax()
{
	printf(" === MODEL MIN/MAX BIOMASS === \n");
	double fMin = this->optimize("min");
	double fMax = this->optimize("max");
	printf("Minimal biomass=%4.2f, Maximum biomass=%4.2f \n", fMin, fMax);

}
*/


		

//void CbModel::setAllZeros(IloNumArray* AllZeros)
//{
//	for(int i=0; i < nRxns; i++)
//	{
//		(*AllZeros)[i] = 0.0;
//	}
//}




