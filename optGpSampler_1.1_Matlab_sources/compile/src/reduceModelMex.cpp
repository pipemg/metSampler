/**
	Flux Measurement Prioritization (FMP)
	HEADER FILE
	@author: WL Megchelenbrink
	@Last update: 03-04-2014
	--------
	@Changelog:
	--------
	03-04-2014:: File created
	--------
*/

#include "reduceModel.h"
#include "mxValidate.h"
#include "CbModelIOMatlab.h"

#if defined(_WIN32) || defined(_WIN64)
    #include <windows.h>
    #include <delayimp.h>
#endif

using namespace arma;
using namespace std;



bool validateReduceModelInput(const mxArray* prhs[])
{
	// Check for valid model
	mxValidate::cbModel(prhs[0], true, true);
	
	char solverName[20];
	mxGetString(prhs[1], solverName, 20);
	

	if(strcmp(solverName,"cplex") != 0 && strcmp(solverName, "gurobi") != 0 && strcmp(solverName, "glpk") != 0)
		mexErrMsgTxt("Invalid solver, choose from {'cplex', 'gurobi', 'glpk'} \n");


	double tolerance = (double)mxGetScalar(prhs[2]);					
	if(tolerance < 0)
		mexErrMsgTxt("Tolerance can not be smaller than 0 \n");

	
	return true;
}



void readMxCellArray(mxArray* cellArray, vector<string>* str)
{
	if(!mxIsCell(cellArray))
		mexErrMsgTxt("Expected cell array input\n");
	
	int nCells = mxGetNumberOfElements(cellArray);
	
	str->resize(nCells);

	for(int i=0; i < nCells; i++)
	{
		mxArray* currentCell = mxGetCell(cellArray, i);
		str->at(i) = mxArrayToString(currentCell);
	}
	
}


void writeMxCellArray(vector<string>* str, mxArray* cellArray)
{
	if(!mxIsCell(cellArray))
		mexErrMsgTxt("Expected cell array output \n");
	
	int nCells = mxGetNumberOfElements(cellArray);
	for(int i=0; i < nCells; i++)
	{
		mxArray* myString = mxCreateString(str->at(i).c_str());
		mxSetCell(cellArray,i, myString);
	}
}





/**
	This is the standard API to Matlab
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{ 
	/////////////
	const int nRequiredInputs = 3;
	if(nrhs != nRequiredInputs)
	{
		mexErrMsgTxt("Invalid number of input arguments \n");
	}

	// Check for valid inputs
	mxArray* mtbS = mxGetField(prhs[0], 0, "S");
	validateReduceModelInput(prhs);
	
	// Read the cbModel
	mxArray* mtbLb = mxGetField(prhs[0], 0, "lb"); 
	mxArray* mtbUb = mxGetField(prhs[0], 0, "ub");
	mxArray* mtbC= mxGetField(prhs[0], 0, "c");
	mxArray* mtbB = mxGetField(prhs[0], 0, "b");
	mxArray* mtbMets = mxGetField(prhs[0], 0, "mets");
	mxArray* mtbRxns = mxGetField(prhs[0], 0, "rxns");


	size_t nMets = mxGetM(mtbS);
	size_t nRxns = mxGetN(mtbS);
	
	sp_mat S;
	CbModelIOMatlab::convertMatrix(mtbS, &S);
	vector<double> lb(mxGetPr(mtbLb), mxGetPr(mtbLb) + nRxns);							// Lower bounds
	vector<double> ub(mxGetPr(mtbUb), mxGetPr(mtbUb) + nRxns);							// Upper bounds
	vector<double> c(mxGetPr(mtbC), mxGetPr(mtbC) + nRxns);								// Objective coefficients
	vector<double> b(mxGetPr(mtbB), mxGetPr(mtbB) + nMets);								// Righthand sides
	vector<string> mets;																// Metabolite names
	vector<string> rxns;																// Reaction names
	
	
	readMxCellArray(mtbMets, &mets);
	readMxCellArray(mtbRxns, &rxns);

	// Parameters
	char solverName[20];
	mxGetString(prhs[1], solverName, 20);
	
	double tolerance = (double)mxGetScalar(prhs[2]);

	
	// Make a lp solver specific model
    CbModelCreator cm(solverName);
    boost::shared_ptr<CbModel> model = cm.create(S, &lb, &ub, &c, &b);

	if(!model)
	{
		char solverErrMsg[100];
		sprintf(solverErrMsg, "Model creation failed, are the %s libraries present? \n", solverName);
		mexErrMsgTxt(solverErrMsg);
	}

	
	/*printf("max biomass before anything ::: %f \n", model->optimize("max"));
	boost::shared_ptr<CbModel> model2 = model->clone();

	printf("max biomass of model clone ::: %f \n", model2->optimize("max"));

	model2->propagate();
	printf("max biomass of model prop ::: %f \n", model2->optimize("max"));

	return;*/
	model->setMetNames(&mets);
	model->setRxnNames(&rxns);

	////////// MODEL REDUCTION ////////////
	vector<int> fixedReactions;		// contains the ids of fixed reaction
	vector<double> fixedFluxes;		// contains the flux (the value) of the fixed reaction
	
	boost::shared_ptr<CbModel> reducedModel;
	bool consistentModels = false;
	try
	{
		// Make original model consistent
		model->propagate();
	
		ReduceModel rm(solverName, tolerance);
		
	/*	for(int i=0; i < model->getNrRxns(); i++)
		{
			if(model->getObjectiveCoef(i) > 0)
				printf("i=%d, objective = %f \n",i, model->getObjectiveCoef(i));
		}*/

		reducedModel = rm.reduce(model, &fixedReactions, &fixedFluxes);
		if(reducedModel)
		{	
			consistentModels = rm.areModelsConsistent(model, reducedModel);
		}
	}
	catch(std::exception& e)
	{
		printf("Error :: %s", e.what());
		mexErrMsgTxt("ReduceModel terminated \n");
	}
	catch(...)
	{
		mexErrMsgTxt("Unexpected exception, ReduceModel terminated \n");
	}

	if(!reducedModel)
	{
		mexErrMsgTxt("Reduced model does not exist \n");
	}

	if(!consistentModels)
	{
		printf("WARNING:: Model objectives are not consistent (max fba orginal = %.4f, max fba reduced model = %.4f)  ..\n", model->optimize("max"), reducedModel->optimize("max"));
	}
	
	////////// RETURN REDUCED MODEL TO MATLAB ////////////
	nMets = reducedModel->getNrMets();
	nRxns = reducedModel->getNrRxns();

	int nzMax = reducedModel->getS()->n_nonzero;

	mxArray* St = mxCreateSparse(0, 0, 0, mxREAL);
	CbModelIOMatlab::convertMatrix(reducedModel->getS(), St);
	reducedModel->getS();
	
	// Create the ouput vectors
	mxArray* lbOut = mxCreateDoubleMatrix(nRxns, 1, mxREAL);
	mxArray* ubOut = mxCreateDoubleMatrix(nRxns, 1, mxREAL);
	mxArray* cOut = mxCreateDoubleMatrix(nRxns, 1, mxREAL);
	mxArray* bOut = mxCreateDoubleMatrix(nMets, 1, mxREAL);
	mxArray* metsOut = mxCreateCellMatrix(nMets, 1);
	mxArray* rxnsOut = mxCreateCellMatrix(nRxns, 1);

	

	// Output metabolite and reaction names
	writeMxCellArray(reducedModel->getMetNames().get(), metsOut);
	writeMxCellArray(reducedModel->getRxnNames().get(), rxnsOut);
	
	// Get pointers to the data
	double* pLb = mxGetPr(lbOut);
	double* pUb = mxGetPr(ubOut);
	double* pC = mxGetPr(cOut);
	double* pB = mxGetPr(bOut);

	sp_mat* Stmp = reducedModel->getS();

	// Fill the data
	for(int i=0; i < nRxns; i++)
	{
		pLb[i] = reducedModel->getLb(i);
		pUb[i] = reducedModel->getUb(i);
		pC[i] = reducedModel->getObjectiveCoef(i);
	}

	for(int i=0; i < nMets; i++)
	{
		pB[i] = reducedModel->getRhs(i);
	}
	
	
	// Put it all together into a Matlab structure
	const int nFields = 7;
	const char* fieldNames[nFields];
	fieldNames[0] = "S";
	fieldNames[1] = "lb";
	fieldNames[2] = "ub";
	fieldNames[3] = "c";
	fieldNames[4] = "b";
	fieldNames[5] = "mets";
	fieldNames[6] = "rxns";
	
	mxArray * modelOut = mxCreateStructMatrix(1, 1, nFields, fieldNames);
	mxSetField(modelOut, 0, "S", St);
	mxSetField(modelOut, 0, "lb", lbOut);
	mxSetField(modelOut, 0, "ub", ubOut);
	mxSetField(modelOut, 0, "c", cOut);
	mxSetField(modelOut, 0, "b", bOut);
	mxSetField(modelOut, 0, "mets", metsOut);
	mxSetField(modelOut, 0, "rxns", rxnsOut);


	int nFixedRxns = fixedReactions.size();

	mxArray* mxFixedRxns = mxCreateNumericMatrix(nFixedRxns, 1, mxINT32_CLASS, mxREAL);
	mxArray* mxFixedFluxes = mxCreateDoubleMatrix(nFixedRxns, 1, mxREAL);
	int* pFixedRxns = (int*)mxGetData(mxFixedRxns);
	double* pFixedFluxes = mxGetPr(mxFixedFluxes);
	
	std::copy(&fixedReactions[0], &fixedReactions[0] + nFixedRxns, pFixedRxns);
	std::copy(&fixedFluxes[0], &fixedFluxes[0] + nFixedRxns, pFixedFluxes);

	// Reduced model
	plhs[0] = modelOut;
	plhs[1] = mxFixedRxns;
	plhs[2] = mxFixedFluxes;
}

