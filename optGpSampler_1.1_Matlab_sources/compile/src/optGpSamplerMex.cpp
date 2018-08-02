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
#include "optGpSampler.h"
#include "mxValidate.h"
//#include "CbModelIO.h"

#if defined(_WIN32) || defined(_WIN64)
    #include <windows.h>
    #include <delayimp.h>
#endif

#include "mex.h"

using namespace arma;
using namespace std;

bool validateOptGpSamplerInput(const mxArray* prhs[], size_t nRxns)
{
	// Check for valid model
	mxValidate::cbModel(prhs[0]);
	
	int nSamples = (int)mxGetScalar(prhs[2]);					
	int nSteps = (int)mxGetScalar(prhs[3]);						
	int nThreads = (int)mxGetScalar(prhs[4]);
	
	char solverName[20];
	mxGetString(prhs[5], solverName, 20);
	
	if(nSamples < 0)
		mexErrMsgTxt("Invalid number of sample points \n");

	if(nSteps < 0)
		mexErrMsgTxt("Invalid number of steps \n");

	if(nThreads < 0 || nThreads > omp_get_num_procs())
	{
		char threadErr[100];
		sprintf(threadErr,"Invalid number of threads (requested = %d, max_threads=%d) \n", nThreads, omp_get_max_threads());
		mexErrMsgTxt(threadErr);
	}

	if(strcasecmp(solverName,"cplex") != 0 && strcasecmp(solverName, "gurobi") != 0 && strcasecmp(solverName, "glpk") != 0)
		mexErrMsgTxt("Invalid solver, choose from {'cplex', 'gurobi', 'glpk'} \n");

	return true;
}




//	This is the standard API to Matlab
/*
int main() 
{ 
	
	sp_mat S;	
	vector<double> lb;
	vector<double> ub;
	vector<double> c;
	vector<double> b;

	string modelName = "Yeast_iND750";
	//string path = "..\\..\\..\\testModels\\csv\\Ecoli_iJO1366\\";
	string path = "D:\\Dropbox\\PHD\\FastCOBRA\\testModels\\csv\\Yeast_iND750\\";
	string solverName = "cplex";


	// Todo later:: improve robustness add more error handling
	CbModelIO cio(modelName, path);
	if(!cio.getStoichiometry(&S))
	{
		std::cerr << "Error reading stoichiometric matrix (S) file " << endl;
		return 0;
	}

	if(!cio.getLb(&lb))
	{
		std::cerr << "Error reading lower bound (lb) file " << endl;
		return 0;
	}
	
	if(!cio.getUb(&ub))
	{
		std::cerr << "Error reading upper bound (ub) file " << endl;
		return 0;
	}
	
	if(!cio.getObjectiveCoefficients(&c))
	{
		std::cerr << "Error reading objective coefficients (c) file " << endl;
		return 0;
	
	}
	
	if(!cio.getRhs(&b))
	{
		std::cerr << "Error reading right hand side (b) file " << endl;
		return 0;
	}


	int nSamples = 10;
	int nSteps = 10;
	int nThreads = 1;
	
	// Make a lp solver specific model
	CbModel* model = CbModelCreator::create(S, &lb, &ub, &c, &b, solverName);
	if(model == NULL)
	{
		char solverErrMsg[100];
		sprintf(solverErrMsg, "Model creation failed, are the %s libraries present? \n", solverName);
		mexErrMsgTxt(solverErrMsg);
	}

	//mat N = zeros(1266, 282);

	CbModel* reducedModel = optGpSampler::reduceModel(model, solverName);
	if (!reducedModel) {
		printf("Failed to create reduced model\n");
	}
	

	//////// SAMPLING ////////////
	optGpSampler ogp(reducedModel, nSamples, nSteps, nThreads);
	time_t start, end;
	time(&start);
	double diff;

	ogp.setWarmupPointsAndFluxBounds();
	time(&end);
	printf ("Warmups took %4.2f seconds\n", difftime (end, start));
	

	printf("so far so good \n");
	
	time(&start);
	if(!ogp.sample())
	{
		mexErrMsgTxt("Error during sampling phase \n");
	};
	time(&end);
	printf ("Sampling took %4.2f seconds\n", difftime(end, start));
	
	return 1;
}	
*/	




/**
	This is the standard API to Matlab
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{ 
	const int nRequiredInputs = 7;
	
	if(nrhs != nRequiredInputs)
	{
		mexErrMsgTxt("Invalid number of input arguments \n");
		return;
	}

	// Check for a valid model
	mxArray* mtbS = mxGetField(prhs[0], 0, "S");
	validateOptGpSamplerInput(prhs,  mxGetN(mtbS));
	
	// Read the cbModel
	mxArray* mtbLb = mxGetField(prhs[0], 0, "lb"); 
	mxArray* mtbUb = mxGetField(prhs[0], 0, "ub");
	mxArray* mtbC= mxGetField(prhs[0], 0, "c");
	mxArray* mtbB = mxGetField(prhs[0], 0, "b");
	
	size_t nMets = mxGetM(mtbS);
	size_t nRxns = mxGetN(mtbS);
	
	sp_mat S;
	CbModelIOMatlab::convertMatrix(mtbS, &S);
	vector<double> lb(mxGetPr(mtbLb), mxGetPr(mtbLb) + nRxns);			// Lower bounds
	vector<double> ub(mxGetPr(mtbUb), mxGetPr(mtbUb) + nRxns);			// Upper bounds
	vector<double> c(mxGetPr(mtbC), mxGetPr(mtbC) + nRxns);				// Objective coefficients
	vector<double> b(mxGetPr(mtbB), mxGetPr(mtbB) + nMets);				// Righthand sides

	// Parameters
	char solverName[20];
	
	double* warmupPts = mxGetPr(prhs[1]);					
	int nSamples = (int)mxGetScalar(prhs[2]);					
	int nSteps = (int)mxGetScalar(prhs[3]);						
	int nThreads = (int)mxGetScalar(prhs[4]);
	mxGetString(prhs[5], solverName, 20);
	bool verbose = (bool)mxGetScalar(prhs[6]);


	// Make a lp solver specific model
    CbModelCreator cm(solverName);
	boost::shared_ptr<CbModel> model = cm.create(S, &lb, &ub, &c, &b);
	
	if(!model)
	{
		char solverErrMsg[100];
		sprintf(solverErrMsg, "Model creation failed, are the %s libraries present? \n", solverName);
		mexErrMsgTxt(solverErrMsg);
	}
	
	int biomassRxnId = 0;
	vector<double>::iterator it = find(c.begin(), c.end(), 1);
	if(it == c.end())
		biomassRxnId = -1;
	else
		biomassRxnId = it - c.begin();
	
	//////// SAMPLING ////////////
	optGpSampler ogp(model, nSamples, nSteps, nThreads, verbose);

	// Warmup phase
	if(mxGetM(prhs[1]) != 0 || mxGetM(prhs[1]) != 0)
	{
		if(mxGetM(prhs[1]) != model->getNrRxns() || mxGetN(prhs[1]) != model->getNrRxns() * 2)
		{
			char txt[100];
			sprintf(txt, "Dimensions of warmup points must correspond to reduced model. Dimensions should be [%d, %d] \n", model->getNrRxns() , 2*model->getNrRxns() * 2);
			mexErrMsgTxt(txt);
		}
		else
		{
			if(verbose)
				printf("Warmup points present, skipping warmup phase \n");
		
			ogp.setWarmupPts(warmupPts);
		}
	}
	else
	{
		ogp.setWarmupPointsAndFluxBounds(biomassRxnId);
	}

	
	// Sampling phase
	bool samplingSucceeded = false;
	int exitcode = 0;
	try
	{
		samplingSucceeded = ogp.sample(&exitcode);
	}
	catch(std::exception& e)
	{
		printf("STD error %s \n", e.what());
	}
	catch(...)
	{
		printf("Unknown error caught \n");
	}

	if(!samplingSucceeded)
	{
		switch(exitcode)
		{
			default:
				char errMsg[100];
				sprintf(errMsg, "Sampling terminated with exit code [%d]\n", exitcode);
				mexErrMsgTxt(errMsg);
			break;

			case 2:
				mexErrMsgTxt("Too many sample points discarded due to nullspace or bound violations\n"
								"Please run optReduceModel with 'fixMatrix=1' \n");
			break;

			case 3:
				mexErrMsgTxt("Sampling terminated by user \n");
			break;

		}

	}

	
	//////// RETURN STUFF TO MATLAB ////////////
	// Copy the warmup points and samples to Matlab
	mxArray* mxWarmups = mxCreateDoubleMatrix(nRxns, nRxns*2, mxREAL);
	mxArray* mxSamples = mxCreateDoubleMatrix(nRxns, nSamples, mxREAL);
		
	const double* warmups = ogp.getWarmupPoints();

	double* mxWarmupsPtr = mxGetPr(mxWarmups);
	double* mxSamplesPtr = mxGetPr(mxSamples);

	std::copy(warmups, warmups + nRxns*nRxns*2, mxWarmupsPtr);
	
	int nSamplesPerThread = (int)ceil( (double)nSamples / nThreads);

	for(int i=0; i < nThreads; i++)
	{
		const double* samples = ogp.getSamplePoints(i);
		int p = nSamplesPerThread * model->getNrRxns();
	
		int nSamplesToCopy = 0;
		if(i < nThreads - 1)
		{
			nSamplesToCopy = model->getNrRxns() * nSamplesPerThread;
		}
		else
		{
			int nSamplesRequired = nSamples - nSamplesPerThread * (nThreads-1);
			nSamplesToCopy = model->getNrRxns() * nSamplesRequired;
		}
		
		std::copy(samples, samples + nSamplesToCopy, mxSamplesPtr + p*i);
	}


	plhs[0] = mxWarmups;
	plhs[1] = mxSamples;
	
	if(verbose)
		ogp.printSummary();
}

