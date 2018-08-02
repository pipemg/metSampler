/* 
	File:			optGpSampler.cpp
	Author:			Wout Megchelenbrink, Radboud University Nijmegen
	Purpose:		Artificial Centering Hit and Run Sampling (ACHR)
	Summary:		Uses parallel ACHR to sample a metabolic network
	Last update:	30/04/2014
	Version:		1.1
	Changes:
	
	/***
	//////// 2014 //////////
	30/04/ CHANGELOG
	- fixed a bug that could cause samples to go out of the nullspace
	- added support for MATLAB mwblas and mwlapack for fast matrix multiplication and SVD
	- added SVD in C++ to get the nullspace
	- added interrupt handler to use control-c termination in Matlab
	- added custom reduceModel
	- added sampling summary
	- added verbose / non-verbose mode
	- added lazy loading of GLPK, GUROBI or CPLEX libraries to allow for easier LP solver plug-in
	- added support to use previously collected warmup points
	- added a warning for reactions that can not carry flux ("dead reactions")	
	- added possibility to only do a warm-up (0 sample points)
	- removed the dependency on ACML libraries (instead uses mwblas and mwlapack that ship with MATLAB)
	- added robustness test results for 11 publicly available (genome-scale) constraint-based metabolic models
	---------------

	//////// 2013 //////////
	03-06-2013
	- Fixed DGEMM for ACML

	18-04-2013:
	 - small bug fix

	04-03-2013
	- can compile against MKL BLAS that ships with Matlab
	

	//////// 2014 //////////
	10-10
	 - added Python bindings
	
	03-08
	- Returns only samples for the list of reactions we're interested in

	02-07
	- Added support for IBM Ilog cplex
	- Uses Mersenne Twister random number generator from Boost C++
	- Made more effective by using less warmups
	- Does a random walk of length X, and then resamples from an initial point
*/

#include "optGpSampler.h"

#ifdef __cplusplus 
    extern "C" bool utIsInterruptPending();
#else
    extern bool utIsInterruptPending();
#endif

#include "blas.h"
#include "lapack.h"
#include "mex.h"

using namespace std;
using namespace arma;

boost::mt19937 random_gen;
boost::variate_generator<boost::mt19937&,boost::uniform_01<> > random_real(random_gen,boost::uniform_01<>());

// Rediculous but UNIX and WINDOW do not have the same function to compare  case-insensitive strings .. 
#if !defined(__unix)
	int strcasecmp(const char *f1, const char *f2)
	{
		return _strcmpi(f1, f2);
	}
#endif

//extern void  mkl_set_num_threads(int number);

/* Construct */
optGpSampler::optGpSampler(const boost::shared_ptr<CbModel>& model, const int nSamples, const int nSteps, const int nThreads, bool verbose)
	: tol(1e-7), model(model), nRxns(model->getNrRxns()), lb(nRxns), ub(nRxns), nSamples(nSamples), nSteps(nSteps), nWarmups(2 * model->getNrRxns()), nThreads(nThreads), verbose(verbose)
{
	warmupTime = 0;
	samplingTime = 0;

	nStepsBeforeProjection = 25;

	// Reduce the model and compute the nullspace of the reduced model
	setNullSpace(model->getS());
	
	// Allocate space for the warmups and samples
	nSamplesPerThread = (int)ceil((double)nSamples / nThreads);		// divide the points over the threads
	warmupPts.resize(nRxns, nWarmups);	
	samplePts.resize(nThreads);

	for(int i=0; i < nThreads; i++)
	{
		samplePts[i].resize(nRxns, (nWarmups + nSamplesPerThread));
	}

	// Seed the random number generator
	random_gen.seed(static_cast<int>(time(0)));

	// Set omp threads to 1 to put BLAS in single thread mode
	//mkl_set_num_threads(1);
	omp_set_num_threads(1);
	
}

// Destruct
optGpSampler::~optGpSampler()
{
	// Reset the available parallel procs
	omp_set_num_threads(omp_get_num_procs());
}


/**
	Use previously sampled warmup points, make sure model are bounds are consistent!
*/
void optGpSampler::setWarmupPts(double* warmupPts)
{
	std::copy(warmupPts, warmupPts + nRxns*nRxns*2, this->warmupPts.memptr());

	lb.resize(nRxns);
	ub.resize(nRxns);

	for(int i=0; i < nRxns; i++)
	{
		lb[i] = model->getLb(i);
		ub[i] = model->getUb(i);
	}
}



/**
	Singular value decomposition by LAPACK
*/
void optGpSampler::svd(mat* U, vec* S, mat* Vt, mat* A)
{
	// Compute singular value decomposition in LAPACK
	ptrdiff_t nRows = model->getNrMets();
	ptrdiff_t nCols = nRxns;
	ptrdiff_t LWORK = min(nRows, nCols) * (6 + 4 * min(nRows, nCols)) + max(nRows, nCols);

	// SVD outputs
	S->resize(min((int)nRows, (int)nCols));
	U->resize((int)nRows, (int)nRows);
	Vt->resize((int)nCols, (int)nCols);

	vector<ptrdiff_t> IWORK(8*min(nRows, nCols));
	ptrdiff_t INFO;
	vector<double> WORK(LWORK);

	// Lapack Singular Value Decomposition (SVD)
	dgesdd("A", &nRows, &nCols, A->memptr(), &nRows, S->memptr(), U->memptr(), &nRows, Vt->memptr(), &nCols, &WORK[0], &LWORK, &IWORK[0], &INFO);
}


/**
	Set the nullspace of the stoichiometric matrix
*/
void optGpSampler::setNullSpace(sp_mat* S)
{
	// Singluar values smaller than svd_tol are considered zero
	const double svd_tol = 1e-12;
	
	mat A((*S));
	
	// SVD outputs
	mat U;
	vec s;
	mat Vt;

	// Do the svd
	svd(&U, &s, &Vt, &A);
	
	uvec nonZeroSingularRows = find(s > svd_tol);
	int rank = nonZeroSingularRows.size();

	// Nullspace; the columns of V corresponding to the zero singular values
	N =  Vt.rows(rank, nRxns-1).t();
	Nt = N.t();
}


/**
	Creates the warmup points by running 2n FBA's
*/
void optGpSampler::setWarmupPointsAndFluxBounds(int biomassRxnId)
{
	double startTime = omp_get_wtime();
	
	vector<double> flux(nRxns);
	vec armaFlux(nRxns);

	for(int i=0; i < nRxns; i++)
	{
		// Change objective to this reaction's flux
		model->setObjectiveCoef(i, 1.0);

		for(int s=0; s <= 1; s++)
		{
			if(s == 0)
			{
				lb[i] = model->optimize("min", &flux);
			}
			else
			{
				ub[i] = model->optimize("max", &flux);
			}
			
			std::copy(flux.begin(), flux.end(), armaFlux.memptr());

			warmupPts.col(2*i + s) = armaFlux;
		}

		if(verbose)
		{
			if(ub[i] - lb[i] < tol && abs(lb[i]) < tol)
			{
				printf("Reaction %d can not carry flux: model reduction recommended for efficiency \n", i);
			}
		}
	}

	if(biomassRxnId > -1)
		model->setObjectiveCoef(biomassRxnId, 1.0);

	double endTime = omp_get_wtime();
	warmupTime = endTime - startTime;
}


//void optGpSampler::printErr(vec* A, int iterNr, double bnd_tolerance)
//{
//	vec prev = (*A);
//
//	if(verbose)
//		{
//			printf("ITERATION = %d \n", iterNr);
//			try
//			{
//				// LB VIOL
//				uvec v1 = find( (*A) < lb - bnd_tolerance);
//				if(v1.size() > 0)
//				{
//					uvec x1 = find( abs( (*A)(v1))  > tol);
//				
//					if(x1.size() > 0)
//					{
//						uvec id = v1(x1);
//						for(int j=0; j < x1.size(); j++)
//						{
//							printf("LB VIOL AT :: j=%d, id=%d, rxnId=%d, FLUX=%.8f (%.8f), lb=%.8f, ub=%.8f \n", j, x1[j], id[j], A->at(id[j]), prev.at(id[j]), lb[id[j]], ub[id[j]]);
//						}
//					}
//				}
//
//
//				// UB VIOL
//				uvec v2 = find( (*A) > ub + bnd_tolerance);
//				if(v2.size() > 0)
//				{
//					uvec x2 = find( abs( (*A)(v2))  > tol);
//				
//					if(x2.size() > 0)
//					{
//						uvec id = v2(x2);
//						for(int j=0; j < x2.size(); j++)
//						{
//							printf("UB VIOL AT ::  j=%d, id=%d, rxnId=%d, FLUX=%.8f (%.8f), lb=%.8f, ub=%.8f \n", j, x2[j], id[j], A->at(id[j]), prev.at(id[j]), lb[id[j]], ub[id[j]]);
//						}
//					}
//				}
//
//
//			}
//			catch(std::exception& e)
//			{
//				printf("not ok ::: %s \n", e.what());
//				return;
//
//			}
//
//			verbose = false;
//		}
//}


/**
	Project a sampled vector onto the nullspace of the stoichiometry matrix
*/
bool optGpSampler::projectOntoNullspace(vec* A, int iterNr)
{
	// Do matrix multiplication with BLAS DGEMM
	vec tmp(N.n_cols, 1);
	matrixMultiply(&tmp, &Nt, A);
	matrixMultiply(A, &N, &tmp);
	
	// Check if the lower- or upper bounds were violated
	const double bnd_tolerance = 1e-4; 
	int i = 0;
	bool accept = true;

	while(any((*A) < lb - bnd_tolerance) || any((*A) > ub + bnd_tolerance))
	{
		// Fix bound violations
		uvec lbViol = find( (*A) < lb);
		if(lbViol.size()  > 0)
			(*A)(lbViol) = lb(lbViol);

		uvec ubViol = find( (*A) > ub);
		if(ubViol.size()  > 0)	
			(*A)(ubViol) = ub(ubViol);

		// Do matrix multiplication with BLAS DGEMM
		vec tmp(N.n_cols, 1);
		matrixMultiply(&tmp, &Nt, A);
		matrixMultiply(A, &N, &tmp);

		if(i > 10) 
		{
			nStepsBeforeProjection = max(25, nStepsBeforeProjection - 100);
			accept = false;
			break;
		}

		i++;
	}
	
	if(i == 0)	
		nStepsBeforeProjection = min(25, nStepsBeforeProjection + 25);
	
	return accept;
}



/*
 * Exit codes :: 0=uninitialized, 1=success, 2=too many samples discarded
*/
bool optGpSampler::sample(int* exitcode)
{
	
	if(nSamples == 0)
	{
		(*exitcode) = 1;
		return true;
	}
	
	uvec validSampling(nThreads);
	validSampling.ones();

	try
	{
		double startTime = omp_get_wtime();
		nPointsDiscarded = 0;
		
		// Move warmups to the center
		vec cp = mean(warmupPts,1);
		warmupPts = warmupPts * 0.33 + 0.67 * cp * ones(1, nWarmups);

		// Add the warmup points temporarily to the samples points
		for(int i=0; i < nThreads; i++)
		{
			samplePts[i].cols(0, nWarmups-1) = warmupPts;
		}
		
		mat* pSamplePts;
		vec curPoint(nRxns);
	    
		int tid;
		int r, i, k;
		bool inNullSpace = false;
		
		vector<int> discarded(nThreads, 0);
		
		#pragma omp parallel private(tid, r, i, k, inNullSpace, pSamplePts) firstprivate(cp, curPoint)  num_threads(nThreads)
		{
			// Get current thread id	
			tid = omp_get_thread_num();

			pSamplePts =& samplePts[tid];
			
			i = nWarmups;
			while(i < nWarmups + nSamplesPerThread)
			{
				// Issue a warning if many samples are discarded (probably near zero values in S)
				if(discarded[tid] >= 100 && discarded[tid] % 100 == 0 && tid == 0)
					printf("WARNING :: %d samples discarded (i=%d) \n");
			

				if(!inNullSpace)
				{
					// Pick different point
					#pragma omp critical
					{
						r = static_cast<int>(floor((random_real() * (i-1) ))); 
						curPoint = pSamplePts->col(r); 
						inNullSpace = true;
					}
				}

				
				// Check for control-C in Matlab
				if(utIsInterruptPending())
				{
					i = nWarmups + nSamplesPerThread;
					(*exitcode) = 3;
					//printf("Operation terminated by user during optGpSampler::sample \n");
					validSampling[tid] = 0;
					continue;
				}
			

				for(k=1; k <= nSteps; k++)
				{
					// Get the next point, in case of error, there is probably something wrong in the matrix S
					if(!getPoint(&curPoint, i, &cp, pSamplePts))
					{
						discarded[tid]++;
						k = nSteps;
						inNullSpace = false;
						continue;
					}
					
					// Correct floating point errors by projecting onto N
					if(k % nStepsBeforeProjection == 0 && inNullSpace)
					{
						if(!projectOntoNullspace(&curPoint))		// discard the point if it is outside the nullspace
						{
							discarded[tid]++;
							inNullSpace = false;
							continue;
						}
					}
				}
				
				
				if(!inNullSpace) 
				{
					continue;
				}

				// Final reprojection
				if(nSteps % nStepsBeforeProjection !=0)
				{
  					if(!projectOntoNullspace(&curPoint) || !inNullSpace)
					{
						discarded[tid]++;
						inNullSpace = false;
						continue;
					}
				}
	
				// Update centerpoint and store current point
				if(i % 25 == 0)
					projectOntoNullspace(&cp);

				cp = (cp*i + curPoint) / (i+1);
				pSamplePts->col(i) = curPoint;
				
				if(i == nWarmups + nSamplesPerThread -1)
					validSampling[tid] = 1;
				
				i++;

			}
		}

		
		for(int i=0; i < nThreads; i++)
		{
			nPointsDiscarded+= discarded[i];
		}

		double endTime = omp_get_wtime();
		samplingTime = endTime - startTime;
		
	}
	catch(std::exception &e)
	{
		printf("Exception caught :: %s \n", e.what());
		return false;
	
	}
	catch(...)
	{
		printf("Unexpected error caught in function createSamplePts. \n");
		return false;
	}

	// Set the exit code to successful sampling
	if((*exitcode) == 0)
		(*exitcode) = 1;

	return arma::any(validSampling==0) ? false : true;
}



bool optGpSampler::getPoint(vec* currentPoint, int i, vec* cp, mat* samplePts, int iterNr)
{
	try
	{
		/***  STEP 2:: PICK A RANDOM DIRECTION  ****/
	
		// Random nr between 0 and i
		int r = 0;
		double randNr = 0.0;
		#pragma omp critical
		{
			r = static_cast<int>(floor((random_real() * (i-1) ))); // take a direction to a random point .. not necessarily a random warmup
			randNr = random_real();
		}
	
		// Get a direction from the center point to the warmup point
		vec u = (samplePts->col(r) - (*cp));
		u = u / norm(u,2);

		// Check if we can move forward or backward
		uvec posDir = find(u > tol);
		uvec negDir = find(u < -tol);
		uvec movable = arma::join_cols(posDir, negDir); 
		
		/** STEP 3:: GET MAXIMUM STEP SIZE **/
		// Get distance to the lower and upperbounds for the movable directions
		vec minStepVec = arma::join_cols( -((*currentPoint)(posDir) - lb(posDir)) / u(posDir), (ub(negDir) - (*currentPoint)(negDir) ) / u(negDir));
		vec maxStepVec = arma::join_cols( (ub(posDir) - (*currentPoint)(posDir)) / u(posDir), -((*currentPoint)(negDir) - lb(negDir)) / u(negDir));
		
		if(minStepVec.size() == 0 || maxStepVec.size() == 0)
			return false;
		
		double maxStep = std::max(as_scalar(arma::min(maxStepVec)), 0.0);
		double minStep = std::min(as_scalar(arma::max(minStepVec)), 0.0);
		double stepSize = as_scalar(minStep + randNr *(maxStep - minStep));
		
		//// STEP 4:: Make the step (stepsize * direction)
		(*currentPoint)(movable) += stepSize * u(movable);
	}
	catch(std::exception& e)
	{
		printf("Exception caught:: %s \n", e.what());
		return false;
	}

	return true;
}


double* optGpSampler::getWarmupPoints()
{
	return warmupPts.memptr();
}

/** 
	Check in what order this is returned
*/
double* optGpSampler::getSamplePoints(int threadId)
{
	return samplePts[threadId].memptr() + nRxns*nRxns*2;
}


// Perform matrix multiplication using the single threaded dgemm blas functions from ACML or MKL
void optGpSampler::matrixMultiply(vec* result, mat* A, vec* B)
{
    // Transpose or Not
    char transp = 'N';

    // Scalar values to use in dgemm
    double alpha = 1.0;
    double beta = 0.0;
            
	// Dimensions of input matrices
	ptrdiff_t m = A->n_rows;
	ptrdiff_t k = A->n_cols;
	ptrdiff_t n = B->n_cols;
    
	// This version of DGEMM from ACML
	dgemm(&transp, &transp, &m, &n, &k, &alpha, A->memptr(), &m, B->memptr(), &k, &beta, result->memptr(), &m); 
}


/* Get random objective vector */
void optGpSampler::setRandomVector(double* vec, int nElems)
{
	for(int i=0; i<nElems; i++)
	{
		vec[i] = random_real() - 0.5;
	}
}


void optGpSampler::printSummary()
{
	printf(" ========================= SAMPLING SUMMARY =======================\n");
	printf("Model contains %d reactions and %d metabolites \n", model->getNrRxns(),  model->getNrMets());
	printf("Warmup phase took %4.2f seconds \n", warmupTime);
	printf("Sampling phase took %4.2f seconds \n", samplingTime);
	printf("Nr steps before projection on the nullspace: minimum=25, maximum=%d\n", nStepsBeforeProjection);
	printf("%d points outside the nullspace have been discarded and were resampled \n", nPointsDiscarded);
	printf("%d sample points successfully sampled with step size=%d \n", nSamples, nSteps);
	printf(" ==================================================================\n");
}









	
