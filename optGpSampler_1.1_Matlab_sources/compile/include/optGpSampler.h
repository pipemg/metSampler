#include "CbModelCreator.h"

/** Boost random number libraries */
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>

#include <omp.h>
#include "CbModelIOMatlab.h"



extern boost::variate_generator<boost::mt19937&,boost::uniform_01<> > random_real;

#if !defined(__unix)
	int strcasecmp(const char *f1, const char *f2);
#endif


class optGpSampler
{

	protected:

		const boost::shared_ptr<CbModel>& model;
		const int nSamples, nSteps, nThreads;
		int nRxns, 	nWarmups, nPointsDiscarded, nSamplesPerThread; 
		const double tol;
		bool verbose;

		int nStepsBeforeProjection;
		arma::vec lb, ub;
		arma::mat N, Nt, warmupPts;
		
		std::vector<arma::mat> samplePts;
		double warmupTime, samplingTime;
	

		///////// *** PROTECTED METHODS *** ////////////////
		void setNullSpace(arma::sp_mat* S);
		void svd(arma::mat* U, arma::vec* S, arma::mat* V,  arma::mat* A);
	
		void setRandomVector(double* vec, int nElems);
		bool projectOntoNullspace(arma::vec* A, int iterNr=0);
		bool getPoint(arma::vec* currentPoint, int i, arma::vec* centerPoint, arma::mat* samplePt, int iterationNr=0);
		void matrixMultiply(arma::vec* result, arma::mat* A, arma::vec* B);
		

	public:
	
		void printErr(arma::vec* A, int iterNr, double bnd_tolerance);

		///////// *** PUBLIC METHODS *** ////////////////
		optGpSampler(const boost::shared_ptr<CbModel>& model, const int nSamples, const int nSteps, const int nThreads, bool verbose=true);
		~optGpSampler();

		void setWarmupPts(double* warmupPts);
		void setWarmupPointsAndFluxBounds(int biomassRxnId);
		void printSummary(); 

		// Sample return an exit code
		bool sample(int* exitcode);

		double* getWarmupPoints();
		double* getSamplePoints(int threadId); 
};


