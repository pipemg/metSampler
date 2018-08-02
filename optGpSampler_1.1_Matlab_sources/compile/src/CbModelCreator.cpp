#include "CbModelCreator.h"

// Supported solvers
#include "CbModelCplex.h"
#include "CbModelGurobi.h"
#include "CbModelGlpk.h"

#if defined(_WIN32) || defined(_WIN64)
	#include <Windows.h>
	#include <delayimp.h>
#elif defined(__unix)
	#include <dlfcn.h>
#endif

using namespace std;
using namespace arma;

// Only for debug info
#include "mex.h"

int CbModelCreator::cnt = 0;
void* CbModelCreator::handle = NULL;

#if defined(_WIN32) || defined(_WIN64)
	CbModel* loadCbModel(const sp_mat &S, vector<double>* lb, vector<double>* ub, vector<double>* c, vector<double>* b, int solverType)
	{
		switch(solverType)
		{
			case 1:
				 return new CbModelCplex(S, lb, ub, c, b);
			break;

			case 2:
				return new CbModelGurobi(S, lb, ub, c, b);
			break;

			case 3:
				return new CbModelGlpk(S, lb, ub, c, b);
			break;
		}
	}
#endif


CbModelCreator::CbModelCreator(string solverName)
{
    handle = NULL;
	this->solverName = solverName;

	// Load libraries when needed (the UNIX way)
	#if defined(__unix)

		if(!handle)
		{
			solverName[0] = toupper(solverName[0]);
			std::string libName = "./linux_lib/libCbModel" + solverName + ".so";

			// Lazy load the library
			handle = dlopen(libName.c_str(), RTLD_LAZY);
		}

		if (!handle)
		{
			cerr << "Cannot load library " << dlerror() << endl;
			return;
		}
	#endif

	// Increase count
	cnt++;
}

CbModelCreator::~CbModelCreator()
{
	// Decrease count
	cnt--;
	closeDynamicLibrary();
}

CbModel* CbModelCreator::createCbModel(const sp_mat& S, vector<double>* lb, vector<double>* ub, vector<double>* c, vector<double>* b, int solverType)
{
	#if defined(_WIN32) || defined(_WIN64)
		__try
		{
			return loadCbModel(S, lb, ub, c, b, solverType);
		}
		__except(EXCEPTION_EXECUTE_HANDLER)
		{
			return NULL;
		}
    #endif

	return NULL;
}


boost::shared_ptr<CbModel> CbModelCreator::create(const arma::sp_mat S, vector<double>* lb, vector<double>* ub, vector<double>* c, vector<double>* b)
{
    // Initialize empty model
	boost::shared_ptr<CbModel> model;

	// Library already loaded by constructor
	// Now load the function to create a model
	#if defined(__unix)

		// If the library did not load, return empty model
		if(!handle)
			return model;

		dlerror();
		const char* dlsym_error;

		// Get access to the create function
		solverName[0] = toupper(solverName[0]);
		std::string symbolName = "create" + solverName + "Model";

		try
		{
			create_t* create_model = (create_t*) dlsym(handle, symbolName.c_str());
			dlsym_error = dlerror();


			if (dlsym_error)
			{
				cerr << "Cannot load function create: " <<  dlsym_error  << endl;
				return model;
			}

			// Create an instance of the class
			model.reset(create_model(S, lb, ub, c, b));
		}
		catch(...)
		{
			printf("Unexpected error while loading:: %s\n", symbolName.c_str());
			return model;
		}
	// Load libraries when needed (the WINDOWS way)
	#elif defined(_WIN32) || defined(_WIN64)
		try
		{
			 std::map<std::string,int> solvers;

			 solvers["cplex"]	= 1;
			 solvers["gurobi"]	= 2;
			 solvers["glpk"]	= 3;

			 model.reset(createCbModel(S, lb, ub, c, b, solvers[solverName]));
		}
		catch(IloException& e)
		{
			printf("Error creating CPLEX model:: %s  \n", e.getMessage());
		}
		catch(std::exception& e)
		{
			printf("Standard library error %s \n", e.what());
		}
		catch(...)
		{
			printf("Unexpected error during %s model creation \n", solverName.c_str());
		}

		if (!model)
		{
			printf("Failed to create %s model\n", solverName.c_str());
		}
	#endif

	
	return model;
}



void CbModelCreator::closeDynamicLibrary()
{
	if(handle && cnt <= 0)
	{
		#if defined(__unix)
			// Close the dynamically loaded library on UNIX systems
			if(strcasecmp(solverName.c_str(),"glpk") == 0)
			{
				free_t* freeGlpkEnvironment = (free_t*) dlsym(handle, "freeGlpkEnvironment");
				freeGlpkEnvironment();
			}		
			
			if(dlclose(handle) != 0)
			{
				printf("Error while closing %s library \n", solverName.c_str());
			}

		#elif defined(_WIN32) || defined(_WIN64)
			// Close the dynamically loaded library on WINDOWS systems
			char loadedDll[25];
			solverName[0] = toupper(solverName[0]);
			sprintf(loadedDll, "CbModel%s.dll", solverName.c_str());
			__FUnloadDelayLoadedDLL2(loadedDll);
		#endif

		// Reset count and handle; there are no active models or model constructors when the library closed
		handle = NULL;
		cnt = 0;
	}
}




