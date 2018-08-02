#include "CbModel.h"
#include <boost/smart_ptr.hpp>

class CbModelCreator
{
	protected:
		typedef CbModel* create_t(const arma::sp_mat S, std::vector<double>* lb, std::vector<double>* ub, std::vector<double>* c, std::vector<double>* b);
		typedef void free_t();
		typedef void destroy_t(CbModel*);

		std::string solverName;

		CbModel* createCbModel(const arma::sp_mat& S, std::vector<double>* lb, std::vector<double>* ub, std::vector<double>* c, std::vector<double>* b, int solverType);
		void closeDynamicLibrary();
		static void *handle;
		static int cnt;

	public:
		CbModelCreator(std::string solverName);
		~CbModelCreator();
		boost::shared_ptr<CbModel> create(const arma::sp_mat S, std::vector<double>* lb, std::vector<double>* ub, std::vector<double>* c, std::vector<double>* b);
};


