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

#include "mex.h"
#include <vector>

class mxValidate  
{
	protected:

	public:
		static bool validReactionIdsInList(std::vector<int>* list, int nRxns);
		static bool cbModel(const mxArray* model, bool checkMetNames=false,  bool checkRxnNames=false);


};
