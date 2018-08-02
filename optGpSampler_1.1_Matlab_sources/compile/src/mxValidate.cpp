#include "mxValidate.h"

bool mxValidate::validReactionIdsInList(std::vector<int>* list, int nRxns)
{
	for(unsigned int i=0; i < list->size(); i++)
	{
		if(list->at(i) < 0 || list->at(i) >= nRxns)
		{
			printf("Warning rxnId=%d exceeds the number of %d total reactions \n", list->at(i)+1, nRxns);
			return false;
		}
	}

	return true;
}


bool mxValidate::cbModel(const mxArray* model, bool checkMetNames, bool checkRxnNames)
{
	// Check for required input fields
	mxArray* mtbS = mxGetField(model, 0, "S");
	if(mtbS == NULL)
		mexErrMsgTxt("model.S not found, please check for stoichiometric matrix \n");
	
	mxArray* mtbLb = mxGetField(model, 0, "lb"); 
	if(mtbLb == NULL)
		mexErrMsgTxt("model.lb not found, please check for lower bound vector \n");

	mxArray* mtbUb = mxGetField(model, 0, "ub");
	if(mtbUb == NULL)
		mexErrMsgTxt("model.ub not found, please check for upper bound vector \n");
	
	mxArray* mtbC= mxGetField(model, 0, "c");
	if(mtbC == NULL)	
		mexErrMsgTxt("model.c not found, please check for objective vector \n");
	
	mxArray* mtbB = mxGetField(model, 0, "b");
	if(mtbB == NULL)	
		mexErrMsgTxt("model.b not found, please check for right hand side vector \n");
	
	if(checkMetNames)
	{
		if(mxGetField(model, 0, "mets") == NULL)	
			mexErrMsgTxt("model.mets not found, please check for metabolite names \n");
	}

	if(checkRxnNames)
	{
		mxArray* mtbRxns = mxGetField(model, 0, "rxns");
		if(mtbRxns == NULL)	
			mexErrMsgTxt("model.rxns not found, please check for reaction names \n");
	}


	// Check length and orientation of vectors
	if(mxGetN(mtbLb) > 1 || mxGetN(mtbUb) > 1 || mxGetN(mtbC) > 1 || mxGetN(mtbB) > 1)
		mexErrMsgTxt("All vectors must be  column vectors \n");

	if(mxGetN(mtbS) != mxGetM(mtbLb))
		mexErrMsgTxt("Length lower bound vector  (model.lb) does not correspond to the number of columns in model.S \n");
	
	if(mxGetN(mtbS) != mxGetM(mtbUb))
		mexErrMsgTxt("Length upper bound vector  (model.ub) does not correspond to the number of columns in model.S \n");
	
	if(mxGetN(mtbS) != mxGetM(mtbC))
		mexErrMsgTxt("Length objective coefficients vector (model.c) does not correspond to the number of columns in model.S \n");
	
	if(mxGetM(mtbS) != mxGetM(mtbB))
		mexErrMsgTxt("Length right hand side vector (model.b) does not correspond to the number of rows in model.S \n");
	

	if(checkMetNames)
	{
		//mxArray* mtbMets = mxGetField(model, 0, "mets");
		//printf("nr mets == %d or nr mets = %d, %d \n", mxGetM(mtbMets), mxGetNumberOfElements(mtbMets), mxGetM(mxGetField(model, 0, "mets")));
		
		if(mxGetM(mtbS) != mxGetM(mxGetField(model, 0, "mets")))
		{
			mexErrMsgTxt("Length metabolite names vector (model.mets) does not correspond to the number of rows in model.S \n");
		}
	}

	if(checkRxnNames)
	{
		if(mxGetN(mtbS) != mxGetM(mxGetField(model, 0, "rxns")))
		{
			mexErrMsgTxt("Length reaction names vector (model.rxns) does not correspond to the number of columns in model.S \n");
		}
	}

	return true;

}


