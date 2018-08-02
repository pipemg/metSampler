/**
	Contains input and output functions for files in MATLAB. Contains many functions
	to convert matlab matrices to other formats (sparse to dense) or to ARMADILLO

	@author: WL Megchelenbrink
	@Last update: 28-03-2014
	--------
	@Changelog:
	--------
	28-03-2014:: File created

	--------
*/



#include "CbModelIOMatlab.h"

using namespace arma;
using namespace std;

#include "mex.h"

/** PUBLIC FUNCTIONS **/

/**
	Converts a sparse mxArray matrix to a dense mxArray matrix or vice versa
*/
void CbModelIOMatlab::convertMatrix(const mxArray* mtbMatrixIn, mxArray* mtbMatrixOut)
{
	if(mxIsSparse(mtbMatrixIn) && !mxIsSparse(mtbMatrixOut))
	{
		CbModelIOMatlab::convertSparseMatrixToDense(mtbMatrixIn, mtbMatrixOut);
	}
	else if(!mxIsSparse(mtbMatrixIn) && mxIsSparse(mtbMatrixOut))
	{
		CbModelIOMatlab::convertDenseMatrixToSparse(mtbMatrixIn, mtbMatrixOut);
	}
}
		


/**
	Converts a sparse or dense mxArray matrix to a sparse armadillo matrix
*/
void CbModelIOMatlab::convertMatrix(const mxArray* mtbMatrixIn, arma::sp_mat* armaMatrixOut)
{
	if(!mxIsSparse(mtbMatrixIn))
	{
		mxArray* tmp = mxCreateSparse(0, 0, 0, mxREAL);
		CbModelIOMatlab::convertDenseMatrixToSparse(mtbMatrixIn, tmp);
		CbModelIOMatlab::convertSparseMatrixMatlabToArmadillo(tmp, armaMatrixOut);
		mxDestroyArray(tmp);
	}
	else
	{
		CbModelIOMatlab::convertSparseMatrixMatlabToArmadillo(mtbMatrixIn, armaMatrixOut);
	}
}


/**
	Converts a sparse or dense mxArray matrix to a dense armadillo matrix
*/
void CbModelIOMatlab::convertMatrix(const mxArray* mtbMatrixIn, arma::mat* armaMatrixOut)
{
	if(mxIsSparse(mtbMatrixIn))
	{
		mxArray* tmp = mxCreateDoubleMatrix(0, 0, mxREAL);
		CbModelIOMatlab::convertSparseMatrixToDense(mtbMatrixIn, tmp);
		CbModelIOMatlab::convertDenseMatrixMatlabToArmadillo(tmp, armaMatrixOut);
		mxDestroyArray(tmp);
	}
	else
	{
		CbModelIOMatlab::convertDenseMatrixMatlabToArmadillo(mtbMatrixIn, armaMatrixOut);
	}
}


/**
 Converts an armadillo matrix to Matlab sparse or dense
*/
void CbModelIOMatlab::convertMatrix(const arma::mat* armaMatrixIn, mxArray* mtbMatrixOut)
{
	if(mxIsSparse(mtbMatrixOut))
	{
		sp_mat tmp; 
		CbModelIOMatlab::convertDenseMatrixToSparse(armaMatrixIn, &tmp);
		CbModelIOMatlab::convertSparseMatrixArmadilloToMatlab(&tmp, mtbMatrixOut);
	}
	else
	{
		CbModelIOMatlab::convertDenseMatrixArmadilloToMatlab(armaMatrixIn, mtbMatrixOut);
	}
}


/**
	Converts a sparse armadillo matrix to Matlab sparse or dense
*/
void CbModelIOMatlab::convertMatrix(const arma::sp_mat* armaMatrixIn, mxArray* mtbMatrixOut)
{
	if(mxIsSparse(mtbMatrixOut))
	{
		CbModelIOMatlab::convertSparseMatrixArmadilloToMatlab(armaMatrixIn, mtbMatrixOut);
	}
	else
	{
		arma::mat tmp; 
		CbModelIOMatlab::convertSparseMatrixToDense(armaMatrixIn, &tmp);
		CbModelIOMatlab::convertDenseMatrixArmadilloToMatlab(&tmp, mtbMatrixOut);
	}
}

/*
	Converts a sparse armadillo matrix to a dense matrix
*/
void CbModelIOMatlab::convertMatrix(const arma::sp_mat* armaMatrixIn, arma::mat* armaMatrixOut)
{
	CbModelIOMatlab::convertSparseMatrixToDense(armaMatrixIn, armaMatrixOut);
}


/*
	Converts a dense armadillo matrix to a sparse matrix
*/
void CbModelIOMatlab::convertMatrix(const arma::mat* armaMatrixIn, arma::sp_mat* armaMatrixOut) 
{
	CbModelIOMatlab::convertDenseMatrixToSparse(armaMatrixIn, armaMatrixOut);
}




/** PROTECTED FUNCTIONS **/

/**
	Converts a sparse Matlab matrix to a sparse Armadillo matrix 
	OK
*/
void CbModelIOMatlab::convertSparseMatrixMatlabToArmadillo(const mxArray* mtlbMatrix, sp_mat* armaMatrix)
{
	if(!mxIsSparse(mtlbMatrix))
	{
		mexErrMsgTxt("The supplied matrix is not in sparse format \n");
		return;
	}
	
	/*
	// Get Matlab sparse structure
	mwIndex* rows = mxGetIr(mtlbMatrix);
	mwIndex* cols = mxGetJc(mtlbMatrix);	
	double* data = mxGetPr(mtlbMatrix);
	*/

	// Set matrix size
	arma::mat tmp; //(nMets, nRxns);
	CbModelIOMatlab::convertMatrix(mtlbMatrix, &tmp);

	(*armaMatrix) = tmp;

//	printf("nnz :: %d \n", armaMatrix->n_nonzero);

	/*
	armaMatrix->set_size(nMets, nRxns);

	//double* vals = armaMatrix->values;
	const uword* ir = armaMatrix->row_indices;
	uword* jc = armaMatrix->col_ptrs;




	for(mwIndex i=0; i < nRxns; i++)	
	{
		for(mwIndex j=cols[i]; j < cols[i+1]; j++)	// loop over the non-zeros in each column
		{
			// rows[j] gives the row index for the jth stored non-zero value, data[j] is the corresponding value
			armaMatrix->at((const uword)rows[j], i) = data[j];
		}
	}
	*/
}


/**
	Converts a dense Matlab matrix to a dense Armadillo matrix 
	OKish
*/
void CbModelIOMatlab::convertDenseMatrixMatlabToArmadillo(const mxArray* mtlbMatrix, mat* armaMatrix)
{
	// Things from Matlab (aka source)
	int nRows = (int)mxGetM(mtlbMatrix);
	int nCols = (int)mxGetN(mtlbMatrix);
	double* srcPtr = mxGetPr(mtlbMatrix);

	if(mxIsSparse(mtlbMatrix))
	{
		mexErrMsgTxt("The supplied matrix is in sparse format, must be dense \n");
		return;
	}

	// Resize arma matrix and copy data
	armaMatrix->resize(nRows, nCols); 

	double* destPtr = armaMatrix->memptr();
	std::copy(srcPtr, srcPtr + nRows*nCols, destPtr);
}



/*
	Converts a sparse Armadillo matrix to a sparse Matlab matrix 
	OK
*/
void CbModelIOMatlab::convertSparseMatrixArmadilloToMatlab(const sp_mat* armaMatrix, mxArray* mlbMatrix)
{
	const uword nMets = armaMatrix->n_rows;
	const uword nRxns = armaMatrix->n_cols;
	const uword nnz = armaMatrix->n_nonzero;

	// Re-allocate space 
	mxSetM(mlbMatrix, nMets);
	mxSetN(mlbMatrix, nRxns);
	mxSetNzmax(mlbMatrix, nnz);
	
	// Allocate memory for row and column pointers	
	mwIndex* jc = (mwIndex*)mxCalloc(nRxns + 1, sizeof(mwIndex));
	mwIndex* ir = (mwIndex*)mxCalloc(nnz, sizeof(mwIndex));
	double* pr = (double*)mxCalloc(nnz, sizeof(double));
	
	// Copy values
	std::copy(armaMatrix->col_ptrs, armaMatrix->col_ptrs + nRxns+1, jc);
	std::copy(armaMatrix->row_indices, armaMatrix->row_indices + nnz, ir);
	std::copy(armaMatrix->values, armaMatrix->values + nnz, pr);

	// Store sparse info in sparse mxArray
	mxSetJc(mlbMatrix, jc);
	mxSetIr(mlbMatrix, ir);
	mxSetPr(mlbMatrix, pr);
}

/*
	Converts a dense Armadillo matrix to a dense Matlab matrix 
	OK
*/
void CbModelIOMatlab::convertDenseMatrixArmadilloToMatlab(const mat* armaMatrix, mxArray* mlbMatrix)
{
	// Get pointers
	const double* armaPtr = armaMatrix->memptr();
	double* mlbPtr = (double*)mxMalloc(armaMatrix->n_elem * sizeof(double));

	// Copy data
	std::copy(armaPtr, armaPtr + armaMatrix->n_elem, mlbPtr);
	mxSetM(mlbMatrix, armaMatrix->n_rows);
	mxSetN(mlbMatrix, armaMatrix->n_cols);
	mxSetPr(mlbMatrix, mlbPtr);
}



/**
	Converts a sparse armadillo matrix to a dense matrix
	OK
*/
void CbModelIOMatlab::convertSparseMatrixToDense(const arma::sp_mat* sparseMatrix, arma::mat* denseMatrix)
{
	(*denseMatrix) = (*sparseMatrix);
}


/**
	Converts a sparse matlab matrix to a dense matrix
	OK
*/
void CbModelIOMatlab::convertSparseMatrixToDense(const mxArray* sparseMatrix, mxArray* denseMatrix)
{
	// Resize and fill with zeros
	int nRows = mxGetM(sparseMatrix);
	int nCols = mxGetN(sparseMatrix);

	// Pointers to sparse ingoing matrix
	mwIndex* jc = mxGetJc(sparseMatrix);
	mwIndex* ir = mxGetIr(sparseMatrix);
	double* prSparse = mxGetPr(sparseMatrix);

	// Allocate memory
	double* prDense = (double*)mxCalloc(nRows * nCols, sizeof(double));

	mwIndex k=0;
	for(mwIndex i=0; i < (mwIndex)nCols; i++)
	{
		for(mwIndex j=jc[i]; j < jc[i+1]; j++)
		{
			// Replace non-zero entries
			prDense[nRows*i + ir[j]] = prSparse[k];
			k++;
		}
	}

	// Pointers to dense outgoing matrix
	mxSetM(denseMatrix, nRows);
	mxSetN(denseMatrix, nCols);
	mxSetPr(denseMatrix, prDense);
	

}


/*
	Converts a dense matlab matrix to a sparse matrix
	OK
*/
void CbModelIOMatlab::convertDenseMatrixToSparse(const mxArray* denseMatrix, mxArray* sparseMatrix)
{
	// Get the matrix dimensions
	int nRows = mxGetM(denseMatrix);
	int nCols = mxGetN(denseMatrix);

	// Get a pointer to the data
	double* prIn = mxGetPr(denseMatrix);

	// Allocate space and temporary vectors
	mwIndex* jcOut = (mwIndex*)mxCalloc(nCols+1, sizeof(mwIndex));
	vector<double> prTmp;
	vector<mwIndex> irTmp;
	
	// Get the non-zero values
	mwIndex nnz = 0;
	for(mwIndex i=0; i < (mwIndex)nCols; i++)
	{
		for(mwIndex j=0; j < (mwIndex)nRows; j++)
		{
			mwIndex index = i*nRows + j;
			if(prIn[index] !=0)
			{
				// Add only non-zero elements
				prTmp.push_back(prIn[index]); 
				irTmp.push_back(j);
				nnz++;
			}
		}

		// Update jc array
		jcOut[i+1] = nnz;
	}

	// Allocate space
	mwIndex* irOut = (mwIndex*)mxCalloc(nnz, sizeof(mwIndex));
	double* prOut = (double*)mxCalloc(nnz, sizeof(double) );
	
	// Copy
	std::copy(&irTmp[0], &irTmp[0] + nnz, irOut);
	std::copy(&prTmp[0], &prTmp[0] + nnz, prOut);

	// Set matrix dimensions and data
	mxSetM(sparseMatrix, nRows);
	mxSetN(sparseMatrix, nCols);
	mxSetJc(sparseMatrix, jcOut);
	mxSetIr(sparseMatrix, irOut);
	mxSetPr(sparseMatrix, prOut);
}


/*
	Converts a dense armadillo matrix to a sparse matrix
*/
void CbModelIOMatlab::convertDenseMatrixToSparse(const arma::mat* denseMatrix, arma::sp_mat* sparseMatrix)
{
	(*sparseMatrix) = (*denseMatrix);
}
