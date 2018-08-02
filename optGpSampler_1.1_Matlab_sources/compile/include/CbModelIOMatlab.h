/**
	HEADER FILE
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


//#include "CbModelIO.h"
#include "mex.h"
#include <armadillo>

class CbModelIOMatlab
{	
	protected:
		static void convertDenseMatrixMatlabToArmadillo(const mxArray* mlbMatrix, arma::mat* armMatrix);
		static void convertSparseMatrixMatlabToArmadillo(const mxArray* mlbMatrix, arma::sp_mat* armMatrix);

		static void convertDenseMatrixArmadilloToMatlab(const arma::mat* armMatrix, mxArray* mlbMatrix);
		static void convertSparseMatrixArmadilloToMatlab(const arma::sp_mat* armMatrix, mxArray* mlbMatrix);

		static void convertSparseMatrixToDense(const arma::sp_mat* sparseMatrix, arma::mat* denseMatrix);
		static void convertSparseMatrixToDense(const mxArray* sparseMatrix, mxArray* denseMatrix);
		static void convertDenseMatrixToSparse(const arma::mat* denseMatrix, arma::sp_mat* sparseMatrix);
		static void convertDenseMatrixToSparse(const mxArray* denseMatrix, mxArray* sparseMatrix);
	

	public:
		static void convertMatrix(const mxArray* mtbMatrixIn, mxArray* mtbMatrixOut);
		static void convertMatrix(const mxArray* mtbMatrixIn, arma::mat* armaMatrixOut);
		static void convertMatrix(const mxArray* mtbMatrixIn, arma::sp_mat* armaMatrixOut);

		static void convertMatrix(const arma::mat* armabMatrixIn, mxArray* mtbMatrixOut);
		static void convertMatrix(const arma::sp_mat* armabMatrixIn, mxArray* mtbMatrixOut);
		static void convertMatrix(const arma::sp_mat* armabMatrixIn, arma::mat* armaMatrixOut);
		static void convertMatrix(const arma::mat* armabMatrixIn, arma::sp_mat* armaMatrixOut);

};
