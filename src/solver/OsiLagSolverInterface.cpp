#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cstdlib>
#include <numeric>
#include <cassert>
#include <cmath>

#include "OsiLagSolverInterface.h"



/**************************************************************************************************/
/*                                     Private helper methods        	    		              */			      
/**************************************************************************************************/

void OsiLagSolverInterface::updateRowMatrix_() const{
    if(!rowMatrixCurrent_){
        rowMatrix_.reverseOrderedCopyOf(colMatrix_);
        rowMatrixCurrent_ = true;
    }
}

void OsiLagSolverInterface::updateColMatrix_() const{
    if(!colMatrixCurrent_) {
        colMatrix_.reverseOrderedCopyOf(rowMatrix_);
        colMatrixCurrent_ = true;
    }   
}

//-----------------------------------------------------------------------------

void OsiLagSolverInterface::checkData_() const{
    int i;
    for(i = getNumRows() - 1; i >= 0; --i){
        if (rowlower_[i] > -1.0e20 && rowupper_[i] < 1.0e20 && rowlower_[i] != rowupper_[i]){
           //throw CoinError("Volume algorithm is unable to handle ranged rows","checkData_", "OsiVolSolverInterface");
        }
    }

    for(i = getNumCols() - 1; i >= 0; --i){
        if (collower_[i] < -1.0e20 || colupper_[i] > 1.0e20) 
            throw CoinError("Volume algorithm is unable to handle infinite bounds","checkData_", "OsiVolSolverInterface");
    }
}

//-----------------------------------------------------------------------------

/* It only works if in the flowForm.ccp we also change the signal of the <= constraints,
so they will be >= constraints in a minimization problem */
/* We need to do that because the AbstractLagrangian does that to compute the dual variables. */
void OsiLagSolverInterface::compute_rc_(const double* u, double* rc) const {
 
    rowMatrix_.transposeTimes(u, rc);

    const int psize = getNumCols();
    std::transform(rc, rc+psize, objcoeffs_, rc, std::minus<double>());
    std::transform(rc, rc+psize, rc, std::negate<double>());
}

//-----------------------------------------------------------------------------

void OsiLagSolverInterface::gutsOfDestructor_(){
    rowMatrix_.clear();
    colMatrix_.clear();
    rowMatrixCurrent_ = true;
    colMatrixCurrent_ = true;

    delete[] colupper_;   colupper_ = 0;  
    delete[] collower_;	collower_ = 0;  
    delete[] continuous_; continuous_ = 0;
    delete[] rowupper_;	rowupper_ = 0; 
    delete[] rowlower_;	rowlower_ = 0; 
    delete[] rowsense_;	rowsense_ = 0; 
    delete[] rhs_;	rhs_ = 0;      
    delete[] rowrange_;	rowrange_ = 0; 
    delete[] objcoeffs_;	objcoeffs_ = 0;

    delete[] colsol_;	        colsol_ = 0;   
    delete[] rowprice_;	        rowprice_ = 0;
    delete[] rowpriceHotStart_;	rowpriceHotStart_ = 0;
    delete[] rc_;       	        rc_ = 0;
    delete[] lhs_;       	        lhs_ = 0;

    lagrangeanCost_ = 0.0;

    maxNumrows_ = 0;
    maxNumcols_ = 0;
}

//-----------------------------------------------------------------------------

void OsiLagSolverInterface::rowRimAllocator_(){
    rowupper_ = new double[maxNumrows_];
    rowlower_ = new double[maxNumrows_];
    rowsense_ = new char[maxNumrows_];
    rhs_      = new double[maxNumrows_];
    rowrange_ = new double[maxNumrows_];
    rowprice_ = new double[maxNumrows_];
    lhs_      = new double[maxNumrows_];
}

//-----------------------------------------------------------------------------

void OsiLagSolverInterface::colRimAllocator_(){
    colupper_  = new double[maxNumcols_];
    collower_  = new double[maxNumcols_];
    continuous_ = new bool[maxNumcols_];
    objcoeffs_ = new double[maxNumcols_];
    colsol_    = new double[maxNumcols_];
    rc_        = new double[maxNumcols_];
}

//-----------------------------------------------------------------------------

void OsiLagSolverInterface::rowRimResize_(const int newSize){
    if (newSize > maxNumrows_) {
        double* rub   = rowupper_;
        double* rlb   = rowlower_;
        char*   sense = rowsense_;
        double* right = rhs_;
        double* range = rowrange_;
        double* dual  = rowprice_;
        double* left  = lhs_;
        maxNumrows_ = CoinMax(1000, (newSize * 5) / 4);
        rowRimAllocator_();
        const int rownum = getNumRows();
        CoinDisjointCopyN(rub  , rownum, rowupper_);
        CoinDisjointCopyN(rlb  , rownum, rowlower_);
        CoinDisjointCopyN(sense, rownum, rowsense_);
        CoinDisjointCopyN(right, rownum, rhs_);
        CoinDisjointCopyN(range, rownum, rowrange_);
        CoinDisjointCopyN(dual , rownum, rowprice_);
        CoinDisjointCopyN(left , rownum, lhs_);
        delete[] rub;
        delete[] rlb;
        delete[] sense;
        delete[] right;
        delete[] range;
        delete[] dual;
        delete[] left;
    }
}

//-----------------------------------------------------------------------------

void OsiLagSolverInterface::colRimResize_(const int newSize){
    if (newSize > maxNumcols_) {
        double* cub = colupper_;
        double* clb = collower_;
        bool* cont  = continuous_;
        double* obj = objcoeffs_;
        double* sol = colsol_;
        double* rc  = rc_;
        maxNumcols_ = CoinMax(1000, (newSize * 5) / 4);
        colRimAllocator_();
        const int colnum = getNumCols();
        CoinDisjointCopyN(cub , colnum, colupper_);
        CoinDisjointCopyN(clb , colnum, collower_);
        CoinDisjointCopyN(cont, colnum, continuous_);
        CoinDisjointCopyN(obj , colnum, objcoeffs_);
        CoinDisjointCopyN(sol , colnum, colsol_);
        CoinDisjointCopyN(rc  , colnum, rc_);
        delete[] cub;
        delete[] clb;
        delete[] cont;
        delete[] obj;
        delete[] sol;
        delete[] rc;
    }
}

//-----------------------------------------------------------------------------

void OsiLagSolverInterface::convertBoundsToSenses_(){
    for (int i = getNumRows() - 1; i >= 0; --i ) {
        convertBoundToSense(rowlower_[i], rowupper_[i],rowsense_[i], rhs_[i], rowrange_[i]);
    }
}

//-----------------------------------------------------------------------------

void OsiLagSolverInterface::convertSensesToBounds_(){
    for (int i = getNumRows() - 1; i >= 0; --i) {
        convertSenseToBound(rowsense_[i], rhs_[i], rowrange_[i],rowlower_[i], rowupper_[i]);
    }
}

/**************************************************************************************************/
/*                                       WarmStart related methods       	                      */      
/**************************************************************************************************/


CoinWarmStart* OsiLagSolverInterface::getEmptyWarmStart () const{ 
    return (dynamic_cast<CoinWarmStart *>(new CoinWarmStartDual())); 
}

CoinWarmStart* OsiLagSolverInterface::getWarmStart() const{
    return new CoinWarmStartDual(getNumRows(), rowprice_);
}

//-----------------------------------------------------------------------------

bool OsiLagSolverInterface::setWarmStart(const CoinWarmStart* warmstart){
    const CoinWarmStartDual* ws = dynamic_cast<const CoinWarmStartDual*>(warmstart);

    if (! ws)
        return false;

    const int ws_size = ws->size();
    if (ws_size != getNumRows() && ws_size != 0) {
        throw CoinError("wrong dual warmstart size", "setWarmStart",
            "OsiVolSolverInterface");
    }

    CoinDisjointCopyN(ws->dual(), ws_size, rowprice_);
    return true;
}

/**************************************************************************************************/
/*                                       HotStart related methods - changed!!!      	                      */      
/**************************************************************************************************/

void OsiLagSolverInterface::markHotStart(){
  delete[] rowpriceHotStart_;
  rowpriceHotStart_ = new double[getNumRows()];
  CoinDisjointCopyN(rowprice_, getNumRows(), rowpriceHotStart_);
}

void OsiLagSolverInterface::solveFromHotStart(){
   int itlimOrig = lagrangianSolver->getNbMaxIterations();
   getIntParam(OsiMaxNumIterationHotStart, lagrangianSolver->getNbMaxIterations());
   CoinDisjointCopyN(rowpriceHotStart_, getNumRows(), rowprice_);
   resolve();
   lagrangianSolver->setNbMaxIterations(itlimOrig);
}

void OsiLagSolverInterface::unmarkHotStart(){
  delete[] rowpriceHotStart_;
  rowpriceHotStart_ = NULL;
}

/**************************************************************************************************/
/*                           Problem information methods (original data)       	                  */      
/**************************************************************************************************/

bool OsiLagSolverInterface::isContinuous(int colNumber) const{
  assert( continuous_!=NULL );
  if ( continuous_[colNumber] ) return true;
  return false;
}

//-----------------------------------------------------------------------------

const CoinPackedMatrix *OsiLagSolverInterface::getMatrixByRow() const {
   updateRowMatrix_();
   return &rowMatrix_;
}

//-----------------------------------------------------------------------

const CoinPackedMatrix *OsiLagSolverInterface::getMatrixByCol() const {
   updateColMatrix_();
   return &colMatrix_;
}

/**************************************************************************************************/
/*                               Problem information methods (results)       	                  */      
/**************************************************************************************************/

std::vector<double*> OsiLagSolverInterface::getDualRays(int /*maxNumRays*/,bool /*fullRay*/) const{
    // *FIXME* : must write the method -LL
    throw CoinError("method is not yet written", "getDualRays",
            "OsiVolSolverInterface");
    return std::vector<double*>();
}

//------------------------------------------------------------------
std::vector<double*> OsiLagSolverInterface::getPrimalRays(int /*maxNumRays*/) const{
    // *FIXME* : must write the method -LL
    throw CoinError("method is not yet written", "getPrimalRays",
            "OsiVolSolverInterface");
    return std::vector<double*>();
}

/**************************************************************************************************/
/*                               Problem modifying methods (rim vectors)       	                  */      
/**************************************************************************************************/

void OsiLagSolverInterface::setColSetBounds(const int* indexFirst,const int* indexLast,const double* boundList){
    while (indexFirst < indexLast) {
        const int ind = *indexFirst;
        collower_[ind] = boundList[0];
        colupper_[ind] = boundList[1];
        ++indexFirst;
        boundList += 2;
    }
}

//-----------------------------------------------------------------------------
void OsiLagSolverInterface::setRowSetBounds(const int* indexFirst,const int* indexLast,const double* boundList){
    if (indexLast - indexFirst < getNumRows() / 3) {
        while (indexFirst < indexLast) {
        setRowBounds(*indexFirst, boundList[0], boundList[1]);
        ++indexFirst;
        boundList += 2;
        }
    } else {
        // it's better to convert everything at once
        while (indexFirst < indexLast) {
        const int ind = *indexFirst;
        rowlower_[ind] = boundList[0];
        rowupper_[ind] = boundList[1];
        ++indexFirst;
        boundList += 2;
        }
        convertBoundsToSenses_();
    }
}

//-----------------------------------------------------------------------------
void OsiLagSolverInterface::setRowSetTypes(const int* indexFirst,const int* indexLast,const char* senseList,const double* rhsList,const double* rangeList){
    if (indexLast - indexFirst < getNumRows() / 3) {
        while (indexFirst < indexLast) {
        setRowType(*indexFirst++, *senseList++, *rhsList++, *rangeList++);
        }
    } else {
        // it's better to convert everything at once
        while (indexFirst < indexLast) {
        const int ind = *indexFirst++;
        rowsense_[ind] = *senseList++;
        rhs_[ind] = *rhsList++;
        rowrange_[ind] = *rangeList++;
        }
        convertSensesToBounds_();
    }
}

//-----------------------------------------------------------------------------

void OsiLagSolverInterface::setContinuous(int index){
    assert(continuous_ != NULL);
    if (index < 0 || index > getNumCols()) {
        throw CoinError("Index out of bound.", "setContinuous",
            "OsiVolSolverInterface");
    }
    continuous_[index] = true;
}

//-----------------------------------------------------------------------

void OsiLagSolverInterface::setInteger(int index){
    assert(continuous_ != NULL);
    if (index < 0 || index > getNumCols()) {
        throw CoinError("Index out of bound.", "setContinuous",
            "OsiVolSolverInterface");
    }
    continuous_[index] = false;
}

//-----------------------------------------------------------------------

void OsiLagSolverInterface::setContinuous(const int* indices, int len){
    assert(continuous_ != NULL);
    const int colnum = getNumCols();
    int i;

    for (i = len - 1; i >= 0; --i) {
        if (indices[i] < 0 || indices[i] > colnum) {
        throw CoinError("Index out of bound.", "setContinuous",
                "OsiVolSolverInterface");
        }
    }
    
    for (i = len - 1; i >= 0; --i) {
        continuous_[indices[i]] = true;
    }
}

//-----------------------------------------------------------------------

void OsiLagSolverInterface::setInteger(const int* indices, int len){
    assert(continuous_ != NULL);
    const int colnum = getNumCols();
    int i;

    for (i = len - 1; i >= 0; --i) {
        if (indices[i] < 0 || indices[i] > colnum) {
        throw CoinError("Index out of bound.", "setContinuous",
                "OsiVolSolverInterface");
        }
    }
    
    for (i = len - 1; i >= 0; --i) {
        continuous_[indices[i]] = false;
    }
}

//-----------------------------------------------------------------------

void OsiLagSolverInterface::setColSolution(const double *colsol){
    CoinDisjointCopyN(colsol, getNumCols(), colsol_);
    // Compute the left hand side (row activity levels)
    colMatrix_.times(colsol_, lhs_);
}

//-----------------------------------------------------------------------

void OsiLagSolverInterface::setRowPrice(const double *rowprice){
    CoinDisjointCopyN(rowprice, getNumRows(), rowprice_);
    compute_rc_(rowprice_, rc_);
}

/**************************************************************************************************/
/*                                 Problem modifying methods (matrix)       	                  */      
/**************************************************************************************************/

void OsiLagSolverInterface::addCol(const CoinPackedVectorBase& vec,const double collb, const double colub, const double obj){
    const int colnum = getNumCols();
    colRimResize_(colnum + 1);
    collower_[colnum]   = collb;
    colupper_[colnum]   = colub;
    objcoeffs_[colnum]  = obj;
    continuous_[colnum] = true;
    colsol_[colnum]     = fabs(collb)<fabs(colub) ? collb : colub;
    rc_[colnum]         = 0.0;

    updateColMatrix_();
    colMatrix_.appendCol(vec);
    rowMatrixCurrent_ = false;
}

//-----------------------------------------------------------------------------

void OsiLagSolverInterface::addCols(const int numcols,const CoinPackedVectorBase * const * cols,const double* collb, const double* colub,   const double* obj){
    if (numcols > 0) {
        const int colnum = getNumCols();
        colRimResize_(colnum + numcols);
        CoinDisjointCopyN(collb, numcols, collower_ + colnum);
        CoinDisjointCopyN(colub, numcols, colupper_ + colnum);
        CoinDisjointCopyN(obj, numcols, objcoeffs_ + colnum);
        CoinFillN(continuous_ + colnum, numcols, true);
        int c;
        for ( c=0; c<numcols; c++ ) {
        if ( fabs(collb[c]) < fabs(colub[c]) ) {
            colsol_[colnum+c] = collb[c];
        }
        else {
            colsol_[colnum+c] = colub[c];
        }
        }
        //CoinFillN(colsol_     + colnum, numcols, 0.0);
        CoinFillN(rc_         + colnum, numcols, 0.0);

        updateColMatrix_();
        colMatrix_.appendCols(numcols, cols);
        rowMatrixCurrent_ = false;
    }
}

//-----------------------------------------------------------------------------

void OsiLagSolverInterface::deleteCols(const int num, const int * columnIndices){
    if (num > 0) {
        int * delPos = new int[num];
        CoinDisjointCopyN(columnIndices, num, delPos);
        std::sort(delPos, delPos + num);
        const int delNum = std::unique(delPos, delPos + num) - delPos;

        const int colnum = getNumCols();
        CoinDeleteEntriesFromArray(collower_, collower_ + colnum,delPos, delPos + delNum);
        CoinDeleteEntriesFromArray(colupper_, colupper_ + colnum,delPos, delPos + delNum);
        CoinDeleteEntriesFromArray(objcoeffs_, objcoeffs_ + colnum,delPos, delPos + delNum);
        CoinDeleteEntriesFromArray(continuous_, continuous_ + colnum,delPos, delPos + delNum);
        CoinDeleteEntriesFromArray(colsol_, colsol_ + colnum,delPos, delPos + delNum);
        CoinDeleteEntriesFromArray(rc_, rc_ + colnum,delPos, delPos + delNum);

        updateColMatrix_();
        colMatrix_.deleteCols(delNum, delPos);
        rowMatrixCurrent_ = false;
    }
}

//-----------------------------------------------------------------------------

void OsiLagSolverInterface::addRow(const CoinPackedVectorBase& vec,const double rowlb, const double rowub){
    const int rownum = getNumRows();
    rowRimResize_(rownum + 1);
    rowlower_[rownum] = rowlb;
    rowupper_[rownum] = rowub;
    convertBoundToSense(rowlb, rowub,rowsense_[rownum], rhs_[rownum], rowrange_[rownum]);
    rowprice_[rownum] = 0.0;
    lhs_[rownum] = 0.0;

    updateRowMatrix_();
    rowMatrix_.appendRow(vec);
    colMatrixCurrent_ = false;
}

//-----------------------------------------------------------------------------

void OsiLagSolverInterface::addRow(const CoinPackedVectorBase& vec,const char rowsen, const double rowrhs,const double rowrng){
    const int rownum = getNumRows();
    rowRimResize_(rownum + 1);
    rowsense_[rownum] = rowsen;
    rhs_[rownum] = rowrhs;
    rowrange_[rownum] = rowrng;
    convertSenseToBound(rowsen, rowrhs, rowrng,rowlower_[rownum], rowupper_[rownum]);
    rowprice_[rownum] = 0.0;
    lhs_[rownum] = 0.0;

    updateRowMatrix_();
    rowMatrix_.appendRow(vec);
    colMatrixCurrent_ = false;
}

//-----------------------------------------------------------------------------

void OsiLagSolverInterface::addRows(const int numrows, const CoinPackedVectorBase * const * rows, const double* rowlb, const double* rowub){
    if (numrows > 0) {
        const int rownum = getNumRows();
        rowRimResize_(rownum + numrows);
        CoinDisjointCopyN(rowlb, numrows, rowlower_ + rownum);
        CoinDisjointCopyN(rowub, numrows, rowupper_ + rownum);
        for (int i = rownum + numrows - 1; i >= rownum; --i) {
        convertBoundToSense(rowlower_[i], rowupper_[i],
                rowsense_[i], rhs_[i], rowrange_[i]);
        }
        CoinFillN(rowprice_ + rownum, numrows, 0.0);
        CoinFillN(lhs_      + rownum, numrows, 0.0);

        updateRowMatrix_();
        rowMatrix_.appendRows(numrows, rows);
        colMatrixCurrent_ = false;
    }
}

//-----------------------------------------------------------------------------

void OsiLagSolverInterface::addRows(const int numrows,const CoinPackedVectorBase * const * rows,const char* rowsen, const double* rowrhs,const double* rowrng){
  if (numrows > 0) {
    const int rownum = getNumRows();
    rowRimResize_(rownum + numrows);
    CoinDisjointCopyN(rowsen, numrows, rowsense_ + rownum);
    CoinDisjointCopyN(rowrhs, numrows, rhs_ + rownum);
    CoinDisjointCopyN(rowrng, numrows, rowrange_ + rownum);
    for (int i = rownum + numrows - 1; i >= rownum; --i) {
      convertSenseToBound(rowsense_[i], rhs_[i], rowrange_[i], rowlower_[i], rowupper_[i]);
    }
    CoinFillN(rowprice_ + rownum, numrows, 0.0);
    CoinFillN(lhs_      + rownum, numrows, 0.0);

    updateRowMatrix_();
    rowMatrix_.appendRows(numrows, rows);
    colMatrixCurrent_ = false;
  }
}

//-----------------------------------------------------------------------------

void  OsiLagSolverInterface::deleteRows(const int num, const int * rowIndices){
    if (num > 0) {
        int * delPos = new int[num];
        CoinDisjointCopyN(rowIndices, num, delPos);
        std::sort(delPos, delPos + num);
        const int delNum = std::unique(delPos, delPos + num) - delPos;

        const int rownum = getNumRows();
        CoinDeleteEntriesFromArray(rowlower_, rowlower_ + rownum, delPos, delPos + delNum);
        CoinDeleteEntriesFromArray(rowupper_, rowupper_ + rownum, delPos, delPos + delNum);
        CoinDeleteEntriesFromArray(rowsense_, rowsense_ + rownum, delPos, delPos + delNum);
        CoinDeleteEntriesFromArray(rowrange_, rowrange_ + rownum, delPos, delPos + delNum);
        CoinDeleteEntriesFromArray(rhs_, rhs_ + rownum, delPos, delPos + delNum);
        CoinDeleteEntriesFromArray(rowprice_, rowprice_ + rownum, delPos, delPos + delNum);
        CoinDeleteEntriesFromArray(lhs_, lhs_ + rownum, delPos, delPos + delNum);

        updateRowMatrix_();
        rowMatrix_.deleteRows(delNum, delPos);
        colMatrixCurrent_ = false;

        delete[] delPos;
    }
}

/**************************************************************************************************/
/*                                          Applying cuts        	    		          */			      
/**************************************************************************************************/

void OsiLagSolverInterface::applyRowCut(const OsiRowCut& rc){
    const int rownum = getNumRows();
    const double lb = rc.lb();
    const double ub = rc.ub();
    rowRimResize_(rownum + 1);
    rowprice_[rownum] = 0.0;
    rowlower_[rownum] = lb;
    rowupper_[rownum] = ub;
    convertBoundToSense(lb, ub,rowsense_[rownum], rhs_[rownum], rowrange_[rownum]);
    updateRowMatrix_();
    rowMatrix_.appendRow(rc.row());
    colMatrixCurrent_ = false;
}

//-----------------------------------------------------------------------

void OsiLagSolverInterface::applyColCut(const OsiColCut& cc){
    int i;

    const double* lb_elem = cc.lbs().getElements();
    const int* lb_ind = cc.lbs().getIndices();
    for (i = cc.lbs().getNumElements() - 1; i >= 0; --i) {
        collower_[lb_ind[i]] = CoinMax(collower_[lb_ind[i]], lb_elem[i]);
    }
    
    const double* ub_elem = cc.ubs().getElements();
    const int* ub_ind = cc.ubs().getIndices();
    for (i = cc.ubs().getNumElements() - 1; i >= 0; --i) {
        colupper_[ub_ind[i]] = CoinMin(colupper_[ub_ind[i]], ub_elem[i]);
    }
}

/**************************************************************************************************/
/*                                             ADDED        	    		                      */			      
/**************************************************************************************************/

void OsiLagSolverInterface::loadProblem(const CoinPackedMatrix& matrix,const double* collb, const double* colub,const double* obj,const double* rowlb, const double* rowub){}			    
void OsiLagSolverInterface::assignProblem(CoinPackedMatrix*& matrix,double*& collb, double*& colub, double*& obj,double*& rowlb, double*& rowub){}   			    
void OsiLagSolverInterface::loadProblem(const CoinPackedMatrix& matrix,const double* collb, const double* colub,const double* obj,const char* rowsen, const double* rowrhs,   const double* rowrng){}
void OsiLagSolverInterface::assignProblem(CoinPackedMatrix*& matrix,double*& collb, double*& colub, double*& obj,char*& rowsen, double*& rowrhs,double*& rowrng){}
void OsiLagSolverInterface::loadProblem(const int numcols, const int numrows,const int* start, const int* index,const double* value,const double* collb, const double* colub,   const double* obj,const double* rowlb, const double* rowub){}
void OsiLagSolverInterface::loadProblem(const int numcols, const int numrows,const int* start, const int* index,const double* value,const double* collb, const double* colub,   const double* obj, const char* rowsen, const double* rowrhs,const double* rowrng){}
int OsiLagSolverInterface::readMps(const char *filename,const char *extension ){ return 0;}
void OsiLagSolverInterface::writeMps(const char *filename,const char *extension ,double objSense) const{}

/**************************************************************************************************/
/*                                             CHANGED        	    		                      */			      
/**************************************************************************************************/

/**************************************************************************************************/
/*                                     Parameter related methods        	    		          */			      
/**************************************************************************************************/

bool OsiLagSolverInterface::setIntParam(OsiIntParam key, int value){
    switch (key) {
        case OsiMaxNumIteration:
            if (value < 0)
                return false;
            lagrangianSolver->setNbMaxIterations(value);
            break;
        case OsiMaxNumIterationHotStart:
            if (value < 0)
                return false;
            OsiSolverInterface::setIntParam(key, value);
            break;
        case OsiLastIntParam:
            return false;
        default:
            return false;
    }
    return true;
}

//-----------------------------------------------------------------------------

bool OsiLagSolverInterface::setDblParam(OsiDblParam key, double value){
    switch (key) {
        case OsiDualObjectiveLimit:
            lagrangianSolver->setUBInit(value);
            break;
        case OsiPrimalObjectiveLimit: // not applicable
            return false;
        case OsiDualTolerance: // only ~0 is applicable, so accept only 1e-50 ...
            return (value == 1e-50);
        case OsiPrimalTolerance:
            if (value < 1e-04 || value > 1e-1)
                return false;
            lagrangianSolver->setPrimalAbsPrecision(value);
            break;
        case OsiObjOffset: 
            return OsiSolverInterface::setDblParam(key, value);
        case OsiLastDblParam:
            return false;
        default:
            return false;
        }
    return true;
}

//-----------------------------------------------------------------------------

bool OsiLagSolverInterface::setStrParam(OsiStrParam key, const std::string & value){
    bool retval=false;
    switch (key) {
        case OsiSolverName:
            return false;
        case OsiProbName:
            OsiSolverInterface::setStrParam(key,value);
            return retval = true;
        case OsiLastStrParam:
            return false;
        default:
            return false;
    }
    return false;
}

//-----------------------------------------------------------------------------

bool OsiLagSolverInterface::getIntParam(OsiIntParam key, int& value) const{
    switch (key) {
        case OsiMaxNumIteration:
            value = lagrangianSolver->getNbMaxIterations();
            break;
        case OsiMaxNumIterationHotStart:
            OsiSolverInterface::getIntParam(key, value);
            break;
        case OsiLastIntParam:
            return false;
        default:
            return false;
    }
    return true;
}

//-----------------------------------------------------------------------------

bool OsiLagSolverInterface::getDblParam(OsiDblParam key, double& value) const{
    switch (key) {
        case OsiDualObjectiveLimit:
            value = lagrangianSolver->getUBInit();
            break;
        case OsiPrimalObjectiveLimit: // not applicable
            return false;
        case OsiDualTolerance: // not applicable, but must return almost 0
            value = 1e-50;
            break;
        case OsiPrimalTolerance:
            value = lagrangianSolver->getPrimalAbsPrecision();
            break;
        case OsiObjOffset:
            OsiSolverInterface::getDblParam(key, value);
            break;
        case OsiLastDblParam:
            return false;
        default:
            return false;
    }
    return true;
}

//-----------------------------------------------------------------------------

bool OsiLagSolverInterface::getStrParam(OsiStrParam key, std::string & value) const{
    switch (key) {
        case OsiProbName:
            OsiSolverInterface::getStrParam(key, value);
            return true;
        case OsiSolverName:
            value = "lag";
            return true;
        case OsiLastStrParam:
            return false;
        default:
            return false;
    }
    return false;
}

/**************************************************************************************************/
/*                Methods returning info on how the solution process terminated       	          */			      
/**************************************************************************************************/

bool OsiLagSolverInterface::isAbandoned() const{
    return false;
}

bool OsiLagSolverInterface::isProvenOptimal() const{
    return (!isDualObjectiveLimitReached() && lagrangianSolver->getStatus()==AbstractLagSolver::STATUS_OPTIMAL);
}

bool OsiLagSolverInterface::isProvenPrimalInfeasible() const{
    return (lagrangianSolver->getStatus()==AbstractLagSolver::STATUS_INFEASIBLE);
}

bool OsiLagSolverInterface::isProvenDualInfeasible() const{
    // LL: *FIXME* : at the moment the volume assumes dual feasibility...
    return false;
}

bool OsiLagSolverInterface::isPrimalObjectiveLimitReached() const{
    // The volume algorithm doesn't know anything about the primal; only the
    // dual is monotone
    return false;
}

bool OsiLagSolverInterface::isDualObjectiveLimitReached() const{
    return (lagrangianSolver->getLB()  >= lagrangianSolver->getUBInit() - DBL_EPSILON);
}

bool OsiLagSolverInterface::isIterationLimitReached() const{
    return (lagrangianSolver->getStatus()==AbstractLagSolver::STATUS_MAX_IT);
}

void OsiLagSolverInterface::initialSolve(){
    // set every entry to 0.0 in the dual solution
    CoinFillN(rowprice_, getNumRows(), 0.0);

    int i;
    
    checkData_();

    /* Only one of these can do any work. */
    updateRowMatrix_();
    updateColMatrix_();

    const int dsize = getNumRows();
    const int psize = getNumCols();

    /* Negate the objective coefficients if necessary. */
    if (objsense_ < 0) {
        std::transform(objcoeffs_, objcoeffs_+psize, objcoeffs_,std::negate<double>());
    }

    lagrangianSolver->run();

    /* extract the solution */
    extractSolution();
}

//-----------------------------------------------------------------------------

/* Resolve an LP relaxation after problem modification */
void OsiLagSolverInterface::resolve(){
    int i;
    
    checkData_();

    /* Only one of these can do any work. */
    updateRowMatrix_();
    updateColMatrix_();

    const int dsize = getNumRows();
    const int psize = getNumCols();

    /* Negate the objective coefficients if necessary. */
    if (objsense_ < 0) {
        std::transform(objcoeffs_, objcoeffs_+psize, objcoeffs_,std::negate<double>());
    }

    /* Set the dual starting point */ /* it is missing to multiply by objsense_ */
    lagrangianSolver->getLagrangianFormulation()->startMultipliers(rowprice_,dsize,objsense_);
    lagrangianSolver->getLagrangianFormulation()->updateLowerUpperBound(collower_,colupper_);
    
    /* Solves the problem */
    lagrangianSolver->run(false,true);

    /* extract the solution */
    extractSolution();
    
}

void OsiLagSolverInterface::extractSolution(){
    const int dsize = getNumRows();
    const int psize = getNumCols();

    /* the lower bound on the objective value */
    lagrangeanCost_ = objsense_ * lagrangianSolver->getLB();

    /* the primal solution. */
    lagrangianSolver->getSolution(colsol_);

    // Reset the objective coefficients if necessary
    if (objsense_ < 0) {
        std::transform(objcoeffs_, objcoeffs_ + psize, objcoeffs_,std::negate<double>());
        // also, multiply the dual solution by -1
        lagrangianSolver->getDualSolution(rowprice_);
        std::transform(rowprice_, rowprice_+dsize, rowprice_,std::negate<double>());
    } else {
        // now we just have to copy the dual
        lagrangianSolver->getDualSolution(rowprice_);
    }

    /* Compute the reduced costs */
    compute_rc_(rowprice_, rc_);

    /* Compute the left hand side (row activity levels). */
    colMatrix_.times(colsol_, lhs_);
}

/**************************************************************************************************/
/*                                  CONSTRUCTORS AND COPY                                         */			      
/**************************************************************************************************/

OsiLagSolverInterface::OsiLagSolverInterface(const Instance &instance): OsiSolverInterface(),
                                            rowMatrixCurrent_(true), rowMatrix_(), colMatrixCurrent_(true),colMatrix_(),
                                            colupper_(0),collower_(0),continuous_(0),rowupper_(0),rowlower_(0),rowsense_(0),
                                            rhs_(0),rowrange_(0),objcoeffs_(0),objsense_(1.0),colsol_(0),rowprice_(0),rc_(0),
                                            lhs_(0),lagrangeanCost_(0.0),rowpriceHotStart_(0),maxNumrows_(0),maxNumcols_(0){
    lagSolverFactory factory;
    lagrangianSolver = factory.createSolver(instance);
}

OsiLagSolverInterface::OsiLagSolverInterface(const OsiLagSolverInterface& x) : OsiSolverInterface(x),
                                            rowMatrixCurrent_(true),rowMatrix_(),colMatrixCurrent_(true),colMatrix_(),
                                            colupper_(0),collower_(0),continuous_(0),rowupper_(0),rowlower_(0),rowsense_(0),
                                            rhs_(0),rowrange_(0),objcoeffs_(0),objsense_(1.0),colsol_(0),rowprice_(0),rc_(0),
                                            lhs_(0),lagrangeanCost_(0.0),rowpriceHotStart_(0),maxNumrows_(0),maxNumcols_(0){
    operator=(x);
    lagSolverFactory factory;
    Instance instance(x.lagrangianSolver->getLagrangianFormulation()->getInstance());
    lagrangianSolver = factory.createSolver(instance);
    lagrangianSolver->getLagrangianFormulation()->copyVariables(x.lagrangianSolver->getLagrangianFormulation()->getMatrixX(),x.lagrangianSolver->getLagrangianFormulation()->getMaxSliceOverall());
}

OsiLagSolverInterface& OsiLagSolverInterface::operator=(const OsiLagSolverInterface& rhs){
    if (&rhs == this)
        return *this;

    OsiSolverInterface::operator=(rhs);
    gutsOfDestructor_();

    rowMatrixCurrent_ = rhs.rowMatrixCurrent_;
    if (rowMatrixCurrent_)
        rowMatrix_ = rhs.rowMatrix_;
    colMatrixCurrent_ = rhs.colMatrixCurrent_;
    if (colMatrixCurrent_)
        colMatrix_ = rhs.colMatrix_;

    if(rhs.maxNumrows_){
        maxNumrows_ = rhs.maxNumrows_;
        rowRimAllocator_();
        const int rownum = getNumRows();
        CoinDisjointCopyN(rhs.rowupper_, rownum, rowupper_);
        CoinDisjointCopyN(rhs.rowlower_, rownum, rowlower_);
        CoinDisjointCopyN(rhs.rowsense_, rownum, rowsense_);
        CoinDisjointCopyN(rhs.rhs_, rownum, rhs_);
        CoinDisjointCopyN(rhs.rowrange_, rownum, rowrange_);
        CoinDisjointCopyN(rhs.rowprice_, rownum, rowprice_);
        CoinDisjointCopyN(rhs.lhs_, rownum, lhs_);
    }
    if(rhs.maxNumcols_){
        maxNumcols_ = rhs.maxNumcols_;
        colRimAllocator_();
        const int colnum = getNumCols();
        CoinDisjointCopyN(rhs.colupper_, colnum, colupper_);
        CoinDisjointCopyN(rhs.collower_, colnum, collower_);
        CoinDisjointCopyN(rhs.continuous_, colnum, continuous_);
        CoinDisjointCopyN(rhs.objcoeffs_, colnum, objcoeffs_);
        CoinDisjointCopyN(rhs.colsol_, colnum, colsol_);
        CoinDisjointCopyN(rhs.rc_, colnum, rc_);
    }

    lagSolverFactory factory;
    Instance instance(rhs.lagrangianSolver->getLagrangianFormulation()->getInstance());
    lagrangianSolver = factory.createSolver(instance);
    lagrangianSolver->getLagrangianFormulation()->copyVariables(rhs.lagrangianSolver->getLagrangianFormulation()->getMatrixX(),rhs.lagrangianSolver->getLagrangianFormulation()->getMaxSliceOverall());

    return *this;
}

OsiSolverInterface * OsiLagSolverInterface::clone(bool copyData) const {
  return copyData ?
    new OsiLagSolverInterface(*this) :
    new OsiLagSolverInterface(this->lagrangianSolver->getLagrangianFormulation()->getInstance());
}

/**************************************************************************************************/
/*                                   LOAD MODEL FORMULATION                                       */			      
/**************************************************************************************************/

void OsiLagSolverInterface::loadModelFormulation(){
    lagrangianSolver->initLagFormulation();
    AbstractLagFormulation* lagformulation = lagrangianSolver->getLagrangianFormulation();
    // add variables.
    setVariables(lagformulation->getLagVariables());
    // add constraints.
    setConstraints(lagformulation->getConstraints());
    // add the first objective function.
    setObjective(lagformulation->getObjFunction(0));
}

void OsiLagSolverInterface::freeFormulation(){
    lagrangianSolver->getLagrangianFormulation()->clearConstraints();
}

/* Defines the variables needed in the MIP formulation. */
void OsiLagSolverInterface::setVariables(const std::vector<Variable> &myVars){
    for (unsigned int i = 0; i < myVars.size(); i++){ 
        addCol(0,NULL,NULL, myVars[i].getLb(), myVars[i].getUb(), 0, myVars[i].getName());
        // std::cout << "Created variable: " << var[d][arc].getName() << std::endl;
    }
    for (unsigned int i = 0; i < myVars.size(); i++){ 
        int pos = myVars[i].getId();
        switch (myVars[i].getType()){
            case Variable::TYPE_BOOLEAN:
                setInteger(pos);
                break;
            case Variable::TYPE_INTEGER:
                setInteger(pos);
                break;
            case Variable::TYPE_REAL:
                break;
            default:
                std::cout << "ERROR: Variable type has not been recognized." << std::endl;
                exit(0);
                break;
        }
    }
    std::cout << "CBC variables have been defined..." << std::endl;
}

/* Defines the constraints needed in the MIP formulation. */
void OsiLagSolverInterface::setConstraints(const std::vector<Constraint> &myConstraints){
    for (unsigned int i = 0; i < myConstraints.size(); i++){ 
        CoinPackedVector constraint;
        for (unsigned int j = 0; j < myConstraints[i].getExpression().getTerms().size(); j++){
            int index = myConstraints[i].getExpression().getTerm_i(j).getVar().getId();
            double coefficient = myConstraints[i].getExpression().getTerm_i(j).getCoeff();
            constraint.insert(index, coefficient);
        }
        addRow(constraint, myConstraints[i].getLb(), myConstraints[i].getUb(), myConstraints[i].getName());
    }
    std::cout << "CBC constraints have been defined..." << std::endl;
}

/** Defines the objective function. **/
void OsiLagSolverInterface::setObjective(const ObjectiveFunction &myObjective){
    
    // Define objective sense: 1 for minimize; -1 for maximize.
    int objSense = 1;
    setObjSense(1);

    // Fill objective coefficients.
    for (unsigned int i = 0; i < myObjective.getExpression().getTerms().size(); i++){
        int index = myObjective.getExpression().getTerm_i(i).getVar().getId();
        double coefficient = myObjective.getExpression().getTerm_i(i).getCoeff();
        setObjCoeff(index, coefficient);
    }
    std::cout << "CBC objective has been defined..." << std::endl;
}

/**************************************************************************************************/
/*                                          DESTRUCTOR                                            */			      
/**************************************************************************************************/

OsiLagSolverInterface::~OsiLagSolverInterface(){
    gutsOfDestructor_();
    delete lagrangianSolver;
}