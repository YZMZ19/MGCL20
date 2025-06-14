/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGNDDArray_HH_
#define _MGNDDArray_HH_
/** @file */
/** @addtogroup BASE
 *  @{
 */
#include "mg/MGCL.h"

// MGNDDArray.h
//

// Forward Declaration
class MGBPointSeq;
class MGKnotVector;
class MGKnot;
class MGIfstream;
class MGOfstream;

///Defines non-decreasing double data array.
///Used for data point abscissa, or knot vector, etc.
///MGNDDArray has size and effective data length.
class MG_DLL_DECLR MGNDDArray{

public:

///Friend Function
MG_DLL_DECLR friend MGNDDArray operator* (double scale, const MGNDDArray& nd);

///String stream Function
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream&, const MGNDDArray&);


////////Special member functions/////////
MGNDDArray():m_capacity(0),m_length(0),m_current(0),m_element(nullptr){;};///void constructor.
virtual ~MGNDDArray(){if(m_element) delete[] m_element;};
MGNDDArray(const MGNDDArray& rhs);//Copy constructor.
MGNDDArray& operator=(const MGNDDArray& rhs);//Copy assignment.
MGNDDArray(MGNDDArray&& rhs);//Move constructor.
MGNDDArray& operator=(MGNDDArray&& rhs);//Move assignment.


/// Constructor MGNDDArray of size n and lenght=n.
/// When data is given, data is an array of length n, and
/// data construct the element data of this MGNDDArray.
/// data[i] must be <= data[i+1].
explicit MGNDDArray(
	int n,		///<size of this.
	const double* data=0///<data array of length n if data!=NULL.
);

/// Construct MGNDDArray so that initial data is init and
/// is incremented by increment.
MGNDDArray(int n, double init, double increment=1.0);

/// From Data Point ordinate, obtain data point seq abscissa.
explicit MGNDDArray(const MGBPointSeq&);

/// From Data Point ordinate with End condition, obtain data point
/// seq abscissa.
MGNDDArray(MGENDCOND begin, MGENDCOND end, const MGBPointSeq&);

///Construct data point seq abscissa from knot vector with End conditions.
///This data point can be input of MGLBRep::buildByInterpolationECWithKTV,
/// or MGSBRep::buildByInterpolationECWithKTV.
///That is, generate data points tau(utaui or vtaui) from knot vector t(tu or tv),
///taking endc into account.
///The order is assumed to be 4.
MGNDDArray(MGENDCOND begin, MGENDCOND end, const MGKnotVector& t);

///From data point, obtain data point of a updated number(nnew).
///If original data point has multiplicities and nnew>=length(),
///original data point parameters and the multiplicities are preserved.
///Same as change_number().
MGNDDArray(const MGNDDArray&, int nnew);

///From data point, obtain data point of updated value range.
///Update so that (*this)(0)=ts, and (*this)(lenght()-1)=te.
///Must be ts<te.
MGNDDArray(const MGNDDArray&, double ts, double te);

/// Construct by extracting sub interval of array2.
MGNDDArray(
	int start_id,			///<Start id of array2(from 0).
	int num,					///<new array length.
	const MGNDDArray& array2	///<Original NDDArray.
);

///Construct by mixing two arrays.
///Mixing is so done as that no too close points are included.
///Although data point multiplicities of array1 are preserved,
///multiplicities of array2 are not.
///DATA POINT MULTIPLICITY IS ALLOWED ONLY IN array1.
MGNDDArray(
	int id1,			///<Start id of array1(from 0).
	int num1,		///<new array length to use of array1.
	const MGNDDArray& array1,///<Original NDDArray1.
	int id2,			///<Start id of array2(from 0).
	int num2,		///<new array length to use of array2..
	const MGNDDArray& array2///<Original NDDArray2.
);

//////////// Operator Overload ////////////
	
///Access to i-th element.
double operator[] (int i) const{return m_element[i];};
double operator() (int i) const{return m_element[i];};

///Access to i-th element.
double& operator[] (int i){return m_element[i];};
double& operator() (int i){return m_element[i];};

///Addition and subtraction of real number.
///All of the elements will be added or subtracted.
MGNDDArray operator+(double) const;
MGNDDArray& operator+=(double);
MGNDDArray operator-(double) const;
MGNDDArray& operator-=(double);

///Unary minus.
///Reverse the ordering of elements by changing all of the signs.
MGNDDArray operator-() const;

///Scaling
MGNDDArray operator*(double scale) const;
MGNDDArray& operator*= (double scale);

///Copmarison operator.
virtual bool operator== (const MGNDDArray& t2) const;
bool operator!= (const MGNDDArray& t2) const{return !operator==(t2);}

//////////// Member Function ////////////

///Add data of multiplicity 1 into data points. mult_max is the maximum
///multiplicity allowed for NDDArray.
///Return value is number of data actually added.
int add_data(double value, int mult_max=1);

///Add data with multiplicity into data points. mult_max is the maximum
///multiplicity allowed.	Return value is number of data actually added.
int add_data(const MGKnot& knot, int mult_max);

///Copy data points from array by removing the mltiple knot.
void copy_removing_multi(
	int start_id,			///<Start id of array2(from 0).
	int num,					///<new array length from start_id
	const MGNDDArray& array2	///<Original NDDArray.
);

///Return a pointer to raw data of MGNDDArray.
const double* data(int i=0) const{ return &m_element[i]; };

///Return a pointer to raw data of MGNDDArray.
double* data(int i=0) { return &m_element[i]; };

///Delete one data at index.
///Return value is new total number of data in the array, generally
///is original_length()-1.
int del_data(int index);

///Test if this is null.
bool is_null()const{return m_length==0;};

/// Return the length of MGNDDArray.
int length() const{ return m_length; };

///Set the length of effective data.
void set_length(int length);

///Set this as a null NDDArray.
virtual void set_null();

/// Store data at the position i of this array.
/// tau(i) must be => tau(i-1).
/// When i>=capacity(), reshape will take place.
/// The data validity after i is not checked.
/// The length is set same as the capacity. User must set the length by set_length().
void store_with_capacityCheck(
	int i,	///<the postion to store at.
	double data///<the data to store.
);

/// Finds index where tau is located in MGNDDArray as an index of knot.
/// 1) index=-1: tau < (*this)(0).
/// 2) 0<= index <n-1: (*this)(0) <= tau< (*this)(n-1) and
///                (*this)(index) <= tau < (*this)(index+1).
/// 3) index=n-1: (*this)(n-1) <= tau.
///Here n=lenght().
virtual int locate(double tau) const;

/// Locate where data of multiplicity of multi is after start.
/// index is the starting point index of this found first after start.
/// index>=start if index>=0.
/// Function's return value locate_multi is actual multiplicity at the
/// index, i.e. locate_multi>=multi.
/// If position of the multiplicity is not found to the end,
/// index=(lenght()-1) (index of the last element) and locate_multi=0
/// will be returned.
/// multi must be >=1.
virtual int locate_multi(int start, int multi, int& index) const;

///Update so that (*this)(0)=ts, and (*this)(lenght()-1)=te.
///Must be ts<te.
virtual void change_range(double ts, double te);

///Update array length.
///Updated array is so generated that the original proportions of
///neighbors hold as much as possible. 
///If original data point has multiplicities and nnew>=length(),
///original data point parameters and the multiplicities are preserved.
MGNDDArray& change_number(int nnew);

///Reference to i-th element.
double ref(int i) const{ return m_element[i];};

///Remove too near data points.
void remove_too_near(
	bool allow_multi=false,	///<indicates if multiple data point is allowed or not,
							///<when allow_multi=false, multiple data points will be removed,
							///<when allow_multi=true, will not be removed.
	double ratio=6.	///<maximum ratio allowed for neighboring span,
			///<let ti=(*this)[i], then
			///<if (t(i+1)-ti)/(ti-t(i-1))>ratio or (t(i+1)-ti)/(ti-t(i-1))<1/ratio,
			///<a data point will be removed(along with the ordinates).
);

///Remove too near data points. Removal will be done with the ordinates.
void remove_too_near(
	MGBPointSeq& ordinates,	///<ordinate,
							///<ordinates.length() must be equal to this->length().
	bool allow_multi=false,	///<indicates if multiple data point is allowed or not,
							///<when allow_multi=false, multiple data points will be removed,
							///<when allow_multi=true, will not be removed.
	double ratio=6.	///<maximum ratio allowed for neighboring span,
			///<let ti=(*this)[i], then
			///<if (t(i+1)-ti)/(ti-t(i-1))>ratio or (t(i+1)-ti)/(ti-t(i-1))<1/ratio,
			///<a data point will be removed(along with the ordinates).
);

/// Reserve capacity.
/// When newCapacity>=1 is input, adequate new capacity() is made to guarantee capacity() >= newCapacity.
/// When newCapacity<=0 is input, the most adequate new capacity() is made .
/// The length() and the data stored are unchanged.
void reserve(int newCapacity = 0);

/// Reserve length n area after the index id.
/// On return, length() will be increased by n.
/// Area from tau[id] to [id+n-1] will contain garbage.
/// The old data fromta[id+1] to [id+lenOld-1] will be transfered to [id+n] ...
void reserveInsertArea(int n, int id);

/// Change the size.
/// start is to indicate from which location of new area to start storing.
/// Although size can be less than original length, some of end data will be
/// lost in this case. When 'start'>0, first 'start' data will be garbage.
/// Original data will be held as long as the storage permits.
/// reshape does update the effective length to size.
void reshape(int size, int start=0);

///Resize the array. Result will contain garbages.
///length() and size() will have nsize.
void resize(int nsize);

/// Return the size of MGNDDArray.
int capacity() const { return m_capacity; };

///Obtain data point from KnotVector t and replace own with it.
MGNDDArray& buildByKnotVector(const MGKnotVector& t);

///Dump Functions.
///Calculate dump size
virtual int dump_size() const;

///Dump Function
virtual int dump(MGOfstream&)const;

///Restore Function
virtual int restore(MGIfstream& );
	
//////////// Member data /////////
protected:
	mutable int m_current;///< current interval of the array is held.

private:
	int m_capacity;	///< Size of m_element.
	int m_length;	///< Length of effective data in m_element.
	double* m_element;

};

/** @} */ // end of BASE group
#endif
