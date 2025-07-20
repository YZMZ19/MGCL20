/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGTransf_HH_
#define _MGTransf_HH_
/** @file */
/** @addtogroup BASE
 *  @{
 */

#include <glm/glm.hpp>
#include "mg/Default.h"
#include "mg/Vector.h"
#include "mg/Matrix.h"

//  Defines class MGTransf

//  Forward Declarations
class MGUnit_vector;
class MGPosition;
class MGIfstream;
class MGOfstream;
class MGIgesOfstream;
class MGObject;

///MGTransf represents a transformation of a space dimension.

///While MGMatrix is for the transformation about the origin(0,0,0),
///MGTransf is a general transformation.
///Transformation consists of a matrix M and a vector V. The matrix defines
///an affine transformation and the vector defines a translation.
///Let T be a MGTransf and A be an object. Then MGTransf transformation of the
///object A is defined as:  A*T (Not T*A), and A*T=A*M+A+V.
///T can be regarded (d+1) space dimension matrix whose last column is
///(0,...,0,1) (of dimension d+1) and the last row is (V,1.) (of dimension d+1).
class MG_DLL_DECLR MGTransf {

public:

///Scaling of the transf.
MG_DLL_DECLR friend MGTransf operator* (double scale, const MGTransf&);

///Genarate a Transf by multiplication of a vector and transf.
//MG_DLL_DECLR friend MGTransf operator* (const MGVector& v, const MGTransf& tr);

///String stream function.
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream&, const MGTransf&);

//////////// Constructor //////////

/// General Space dimension Constructor.
///When dim>=1, MGMatrix is initialized as a unit matrix.
explicit MGTransf(int dim=0);

///Construct Transf from a matrix and a vector.
///***** This is the fundamental constructor.*****
MGTransf(const MGMatrix&, const MGVector&);
MGTransf(MGMatrix&&, const MGVector&);

///Transf from same scaling factor for all coordinates and a translation.
MGTransf(double scale, const MGVector&);

///Construct Transf by copying old Transf, changing space dimension,
///and ordering of coordinates.
MGTransf (int sdim, const MGTransf& transf,
			int start1=0, int start2=0);

///////// 2D Constructor.////////////

///2D Transf from different scaling factors for x and y coordinates.
///No translation.
MGTransf(double scalex, double scaley);

/// Construct 2D space Transf to transform data for 'origin' to be
///origin(0,0,0) and for 'unit' to be x-coordimate.
MGTransf(const MGUnit_vector& unit,  ///< unit vector to be x-coordinate.
         const MGPosition& origin	 ///<origin to be origin.
);

//////////// 3D Constructor.////////////

///3D Transf from different scaling factors for x, y and z coordinates.
///No translation.
MGTransf(double scalex, double scaley, double scalez);

///3D Transf to transform P to be origin(0,0,0), uvecx to be x-axis
///and uvecy to be y-axis.

///If uvecx and uvecy does not cross at right angle, uvecy will be
///transformed.
MGTransf(const MGUnit_vector& uvecx, ///<1st unit vector.
         const MGUnit_vector& uvecy, ///<2nd unit vector
         const MGPosition& P		///<origin
 );	

///Transf to transform the line segment (P0, P1) to the line segment(Q0, Q1)
///P0 is transformed to Q0 and P1 is transformed to Q1.
///Space dimension of each points can be any number more than 1.
MGTransf(
	const MGPosition& P0, const MGPosition& P1,
    const MGPosition& Q0, const MGPosition& Q1);

//////////// Operator overload. /////////

///Reference (i,j)-th element.
double operator() (int i, int j) const{return ref(i,j);};

///Access to (i,j)-th element.
double& operator() (int i, int j) ;

///Functor to apply this transform to the object
///Function's return is the reference to object.
MGObject& operator()(MGObject& object);
MGObject* operator()(MGObject* object);

///Translation of Transf.
MGTransf operator+ (const MGVector&) const;

///Translation of Transf.
MGTransf& operator+= (const MGVector&);

///Translation of Transf.
MGTransf operator- (const MGVector&) const;

///Translation of Transf.
MGTransf& operator-= (const MGVector&);

///Scaling of the transf.
MGTransf operator* (double scale) const;

///Multiplication of transf and matrix.
MGTransf operator* (const MGMatrix&) const;

///Multiplication of two transfs.
MGTransf operator* (const MGTransf& ) const;

///Scaling of the transf.
MGTransf& operator*= (double scale);

///Multiplication of transf and matrix.
MGTransf& operator*= (const MGMatrix&);

///Multiplication of two transfs.
MGTransf& operator*= (const MGTransf&);

///Boolean operation.

///Equal comparison of two transf.
bool operator== (const MGTransf&) const;
bool operator!= (const MGTransf&) const;

//////////// Member Function /////////

///Return affine matrix of the transformation.
const MGMatrix& affine() const{return m_affine;};

///Convert this transf matrix to OpenGL Matrix.
void convert_to_glMatrix(
	glm::mat4& glMatI//double glMat[16]	///OpenGL Matrix will be output.
)const;

///Test if this is null.
bool is_null()const{return m_affine.is_null();};

///PD124=Transformation matrix.

///Function's return value is the directory entry id created.
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

///Reference (i,j)-th element.
double ref(int i, int j) const;

///Change the space dimension of this MGTransf to the new sdim.
void resize(int sdim);

///Obtain the scaling factor of this transf.
double scale()const{return m_affine.scale();};

///  Return space dimension
int sdim() const{return m_affine.sdim();};

/// Construct a transf to do the transformation of matrix around input point 
/// instead of origin of matrix, and replace own transf with it.
MGTransf& set_matrix(const MGMatrix& mat, const MGPosition& point);
MGTransf& set_matrix(MGMatrix&& mat, const MGPosition& point);
   
///Set this as a null transf.
void set_null();

///2D mirror reflection transf about a line whose direction
///is vec and passes through point P.
MGTransf& set_reflect_2D(	const MGVector& vec,
							const MGPosition& P= mgORIGIN_2D);

///Rotation 2D matrix around point P by angle.
MGTransf& set_rotate_2D(	double angle,
							const MGPosition& P= mgORIGIN_2D); 

///3D mirror reflection transf about a plane whose normal
///is vec and passes through point P.
MGTransf& set_reflect_3D(
	const MGVector&,
	const MGPosition& P= mgORIGIN);

///3D rotation matrix around the straight line whose direction is vec
///and that passes through point P.
MGTransf& set_rotate_3D(
	const MGVector& vec,
	double angle,                                      
	const MGPosition& P= mgORIGIN); 

///Update to scaling transf of same scaling value for each coordinate.
///Not change space dimension.
///No translation.
MGTransf& set_scale(double scale);

///Update to scaling trans of different scaling values for each coordinate.

///Not change space dimension.
///No translation.
///Scales are input through scale[sdim()].
MGTransf& set_diff_scale(double* scale);

///Set up this transf from OpenGL matrix.
MGTransf& set_glMatrix(const double glMat[16]);

///Return translation part of the transf.
const MGVector& translation() const{return m_translation;};

///Dump Functions
int dump_size() const;

///Dump Function
int dump(MGOfstream& ) const;

///Restore Function
int restore(MGIfstream& );

private:

//////////// Member Data /////////
  MGMatrix m_affine;		///< Affine transformation part.
  MGVector m_translation;	///< Translation part.

};


/** @} */ // end of BASE group
#endif
