/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGFFisect_HH_
#define _MGFFisect_HH_
/** @addtogroup IsectContainer
 *  @{
 */

#include "mg/MGCL.h"
#include "mg/isect.h"
#include "mg/Curve.h"
#include "mg/FPline.h"
#include "topo/Face.h"

// MGFFisect.h
// Header for MGFFisect

//Forward Declaration
class MGObject;
class MGCurve;
class MGFace;

///MGFFisect represents one intersection line of a MGFace and MGFace or MGSurface.

class MG_DLL_DECLR MGFFisect:public MGisect{
private:	
	std::unique_ptr<MGCurve> m_curve;	///<world coordinate curve expression. A newed object.
	MGFPline m_face1line;///<1st Face's parameter line expression.
	MGFPline m_face2line;///<2nd Face's parameter line expression,
						///<may be null when the objective is not a face but a surface.

public:

////////Special member functions/////////
MGFFisect(){;};	///void constructor.
~MGFFisect()=default;			///Destructor.
MGFFisect(const MGFFisect&)=delete;///Copy constructor.
MGFFisect& operator= (const MGFFisect&)=delete;///Copy assignment.
MGFFisect(MGFFisect&&)=default;		///Move constructor.
MGFFisect& operator= (MGFFisect&&)=default;///Move assignment.

///Construct providing all the raw data.
///The ownership of iline, face1, and face2 are all transfered to MGFFisect.
///All of these objects must be newed ones.
MGFFisect(
	MGCurve* iline,	///<Pointer of a newed object.
	const MGFace* face1,///<1st face.
	MGCurve* param1,///<Pointer of a newed object.
	const MGFace* face2,///<2nd face.
	MGCurve* param2///<Pointer of a newed object.
):m_curve(iline),m_face1line(face1,param1),m_face2line(face2,param2){;};

///Construct providing all the raw data.
///Copy version. Copy of the three curves will take place.
MGFFisect(
	const MGCurve& iline,///<Reference to intersection curve.
	const MGFace* face1,///<Face 1.
	const MGCurve& param1,///<1st Face's parameter line expression.
	const MGFace* face2,///<Face 2.
	const MGCurve& param2///<2nd Face's parameter line expression.
);

///Construct providing all the raw data.
MGFFisect(
	MGCurve* iline,	///<Pointer of a newed object.
	MGFPline&& face1uv,///<Face1's FPline.
	MGFPline&& face2uv///<Face2's FPline.
);

///Comparison operator.
bool operator< (const MGFFisect& ssi2)const;
bool operator> (const MGFFisect& ssi2)const{return ssi2<(*this);};
bool operator<= (const MGFFisect& ssi2)const{return !(ssi2<(*this));};
bool operator>= (const MGFFisect& ssi2)const{return !((*this)<ssi2);};
bool operator== (const MGFFisect& ssi2)const;
bool operator!= (const MGFFisect& ssi2)const{return !operator==(ssi2);};

///Ordering functions.
bool operator< (const MGisect& is)const;
bool operator< (const MGCCisect& is)const{return false;};
bool operator< (const MGCSisect& is)const{return false;};
bool operator< (const MGCFisect& is)const{return false;};
bool operator< (const MGSSisect& is)const{return false;};
bool operator== (const MGisect& is)const;

//////////// Memeber Function ////////////

///Return the object of the intersection(world coordinates representation).
const MGObject& isect()const{return *m_curve;};

///Return the 1st object's parameter value curve of the intersection.
const MGCurve* isect1_param1()const{return &(m_face1line.uvline());};

///Return the 2nd object's parameter value curve of the intersection.
const MGCurve* isect1_param2()const{return &(m_face2line.uvline());};

///Return the manifold dimension of the intersection, i.e.
///0: when the intersection is a point,
///1: when                  is a curve,
///2: when                  is a surface.
int manifold_dimension()const{return 1;};

/// Output function.
std::ostream& toString(std::ostream& ostrm)const;

///Exchange 1st and 2nd order of the parameter line representation.
void exchange12() override;

protected:

	///Get the object1 pointer.
	const MGObject* object1(const MGObject* obj)const{
		return m_face1line.face()->object_pointer();
	};

	///Get the object2 pointer.
	const MGObject* object2(const MGObject* obj)const{
		return m_face2line.face()->object_pointer();
	};
};

/** @} */ // end of IsectContainer group
#endif
