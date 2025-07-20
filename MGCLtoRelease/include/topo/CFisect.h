/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGCFisect_HH_
#define _MGCFisect_HH_
/** @addtogroup IsectContainer
 *  @{
 */

#include "mg/MGCL.h"
#include "mg/CSisect.h"

class MGPosition;
class MGFSurface;

//
//Define MGCFisect Class.

///MGCFisect is to represent an intersection of a fsurface and a curve.

///(MGCSisect csi, MGFSurface* f) where csi consists of world point, curve parameter,
///and face(surface) parameter, and f is a fsurface pointer.
class MG_DLL_DECLR MGCFisect:public MGCSisect{

public:
/////////Constructor/////////

///void constructor.
MGCFisect():m_face(nullptr){;};

///Construct from all the necessary data.
MGCFisect(
	const MGCSisect& csi,	///<isect data (point, curve parameter value,
							///<            surface parameter value)
	const MGFSurface& face		///<face.
):MGCSisect(csi), m_face(&face){;};

///Construct from all the necessary data.
MGCFisect(
	const MGPosition& point,	///<World coordinate point data of the isect.
	const double& t,			///<curve parameter value of the isect.
	const MGPosition& uv,		///<Face(Surface) parameter value of the isect.
	const MGFSurface& face		///<face.
);

/////////Operator oveload/////////

bool operator< (const MGCFisect& fp)const;
bool operator== (const MGCFisect& fp)const;

///Ordering functions.
bool operator< (const MGisect& is)const;
bool operator< (const MGCCisect& is)const{return false;};
bool operator< (const MGCSisect& is)const{return false;};
bool operator< (const MGSSisect& is)const{return true;};
//bool operator< (const MGFFisect& is)const{return true;};
bool operator== (const MGisect& is)const;

/////////Member function/////////

///Return isect data.
const MGCSisect& csi()const{return *this;};

///Exchange 1st and 2nd order of the parameter line representation.
void exchange12(){;};

///return the face.
const MGFSurface& face()const{return *m_face;};

///Return the manifold dimension of the intersection, i.e.
///0: when the intersection is a point,
///1: when                  is a curve,
///2: when                  is a surface.
int manifold_dimension()const{return 0;};

/// Output virtual function.
std::ostream& toString(std::ostream& ostrm)const;

///Return the parameter value of the surface.
/// 交点の Surface のパラメータ値を返却する。
const MGPosition& param_face() const{return MGCSisect::param_surface();};

private:
	const MGFSurface* m_face;	///<face pointer of the intersection point.

};

/** @} */ // end of IsectContainer group
#endif
