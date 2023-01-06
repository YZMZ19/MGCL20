/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _mgTLInputParam_HH_
#define _mgTLInputParam_HH_

////////////

class MGIfstream;
class MGOfstream;
#include "mg/MGCL.h"
#include "mg/Object.h"
#include "mgGL/VBO.h"

/** @file */
/** @addtogroup UseTessellation
 *  @{
 */

///A class that contains all the necessary input parameters to make tessellation.

///This is used to construct mgTLData(the tessellation of a surface), or for
///other parameter for tessellation.
class MG_DLL_DECLR mgTLInputParam{
public:

friend std::ostream& operator<< (std::ostream& out, const mgTLInputParam& para);

/// Serialization fucntion.
friend MGOfstream& operator<< (MGOfstream& buf, const mgTLInputParam& para);
friend MGIfstream& operator>> (MGIfstream& buf, mgTLInputParam& para);

mgTLInputParam(
	double crvTol=.15,			///<バウンダリのトレランス
	double surfTol=.2,			///<平面とみなすトレランス
	double max_ratio=2.,	///<最大アスペクト比
	MGCL::fan_kind fk=MGCL::MULTIPLE_TRIANGLES,
		///<fk=SINGLE_TRIANGLE:   1 triangle/FAN
		///<fk=MULTIPLE_TRIANGLES: as many triangles as possible/FAN
	int minimum_tri=8,	///<Specify minimum number of triangles.
	double max_edge_len=-1.	///<when max_edge_len<=0, this means no limits on an edge length.
);

///Construct from the object box data and the span length to draw object.
///span_length=MGOpenGLView::span_length().
mgTLInputParam(
	const MGObject& obj,
	double span_length
);

///////////// Operator overload./////////

////////// member function /////////////
double crvTol()const{return m_crvTol;};
void set_crvTol(double new_tol){ m_crvTol = new_tol;}

double surfTol()const{return m_surfTol;};
void set_surfTol(double new_tol){ m_surfTol = new_tol;}

double max_ratio()const{return m_max_ratio;};
void set_max_ratio(double new_ratio){ m_max_ratio = new_ratio;}

MGCL::fan_kind fanKind()const{return m_fk;};
void set_fanKind(MGCL::fan_kind new_fan){ m_fk = new_fan;}

int minimum_tri()const{return m_minimum_tri;};
void set_minimum_tri(int new_tri){ m_minimum_tri = new_tri;}

double max_edge_len()const{return m_max_edge_len;};
void set_max_edge_len(double new_len){ m_max_edge_len = new_len;}

private:
	double m_crvTol;	///<バウンダリのトレランス
	double m_surfTol;	///<平面とみなすトレランス
	double m_max_ratio;	///<最大アスペクト比
	MGCL::fan_kind m_fk;
		///< =SINGLE_TRIANGLE,
		///<	1 triangle/FAN(default) and STRIP for as many as posible triangles.
		///<	STRIP triangles may cover multiple rectangles.
		///< =MULTIPLE_TRIANGLES,
		///<	as many triangles as possible/FAN and STRIP for as many as posible triangles.
		///<	STRIP triangles may cover multiple rectangles.
		///< =SINGLE_TRIANGLE_NO_STRIP,
		///<	SINGLE_TRIANGLE, but STRIP triangles cover only one tessellated rectagle.
		///< =MULTIPLE_TRIANGLES_NO_STRIP,
		///<	MULTIPLE_TRIANGLES, but STRIP triangles cover only one tessellated rectagle.
	int m_minimum_tri;	///<Specify minimum number of triangles.
	double m_max_edge_len;	///<when max_edge_len<=0, this means no limits on an edge length.
};

/** @} */ // end of UseTessellation group
#endif
