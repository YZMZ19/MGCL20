#ifndef _mgTL2Triangle_HH_
#define _mgTL2Triangle_HH_

#include <vector>
#include "mg/MGCL.h"
#include "mg/Position.h"
#include "Tl2/TL2parameter.h"

/****************************************************************/
/*   Copyright (c) 2019 by System fugen G.K.                */
/*                       All rights reserved.                   */
/****************************************************************/
class MGSurface;
class MGFace;
class mgTL2Face;

/** @file */
/** @addtogroup UseTessellation
 *  @{
 */


///mgTL2Triangle holds (multiple) triangles data, which are a fan or a strip.

///The space dimension of all the elements(MGPosition) of mgTL2Triangle is 2, 3 or 6.
///When 2, it is (u,v) that is a parameter value of the target surface.
///When 3, it is (x,y,z) and when 6, it is (x,y,z,xn,yn,zn). Here (xn,yn,zn) is the unit normal
///at the position (x,y,z) of the surface.
///See OpenGL for fan or strip.
class mgTL2Triangle{

public:
	typedef std::vector<MGPosition>::iterator iterator;
	typedef std::vector<MGPosition>::const_iterator const_iterator;

friend std::ostream& operator<< (std::ostream& out, const mgTL2Triangle& triangle);

//////////// constructor ///////////////
mgTL2Triangle(
	mgTESTRIANG type = mgTESTRIANG::mgTESTRIANG_FAN	///<triangles type, mgTESTRIANG_FAN or mgTESTRIANG_STRIP.
):m_type(type){;};

mgTL2Triangle(
	int n,	///<array length. All of the result MGPosition data included will be null.
	mgTESTRIANG type = mgTESTRIANG::mgTESTRIANG_FAN	///<triangles type, mgTESTRIANG_FAN or mgTESTRIANG_STRIP.
);

//////////// Operator overload ///////////////

const MGPosition& operator[](int i)const{return m_xyzs[i];};
MGPosition& operator[](int i){return m_xyzs[i];};

//////////// Member Function ///////////////

iterator begin(){return m_xyzs.begin();};
iterator end(){return m_xyzs.end();};
const_iterator begin()const{return m_xyzs.begin();};
const_iterator end()const{return m_xyzs.end();};

///タイプを返却する mgTESTRIANG_FAN mgTESTRIANG_STRIP
mgTESTRIANG getGeometryType()const{return m_type;};

void push_back(const MGPosition& xyz){m_xyzs.push_back(xyz);};
void pop_back(){m_xyzs.pop_back();};

///タイプを設定する
void setGeometryType(mgTESTRIANG type){m_type = type;};

///Obtain the number of points included.
int size()const{return int(m_xyzs.size());};

private:
	mgTESTRIANG m_type;		///<mgTESTRIANG_FAN, or mgTESTRIANG_STRIP
	std::vector<MGPosition> m_xyzs;	///<A vector of coordinates data that makes a fan or a strip as:
		///<The space dimension of all the elements(MGPosition) of mgTL2Triangle is 2, 3 or 6.
		///<When 2, it is (u,v) that is a parameter value of the target surface.
		///<When 3, it is (x,y,z) and when 6, it is (x,y,z,xn,yn,zn). Here (xn,yn,zn) is the unit normal
		///<at the position (x,y,z) of the surface.
};

/** @} */ // end of UseTessellation group
#endif
