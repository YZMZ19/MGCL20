#ifndef _mgTL2Triangles_HH_
#define _mgTL2Triangles_HH_

#include "Tl2/TL2Triangle.h"

/****************************************************************/
/*   Copyright (c) 2019 by System fugen G.K.                */
/*                       All rights reserved.                   */
/****************************************************************/

class MGFace;
class mgTL2Face;
class mgTL2Fans;
class mgTL2LPline;

/** @addtogroup UseTessellation
 *  @{
 */

///A vector of mgTL2Triangle's.

///mgTL2Triangle holds multiple triangles, and mgTL2Triangles holds vector of mgTriangle's.
///ポリゴンのベクトル保持クラス
class MG_DLL_DECLR mgTL2Triangles{

public:
	using MYELM = std::unique_ptr<mgTL2Triangle>;
	using MYVEC = std::vector<MYELM>;
	typedef MYVEC::const_iterator const_iterator;
	typedef MYVEC::iterator iterator;

friend std::ostream& operator<< (std::ostream& out, const mgTL2Triangles& tlTriangles);

//////////// constructor ///////////////

mgTL2Triangles(
	MGCL::TL_DATA_KIND dkind=MGCL::UV,
	///< Indicates if triangles's data is (u,v) data, (x,y,z) or (x,y,z) with the normal.
	///< =UV(0): (u,v) data. The space dimension of all the elements(MGPosition) of mgTL2Triangle is 2.
	///< =XYZ(1): (x,y,z) data. The space dimension of all the elements(MGPosition) of mgTL2Triangle is 3.
	///< =XYZNormal(2): (x,y,z) with (xn,yn,zn)=normal data. The space dimension of all the elements(MGPosition)
	///<     of mgTL2Triangle is 6 and includes (x,y,z,xn,yn,zn). Here (xn,yn,zn) is the unit normal
	///<     at the position (x,y,z) of the surface.
	const MGSurface* surf=0
);

///Move constructor.

////////Special member functions/////////
virtual ~mgTL2Triangles() = default;
mgTL2Triangles(const mgTL2Triangles&) = delete;
mgTL2Triangles(mgTL2Triangles&& tris) = default;
mgTL2Triangles& operator=(const mgTL2Triangles&) = delete;
mgTL2Triangles& operator=(mgTL2Triangles&&) = default;

//////////// operator overload ////////////

///Return the i-th mgTL2Triangle.
const mgTL2Triangle& operator[](int i)const;
mgTL2Triangle& operator[](int i);

////////// member function /////////////
iterator begin(){return m_triangles.begin();}
iterator end(){return m_triangles.end();}
const_iterator begin()const{return m_triangles.begin();}
const_iterator end()const{return m_triangles.end();}

MYELM& front(){return m_triangles.front();};
MYELM& back(){return m_triangles.back();};

const MYELM& front()const{return m_triangles.front();};
const MYELM& back()const{return m_triangles.back();};

///Return true if m_triangles' data is (u,v).
bool is_uv()const{return m_kind==MGCL::UV;};

///Make fan data from plinen and pline0's 1st or end point that is pivot.
///First or end is designated by FromEndOfPline0. if true from 1st point.
///The number of vertices of plinen is 2 or greater than 2.
///The fan is made from the pline0's start(end) point to n vertices of plinen.
///Only 1st(end) point of pline0 is accessed if need2ndPointOfpline0=flase.
///When need2ndPointOfpline0=true, 2nd point(neibor of 1st or end) of pline0
///is added to make an additional triangle of to the 1st of this fan:
///FromEndOfPline0=true:(pivot, 2nd point, 1st of plinen),
///FromEndOfPline0=false:(pivot, 1st of plinen, ...., n-2 of pline0).
void makeFan(
	const mgTL2LPline& plinen,
	const mgTL2LPline& pline0,
	bool need2ndPointOfpline0=false,
	bool FromEndOfPline0=false
);

///Return true if m_triangles' data is (u,v).
bool need_normal()const{return m_kind==MGCL::XYZNormal;};

///return the data kind of this.
MGCL::TL_DATA_KIND get_kind()const{return m_kind;};

void push_back(mgTL2Triangle* pTlTriangle){m_triangles.emplace_back(pTlTriangle);};
void push_back(mgTL2Triangles&& tris){
	std::move(tris.begin(), tris.end(),std::back_inserter(m_triangles));
};

///Push back a triangel of type mgTESTRIANG_FAN.
///The veritces of the triangle is verticesIDs which are vertex id of fans.
void push_back(
	const mgTL2Fans& fans,///<Target to push back.
	int nvTri,///<Number of vertices in verticesIDs. generally nvTri!=verticesIDs.size().
	const std::vector<int>& verticesIDs///< Vertex ids.
);

int size()const{return int(m_triangles.size());};

void set_surface(const MGSurface* surf){m_surface=surf;};
const MGSurface* surface()const{return m_surface;};


private:
	const MGSurface* m_surface;///<Surface pointer of this triangles.
	MGCL::TL_DATA_KIND m_kind;///< Indicates if m_triangles's data is (u,v) data,
				///< (x,y,z), or (x,y,z) with the normal.
	///< =UV(0): (u,v) data. The space dimension of all the elements(MGPosition) of mgTL2Triangle is 2.
	///< =XYZ(1): (x,y,z) data. The space dimension of all the elements(MGPosition) of mgTL2Triangle is 3.
	///< =XYZNormal(2): (x,y,z) with (xn,yn,zn)=normal data. The space dimension of all the elements(MGPosition)
	///<     of mgTL2Triangle is 6 and includes (x,y,z,xn,yn,zn). Here (xn,yn,zn) is the unit normal
	///<     at the position (x,y,z) of the surface.

	MYVEC m_triangles;	///<ポリゴンのベクトル

};

/** @} */ // end of UseTessellation group
#endif
