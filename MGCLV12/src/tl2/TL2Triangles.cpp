#include "StdAfx.h"
#include "mg/Position.h"
#include "Tl2/TL2LPline.h"
#include "Tl2/TL2Triangles.h"
#include "Tl2/TL2Fans.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//Constructor.
mgTL2Triangles::mgTL2Triangles(
	MGCL::TL_DATA_KIND dkind,
	const MGSurface* surf
):m_kind(dkind),m_surface(surf),m_triangles(){
}

//Return the i-th mgTL2Triangle.
const mgTL2Triangle& mgTL2Triangles::operator[](int i)const{
	return *(m_triangles[i]);
}
mgTL2Triangle& mgTL2Triangles::operator[](int i){
	return *(m_triangles[i]);
}

static const std::string kind_string[3]={"(u,v)", "(x,y,z)", "(x,y,z) with normal"};
std::ostream& operator<< (std::ostream& out, const mgTL2Triangles& tris){
	out<<"TLTriangles="<<(&tris)<<", surface="<<"tris.m_surface"<<",num of triangles="
		<<tris.m_triangles.size()<<", kind="<<kind_string[tris.get_kind()]<<std::endl;
	mgTL2Triangles::const_iterator iter = tris.begin();
	for(; iter != tris.end(); iter++)
		out<<**iter;
	out<<std::endl;
	return out;
}

///Push back a triangel of type mgTESTRIANG_FAN.
///The veritces of the triangle is verticesIDs which are vertex id of fans.
///The vertex ids are converted to (x,y,x) or (u,v) according to triangles.is_uv().
void mgTL2Triangles::push_back(
	const mgTL2Fans& fans,
	int nvTri,//Number of vertices in verticesIDs. generally nvTri!=verticesIDs.size().
	const std::vector<int>& verticesIDs
){
	if(nvTri<=2)
		return;

	bool uv_is_required=is_uv();
	bool normal_is_required=need_normal();
	mgTL2Triangle* cfan=new mgTL2Triangle(nvTri, mgTESTRIANG::mgTESTRIANG_FAN);
	for(int i=0; i<nvTri; i++){
		int id=verticesIDs[i];
		(*cfan)[i]= uv_is_required ? fans.uv(id) : fans.xyz(id,normal_is_required);
	}
	push_back(cfan);
}

///Make fan data from plinen and pline0's 1st or end point that is pivot.
///First or end is designated by FromEndOfPline0. if true from 1st point.
///The number of vertices of plinen is 2 or greater than 2.
///The fan is made from the pline0's start(end) point to n vertices of plinen.
///Only 1st(end) point of pline0 is accessed if need2ndPointOfpline0=flase.
///When need2ndPointOfpline0=true, 2nd point(neibor of 1st or end) of pline0
///is added to make an additional triangle of to the 1st of this fan:
///FromEndOfPline0=true:(pivot, 2nd point, 1st of plinen),
///FromEndOfPline0=false:(pivot, 1st of plinen, ...., n-2 of pline0).
void mgTL2Triangles::makeFan(
	const mgTL2LPline& plinen,
	const mgTL2LPline& pline0,
	bool need2ndPointOfpline0,
	bool FromEndOfPline0
){
	int nlpn=plinen.number_of_points();
	int nlp0=pline0.number_of_points();
	int np=need2ndPointOfpline0 ? nlpn+2:nlpn+1;
	mgTL2Triangle* fanP=new mgTL2Triangle(np, mgTESTRIANG::mgTESTRIANG_FAN);
	mgTL2Triangle& fan=*fanP;
	bool uvRequired = is_uv(), normal_is_required = need_normal();
	auto PonEdge = [uvRequired, normal_is_required](const mgTL2LPline& edge, int i) {
		return uvRequired ? edge.uv(i) : edge.xyz(i, normal_is_required);
	};

	int idStart=FromEndOfPline0 ? nlp0-1:0;
	fan[0]= PonEdge(pline0,idStart);
	int i=1;
	if(!FromEndOfPline0 && need2ndPointOfpline0)
		fan[i++]= PonEdge(pline0,1);
	for(int j=0; j<nlpn; j++)
		fan[i++]= PonEdge(plinen,j);
	if(FromEndOfPline0 && need2ndPointOfpline0)
		fan[i++]= PonEdge(pline0,nlp0-2);
	push_back(fanP);
}
