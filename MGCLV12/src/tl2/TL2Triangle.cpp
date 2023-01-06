#include "StdAfx.h"
#include "mg/Position.h"
#include "mg/Surface.h"
#include "topo/Face.h"
#include "Tl2/TL2Triangle.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif


mgTL2Triangle::mgTL2Triangle(
	int n,	//array length. All of the result MGPosition data included will be null.
	mgTESTRIANG typ	///<triangles type, mgTESTRIANG_FAN or mgTESTRIANG_STRIP.
):m_type(typ),m_xyzs(n){;}

static const char* type2[3]={"mgTESTRIANG_UNKNOWN","mgTESTRIANG_FAN","mgTESTRIANG_STRIP"};
std::ostream& operator<< (std::ostream& out, const mgTL2Triangle& triangle){
	out<<"TL2Triangle="<<(&triangle)<<" num of points="<<triangle.m_xyzs.size()
		<<",type="<<type2[static_cast<int>(triangle.getGeometryType())]
		<<":: "<<std::endl;

	int nm1=triangle.size()-1;
	mgTL2Triangle::const_iterator iter = triangle.begin(), itend=triangle.end();
	for(int i=0; iter!=itend; iter++, i++){
		out<<i<<":"<<*iter;
		if(i<nm1)
			out<<",";
	}
	out<<std::endl;
	return out;
}
