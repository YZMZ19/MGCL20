/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "StdAfx.h"
#include "mg/tolerance.h"
#include "mg/Unit_vector.h"
#include "mg/MGStl.h"
#include "Tl2/TL2Triangle.h"
#include "Tl2/TL2Triangles.h"

//Construct a surface or a face whose data are (u,v).
MGStl::MGStl(
	const mgTL2Triangles& tris///<mgTL2Triangles whose data depend on tris.get_kind();
){
	//IdentifyPosition使用用
	//三角形の頂点の座標、頂点の番号を格納するmap
	triangleMap VertexMap;
	AddTL2Data(tris,VertexMap);
}

MGStl::MGStl(
	double error,	//Error to regard two points are the same.
	const mgTL2Triangles& tris///<mgTL2Triangles whose data depend on tris.get_kind();
){
	triangleMap VertexMap;
	AddTL2Data(error,tris,VertexMap);
}
	
///conversion constructor from tessellation data.
MGStl::MGStl(const std::vector<mgTL2Triangles>& tlDataVector){
	triangleMap VertexMap;
	std::vector<mgTL2Triangles>::const_iterator
		i=tlDataVector.begin(), iend=tlDataVector.end();
	for(;i!=iend; i++)
		AddTL2Data(*i,VertexMap);
}

///mgTL2TrianglesをMGStlに追加する
void MGStl::AddTL2Data(
	double error,	///<Error to regard two points are the same.
	const mgTL2Triangles& tris,///<mgTL2Triangles whose data depend on tris.get_kind();
	triangleMap& VertexMap
){
	mgTolSetWCZero wczeroSet(error);//Set&save the error.

	const MGSurface& surf=*(tris.surface());
	bool trisData_is_uv=tris.get_kind()==MGCL::UV;
	int ntri=tris.size();
	for(int itri=0; itri<ntri; itri++){
		const mgTL2Triangle& tri=tris[itri];
		mgTESTRIANG geoType = tri.getGeometryType();//ポリゴンの形状(ファンもしくはストリップ)
		int nVert=tri.size();
		int nVm2 = nVert-2;

		MGPosition xyz0=trisData_is_uv ? MGPosition(surf.eval(tri[0])) : tri[0];
		int id0=IdentifyPosition(xyz0, VertexMap);

		MGPosition xyz1=trisData_is_uv ? MGPosition(surf.eval(tri[1])) : tri[1];
		int id1=IdentifyPosition(xyz1,VertexMap);

		// ポリゴンごとにループ(indicesを取得する)
		for(int i=0; i<nVm2; i++){
			MGPosition xyz2=trisData_is_uv ? MGPosition(surf.eval(tri[i+2])) : tri[i+2];
			int id2=IdentifyPosition(xyz2,VertexMap);

			MGUnit_vector N;
			m_indices.push_back(id0);
			if(geoType== mgTESTRIANG::mgTESTRIANG_FAN || !(i%2)){
				m_indices.push_back(id1);
				m_indices.push_back(id2);
				N=UnitNormal(xyz0,xyz1,xyz2);
			}else{
				m_indices.push_back(id2);
				m_indices.push_back(id1);
				N=UnitNormal(xyz0,xyz2,xyz1);
			}
			// 三角形の法線ベクトルを求める
			m_vecNormlTriang.push_back(N);
			if(geoType== mgTESTRIANG::mgTESTRIANG_STRIP){
				id0=id1; xyz0=xyz1;
			}
			id1=id2; xyz1=xyz2;
		}
	}
}

//mgTL2TrianglesをMGStlに追加する
void MGStl::AddTL2Data(
	const mgTL2Triangles& tris,///<mgTL2Triangles whose data depend on tris.get_kind();
	triangleMap& VertexMap
){
	const MGSurface& surf=*(tris.surface());
	double error=surf.get_box().len()*MGTolerance::rc_zero();
	AddTL2Data(error,tris,VertexMap);
}
