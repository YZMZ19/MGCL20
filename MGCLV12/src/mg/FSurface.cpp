/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Tolerance.h"
#include "mg/Straight.h"
#include "mg/LBRep.h"
#include "mg/TrimmedCurve.h"
#include "mg/SurfCurve.h"
#include "mg/CParam_list.h"
#include "mg/Plane.h"
#include "mg/SBRep.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "topo/Shell.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

//Comparison operator.
bool MGFSurface::operator< (const MGFSurface& f2) const{
	const MGFace* face1=dynamic_cast<const MGFace*>(this);
	const MGFace* face2=dynamic_cast<const MGFace*>(&f2);
	if(face1 && face2)//If both are Face.
		return (*face1)<(*face2);

	if((!face1) && (!face2))//If both are surface.
		return *(get_surface_pointer())< *(f2.get_surface_pointer());

	const MGObject* obj1 = object_pointer();
	const MGObject* obj2 = f2.object_pointer();
	return obj1->identify_type()<obj2->identify_type();
}

//Get the box of the object.
const MGBox& MGFSurface::get_box() const{
	const MGObject* obj=object_pointer();
	return obj->box();
}

/**
 *  @brief eval_discrete_deviationの下請け関数
 *  @param face1 一方の面データ
 *  @param face2 もう一方の面データ
 *  @param wcrv1 face1側エッジのワールドカーブ.
 *  @param wcrv2 face2側エッジのワールドカーブ.
 *  @param pcrv1 face1側エッジの面上パラメータカーブ.
 *  @param pcrv2 face2側エッジの面上パラメータカーブ.
 *
 *  wcrv1とwcrv2はほぼ共線となっている。
 */
void deviation(
	const MGSurface& srf1,///<The target 1st surface.
	const MGSurface& srf2,///<The target 2nd surface.
	const MGCurve&    wcrv1,///<The target 1st world curve of face1.
	const MGCurve&    wcrv2,///<The target 2nd world curve of face2.
	const MGCurve&    pcrv1,///<The target 1st parameter curve of face1.
	const MGCurve&    pcrv2,///<The target 2nd parameter curve of face2.
	double commonSpan[4],//parameter range of wcrv1 and wcrv2
	double commonSpanPara[4],//parameter range of pcrv1 and pcrv2
	int npoint,		///<num of inner discrete points.
	std::vector<MGPosition>& uvuvs///<face1とface2のそれぞれのパラメータ値(u1,v1), (u2,v2)が4次元の点として返る。
){
	MGPosition uvuv(4);
	double pt1s = commonSpanPara[0], pt2s = commonSpanPara[2];
	double pt1e = commonSpanPara[1], pt2e = commonSpanPara[3];
	uvuv.store_at(0, pcrv1.eval(pt1s));
	uvuv.store_at(2, pcrv2.eval(pt2s));
	uvuvs.push_back(uvuv);

	double wt1s = commonSpan[0], wt2s=commonSpan[2];
	double wt1e = commonSpan[1], wt2e = commonSpan[3];
	if(npoint>0){
		int np1 = npoint+1;
		double dwt1 = (wt1e-wt1s)/np1;
		double dwt2 = (wt2e-wt2s)/np1;
		double dpt1 = (pt1e-pt1s)/np1;
		double dpt2 = (pt2e-pt2s)/np1;

		// pos1, pos2の、それぞれの最初と最後の要素以外のすべての要素を計算する
		for(int i = 1; i < np1; i++){
			double di = double(i);
			double wt1 = wt1s+dwt1*di;
			double wt2 = wt2s+dwt2*di;
			MGPosition P1 = wcrv1.eval(wt1);
			if(!wcrv2.perp_guess(wt2s, wt2e, P1, wt2, wt2)){
				wt2 = wcrv2.closest(P1);
			}

			double pt1 = pt1s+dpt1*di;
			double pt2 = pt2s+dpt2*di;
			pt1 = srf1.param_of_pcurve(wt1, wcrv1, pcrv1, &pt1);
			pt2 = srf2.param_of_pcurve(wt2, wcrv2, pcrv2, &pt2);
			uvuv.store_at(0, pcrv1.eval(pt1), 0, 2);
			uvuv.store_at(2, pcrv2.eval(pt2), 0, 2);
			uvuvs.push_back(uvuv);
		}
	}

	// pos1, pos2それぞれの最後の点を決定する
	uvuv.store_at(0, pcrv1.eval(commonSpanPara[1]));
	uvuv.store_at(2, pcrv2.eval(commonSpanPara[3]));
	uvuvs.push_back(std::move(uvuv));
}

///Evaluate deviations of two faces(this and face2) at npoint
///discrete points.
///(1)Search the common edges which have the distance within tolerance.
///(2)Compute the nearest points from npoint discrete points of this to face2.
///Let uvuvi=uvuvs[i], then
///uvuvi[0], [1] are this face's parameter value(u1,v1), and uvuvi[2], [3] are
///parameter value(u2,v2) of face2 which is the nearest point from the point (u1, v1).
void MGFSurface::eval_discrete_deviation(
	const MGFSurface& face2,///< The target 2nd surface.
	std::vector<MGPosition>& uvuvs,///< the nearest points pairs are set.
	int npoint,		///<indicates how many discrete points be obtained.
	double tolerance///<tolerance to get two edge to compute deviation.
)const{

	const MGFSurface& face1=*this;
	const MGSurface& srf1=*(face1.get_surface_pointer());
	const MGSurface& srf2=*(face2.get_surface_pointer());

	// 1. faceの境界線を取得する。
	std::vector<UniqueCurve> outer1 = face1.outer_boundary();
	std::vector<UniqueCurve> outer2 = face2.outer_boundary();
	std::vector<UniqueCurve> outerp1 = face1.outer_boundary_param();
	std::vector<UniqueCurve> outerp2 = face2.outer_boundary_param();
	//Note that outer1 and outer01 have the same direction and the same parameter range.

	// 2. エッジ同士の離れをチェックし、近いようであれば離れ計算の対象にする。
	// ループスタート
	for(size_t i=0, n= outer1.size(); i<n; i++){
		// face1側のエッジのworld curve
		const MGCurve& wcrv1 = *outer1[i];
		const MGCurve& pcrv1 = *outerp1[i];
		MGPosition crv1s = wcrv1.start_point();
		MGPosition crv1e = wcrv1.end_point();
		for(size_t j=0; j<outer2.size(); j++){
			// face2側のエッジのworld curve
			const MGCurve& wcrv2 = *outer2[j];
			const MGCurve& pcrv2 = *outerp2[j];

			std::vector<double> spans;
			// 一時的にトレランスを緩くしてからMGCurve::common関数を呼ぶ。
			mgTolSetWCLineZero wlzeroSet(tolerance, tolerance);
			int cmnret = wcrv1.common(wcrv2, spans);
			wlzeroSet.restore();

			// 共線でなければ次のエッジの組み合わせへジャンプ。
			if(cmnret <= 0)
				continue;

			size_t k, nspan=spans.size()/4;
			double totalLen=fabs(spans[1]-spans[0]);
			for(k=1; k<nspan; k++)
				totalLen+=fabs(spans[4*k+1]-spans[4*k]);

			for(k=0; k<nspan; k++){
				double commonSpanPara[4];
				size_t fk=4*k;
				// face1側エッジの共線部分曲線を取り出す
				double* commonSpan = &(spans[fk]);
				double sp0 = spans[fk], sp1 = spans[fk+1];
				double& spara1= commonSpanPara[0]=srf1.param_of_pcurve(sp0,wcrv1,pcrv1,&sp0);
				double& epara1= commonSpanPara[1] = srf1.param_of_pcurve(sp1,wcrv1,pcrv1,&sp1);

				// face2側エッジの共線部分曲線を取り出す
				double sp2 = spans[fk+2], sp3 = spans[fk+3];
				double& spara2= commonSpanPara[2] = srf2.param_of_pcurve(sp2,wcrv2,pcrv2,&sp1);
				double& epara2= commonSpanPara[3] = srf2.param_of_pcurve(sp3,wcrv2,pcrv2,&sp3);

				// 補助関数に丸投げ
				int npointk=int(npoint*fabs(sp1-sp0)/totalLen)+1;
				deviation(srf1, srf2, wcrv1, wcrv2, pcrv1, pcrv2,commonSpan,commonSpanPara,npointk,uvuvs);
			}
		}
	}
}

//Obtain all the boundaries(i.e., outer boundary and all the inner boundaries)
std::vector<UniqueCurve> MGFSurface::get_all_boundaries(void)const{
	std::vector<UniqueCurve> crvs=outer_boundary();
	int n=number_of_inner_boundaries();
	for(int i=0; i<n; i++){
		std::vector<UniqueCurve>inneri = inner_boundary(i);
		std::move(inneri.begin(), inneri.end(), std::back_inserter(crvs));
	}
	return crvs;
}

//Obtain parameter space error.
double MGFSurface::param_error() const{
	const MGSurface& f=*(get_surface_pointer());
	return f.param_error();
}
double MGFSurface::param_error_u() const{
	const MGSurface& f=*(get_surface_pointer());
	return f.param_error_u();
}
double MGFSurface::param_error_v() const{
	const MGSurface& f=*(get_surface_pointer());
	return f.param_error_v();
}

///Return parameter value of the middle point of the surface.
///The middle point is defined as the parameter (u,v) where
///u=(param_s_u()+param_e_u())/2, and v likewise.
MGPosition MGFSurface::param_mid()const{
	MGBox uvb=param_range();
	return uvb.mid();
}

//指定点から最も近い、垂線の足とパラメータ値を返す。
//Return the foot of the perpendicular straight line from p that is 
//nearest to point P.
//Function's return value is whether point is obtained(>0) or not(0)
int MGFSurface::perp_one(
	const MGPosition& P, // 指定点(point)
	MGPosition& uv 		//Parameter value of the surface will be returned.
)const{
	MGPosition_list list=perps(P);
	int n=list.entries();
	if(n){	//Compute the nearest point.
		MGPosition_list::iterator i=list.begin();
		uv=*i++;
		MGPosition uv2;
		double dist=(eval(uv)-P).len(), dist2;
		for(;i!=list.end();i++){
			uv2=(*i);
			dist2=(eval(uv2)-P).len();
			if(dist2<dist){uv=uv2; dist=dist2;}
		}
	}
	return n;
}

//Obtain main parameter lines of the FSurface without boundaries.
//inner_skeleton includes only inner parameter lines without boundaries.
//density indicates how many inner parameter lines are necessary
//for both u and v directions.
std::vector<UniqueCurve> MGFSurface::inner_skeleton(int density) const{
	std::vector<UniqueCurve> crv_list;
	if(density>0){
		if(density>10)
			density=10;
		MGBox prange = param_range();
		const MGInterval& uspan=prange[0];
		const MGInterval& vspan=prange[1];
		double u0=uspan[0], u1=uspan[1];
		double ulen=u1-u0;
		double v0=vspan[0], v1=vspan[1];
		double vlen=v1-v0;
		double divider=double(density+1);
		for(int i=1; i<=density; i++){
			double ith=double(i)/divider;
			double u=u0+ulen*ith;
			std::vector<UniqueCurve> pcrvU= parameter_curves(1, u);
			std::move(pcrvU.begin(), pcrvU.end(), std::back_inserter(crv_list));
			double v=v0+vlen*ith;
			std::vector<UniqueCurve> pcrvV = parameter_curves(0, v);
			std::move(pcrvV.begin(), pcrvV.end(), std::back_inserter(crv_list));
		}
	}
	return crv_list;
}

//Obtain boundary and main parameter lines of the FSurface.
//skeleton includes boundary() and inner parameter lines.
//density indicates how many inner parameter lines are necessary
//for both u and v directions.
std::vector<UniqueCurve> MGFSurface::skeleton(int density) const{
	std::vector<UniqueCurve> crv_list=get_all_boundaries();
	if(density<0)
		density=0;

	std::vector<UniqueCurve> Inner = inner_skeleton(density);
	std::move(Inner.begin(), Inner.end(), std::back_inserter(crv_list));
	return crv_list;
}

//Obtain all the parameter curves at knots of u and v knot vector.
std::vector<UniqueCurve> MGFSurface::skeleton_at_knots()const{
	const MGSurface* srf=get_surface_pointer();
	const MGBox& pbox=box_param2();
	double u0=pbox[0].low_point(), u1=pbox[0].high_point();
	double v0=pbox[1].low_point(), v1=pbox[1].high_point();

	//all the parameter lines are necessary.
	const MGKnotVector& tu=srf->knot_vector_u();
	std::vector<UniqueCurve> crvs = get_all_boundaries();

	if(tu!=mgNULL_KNOT_VECTOR){
		int num1=tu.bdim()-1;
		int ku=tu.order();
		double uold=tu[ku]-1.;
		for(int i=ku; i<=num1; i++){
			double u=tu[i];
			if(u==uold)
				continue;
			if(u0<u && u<u1){
				std::vector<UniqueCurve> pcrvsi=parameter_curves(1,u);
				std::move(pcrvsi.begin(), pcrvsi.end(), std::back_inserter(crvs));
				uold=u;
			}
		}
	}
	
	const MGKnotVector& tv=srf->knot_vector_v();
	if(tv!=mgNULL_KNOT_VECTOR){
		int nvm1=tv.bdim()-1;
		int kv=tv.order();
		double vold=tv[kv]-1.;
		for(int i=kv; i<=nvm1; i++){
			double v=tv[i];
			if(v==vold)
				continue;
			if(v0<v && v<v1){
				std::vector<UniqueCurve> pcrvsi=parameter_curves(0,v);
				std::move(pcrvsi.begin(), pcrvsi.end(), std::back_inserter(crvs));
				vold=v;
			}
		}
	}
	return crvs;
}

void trimProject(
	const std::vector<const MGCurve*>& trimmers,//Trimmer curves
	const MGVector&         dir, // dir == mgNULL_VEC means pull-projection
	const MGFSurface&       surf,
	std::vector<UniqueCurve>&     result_world,
	std::vector<UniqueCurve>&     result_param
){
	std::vector<const MGCurve*>::const_iterator i=trimmers.begin(), iend=trimmers.end();
	for(; i != iend; ++i){
		const MGCurve* crvi=*i;
		if(crvi->order()==2 && crvi->bdim()>2){
			//crvi is a polyline that has multiple segments. Subdivide the polyline into one line segment.
			const MGLBRep* lbi=dynamic_cast<const MGLBRep*>(crvi);
			if(!lbi)
				continue;//This must not take place.
			int n=lbi->bdim();
			for(int j=1; j<n; j++){
				MGLBRep lineJ;
				lbi->shrinkToKnots(j, j+1, lineJ);
				surf.project(lineJ, result_param, result_world, dir);
			}
		}else
			surf.project(*crvi, result_param, result_world, dir);
	}
}

//Exclusiv program to sort network loops of build_networks() function.
bool networkcomp(
	const UniqueLoop& netw1,
	const UniqueLoop& netw2
){
	const MGBox b1=netw1->box();
	double l1=b1.len();
	const MGBox b2=netw2->box();
	double l2=b2.len();
	return l1>l2;
}

//Build networks of surf, given parameter curves vector.
void build_networks(
	const MGFSurface& surf,		//The objective surface
	const std::vector<UniqueCurve>& pcurves,//(u,v) 2D parameter curves of surf.
	std::vector<UniqueLoop>& networks	///<Built networks
){
	double err=surf.param_error()*8.;
	size_t n=pcurves.size();
	std::vector<const MGCurve*> pcurves2(n);
	for (size_t i = 0; i < n; i++) {
		pcurves2[i] = pcurves[i].get();
	}

	while(true){
		size_t npcurves=pcurves2.size();
		if(!npcurves)
			break;;

		//Search non processed curve.
		const MGCurve* curve=0;
		while(curve==0 && npcurves){
			npcurves--;
			curve=pcurves2[npcurves];pcurves2.pop_back();
		}
		if(!curve)
			break;//If non processed curves not found.

		
		int j, nnet=(int)networks.size();
		for(j=nnet-1; j>=0; j--){
			MGLoop& netj=*networks[j];
			if(netj.merge_network(*curve,err)){
				break;
			}
		}
		if(j>=0)
			continue;

		std::unique_ptr<MGLoop> networki(new MGLoop(new MGEdge(curve->clone())));
		if(curve->is_closedWithError(err))
			networki->make_close();
		npcurves=pcurves2.size();
		bool merged=true;
		while(merged){
			merged=false;
			for(int i=(int)(npcurves-1); i>=0; i--){
				const MGCurve* pcurve2i=pcurves2[i];
				if(pcurve2i==0)
					continue;
				bool mergedi=networki->merge_network(*pcurve2i,err);
				if(mergedi){
					pcurves2[i]=0;
					merged=true;
				}
			}
		}
		networks.emplace_back(networki.release());
	}

	int nnet=(int)networks.size();
	for(int j=nnet-1; j>=0; j--){
		MGLoop& netj=*networks[j];
		bool deleted;
		do{
			deleted=netj.remove_pendent_edge(surf);
		}while(deleted);
		netj.remove_garbage_edge(err);
		if(netj.number_of_edges()==0)
			networks.erase(networks.begin()+j);
	}
	std::sort(networks.begin(),networks.end(),networkcomp);
}

//Trim this fsurface with trimmers. trimmers are 3D curves and will be projected
//onto this surface tword the direction dir. If dir is null vector, surface normal
//prjection will be performed. Trimming is so performed that the smallest region enclosed
//by trimmers that includes the surface point uv will be removed. 
void MGFSurface::trim(
	const std::vector<const MGCurve*>& trimmers,//Trimmer curves
	const MGVector&  dir,	//trimmers projection direction.
	const MGPosition& uv,	//surface parameter (u,v) that indicates the region to remove.
							//The smallest region that inclued uv will be removed.
	std::vector<UniqueFace>& faces//Result trimmed face(s) will be appended.
			//If no trimming was performed, no faces will be appended.
)const{
	// get surf-curve by project
	std::vector<UniqueCurve> wcurves, pcurves;
	trimProject(trimmers,dir,*this,wcurves,pcurves);

	std::vector<UniqueLoop> networks;	//Built networks
	build_networks(*this,pcurves,networks);

	MGFace face(*this);
	face.remove_inactive_loops();
	face.make_outer_boundary();
	face.trim(networks,uv,faces);
}
void MGFSurface::trim(
	const std::vector<UniqueCurve>& trimmers,	///<Trimmer curves
	const MGVector&  dir,	///<trimmers projection direction.
	const MGPosition& uv,	///<surface parameter (u,v) that indicates the region to remove,
							///<The smallest region that inclued uv will be removed.
	std::vector<UniqueFace>& faces///<Result trimmed face(s) will be appended,
			///<If no trimming was performed, no faces will be appended.
)const {
	std::vector<const MGCurve*> trimmers2;	//Trimmer curves
	trimmers2.resize(trimmers.size());
	extractConstPointerVec(trimmers.begin(),trimmers.end(), trimmers2.begin());
	trim(trimmers2, dir, uv, faces);
}

//Extract a sub surface with trimmers. trimmers are 3D curves and will be projected
//onto this surface toword the direction dir. If dir is null vector, surface normal
//projection will be performed. Extraction is so performed that the smallest region
//enclosed by trimmers that includes the surface point uv is extracted. 
void MGFSurface::extract(
	const std::vector<const MGCurve*>& trimmers,	//Trimmer curves
	const MGVector&  dir,	//trimmers projection direction.
	const MGPosition& uv,	//surface parameter (u,v) that indicates the region to extract.
							//The smallest region that inclued uv will be extracted.
	std::unique_ptr<MGFace>& eface//Result extracted face will be output.
)const{
	MGFace face(*this);
	face.remove_inactive_loops();
	face.make_outer_boundary();
	double errSave=MGTolerance::wc_zero();

	// 1 - get surf-curve by project
	std::vector<UniqueCurve> wcurves, pcurves;
	trimProject(trimmers,dir,face,wcurves,pcurves);

	std::vector<UniqueLoop> networks;	//Built networks
	build_networks(face,pcurves,networks);
	face.extract_sub_face(networks,uv,eface);
}
void MGFSurface::extract(
	const std::vector<UniqueCurve>& trimmers,	///<Trimmer curves
	const MGVector&  dir,	///<trimmers projection direction.
	const MGPosition& uv,	///<surface parameter (u,v) that indicates the region to extract.
							///<The smallest region that inclued uv will be extracted.
	std::unique_ptr<MGFace>& eface///<Result extracted face will be output.
)const {
	std::vector<const MGCurve*> trimmers2;	//Trimmer curves
	trimmers2.resize(trimmers.size());
	extractConstPointerVec(trimmers.begin(), trimmers.end(), trimmers2.begin());
	extract(trimmers2, dir, uv, eface);
}

///split this fsurface with splitters. splitters are 3D (x,y,z) curves that may not
///lie on the surface.
void MGFSurface::split(
	const std::vector<const MGCurve*>& splitters,	//splitter world curves.
	const MGVector&  dir,	//splitter projection direction.
							//If dir.is_null(), normal projection will be performed.
	std::vector<UniqueFace>& faces//Result splitted face(s) will be appended.
			//If no splitting was performed, no faces will be appended.
)const{
	// 1 - get surf-curve by project
	std::vector<UniqueCurve> wcurves, pcurves;
	trimProject(splitters,dir,*this,wcurves,pcurves);
	split(pcurves,faces);
}

///split this fsurface with splitters. splitters are 2D (u,v) surfaces's parameter curves.
void MGFSurface::split(
	const std::vector<UniqueCurve>& splitters,//splitter (u,v) curves.
	std::vector<UniqueFace>& faces//Result splitted face(s) will be appended.
			//If no splitting was performed, no faces will be appended.
)const{
	MGFace face(*this);
	face.remove_inactive_loops();
	face.make_outer_boundary();

	std::vector<UniqueLoop> networks;	//Built networks
	build_networks(face,splitters,networks);
	face.split(networks,faces);
}

MGisects MGFSurface::isectFS(const MGObject & obj2) const{
	auto curve = dynamic_cast<const MGCurve*>(&obj2);
	if(curve)
		return isectFS(*curve);
	auto surf2 = dynamic_cast<const MGSurface*>(&obj2);
	if(surf2)
		return isectFS(*surf2);
	auto face2 = dynamic_cast<const MGFace*>(&obj2);
	if(face2)
		return isectFS(*face2);
	auto shell2 = dynamic_cast<const MGShell*>(&obj2);
	if(shell2)
		return isectFS(*shell2);
	return MGisects();
}

///Intersection.
MGHHisects MGFSurface::isectFS(const MGShell& shell2) const{
	MGHHisects is = shell2.isectFS(*this);
	is.exchange12();
	return is;
}
MGSSisects MGFSurface::isectFS(const MGFSurface& fsurf) const{
	auto surf2 = dynamic_cast<const MGSurface*>(&fsurf);
	if(surf2)
		return isectFS(*surf2);
	auto face2 = dynamic_cast<const MGFace*>(&fsurf);
	if(face2)
		return isectFS(*face2);
	return MGSSisects();
}
