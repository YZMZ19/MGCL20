/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Group.h"
#include "mg/CParam_list.h"
#include "mg/Position_list.h"
#include "mg/CSisects.h"
#include "mg/SSisects.h"
#include "mg/Straight.h"
#include "mg/SurfCurve.h"
#include "mg/Curve.h"
#include "mg/LBRep.h"
#include "mg/Plane.h"
#include "mg/FSurface.h"
#include "mg/SBRep.h"
#include "mg/isects.h"
#include "mg/Tolerance.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "topo/Shell.h"
#include "topo/HHisects.h"
#include "topo/FOuterCurve.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

//
//Implements MGFace Class.
//MGFace is an instance of MGCell.

///////Member Function///////

//Compute intersection points of an inner parameter line of this face and face2.
//The intersection point is used to compute surface to surface intersection lines.
//Function's return value is at most one intersection point in uvuv_list.
//One member of uvuv_list is (u1,v1,u2,v2), where (u1,v1) is a parameter of
//this surface and (u2,v2) is a parameter of surf.
MGPosition_list MGFace::intersectInner(
	const MGFace& face2		//The second surface.
) const{
	MGPosition_list uvuv_list;
	const MGBox& pbox1=box_param(); const MGBox& pbox2=face2.box_param();

	double u10=pbox1[0].low_point(), u11=pbox1[0].high_point();
	double v10=pbox1[1].low_point(), v11=pbox1[1].high_point();
	const MGKnotVector& t1u=knot_vector_u();
	const MGKnotVector& t1v=knot_vector_v();
	int nspan1u=t1u.locate(u11,1)+1-t1u.locate(u10);
	int nspan1v=t1v.locate(v11,1)+1-t1v.locate(v10);

	double u20=pbox2[0].low_point(), u21=pbox2[0].high_point();
	double v20=pbox2[1].low_point(), v21=pbox2[1].high_point();
	const MGKnotVector& t2u=face2.knot_vector_u();
	const MGKnotVector& t2v=face2.knot_vector_v();
	int nspan2u=t2u.locate(u21,1)+1-t2u.locate(u20);
	int nspan2v=t2v.locate(v21,1)+1-t2v.locate(v20);

	int maximum;
	if(nspan1u<nspan1v){
		if(nspan1v<nspan2u){
			if(nspan2u<nspan2v) maximum=3; else maximum=2;
		}else{
			if(nspan1v<nspan2v) maximum=3; else maximum=1;
		}
	}else{
		if(nspan1u<nspan2u){
			if(nspan2u<nspan2v) maximum=3; else maximum=2;
		}else{
			if(nspan1u<nspan2v) maximum=3; else maximum=0;
		}
	}
	const MGFace *f1, *f2;
	int isU=1; if(maximum%2) isU=0;
	double t0,t1;
	int nspan;
	if(maximum<=1){
		f1=this; f2=&face2;
		if(isU){ t0=u10; t1=u11; nspan=nspan1u;}
		else{ t0=v10; t1=v11; nspan=nspan1v;}
	}else{
		f2=this; f1=&face2;
		if(isU){ t0=u20; t1=u21; nspan=nspan2u;}
		else{ t0=v20; t1=v21; nspan=nspan2v;}
	}
	double tau=(t0+t1)*.5;
	double delta=(t1-t0)/double(nspan);
	int sign=-1;
	MGCSisects csiList;
	for(int i=0; i<nspan; i++){
		tau+=double(i*sign)*delta;
		std::vector<UniqueCurve> crvs;
		if(tau>t0 && tau<t1)
			crvs=f1->parameter_curves(isU,tau);
		for(size_t j=0; j<crvs.size(); j++){
			csiList=crvs[j]->isect(*f2);
			if(csiList.size()) break;
		}
		if(csiList.size()) break;
		sign*=-1;
	}
	if(!csiList.size())
		return uvuv_list;

	double u1,v1;
	auto& csi=isectCast<MGCSisect>(csiList.begin());
	u1=tau; v1=csi.param_curve();
	if(!isU){ u1=v1; v1=tau;}
	const MGPosition& uv2=csi.param_surface();
	MGPosition uvuv(4);
	if(f1==this){
		uvuv(0)=u1; uvuv(1)=v1; uvuv.store_at(2,uv2,0,2);
	}else{
		uvuv(2)=u1; uvuv(3)=v1; uvuv.store_at(0,uv2,0,2);
	}
	uvuv_list.append(uvuv);
	return uvuv_list;
}

//Compute intersection points of an inner parameter line of this face and sf2.
//The intersection point is used to compute surface to surface intersection lines.
//Function's return value is at most one intersection point in uvuv_list.
//One member of uvuv_list is (u1,v1,u2,v2), where (u1,v1) is a parameter of
//this surface and (u2,v2) is a parameter of surf.
MGPosition_list MGFace::intersectInner(
	const MGSurface& sf2		//The second surface.
) const{
	MGPosition_list uvuv_list;
	const MGBox& pbox1=box_param();

	double u10=pbox1[0].low_point(), u11=pbox1[0].high_point();
	double v10=pbox1[1].low_point(), v11=pbox1[1].high_point();
	const MGKnotVector& t1u=knot_vector_u();
	const MGKnotVector& t1v=knot_vector_v();
	int nspan1u=t1u.locate(u11,1)+1-t1u.locate(u10);
	int nspan1v=t1v.locate(v11,1)+1-t1v.locate(v10);

	double u20=sf2.param_s_u(), u21=sf2.param_e_u();
	double v20=sf2.param_s_v(), v21=sf2.param_e_v();
	const MGKnotVector& t2u=sf2.knot_vector_u();
	const MGKnotVector& t2v=sf2.knot_vector_v();
	int nspan2u=t2u.locate(u21,1)+1-t2u.locate(u20);
	int nspan2v=t2v.locate(v21,1)+1-t2v.locate(v20);

	int maximum;
	if(nspan1u<nspan1v){
		if(nspan1v<nspan2u){
			if(nspan2u<nspan2v) maximum=3; else maximum=2;
		}else{
			if(nspan1v<nspan2v) maximum=3; else maximum=1;
		}
	}else{
		if(nspan1u<nspan2u){
			if(nspan2u<nspan2v) maximum=3; else maximum=2;
		}else{
			if(nspan1u<nspan2v) maximum=3; else maximum=0;
		}
	}
	int isU=1; if(maximum%2) isU=0;
	double t0,t1,tau;
	int nspan;
	int sign=-1;
	MGCSisects csiList;
	if(maximum<=1){
		if(isU){ t0=u10; t1=u11; nspan=nspan1u;}
		else{ t0=v10; t1=v11; nspan=nspan1v;}
		tau=(t0+t1)*.5;
		double delta=(t1-t0)/double(nspan);
		for(int i=0; i<nspan; i++){
			tau+=double(i*sign)*delta;
			std::vector<UniqueCurve> crvs;
			if(tau>t0 && tau<t1)
				crvs=parameter_curves(isU,tau);
			for(size_t j=0; j<crvs.size(); j++){
				csiList=sf2.isect(*(crvs[j]));
				if(csiList.size()) break;
			}
			if(csiList.size()) break;
			sign*=-1;
		}
	}else{
		if(isU){ t0=u20; t1=u21; nspan=nspan2u;}
		else{ t0=v20; t1=v21; nspan=nspan2v;}
		tau=(t0+t1)*.5;
		double delta=(t1-t0)/double(nspan);
		for(int i=0; i<nspan; i++){
			tau+=double(i*sign)*delta;
			if(tau>t0 && tau<t1){
				MGCurve* crv=sf2.parameter_curve(isU,tau);
				csiList=isect(*crv); delete crv;
				if(csiList.size()) break;
			}
			if(csiList.size()) break;
			sign*=-1;
		}
	}
	if(!csiList.size()) return uvuv_list;

	double u1,v1;
	auto& csi=isectCast<MGCSisect>(csiList.begin());
	u1=tau; v1=csi.param_curve();
	if(!isU){ u1=v1; v1=tau;}
	const MGPosition& uv2=csi.param_surface();
	MGPosition uvuv(4);
	if(maximum<=1){
		uvuv(0)=u1; uvuv(1)=v1; uvuv.store_at(2,uv2,0,2);
	}else{
		uvuv(2)=u1; uvuv(3)=v1; uvuv.store_at(0,uv2,0,2);
	}
	uvuv_list.append(uvuv);
	return uvuv_list;
}

///Compute the intersections of two objects.
///Intersections are obtained from two objects, which are known using
///the MGisects::object1() and object2().
///****NOTE****
///When two objects' manifold dimension are the same, object1 is this object
///at the invocation of MGObject::intersection(), and object2 is the argument
///object.
///However, their manifold dimension are not the same, object1 is always
///the lower dimension's object and object2 is the higer dimension's object.
MGisects MGFace::isect(const MGObject& obj2)const{
	auto curve = dynamic_cast<const MGCurve*>(&obj2);
	if(curve)
		return curve->isect(*this);
	auto surf2 = dynamic_cast<const MGSurface*>(&obj2);
	if(surf2)
		return intersect(*surf2);
	auto face2 = dynamic_cast<const MGFace*>(&obj2);
	if(face2)
		return intersect(*face2);
	auto shell2 = dynamic_cast<const MGShell*>(&obj2);
	if(shell2){
		const MGFSurface* fs = static_cast<const MGFSurface*>(this);
		MGHHisects is=shell2->isectFS(*fs);
		is.exchange12();
		return std::move(is);
	}
	return MGisects();
}

//Curve to face intersection.
MGCSisects MGFace::isect(const MGCurve& curv)const{
	return curv.isect(*this);
}

MGSSisects MGFace::isect(const MGFSurface& surf2)const{
	return intersect(surf2);
}
MGSSisects MGFace::isect(const MGFace& face2)const{
	return intersect(face2);
}
MGSSisects MGFace::isect(const MGSurface& surf2)const{
	return intersect(surf2);
}

//Face to Face intersection.
MGSSisects MGFace::intersect(const MGFSurface& face2) const{
	if(!face2.has_commonFS(*this))
		return MGSSisects(this,&face2);

	MGPosition_list uvuv_list;

	//Compute intersection points of this face's boundary and face2.
	//Will be stored in uvuv_list. The one member is uvuv(7):
	//0,1:(u,v) of this, 2,3:(u,v) of face2, 4-6:direction of isect.
	intersect12Boundary(face2,uvuv_list);
	return isectEP(uvuv_list,face2);
}

//Obtain parameter u(kcod=1), or v(kcod=0) of the intersection point of
//v=x(const) line(kcod=1), or u=x(const) line (kcod=0) with
//all the face boundaries.
std::vector<double> MGFace::isect1D_with_boundaries(
	double x,							//coordinate value of kcod.
	int kcod)						//Coordinate kind, =0:u, =1: v.
const{
	MGPosition_list plist;
	mgTolSetWCZero wczeroSet(parameter_error());//Set&save the error.

	//1. outer boundaries.
	std::vector<UniqueCurve> curves=outer_boundary_param();
	int n=(int)curves.size(),i;
	for(i=0; i<n; i++){
		MGCParam_list list=curves[i]->isect_1D(x,kcod);
		MGCParam_list::iterator ci=list.begin(), ciend=list.end();
		while(ci!=ciend){
			MGPosition uv=curves[i]->eval(*ci++);
			plist.append(*this, uv);
		}
	}

	//2. inner boundaries.
	int n_inner=number_of_inner_boundaries();
	for(int j=0; j<n_inner; j++){
		std::vector<UniqueCurve> curvesj=inner_boundary_param(j);
		n=(int)curvesj.size();
		for(i=0; i<n; i++){
			MGCParam_list list=curvesj[i]->isect_1D(x,kcod);
			MGCParam_list::iterator ci=list.begin(), ciend=list.end();
			while(ci!=ciend){
				MGPosition uv=curvesj[i]->eval(*ci++);
				plist.append(*this, uv);
			}
		}
	}

	int other_cod=(kcod+1)%2;
	int np=(int)plist.size();
	std::vector<double> ivec(np);
	MGPosition_list::iterator pi=plist.begin();
	for(i=0; i<np; i++)	ivec[i]=(*pi++).ref(other_cod);
	std::sort(ivec.begin(),ivec.end());

	return ivec;
}

//Compute intersection points of this face's boundary(outer and inners) with
//face2. If intersection points are found and the boundary is a loop,
//the point's edge pointer(of this) will be stored in a member uvuv of uvuvs.
//uvuv[7] is the edge pointer. If the boundary is not a loop(that is, a perimeter of
//Surfaces), uvuv.sdim()==7 and an edge pointer is not returned.
//When uvuv.sdim()==8, the edge pointer of uvuv[7] is accessed through union mgEdgeP.
//uvuvs[i] is i-th intersection points.
int MGFace::isect_boundary(
	const MGFSurface& face2,
	MGPosition_list& uvuvs,
	//id1 and id2 are the ids of uvuv where this face's and f2's parameters
	//are to be stored in a member of uvuvs.
	//This face's (u,v) is stored in uvuv(id1) and (id1+1).
	//f2's (u,v) is stored in uvuv(id2) and (id2+1).
	//id2=0 if id1=2, and id2=2 if id1=0.
	int id1
)const{
	int inum=0;

	//Compute intersection points of the outer boundary of this face and face2.
	//Will be stored in uvuvs. The one member is uvuv(8):
	//0,1:(u,v) of face, 2,3:(u,v) of surf, 4-6:direction of isect.
	//7:for edge pointer at the intersection point.
	inum+=isect_outcurves(face2,uvuvs,id1);

	//Compute intersection points of inner boundaries of face fi and surf.
	int ibStrt;
	int nin=number_of_inner_boundaries(ibStrt);
	for(int k=0; k<nin; k++)
		inum+=isect_lpcurves(face2,ibStrt+k,uvuvs,id1);
	return inum;
}

//Compute intersection points between the boundary of iid-th inner boundary
//of this face and face2 to compute intersections of face with face2.
//Function's return value is the number of ip's obtained before appending
//into uvuv_list, may not be equal to the enlarged size of uvuv_list.
int MGFace::isect_incurves(
	const MGFSurface& face2,
	int iid,				//Loop id of this face.
	MGPosition_list& uvuv_list,	//intersection points will be appended.
		//One member in the list is of sdim 8,
		//(4,5,6) is the direction vector, and (7) is Edge pointer of the point.
	int id1			//id of uvuv(a member of uvuv_list).
		//uvuv(id1) for this face parameter uvuv(id2) for face2 parameter.
		//id2=0 if id1=2, and id2=2 if id1=0.
)const{
	int start_id;
	int num=number_of_inner_boundaries(start_id);
	return isect_lpcurves(face2,start_id+iid,uvuv_list,id1);
}

//"isect1_incr_pline" is a dedicated function of isect_start_incr, will get
// shortest parameter line necessary to compute intersection.
MGCurve* MGFace::isect_incr_pline(
	const MGPosition& uv,	//last intersection point.
	int kdt,				//Input if u=const v-parameter line or not.
							// true:u=const, false:v=const.
	double du, double dv,	//Incremental parameter length.
	double& u,				//next u value will be output
	double& v,				//next v value will be output
	int incr			//Incremental valuse of B-coef's id.
)const{
	const MGSurface* surf=surface();
	return surf->isect_incr_pline(uv,kdt,du,dv,u,v,incr);
}

//Compute intersection points between loop lp of this face and face2
//to compute intersections of face with face2.
//Function's return value is the number of ip's obtained before appending
//into uvuv_list, may not be equal to the enlarged size of uvuv_list.
int MGFace::isect_lpcurves(
	const MGFSurface& face2,		//srf!=null, face2=null.
	const MGLoop& lp,				//Loop id of this face.
	MGPosition_list& uvuv_list,	//intersection points will be appended.
		//One member in the list is of sdim 8,
		//(4,5,6) is the direction vector, and (7) is Edge pointer of the point.
	int id1			//id of uvuv(a member of uvuv_list).
		//uvuv(id1) for this face parameter uvuv(id2) for face2 parameter.
		//id2=0 if id1=2, and id2=2 if id1=0.
)const{
	int id2=2; if(id1==2) id2=0;
	double zero_angle=MGTolerance::angle_zero();

	mgEdgeP edgeP;//Union to access MGEdge* in MGPostion.
	MGLSPoint_vector lsps=lp.isect_binder(face2);
	int m=lsps.entries();
	for(int j=0; j<m; j++){
		const MGLSPoint& lsp=lsps[j];
		const MGPosition& suv=lsp.surface_param();
		MGPosition uvuv(8,suv,id2,0);
		double tb=lsp.binder_param();	//tb is binder edge's curve parameter value.
		const MGEdge* pedge=lsp.parameter_edge();
		double tp=pedge->param_pcell(tb);//tp is parameter edge's parameter.
		MGPosition uv=pedge->eval(tp);
		uvuv.store_at(id1,uv,0,2);
		MGVector N2=face2.normal(suv);
		MGVector N=normal(uv), V=pedge->eval_star(tp,1);
		double vn2angle=V.cangle(N2); if(vn2angle<.0) vn2angle*=-1.;
		if(vn2angle<=zero_angle){
			//vn2angle<=zero_angle means V and srf(or face2) are almost parallel.
			//double tps=pedge->param_s(), tpe=pedge->param_e();
			//if(tp>(tps+tpe)*.5) V*=-1.;
			//uvuv.store_at(4,V,0,3);
			uvuv.store_at(4,MGDefault::zero_vector(),0,3);
		}else uvuv.store_at(4,N2*(N*V)*N2,0,3);
		edgeP.pointer=pedge; uvuv(7)=edgeP.doubleV;
		uvuv_list.append(*this,face2,uvuv);
	}
	return m;
}

//Compute intersection points between loop lpid of this face and face2
//to compute intersections of face with face2.
//Function's return value is the number of ip's obtained before appending
//into uvuv_list, may not be equal to the enlarged size of uvuv_list.
int MGFace::isect_lpcurves(
	const MGFSurface& face2,		//srf!=null, face2=null.
	int lpid,				//Loop id of this face.
	MGPosition_list& uvuv_list,	//intersection points will be appended.
		//One member in the list is of sdim 8,
		//(4,5,6) is the direction vector, and (7) is Edge pointer of the point.
	int id1			//id of uvuv(a member of uvuv_list).
		//uvuv(id1) for this face parameter uvuv(id2) for face2 parameter.
		//id2=0 if id1=2, and id2=2 if id1=0.
)const{
	const MGLoop& lp=*(loop(lpid));
	return isect_lpcurves(face2,lp,uvuv_list,id1);
}

//Compute intersection points of outer boundary curves of this face 
//with face2 to compute intersections.
//Function's return value is the number of ip's obtained(appended)
//into uvuv_list, may not be equal to the enlarged size of uvuv_list.
int MGFace::isect_outcurves(
	const MGFSurface& face2,
	MGPosition_list& uvuv_list,	//intersection points will be appended.
		//One member in the list is of sdim 7,
		//and the last three elements are the ip direction vector.
	int id1			//id of uvuv(a member of uvuv_list).
		//uvuv(id1) for this face parameter uvuv(id2) for srf or face2 parameter.
		//id2=0 if id1=2, and id2=2 if id1=0.
)const{
	int id2=2; if(id1==2) id2=0;
	const MGSurface* srf1=surface();
	std::vector<UniqueCurve> perim(srf1->perimeter_num());
		//To save whole perimeter curve construction
		//at multiple processes of same perimeters.

	std::vector<MGFOuterCurve> cid=outer_curve();
	
	//Compute intersection points.
	int inum=0;
	std::vector<MGFOuterCurve>::iterator i=cid.begin(), ie=cid.end();
	for(; i!=ie; i++){
		if(i->is_loop()){
			//When loop.
			const MGLoop& lpi=*(i->loop());
			inum+=isect_lpcurves(face2,lpi,uvuv_list,id1);
		}else{
			//When perimeter.
			int pid=i->perimeter_id();
			if(perim[pid]==0) perim[pid].reset(srf1->perimeter_curve(pid));
			double t0,t1; i->range(t0,t1);
			MGTrimmedCurve crvi(*(perim[pid]), t0,t1);
			double error=MGTolerance::wc_zero();
			mgTolSetWCZero wczeroSet(error*.1);//Set&save the error.
			MGCSisects csilist=face2.isectFS(crvi);
			wczeroSet.restore();
			int m=csilist.entries(); inum+=m;
			for(int j=0; j<m; j++){
				auto cs=csilist.removeFirst();
				MGPosition uvuv(4,cs->param_surface(),id2,0);
				double t=cs->param_curve();	//Perimeter parameter.
				MGPosition uv=srf1->perimeter_uv(pid,t);
				uvuv.store_at(id1,uv,0,2);
				uvuv_list.prepend(*this,face2,uvuv);
			}
		}
	}
	return inum;
}

//Compute intersection lines, given end points of the i.l.
MGSSisects MGFace::isectEP(
	MGPosition_list& uvuv_list,	//End points list of the intersection.
		//On return, uvuv_list.size() will be 0.
	const MGFSurface& fsrf2		//The second surface.
) const{
	MGSSisects lst(this,&fsrf2);
	if(uvuv_list.size()==0) return lst;

	MGSSisect ssi;
	MGPosition_list::iterator uvuv_id;
	int obtained;
	const MGPlane* pl1=dynamic_cast<const MGPlane*>(surface());
	const MGPlane* pl2=dynamic_cast<const MGPlane*>(fsrf2.get_surface_pointer());
	if(pl1 && pl2){//When both are planes.
		while(uvuv_list.entries()){
			MGPosition uvuvS=uvuv_list.removeFirst();
			if(obtained=pl1->isect_startHPL(uvuvS,uvuv_list,*pl2,ssi,uvuv_id)){
				lst.append(std::move(ssi));
				if(obtained==3) uvuv_list.removeAt(uvuv_id);
			}
		}
	}else if(pl2){//When this is not a plane and srf2 is a plane.
		//Compute intersection lines using isect_with_plane.
		lst=isect_with_plane(uvuv_list,*pl2,fsrf2);
	}else if(pl1){//When this is a plane and srf2 is not a plane.
		MGPosition_list uvuvlist2;
		MGPosition_list::iterator i=uvuv_list.begin(), ie=uvuv_list.end();
		for(; i!=ie; i++){
			MGPosition uvuv2(*i);
			uvuv2(0)=(*i)[2]; uvuv2(1)=(*i)[3];
			uvuv2(2)=(*i)[0]; uvuv2(3)=(*i)[1];
			uvuvlist2.append(uvuv2);
		}
		lst=fsrf2.isect_with_plane(uvuvlist2,*pl1,*this);
		lst.exchange12();
		uvuv_list.clear();
	}else{//When both are not planes.
	//Compute intersection line using isect_with_surf.
		lst=isect_with_surf(uvuv_list,fsrf2);
	}
	return lst;
}

///Define curve division number when a curve crv be projected onto this MGFSurface.
///The result is used in prj2GetParamRange().
int MGFace::get_proj_divnum(const MGCurve& crv)const{
	int divnum=crv.intersect_dnum()+2;
	int first;
	int niloop=number_of_inner_boundaries(first);
	if(!niloop)
		return divnum;

	const MGBox& bpara=box_param();
	const MGInterval& urange=bpara[0];
	double u0=urange.low_point(), u1=urange.high_point();
	const MGInterval& vrange=bpara[1];
	double v0=vrange.low_point(), v1=vrange.high_point();
	double uspan=u1-u0, vspan=v1-v0;
	double uspan2=uspan, vspan2=vspan;
	for(int i=0; i<niloop; i++){
		const MGBox& ilbox=loop(first+i)->box();
		const MGInterval& urangei=ilbox[0];
		double u0i=urangei.low_point(), u1i=urangei.high_point();
		const MGInterval& vrangei=ilbox[1];
		double v0i=vrangei.low_point(), v1i=vrangei.high_point();
		double usi=u0i-u0;
		if(uspan>usi)
			uspan=usi;
		usi=u1-u1i;
		if(uspan>usi)
			uspan=usi;
		double vsi=v0i-v0;
		if(vspan>vsi)
			vspan=vsi;
		usi=u1-u1i;
		if(vspan>vsi)
			vspan=vsi;
	}
	int divn2=int(uspan2/uspan);
	int divn3=int(vspan2/vspan);
	if(divnum<divn2)
		divnum=divn2;
	if(divnum<divn3)
		divnum=divn3;
	return int(divnum*1.4);
}

///Execute polar-scaling to all the MGCurve and MGFace of this group.
///curve's (x,y) are updated. No other coordinates are unchanged.
///The updated result curve is always MGLBRep.
///For MGFace, the boundaries are polar-scaled.
///
///Rotation is performed from the angle range (angleBase,angle1) to
///(angleBase,angle2).
///That is, when angle1=angle2, no change is done.
///When angle2 is angleBase, all the data will lie on the straight of from origin to
///(cos(angleBase), sin(angleBase)).
///angle1-angleBase must be >MGTolerance::angle_zero().
///IF a member gel is not MGCurve nor MGFace, it is unchanged.
void MGFace::scalePolar(
	double angleBase,	///<base angle.
	double angle1,		
	double angle2
){
	make_outer_boundary();
	int nloop=number_of_loops();
	for(int i=0; i<nloop; i++){
		UniqueLoop& loopi=loop(i);
		int nedge=loopi->number_of_edges();
		for(int j=0; j<nedge; j++){
			MGEdge* edgej=loopi->edge(j);
			MGEdge* bedgej=edgej->binder_edge();
			if(bedgej){
				bedgej->set_extent_as_null();//Clear the old extent.
			}
			MGCurve* crvOld=edgej->base_curve();
			std::unique_ptr<MGLBRep> crvNew=crvOld->scalePolar(angleBase,angle1,angle2);
			edgej->set_extent(std::move(crvNew));
		}
	}
	invalidateBox();
}

///Execute polar-scaling to all the MGCurve and MGFace of this group.
///curve's (x,y) are updated. No other coordinates are unchanged.
///The updated result curve is always MGLBRep.
///For MGFace, the boundaries are polar-scaled.
///
///Rotation is performed from the angle range (angleBase,angle1) to
///(angleBase,angle2).
///That is, when angle1=angle2, no change is done.
///When angle2 is angleBase, all the data will lie on the straight of from origin to
///(cos(angleBase), sin(angleBase)).
///angle1-angleBase must be >MGTolerance::angle_zero().
///IF a member gel is not MGCurve nor MGFace, it is unchanged.
void MGGroup::scalePolar(
	double angleBase,	///<base angle.
	double angle1,		
	double angle2
){
	iterator i=begin(), ie=end(), ip;	
	for(; i!=ie; i=ip){
		ip=i;ip++;//Save i++ into ip.
		MGGel* gel = i->get();
		MGCurve* curve = dynamic_cast<MGCurve*>(gel);
		if(curve){
			std::unique_ptr<MGLBRep> lbi=curve->scalePolar(angleBase,angle1,angle2);
			iterator j=erase(i);
			insert(j,std::move(lbi));
			delete curve;
			continue;
		}
		MGFace* f=dynamic_cast<MGFace*>(gel);
		if(f){
			f->scalePolar(angleBase,angle1,angle2);
			continue;
		}
		MGGroup* g = dynamic_cast<MGGroup*>(gel);
		if(g)
			g->scalePolar(angleBase,angle1,angle2);
	}
}

///Rotate only the boundary of this face, but do not rotate the base surface.
///This is designed for the face of scalePolar().
void MGFace::rotateBoundary(const MGMatrix& mat){
	int nloop=number_of_loops();
	for(int i=0; i<nloop; i++){
		UniqueLoop& loopi=loop(i);
		int nedge=loopi->number_of_edges();
		for(int j=0; j<nedge; j++){
			MGEdge* edgej=loopi->edge(j);
			MGEdge* bedgej=edgej->binder_edge();
			if(bedgej){
				bedgej->set_extent_as_null();//Clear the old extent.
			}
			MGCurve& crvj=*(edgej->base_curve());
			crvj*=mat;
		}
	}
	invalidateBox();
}
