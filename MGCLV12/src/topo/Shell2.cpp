/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include <iterator>
#include "mg/Box.h"
#include "mg/Position_list.h"
#include "mg/SSisect.h"
#include "mg/isects.h"
#include "mg/Tolerance.h"
#include "topo/LEPoint.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "topo/Shell.h"
#include "topo/HHisects.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGShell Class.

///////Member Function///////

//Compute the closest point from a point to this shell.
MGFPoint MGShell::closest(const MGPosition& point) const{
	int n=number_of_pcells(); if(!n) return MGFPoint();
	double* dist=new double[n], min, disti;
	int j;

	const_iterator i=pcell_begin(), ie=pcell_end();
	const MGFace* f=0;
	for(j=0; j<n; i++, j++){
		const MGFace* fi=dynamic_cast<const MGFace*>(i->get());
		if(!fi)
			continue;
		MGBox boxi=fi->box();dist[j]=disti=boxi.distance(point);
		if(!f){//For the first fi.
			f=fi; min=disti;
		}else{
			if(disti<min){
				f=fi; min=disti;
			}
		}
	}

	MGPosition uv=f->closest(point);
	min=(f->eval(uv)-point).len();
	const MGFace* fmin=f;
	if(!MGAZero(min)){ 

	for(i=pcell_begin(), j=0; j<n; i++, j++){
		const MGFace* fi=dynamic_cast<const MGFace*>(i->get());
		if(!fi || fi==f)
			continue;
		if(min<dist[j])
			continue;
		MGPosition uvi=fi->closest(point);
		disti=(fi->eval(uvi)-point).len();
		if(disti<min){
			min=disti; fmin=fi; uv=uvi;
			if(MGAZero(min))
				break;
		}
	}

	}

	delete[] dist;
	return MGFPoint(*fmin, uv);
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
MGisects MGShell::isect(const MGObject& obj2)const{
	auto curve = dynamic_cast<const MGCurve*>(&obj2);
	if(curve)
		return isect(*curve);
	auto surf2 = dynamic_cast<const MGFSurface*>(&obj2);
	if(surf2)
		return isectFS(*surf2);
	auto shell2 = dynamic_cast<const MGShell*>(&obj2);
	if(shell2)
		return isect(*shell2);
	return MGisects();
}

//Intersection of a shell and a curve.
MGCFisects MGShell::isect(const MGCurve& curve) const{
	MGCFisects is(&curve);
	const_iterator i=pcell_begin(), ie=pcell_end();
	for(; i!=ie; i++){
		const MGFace* fi=dynamic_cast<const MGFace*>(i->get());
		if(!fi)
			continue;
		MGCSisects list=fi->isect(curve);
		while(!list.empty())
			is.emplace_back(*(list.removeFirst()),*fi);
	}
	return is;
}

MGHHisects MGFace::isect(const MGShell& shell) const{
	MGHHisects hhi=shell.isectFS(*this);
	hhi.exchange12();
	return hhi;
}
MGHHisects MGSurface::isect(const MGShell& shell) const{
	MGHHisects hhi=shell.isectFS(*this);
	hhi.exchange12();
	return hhi;
}

//Intersection of two shells.
MGHHisects MGShell::isect(const MGShell& shell2)const{
	MGHHisects vec;
	int nf1=number_of_faces(), nf2=shell2.number_of_faces();
		//nf1,2=num of faces included in shell1,2.
	if(!nf1 || !nf2) return vec;
	if(!has_common(shell2)) return vec;

	MGHHisects lines2;
	const_iterator i=pcell_begin(), ie=pcell_end();
	for(; i!=ie; i++){
		const MGFace* fi=face(i); if(!fi) continue;
		const_iterator i2=shell2.pcell_begin(), ie2=shell2.pcell_end();
		for(; i2!=ie2; i2++){
			const MGFace* fi2=shell2.face(i2); if(!fi2) continue;
			MGSSisects ilst=fi->isect(*fi2);
			size_t ncrvs=ilst.size();
			for(size_t j=0; j<ncrvs; j++){
				auto ssi=ilst.releaseFront<MGSSisect>();
				lines2.append(fi,fi2,std::move(ssi));
			}
		}
	}

	for(MGHHisects::iterator hhi=lines2.begin(); hhi!=lines2.end(); hhi++){
		MGHHisect& is=isectCast<MGHHisect>(hhi);
		if(is.is_null()) continue;
		is.build_one(lines2);
		if(is.is_null()) continue;
		vec.emplace_back(hhi->release());
	}
	return vec;
}

//Intersection of a shell and a face.
//This shell's face is face1 in HHisect and face2 is face.
MGHHisects MGShell::isectFS(const MGFSurface& f) const{
	MGHHisects vec;

	int nf=number_of_pcells();//nf=num of faces included.
	if(!nf) return vec;
	if(!has_common(*f.object_pointer()))
		return vec;

	MGHHisects lines2;
	const_iterator i=pcell_begin(), ie=pcell_end();
	for(; i!=ie; i++){
		const MGFace* fi=face(i);
		if(!fi)
			continue;
		MGSSisects ilst=fi->isect(f);
		size_t ncrvs=ilst.size();
		for(size_t j=0; j<ncrvs; j++){
			std::unique_ptr<MGSSisect> ssi=ilst.releaseFront<MGSSisect>();
			lines2.append(fi,0,std::move(ssi));
		}
	}

	size_t nhhi = lines2.size();
	for(size_t j = 0; j<nhhi; j++){
		std::unique_ptr<MGHHisect> is = lines2.releaseFront<MGHHisect>();
		if(is->is_null()) continue;
		is->build_one(lines2);
		if(is->is_null()) continue;
		vec.push_back(std::move(is));
	}
	return vec;
}

//Get the partner point uvuv_id of the uvuv_list[j].
//If found, return true, if not, return false.
//Shell*(face or surface) version.
//The other one must not be shell, a surface or a face.
bool MGShell::isect_partner(
	const MGPosition& uvuv,//The point to get the partner.
	MGPosition_list* uvuv_list,
	int& j,		//Partner edge's face id(of uvuv_list) will be returned.
	MGPosition_list::iterator& uvuvItr
		//Partner point's iterator of uvuv_list[j] will be returned.
)const{
	if(uvuv.sdim()<8)
		return false;

	mgEdgeP edgP;	//Work area to access to MGEdge* in uvuv_list.
	const MGEdge* e; edgP.doubleV=uvuv[7]; e=edgP.pointer;
	const MGEdge* pe=e->first_partner();
	if(!pe) return false;//If e had no partnres.

	const MGCell* pf=pe->star();
	auto equalTopf = [pf](const UniqueCell& pcel){return pcel.get()==pf; };
	const_iterator fitr=std::find_if(pcell_begin(), pcell_end(), equalTopf);
	if(fitr==pcell_end())
		return false;

	j=(int)std::distance(pcell_begin(), fitr);
		//j is the id of uvuv_list of face pf.		

	MGPosition P=e->face()->eval(uvuv[0],uvuv[1]);

	//Find uvuv_list[j]'s uvuv that has the same edge and nearest to P.
	MGPosition_list::iterator itr=uvuv_list[j].begin(), itre=uvuv_list[j].end();
	double dlen=-1.;
	for(; itr!=itre; itr++){
		if(itr->sdim()<8) break;
			//We can break since elements of sdim()<8 are set last in the array.
		const MGEdge* et; edgP.doubleV=itr->ref(7); et=edgP.pointer;
		if(et==pe){
			MGPosition Q=et->face()->eval((*itr)[0],(*itr)[1]);
			if(dlen<0.){
				dlen=(P-Q).len(); uvuvItr=itr;
			}else{
				double dlen2=(P-Q).len();
				if(dlen>dlen2){
					dlen=dlen2; uvuvItr=itr;
				}
			}
		}
	}
	if(dlen<0.) return false;
	return true;
}

//Get the partner point uvuv_id of the uvuv_list[j].
//If found, return true, if not, return false.
//Shell*shell intersection version.
bool MGShell::isect_partner(
	const MGPosition& uvuv,//The point to get the partner.
	MGPosition_list* uvuv_list,
	int& i,
	int& j,		//Partner edge's face id(of uvuv_list) will be returned.
	MGPosition_list::iterator& uvuvItr
		//Partner point's iterator of uvuv_list[i+nf1*j] will be returned.
)const{
	if(uvuv.sdim()<8)
		return false;

	mgEdgeP edgP;	//Work area to access to MGEdge* in uvuv_list.
	const MGEdge* e; edgP.doubleV=uvuv[7]; e=edgP.pointer;
	const MGEdge* pe=e->first_partner();
	if(!pe) return false;//If e had no partnres.
	const MGCell* pf=pe->star();
	const MGComplex* shl=pf->parent_complex();
		//Shell of pf and e->star() is the same.
	auto equalTopf = [pf](const UniqueCell& pcel){return pcel.get()==pf; };
	const_iterator fitr=std::find_if(shl->pcell_begin(), shl->pcell_end(), equalTopf);
	if(fitr==shl->pcell_end()) return false;

	int k, kp1;
		//Id of uvuv to evaluate the position(this shell or the other shell).
	if(shl==this){
		i=(int)std::distance(pcell_begin(), fitr);
		k=0; 
	}else{
		j=(int)std::distance(shl->pcell_begin(), fitr);
		k=2; 
	}
	kp1=k+1;
	MGPosition P=e->face()->eval(uvuv[k],uvuv[kp1]);
	
	//i+j*number_faces() is the id of uvuv_list of face pf.
	int id=i+j*number_of_faces();
	//Find uvuv_list[j]'s uvuv that has the same edge and nearest to P.
	MGPosition_list::iterator itr=uvuv_list[id].begin(), itre=uvuv_list[id].end();
	double dlen=-1.;
	for(; itr!=itre; itr++){
		const MGEdge* et; edgP.doubleV=itr->ref(7); et=edgP.pointer;
		if(et==pe){
			MGPosition Q=et->face()->eval((*itr)[k],(*itr)[kp1]);
			if(dlen<0.){
				dlen=(P-Q).len(); uvuvItr=itr;
			}else{
				double dlen2=(P-Q).len();
				if(dlen>dlen2){
					dlen=dlen2; uvuvItr=itr;
				}
			}
		}
	}
	if(dlen<0.) return false;
	return true;
}

//Test if a point is on the shell or not.
//Function's return value is true if the point is on the shell, and false if not.
//The point parameter of the shell is returned in fp if true is returned.
//If false is returned, the closest point of the shell will be returned in fp.
bool MGShell::on(const MGPosition& point,
	MGFPoint& fp				//Shell's point parameter value.
)const{
	fp=closest(point);
	return MGAZero((fp.eval()-point).len())!=0;
}

//Obtain perpendicular points of a shell from a point.
std::vector<MGFPoint> MGShell::perps(const MGPosition& point) const{
	double error=MGTolerance::wc_zero(); error*=4.;
	const_iterator i=pcell_begin(), ie=pcell_end();
	MGPosition_list Ps;
		//This is used for the same points to be not added in the answer. 
	std::vector<MGFPoint> fps;//return value.
	for(; i!=ie; i++){
		const MGFace* fi=face(i); if(!fi) continue;
		MGPosition_list uvs=fi->perps(point);
		MGPosition_list::iterator j=uvs.begin(), je=uvs.end(), jwk;
		for(;j!=je; j++){
			MGPosition P=fi->eval(*j);
			MGBox box(P,error);
			if(!Ps.in(box,jwk)){
				Ps.push_back(P);
				fps.push_back(MGFPoint(*fi,*j));
			}
		}
	}
	return fps;
}

//Obtain the projected curves of a curve onto the shell.
//The direction of the projection is along the vector vec if the vec is not
//NULL, and normal to the shell if the vec is NULL.
//Output of 'project' is two kind of curves:
//one is general world coordinate curves(iline() of the MGHHisect members of lines),
//and the other is (u,v) curves of the parameter space of the surfaces
//(uvline1() of the MGHHisect members of lines ).
//*** uvline2() of the MGHHisect members of lines is a deque of length zero
//(not used).
//Function's return value is the number of curves obtained.
//When <0 is returned, some internal error occured.
int MGShell::project(
	const MGCurve& crv,		//given world coordinate curve to project.
	MGHHisects& lines,
			//World coordinates (x,y,z) lines and (u,v) lines of the projected
			//curves will be returned in lines.
	const MGVector& vec	//projection direction.
			//if vec = NULL, then projection that is normal to the shell.
)const{
	MGHHisects lines2;
	const_iterator i=pcell_begin(), ie=pcell_end();
	for(; i!=ie; i++){
		const MGFace* fi=face(i); if(!fi) continue;
		std::vector<UniqueCurve> crv_uvs;		//uv projection curve.
		std::vector<UniqueCurve> crvs;		//projection curve.
		int ncrvs=fi->project(crv,crv_uvs,crvs,vec);
		if(ncrvs<0)
			continue;//since error detected.
		for(int j = 0; j < ncrvs; j++){
			lines2.emplace_back<MGHHisect>(crvs[j].release(), fi, crv_uvs[j].release());
		}
	}
	size_t nlines=lines2.size();
	if(!nlines) return 0;

	size_t k, nlines2=0;

	//Search hhisect closest to the start point of crv.
	MGPosition PS=crv.start_point();
	double dist=-1.;
	MGHHisects::iterator hhi=lines2.begin(), minhhi;
	for(k=0; k<nlines; k++, hhi++){
		auto& hhii=isectCast<MGHHisect>(hhi);
		MGVector diff=hhii.iline().start_point()-PS;
		double dist2=diff%diff;
		if(dist<0. || dist2<dist){
			minhhi=hhi;
			dist=dist2;
		}
	}

	MGHHisect is(std::move(isectCast<MGHHisect>(minhhi)));
	hhi=lines2.begin();
	for(k=0; k<nlines; k++,hhi++){
		auto& hhii=isectCast<MGHHisect>(hhi);
		if(k)
			is=std::move(hhii);
		if(is.is_null()) continue;
		is.build_one(lines2);
		if(is.is_null()) continue;
		lines.emplace_back(new MGHHisect(std::move(is)));
		nlines2++;
	}
	return (int)nlines2;
}
