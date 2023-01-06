/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#include "StdAfx.h"
#include "mg/CParam_list.h"
#include "mg/Position_list.h"
#include "mg/CSisects.h"
#include "mg/SSisects.h"
#include "mg/Straight.h"
#include "mg/SurfCurve.h"
#include "mg/Curve.h"
#include "mg/LBRep.h"
#include "mg/TrimmedCurve.h"
#include "mg/Tolerance.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "topo/LEPoint.h"
#include "topo/LPoint.h"
#include "topo/Face.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
using namespace std;

//
//Implements MGFace Class.
//MGFace is an instance of MGCell.

/////// Constructor ///////

//Fundamental constructor.
//Construct a face from geometry of manifold dimension 2
//and the boundaries.
//The constructor takes the ownership of geo and MGLoop in boundaries.
//boundaries must be loops.
MGFace::MGFace(
	MGSurface* geo,
	const std::vector<UniqueLoop>& bndries
):MGCell(geo){
	assert(geo->manifold_dimension()==2);
	for(auto& li:bndries)
		m_boundaries.emplace_back(li->cloneWithParent(*this));
}
MGFace::MGFace(
	MGSurface* geo,
	std::vector<UniqueLoop>&& bndries
):MGCell(geo), m_boundaries(std::move(bndries)){
	assert(geo->manifold_dimension()==2);
	for(auto& loop : boundaries()){
		loop->set_parent(*this);
	}
}

//Conversion constructor from MGFSurface to MGFace.
MGFace::MGFace(const MGFSurface& surf)
:MGCell(surf.get_surface_pointer()->clone())
, m_box_param(surf.param_range()){
	copy_appearance(*(surf.object_pointer()));
	const MGFace* f=dynamic_cast<const MGFace*>(&surf);
	if(f)
		copy_all_boundaries(*f);
}

//Construct a face by copying boundaries(only parameter rep of the boundary)
//from argument boundaries.
MGFace::MGFace(
	const MGSurface& surf,
	const std::vector<MGLoop*>& boundaries
):MGCell(surf){
	std::vector<MGLoop*>::const_iterator i=boundaries.begin(),
		iend=boundaries.end();
	for(;i!=iend; i++){
		append_boundary((*i)->clone());
	}
}

//This form is to input a newed surface. The constructor takes the ownership
//of the surf.
MGFace::MGFace(
	UniqueSurface&& surf,
	std::vector<UniqueLoop>&& boundaries
):MGCell(surf.release()){
	for (auto& loopi : boundaries)
		append_boundary(loopi.release());
}

//Copy constructor.
MGFace::MGFace(const MGFace& face):MGFace(face,true){
}

//Move constructor.
MGFace::MGFace(MGFace&& face):MGCell(std::move(face)),
m_box_param(std::move(face.m_box_param)), m_boundaries(std::move(face.m_boundaries)){
	for(auto& loop : boundaries()){
		loop->set_parent(*this);
	}
}

//Constructor.
MGFace::MGFace(const MGFace& face, bool boundaryIsNecessary)
:MGCell(face), m_box_param(face.m_box_param){
	if(boundaryIsNecessary)
		copy_all_boundaries(face);
}

MGFace::~MGFace() {
	free_neighbours();
	for (UniqueLoop& li : m_boundaries)
		li.reset();
}

//Make a clone of the cell.
MGFace* MGFace::clone() const{
	return new MGFace(*this, true);
}
MGFace* MGFace::cloneWithoutBoundary() const{
	return new MGFace(*this, false);
}
/////// operator overload///////

// Faceに平行移動を行ないオブジェクトを生成する。
//Translation of the Face
MGFace MGFace::operator+ (const MGVector& v) const{
	MGFace face(*this);
	face+=v;
	return face;
}

// Faceに逆方向の平行移動を行ないオブジェクトを生成する。
//Translation of the Face
MGFace MGFace::operator- (const MGVector& v) const{
	MGFace face(*this);
	face-=v;
	return face;
}

//Faceのスケーリングを行い，Faceを作成する。
//Scaling of the Face by a double.
MGFace MGFace::operator* (double s) const{
	MGFace face(*this);
	face*=s;
	return face;
}

//Faceのスケーリングを行い，Faceを作成する。
//Scaling of the Face by a double.
MGFace operator* (double s, const MGFace& face){
	MGFace tface(face);
	tface*=s;
	return tface;
}

// 与えられた変換でFaceの変換を行い，Faceを作成する。
//Transformation of the Face by a matrix.
MGFace MGFace::operator* (const MGMatrix& mat) const{
	MGFace face(*this);
	face*=mat;
	return face;
}

// 与えられた変換によってトランスフォームをおこないFaceを生成する。
//Transformation of the Face by a MGTransf.
MGFace MGFace::operator* (const MGTransf& tr) const{
	MGFace face(*this);
	face*=tr;
	return face;
}

//Assignment.
//When the leaf object of this and cell2 are not equal, this assignment
//does nothing.
MGFace& MGFace::operator=(const MGFace& face){
	if(this==&face)
		return *this;

	MGCell::operator=(face);
	m_boundaries.clear();
	copy_all_boundaries(face);
	m_box_param=face.m_box_param;
	return *this;
}
MGFace& MGFace::operator=(MGFace && f2){
	if(this==&f2)
		return *this;

	MGCell::operator=(std::move(f2));
	m_box_param = std::move(f2.m_box_param);
	m_boundaries = std::move(f2.m_boundaries);
	for(auto& loop : m_boundaries){
		loop->set_parent(*this);
	}
	return *this;
}
MGFace& MGFace::operator=(const MGGel& gel2){
	const MGFace* gel2_is_this=dynamic_cast<const MGFace*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}
bool MGFace::operator<(const MGGel& gel2)const{
	const MGFace* f2=dynamic_cast<const MGFace*>(&gel2);
	if(f2)
		return operator<(*f2);
	return identify_type() < gel2.identify_type();
}

std::ostream& MGFace::toString(std::ostream& ostrm) const{
	ostrm<<"<<Face="<<(const MGGel*)this;
	ostrm<<", box_param="<<m_box_param;
	ostrm<<std::endl<<", boundaries(Loops) = "<<m_boundaries.size();
	const_iterator bs = m_boundaries.begin(), be = m_boundaries.end();
	int n = 0;
	if (bs == be) {
		ostrm<<std::endl;
	}
	else{
		ostrm<<"::" <<std::endl;
		ostrm<<(**bs).whoami()<<n++<<"=";
		ostrm<<(**bs);
		for(bs++; bs!=be; bs++){
			ostrm<<std::endl;
			ostrm<<n++<<"="<<(**bs);
		}
	}
	MGCell::toString(ostrm);
	ostrm<<"=Face>>"<<endl;
	return ostrm;
}

///////Member Function///////

//Generate arrow data of the tangent along u and v and the normal
//at the parameter value (u,v) of the surface.
//data[0] is the origin of the u-tangent arrow, data[1] is the top of the u-tangent arrow,
//data[2], [3] are two bottoms of u-tangent arrowhead.
//data[0], [4], [5], [6] are the points of v-tangent arrow.
//data[0], [7], [8], [9] are the points of v-tangent arrow.
void MGFace::arrow(double u,double v, MGPosition data[10])const{
	const MGSurface* f=surface();
	f->arrow(box(),u,v,data);
}

//Add a new loop to this face as aboundary.
//When the old inner loops that are outside the nloop will be removed from this.
//nloop can be inner or outer.
///nloop must be a newed MGLoop, and the ownership is transfered to this.
void MGFace::add_boundary(MGLoop* nloop){
	if(nloop->is_outer_boundary()){
		if(hasOuterBoundaryLoop())
			erase_boundary(int(0));
	}

	int nlps=number_of_loops();
	for(int i=nlps-1; i>=0; i--){
		if(!nloop->inside(loop(i)->start_point()))
			erase_boundary(i);
	}

	m_box_param.set_null();
	if(nloop->is_outer_boundary())
		prepend_boundary(nloop);
	else
		append_boundary(nloop);

	shrink_base_surface_to_knot();
}

//Obtain iterator of m_boundaries.
MGFace::const_iterator MGFace::boundaryIterator(const MGLoop* bnd) const{
	const_iterator itr=m_boundaries.begin(), itre=m_boundaries.end();
	for(; itr!=itre; itr++){ if((*itr).get()==bnd) break;}
	return itr;
}
MGFace::iterator MGFace::boundaryIterator(MGLoop* bnd){
	iterator itr=m_boundaries.begin(), itre=m_boundaries.end();
	for(; itr!=itre; itr++) {if((*itr).get()==bnd) break;}
	return itr;
}

//Append new one boundary to boundary vectors.
//Returned is the number of boudaries after appending.
//bound must be a newed MGLoop, and the ownership is transfered to this.
//*** append_boundary does not check validity with other loops
//(e.x. already existed loops will be outside the new boudanry bound).
//If the validity check is necessary, use add_boudanry().
int MGFace::append_boundary(MGLoop* bound){
	m_box_param.set_null();
	if(bound){
		m_boundaries.emplace_back(bound);
		bound->set_parent(*this);
	}
	invalidateBox();
	return (int)m_boundaries.size();
}

//Prepend new one boundary to boundary vectors.
//Returned is the number of boudaries after prepending.
int MGFace::prepend_boundary(MGLoop* bound){
	m_boundaries.emplace(m_boundaries.begin(), bound);
	bound->set_parent(*this);
	invalidateBox();
	return (int)m_boundaries.size();
}

///Copy all boundaries of cellin into this(but does not copy own binder cell relation),
///and the copied boundarys' MGPCell* association  of the original(key) and the new one into cmap.
//1st(key) is original MGPCell*(MGEdge) of boundaries(m_boundaries(MGLoop's)).
//2nd is copyied new.
void MGFace::copy_all_boundaries(
	const MGCell& cellin,
	std::map<const MGPCell*, MGPCell*>* cmap
){
	assert(dynamic_cast<const MGFace*>(&cellin));
	const MGFace* f = dynamic_cast<const MGFace*>(&cellin);
	for(auto& loopi:f->m_boundaries){
		m_boundaries.push_back(loopi->cloneWithParentAndMap(*this, *cmap));
	}
	invalidateBox();
}

//Erase i-th boundary.
//erase_boundary remove from this cell's bounary and destruct the boundary.
void MGFace::erase_boundary(int i){
	assert(i<number_of_boundaries());
	iterator itr = m_boundaries.begin();
	std::advance(itr, i);
	erase_boundary(itr);
}
void MGFace::erase_boundary(iterator i){
	m_boundaries.erase(i);
	invalidateBox();
}

//Erase specified boundary.
//erase_boundary remove from this cell's bounary and destruct the boundary.
void MGFace::erase_boundary(MGLoop* bnd){
	iterator itr = m_boundaries.begin();
	iterator itre = m_boundaries.end();
	for(; itr!=itre; itr++){
		if(bnd==(*itr).get()){
			m_boundaries.erase(itr);
			break;
		}
	}
	invalidateBox();
}

//Obtain the center parameter value of this cell.
MGPosition MGFace::center_param() const{
	MGPosition cntr;
	const_iterator bb = m_boundaries.begin(), be = m_boundaries.end();
	if(bb==be) cntr = m_extent->center_param();//If does not have boundary.
	else{
		//This manifold dimension is always greater than 1.
		const MGComplex* bnd = dynamic_cast<const MGComplex*>(bb->get());
		//Barycenter of the vertices of the first boundary is the center.
		cntr = bnd->center();
	}
	return cntr;
}

//Get the clone of this as a MGFace.
//If this is MGSurface, it is converted to MGFace.
MGFace* MGSurface::clone_as_face()const{
	return new MGFace(copy_surface());
}

///Free specified boundary(bound) from a member of parent cell's boundaries.
///Return MGComplex* bound if freed normally.
///If bound was not a member of the boundaries, return 0.
///Only free, does not destruct the boundary.
MGComplex* MGFace::free_boundary(const MGComplex* bound){
	iterator itr = m_boundaries.begin();
	iterator itre = m_boundaries.end();
	for(; itr!=itre; itr++){
		if(bound==(*itr).get()){
			MGLoop* bn = itr->release();
			m_boundaries.erase(itr);
			invalidateBox();
			return bn;
		}
	}
	return 0;
}

///If this had boundary binders, free them. As the result this
///will have no neighbours.
void MGFace::free_neighbours(){
	for(auto& loopi: m_boundaries){
		for(auto& edgej:loopi->pcells()){
			MGEdge* e = dynamic_cast<MGEdge*>(edgej.get());
			e->resetBinder();
		}
	}
}

//Obtain all the boundary curves(world coordinates representation)
//of the face.
//That is, all of the outer boundaries and all of the inner boundaries.
std::vector<UniqueCurve> MGFace::face_boundaries()const{
	std::vector<UniqueCurve> crvs=outer_boundary();
	int n=number_of_inner_boundaries(),i;
	for(i=0; i<n; i++){
		std::vector<UniqueCurve> crvsi=inner_boundary(i);
		std::move(crvsi.begin(), crvsi.end(), std::back_inserter(crvs));
	}
	return crvs;
}

///Get all the MGPCell* of the all the boundaries of this.
std::vector<const MGPCell*> MGFace::getBoundaryPcells()const{
	std::vector<const MGPCell*> pcells;
	for(auto& loopi: m_boundaries){
		for(auto& edgej:loopi->pcells()){
			const MGEdge* e = dynamic_cast<const MGEdge*>(edgej.get());
			pcells.push_back(e);
		}
	}
	return pcells;
}

//Return box of the parameter space of the face.
//After trimmed one.
const MGBox& MGFace::box_param() const{
	if(m_box_param.is_null())
		compute_box_param();
	return m_box_param;
}

//Compute closest point from a point.
//Returned is the parameter value of the face that is closest to point.
MGPosition MGFace::closest(const MGPosition& point) const{
	const MGSurface* srf=surface();
	MGPosition P,P1; double dist,dist1; MGPosition uv;
	MGVector dif;
	//(P,uv,dist) will be the colsest pair of the face to point.
	//uv:parameter value of P, dist:distance between P and point.

	//First, compute closest candidates on the surface.
	int n=perp_one(point,uv);

	std::vector<UniqueCurve> prmtr=outer_boundary();
	if(n==0){P=prmtr[0]->start_point(); uv=MGPosition();}
	else P=eval(uv);
	dif=(P-point); dist=dif%dif;

	//Second, compute closest candidates on the boundaries.
	//	2.1 Outer boundary.
	int pnum,i;
	pnum=(int)prmtr.size();
	for(i=0; i<pnum; i++){
		 double t=prmtr[i]->closest(point);
		 P1=prmtr[i]->eval(t); dif=P1-point; dist1=dif%dif;
		 if(dist1<dist){dist=dist1; P=P1; uv=MGPosition();}
	}

	//	2.2 Inner boundaries.
	int m=number_of_inner_boundaries();
	for(int j=0; j<m; j++){
		std::vector<UniqueCurve> prmtrj=inner_boundary(j);
		pnum=(int)prmtrj.size();
		for(i=0; i<pnum; i++){
			 double t=prmtrj[i]->closest(point);
			 P1=prmtrj[i]->eval(t); dif=P1-point; dist1=dif%dif;
			 if(dist1<dist){dist=dist1; P=P1; uv=MGPosition();}
		}
	}

	if(uv.sdim()==0) uv=param(P);
	return uv;
}

//Compute the closest point on all the perimeters of the surface.
//The point is returned as the parameter value (u,v) of this surface.
MGPosition MGFace::closest_on_boundary(const MGStraight& sl)const{
	int nloop=number_of_loops();
	if(!nloop)
		return surface()->closest_on_boundary(sl);

	double distance=-1.;
	MGLEPoint closest;
	for(int i=0; i<nloop; i++){
		const UniqueLoop& li=loop(i);
		li->closest_world(sl,closest,distance);
	}
	return closest.eval();
}

//compute box of the cell in m_box.
//Currently this does not compute correct box, compute m_extent box.
void MGFace::compute_box(MGBox& bx) const{
	const MGSurface* srf=surface();
	if(srf)
		bx =MGBox(srf->box_limitted(box_param()));
	else{
		bx = MGBox();
	}
}

//Compute parameter range box.
void MGFace::compute_box_param() const{
	if(hasOuterBoundaryLoop()){
		const UniqueLoop& lp=loop(int(0));
		m_box_param=lp->box();
	}else{
		const MGSurface* srf=surface();
		if(srf)
			m_box_param=srf->param_range();//Compute parameter range box.
		else{
			int n=number_of_boundaries();
			if(n){
				const UniqueLoop& lp=loop(int(0));
				m_box_param=lp->box();
			}
		}
	}
}

//set box as null(to set the box as initial)
void MGFace::invalidateBox()const{
	MGCell::invalidateBox();
	m_box_param.set_null();
}

//Test if directions of parameter curve and world curve of the face boundary
//is equal or not. This function can be used to test the pair of
//the output of outer_boundary() and outer_boundary_param(), or the pair of
//inner_boundary() and inner_boundary_param().
//Return is:
//true if equal direction, false if opposite direction.
bool MGFace::equal_direction(
	const std::vector<UniqueCurve>& wcurves,
		//output of outer_boundary() or inner_boundary().
	const std::vector<UniqueCurve>& pcurves,
		//output of outer_boundary_param() or inner_boundary_param().
	int i
		//id of the curve in wcurves and pcurves to test the direction.
		) const{
	assert(size_t(i)<wcurves.size() && wcurves.size()==pcurves.size());

	MGSurfCurve scrv(*(surface()),*(pcurves[i]));
	double t1=(scrv.param_s()+scrv.param_e())*0.5;
	double t2=(wcurves[i]->param_s()+wcurves[i]->param_e())*0.5;
	return (scrv.direction(t1))%(wcurves[i]->direction(t2))>0.;
}

//Evaluate.
//Input parameter value is not checked if it is in_range() or not.
//Even if it is not in_range(), surface evaluation will be executed.
MGVector MGFace::eval(double u, double v,	//Face parameter value(u,v)
			  int ndu, int ndv) const//Order of derivative.
{
	const MGSurface* srf=surface();
	return srf->eval(u,v,ndu,ndv);
}
MGVector MGFace::eval(const MGPosition& uv,	//Face parameter value(u,v)
			  int ndu, int ndv) const//Order of derivative.
{
	const MGSurface* srf=surface();
	return srf->eval(uv,ndu,ndv);
}

//Get inner_aboundary loops included in the input box.
std::vector<const MGLoop*> MGFace::get_inner_boundary_loops(
	const MGBox& uvbox
)const{
	std::vector<const MGLoop*> lps;
	int is;
	int n=number_of_inner_boundaries(is);
	if(n==0) return lps;
	double u0=uvbox[0].low_point(), u1=uvbox[0].high_point();
	double v0=uvbox[1].low_point(), v1=uvbox[1].high_point();
	for(int i=0; i<n; i++){
		const MGLoop* lpi=loop(i+is).get();
		const MGBox& bx=lpi->box();
		const MGInterval& urange=bx[0];
		if(urange<u0) continue;
		if(urange>u1) continue;
		const MGInterval& vrange=bx[1];
		if(vrange<v0) continue;
		if(vrange>v1) continue;
		lps.push_back(lpi);
	}
	return lps;
}

//Test if this face has boundary loops or not in the specified box.
//If this has one, return true.
bool MGFace::hasLoop(const MGBox& uvbox) const{
	int n=number_of_loops();
	if(n==0) return false;
	double u0=uvbox[0].low_point(), u1=uvbox[0].high_point();
	double v0=uvbox[1].low_point(), v1=uvbox[1].high_point();
	for(int i=0; i<n; i++){
		const MGBox& bx=loop(i)->box();
		const MGInterval& urange=bx[0];
		const MGInterval& vrange=bx[1];
		if(urange>=u0 && urange<=u1 && vrange>=v0 && vrange<=v1) return true;
	}
	return false;
}

//Test if this face has the outer boundary loop instead of perimeter boundary
//loops. If this has the outer boundary loop and has not perimeter boundary
//loops, return true.
bool MGFace::hasOuterBoundaryLoop() const{
	int n=number_of_boundaries();
	if(!n)
		return false;
	if(loop(int(0))->is_outer_boundary())
		return true;
	return false;
}

//Test if this face has perimeter boundary loops or not.
//If this has one, return true.
bool MGFace::hasPerimeterBoundaryLoop() const{
	int n=number_of_boundaries();
	if(!n) return false;
	if(loop(int(0))->is_perimeter_boundary()) return true;
	return false;
}

//Obtain i-th inner_boundary curves(world coordinates representation)
//of the face. Let the output of inner_boundary(i) be wcurves and
//of inner_boundary_param(i) be pcurves, then wcurves[j] corresponds
//to pcurves[j] one to one. Number of inner_boundary can be obtained
//by the function number_of_inner_boundaries().
std::vector<UniqueCurve> MGFace::inner_boundary(int i)const{
	int start;
	int nb=number_of_inner_boundaries(start);
	std::vector<UniqueCurve> crvs;
	if(nb && (i<nb)){
		int j=start+i;
		const UniqueLoop& lp=loop(j);
		crvs=lp->curves_world();
	}
	return crvs;
}

//Obtain i-th inner_boundary curves(parameter space representation)
//of the face. Let the output of inner_boundary(i) be wcurves and
//of inner_boundary_param(i) be pcurves, then wcurves[j] corresponds
//to pcurves[j] one to one. Number of inner_boundary can be obtained
//by the function number_of_inner_boundaries().
std::vector<UniqueCurve> MGFace::inner_boundary_param(int i)const{
	int start;
	std::vector<UniqueCurve> crvs;
	int nb=number_of_inner_boundaries(start);
	if(nb && (i<nb)){
		int j=start+i;
		const UniqueLoop& lp=loop(j);
		crvs = lp->curves();
	}
	return crvs;
}

//Test if (u,v) is inside inner boundary. inside means not on the
//boundary and not included inside the face.
//If true is returned, the id of m_boundaies is returned.
//Function's return value is:
//  0:outside of all the inner loops(not on the loop)
//  1:unknown
//  2:inside an inner loop(not on the loop), and the loop id is returned in id.
// otherwise:on the loop(int(MGEdge* of parameter edge))will be returned, 
//			 and the loop id is returned in id.
size_t MGFace::inside_inner_boundary(
	const MGPosition& uv, int& id
)const{
	int start;
	int n=number_of_inner_boundaries(start);
	for(int i=0; i<n; i++){
		id=start+i;
		int in_the_loop=loop(id)->inside(uv);
		if(!in_the_loop)
			return 2;
		if(in_the_loop>2)
			return in_the_loop;
	}
	return 0;
}

//Test if (u,v) is inside the outer boundary.
//Inside the outer boundary means that inside outer_boudary_param() or not.
//***Caution***
//(1)This must not be used for faces that do not have perimeter or outer boundary
//loop.
//(2)inside_outer_boundary does not check about the inner loops.
//
//Function's return value is:
//  0:outside the outer boundary(not on a loop)
//  1:unknown
//  2:inside the outer boundary(not on a loop)
// otherwise:on the outer boundary loop(including perimeter loops)
size_t MGFace::inside_outer_boundary(
	const MGPosition& uv
)const{
	std::vector<const MGLoop*> loops;
	extract_loops(loops);
	return inside_outer_loop(uv,loops,surface());
}

//Test if parameter value (u,v) is in the range of the face parameter.
bool MGFace::in_range(double u, double v)const{
	return in_range(MGPosition(u,v));
}
bool MGFace::in_range(const MGPosition& uv)const{
	if(box_param()<<uv) return false;
	else if(!number_of_boundaries())
		return true;
	else{
		int dummy;
		size_t in=inside_inner_boundary(uv,dummy);
		if(in>2)
			return true;
		else if(in==2)
			return false;
		else if(no_outer_boundaries())
			return true;
		else
			return inside_outer_boundary(uv)>=2;
	}
}

//Test if (u,v) is inside the face.
//Function's return value is:
//  0:outside the face.
//  <0:(u,v) is on an inner boundary, and abs(return code) is the loop id.
//  1:unknown.
//  2:inside the face, not on a boundary.
//  4:(u,v) is on the outer boundary(including perimeter loops).
//  >=10: (u,v) is on a perimeter, (10+perimeter number) will be returned.
int MGFace::in_range_with_on(
	const MGPosition& uv
)const{
	if(box_param()<<uv)
		return 0;

	int loop_id;
	size_t in=inside_inner_boundary(uv,loop_id);
	if(in>2)
		return -int(loop_id);
	else if(in==2)
		return 0;

	if(no_outer_boundaries()){
		const MGSurface* srf=surface();
		if(srf)
			return srf->in_range_with_on(uv);
		else
			return  2;
	}

	in=inside_outer_boundary(uv);
	if(in>2)
		return 4;
	else
		return int(in);
}

//Negate the face.
void MGFace::negate(){
	MGCell::negate();
	if(!m_box_param.is_null())
		m_box_param=MGBox(2,m_box_param,1,0);
		//Exchange parameter range of(u,v).
}

//Negate the boundaries.
void MGFace::negate_boundary(){
	//Negate each boudary.
	iterator bitr = m_boundaries.begin(), bitrend = m_boundaries.end();
	for(; bitr!=bitrend; bitr++)
		(*bitr)->negate_as_boundary(this);

	//Sort the boundaries.
	auto compareComplx=[](const UniqueLoop& b1, const UniqueLoop& b2){return (*b1)<(*b2); };
	std::sort(m_boundaries.begin(), m_boundaries.end(), compareComplx);
}

//Test if no outer boundary except the surface perimeters.
//That is, test if the following two conditions are satisfied:
//         1. no perimeter boundaries.
//         2. no outer boundary.
bool MGFace::no_outer_boundaries()const{
	int n=number_of_boundaries();
	if(!n) return true;
	const UniqueLoop& lp_first=loop(int(0));
	if(lp_first->is_perimeter_boundary()) return false;
	if(lp_first->is_outer_boundary()) return false;
	return true;
}

//Get number of inner boundaries.
//Returned i is the id of the first inner boundary loop if inner boundaries
//exist.
int MGFace::number_of_inner_boundaries(int& i)const{
	int n=number_of_boundaries();
	i=0;
	while(i<n && loop(i)->is_perimeter_boundary()) i++;
	while(i<n && loop(i)->is_outer_boundary()) i++;
	int j=i;
	while(j<n && loop(j)->is_inner_boundary()) j++;
	return j-i;
}

//Compute number of active loops.
int MGFace::number_of_loops()const{
	int n=number_of_boundaries();
	int i=0;
	while(i<n && loop(i)->is_perimeter_boundary()) i++;
	while(i<n && loop(i)->is_outer_boundary()) i++;
	while(i<n && loop(i)->is_inner_boundary()) i++;
	return i;
}

//Get number of perimeter boundary loop.
int  MGFace::number_of_perimeter_boundaries()const{
	int n=number_of_boundaries();
	int i=0;
	while(i<n && loop(i)->is_perimeter_boundary()) i++;
	return i;
}

//Test if a point P is on the face.
//Returned is true if the point P is on the face.
//false if P was not on the face.
bool MGFace::on(const MGPosition& P,
		MGPosition& uv	//Parameter value of the face is returrned.
						//Even if P is not on the face, nearest point
						//parameter value will be returned.
		) const{
	uv=closest(P);
	return P==eval(uv);
}

//Test if input (u,v) is parameter value on a perimeter of the base surface.
//If u or v is on a perimeter, they will be updated to the perimeter value.
bool MGFace::on_a_perimeter(
	double& u, double& v,		//Surface parameter (u,v)
	int& perim_num	//if function returns true,
						//the perimete number is output.
)const{
	return surface()->on_a_perimeter(u,v,perim_num);
}

//Obtain parameter value of the face whose world coordinates are P.
MGPosition MGFace::param(const MGPosition& P)const{
	const MGSurface* srf=surface();
	MGPosition uv=srf->param(P);
	return range(uv);
}

// Return ending parameter value.
double MGFace::param_e_u()const{
	const MGBox& uvbox=box_param();
	return uvbox[0].high_point();
}
double MGFace::param_e_v()const{
	const MGBox& uvbox=box_param();
	return uvbox[1].high_point();
}

// Return starting parameter value of the base surface.
double MGFace::param_s_u()const{
	const MGBox& uvbox=box_param();
	return uvbox[0].low_point();
}
double MGFace::param_s_v()const{
	const MGBox& uvbox=box_param();
	return uvbox[1].low_point();
}

//Return the foot of the perpendicular straight line from P.
//Computation is done from the guess parameter value.
//Function's return value is whether point is obtained(true) or not(false).
bool MGFace::perp_guess(
	const MGPosition& P,		//Point
	const MGPosition& uvguess,	// guess parameter value of the shell
	MGPosition& uv				// Parameter value will be returned.
) const{
	const MGBox& pbox=box_param();	//Parameter range of the face.
	MGPosition uv0(pbox[0].low_point(), pbox[1].low_point());
	MGPosition uv1(pbox[0].high_point(), pbox[1].high_point());

	bool obtained=false;
	if(surface()->perp_guess(uv0,uv1,P,uvguess,uv)) obtained=in_range(uv);
	return obtained;
}

//Compute perpendicular points of a curve and the face, given
//guess starting paramter values.
//Function's return value is:
//   perp_guess=true if perpendicular points obtained,
//   perp_guess=false if perpendicular points not obtained,
bool MGFace::perp_guess(
	const MGCurve& curve,	//curve.
	const MGPosition& uvguess,	//Guess parameter value of the face.
	double tguess,			//Guess parameter value of the curve.
	MGPosition& uv,			//perpendicular point's parameter values of the shell
	double& t				//will be output.
) const{
	const MGBox& pbox=box_param();	//Parameter range of the face.
	MGPosition uv0(pbox[0].low_point(), pbox[1].low_point());
	MGPosition uv1(pbox[0].high_point(), pbox[1].high_point());

	bool obtained=false;
	MGPosition tuvguess(3), tuv;
	tuvguess(0)=tguess;tuvguess(1)=uvguess[0];tuvguess(2)=uvguess[1];
	if(surface()->perp_guess(uv0,uv1,curve,1.,-1.,tuvguess,tuv)){
		t=tuv[0]; uv=MGPosition(2,tuv,0,1);
		obtained=in_range(uv);
	}
	return obtained;
}

//指定点から最も近い、垂線の足とパラメータ値を返す。
//Return the foot of the perpendicular straight line from p that is 
//nearest to point p.
// Function's return value is whether point is obtained(1) or not(0)
int MGFace::perp_point (
	const MGPosition& p,		// 指定点(point)
	MGPosition& uv,		//Parameter value of the surface will be returned.
	const MGPosition* uvguess	// guess parameter value of surface
	) const
{
	MGPosition uv0(1.,1.), uv1(0.,0.);
	if(uvguess) return perp_guess(p,*uvguess,uv);
	else        return perp_one(p,uv);
}

//Compute perpendicular points on the face from a point P((x,y,z)).
//MGPosition uv in the MGPosition_list is:
//uv(0): u parameter, and uv(1): v parameter of the face.
//Generally number of uv are more than one.
MGPosition_list MGFace::perps(const MGPosition& P) const{
	MGPosition_list uvs;
	const MGSurface* srf=surface(); if(!srf) return uvs;
	uvs=srf->perps(P);
	remove_outside_param(uvs);
	return uvs;
}

//Compute the parameter value of the closest point from the straight to
//this object.
//sl is the eye projection line whose direction is from yon to hither, and if
//sl had multiple intersection points, The closest point to the eye will be selected.
MGPosition MGFace::pick_closest(const MGStraight& sl)const{
	const MGSurface* srf=surface();
	if(no_outer_boundaries())
		return srf->pick_closest(sl);

	MGCSisects ises=srf->isectSl(sl);
	double t=.0; MGPosition uv;
	MGCSisects::iterator is=ises.begin(), ie=ises.end();
	for(; is!=ie; is++){
		MGCSisect& csi=isectCast<MGCSisect>(is);
		const MGPosition& uv2=csi.param_surface();
		if(in_range(uv2)){
			double t2=csi.param_curve();
			if(uv.is_null()){
				t=t2; uv=uv2;
			}else{
				if(t2>t){
					t=t2; uv=uv2;
				}
			}
		}
	}
	if(!uv.is_null())
		return uv;

	//Second, compute closest candidates on the boundaries.
	return closest_on_boundary(sl);
}

//Dedicated function of range.
//Will check if point (u,v) is inside inner boundary or not.
//If inside an inner boundary, obtain the closest point to the boundary.
MGPosition MGFace::range_check_inner_boundary(
	const MGPosition& uv
)const{
	int loop_id;
	double len;
	size_t in=inside_inner_boundary(uv,loop_id);
	if(in>=2){
		const UniqueLoop& cloop=loop(loop_id);
		MGLEPoint lp=cloop->closest(uv,len);
		return cloop->eval(lp);
	}else{
		return uv;
	}
}

//Round the input parameter (u,v) of the face to the nearest point of
//the face parameter range.
//**********This algorithm should be improved later.**********
//(Especially about using MGClosest_to_curves)
MGPosition MGFace::range(const MGPosition& uv) const{
	const MGSurface* srf=surface();
	if(!srf)
		return uv;

	mgTolSetWCZero wczeroSet(parameter_error());//Set&save the error.
	MGPosition uv_new;
	if(!srf->in_range(uv)){
		if(no_outer_boundaries())
			uv_new=srf->range(uv);
		else
			uv_new=MGClosest_to_curves(uv,outer_boundary_param());
	}else if(!number_of_boundaries()) {
	//Now uv is in the parameter range of surface.	
		uv_new=uv;	//When no boundaries.
	}else if(no_outer_boundaries()){
		uv_new=range_check_inner_boundary(uv);
	}else if (inside_outer_boundary(uv))
	//When there exists perimeter boundary.	
		uv_new=range_check_inner_boundary(uv);
	else
		uv_new=MGClosest_to_curves(uv,outer_boundary_param());

	return uv_new;
}

//Test if this face has an inactive loop.
//If this has one, return true.
bool MGFace::hasInactiveLoop()const{
	int n=number_of_boundaries();
	for(int i=n-1; i>=0; i--){
		if(loop(i)->is_inactive(this))
			return true;
	}
	return false;
}

//Remove inactive loops from this face.
void MGFace::remove_inactive_loops(){
	int n=number_of_boundaries();
	for(int i=n-1; i>=0; i--){
		if(loop(i)->is_inactive(this))
			erase_boundary(i);
	}
}

//Remove parameter uv from uvs that is outside face parameter range.
void MGFace::remove_outside_param(MGPosition_list& uvs)const{
	MGPosition_list::iterator i=uvs.begin(), iend=uvs.end(), i1;
	while(i!=iend){
		i1=i; i1++;
		if(!in_range(*i)) uvs.removeAt(i);
		i=i1;
	}
}

//Get surface pointer.
MGSurface* MGFace::surface(){
	UniqueGeometry& cell=extent();
	return dynamic_cast<MGSurface*>(cell.get());
}
const MGSurface* MGFace::surface() const{
	const UniqueGeometry& cell=extent();
	return dynamic_cast<const MGSurface*>(cell.get());
}

//Obtain the closest point from point uv to vector of curves.
//MGClosest_to_curves does not change wc_zero, and so calling program of
//MGClosest_to_curves should change it if necessary.
MGPosition MGClosest_to_curves(
	const MGPosition& uv,				//Point.
	const std::vector<UniqueCurve>& curves	//vector of curves.
){
	int n=(int)curves.size(); if(n==0) return uv;
	double t=curves[0]->closest(uv);
	MGPosition uv_min=curves[0]->eval(t), uv1;
	MGVector dif=uv-uv_min;
	double lenmin=dif%dif, len1;
	for(int i=1; i<n; i++){
		t=curves[i]->closest(uv);
		uv1=curves[i]->eval(t);
		dif=uv-uv1; len1=dif%dif;
		if(len1<lenmin){lenmin=len1; uv_min=uv1;}
	}
	return uv_min;
}
