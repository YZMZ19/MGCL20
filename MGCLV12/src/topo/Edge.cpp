/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/Box.h"
#include "mg/Curve.h"
#include "mg/Straight.h"
#include "mg/LBRep.h"
#include "mg/TrimmedCurve.h"
#include "mg/Surface.h"
#include "mg/SurfCurve.h"
#include "mg/Tolerance.h"
#include "topo/PVertex.h"
#include "topo/BVertex.h"
#include "topo/Loop.h"
#include "topo/Edge.h"
#include "topo/Face.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

//Implements MGEdge Class.
//MGEdge is a 1D minifold instance of MGCell.

/////// Constructor ///////

//void constructor.
MGEdge::MGEdge():m_equal_to_binder(0){
}

MGEdge::MGEdge(const MGEdge& e, bool boundaryIsNecessary)
: MGCell(e), m_equal_to_binder(0) {//We do not copy PCell and BCell data.
	if(boundaryIsNecessary)
		copy_all_boundaries(e);
}
//Copy constructor.
MGEdge::MGEdge(const MGEdge& e)
: MGEdge(e, true){
}

//Move constructor.
MGEdge::MGEdge(MGEdge&& eFrom)noexcept
:MGCell(std::move(eFrom)), m_equal_to_binder(0) {//We do not copy PCell and BCell data.
	for(int i = 0; i<2; i++){
		UniquePVertex& toV = m_vertex[i];
		toV = std::move(eFrom.m_vertex[i]);
		if(toV){
			toV->set_edge(this);
		}
	}
}

//Fundamental constructor.
//Construct an edge from geometry of manifold dimension 1
//and the boundaries.
//The constructor takes the ownership of geo and boundaries.
MGEdge::MGEdge(
	MGCurve* geo,
	MGPVertex* boundaries[2]
):MGCell(geo), m_equal_to_binder(0){
	for(int i=0; i<2; i++){
		m_vertex[i].reset(boundaries[i]);
		if(m_vertex[i])
			m_vertex[i]->set_edge(this);
	}
	assert(geo->manifold_dimension()==1);
}

//Make an edge with a binder of a boundary that has active start and end vertex
//(if the curve is not infinite straight line).
//The second form that input MGCurve* takes the ownership of the crv
//into the MGEdge, must not delete the object and the object must be
//newed one.
MGEdge::MGEdge(const MGCurve& crv):MGCell(crv), m_equal_to_binder(0){
}
MGEdge::MGEdge(MGCurve* crv):MGCell(crv), m_equal_to_binder(0){

}
MGEdge::~MGEdge() {
	free_neighbours();
}

///Set both end parameter vertices.
void MGEdge::setBothEnds(const MGInterval& range){
	const MGCurve& crv = *base_curve();
	MGInterval prange = crv.param_range();
	prange &= range;
	if(prange.finite_below())
		set_start(prange.low_point());
	if(prange.finite_above())
		set_end(prange.high_point());
}

//Make an edge of a boundary that has active start and end vertex.
//range is the parameter range of crv.
//The second form that input MGCurve* takes the ownership of the crv
//into the MGEdge, must not delete the object and the object must be
//newed one.
MGEdge::MGEdge(const MGCurve& crv, const MGInterval& range)
:MGCell(crv), m_equal_to_binder(0){
	setBothEnds(range);
}
MGEdge::MGEdge(MGCurve* crv, const MGInterval& range)
:MGCell(crv), m_equal_to_binder(0){
	setBothEnds(range);
}

//Make an edge with a binder of a boundary that has active start and end vertex.
MGEdge::MGEdge(
	const MGSurface&surf,//Parent surface of which this edge makes a boundary
	const MGCurve& pcrv, //Parameter curve of the surface surf.
	const MGInterval& prange,//param range of pcrv.
	const MGCurve& wcrv  //World coordinate curve of the surface surf.
						//wcrv will be trimmed by prange of pcrv.
):MGEdge(pcrv, prange){
	double t1=param_s(), t2=param_e(), s1,s2;
	wcrv.on(surf.eval(pcrv.eval(t1)),s1);
	wcrv.on(surf.eval(pcrv.eval(t2)),s2);
	resetBinder(new MGEdge(wcrv,MGInterval(MGInterval(s1),s2)));
}

//Make a clone of this edge.
//clone() does not copy the binder cell relation.
MGEdge* MGEdge::clone() const{
	MGEdge* e= new MGEdge(*this,true);
	return e;
}
MGEdge* MGEdge::cloneWithoutBoundary() const{
	MGEdge* e= new MGEdge(*this, false);
	return e;
}

///Clone this BCell, building the new partner membership of partners.
SharedBCell MGEdge::cloneWithPartnerMembers(
	std::vector<const MGPCell*>& newPartners
)const{
	SharedBCell binderEdge(new MGEdge(*this));
	binderEdge->m_partners = std::make_unique< PCellVec>(newPartners);
	return binderEdge;
}

//Assignment.
//does not change binder and partner relation,
//does not change parent complex.
MGEdge& MGEdge::operator=(const MGEdge& e){
	if(this==&e)
		return *this;

	MGCell::operator=(e);
	for(int i=0; i<2; i++){
		if(e.m_vertex[i]){
			if(m_vertex[i])
				m_vertex[i]->set_t(e.m_vertex[i]->t());
			else{
				m_vertex[i].reset(new MGPVertex(*(e.m_vertex[i]),this));
			}
		}else{
			m_vertex[i].reset();
		}
	}
	m_equal_to_binder=0;
	return *this;
}
MGEdge& MGEdge::operator=(const MGGel& gel2){
	const MGEdge* gel2_is_this=dynamic_cast<const MGEdge*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGEdge::operator==(const MGEdge& e2)const{
	return this==&e2;
}
bool MGEdge::operator==(const MGGel& gel2)const{
	const MGEdge* gel2_is_this=dynamic_cast<const MGEdge*>(&gel2);
	if(gel2_is_this)
		return operator==(*gel2_is_this);
	return false;
}
bool MGEdge::operator<(const MGEdge& e2)const{
	if(loop()== e2.loop())
		return MGCell::operator<(e2);

	MGTrimmedCurve crv1=trimmed_curve();
	MGTrimmedCurve crv2=trimmed_curve();
	return crv1<crv2;
}
bool MGEdge::operator<(const MGGel& gel2)const{
	const MGEdge* gel2_is_this=dynamic_cast<const MGEdge*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return identify_type()<gel2.identify_type();
}
bool MGEdge::operator!=(const MGGel& gel2)const{
	return !((*this)==gel2);
}
bool MGEdge::operator!=(const MGEdge& e2)const{
	return !((*this)== e2);
}

// Edge に平行移動を行ないオブジェクトを生成する。
//Translation of the Edge
MGEdge MGEdge::operator+ (const MGVector& v) const{
	MGEdge nEdge(*this);
	nEdge+=v;
	return nEdge;
}

// Edgeに逆方向の平行移動を行ないオブジェクトを生成する。
//Translation of the Edge
MGEdge MGEdge::operator- (const MGVector& v) const{
	MGEdge nEdge(*this);
	nEdge-=v;
	return nEdge;
}

//Edgeのスケーリングを行い，Edgeを作成する。
//Scaling of the Edge by a double.
MGEdge MGEdge::operator* (double s) const{
	MGEdge nEdge(*this);
	nEdge*=s;
	return nEdge;
}

//Edgeのスケーリングを行い，Edgeを作成する。
//Scaling of the Edge by a double.
MGEdge operator* (double s, const MGEdge& Edge){
	MGEdge nEdge(Edge);
	nEdge*=s;
	return nEdge;
}

// 与えられた変換でEdgeの変換を行い，Edgeを作成する。
//Transformation of the Edge by a matrix.
MGEdge MGEdge::operator* (const MGMatrix& mat) const{
	MGEdge nEdge(*this);
	nEdge*=mat;
	return nEdge;
}

// 与えられた変換によってトランスフォームをおこないEdgeを生成する。
//Transformation of the Edge by a MGTransf.
MGEdge MGEdge::operator* (const MGTransf& tr) const{
	MGEdge nEdge(*this);
	nEdge*=tr;
	return nEdge;
}

//Object transformation.
MGEdge& MGEdge::operator+=(const MGVector& v){
	MGCell::operator+=(v);
	m_equal_to_binder=0;
	return *this;
}
MGEdge& MGEdge::operator*=(double scale){
	MGCell::operator*=(scale);
	m_equal_to_binder=0;
	return *this;
}
MGEdge& MGEdge::operator*=(const MGMatrix& mat){
	MGCell::operator*=(mat);
	m_equal_to_binder=0;
	return *this;
}
MGEdge& MGEdge::operator*=(const MGTransf& tr){
	MGCell::operator*=(tr);
	m_equal_to_binder=0;
	return *this;
}

//Get after edge in the loop sequence.
//The aft_edge is the first neighbour edge.
const MGEdge* MGEdge::aft_edge(bool at_end,int* vertexID) const{
	const MGEdge* aft=0;
	int vid= at_end ? 1:0;

	if(m_vertex[vid]){
		std::vector<const MGPCell*> ptrs=m_vertex[vid]->partners();
		if(ptrs.size()){
			const MGPVertex* pv=static_cast<const MGPVertex*>(ptrs[0]);
			aft = dynamic_cast<const MGEdge*>(pv->star());
			if(vertexID){
				*vertexID= (pv->is_start_vertex()) ? 0:1;
			}
		}
	}
	return aft;
}
MGEdge* MGEdge::aft_edge(bool at_end,int* vertexID){
	const MGEdge* cthis=this;
	return const_cast<MGEdge*>(cthis->aft_edge(at_end,vertexID));
}

//Obtain binder edge pointer.
MGEdge* MGEdge::binder_edge() const{
	return dynamic_cast<MGEdge*>(binder().get());
}

//Obtain the center parameter value of this cell.
MGPosition MGEdge::center_param() const{
	double t=param_s();
	t+=param_e();
	t*=0.5;
	return MGPosition(1,&t);
}

//Bind this edge  to edge2(is a MGEdge). Both edges are parameter edges of faces.
//This cell is a pcell of a boundary of a face A,
//and cell2 is also is a pcell of a face B.
//That is, this cell is a part of a boundary of face A,
//and cell2 is a part of a boundary of face B.
void MGEdge::bind(MGEdge& edge2){
	m_equal_to_binder=edge2.m_equal_to_binder=0;
	MGPCell::bind(edge2);
}

//compute box of the cell in m_box.
void MGEdge::compute_box(MGBox& bx) const{
	MGTrimmedCurve tc = trimmed_curve();
	tc.compute_box(bx);
}

///Connect the start(id=0) or end(id=1) of this to the pvert's edge at pvert.
void MGEdge::connect(int id, MGPVertex& pvert){
	MGEdge* Eofpvert=pvert.starEdge();
	assert(Eofpvert);
	int StartEnd=pvert.is_start_vertex() ? 0:1;
	connect_at_id(id,Eofpvert,StartEnd);
}

///Test if this has non null vertex at start or end(when atStart=false),
///and if does not have, make one.
void MGEdge::ensureHasVertex(int id) {
	assert(id == 1 || id == 0);
	if (!m_vertex[id]) {
		const MGCurve* bcrv = base_curve();
		double t = id ? bcrv->param_e():bcrv->param_s();
		setStartOrEnd(id, t);
	}
}

///Test if this has non null vertex at start and end,
///and if does not have, make ones.
void MGEdge::ensureHasVerticesAtBothEnds() {
	ensureHasVertex(0);
	ensureHasVertex(1);
}

//Connect the start(id1=0) or end(id1=1) of this to the start(id2=0) or
// the end(id2=1) of e2.
//If both edges of this and e2 are members of a complex, they must be
//the same.
void MGEdge::connect_at_id(int id1, MGEdge* e2, int id2){
	assert(e2->is_pcell());
	assert(id1<=1 && id2<=1);

	ensureHasVertex(id1);
	e2->ensureHasVertex(id2);

	MGComplex* complex=parent_complex();//is MGLoop.
	MGComplex* complex2=e2->parent_complex();
	assert(!complex || !complex2 || complex==complex2);
	//If both has complexes, they must be the same.

	if(complex && !complex2){
		if(id1==0)
			complex->prepend_pcell(e2);
		else
			complex->append_pcell(e2);
	}else if(!complex && complex2){
		if(id1==0)
			complex2->prepend_pcell(this);
		else
			complex2->append_pcell(this);
	}

	merge_2bcells(*m_vertex[id1], *e2->m_vertex[id2]);
}

///Copy all boundaries of cellin into this(but does not copy own binder cell relation),
//and register cellin's boundary MGPCell* association into cmap.
//1st(key) is original MGPCell* of boundaries of this,
//2nd is copyied new.
void MGEdge::copy_all_boundaries(
	const MGCell& cellin, 
	std::map<const MGPCell*, MGPCell*>* cmap//1st(key) is original and 2nd is copyied new.
){
	assert(dynamic_cast<const MGEdge*>(&cellin));
	const MGEdge* edge= dynamic_cast<const MGEdge*>(&cellin);
	for(int i=0; i<2; i++){
		UniquePVertex& toV = m_vertex[i];
		if (toV)
			toV.reset();
		const UniquePVertex& fromV = edge->m_vertex[i];
		if(fromV){
			toV = std::make_unique<MGPVertex>(*fromV, this);
			if(cmap)
				cmap->insert(std::make_pair(fromV.get(),toV.get()));
		}
	}
	invalidateBox();
}

//Return curve pointer of this edge.
MGCurve* MGEdge::base_curve(){
	UniqueGeometry& geo=extent();
	return dynamic_cast<MGCurve*>(geo.get());
}
const MGCurve* MGEdge::base_curve() const{
	const UniqueGeometry& geo=extent();
	return dynamic_cast<const MGCurve*>(geo.get());
}

//Return curve pointer cut by start and end parameter range.
//Output is newed curve object, must be deleted.
MGCurve* MGEdge::curve_limitted() const{
	const MGCurve* crv=base_curve();
	if(crv) return crv->copy_limitted(range());
	else return 0;
}

//Obtain this edge's parent loop's edge iterator as a member of the parent loop.
//This must have the parent loop, i.e. loop() must not null.
MGComplex::const_iterator MGEdge::edge_iterator()const{
	const MGLoop* lp=loop();
	MGComplex::const_iterator i=lp->pcell_begin(), ie=lp->pcell_end();
	const MGCell* enb=static_cast<const MGCell*>(this);
	for(; i!=ie; i++)
		if(i->get()==enb)
			break;
	return i;
}

MGComplex::iterator MGEdge::edge_iterator(){
	MGLoop* lp=loop();
	MGComplex::iterator i=lp->pcell_begin(), ie=lp->pcell_end();
	MGCell* enb=static_cast<MGCell*>(this);
	for(; i!=ie; i++)
		if(i->get()==enb)
			break;
	return i;
}

//Obtain edge number as a member of the parent loop.
//This must have the parent loop, i.e. loop() must not null.
int MGEdge::edge_num()const{
	const MGLoop* lp=loop();
	return lp->edge_num(this);
}

//Test if SurfCurve of the edge has equal direction to binder edge's direction.
//Returned is true if eaual, false if not.
bool MGEdge::equal_direction_to_binder()const{
	if(m_equal_to_binder)
		return (m_equal_to_binder==1);

	const MGEdge* bedge=binder_edge();
	if(!bedge)
		return false;	//If does not have binder.
	const MGSurface* srf=star_surface();
	if(!srf)
		return false;		//If does not have star surface. 

	m_equal_to_binder=srf->equal_direction(trimmed_curve(), bedge->trimmed_curve());
	return (m_equal_to_binder==1);
}

//Evaluate the nderiv's derivative at parameter t.
//Evaluate of the curve's data.
MGVector MGEdge::eval(double t, int nderiv)const{
	assert(base_curve());
	return base_curve()->eval(t,nderiv);
}

//Evaluation of the star curves of the edge at the point t.
//When nderi=0, get a position of the surface at the boundary point t.
//The star curve is SurfCurve(face's surface, edge's curve).
//(The star curve has the same world coordinate with the binder curve's, but
//their direction may be opposite. The star curve has always the same direction
//as the loop.)
MGVector MGEdge::eval_star(
	double t,		//Parameter value of this parameter edge's curve.
	int nderi	//Order of derivative.
)const{
	assert(face()->surface());//This edge must be a boundary member of a face.
	const MGCurve* crv=base_curve();
	if(crv){
		MGSurfCurve starCurve(*(face()->surface()),*crv);
		return starCurve.eval(t, nderi);
	}
	else return MGVector();
}

///Get extent geometry, may be null if this does not have extent.
const MGCurve* MGEdge::extentBC() const{
	return dynamic_cast<const MGCurve*>( extent().get());
}
MGCurve* MGEdge::extentBC(){
	return dynamic_cast<MGCurve*>(extent().get());
}

//Get the star face pointer.
const MGFace* MGEdge::face() const{
	return dynamic_cast<const MGFace*>(star());
}
MGFace* MGEdge::face(){
	return dynamic_cast<MGFace*>(star());
}

//Get the 1st partner edge of this edge.
const MGEdge* MGEdge::first_partner() const{
	const MGEdge* pedge=0;
	std::vector<const MGPCell*> cell=partners();
	if(cell.size()){
		pedge=dynamic_cast<const MGEdge*>(cell[0]);
	}
	return pedge;
}

//free start(i=0) or end(i=1) boundary's bindness.
void MGEdge::free_vertex(int i) {
	assert(i == 0 || i == 1);
	if (m_vertex[i])
		m_vertex[i]->resetBinder();

}

///If this had boundary binders, free them. As the result this
///will have no neighbours.
void MGEdge::free_neighbours(){
	free_vertex(0);
	free_vertex(1);
}

///Comparison of two MGEdges as MGPCells.
///is_less_than() defines the order of partners to store in MGBCell.
/// *****Currently definite MGEdged order not defined for MGBCell.******
bool MGEdge::is_less_than(const MGEdge& pe2)const {
	return operator<(pe2);
}
bool MGEdge::is_less_than(const MGPCell& pcel2)const{
	const MGEdge* e2 = dynamic_cast<const MGEdge*>(&pcel2);
	if(e2)
		return operator<(*e2);
	return identify_type()<pcel2.identify_type();
}

///test if parameter t is the one of the end point of the loop.
bool MGEdge::is_end_point(double t)const{
	double te=param_e();
	double error=parameter_error();
	return (te-error<=t && t<=te+error);
}

///test if parameter t is the one of the start point of the loop.
bool MGEdge::is_start_point(double t)const{
	double ts=param_s();
	double error=parameter_error();
	return (ts-error<=t && t<=ts+error);
}

//Connect this and e2.
//If start==true, start of this edge to end of e2;
//If start==false, end of this edge to start of e2;
//e2 must be a newed object, and the ownership is transfered to the system.
void MGEdge::join(bool start, MGEdge* e2){
	int i=1,j=0;
	if(start){ i=0; j=1;}
	connect_at_id(i,e2,j);
}

//Return parent loop pointer.
const MGLoop* MGEdge::loop()const{
	const MGComplex*  cmp=parent_complex();
	return dynamic_cast<const MGLoop*>(cmp);
}
MGLoop* MGEdge::loop(){
	MGComplex* cmp=parent_complex();
	return dynamic_cast<MGLoop*>(cmp);
}

//Negate the boundary.
void MGEdge::negate_boundary(){
	m_vertex[0].swap(m_vertex[1]);
	MGCurve* crv=base_curve();
	if(!crv)
		return;

	for(int i=0; i<2; i++){
		UniquePVertex& v=m_vertex[i];
		if(v){
			double param=v->t();
			v->set_t(crv->negate_param(param));
		}
	}
}

//Compute the mid point of this edge.
//Mid point is the point of the paramete mid=(param_s()+param_e())*.5
MGPosition MGEdge::mid_point()const{
	double t=(param_s()+param_e())*.5;
	return eval(t);
}

//Negate the direction of the cell.
void MGEdge::negate(){
	MGCell::negate();
	m_equal_to_binder*=-1;
}

//Obtain all the neighbours.
//The neighbours do not contain this cell except when this cell is
//connected to this cell itself(closed cell).
std::vector<const MGCell*> MGEdge::neighbours() const{
	std::vector<const MGCell*> nbrs;
	for(int i=0; i<2; i++){	//Loop for the two m_vertex[].
		std::vector<const MGPCell*> prtnrs;
		if(m_vertex[i])
			prtnrs=m_vertex[i]->partners();
		std::vector<const MGPCell*>::iterator
			prtnr_now=prtnrs.begin(), prtnr_end=prtnrs.end();
		for(;prtnr_now!=prtnr_end; prtnr_now++){
			const MGCell* nbr=(*prtnr_now)->star();
			if(nbr){
				if(std::find(nbrs.begin(),nbrs.end(),nbr)==nbrs.end())
					nbrs.push_back(nbr);
			}
		}
	}
	return nbrs;
}

//Get the perimeter number where this edge is on.
//If this is not on any perimeter, -1 will be returned.
int MGEdge::surface_perimeter() const{
	const MGFace* f=face(); if(!f) return -1;
	const MGSurface* sf=f->surface(); if(!sf) return -1;
	return surface_perimeter(*sf);
}
int MGEdge::surface_perimeter(const MGFace& face) const{
	const MGSurface* sf=face.surface(); if(!sf) return -1;
	return surface_perimeter(*sf);
}
int MGEdge::surface_perimeter(const MGSurface& sf) const{
	int pnum;
	if(sf.on_perimeter(*(base_curve()),pnum)) return pnum;
	return -1;
}

ostream& MGEdge::toString(ostream& ostrm) const{
	std::string me = whoami();
	ostrm<<"<<"<<me<<"="<< (const MGGel*)this;
	ostrm<<", equal_to_binder="<<m_equal_to_binder;
	is_bcell()?MGBCell::toString(ostrm):MGPCell::toString(ostrm);
	ostrm << std::endl;
	MGCell::toString(ostrm);

	const UniquePVertex& pvS = m_vertex[0];
	const UniquePVertex& pvE = m_vertex[1];
	if (!pvS || !pvE) {
		ostrm << endl << ",vertex[";
		if (pvS)
			ostrm << (const MGGel*)pvS.get();
		else
			ostrm << "Null";
		ostrm << ",";
		if (pvS)
			ostrm << (const MGGel*)pvE.get();
		else
			ostrm << "Null";
		ostrm<<"],";
	}
	if(pvS)
		ostrm<<endl<<",Vertex0="<<*pvS;
	if(pvE)
		ostrm<<endl<<",Vertex1="<<*(pvE);
	if(pvS || pvE){
		MGBVertex* bvS = pvS ? pvS->binder_vertex():0 ;
		MGBVertex* bvE = pvE ? pvE->binder_vertex():0 ;
		if(bvS) ostrm<<std::endl<<*bvS;
		if(bvE) ostrm<<std::endl<<*bvE;
	}
	ostrm<<endl;
	ostrm<<"="<<me<<" >> "<<endl;
	return ostrm;
}

//Obtain the parameter of the binder edge's curve that represent
//the same point as sp. sp is a parameter value of this parameter edge.
//Let S() is the star(surface) of this edge, and fp() is the curve of this cell
//which is a boundary of S(). And fb() is the binder curve of this edge.
//Then S(fp(sp))=fb(param_bcell(sp)).
//This is a parameter edge and have the binder, and the parameter sp is a parameter
//of this cell's curve. If this does not have a binder, return -1.
double MGEdge::param_bcell(double sp, const double* guess)const{
	const MGEdge* bdr=binder_edge();
	const MGSurface* srf=star_surface();
	if(!bdr || !srf) return -1.;

	const MGCurve& bcrv=*(bdr->base_curve());
	MGPosition P=srf->eval(eval(sp));
	double tb;
	if(guess){
		double tguess=*guess;
		tb=tguess;
		if(!bcrv.perp_guess(1.,0.,P,tguess,tb))
			bcrv.on(P,tb);
	}else{
		bcrv.on(P,tb);
	}
	return tb;
}

//This must be a parameter edge.
//Obtain the parameter of this parameter edge's curve that represent the same
//point as the binder edge's paramter tb.
//Let S() is the star(surface) of this edge, and fp() is the curve of this cell
//which is a boundary of S(). And fb() is the binder curve.
//Then S(fp(param_pcell(tb)))=fb(tb).
//This edge must have the binder edge, and the parameter tb is the parameter
//of the binder edge's curve. If this does not have a binder, return -1.
double MGEdge::param_pcell(double tb, const double* guess)const{
	const MGEdge* bdr=binder_edge();
	const MGCurve& pcrv=trimmed_curve();
	const MGSurface* srf=star_surface();
	if(!bdr || !srf)
		return -1.;

	return srf->param_of_pcurve(tb,bdr->trimmed_curve(),pcrv,guess);
}

//Obtain start or end parameter value of the edge.
double MGEdge::param_e()const	//End parameter value.
{
	if(active_end())
		return m_vertex[1]->t();

	double t=mgInfiniteVal;
	const MGCurve* crv=base_curve();
	if(crv)
		t=crv->param_e();
	return t;
}
double MGEdge::param_s()const	//Start parameter value.
{
	if(active_start())
		return m_vertex[0]->t();

	double t=-mgInfiniteVal;
	const MGCurve* crv=base_curve();
	if(crv)
		t=crv->param_s();
	return t;
}

//Compute the parameter value of the closest point from the straight to
//this object.
//sl is the eye projection line whose direction is from yon to hither, and if
//sl had multiple intersection points, The closest point to the eye will be selected.
MGPosition MGEdge::pick_closest(const MGStraight& sl)const{
	MGCurve* crv=curve_limitted();
	MGPosition prm=crv->pick_closest(sl);
	delete crv;
	return prm;
}

//Obtain partner edges.
//Partners represent same world's(same cell's parameter) coordinates.
//Parameter edges' partners are parameter edges.
//Binder edges' partners are binder edges.
//The partners do not include this edge except when star cell is
//connected to the star cell itself(closed only by the star cell).
std::vector<const MGEdge*> MGEdge::partner_edges() const{
	assert(is_pcell());
	std::vector<const MGPCell*> cells=partners();
	std::vector<const MGEdge*> pedges;
	int ncells=(int)cells.size();
	for(int i=0; i<ncells; i++){
		const MGEdge* edgei= dynamic_cast<const MGEdge*>(cells[i]);
		pedges.push_back(edgei);
	}
	return pedges;
}

//Get previous edge in the loop sequence.
//The pre_edge is the first neighbour edge.
const MGEdge* MGEdge::pre_edge(bool at_start) const{
	const MGEdge* pre=0;
	int vid = at_start ? 0:1;
	if(m_vertex[vid]){
		std::vector<const MGPCell*> ptrs=m_vertex[vid]->partners();
		size_t m=ptrs.size();
		if(m){
			pre= dynamic_cast<const MGEdge*>(ptrs[m-1]->star());
		}
	}
	return pre;
}
MGEdge* MGEdge::pre_edge(bool at_start){
	const MGEdge* cthis=this;
	return const_cast<MGEdge*>(cthis->pre_edge(at_start));
}

//Get parameter range of the edge.
MGInterval MGEdge::range()const{
	MGEReal t0(MGINFINITE_TYPE::MGINFINITE_MINUS), t1(MGINFINITE_TYPE::MGINFINITE_PLUS);
	if(active_start())
		t0=m_vertex[0]->t();
	if(active_end())
		t1=m_vertex[1]->t();
	const MGCurve* crv=base_curve();
	if(MGREqual_base(t0,t1, MGEReal(crv->param_span()))){
		double error=crv->param_error();error*=.5;
		double tmid=(t0.value()+t1.value())*.5;
		return MGInterval(tmid+error, tmid-error);//This is empty interval.
	}
	return MGInterval(t0,t1);
}

///Release the extent of this binder cell.
MGGeometry* MGEdge::release_extentBC(){
	invalidateBox();
	return m_extent.release();
}

///Set extent of this binder cell.
void MGEdge::set_extentBC(std::unique_ptr<MGGeometry>&& extent){
	assert(dynamic_cast<MGCurve*>(extent.get()));
	set_extent(std::move(extent));
}
void MGEdge::set_extent_as_nullBC(){
	set_extent_as_null();
}

//Set binder cell edge to this parameter cell.
//This curve's coordinates are parameters of a face. And input wcrv's
//coordinates are world coordinate of the face.
//Parameter range of the curve is from start to end of the wcrv when no range
//is specified.
MGEdge* MGEdge::set_binder_edge(const MGCurve& wcrv)const{
	MGEdge* bndr=new MGEdge(wcrv);
	resetBinder(bndr);
	return bndr;
}
MGEdge* MGEdge::set_binder_edge(const MGCurve& wcrv, const MGInterval& range)
const{
	MGEdge* bndr=new MGEdge(wcrv,range);
	resetBinder(bndr);
	return bndr;
}
MGEdge* MGEdge::set_binder_edge(MGCurve* wcrv)const{
	MGEdge* bndr=new MGEdge(wcrv);
	resetBinder(bndr);
	return bndr;
}
MGEdge* MGEdge::set_binder_edge(MGCurve* wcrv, const MGInterval& range)const{
	MGEdge* bndr=new MGEdge(wcrv,range);
	resetBinder(bndr);
	return bndr;
}

//Set(make) start or end vertex at the parameter t.
//When id=0, start. id1, end.
void MGEdge::setStartOrEnd(int id, double t){
	assert(0<=id && id<=1);
	if(m_vertex[id]){
		m_vertex[id]->resetBinder();
		m_vertex[id]->set_t(t);
	} else{
		m_vertex[id].reset(new MGPVertex(t, this));
	}
	id ? base_curve()->trim_end(t):	base_curve()->trim_start(t);
}

///Set only parameter range of this edge.
///Does not change the edge connection like set_start or set_end.
void MGEdge::set_only_param_range(double ts, double te){
	if(m_vertex[0]){
		m_vertex[0]->set_t(ts);
	}else{
		m_vertex[0].reset(new MGPVertex(ts,this));
	}
	if(m_vertex[1]){
		m_vertex[1]->set_t(te);
	}else{
		m_vertex[1].reset(new MGPVertex(te,this));
	}
}

///Obtain star cells.
const MGCell* MGEdge::star() const{
	assert(is_pcell());
	const MGCell* astar=nullptr;
	const MGComplex* boundry = parent_complex();
	if(boundry) astar= boundry->star();
	return astar;
}
MGCell* MGEdge::star(){
	const MGCell* cthis=const_cast<const MGEdge*>(this)->star();
	return const_cast<MGCell*>(cthis);
}

//Obtain star surface.
//Star cell of this must be a face. If not, return null.
//If does not have star surface, returns null.
const MGSurface* MGEdge::star_surface()const{
	const MGCell* str=star();
	if(str)
		return dynamic_cast<const MGSurface*>(str->extent().get());
	else return 0;
}

//Trim the edge at parameter t.
//When start=true, trim start, and the result is from t to end.
//When start=false, trim end, and the result is from start to t.
void MGEdge::trim(double t, bool start){
	//Set pcell parameter.
	if(start){
		assert((param_e()-t)>parameter_error());//t must not end point.
		free_vertex(0);
		set_start(t);
	}else{
		assert((t-param_s())>parameter_error());//t must not start point.
		free_vertex(1);
		set_end(t);
	}

	//Set binder cell parameter.
	MGCurve* wcrv=world_curve();
	if(wcrv){
		//When binder curve exist, we change the parameter range.
		const MGSurface* srf=star_surface();
		if(srf){
			MGPosition uv(eval(t));
			MGPosition P(srf->eval(uv));
			double t_w; wcrv->on(P, t_w);
			MGEdge* ew=binder_edge();
			if(equal_direction_to_binder()){
				if(start)
					ew->set_start(t_w);
				else
					ew->set_end(t_w);
			}else{
				if(start)
					ew->set_end(t_w);
				else
					ew->set_start(t_w);
			}
		}else{
			//In this case, we cannot modidfy world_curve of this edge.
			//All we can do is to free(delete) the world_curve.
			resetBinder();
		}
	}
	m_box.set_null();
}

//Get trimmed curve representation of the edge.
MGTrimmedCurve MGEdge::trimmed_curve() const{
	return MGTrimmedCurve(*(base_curve()),range());
}

//Return world curve pointer of this edge. That is, curve pointer
//of this edge's binder edge. May be null when no binder, or the binder
//does not have an extent.
MGCurve* MGEdge::world_curve(){
	MGCurve* crv=0;
	MGEdge* bndr=binder_edge();
	if(bndr)
		crv=bndr->base_curve();
	return crv;
}
const MGCurve* MGEdge::world_curve() const{
	const MGCurve* crv=0;
	MGEdge* bndr=binder_edge();
	if(bndr)
		crv=bndr->base_curve();
	return crv;
}


///Get all the MGPCell* of the all the boundaries of this.
std::vector<const MGPCell*> MGEdge::getBoundaryPcells()const{
	std::vector<const MGPCell*> pcells;
	for(int i = 0; i<2; i++){
		if(m_vertex[i])
			pcells.push_back(m_vertex[i].get());
	}
	return pcells;
}

///Get the name of the class.
std::string MGEdge::whoami()const {
	return is_pcell() ? std::string("PEdge"): std::string("BEdge");
}
