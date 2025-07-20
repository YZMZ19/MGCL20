/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGEdge_HH_
#define _MGEdge_HH_

#include <map>
#include "mg/Box.h"
#include "mg/TrimmedCurve.h"
#include "topo/Complex.h"
#include "topo/BCell.h"
#include "topo/PVertex.h"
#include "topo/Cell.h"

class MGStraight;
class MGSurface;
class MGPVertex;
class MGBVertex;
class MGLoop;
class MGFace;

//
//Define MGEdge Class.

/** @addtogroup TOPO
 *  @{
 */

///MGEdge is an instance of MGCell, represents a boundary element of 2D manifold.

///MGEdge constitues an MGLoop that is a boundary of MGFace.
///MGEdge can be a parameter cell or a binder cell. The coordinates of a parameter cell
///MGEdge is (u,v) surface parameter, and the ones of binder cell MGEdge is
///(x,y,z) of the world.
class MG_DLL_DECLR MGEdge: public MGCell, public MGPCell, public MGBCell{

private:

	std::unique_ptr<MGPVertex> m_vertex[2];	///<[0]:for start, [1] for end point.
		///<When m_vertex[.]=0, the side has no boundary
		///<and the boundary is the same as the curve's(geometry of MGCell).

	mutable int m_equal_to_binder;
	///<Flag if this curve's direction is equal to the binder or not.
	///<=0:unknown, =1: equal, =-1:opposite direction.

public:

///Edgeのスケーリングを行い，Edgeを作成する。
///Scaling of the Edge by a double.
MG_DLL_DECLR friend MGEdge operator* (double s, const MGEdge& e);

///////// Constructor /////////

///Special constructors.
MGEdge();
virtual ~MGEdge();

MGEdge(const MGEdge& bv);
MGEdge(MGEdge&& rhs)noexcept;//Move constructor.

MGEdge(const MGEdge& e, bool boundaryIsNecessary);
MGEdge& operator=(const MGEdge& rhs);
MGEdge& operator=(MGEdge&& rhs) = default;//Move assignment.

///Fundamental constructor.
///Construct an edge from geometry of manifold dimension 1.
///The constructor takes the ownership of geo and MGPVertex* in boundaries.
MGEdge(
	MGCurve* geo,
	MGPVertex* boundaries[2]
);

///Make an edge of a boundary that has active start and
///end vertex if the curve is not infinite straight line.
///The second form that input MGCurve* takes the ownership of the crv
///into the MGEdge, must not delete the object and the object must be
///newed one.
explicit MGEdge(const MGCurve& crv);
explicit MGEdge(MGCurve* crv);

///Make an edge of a boundary(MGCell that has active start and end vertex).

///range is the parameter range of crv.
///The second form that input MGCurve* takes the ownership of the crv
///into the MGEdge, must not delete the object and the object must be
///newed one.
MGEdge(const MGCurve& crv, const MGInterval& range);
MGEdge(MGCurve* crv, const MGInterval& range);

///Make an edge with a binder of a boundary
///(MGCell that has active start and end vertex).
MGEdge(
	const MGSurface&surf,///<Parent surface of which this edge makes a boundary
	const MGCurve& pcrv, ///<Parameter curve of the surface surf.
	const MGInterval& prange,///<param range of pcrv.
	const MGCurve& wcrv	///<World coordinate curve of the surface surf.
						///<wcrv will be trimmed by prange of pcrv.
);

///Make a clone of the cell.
///clone() does not copy the binder cell of this.
MGEdge* clone() const override;

///Make a clone of the cell that does not have boundaries.
///Does not copy the binder cell of this.
MGEdge* cloneWithoutBoundary() const override;

///////// operator overload/////////

///Assignment.
///When the leaf object of this and cell2 are not equal, this assignment
///does nothing.
///does not change binder and partner relation,
///does not change parent complex.
MGEdge& operator=(const MGGel& gel2);

/// Edge に平行移動を行ないオブジェクトを生成する。
///Translation of the Edge
MGEdge operator+ (const MGVector& v) const;

/// Edgeに逆方向の平行移動を行ないオブジェクトを生成する。
///Translation of the Edge
MGEdge operator- (const MGVector& v) const;

///Edgeのスケーリングを行い，Edgeを作成する。
///Scaling of the Edge by a double.
MGEdge operator* (double s) const;

/// 与えられた変換でEdgeの変換を行い，Edgeを作成する。
///Transformation of the Edge by a matrix.
MGEdge operator* (const MGMatrix& mat) const;

/// 与えられた変換によってトランスフォームをおこないEdgeを生成する。
///Transformation of the Edge by a MGTransf.
MGEdge operator* (const MGTransf& tr) const;

///Complexのスケーリングを行い，Complexを作成する。
///Scaling of the Complex by a double.
MGEdge operator/ (double s) const{return (*this)*(1./s);};

///Comparison of two objects.
bool operator==(const MGEdge& e2)const;
std::partial_ordering operator<=>(const MGEdge& gel2)const;

//gel2 must be the same class as this.
bool equal_test(const MGGel& gel2)const override;

//gel2 must be the same class as this.
std::partial_ordering ordering_test(const MGGel& gel2)const override;

///Object transformation.
MGEdge& operator+=(const MGVector& v)override;
MGEdge& operator-=(const MGVector& v)override{return operator+=(-v);};
MGEdge& operator*=(double scale)override;
MGEdge& operator*=(const MGMatrix& mat)override;
MGEdge& operator*=(const MGTransf& tr)override;

///Transform the binder.
void binder_tr(const MGVector& v) override{ *this+=v; };//Translation.
void binder_tr(double s) override{ *this*=s; };//Scaling.
void binder_tr(const MGMatrix& mat) override{ *this*=mat; };//Matrix transform.
void binder_tr(const MGTransf& tr) override{ *this*=tr; };//Trans transform.

///Write out the edge data to ostream.
std::ostream& toString(std::ostream& ostrm) const;

/////////Member Function/////////

///Test if active at start or end.
bool active_end() const{ return m_vertex[1]!=0;}
bool active_start() const{ return m_vertex[0]!=0;}

///Get after edge in the loop sequence.
///The aft_edge is the first neighbour edge.
const MGEdge* aft_edge(bool at_end=true, int* vertexID=0)const;
MGEdge* aft_edge(bool at_end=true, int* vertexID=0);

///Obtain binder edge pointer.
///Null when this does not have binder.
MGEdge* binder_edge() const;

///Obtain binder of this if this is a binder(null when not).
const MGBCell* binderCell()const override{return binder_edge();};
MGBCell* binderCell()override{return binder_edge();};

///Obtain the center parameter value of this cell.
MGPosition center_param() const;

///Compute the continuities between this edge(edge1) and the edge2.

///This edge and edge2 must be parameter edges of each face.
///In distance, tangent, and normal, the following output will be set:
///distance[0-6] as:
///	[0] edge1's curve parameter that has the maximum distance with edge2.
///	[1] edge2's curve parameter that has the maximum distance with edge1.
///  [2] the evaluated maximum distance between edge1 and edge2 at distance[0] and [1]
///	[3] edge1's curve parameter that has the minimum distance with edge2.
///	[4] edge2's curve parameter that has the minimum distance with edge1.
///  [5] the evaluated minimum distance between edge1 and edge2 at distance[3] and [4]
///	[6] mean distance between edge1 and edge2.
///tangent[0-3] as:
///	[0] edge1's curve parameter that has the maximum tangent difference with edge2.
///	[1] edge2's curve parameter that has the maximum tangent difference with edge1.
///  [2] the evaluated maximum tangent difference between edge1 and edge2 at tangent[0] and [1].
///	[3] mean tangent difference between edge1 and edge2.
///normal[0-3] as:
///	[0] edge1's curve parameter that has the maximum normal difference with edge2.
///	[1] edge2's curve parameter that has the maximum normal difference with edge1.
///  [2] the evaluated maximum normal difference between edge1 and edge2 at normal[0] and [1].
///	[3] mean normal difference between edge1 and edge2.
void compute_continuity(
	const MGEdge& edge2,
	double diatance[7],
	double tangent[4],
	double normal[4]
)const;

///Connect this edge  to cell2(is an MGEdge). Both edges are parameter edges of faces.

///This cell is a boundary of an MGFace A,
///and cell2 is also is a boundary of another MGFace B.
///connet() binds two faces A and B through this edge and the cell2 edge.
void bind(MGEdge& cell2);

///Connect the start(id1=0) or end(id1=1) of this to the start(id2=0) or
/// the end(id2=1) of e2.

///If both edges of this and e2 are members of a complex, they must be the same.
///e2 must be a newed object, and the owneship is transfered to the system.
void connect_at_id(int id1, MGEdge* e2, int id2);

///Connect the start of this to the pvert's edge at pvert.
void connect_at_start(MGPVertex& pvert){connect(0,pvert);};

///Connect the end of this to the pvert's edge at pvert.
void connect_at_end(MGPVertex& pvert){connect(1,pvert);};

///Return curve pointer of this edge.

///Null when this does not have geometry.
///The expression is of parameter space of face.
MGCurve* base_curve();
const MGCurve* base_curve() const;

///Compute the box from the scratch.
void compute_box(MGBox& bx) const;

///Return curve pointer cut by start and end parameter range.

///Output is newed curve object, must be deleted.
///Null when this does not have geometry.
///The expression is of parameter space of face if this is parameter edge of a face.
///curve_limitted() does not return MGTrimmedCurve, returns real curve.
MGCurve* curve_limitted() const;

///Draw curve in this coordinates.
void drawWire(
	mgVBO& vbo,///<The target graphic object.
	int line_density=1	///<line density to draw a surface in wire mode.
)const;

//Draw 3D curve in the star face world coordinates.
///This is a boundary edge(MGPCell) of a face, and the curves are extracted
///from the binder edge and drawn.
void drawWire_in_star(
	mgVBO& vbo,
	int line_density	//line density to draw a surface in wire mode.
)const;

///Draw 3D point(vertex) in world coordinates.
///The object is converted to point(s) and is drawn.
///This is valid only for topology objects or MGPoint.
void drawVertex(
	mgVBO& vbo///<The target graphic object.
)const;

//Obtain this edge's parent loop's edge iterator as a member of the parent loop.
//This must have the parent loop, i.e. loop() must not null.
MGComplex::const_iterator edge_iterator()const;
MGComplex::iterator edge_iterator();

//Obtain this edge number as a member of the parent loop.
//This must have the parent loop, i.e. loop() must not null.
int edge_num()const;

///Obtain the end point of the edge.
MGPosition end_point()const{return eval(param_e());};

///Test if this has non null vertex at start(id=0) or end(id=1),
///and if does not have, make one.
void ensureHasVertex(int id);

///Test if this has non null vertex at start and end,
///and if does not have, make ones.
void ensureHasVerticesAtBothEnds();

///Evaluate the nderiv's derivative at parameter t.
///Evaluate of the curve's data.
MGVector eval(double t, int nderiv=0)const;

///Evaluation of the star curves of the edge at the point t.
///When nderi=0, get a position of the surface at the boundary point t.
///The star curve is SurfCurve(face's surface, edge's curve).
///(The star curve has the same world coordinate with the binder curve's, but
///their direction may be opposite. The star curve has always the same direction
///as the loop.)
MGVector eval_star(
	double t,		///<Parameter value of this parameter edge's curve.
	int nderi=0	///<Order of derivative.
)const;

///Test if SurfCurve of the edge has equal direction to binder edge's direction.
///Returned is true if eaual, false if not.
bool equal_direction_to_binder()const;

///Get extent geometry, may be null if this does not have extent.
const MGCurve* extentBC() const override;
MGCurve* extentBC()  override;

///Get the star face pointer.
const MGFace* face() const;
MGFace* face();

///Obtain star cells.
const MGCell* star() const override;
MGCell* star() override;

///Get the 1st partner edge of this edge.
///This is a parameter cell.
const MGEdge* first_partner() const;

///If this had boundary binders, free them. As the result, this
///will have no neighbours.
void free_neighbours() override;

///Return Object's type ID (TID)
long identify_type()const;
long identify_typeBC()const{ return identify_type(); };
long identify_typePC()const{ return identify_type(); };

///Ask if this is binder cell.
bool is_bcell() const{ return number_of_partner_members()>=1; };
bool is_pcell() const{ return number_of_partner_members()==0; };

///Test if this edge's start point(when start=true) and edge2 is connected
///and their directions are the same.
///When start=false, this edge's end point is tested.
bool is_connected_and_same_direction(
	bool start,
	const MGEdge& edge2
)const;

///Comparison of two MGEdges as MGPCells.
///is_less_than() defines the order of partners to store in MGBCell.
bool is_less_than(const MGPCell& pcel2)const override;
bool is_less_than(const MGEdge& pe2)const;

///test if parameter t is the one of the end point of the loop.
bool is_end_point(double t)const;

///test if parameter t is the one of the start point of the loop.
bool is_start_point(double t)const;

///Test if this is a free edge.
///Free edges are ones that do not have partner edges.
bool is_free()const{ return number_of_partners()==0;};

///Intersection of a shell and a curve.
MGisects isect(const MGObject& obj2)const override{
	return MGisects();
};

///Connect this and e2.
///If start==true, start of this edge to end of e2;
///If start==false, end of this edge to start of e2;
///e2 must be a newed object, and the ownership is transfered to the system.
void join(bool start, MGEdge* e2);

///Return parent loop pointer.
const MGLoop* loop() const;
MGLoop* loop();

///Make a binder cell of this parameter cell.
///Returned is the binder reference.
///The binder may have no geometry, and only has binder and parameter cell relationship.
std::shared_ptr<MGBCell>& make_binder() const override;

///Make a binder associated with the world curve rep.
///Returned is the binder edge pointer.
///If the parameter edge had already the binder,
///make_binder_with_curve only returns the pointer.
///*** This edge must be a member of a loop that is a boundary of a face.
MGEdge* make_binder_with_curve()const;

///Make sure that this has an extent expression.
///When this did not have an extent, make the extent from the partner
///member's parameter expression and the star cell.
///This must be a binder cell that has partner members that are
///boundaries. When this is not the case or this had an extent already,
///it does nothing.
void make_extent() const override;

///Obtain manifold dimension.
int manifold_dimension() const{return 1;};

///Obtain the i-th member partner edge. This must be a binder edge.
const MGEdge* partner_member_edge(int i)const;

///Compute the mid point of this edge.
///Mid point is the point of the paramete mid=(param_s()+param_e())*.5
MGPosition mid_point()const;

///Negate the direction of the cell.
void negate();

///Obtain all the neighbours.
///Neighbours are 
///The neighbours do not contain this cell except when this cell is
///connected to this cell itself(closed cell).
std::vector<const MGCell*> neighbours() const;

///Test if the edge is a part of a surface perimeter.
bool on_surface_perimeter() const{return surface_perimeter()>=0;};
bool on_surface_perimeter(const MGFace& f) const{return surface_perimeter(f)>=0;};
bool on_surface_perimeter(const MGSurface& sf) const{return surface_perimeter(sf)>=0;};

///Obtain the parameter of the binder edge's curve that represent
///the same point as sp.

///sp is a parameter value of this parameter edge.
///Let S() is the star(surface) of this edge, and fp() is the curve of this cell
///which is a boundary of S(). And fb() is the binder curve of this edge.
///Then S(fp(sp))=fb(param_bcell(sp)).
///This is a parameter edge and have the binder, and the parameter sp is a parameter
///of this cell's curve. If this does not have a binder, return -1.
double param_bcell(double sp, const double* guess=0)const;

///Obtain the parameter of this parameter edge's curve that represent the same
///point as the binder edge's paramter tb.

///This must be a parameter edge.
///Let S() is the star(surface) of this edge, and fp() is the curve of this cell
///which is a boundary of S(). And fb() is the binder curve.
///Then S(fp(param_pcell(tb)))=fb(tb).
///This edge must have the binder edge, and the parameter tb is the parameter
///of the binder edge's curve. If this does not have a binder, return -1.
double param_pcell(double tb, const double* guess=0)const;

///Obtain end parameter value of the edge.
double param_e()const;

///Obtain start parameter value of the edge.
double param_s()const;

///Obtain parameter span of this edge.
double param_span()const{return param_e()-param_s();};

///Obtain partner edges.

///Partners represent same world's(same cell's parameter) coordinates.
///Parameter edges' partners are parameter edges.
///Binder edges' partners are binder edges.
///The partners do not include this edge except when star cell is
///connected to the star cell itself(closed only by the star cell).
std::vector<const MGEdge*> partner_edges() const;

///Compute the parameter value of the closest point from the straight to
///this object.

///sl is the eye projection line whose direction is from yon to hither, and if
///sl had multiple intersection points, The closest point to the eye will be selected.
MGPosition pick_closest(const MGStraight& sl)const;

///Approximate the parameter edge by a polyline and replace this edge
///expression by the polyline.

///Polyline approximation is so done that
///the correspoinding binder edge can be appximated by the polyline connecting
///each binder edge's point that corresponds to the each this edge's point.
///(1) This must be a parameter cell edge.
///(2) This edge must be a member of a loop which is a boundary of a face.
///(3) If this edge did not have a binder edge, polygonize generates the binder edge.
///(The tolerance used to generate the binder is MGTolerance::line_zero(),
/// not input error.)
///Input error is tolerance allowed between the polygon and the original curve.
void polygonize(double error);

///Get previous edge in the loop sequence.

///The pre_edge is the first neighbour edge.
const MGEdge* pre_edge(bool at_start=true) const;
MGEdge* pre_edge(bool at_start=true);

///Get parameter range of the edge.
MGInterval range()const;

///Release the extent of this binder cell.
MGGeometry* release_extentBC() override;

///Set binder cell edge to this parameter cell.

///Set extent of this binder cell.
void set_extentBC(std::unique_ptr<MGGeometry>&& extent)override;
void set_extent_as_nullBC() override;

///This curve's coordinates are of parameter space of a face. And input crv's
///coordinates are world coordinate of the face.
///range is the parameter range of wcrv.
///Parameter range of the wcrv is from start to end of the wcrv when no range
///is specified.
///Function return value is the binder's pointer generated.
MGEdge* set_binder_edge(const MGCurve& wcrv)const;
MGEdge* set_binder_edge(const MGCurve& wcrv, const MGInterval& range)const;

///These forms give the ownership of wcrv to the edge.
///That is, wcrv must be newed one and users must not delete it.
///Others are same as above "set_binder_edge(const MGCurve& wcrv)" form.
///Function return value is the binder's pointer generated.
MGEdge* set_binder_edge(MGCurve* wcrv)const;
MGEdge* set_binder_edge(MGCurve* wcrv, const MGInterval& range)const;

///Set start point(boundary) data.
///If this is connected to other edges at end, the connectin will be freed.
void set_end(double t){	setStartOrEnd(1, t);};///Parameter value of the start point.

///Set start point(boundary) data.
///If this is connected to other edges at start, the connectin will be freed.
void set_start(double t){ setStartOrEnd(0, t); };///Parameter value of the start point.

///Set only parameter range of this edge.
///Does not change the edge connection like set_start or set_end.
void set_only_param_range(double ts, double te);

///Obtain star surface.
///Star cell of this must be a face. If not, return null.
///If does not have star surface, returns null.
const MGSurface* star_surface()const;

///Obtain the end point of the edge.
MGPosition start_point()const{return eval(param_s());};

///Get the perimeter number where this edge is on.
///If this is not on any perimeter, -1 will be returned.
int surface_perimeter() const;
int surface_perimeter(const MGSurface& sf) const;
int surface_perimeter(const MGFace& face) const;

///Trim the loop. Result is from start to t1.
void trim_end(double t){trim(t,false);};

///Trim the loop. Result is from t1 to end.
void trim_start(double t){trim(t,true);};

///Get trimmed curve representation of the edge.
MGTrimmedCurve trimmed_curve() const;

///Get the vertex at the start or end.
const std::unique_ptr<MGPVertex>& vertex(int id)const{return m_vertex[id];};
const std::unique_ptr<MGPVertex>& vertex_start()const{return m_vertex[0];};
const std::unique_ptr<MGPVertex>& vertex_end()const{return m_vertex[1];};
std::unique_ptr<MGPVertex>& vertex(int id){return m_vertex[id];};
std::unique_ptr<MGPVertex>& vertex_start(){return m_vertex[0];};
std::unique_ptr<MGPVertex>& vertex_end(){return m_vertex[1];};

///Return world curve pointer of this edge. That is, curve pointer
///of this edge's binder edge.
///May be null when no binder, or the binder does not have an extent.
MGCurve* world_curve();
const MGCurve* world_curve() const;

///Get the name of the class.
std::string whoami()const override;

protected:

///Read Object's member data.
void ReadMembers(MGIfstream& buf);

///Write Object's Member Data
void WriteMembers(MGOfstream& buf) const;

///Clone this BCell, building the new partner membership of partners.
SharedBCell cloneWithPartnerMembers(std::vector<const MGPCell*>& newPartners)const override;

//free start(i=0) or end(i=1) boundary's bindness.
void free_vertex(int i);

///Free specified boundary(bound) from a member of parent cell's boundaries.
///Return MGComplex* bound if freed normally.
///If bound was not a member of the boundaries, return 0.
///Only free, does not destruct the boundary.
MGComplex* free_boundary(const MGComplex* bound)override{return 0;};

///Get all the MGPCell* of the all the boundaries of this.
std::vector<const MGPCell*> getBoundaryPcells()const override;

private:

///Connect the start(id=0) or end(id=1) of this to the pvert's edge at pvert.
void connect(int id, MGPVertex& pvert);

///Compute continuity, given the evaluation interval and the division number.
///This is a parameter edge that has the star face.
void compute_continuity2(
	const MGInterval& span,///<this edge's binder edge's parameter span
	int npoint,				///<division number of this edge's interval sspan
	const MGEdge& edge2,	///<the second parameter edge that has the star face.
	double distance[7],		///<evaluated data will be set in distance, tangent, and normal.
	double tangent[4],		///<See compute_continuity
	double normal[4]
)const;

///Copy all boundaries of cellin into this(but does not copy own binder cell relation),
//and register cellin's boundary MGPCell* association into cmap.
//1st(key) is original MGPCell*(MGPVertex) of boundaries(m_vertex[]).
//2nd is copyied new.
void copy_all_boundaries(
	const MGCell& cellin,
	std::map<const MGPCell*, MGPCell*>* cmap=0
) override;

///Make this cell's binder cell's extent expression.
///Returned is a MGGeometry pointer generated by new.
///When this cell does not have star cell, null pointer will be returned.
///make_binder_extent() only makes the expression, and does nothing to
///the topology structure.
std::unique_ptr<MGGeometry> make_binder_extent() const override;

///Negate the boundary.
void negate_boundary()override;

//Set(make) start or end vertex at the parameter t.
//When id=0, start. id1, end.
///If this is connected to other edges at id, the connectin will be freed.
void setStartOrEnd(int id, double t);

///Set both end parameter vertices.
void setBothEnds(const MGInterval& range);

///Trim the edge at parameter t.
///When start=true, trim start, and the result is from t to end.
///When start=false, trim end, and the result is from start to t.
void trim(double t, bool start);

friend class MGFace;
friend class MGLoop;
friend class MGIgesIfstream;

};

/** @} */ // end of TOPO group
#endif
