/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGLoop_HH_
#define _MGLoop_HH_

#include "mg/Curve.h"
#include "topo/LSPoint_vector.h"
#include "topo/Complex.h"
#include "topo/TrimLoop.h"

class MGInterval;
class MGPosition;
class MGCurve;
class MGStraight;
class MGSurface;
class MGLCisect_vector;
class MGLLisect_vector;
class MGPVertex;
class MGLEPoint;
class MGLPoint;
class MGLoop;
class MGEdge;
class MGFSurface;
class MGFace;
/** @file */
/** @addtogroup TOPO
 *  @{
 */

//
//Define MGLoop Class.

///MGLoop is a boundary of a face, a boundary of 2D manifold cell.

///MGLoop accepts parameter space curve and world space curve of
///of a boundary curve, and constructs a boundary of a face from the
///two types of curves.
///Input curves direction indicate which part of the face will be target
///part after trimed by the boundary. In 2D space (u,v) of the parameter
///space, LEFT side of the parameter curve along the curve's direction
///is the target part of face.
class MG_DLL_DECLR MGLoop:public MGComplex{

public:
	enum LoopKind{
		UNDEFINED = -1,
		INACTIVE = 0,
		PERIMITER_LOOP = 1,
		OUTER_LOOP = 2,	///<must be cloded.
		INNER_LOOP = 3,	///<must be cloded.
		NETWORK = 4
	};

private:
	mutable LoopKind m_kind;	///<Kind of loop.

	mutable double m_area;		///<When closed, loop area will be stored.
	///<m_area>0. means inside of this closed loop is active(outer boundary).
	///<m_area<0. means outside of this closed loop is active(inner boundary).

	mutable int m_perim_num_s; ///<Start point perimeter number.
	mutable int m_perim_num_e; ///<End point perimeter number.
	///<valid only when m_kind==PERIMITER_LOOP.
	///<Hold the perimeter number start or end point of the loop lies on.

public:

///Get edge pointer from its iterator in MGComplex(MGLoop).
MG_DLL_DECLR friend MGEdge* edge_from_iterator(iterator i);
MG_DLL_DECLR friend const MGEdge* edge_from_iterator(const_iterator i);

///Evaluation of the loop at the point t.

///When nderi=0, get a parameter (u,v) of the surface at the boundary point.
MG_DLL_DECLR friend MGVector eval(const MGLEPoint& t, int nderi);

///Test if (u,v) is inside the outer boundary(of std::vector<MGLoop*>& boundaries).

///Inside the outer boundary means that inside outer_boudary_param() or not.
///This must not be used for faces that do not have perimeter or outer boundary
///loop.
///Function's return value is:
///  0:outside the outer boundary(not on a loop)
///  1:unknown
///  2:inside the outer boundary(not on a loop)
/// otherwise:on the outer boundary loop 
MG_DLL_DECLR friend size_t inside_outer_loop(
	const MGPosition& uv,
	const std::vector<const MGLoop*>& loop,
	const MGSurface* surf
);

/////////Constructor/////////

///Void constructor.
MGLoop();
~MGLoop();

///Copy constructor.
MGLoop(const MGLoop& loop2)=default;
MGLoop(MGLoop&& loop2) = default;

///Construct a loop of one edge.
explicit MGLoop(MGEdge* edge);

///Construct a Loop of one edge of one curve cell.
///param_curve is parameter space representation of the face
///of which this loop will be a boundary, will make parameter cell.
///world_curve is world coordinate representation of the face
///of which this loop will be a boundary, will make binder cell.
///range1 is parameter range of the curve param_curve,
///range2 is parameter range of the curve world_curve.
///When range1,2 are not specified, the start and the end of the curve are
///treated as their ranges.
///***param_curve and world_curve must have the same direction.
MGLoop(const MGCurve& param_curve, const MGCurve& world_curve);
MGLoop(const MGCurve& param_curve, const MGInterval& range1,
	   const MGCurve& world_curve, const MGInterval& range2);
MGLoop(std::unique_ptr<MGCurve>& param_curve, std::unique_ptr<MGCurve>& world_curve);

///Assignment.
///When the leaf object of this and bnd2 are not equal, this assignment
///does nothing.
MGLoop& operator=(const MGGel& gel2);
MGLoop& operator=(const MGLoop& gel2);
MGLoop& operator=(MGLoop&& gel2);

///Object transformation.
MGLoop& operator+=(const MGVector& v);
MGLoop& operator-=(const MGVector& v);
MGLoop& operator*=(double scale);
MGLoop& operator*=(const MGMatrix& mat);
MGLoop& operator*=(const MGTransf& tr);

/// Complex�ɕ��s�ړ����s�Ȃ��I�u�W�F�N�g�𐶐�����B
///Translation of the Complex
MGLoop operator+ (const MGVector& v) const;

/// Complex�ɋt�����̕��s�ړ����s�Ȃ��I�u�W�F�N�g�𐶐�����B
///Translation of the Complex
MGLoop operator- (const MGVector& v) const{ return operator+(-v); };

///Complex�̃X�P�[�����O���s���CComplex���쐬����B
///Scaling of the Complex by a double.
MGLoop operator* (double s) const;

/// �^����ꂽ�ϊ���Complex�̕ϊ����s���CComplex���쐬����B
///Transformation of the Complex by a matrix.
MGLoop operator* (const MGMatrix& mat) const;

/// �^����ꂽ�ϊ��ɂ���ăg�����X�t�H�[���������Ȃ�Complex�𐶐�����B
///Transformation of the Complex by a MGTransf.
MGLoop operator* (const MGTransf& tr) const;

///Complex�̃X�P�[�����O���s���CComplex���쐬����B
///Scaling of the Complex by a double.
MGLoop operator/ (double s) const{ return (*this)*(1./s); };

///This operator is to sort loops in the order:
///  1. Perimeter boundary.
///  2. Outer boundary.
///  3. Inner boundary.
///  4. Inactive loop.
std::partial_ordering operator<=>(const MGLoop& gel2)const;
bool operator==(const MGLoop& lp2)const; //we do not define.

//gel2 must be the same class as this.
bool equal_test(const MGGel& gel2)const override;

//gel2 must be the same class as this.
std::partial_ordering ordering_test(const MGGel& gel2)const override;

///Test if this is active boundary.
bool active() const;

///Append edge to the end of loop.
///"append" connects the edge's start to the end of the loop.
void append(MGEdge* edge);

///Build one edge of srf from the curve wcrv on srf and common information
///pspan and peri_num, which are a perimeter peri_num's parameter spans(psapn).
///wcrv must not be a MGCompositeCurve.
///One edge is generated and append to this loop.
///Function's return value is MGEdge* appended.
MGEdge* append_edge_from_crvWorld(
	const MGSurface& srf,///<The target surface.
	const MGCurve& wcrv,///<curve of world coordinates on this face that may
		///<coincide to a perimeter of the surface of this face.
	double& tLast,///<wcrv's parameter value to start is input and the end param value
		///<of the last edge generated will be output.
	double terror,///<wcrv's parameter space error.
	const std::vector<double>& pspan,///<common parameter value of a perimeter peri_num.
	int peri_num,///<pspan and peri_num are output of getPerimeterCommon(). Refer to it.
	bool orientation_is_opposite=false///<orientation flag of wcrv to edge. True if opposite.
);

///Compute curvilinear integral of the parameter space of the area
///sorrounded by the loop.
double area()const;

///Test if loop is a perimeter boundary or not.
///If yes, get the perimeter ids.
///pid_s, e are valid only when both_end_on_perimeter() is true,
///contain perimeter numbers of the start or end of the loop.
bool both_end_on_perimeter(int& pid_s, int& pid_e,const MGFSurface* srf=0) const;

///Make a clone.
///Output is a newed object, must be deleted by calling program.
MGLoop* cloneWithParent(MGCell& parent)const;
MGLoop* clone()const;

///Test if this is closed boundary.
bool closed() const;

///Compute closest point from the point P to this loop.
///Returned is the loop's point, and distance is the length of distance
///between P and MGLEPoint.
///All the coordinate values are of parameter space of the face's surface.
MGLEPoint closest(const MGPosition& P, double& distance) const;

///Compute closest point closest(P in the world) of the loop(MGLEPoint)
///from the straight sl to this loop's world coordinate rep.

///The function's return value is:
///true: if a point whose distance is smaller than dist is found,
///false: if not.
///Returned to closest is the loop's point, and the distance between
///P and MGLEPoint's world rep is set into dist if found.
///The coordinates of P are of world coordinate space of the loops
///world coordinate rep. This must be a boundary of a surface(that is,
///face() must not be null.
bool closest_world(
	const MGStraight& sl,
	MGLEPoint& closest,
	double& dist
)const;

///Compute common range of two loops, this and loop2.

///Function's return value is number of ranges obtained.
///In variable ranges, common ranges are output as:
///  Let n be output of the function, then ranges1.size()=ranges2.size()=2*n.
///  ranges1[2*i+0] and ranges1[2*i+1] are parameter values of this loop, and
///  ranges2[2*i+0] and ranges2[2*i+1] are parameter values of loop2.
///  Although ranges1[2*i+0] < ranges1[2*i+1] always holds, 
///  ranges2[2*i+0] < ranges2[2*i+1], or ranges2[2*i+0] > ranges2[2*i+1].
///  Let f1() be this loop, and f2() be loop2, then
///  f1(ranges1[j]) and f2(ranges2[j]) represent the same point in star
///  Face world for 0<=j<n*2 .
///This and loop2 must have each star faces.
int common(
	const MGLoop& loop2,
	std::vector<MGLEPoint>& pranges1,
	std::vector<double>& branges1,
	std::vector<MGLEPoint>& pranges2,
	std::vector<double>& branges2
)const;

///Compute curvilinear integral of the loop.
///Computation is done in the parameter space of the face.
double compute_area()const;

///Copy loop data into this.
///This boundary data is cleared and loop's boundary is copied into this.
void copy_boundary(const MGComplex& loop);

///Obtain vector of curves(TrimmedCurve) of the loop.
///The curves are of parameter space expression.
///Let crvs be the output of curves() and wcrvs of curves_world(),
///then crvs[i] corresponds to wcrvs[i] one by one.
std::vector<UniqueCurve> curves()const;

///Obtain vector of curves(world coordinate expression) of the loop.
///Output curves are MGTrimmedCurve of edge's binders.
///When some of the edges do not have binders, they will be created.
///Let crvs be the output of curves() and wcrvs of curves_world(),
///then crvs[i] corresponds to wcrvs[i] one by one.
std::vector<UniqueCurve> curves_world()const;

///Draw vertices of this loop. The vertices coordinate kind is the same as
///the edges. That is, the vertices are MGBVertex.
void drawVertex(
	mgVBO& vbo
)const;

///Draw 3D curve in the star face world coordinates.
///This is a boundary loop(MGPCell's) of a face, and the curves are extracted
///from the binder edges of this and drawn.
virtual void drawWire_in_star(
	mgVBO& vbo,///<Target graphic object.
	int line_density = 1	///<line density to draw a surface in wire mode.
)const;

///Return i-th edge pointer.
MGEdge* edge(int i);
const MGEdge* edge(int i) const;

///Test if at least one edge is included in this loop.
bool edge_exist() const{return pcell_exist();};

///Get edge number in this loop.
///If e is not a member of this loop, 0 will  be returned.
int edge_num(const MGEdge* e)const;

///Return end point of this loop as MGPosition.
///The point is of parameter space.
MGPosition end_point() const;

///Return end point of this loop as MGLEPoint.
///loop must include at least one edge, or this output is undefined.
MGLEPoint end_LPoint() const;

///Get error of this loop.
///This is obtained from the parent surface. If parent surface did not
/// exist, error=the box of the loop by relative zero.
double error()const;

///Evaluation of the loop at the point t.
MGVector eval(const MGLPoint& t, int nderi=0)const;

///Evaluation of the loop at i-th edge's parameter t.
MGVector eval(int i, double t, int nderi=0)const;

///Return pointer of the face. If the loop is not a boundary of any face,
///null will be returned.
const MGFace* face() const;
MGFace* face();

///Return pointer of the first edge.
const MGEdge* first_edge() const;
MGEdge* first_edge();

///Get the loop id of this loop in the star face baoundary.
///Let face=face(), then face->loop(get_loop_id_in_face())=this;
///When this does not have star face, or this is not a boundary of a face,
///-1 will be returned.
int get_loop_id_in_face()const;

///Test if parameter value (u,v) is inside this loop or not.

///inside means inside face, that is,
///if the loop is inner, outside inner loops and
///if the loop is outer boundary loop, inside the outer boundary loop.
///This can be used for perimeter boundary loops.
///Returned is:
///  0:outside(not on the loop)
///  1:unknown
///  2:inside(not on the loop)
/// otherwise:on the loop(int(parameter edge id +100))will be returned.
int inside(double u, double v) const;
int inside(const MGPosition& uv) const;

///Return Object's type ID (TID)
long identify_type()const;

///Intersection of a shell and a curve.
MGisects isect(const MGObject& obj2)const override{
	return MGisects();
};

///Compute intersections of this loop(parameter rep of a face)
///and param_curve.
MGLCisect_vector isectLoop(const MGCurve& param_curve) const;

///Compute intersection points of 1D sub curves of the original loop.
///Parameter values of intersection points(MGLEPoint's) will be returned.
std::vector<MGLEPoint> isect_1D(						
	double f,			///< Coordinate value
	int coordinate=0	///< Coordinate kind of the data f(from 0).
) const;	

///Compute intersections of this loop(parameter rep of the face) and param_curve.

///isect_with_endpoints() includes endpoints of both if they are close enough
///(within tolerance).
MGLCisect_vector isect_with_endpoints
(const MGCurve& param_curve) const;

///Test if this loop is inactive or not.
bool is_inactive(const MGFSurface* srf=0) const;

///Test if this loop is inner boundary.

///Inner boundary is:
///  (1) closed loop. (2) the direction is clockwise.
bool is_inner_boundary(const MGFSurface* srf=0) const;

///Test if this loop is outer boundary.

///outer boundary is:
///  (1) closed loop. (2) the direction is anti-clockwise.
bool is_outer_boundary(const MGFSurface* srf=0) const;

///Test if this loop is perimeter boundary.

///Perimeter boundary is:
///  both ends are on the surface perimeter.
bool is_perimeter_boundary(const MGFSurface* srf=0) const;

///Test if this loop is network.
bool is_network(const MGFSurface* srf=0) const{return m_kind==NETWORK;};

///Compute intersections of two loops.
MGLLisect_vector isect(const MGLoop& loop2)const;

///Join two loops.

///start indicates which end of this loop loop2 should be connected to.
///start=true: connect start of this and loop2's end.
///start=false: connect end of this and loop2's start.
///The 2nd form takes the ownership of loop2, will delete loop2.
void join(bool start, const MGLoop& loop2);
void join(bool start, MGLoop* loop2);
void join(bool start, std::unique_ptr<MGLoop>& loop2);

///Join adjacent two C1 continuous edges to one edge.
///When the number of edges is equal to or less than 4, join is not executed.
void join_C1_edges(
	int parameter_normalization=2
		///Specify how rebuilt edes's parameterization be done.
		//=0: no parameter normalization.
		//=1: normalize to range=(n, n+1) for the n-th edge.
		//=2: normalize to make the average length of the 1st derivative 
		//    is as equal to 1. as possible. Let the n-th edge's parameter range
		//    is (TSn, TEn), then, TSn=TEn-1.
);

///Return edge pointer of the last edge.
const MGEdge* last_edge() const;
MGEdge* last_edge();

///Make this loop as closed.

///This loop's 1st edge's start point must be the same as the last edge's end point.
///However, this is not tested in make_close.
void make_close();

///Make a vertex at lp and subdivide the edge into two edges.
///Returned is true if subdivision is done and false if no subdivision
///is done since lp was one of existed vertex.
///When function's return value is true, pre is always the same edge as
///lp's edge and the iterator is unchanged.
bool make_vertex(
	const MGLEPoint& lp,///<point to subdivide of this loop.
	MGEdge*& pre,		///<pre and aft-edge of the lp will be output
	MGEdge*& aft,		///<after make_vertex's execution, may be null.
	MGPVertex** pvertex=0 ///MGPVertex of the made vertex which is the boundary of lp.
);

///Get manifold dimension.
int manifold_dimension() const{return 1;};

///Merge with param_curve as a network loop.

///Function's return value is:
///  true: merge is done and param_curve is processed into this loop as a network.
///  false: merge is not performed since no intersection with this loop were found.
///When false is returned, this loop does not have the input param_curve information
///as an edge information.
///This loop must not be empty loop, ant the kind is always changed to NETWORK.
bool merge_network(
	const MGCurve& param_curve,///<Curve to merge.
	double error	///<Error to judge two param_curve' points are equal.
);

///Merge loop2 or param_curve to the existing loop data, and build a loop.

///Returned is if merge was done(true) or not(false).
///When no intersection was found, merge is not executed.
///param_curve is parameter space representation of the face
///of which this loop will be a boundary, will make a parameter cell.
///world_curve is world coordinate representation of the face
///of which this loop will be a boundary, will make binder cell.
///range1 is parameter range of the curve param_curve.
///When world_curve is input, it will be trimmed accordin to the range
///of param_curve.
///param_curve and world_curve may be opposite direction.
///When more than one closed loops are detected, first one from
///the old loop start point is employed, and other loops are discarded.
bool merge_trim(const MGCurve& param_curve);
bool merge_trim(const MGCurve& param_curve, const MGInterval& range1);
bool merge_trim(const MGCurve& param_curve, const MGCurve& world_curve);
bool merge_trim(const MGLoop& loop2);

///Compute the mid point of this loop.

///Mid point is the point of the mid point of m-th edge,
///where m=number_edges()/2.
MGPosition mid_point()const;

///Reverse the direction of the boundary.
///(Coordinate transformation is not performed.)
void negate();

///Negate the boundary according to the parent cell negation.

///1. Transform the coordinates of the bondary cell.
///(This transfromation depends on how the parent cell is transformed
///when negate() is invoked. So, the member cells of this boundary
///are transformed by negate_transoform of the parent cell.)
///2. Reverse the direction of the parameter cells(negate each cell).
///3. Reverse the ordering of the parameter cells.
///4. Negate the binders.
void negate_as_boundary(const MGCell* parent=0);

///Get the number of edge included.
int number_of_edges()const{return number_of_pcells();};

///Test if start or end point of the loop is on perimeter of the surface.

///Returned is true when on perimeter, false if not.
///When on perimeter, perimeter number of the surface is returned in pid_x.
bool on_perimeter_end(int& pid_e,const MGFSurface* surf=0)const;
bool on_perimeter_start(int& pid_s,const MGFSurface* surf=0)const;

///Test if all the edges included are on a surface perimeter.
bool on_surface_perimeter(const MGFace& f)const;

///Output function.
std::ostream& toString(std::ostream&) const;

///Prepend edge to the start of the loop.
void prepend(MGEdge* e);

///Remove pendent edge.
///Function's return value is true if an edge was removed.
bool remove_pendent_edge(const MGFSurface& face);

///Remove garbage edge.
///Error is the parameter error allowed of the edge length to judge removal.
void remove_garbage_edge(double error);

///Obtain parent surface pointer.
const MGSurface* surface()const;

///Return start point of this loop as MGLEPoint.
///loop must include at least one edge, or this output is undefined.
MGLEPoint start_LPoint() const;

///Return start point of this loop as MGPosition.
///The point is of parameter space.
MGPosition start_point() const;

///Subdivide this loop so that one parameter range in ranges becomes one edge.

///In ranges, parameter range of this loop is stored as:
///Let n=ranges.size()/2, then span from ranges[2*i] to ranges[2*i+1] is one
///parameter span that is supposed to be one edge for 0<=i<n.
///Returned are new edge pointers that correspond to ranges[2*i] to ranges[2*i+1]
///after subdivided for 0<=i<n.
///****Currently this does not conform to non_manifold model.
///That is, this loop must not have partner edges already.
std::vector<MGEdge*> subdivide(
	const std::vector<MGLEPoint>& ranges
);

///Subdivide this loop so that one parameter range from le1 to le2 becomes one edge.

///le1, and 2 must be of the same edge.
///new edge between le1 and le2 will be output.
MGEdge* subdivide(
	MGLEPoint& le1, MGLEPoint& le2
);	

///Trim the loop. Result loop is from t1 to t2.

///The loop can be closed one. In this case, t1 can be >t2.
///When not closed, t1 must be less than t2.
void trim(const MGLEPoint& t1, const MGLEPoint& t2);

///Trim the loop. Result is from start to t1.
void trim_end(const MGLEPoint& t1);

///Trim the loop. Result is from t1 to end.
void trim_start(const MGLEPoint& t1);

protected:

///Construct from boundary complex(i.e. MGLoop).
///This constructor takes the ownership of MGCell* in boundary.
MGLoop(
	std::list<MGCell*> boundaries	///<Boundary data of the super class MGComplex.
);

///Copy constructor with member edge's original and copied one mapping.
MGLoop(
	const MGLoop& loop,		///<original Loop.
	std::map<const MGPCell*, MGPCell*>& cmap //cellmap to register MGEdge association.
	 //original edge of loop and copied edge association are registered in cmap.
);

///Make a clone.
///This form is to register MGPCell(MGEdge) association of the member cells.
///Association of the original cell(key) and copied new one is set in cmap.
///When parent is specified, clone's parent is set to the parent.
UniqueLoop cloneWithParentAndMap(MGFace& parent, std::map<const MGPCell*, MGPCell*>& cmap) const;

///Get the name of the class.
std::string whoami()const{return "Loop";};

///Read Object's member data.
void ReadMembers(MGIfstream& buf);

///Write Object's Member Data
void WriteMembers(MGOfstream& buf) const;

private:

///Compute loop kind.
void compute_kind(const MGFSurface* srf=0)const;

///Get loop kind of INACTIVE, PERIMITER_LOOP, OUTER_LOOP, or INNER_LOOP.
///get_kind will not get the kind of NETWORK.
void get_kind(const MGFSurface* srf=0) const;

///Get perimeter number where start or end point lies on.
///Function's return value is:
/// true if both end on perimeter,
/// false if any of start or end point not on perimeter.
///pid_s, _e are valid only when true is returned.
bool get_perimeter_num(int& pid_s, int& pid_e,const MGFSurface* srf=0) const;

///Test if (u,v) is in, on, or out of  curves by obtaining intersetion points
///of sl with the curves. The loop has a direction and if this is
///inner, inside means outside of the closed curves.
///Returned is:
///  0:outside(not on the loop)
///  1:unknown
///  2:inside(not on the loop)
/// otherwise:on the loop(int(parameter edge id +100))will be returned.
///When 1 is returned, try again by changing the sample straight line sl.
int in_on_curves(
	const MGPosition& uv,	///<Test point
	const double error[3],
		///<error[0]: Tolerance allowed. Use for the test of on curve.
		///<error[1]: loop's u span length
		///<error[2]: loop's v span length
	const std::vector<UniqueCurve>& crvs,///<Vector of curves that make loop.
		///<This has the direction as a face boundary, and the in-ness is
		///<determined according to the direcion.
	const MGStraight& sl	///<Sample straight line to test the inclusion.
)const;	

///Test if (u,v) is in, on, or out of this loop.
///inside_test() does not use any face data of this loop. That is, this loop
///does not have to have parent face. The loop has a direction and if this is
///inner, inside means outside of the closed curves.
///This loop does not have to make a closed loop.
///returned is:
///  0:outside(not on the loop)
///  1:unknown
///  2:inside(not on the loop)
/// otherwise:on the loop(int(parameter edge id +100))will be returned.
int inside_test(
	const MGPosition& uv	///<Point to test if inside of the loop.
)const;

///Compute intersection points of 1D sub curves of the original loop.
///Parameter values of intersection points(MGLEPoint's) will be returned.
///This is for tessellation and intersection with perimeter boudary edges will
///be excluded.
std::vector<MGLEPoint> isect_1D_tess(						
	double f,			///< Coordinate value
	int coordinate=0	///< Coordinate kind of the data f(from 0).
) const;	

///Compute intersections of this loop's binder edges with face.
///If some parameter edges did not have a binder, isect_binder generates them.
MGLSPoint_vector isect_binder(const MGFSurface& face) const;

///Make loop
void make_loop(std::unique_ptr<MGCurve>& param_curve, std::unique_ptr<MGCurve>& world_curve);

///Make a vertex at lp and subdivide the edge into two edges.
///Returned is true if subdivision is done and false if no subdivision
///is done since lp was one of existed vertex.
///When function's return value is true, pre is always the same edge as
///lp's edge and the iterator is unchanged.
bool makeVertexWithLcis(
	const MGLEPoint& lp,///<point to subdivide of this loop.
	MGLCisect_vector& lcis,///<MGLCisect_vector before makeVertex,
				///will be so updated according to the vertex make.
	MGPVertex*& pvertex ///MGPVertex of the made vertex.
);

///Internal program for merge functions.
bool merge_trim(const MGCurve& param_curve, const MGInterval& range1,
	   const MGCurve* world_curve);

///Internal funtion for merge_trim(loop).
bool merge_loop(const MGLoop& loop2);

///Dedicated function of trim_start, end, will set new parameter
///range of the trimming edge.
///When start=true, set as start point data, else end point data.
void trim_param_set(const MGLEPoint& t, bool start);

friend class MGFace;
friend class MGEdge;

};

///Get edge pointer from its iterator in MGComplex(MGLoop).
MG_DLL_DECLR MGEdge* edge_from_iterator(MGComplex::iterator i);
MG_DLL_DECLR const MGEdge* edge_from_iterator(MGComplex::const_iterator i);

///When nderi=0, get a parameter (u,v) of the surface at the boundary point.
MG_DLL_DECLR MGVector eval(const MGLEPoint& t, int nderi = 0);

MG_DLL_DECLR size_t inside_outer_loop(
	const MGPosition& uv,
	const std::vector<const MGLoop*>& loop,
	const MGSurface* surf = 0
);

///Build networks of surf, given parameter curves vector.
extern MG_DLL_DECLR void build_networks(
	const MGFSurface& surf,		///<The objective surface
	const std::vector<UniqueCurve>& pcurves,///<(u,v) 2D parameter curves of surf.
	std::vector<UniqueLoop>& networks	///<Built networks
);

///Get the clone of ts.loop()(==te.loop()), and trim th loop from ts to te.
///Function's return value is the trimmed loop.
extern MG_DLL_DECLR std::unique_ptr<MGLoop> trim_out_subloop(const MGLEPoint& ts, const MGLEPoint& te);

/** @} */ // end of TOPO group
#endif
