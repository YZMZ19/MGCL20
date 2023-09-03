/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "topo/LEPoint.h"
#include "topo/BVertex.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "topo/Shell.h"
#include "mgGl/Appearance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGShell Class.

///////Constructor////////

//Construct a shell of one face.
//Copy version of face.
MGShell::MGShell(const MGFace& face){
	MGFace* nf=new MGFace(face);
	nf->make_outer_boundary();
	nf->remove_appearance();
	append_pcell(nf);
	copy_appearance(face);
}

//Construct a shell of one face.
//face is a pointer to a newed face, and the ownership
//will be transfered to MGShell.
MGShell::MGShell(MGFace* face){
	copy_appearance(*face);
	face->remove_appearance();
	face->make_outer_boundary();
	append_pcell(face);
}

//Construct a shell of one face.
//Copy version of face.
MGShell::MGShell(const MGFSurface& face){
	copy_appearance(*(face.object_pointer()));
	MGFace* f=face.clone_as_face();
	f->remove_appearance();
	f->make_outer_boundary();
	append_pcell(f);
}

//Construct a shell of one face.
//face is a pointer to newed face, and the ownership
//will be transfered to MGShell.
MGShell::MGShell(MGFSurface* face){
	copy_appearance(*(face->object_pointer()));
	MGFace* f=face->make_face();
	f->remove_appearance();
	f->make_outer_boundary();
	append_pcell(f);
}

//Fundamental constructor.
//Construct from boundary complex(i.e. MGLoop).
//This constructor takes the ownership of MGCell* in boundary.
MGShell::MGShell(
	std::list<MGCell*> faces
		//Boundary data of the super class MGComplex(List of faces).
):MGComplex(faces){;}

///////Destructor///////
MGShell::~MGShell() { ; }

///////operator overload///////


/// Complexに平行移動を行ないオブジェクトを生成する。
///Translation of the Complex
MGShell MGShell::operator+ (const MGVector& v) const{
	MGShell shell2(*this);
	shell2 += v;
	return shell2;
}

//Complexのスケーリングを行い，Complexを作成する。
//Scaling of the Complex by a double.
MGShell MGShell::operator* (double s) const{
	MGShell shell2(*this);
	shell2 *= s;
	return shell2;
}

/// 与えられた変換でComplexの変換を行い，Complexを作成する。
///Transformation of the Complex by a matrix.
MGShell MGShell::operator* (const MGMatrix& mat) const{
	MGShell shell2(*this);
	shell2 *= mat;
	return shell2;
}

/// 与えられた変換によってトランスフォームをおこないComplexを生成する。
///Transformation of the Complex by a MGTransf.
MGShell MGShell::operator* (const MGTransf& tr) const{
	MGShell shell2(*this);
	shell2 *= tr;
	return shell2;
}

///comparison
bool MGShell::operator==(const MGShell& s2)const{
	return this == &s2;
}
std::partial_ordering MGShell::operator<=>(const MGShell& s2)const {
	return this <=> &s2;
}

bool MGShell::equal_test(const MGGel& g2)const {
	auto c = typeCompare(g2);
	return c == 0 ? *this == dynamic_cast<const MGShell&>(g2) : false;
}
std::partial_ordering MGShell::ordering_test(const MGGel& g2)const {
	auto c = typeCompare(g2);
	return c == 0 ? *this <=> dynamic_cast<const MGShell&>(g2) : c;
}

MGShell& MGShell::operator=(const MGGel& gel2){
	const MGShell* gel2_is_this=dynamic_cast<const MGShell*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

///////Member Function///////

//Make a clone.
MGShell* MGShell::cloneWithParent(MGCell& parent) const{
	MGShell* shl=clone();
	shl->set_parent(parent);
	return shl;
}
MGShell* MGShell::clone() const{
	MGShell* shl=new MGShell(*this);
	return shl;
}

//Exlude non free edges or opposite direction's edges from their vectors.
void exclude(
	std::vector<MGLEPoint>& pr1,
	std::vector<MGLEPoint>& pr2
){
	if(pr1.empty())
		return;

	//Exclude non free edges, or opposite direction's edges.
	size_t n=pr1.size();
	for(int i=(int)(n-1); i>=0; i--){
		std::vector<MGLEPoint>::iterator i1=pr1.begin()+i, i2=pr2.begin()+i;
		if( (i1->eval_star(1)%i2->eval_star(1))>0. || !(i1->edge()->is_free())
			|| !(i2->edge()->is_free()) ){
			pr1.erase(i1); pr2.erase(i2);
		}
	}
}

//Get face pointer from its iterator in MGComplex(MGShell).
MGFace* MGShell::face(iterator i){
	return static_cast<MGFace*>(i->get());
}
const MGFace* MGShell::face(const_iterator i)const{
	return static_cast<const MGFace*>(i->get());
}

//Get the face pointer from its pcell id in MGComplex.
MGFace* MGShell::face(int i){
	iterator j=pcell_begin();
	std::advance(j,i);
	return face(j);
}
const MGFace* MGShell::face(int i)const{
	const_iterator j=pcell_begin();
	std::advance(j,i);
	return face(j);
}

//This is a proprietry routine for merge_at_common_edge, will merge
//loops at free common edges.
bool merge_2loops(bool can_negate, MGLoop& lp1, MGLoop& lp2){
	std::vector<MGLEPoint> pr1,pr2;
	std::vector<double> br1,br2;
	if(!lp1.common(lp2,pr1,br1,pr2,br2)) return false;

	if(can_negate && pr1[0].eval_star(1)%pr2[0].eval_star(1)>0.){
		//negate the face.
		lp2.face()->negate();
		lp1.common(lp2,pr1,br1,pr2,br2);
	}
	exclude(pr1,pr2);
	//Exclude pairs of edges that are not the same direction as the first edge.
	if(pr1.empty()) return false;

	std::vector<MGEdge*> e1s=lp1.subdivide(pr1);
	std::vector<MGEdge*> e2s=lp2.subdivide(pr2);
	int ne=(int)e1s.size();
	for(int i=0; i<ne; i++)
		e1s[i]->bind(*(e2s[i]));
	return true;
}

//Merge a face at free edges of this shell.
//Function's return value is
//   false: merge was not done because no common edges were found.
//   true: merge of the shell and the face was done.
//      Case of true includes that merge was done only partialy.
//      This case occurs when after some common edges are merged,
//      some common edges are found to have contradictionary direction.
//      Merge is done only for the first common edges found to have
//      same direction.
//When the function's return value is false, the face will not be added to
//the shell. Merge negates the input face direction when this shell's direction
//is not the same as the input face's one along the common edge(s).
//The second form is a pointer version. The face must be newed object
//and will be destructed after merge when function's return value is 0
//(merge was processed).
bool MGShell::merge_at_common_edge(const MGFace& face){
	MGFace* nf=new MGFace(face);
	bool merged=merge_at_common_edge(nf);
	if(!merged) delete nf;
	return merged;
}
bool MGShell::merge_at_common_edge(const MGFSurface& face){
	MGFSurface* nf=face.clone_fsurface();
	bool merged=merge_at_common_edge(nf);
	if(!merged) delete nf;
	return merged;
}
bool MGShell::merge_at_common_edge(MGFSurface* face){
	MGFace* f=dynamic_cast<MGFace*>(face);
	if(f)
		return merge_at_common_edge(f);

	MGSurface* srf=dynamic_cast<MGSurface*>(face);
	MGFace* f2=new MGFace(srf);
	if(merge_at_common_edge(f2))
		return true;
	delete f2;
	return false;
}

//This form is a pointer version. The face must be newed object
//and will be destructed after merge when function's return value is 0
//(merge was processed).
//This shell may be dummy shell. In this case, the face is added to the 1st
//face in the shell.
bool MGShell::merge_at_common_edge(MGFace* face){
	int i1,i2, i1save,i2save; //id of inner boudary of this shell's loop and face.
	int nin1,nin2; //number of inner boudaries of this shell's loop and face.
	int j1,j2;	//inner boudary loop variable.

	face->make_outer_boundary();
	UniqueLoop& out2=face->loop(int(0));
	nin2=face->number_of_inner_boundaries(i2save);

	iterator ci=pcell_begin(), ce=pcell_end();
	if(ci==ce){
		copy_appearance(*face);
		face->remove_appearance();
		append_pcell(face);
		return true;
	}

	bool merged=false, merged2;
	for(;ci!=ce; ci++){
		MGFace* fi=dynamic_cast<MGFace*>(ci->get());
		if(!fi) continue;
		if(fi==face) continue;

		// 1. outer1 versus outer2.
		UniqueLoop& out1=fi->loop(int(0));
		merged2=merge_2loops(!merged, *out1, *out2);
		if(merged2){
			merged=true;
		}

		// 2. outer1 versus inner2.
		i2=i2save;
		for(j2=0; j2<nin2; j2++, i2++){
			merged2=merge_2loops(!merged,*out1,*(face->loop(i2)));
			if(merged2) merged=true;
		}

		nin1=fi->number_of_inner_boundaries(i1save);
		i1=i1save;
		for(j1=0; j1<nin1; j1++, i1++){
		// 3. inner1 versus outer2.
			UniqueLoop& inner1=fi->loop(i1);
			merged2=merge_2loops(!merged,*inner1, *out2);
			if(merged2) merged=true;

		// 4. inner1 versus inner2.
			i2=i2save;
			for(j2=0; j2<nin2; j2++, i2++){
				merged2=merge_2loops(!merged,*inner1,*(face->loop(i2)));
				if(merged2) merged=true;
			}
		}
	}
	if(merged)
		face->remove_appearance();

	return merged;
}

//Debug Function
std::ostream& MGShell::toString(std::ostream& ostrm) const{
	ostrm<<"<<MGShell="<<(const MGGel*)this;
	MGComplex::toString(ostrm);
	ostrm<<"=MGShell>>"<<std::endl;
	return ostrm;
}

//Obtain boundary and main parameter lines of the FSurface.
//skeleton includes boundary() and inner parameter lines.
//density indicates how many inner parameter lines are necessary
//for both u and v directions.
std::vector<UniqueCurve> MGShell::skeleton(int density) const{
	const_iterator ci=pcell_begin(), ce=pcell_end();
	std::vector<UniqueCurve> crvs;
	for(;ci!=ce; ci++){
		MGFace* fi=static_cast<MGFace*>(ci->get());
		std::vector<UniqueCurve> si=fi->skeleton(density);
		std::move(si.begin(), si.end(), std::back_inserter(crvs));
	}
	return crvs;
}

//Obtain all the parameter curves at knots of u and v knot vector.
std::vector<UniqueCurve> MGShell::skeleton_at_knots()const{
	const_iterator ci=pcell_begin(), ce=pcell_end();
	std::vector<UniqueCurve> crvs;
	for(;ci!=ce; ci++){
		MGFace* fi=dynamic_cast<MGFace*>(ci->get());
		std::vector<UniqueCurve> si=fi->skeleton_at_knots();
		std::move(si.begin(), si.end(), std::back_inserter(crvs));
	}
	return crvs;
}
