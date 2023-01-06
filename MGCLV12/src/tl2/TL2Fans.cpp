#include "StdAfx.h"
#include "mg/Tolerance.h"
#include "mg/Position.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "Tl2/TL2Parameter.h"
#include "Tl2/TL2LPline.h"
#include "Tl2/TL2Fan.h"
#include "Tl2/TL2Fans.h"
#include "Tl2/TL2Triangle.h"
#include "Tl2/TL2Triangles.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

/****************************************************************/
/*   Copyright (c) 2019 by System fugen G.K.                */
/*                       All rights reserved.                   */
/****************************************************************/

//////////// mgTL2Fans///////////

////////// private class //////////

class mgTL2FanEdge{
private:
	int m_start, m_end;//Edge's start and end ids.
public:
	mgTL2FanEdge(){;};
	mgTL2FanEdge(int start, int end):m_start(start),m_end(end){;};
	int start()const {return m_start;};
	int end()const {return m_end;};
};

class mgTL2FanEdges{
private:
	std::deque<mgTL2FanEdge> m_edges;
public:
	typedef std::deque<mgTL2FanEdge>::const_iterator eitr;
	mgTL2FanEdges(){;};
	eitr begin()const{return m_edges.begin();};
	eitr end()const{return m_edges.end();};
	bool empty()const{return m_edges.empty();};
	void pop_back(){m_edges.pop_back();};
	void push_back(int start, int end){
		m_edges.push_back(mgTL2FanEdge(start,end));
	};
	void push_front(int start, int end){
		m_edges.push_front(mgTL2FanEdge(start,end));
	};
	void push_front(const mgTL2FanEdge& edge){m_edges.push_front(edge);};
};

std::ostream& operator<< (std::ostream& out, const mgTL2FanEdges& edges){
	mgTL2FanEdges::eitr i=edges.begin(), ie=edges.end();
	for(int j=0; i!=ie; i++, j++){
		out<<" "<<j<<"("<<(*i).start()<<","<<(*i).end()<<")";
	}
	return out;
}

//ファンを作成する際に使用するステータス
enum mgTRIANG_STATUS{
	UNKNOWN,	//0
	TWOTOUCH,	//1
	NOTOUCH,	//2
	REGULAR,	//3
	RIGHTTOUCH,	//4
	LEFTTOUCH	//5
};

///////////constructor////////////

mgTL2Fans::mgTL2Fans(
	const MGLoop& polygon		///<The target polygon, that is, the outer loop of an MGFace.
){
	int npoly=polygon.number_of_edges();
	m_polylines.resize(npoly);
	for(int j=0; j<npoly; j++){
		const mgTL2Polyline* pline=TL2Polyline(polygon.edge(j));
		m_polylines[j]=mgTL2LPline(pline);
	}
	initialize();
}

mgTL2Fans::mgTL2Fans(
	std::vector<const mgTL2Polyline*>& polylines//Edges of the polyline that are mgTL2Polyline.
){
	size_t npoly=polylines.size();
	m_polylines.resize(npoly);
	for(size_t j=0; j<npoly; j++)
		m_polylines[j]=mgTL2LPline(polylines[j]);
	initialize();
}

mgTL2Fans::mgTL2Fans(
	const mgTL2LPline pline[4]///Four edges that constitute a closed polygon.
){
	int npoly=4;
	m_polylines.resize(npoly);
	for(int j=0; j<npoly; j++)
		m_polylines[j]=pline[j];
	initialize();
}

mgTL2Fans::mgTL2Fans(
	const mgTL2Polyline* pline[4]///Four edges that constitute a closed polygon.
){
	int npoly=4;
	m_polylines.resize(npoly);
	for(int j=0; j<npoly; j++)
		m_polylines[j]=mgTL2LPline(pline[j]);
	initialize();
}


//Construct mgTL2Fans, that is, construct m_fans from m_polylines.
void mgTL2Fans::initialize(){
	size_t nvertices=0;
	size_t npoly=m_polylines.size();
	for(size_t j=0; j<npoly; j++){
		nvertices+=(m_polylines[j].number_of_points()-1);
	}
	if(nvertices<3)
		return;	//一応チェックしておく

	//多角形の辺をスタックに積む
	m_fans.resize(nvertices);
	mgTL2FanEdges edgeStack;		//スタックエッジ
	init_edgeStack(edgeStack);
	if(nvertices==3)
		return;

	//スタックがからになるまで、頂点と周辺の頂点リストのベクトルを求める
	size_t nv2=nvertices*2;
	for(size_t i=0; !edgeStack.empty() && i<nv2;i++){
		//スタックからポップする
		mgTL2FanEdges::eitr endEdgeIter = edgeStack.end();
		const mgTL2FanEdge edge = *(--endEdgeIter); edgeStack.pop_back();
		int alpha=edge.start(), beta=edge.end();
		if(is_boundary(alpha,beta)){
			if(used(alpha,beta)) continue;
		}

		//3点目の頂点を求める
		int status;
		int gamma = find3rdV(alpha, beta, status);
		//求まったステータスに応じた処理を行う
		if(status==UNKNOWN)continue;						//UNKNOWN==0
		switch(status){
			case TWOTOUCH: continue;						//TWOTOUCH==1
			case NOTOUCH: edgeStack.push_front(edge); break;//NOTOUCH==2
			case REGULAR:									//REGULAR==3
				push1Vaft(alpha, beta, gamma);
				push1Vbefore(beta, gamma, alpha);
				push2V(gamma, alpha, beta);
				set_edge_used(alpha, beta);
				set_edge_used(beta, gamma);
				set_edge_used(gamma, alpha);
				break;
			case RIGHTTOUCH:								//RIGHTTOUCH==4
				push1Vbefore(beta, gamma, alpha);
				push1Vaft(gamma, alpha, beta);
				set_edge_used(alpha, beta);
				set_edge_used(beta, gamma);
				break;
			case LEFTTOUCH:									//LEFTTOUCH==5
				push1Vaft(alpha, beta, gamma);
				push1Vbefore(gamma, alpha, beta);
				set_edge_used(alpha, beta);
				set_edge_used(gamma, alpha);
				break;
			default: break;//求まらなかったとき				
		}
		if(!is_boundary(alpha,gamma) && (status!=RIGHTTOUCH)){
			edgeStack.push_back(alpha, gamma);
		}
		if(!is_boundary(gamma,beta) && (status!=LEFTTOUCH)){
			edgeStack.push_back(gamma, beta);
		}
	}
}

//3点目の頂点(id of m_fans)を求める
int mgTL2Fans::find3rdV(
	int		alpha,	//エッジの始点(id of m_fans)
	int		beta,	//エッジの終点(id of m_fans)
	int&		status	//ステータス
){
	const int zcoord=2;
	const MGPosition sPos=uv(alpha);
	const MGPosition ePos=uv(beta);
	MGVector cedge(ePos-sPos);

	//多角形をループして3点目を求める
	bool gamma_is_used, right_touch, left_touch;
	double maxCang=2.0;
	int gamma=0;
	int nvertices=size();
	for(int i=0; i<nvertices; i++){
		if((alpha==i) || (beta==i))
			continue;//始終点と等しいときパスする

		//外積が負ではいけない
		const MGPosition p3 = uv(i);//3点目の座標
		if((cedge*MGVector(p3-sPos))[zcoord] <= 0.0)
			continue;

		bool i_is_used, iright_touch=false, ileft_touch=false;
		double cang((sPos-p3).cangle(ePos-p3));
		if(cang < maxCang){//When the angle at i is larger than the before,
			//check if the edges (i,alpha) and (i,beta) do not have
			//any intersections with existing edges.
			i_is_used=used(i);
			if(i_is_used){
				iright_touch=used(i,alpha);
				if(!iright_touch){//共通エッジがあれば交点がない
					if(has_isect(i,alpha))
						continue;
				}
				ileft_touch=used(beta,i);
				if(!ileft_touch){//共通エッジがあれば交点がない
					if(has_isect(beta,i))
						continue;
				}
			}else{
				if(has_isect(i,alpha))
					continue;
				if(has_isect(beta,i))
					continue;
			}
			maxCang = cang;
			gamma_is_used=i_is_used;
			right_touch = iright_touch;
			left_touch = ileft_touch;
			gamma=i;
		}
	}
	if(maxCang > 1.5){//3点目が見つからなかった
		status=UNKNOWN;
		return gamma;
	}

	//ステータスの更新を行う
	if(!gamma_is_used){//未使用の点のときREGULAR
		status=REGULAR;
	}else if(right_touch){
		if(left_touch) status=TWOTOUCH;
		else status = RIGHTTOUCH;
	}else{
		if(left_touch) status=LEFTTOUCH;
		else status=NOTOUCH;
	}
	return gamma;
}

//線分(v1,v2)が既存のedgeと交点があるかどうかを調べる
//Function returns true if (v1,v2) had an isect.
bool mgTL2Fans::has_isect(
		int 	v1,	//線分1の点
		int 	v2	//線分1の点
)const{
	if(is_boundary(v1,v2))
 return false;

	double error=MGTolerance::rc_zero();
	const MGPosition p1=uv(v1), p2=uv(v2);
	const MGVector dir1(p2-p1);
	double udir1=dir1[0], vdir1=dir1[1];
	int nfan=size()-1;
	for(int i=0; i<nfan; i++){//Loop over fans(m_fans).
	if(v1==i || v2==i)
		continue;

	const MGPosition p3=uv(i);//Center of the fani.
	const mgTL2Fan& fani=*(m_fans[i]);
	int nv=fani.size();
	MGVector dir31(p3-p1), dir32(p3-p2);
	double udir31=dir31[0], vdir31=dir31[1];
	double udir32=dir32[0], vdir32=dir32[1];
	double z13=vdir31*udir1-udir31*vdir1;
		//z value of vector product of dir31 and dir1
	for(int j=0; j<nv; j++){//Loop over vertices on the fani.
		int vj=fani[j];
		if(vj<=i)
			continue;
		if(vj==v1 || vj==v2)
			continue;
//		if(is_boundary(i,vj)) continue;
		if(!is_boundary(i,vj) && !used(i,vj))
			continue;

		const MGPosition p4=uv(vj);
		MGVector dir41(p4-p1);
		double z14=dir41[1]*udir1-dir41[0]*vdir1;
			//z value of vector product of dir41 and dir1
		double z134=z13*z14;
		if(z134>=0.)
			continue;//This means p3 and p4 are located at the same side about the straight line (p1, p2).
//		if(-z134<=error) continue;///////////

		const MGVector dir2(p4 - p3);
		double udir2=dir2[0], vdir2=dir2[1];
		double z23=vdir31*udir2-udir31*vdir2;
		double z24=vdir32*udir2-udir32*vdir2;
		double z234=z23*z24;
		if(z234>=0.)
			continue;//This means p1 and p2 are located at the same side about the
							//straight line (p3, p4).
//		if(-z234<=error) continue;///////////

		return true;//In this case, (p1,p2) and (p3,p4) has and intersection.
	}

	}
	return false;
}

//Fan生成に必要な変数の準備を行う
//edgeStackにedgeのstackを積む
void mgTL2Fans::init_edgeStack(
	mgTL2FanEdges& edgeStack
){
	//多角形の辺をスタックに積む、このとき頂点を未使用にしておく
	//push edges on the stack edgeStack.
	int nvertices=size();
	int im1=nvertices-1;
	for(int i=0; i<nvertices; i++){
		int ip1=i+1;
		if(ip1==nvertices)
			ip1=0;
		m_fans[i].reset(new mgTL2Fan(ip1,im1));
		edgeStack.push_back(i,ip1);
		im1=i;
	}
}

//Test if the edge(alpha, beta) is boundary or not.
bool mgTL2Fans::is_boundary(int alpha, int beta) const{
	if(alpha>beta){
		int temp=alpha; alpha=beta; beta=temp;
	}
	int n=size();
	if(alpha==0 && beta==n-1)
		return true;
	int alpha_n=(alpha+1)%n;
	return alpha_n==beta;
}

//目的：中心点(center)がalphaのmgTL2Fanに対し､頂点gammaを基準betaの
//後ろに追加する
void mgTL2Fans::push1Vaft(
	int	alpha,	//中心点のインデックス
	int	beta,	//基準となる頂点のインデックス
	int	gamma	//追加する頂点のインデックス
){
	mgTL2Fan& fan=*(m_fans[alpha]);
	mgTL2Fan::IndexItr j=fan.find_aft(beta);
	if(j!=fan.end()){//If not found
		j++;
		if(j!=fan.end()&&(*j)==gamma) return;
	}
	fan.insert(j,gamma);
	fan.set_vertex_used();
}

//目的：中心点(center)がalphaのmgTL2Fanに対し､頂点betaを基準gammaの
//前に追加する
void mgTL2Fans::push1Vbefore(
	int	alpha,	//中心点のインデックス
	int	beta,	//追加する頂点のインデックス
	int	gamma	//基準となる頂点のインデックス
){
	mgTL2Fan& fan=*(m_fans[alpha]);
	mgTL2Fan::IndexItr j=fan.find(gamma);
	if(j!=fan.end()){//If not found
		if(j!=fan.begin()){
			mgTL2Fan::IndexItr jm1=j;jm1--;
			if((*jm1)==beta) return;
		}
	}else{
		j=fan.begin();
	}
	fan.insert(j,beta);
	fan.set_vertex_used();
}

//目的：中心点(center)がalphaのfanに頂点(beta,gamma)を新規に作成する
//push2V is invoked only for unused vertices.
void mgTL2Fans::push2V(
	int	alpha,	//γのインデックス
	int	beta,	//αのインデックス
	int	gamma	//βのインデックス
){
	mgTL2Fan& fan=*(m_fans[alpha]);
	mgTL2Fan::IndexItr j=fan.begin();
	if((*j++)!=beta){
		fan.insert(j,beta);
	}
	mgTL2Fan::IndexItr je=fan.end(); j=je;j--;
	if(*j!=gamma) fan.insert(j,gamma);
	fan.set_vertex_used();
}

//Set edge(i,j) as used.
void mgTL2Fans::set_edge_used(int alpha, int beta){
	if(alpha>beta){
		int temp=alpha; alpha=beta; beta=temp;
	}
	MYELM& fan=m_fans[alpha];
	fan->set_edge_used(beta);
}

//check if vertex(alpha) is used or not.
bool mgTL2Fans::used(int alpha) const{
	return m_fans[alpha]->vertex_is_used();
}

//check if edge(alpha, beta) is used or not.
bool mgTL2Fans::used(int alpha, int beta) const{
	if(alpha>beta){
		int temp=alpha; alpha=beta; beta=temp;
	}
	const MYELM& fan=m_fans[alpha];
	return fan->edge_is_used(beta);
}

MGPosition mgTL2Fans::uv(int i)const{
	assert(i<size());
	int npoly=(int)m_polylines.size();
	for(int j=0; j<npoly; j++){
		const mgTL2LPline& pline=m_polylines[j];
		int nj=pline.number_of_points();
		if(i<nj)
			return pline.uv(i);
		i-=nj-1;
	}
	assert(false);
	return MGPosition(2);
}

MGPosition mgTL2Fans::xyz(int i, bool need_normal)const{
	assert(i<size());
	int npoly=(int)m_polylines.size();
	for(int j=0; j<npoly; j++){
		const mgTL2LPline& pline=m_polylines[j];
		int nj=pline.number_of_points();
		if(i<nj)
			return pline.xyz(i,need_normal);
		i-=nj-1;
	}
	assert(false);
	return MGPosition(3);
}

std::ostream& operator<< (std::ostream& out, const mgTL2Fans& fans){
	out<<"mgTL2Fans::num of fans="<<fans.size()<<std::endl;
	size_t npoly=fans.m_polylines.size();
	for(size_t j=0; j<npoly; j++)
		out<<j<<"-th edge:"<<fans.m_polylines[j]<<std::endl;
	mgTL2Fans::const_iterator i=fans.begin(), ie=fans.end();
	for(int j=0; i!=ie; i++,j++){
		out<<j<<"(";
		(**i).print_indices(out);
		out<<")"<<endl;
	}
	return out;
}

//Private class to sort the fans in mgTLFans according to the vertex number.
class mgTL2FanSize{
private:
	int m_center;	//id of mgTLFans(center id).
	int m_vnum;	//number of vertices of the mgTL2Fan.
public:
	int center()const{return m_center;};
	void set(int center, int vnum){m_center=center; m_vnum=vnum;};
	int vnum()const{return m_vnum;};
};

//mgTlTriangleをpolyの長さの順にソートするためのクラス
class mgTl2fansizeSort{
public:
	bool operator()(const mgTL2FanSize& tf1, const mgTL2FanSize& tf2)
		const{return tf1.vnum() > tf2.vnum();}
};

///Triangulate polygon.
///The result will be appended onto triangles.
void mgTL2Fans::triangulate(
	mgTL2Triangles& triangles	///<Triangulated data will be appended.
)const{
	const mgTL2Fans& fans=*this;

	//頂点と周辺の頂点リストのベクトルから3角形FANのベクトルを作成する
	int nfan=fans.size();
	std::vector<mgTL2FanSize> fansizes(nfan);
	for(int i=0; i<nfan; i++){
		fansizes[i].set(i,fans[i]->size());
	}
	std::sort(fansizes.begin(), fansizes.end(), mgTl2fansizeSort());
	std::vector<bool> vused(nfan,false);
		//Flag if the corresponding vertex is already processed to make fan.
		//vused[i] corrsponds to the vertex fans[i] is used or not.

	for(int ifan=0; ifan<nfan; ifan++){
		//Loop over fansizes vector. ifan is an id of fansizes.
		int center=fansizes[ifan].center();
		const mgTL2Fan& fan=*(fans[center]);
		int nvert=fan.size();
		if(nvert<=1)
			continue;//To process the next fan.

		std::vector<int> verticesIDs(nvert+1);//vertices id vector.
		int nvTri=1;
		verticesIDs[0]=center;
		for(int m=0; m<nvert; m++){
			int v1=fan[m];
			//現在頂点v1がすでに作成されたfanの一部であれば元のfanを分割する
			if(vused[v1]){
				triangles.push_back(fans,nvTri,verticesIDs);
				nvTri=1;
			}else{
				verticesIDs[nvTri++]=v1;
			}
		}
		triangles.push_back(fans,nvTri,verticesIDs);
		vused[center]=true;
	}
}

///Triangulate polygon.
///The result will be appended onto triangles.
void triangulate(
	const MGLoop& polygon,///<Target MGLoop to triangulate whose edges' base_curve() must be
		///mgTL2Polyline.
	mgTL2Triangles& triangles///<Triangulated data will be appended.
){
	//各多角形のファンを生成する
	mgTL2Fans fans(polygon);//全体の三角形FANのベクトル
	fans.triangulate(triangles);
}

///Triangulate polygon.
///The result will be appended onto triangles.
void triangulate(
	const mgTL2LPline polygon[4],///<Four edges that constitute a closed polygon.
	mgTL2Triangles& triangles ///<Triangulated data will be appended.
){
	mgTL2Fans fans(polygon);//全体の三角形FANのベクトル
	fans.triangulate(triangles);
}

///Triangulate polygon.
///The result will be appended onto triangles.
void triangulate(
	const mgTL2Polyline* polygon[4],///<Four edges that constitute a closed polygon.
	mgTL2Triangles& triangles///<Triangulated data will be appended.
){
	mgTL2Fans fans(polygon);//全体の三角形FANのベクトル
	fans.triangulate(triangles);
}
