#include "StdAfx.h"
#include "Tl2/TL2Triangle.h"
#include "Tl2/TL2Fan.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
	
//Test if the edge(i,j) is used or not where i is the index of this edge.
bool mgTL2Fan::edge_is_used(int j)const{
	CEUitr l,ls=m_used_edges.begin(), le=m_used_edges.end();
	l=std::find(ls,le, j);
	return l!=le;
}

//頂点周辺の頂点リストからindexを検索する（後ろから検索）
//If found, iterator of the index be returned.
//If not found , end() will be returned.
mgTL2Fan::IndexItr mgTL2Fan::find_aft(int index){
	int n=size();
	for(int j=n-1; j>=0; j--){
		if(m_indices[j]==index) return begin()+j;
	}
	return end();
}

//Print out indices as "|n0,n1,....
void mgTL2Fan::print_indices(std::ostream& out)const{
	int n=size();
	int nm1=n-1;
	for(int i=0; i<n; i++){
		int id=(*this)[i];
		out<<id;
		if(i<nm1)
			out<<",";
	}
}

//Set the edge(i,j) as used where i is the index of this fan's vertex.
void mgTL2Fan::set_edge_used(int j){
	EUitr l,ls=m_used_edges.begin(), le=m_used_edges.end();
	l=std::find(ls,le, j);
	if(l==le)
		m_used_edges.push_back(j);
}

std::ostream& operator<< (std::ostream& out, const mgTL2Fan& fan){
	int n=fan.size();
	out<<"Fan::num of indices="<<n<<"::";
	if(n)
		out<<fan.m_indices[0];
	for(int i=1; i<n; i++)
		out<<","<<fan.m_indices[i];
	return out;
}
