#ifndef _mgTL2Fan_HH_
#define _mgTL2Fan_HH_

#include <algorithm>
#include <vector>
#include <iosfwd>

#include "mg/MGCL.h"

/****************************************************************/
/*   Copyright (c) 2019 by System fugen G.K.                */
/*                       All rights reserved.                   */
/****************************************************************/

///////////// mgTL2Fan /////////////

/** @file */
/** @addtogroup UseTessellation
 *  @{
 */


// private class for tessellation.

///mgTL2Fan is a point list to constitue a fan.

///It is always a member of mgTL2Fans that is temporary data to generate
///mgTL2Triangles, and contains id of vertices of the outer loop of a MGFace.
///mgTL2Fan does not include its start point id since the id in mgTL2Fans is the 1st id.
///Let mgTL2Fan& fani=mgTL2Fans[i], then (i,fani[0], fani[1], ..., fani[n-1]) constitutes
///a fan(all of them are ids of vertices of the outer loop of a MGFaces).
///The actual triangles are:
///(i,fani[0], fani[1]), (i,fani[1], fani[2]), ..., (i,fani[i], fani[i+1]), ...,
///(i, fani[n-2], fnai[n-1]).
///See Les Piegel and Wayne Tiller's paper about this point list.
///"Geometry-based triangulation of trimmed NURBS surfaces", Computer-Aided Desigh,
///Vol.30, No.1, pp.11-16, 1998.
class mgTL2Fan{
	typedef std::deque<int> mgTL2deqIndex;

private:
	mgTL2deqIndex	m_indices;///<includes vertices from the second points to the last.
	std::vector<int> m_used_edges;
					///<Include j's array if the edge(i,j) is a used edge in mgTL2Fans.
					///<j is always greater than i, and the order of j's is undefined.
					///<Here i is this fan's index.				
	bool m_used = false;	///<flag to indicate if this vertex is used or not.

public:
	typedef mgTL2deqIndex::iterator IndexItr;
	typedef mgTL2deqIndex::const_iterator CIndexItr;
	typedef mgTL2deqIndex::reverse_iterator ritr;
	typedef std::vector<int>::iterator EUitr;///<Edge Used iterator.
	typedef std::vector<int>::const_iterator CEUitr;///<Edge Used iterator.

friend std::ostream& operator<< (std::ostream& out, const mgTL2Fan& fan);

///////////////Constructor///////////////

	mgTL2Fan()=default;
	mgTL2Fan(int v1):m_indices(1,v1){;};
	mgTL2Fan(int v1, int v2):m_indices(2){m_indices[0]=v1;m_indices[1]=v2;};
	mgTL2Fan(mgTL2deqIndex& index):m_indices(index){;};

//////////// Operator overload ///////////////
	int operator[](int i)const{return m_indices[i];};

/////////////member functions////////////

	int back()const{return m_indices.back();};
	IndexItr begin(){return m_indices.begin();};
	CIndexItr begin()const{return m_indices.begin();};
	
	///Test if the edge(i,j) is used or not where i is the index of this edge.
	bool edge_is_used(int j)const;

	IndexItr end(){return m_indices.end();};
	CIndexItr end()const{return m_indices.end();};
	void erase(IndexItr iter){m_indices.erase(iter);};

	///頂点周辺の頂点リストからindexを検索する（前から検索）
	CIndexItr find(int index)const{
		return std::find(m_indices.begin(), m_indices.end(), index);
	}
	///頂点周辺の頂点リストからindexを検索する（前から検索）
	IndexItr find(int index){
		return std::find(m_indices.begin(), m_indices.end(), index);
	}

	///頂点周辺の頂点リストからindexを検索する（後ろから検索）
	IndexItr find_aft(int index);

	int front()const{return m_indices.front();};

	///Insert the index before the position iter.
	IndexItr insert(IndexItr iter, int index){
		return m_indices.insert(iter, index);
	};

	const mgTL2deqIndex& indices()const{return m_indices;};
	void push_back(int index){m_indices.push_back(index);};
	void push_front(int index){m_indices.push_front(index);};
	void pop_back(){m_indices.pop_back();};
	void pop_front(){m_indices.pop_front();};

	///Print out indices as "|n0,n1,....
	void print_indices(std::ostream& out)const;

	ritr rbegin(){return m_indices.rbegin();};
	ritr rend(){return m_indices.rend();};

	///Set this vertex as used.
	void set_vertex_used(){m_used=true;};

	///Set the edge(i,j) as used where i is the index of this fan's vertex.
	void set_edge_used(int j);

	int size()const{return (int)m_indices.size();};

	bool vertex_is_used()const{return m_used;};

};

/** @} */ // end of UseTessellation group
#endif
