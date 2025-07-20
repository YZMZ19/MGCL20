#ifndef _mgTL2Fans_HH_
#define _mgTL2Fans_HH_

#include <vector>
#include "Tl2/TL2Fan.h"
#include "Tl2/TL2Triangle.h"

/****************************************************************/
/*   Copyright (c) 2019 by System fugen G.K.                */
/*                       All rights reserved.                   */
/****************************************************************/

///////////// mgTL2Fans /////////////

class MGPosition;
class MGLoop;
class mgTL2parameter;
class mgTL2Polyline;
class mgTL2FanEdges;
class mgTL2LPline;
class mgTL2Triangles;

/** @file */
/** @addtogroup UseTessellation
 *  @{
 */


///Defines a vector mgTL2Fan.

///  private class for tessellation.
class mgTL2Fans{
public:
	using MYELM = std::unique_ptr< mgTL2Fan>;
	using MYVEC = std::vector<MYELM>;
	typedef MYVEC::iterator iterator;
	typedef MYVEC::const_iterator const_iterator;

friend std::ostream& operator<< (std::ostream& out, const mgTL2Fans& fans);

/////////////constructor////////////////////////
mgTL2Fans(
	const MGLoop& polygon		///<The target polygon, that is, a outer loop of a MGFace.
);

mgTL2Fans(
	std::vector<const mgTL2Polyline*>& polylines//Edges of the polyline that are mgTL2Polyline.
);

mgTL2Fans(
	const mgTL2LPline pline[4]///Four edges that constitute a closed polygon.
);

mgTL2Fans(
	const mgTL2Polyline* pline[4]///Four edges that constitute a closed polygon.
);

///Access to i-th element.
const MYELM& operator[](int i)const{return m_fans[i];};
MYELM& operator[](int i){return m_fans[i];};

/////////////member function////////////////////////

	///Begin, end.
	const_iterator begin()const{return m_fans.begin();};
	const_iterator end()const{return m_fans.end();};
	iterator begin(){return m_fans.begin();};
	iterator end(){return m_fans.end();};


	///3�_�ڂ̒��_(id of m_fans)�����߂�(���_�͎g�p�_�ƂȂ�)
	int find3rdV(
		int		alpha,	///<�G�b�W�̎n�_(id of m_fans)
		int		beta,	///<�G�b�W�̏I�_(id of m_fans)
		int&		status	///<�X�e�[�^�X
	);

	///Fan�����ɕK�v�ȕϐ��̏������s��
	///edgeStack��edge��stack��ς�
	void init_edgeStack(
		mgTL2FanEdges& edgeStack
	);

	///Test if the edge(alpha, beta) is boundary or not.
	bool is_boundary(int alpha, int beta) const;

	///����(v1,v2)��������edge�ƌ�_�����邩�ǂ����𒲂ׂ�
	///Function's return: true if had isect.
	bool has_isect(
		int 	v1,	///<����1�̓_
		int 	v2	///<����1�̓_
	)const;

	///Push back a fan.
	void push_back(mgTL2Fan* fan){m_fans.emplace_back(fan);};
	void push_back(mgTL2Fans&& fans){
		std::move(fans.begin(), fans.end(),std::back_inserter(m_fans));};

	///�ړI�F���S�_(center)��alpha��mgTL2Fan�ɑ΂�����_gamma���beta��
	///���ɒǉ�����
	void push1Vaft(
		int	alpha,	///<���S�_�̃C���f�b�N�X
		int	beta,	///<�ǉ����钸�_�̃C���f�b�N�X
		int	gamma	///<��ƂȂ钸�_�̃C���f�b�N�X
	);

	///�ړI�F���S�_(center)��alpha��mgTL2Fan�ɑ΂�����_beta���gamma��
	///�O�ɒǉ�����
	void push1Vbefore(
		int	alpha,	///<���S�_�̃C���f�b�N�X
		int	beta,	///<�ǉ����钸�_�̃C���f�b�N�X
		int	gamma	///<��ƂȂ钸�_�̃C���f�b�N�X
	);

	///�ړI�F���S�_(center)��gamma�Œ��_(alpha,beta)�̂��̂�V�K�ɍ쐬����
	void push2V(
		int	gamma,	///<���̃C���f�b�N�X
		int	alpha,	///<���̃C���f�b�N�X
		int	beta	///<���̃C���f�b�N�X
	);

	///Set edge(alpha,j) as used.
	void set_edge_used(int alpha, int beta);

	///Get the number of fans included.
	int size() const{return (int)m_fans.size();};

	///Triangulate polygon.
	///The result will be appended onto triangles.
	void triangulate(
		mgTL2Triangles& triangles	///<Triangulated data will be appended.
	)const;

	///check if vertex(alpha) is used or not.
	bool used(int alpha) const;

	///check if edge(alpha, beta) is used or not.
	bool used(int alpha, int beta) const;

	///Retrieve surface parameter value (u,v) of the i-th vertex.
	MGPosition uv(int i)const;

	///Retrieve world coordinate value (x,y,z) of the i-th vertex.
	MGPosition xyz(int i, bool need_normal)const;

////////////////member data///////////////
private:
	MYVEC m_fans;
	std::vector<mgTL2LPline> m_polylines;///<Edges of the polyline that are mgTL2LPline.

	//Construct mgTL2Fans, that is, construct m_fans from m_polylines.
	void initialize();
};

///Triangulate polygon.
///The result will be appended onto triangles.
void triangulate(
	const MGLoop& polygon,///<Target MGLoop to triangulate whose edges' base_curve() must be
		///mgTL2Polyline.
	mgTL2Triangles& triangles///<Triangulated data will be appended.
);

///Triangulate polygon.
///The result will be appended onto triangles.
void triangulate(
	const mgTL2LPline polygon[4],///<Four edges that constitute a closed polygon.
	mgTL2Triangles& triangles ///<Triangulated data will be appended.
);

///Triangulate polygon.
///The result will be appended onto triangles.
void triangulate(
	const mgTL2Polyline* polygon[4],///<Four edges that constitute a closed polygon.
	mgTL2Triangles& triangles///<Triangulated data will be appended.
);

/** @} */ // end of UseTessellation group
#endif
