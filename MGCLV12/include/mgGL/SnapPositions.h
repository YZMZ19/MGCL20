/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGPositions_HH_
#define _MGPositions_HH_

#include <vector>
#include <list>
#include <set>
#include "mg/Position.h"
#include "mg/PickObjects.h"
#include "mgGL/VBO.h"
#include "mgGL/Appearance.h"

//
//Define MGSnapPositions Class.

class MGCurve;
class MGGroup;
class MGGel;
class MGFace;
class MGShell;
/** @addtogroup GLAttrib
 *  @{
 */

///MGSnapPositions is a class to store candidate snap positions.

///Positions are stored in an array(vector) of MGPosition's.
///One snap kind candidate position date of all the target objects
///are extracted, and stored in a MGSnapPositions.
///This is to fasten the snap position picking.
class MG_DLL_DECLR MGSnapPositions: public mgVBO{

public:

///Define snap kind enum.
enum snap_kind{
	DELETE_POINT=-1,
	nopos=0,
	endpos,
	knotpos,
	vertexpos,
	nearpos,
	centerpos,
	ON_SURFACE
};

////////////Constructor////////////

///Void constructor.
MGSnapPositions(snap_kind kind=nopos);

///Copy constructor.
MGSnapPositions(const MGSnapPositions& sp);

///Virtual Destructor
virtual ~MGSnapPositions();

////////////Member Function////////////

typedef std::vector<MGPosition>::const_iterator const_iterator;
typedef std::vector<MGPosition>::iterator iterator;
typedef std::pair<const MGObject*, int> obj_num;

///append this positions in m_positions into points.
void append_position_data(std::vector<MGPosition>& points)const;

///Assignment
MGSnapPositions& operator=(const MGSnapPositions& sp);

///Extract candidate position data into this.
void extract(
	const MGCurve& crv	///<the curve to extract.
);
void extract(
	const MGSurface& srf///<the surface to extract.
);
void extract(
	const MGPoint& point	///<the curve to extract.
);
void extract(
	const MGFace& face	///<the curve to extract.
);
void extract(
	const MGShell& shell	///<the surface to extract.
);
void extract(
	const std::list<const MGGel*>& gel_list	///<the list to extract.
);
void extract(
	const MGGel& gel	///<the gel to extract.
);
void extract(
	const MGPickObjects& pobjs	///<array of pick objects to extract.
);

///Prepare for selection draw.
virtual void make_display_list(MGCL::VIEWMODE vmode=MGCL::DONTCARE);

///描画関数selectionDraw()は、Object選択のための表示処理をする。

///通常のdrawとの相違：///Colorとしてm_bufferIDを用い、size処理以外の
///attributesの処理（normal, texture, color)をしない。
// 派生もとのクラスはVirtual宣言されてるが、このクラスから派生クラスを作ることは
// 想定していないので、Virtualをはずしている。
void selectionDraw(MGCL::VIEWMODE viewMode=MGCL::DONTCARE);

/// Inputting output(selected) of pick_to_select_buf, obtains pick data.
void get_pick_data(
	const std::set<unsigned>& selected,///< selected data of pick_to_select_buf.
	MGPosition& point,	///<point data will be output.
	const MGObject*& obj,///<point's object will be rturned.
	MGPosition& t	///<When obj is an MGCurve and
		///<this snap kind is nearpos, endpos, knotpos, or centerpos,
		///<the point's parameter value of the curve be returned: t.sdim()=1.
		///<When obj is an MGFSurface and this snap kind is centerpos, or vertexpos
		///<the point's parameter value(u,v) be returned:t.sdim()=2.
)const;

///Get and set the snap_kind.
snap_kind get_snap_kind()const{return m_snap_kind;};
void set_snap_kind(snap_kind kind){m_snap_kind=kind;};

///Get the object of the position posID of m_positions.
const MGObject* object(int posID)const;

///The container m_positions' access functions.
const MGPosition& back()const{return m_positions.back();};
MGPosition& back(){return m_positions.back();};
iterator begin(){return m_positions.begin();};
const_iterator begin()const{return m_positions.begin();};
void clear();
bool empty()const{return m_positions.empty();};
iterator end(){return m_positions.end();};
const_iterator end()const{return m_positions.end();};
const MGPosition& front()const{return m_positions.front();};
MGPosition& front(){return m_positions.front();};
const MGPosition& operator[](int i)const{return m_positions[i];};
MGPosition& operator[](int i){return m_positions[i];};
void pop_back(){m_positions.pop_back();};
void push_back(const MGPosition& pos){m_positions.push_back(pos);};
size_t size(){return m_positions.size();};

///Get the point data array.
const std::vector<MGPosition>& points()const{return m_positions;};
std::vector<MGPosition>& points(){return m_positions;};

private:
std::vector<obj_num> m_obj_nums;///<obj_num includes how many data are included in
	///<m_positions for an MGObject.
	///<m_positions.size()=sum of(m_obj_nums[i].second) for i=0,...,m_obj_nums.size()-1.
std::vector<MGPosition> m_positions;///<All the position data will be stored.
snap_kind m_snap_kind;


};

/** @} */ // end of GLAttrib group
#endif
