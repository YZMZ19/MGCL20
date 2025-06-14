/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
/// MGPickObject.h : MGPickObject クラスの宣言およびインターフェイスの定義をします。
///
///////////////////////////////////////////////////////////////////////////////////////

#pragma once
#include <vector>
#include "mg/MGCL.h"
#include "mg/Position.h"
#include "mg/GelPosition.h"

class MGGroup;
class MGObject;

/** @addtogroup MGObjectRelated
 *  @{
 */

/// MGPickObject is a class to locate where an object is in a group hierarchy. 

///Generally, A group includes other groups, and the included groups
/// include other groups. In this way, the groups make a group hierachy.
/// MGPickObject represents this hierarcy.
/// m_group is the top MGGroup that includes
/// the object m_object if m_Ghierarchy.size()==0. If m_Ghierarchy.size()>0,
/// m_group includes m_Ghierarchy[0].
///	Let n=m_Ghierarchy.size(), then group m_Ghierarchy[i-1] includes
/// m_Ghierarchy[i] for n=0,...,n-2. m_Ghierarchy[n-1] includes m_object;
/// m_object is the leaf MGObject pointer.
/// Although m_Ghierarchy[i] for i=0,...,n-2 are always MGGroup, m_Ghierarchy[n-1] may be
/// MGShell that includes MGFace. In this case, m_object is the MGFace.
class MG_DLL_DECLR MGPickObject:public MGGelPosition{

public:

////////Special member functions/////////
MGPickObject()=default;
virtual ~MGPickObject()=default;
MGPickObject(const MGPickObject&)=default;///Copy constructor.
MGPickObject& operator= (const MGPickObject&)=default;///Copy assignment.
MGPickObject(MGPickObject&&)=default;		///Move constructor.
MGPickObject& operator= (MGPickObject&&)=default;///Move assignment.

///Constructor of no hierarched group(m_Ghierarchy.size()==0).
explicit MGPickObject(MGGroup* group, MGObject* obj=0)
:MGGelPosition(group,obj){;};

///conversion constructor.
MGPickObject(const MGGelPosition& gelp2):MGGelPosition(gelp2){;};

///Comparison operator.
bool operator<(const MGPickObject& po2)const;
bool operator>(const MGPickObject& po2)const{return po2<(*this);};
bool operator<=(const MGPickObject& po2)const{return !((*this)>po2);};
bool operator>=(const MGPickObject& po2)const{return !(po2>(*this));};

//////////////////オペレーション//////////////

///Generate a newed clone object.
virtual MGPickObject* clone()const;

///Highlightthe object using the display list of this object.
virtual void hilight_using_display_list(
	int line_density	///<line density to draw a surface in wire mode.
)const;

///Get the parameter value of the object at the picked position.
MGPosition& parameter(){return m_parameter;};
const MGPosition& parameter()const{return m_parameter;};

///Set the object parameter value.
void set_parameter(const MGPosition& param){m_parameter=param;};

private:

	// parameter value at the picked position.
	MGPosition m_parameter;
		//m_parameter.sdim()=leaf_object()->manifold_dimension().
};

/** @} */ // end of MGObjectRelated group
