/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#pragma once


//Define MGAppearance Class.
class MGOfstream;
class MGIfstream;
class MGGLAttrib;
class MGLight;
#include <iosfwd>
#include "mg/MGCL.h"
#include "mg/Group.h"
#include "mgGL/color.h"

#ifndef _CONSOLE
#include "mgGL/RenderAttr.h"
#endif //_CONSOLE


/** @addtogroup DisplayHandling
 *  @{
 */

///A container class to hold MGGLAttrib objects.

///MGAppearance acts just like as std::auto_ptr.
///That is, MGAppearance holds newed object pointers of MGGLAttrib,
///and when copy constructor or assignment operator is invoked,
///the pointer ownership is transfered to the new MGAppearance object.
///A list of newed MGGLAttrib object pointer will be stored in parent MGGroup.
///No two same leaf type MGGLAttrib objects are included in this list.
class MG_DLL_DECLR MGAppearance : public MGAttrib {

public:
	typedef MGGroup::iterator iterator;
	typedef MGGroup::const_iterator const_iterator;

protected:
	MGGroup m_glattribs;///<Attribute elements of this appearance.
	bool m_no_display;///<True if not to display, false if to display.

public:

	////////Special member functions/////////
	MGAppearance() :m_no_display(false) { ; };
	~MGAppearance() = default;
	MGAppearance(const MGAppearance& rhs) = delete;
	MGAppearance& operator=(const MGAppearance& rhs) = delete;
	MGAppearance(MGAppearance&& rhs) = default;
	MGAppearance& operator=(MGAppearance&& rhs);

	///comparison
	bool operator<(const MGAppearance& gel2)const;
	bool operator<(const MGGel& gel2)const;

	///To print out the contents.
	std::ostream& toString(
		std::ostream& ostrm
	)const override;

	MGAttrib* back() { return static_cast<MGAttrib*>(m_glattribs.back().get()); };
	const MGAttrib* back()const { return static_cast<const MGAttrib*>(m_glattribs.back().get()); };
	iterator begin() { return m_glattribs.begin(); };
	const_iterator begin()const { return m_glattribs.begin(); };
	void clear() { m_glattribs.clear(); };
	bool empty()const { return m_glattribs.empty(); };
	iterator end() { return m_glattribs.end(); };
	const_iterator end()const { return m_glattribs.end(); };
	MGAttrib* front() { return static_cast<MGAttrib*>(m_glattribs.front().get()); };
	const MGAttrib* front()const { return static_cast<const MGAttrib*>(m_glattribs.front().get()); };
	void pop_back() { m_glattribs.pop_back(); };
	void pop_front() { m_glattribs.pop_front(); };
	size_t size()const { return m_glattribs.size(); };

	///Add a light. light must be a newed object, and the ownership will be
	///transfered to this object.
	///Function's return value is the number of lights after added.
	int add_light(MGLight* light);

	///Test if this MGAppearance can be removed or not.
	bool can_be_removed()const;

	///Generate copied gel of this gel.
	///Returned is a newed object. User must delete the object.
	MGAppearance* clone()const;

	///Obtain display list name of the curren rendering context(MGOpenGLView).
	///0(null) means the current MGOpenGLView=null or MGAppearance..
	virtual mgVBO* dlist_name()const { return 0; };

	///Judge if the display list for vmode is made or not.
	bool displayList_is_made(MGCL::VIEWMODE vmode)const { return true; };

	///Erase the specified attribute.
	///Function's return value is the iterator after the erased data.
	iterator erase(iterator i) { return m_glattribs.erase(i); };

	///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
	int get_draw_attrib_mask()const;
	int get_render_attrib_mask()const;

	/// Return This object's typeID
	long identify_type() const { return MGAPPEARANCE_TID; };

	///Get manifold dimension.
	///MGGroup returns right one, MGGroup return 2, and others return -1.
	int manifold_dimension() const { return -1; };

	///Test this is no_display MGAppearance.
	bool no_display()const { return m_no_display; };

	///Release the attribute of specified type.
	///Function's return value is the MGGLAttrib* that is released.
	std::unique_ptr<MGGLAttrib> release_attrib(long tid);

	///Search the same type MGGLAttrib leaf class object in this list.
	///Function's return value is the iterator found.
	///If not found, end() will be returned.
	iterator search(const MGGLAttrib* atr);
	iterator search_by_id(MGGEL_TID tid);
	const_iterator search(const MGGLAttrib* atr)const;
	const_iterator search_by_id(MGGEL_TID tid)const;

	///Set the attribute in this list. attr must be a newed object, and the
	///ownership will be transfered to this MGAppearance.
	void set_attrib(MGGLAttrib* attr);
	void set_attrib(UniqueGLAttribVec& attrs);

	///Set the attribute in this list. attr must be a newed object, and the
	///ownership will be transfered to this MGAppearance.
	///When the appearance held an attribute, the old one will be returned
	///as the function's return value. Users must delete it.
	MGGLAttrib* set_attrib_with_old(MGGLAttrib* attr);

	///Set color data.
	void set_color(const MGColor& color);
	void set_color(const float color[4]);
	void set_color(float red, float green, float blue, float alpha = 1.);

	///Set display/no display.
	void set_display() { m_no_display = false; };
	void set_no_display() { m_no_display = true; }

	///Set the material.

	///Set Line width.
	void setLineWidth(float width);

	///Get the class name.
	virtual std::string whoami()const { return "Appearance"; };

	///メンバデータを読み出す関数
	void ReadMembers(MGIfstream& buf);

	///メンバデータを書き込む関数
	void WriteMembers(MGOfstream& buf) const;

	///draw GLAttributes process.
	void drawAttrib(
		mgVBO& vbo,///<The target graphic object.
		bool no_color = false	///<if true, color attribute will be neglected.
	)const override;

	///render GLAttributes process.
	void render(mgVBO& vbo)const;

#ifndef _CONSOLE

	///When rs=FRONT_AND_BACK and different material for the back side
	///is used, set_back_material must be invoked after invoking set_material.
	///Else the same material will be appllied for the both sides.
	void set_material(
		MGRenderAttr::RENDERSIDE rs,
		const float ambient[3],
		const float diffuse[3],
		const float specular[3],
		const float emission[3],
		float shininess = 0.,
		float transparency = 0.
	);

	///Set the back side material.

	///Invoking set_back_material means two sided material
	///and setting different material to the back side.
	///Before use of set_back_material, set_material must be invoked first.
	///set_back_material will set two sided material.
	void set_back_material(
		const float ambient[3],
		const float diffuse[3],
		const float specular[3],
		const float emission[3],
		float shininess = 0.,
		float transparency = 0.
	);

#endif //_CONSOLE

	///Set Line stipple.

	///When factor=0 is input, line pattern is disabled. 実線となる
	///When factor<0, the stipple attribute is undefined. This means the attribute
	///is defined by the environment.
	///When factor<=0, pattern is unnecessary.
	void setLineStipple(short int factor, unsigned short pattern);

private:
	iterator insert(iterator i, MGAttrib* atr) {
		return m_glattribs.insert(i, std::unique_ptr<MGGel>(atr));
	};
	void push_back(MGAttrib* atr) { m_glattribs.push_back(std::unique_ptr<MGGel>(atr)); };
	void push_front(MGAttrib* atr) { m_glattribs.push_front(std::unique_ptr<MGGel>(atr)); };
};

/** @} */ // end of DisplayHandling group
