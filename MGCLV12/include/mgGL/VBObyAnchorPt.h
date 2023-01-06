
#if !defined(_VBOBYANCHORPT__INCLUDED_)
#define _VBOBYANCHORPT__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include "mgGL/Color.h"
#include "mgGL/VBO.h"
#include "mgGL/VBOElement.h"

/** @file */
/** @addtogroup DisplayHandling
 *  @{
 */

/////////////////////////////////////////////////////////////////////////////
// VBObyAnchorPt

/// 文字列描画のためのVBOクラス.

/// AnchorPoint(ワールド座標系)からの相対位置(スクリーン座標系)で
/// 文字列の形状を記述しています。
class MG_DLL_DECLR VBObyAnchorPt :public mgVBO{
protected:
	MGPosition m_anchorPointStatic;///The size specified by setStaticAttribAnchorPoint.
		///The anchorPosition of mgVBOLeaf generated is set to this postion.

public:
	///デフォルトコンストラクタ.
	VBObyAnchorPt(mgGLSL::CoordinateType coordinateType);

	///copy constructor.
	VBObyAnchorPt(const VBObyAnchorPt& vbo);

	///Assignment.
	VBObyAnchorPt& operator=(const VBObyAnchorPt& vbo);

	virtual ~VBObyAnchorPt();

	///@cond
	///When this mgVBOPointer, return the vbo pointer referenced.
	// 自分自身がVBOPointerかどうかをきいてる。
	virtual mgVBO* vboPointer(){return NULL;};
	virtual const mgVBOLeaf* leafPointer(){return NULL;};
	///@endcond

	///Set anchor point.
	void setStaticAttribAnchorPoint(const MGPosition& pos);

	///Get anchor point.
	const MGPosition& anchorPointStatic()const;

protected:
	virtual const MGPosition* getAnchorPosition() override;

};


/** @} */ // end of DisplayHandling group
#endif //_VBOBYANCHORPT__INCLUDED_