#include "mg/MGCL.h"
#include "mg/Position.h"
#include "mg/LBRep.h"

///Bezier curveをつまみ変形するためのclass
class MG_DLL_DECLR MGBezierMove{
public:

	MGBezierMove(
		const MGLBRep& bezier,
		bool Divide,
		double t_pivot
	);

	void movePoint(const MGPosition& ToPoint, MGLBRep& mevedBezier);

private:
	const MGLBRep& m_originalBezier;
	double m_pivot;//移動させる場所を示すm_originalBezierのパラメータ値
	MGPosition m_pivotP;
	bool m_closed;
	int m_index;
	bool m_divide;//折れとしてよいかどうか, =true:折れとしてよい。
	bool m_bIsParallel_s, m_bIsParallel_e;//それぞれ始点、終点側が折れているかどうか。
	
	//指定された点の前後の方向線が平行か調べる。
	//return : 平行ではない･前後に点が無い時はfalse, 平行な時はtrue
	bool isParallel(
		const MGBPointSeq& bp//=調べる点列
		, int index//=調べる点のインデックス(0,1,2,･･･,n-1)
	);
};
