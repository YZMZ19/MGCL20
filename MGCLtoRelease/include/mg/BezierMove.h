#include "mg/MGCL.h"
#include "mg/Position.h"
#include "mg/LBRep.h"

///Bezier curve���܂ݕό`���邽�߂�class
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
	double m_pivot;//�ړ�������ꏊ������m_originalBezier�̃p�����[�^�l
	MGPosition m_pivotP;
	bool m_closed;
	int m_index;
	bool m_divide;//�܂�Ƃ��Ă悢���ǂ���, =true:�܂�Ƃ��Ă悢�B
	bool m_bIsParallel_s, m_bIsParallel_e;//���ꂼ��n�_�A�I�_�����܂�Ă��邩�ǂ����B
	
	//�w�肳�ꂽ�_�̑O��̕����������s�����ׂ�B
	//return : ���s�ł͂Ȃ���O��ɓ_����������false, ���s�Ȏ���true
	bool isParallel(
		const MGBPointSeq& bp//=���ׂ�_��
		, int index//=���ׂ�_�̃C���f�b�N�X(0,1,2,���,n-1)
	);
};
