/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGInterval_HH_
#define _MGInterval_HH_
#include "mg/EReal.h"
/** @file */
/** @addtogroup BASE
 *  @{
 */

//  MGInterval.h
//  header for class MGInterval
// Forward Declaration
class MGIfstream;
class MGOfstream;

/// Interval of 1 dimension, i.e. MGInterval is a real line.

///  ��������̋�Ԃ�\���B�Q�� MGEReal �ŕ\�������B
class MG_DLL_DECLR MGInterval{

	///  Two points on real line. ��������̂Q�_ 
	MGEReal m_low;	///<smaller one.
	MGEReal m_high;	///<larger one.

public:

///+operator.
MG_DLL_DECLR friend MGInterval operator+ (double, const MGInterval& );

///Scalar multiplication to the interval.
///The mid_point of the interval does not change. The width of the interval
///is widened(narrowed) by the multiplication of the value a.
MG_DLL_DECLR friend MGInterval operator* (double a, const MGInterval& );

///String stream function
MG_DLL_DECLR friend std::ostream& operator<< (std::ostream&, const MGInterval&);

//////////////Constructor////////////

/// void Constructor
MGInterval();

///Construct interval from interval type and a point.
///When the type is MGINTERVAL_FINITE, start and end is the point.
///MGINTERVAL_TYPE:MGINTERVAL_EMPTY, MGINTERVAL_FINITE,
///MGINTERVAL_FINITE_ABOVE, MGINTERVAL_FINITE_BELOW, or MGINTERVAL_INFINITE.
///MGINTERVAL_FINITE_ABOVE is infinite below and finite above.
///MGINTERVAL_FINITE_BELOW is infinite above and finite below.
/// �C���^�[�o���^�C�v�ƈ�_���w�肵�ăC���^�[�o���𐶐��B
///     MGINTERVAL_FINITE �w�莞 �w�肳�ꂽ��_�̓_�͈͂̃C���^�[�o��
///     MGINTERVAL_FINITE_ABOVE, MGINTERVAL_FINITE_BELOW �w�莞 
///     �L���̈�_���w��
///     MGINTERVAL_INFINITE �w�莞 double �̒l�̎w��͕K�v�Ȃ�
explicit MGInterval(MGINTERVAL_TYPE, double=0.0);

///Construct interval from two MGEReal points.
///t1 is lower limit, and t2 is upper limit.
///When t1> t2, the interval is empry.
/// �Q��MGEReal�̒l����C���^�[�o���𐶐��B
///***** This is the fundamental constructor.
MGInterval(const MGEReal& t1, const MGEReal& t2); 

/// Construct an interval which contains both input interval and a value.
MGInterval(const MGInterval&, const MGEReal&); 

///Construct a interval of one point.
explicit MGInterval(double t);

////////////Operator overload////////////

///Return low or high point.
///0<=i<=1, and for i=0, return low point, for i=1, return high point.
const MGEReal& operator() (int i)const;
MGEReal& operator() (int i);
const MGEReal& operator[] (int i)const;
MGEReal& operator[] (int i);

///Addition of two intervals. Generated interval is the minimum one
///that includes both intervals. The same as operator| .
///  ���g�̃C���^�[�o���Ɨ^����ꂽ�C���^�[�o���̉��Z���s���I�u�W�F�N�g
///  �𐶐�����B
MGInterval operator+ (const MGInterval& ) const;

///Translation of the interval.
///  ���g�̃C���^�[�o���� double �̉��Z���s���u�W�F�N�g��������B
MGInterval operator+ (double) const;

///Addition of two intervals. Updated interval is the minimum one
///that includes both original intervals. The same as operator|=. 
///  ���g�̃C���^�[�o���ɉ��Z�����g�̃C���^�[�o���Ƃ���B
MGInterval& operator+= (const MGInterval& intrvl2){return operator |=(intrvl2);};

///Translation of the interval.
MGInterval& operator+= (double);

///Unary minus.
///  �O�u�P���}�C�i�X�B�I�u�W�F�N�g�𐶐��B
MGInterval operator- () const;

///Subtraction. Subtraction is defined as:
/// "this" minus the common of the two intervals. 
///  ���g�̃C���^�[�o���Ɨ^����ꂽ�C���^�[�o���̌��Z���s���I�u�W�F�N�g
///  �𐶐�����B
MGInterval operator- (const MGInterval& ) const ;

///Translation of the interval.
///  ���g�̃C���^�[�o���� double �̌��Z���s���u�W�F�N�g��������B
MGInterval operator- (double) const;

///Subtraction. Subtraction is defined as:
/// "this" minus the common of the two intervals. 
///  ���g�̃C���^�[�o���Ɍ��Z�����g�̃C���^�[�o���Ƃ���B
MGInterval& operator-= (const MGInterval& );

///Translation of the interval.
MGInterval& operator-= (double);

///Scalar multiplication to the interval.
///The mid_point of the interval does not change. The width of the interval
///is widened(narrowed) by the multiplication of the value a.
///  �X�J���[�̏�Z���s���I�u�W�F�N�g�𐶐�����B
MGInterval operator* (double a) const;

///Scalar multiplication to the interval.
///  �X�J���[�̏�Z���s�����g�̃C���^�[�o���Ƃ���B
///The mid_point of the interval does not change. The width of the interval
///is widened(narrowed) by the multiplication of the value a.
MGInterval& operator*= (double a);

///Scalar division to the interval.
///  �X�J���[�̏��Z���s���I�u�W�F�N�g�𐶐�����B
///The mid_point of the interval does not change. The width of the interval
///is widened(narrowed) by the division of the value a.
MGInterval operator/ (double a) const; 

///Scalar division to the interval.
///  �X�J���[�̏��Z���s�����g�̃C���^�[�o���Ƃ���B
///The mid_point of the interval does not change. The width of the interval
///is widened(narrowed) by the division of the value a.
MGInterval& operator/= (double a);

///Addition of two intervals. Generated interval is the minimum one
///that includes both intervals. The same as operator+ .
///  ���g�̃C���^�[�o���Ɨ^����ꂽ�C���^�[�o����OR���������C���^�[�o����
///  ��������B
MGInterval operator| (const MGInterval&) const ;

///Addition of two intervals. Updated interval is the minimum one
///that includes both original intervals. The same as operator+=. 
///���g�̃C���^�[�o���Ɨ^����ꂽ�C���^�[�o����OR�������Ď��g�̃C���^�[�o��
///�Ƃ���B
MGInterval& operator|= (const MGInterval&);

///And operation. Generate the interval common to the both intervals.
///���g�̃C���^�[�o���Ɨ^����ꂽ�C���^�[�o���̋��ʕ����̃C���^�[�o����
///��������B
MGInterval operator& (const MGInterval&) const ;

///And operation. Update to the interval common to the both intervals.
///  ���g�̃C���^�[�o���Ɨ^����ꂽ�C���^�[�o�����ʕ��������g�̃C���^�[�o��
///  �Ƃ���B
MGInterval& operator&= (const MGInterval&);

///  Boolean ���Z

///Test if the two intervals have common part.
///  ���g�̃C���^�[�o���Ɨ^����ꂽ�C���^�[�o�������ʕ����������Ă��邩
///  �ǂ�����ԋp����B
bool operator&& (const MGInterval&)const;

///Test if this includes intrval2. 
/// ���g�̃C���^�[�o���� empty �̏ꍇ false(0) ��ԋp���A�^����ꂽ�C���^�[�o��
/// �� empty �̏ꍇ true(1) ��ԋp����B�^����ꂽ�C���^�[�o���̉��[�����g
/// �̃C���^�[�o���̉��[���傫���A�^����ꂽ�C���^�[�o���̏�[�����g��
/// �C���^�[�o���̏�[������������ true(1) ��ԋp����B
/// �܂��A�����̃C���^�[�o���� empty �̎� true(1) ��ԋp����B
bool includes(const MGInterval& intrval2) const;

///Comparison
bool operator== (const MGInterval& i2) const;
std::partial_ordering operator<=> (const MGInterval& i2) const;

/// <summary>
/// equivalent of MGInterval i & MGEReal means i.includes(e).
/// </summary>
/// <param name="i2"></param>
/// <returns></returns>
bool operator== (const MGEReal& i2) const;
bool operator== (double t) const;
auto operator<=> (const MGEReal& i2) const->std::partial_ordering;
auto operator<=> (double t) const->std::partial_ordering;

///Chnage the interval range to [t0, t1] if t0<t1,
///                          to [t1, t0] if t1<t0.
void change_range(
	double t0,	///<new range, start.
	double t1	///<new range, end.
);

///Test if this is empty:�C���^�[�o���� empty ���ǂ������肷��B
//Interval of m_low==m_high is not empty when m_low and m_high are finite.
bool empty() const;

///Expand the interval so that this includes the doube val.
void expand(const MGEReal& val);

///Test if finite.
///  �C���^�[�o�����L�����ǂ������肷��B
bool finite() const { return !infinite(); };

///Test if finite above(infinite below).
///  �C���^�[�o��������L�����ǂ������肷��B
bool finite_above() const;
 
///Test if finite below(infinite above).
///  �C���^�[�o���������L�����ǂ������肷��B
bool finite_below() const;

///Get the high value.
const MGEReal& high() const {return m_high;};
MGEReal& high() {return m_high;};

///Return high point of the interval. When infinite_above(),
///mgInfiniteVal will be returned.�@When empty, return 0.
///����l��ԋp����B
double high_point() const { return m_high.value(); };

///Test if this finite interval includes the value,
///taking the relative error of the interval into the account.
bool includes(double t)const;
bool includes(const MGEReal& t)const;

///Test if infinite.
///  �C���^�[�o�����������ǂ������肷��B
bool infinite() const;

///Test if infinite_above.
bool infinite_above() const;

///Test if infinite_below.
bool infinite_below() const;

///Compute interpolated point of the interval. Given parameter t,
///the interpolated point=(1-t)*low_point()+t*high_point().
///t may be t<0 or t>1. Valid only when finite().
///  �C���^�[�o���̕�Ԃ�ԋp����B�^����ꂽ param �ɑ΂��āA
///  (1-t)*low_point()+t*high_point() ��ԋp����B
///  ( 0 > t, 1 < t �̎��͊O�� )
double interpolate(double t) const;

///Test if this is empty interval.
bool is_null() const{return empty();};

///Compute the length of the interval.
///When empty, length()<=0.0.
///��[�Ɖ��[�̍��ق�ԋp����Bempty �͔񐳒l�i�[���܂��͕��̒l�j
///��ԋp����B
MGEReal length() const;

///Get the low value.
const MGEReal& low() const {return m_low;};
MGEReal& low() {return m_low;};

///Return low point of the interval. When infinite_below(),
/// -mgInfiniteVal will be returned. When empty, return 0.
///  �C���^�[�o�����L���̎��L���ŉ����l��ԋp����B��̎� 0.0 ��
///  �ԋp����B
double low_point() const { return m_low.value(); };


///Return mid point of the interval. Valid only when finite().
///  �C���^�[�o�����L���̎��L���Œ��_��ԋp����B
double mid_point() const;

///Compute relative_error*length();
double relative_error() const;

///Round the input value into this interval's range, that is:
///round_into_interval(t)=(*this)[0] if t<=(*this)[0],
///round_into_interval(t)=(*this)[1] if t>=(*this)[1].
double round_into_interval(double t)const;

///Set the interval length to len so that the midpoint does not change.
void set_length(double len);

///Set this as null(empty) interval
void set_null(double t=0.);
void set_empty(double t=0.){set_null(t);};

///Update high point of the interval. When given point is less than
///low point, the interval becomes empty.
///  ���g�̃C���^�[�o���̏���l��ύX����B�^����ꂽ�l�������l��菬����
///  �ꍇ�A���g�̃C���^�[�o���� empty �ɂȂ�A�㉺���l�̓��ꊷ���͂Ȃ��B
void set_high_point(const MGEReal&);
void set_high(const MGEReal& t){set_high_point(t);};

///Update low point of the interval. When given point is greater than
///high point, the interval becomes empty.
///  ���g�̃C���^�[�o���̉����l��ύX����B�^����ꂽ�l������l���傫��
///  �ꍇ�A���g�̃C���^�[�o���� empty �ɂȂ�A�㉺���l�̓��ꊷ���͂Ȃ��B
void set_low_point(const MGEReal&);
void set_low(const MGEReal& t){set_low_point(t);};

///Return interval type.
///  �C���^�[�o���̃^�C�v��ԋp����B
MGINTERVAL_TYPE type() const;

///Dump Functions.
///Calculate dump size
int dump_size() const;

///Dump Function
int dump(MGOfstream& ) const;

///Restore Function
int restore(MGIfstream& );

private:


};

/** @} */ // end of BASE group
#endif
