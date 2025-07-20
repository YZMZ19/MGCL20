/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGTolerance_HH_
#define _MGTolerance_HH_
/** @file */
/** @addtogroup BASE
 *  @{
 */
#include "mg/MGCL.h"
///  Forward Declarations
class MGIfstream;
class MGOfstream;
class MGEReal;

//  MGTolerance.h

#define MG_MAX_TOL_STACK_SIZE 10 ///< stack size to push the current tolerance.

//  Defines tolerance.
//

///MGTolerance is a class to hold various tolerance data used in MGCL.

///Tolerances used are follows:
///(Out of these, wc_zero, line_zero are important. If not familiar, forget 
///rc_zero, mach_zero, m_max_knot_ratio, and angle_zero. Defalut values will be OK.)
///(1) wc_zero(world coordinate zero value):
///Two points distance that should be regarded as coincidence in user's world coordinate.
///Typically, intersection of 2 curves are so computed that the i-points be within
///this tolerance.
///(2) line_zero(Line zero in world coordinates):
///Two curves distance that should be regarded as coincidence.
///line_zero>=wc_zero is recommended.
///(3) rc_zero(normalized relative coordinate zero):
///Two points distance that should be regarded as coincidence in normalized (0.,1.) space.
///This is used when the computation is ecexuted other than world coordinates, like
///in parameter space.
///Generally speaking, when the span of the world interedted is L, L*rc_zero=wc_zero.
///(4) angle_zero: Angle that should be regarded as zero angle in radian.
///Generally speaking, mgDBLPAI*rc_zero=angle_zero.
///(5) m_max_knot_ratio:Maximum ratio of neighboring two spans of a knot vector.
///This is used to stabilize the computation to solve linear equations for B-Spline coefficients. 
///(6) mach_zero(machine zero): the smallest floating value to avoid exception of division zero.
class MG_DLL_DECLR MGTolerance {
private:
	MGTolerance();
	MGTolerance(const MGTolerance&)=delete;
	MGTolerance& operator=(const MGTolerance&)=delete;

public:

	///String stream function.
	MG_DLL_DECLR friend std::ostream& operator << (std::ostream&, const MGTolerance&);

	///Test if data is machine zero. Machine Zero version
	MG_DLL_DECLR friend bool MGMZero(double data);

	///Test if two double is equal in world coordinate.
	MG_DLL_DECLR friend bool MGAEqual(double data1, double data2);

	///Test if data is world coordinate zero.
	MG_DLL_DECLR friend bool MGAZero(double data);

	///Test if difference of two data is less than MGTolerance::rc_zero(). 
	MG_DLL_DECLR friend bool MGREqual(double data1, double data2);

	///Test if difference of two data is less than MGTolerance::rc_zero().

	///Test is performed after changing data1 and data2 proportionally for the bigger one to be 1.
	MG_DLL_DECLR friend bool MGREqual2(double data1, double data2);

	///Test if difference of two data is equal.

	///Comparison is:
	///test if abs(data1-data2)/base_length is less than MGTolerance::rc_zero().
	MG_DLL_DECLR friend bool MGREqual_base(double data1, double data2, double base_length);
	MG_DLL_DECLR friend bool MGREqual_baseEReal(MGEReal data1, MGEReal data2, const MGEReal& base_length);

	///Test if data is less or equal to rc_zero().

	///トレランスを考慮して与えられた値が０か調べる1-----Relative Version
	MG_DLL_DECLR friend bool MGRZero(double data);

	///Test if data is less or equal to rc_zero() compared to base_length.

	///Comparison is done after data and base_length are so changed
	///that base_length is 1.
	///If base_length is zero, MGRZero2 returns always false.
	MG_DLL_DECLR friend bool MGRZero2(double data, double base_length);
	MG_DLL_DECLR friend bool MGRZero2(double data, const MGEReal& base_length);

	///Test if angle is right from the value cos(angle).
	MG_DLL_DECLR friend bool MGRight_angle(double cos_data);

	///Test if data is zero angle taking tolerance into account.

	/// If data is small enough, data is almost equal to sin(data).
	/// This fact says we can use sine value instead of radian as data input.
	MG_DLL_DECLR friend bool MGZero_angle(double data);


	//////////// Destructor. 仮想デストラクタ ////////////
	~MGTolerance();


	///Get the tolerance instance. Returned is a singleton.
	static MGTolerance& instance();


	///Return currently used tolerance stack length.
	int stack_length() const { return m_count; };

	///Return maximum tolerance stack size.
	int max_stack_size() const { return MG_MAX_TOL_STACK_SIZE; };

	///Update m_mach_zero. Returned is old mach_zero used so far.
	static double set_mach_zero(double);

	///Update m_wc_zero. Returned is old wc_zero used so far.
	static double set_wc_zero(double);

	///Update m_rc_zero. Returned is old rc_zero used so far.
	static double set_rc_zero(double);

	///Update m_angle_zero. Returned is old angle_zero used so far.
	static double set_angle_zero(double);

	///Update m_line_zero. Returned is old line_zero used so far.
	static double set_line_zero(double);

	///Update m_max_knot_ratio. Returned is old max_knot_ratio used so far.
	static double set_max_knot_ratio(double);

	///Push all the tolerance to stack.
	static void push();

	///Pop all the tolerance from stack.
	static void pop();

	///Return machine zero.
	static double mach_zero() { return instance().m_mach_zero; };

	///Return world zero.
	static double wc_zero() { return instance().m_wc_zero; };

	///Return square of world zero.
	static double wc_zero_sqr() { return instance().m_wc_zero_sqr; };

	///Return relative zero.
	static double rc_zero() { return instance().m_rc_zero; };

	///Return square of relative zero.
	static double rc_zero_sqr() { return instance().m_rc_zero_sqr; };

	///Return angle zero.
	static double angle_zero() { return instance().m_angle_zero; };

	///Return line zero.
	static double line_zero() { return instance().m_line_zero; };

	///Return maximum knot ratio.
	static double max_knot_ratio() { return instance().m_max_knot_ratio; };

	///Dump Functions.
	int dump_size() const;

	///Dump Function.
	int dump(MGOfstream&) const;

	///Restore Function.
	int restore(MGIfstream&);

private:

	//////////// Memeber data. メンバデータ //////////

	double m_mach_zero;	///< 有効値として割り算ができる最小の数値
				///<Machine zero
	double m_wc_zero;	///< 等しいとみなす２点間の距離(ワールド座標）
				///<Two points distance that should be regarded as coincidence
				///<in user's world coordinate.
	double m_wc_zero_sqr;///< m_wc_zero の２乗
				///<Square of m_wc_zero.
	double m_rc_zero;	///< 等しいとみなす２点間の距離(相対座標）
				///<Two points distance that should be regarded as coincidence
				///<in normalized (0.,1.) space.
	double m_rc_zero_sqr;///< m_wc_zero の２乗
				///<Square of m_rc_zero.
	double m_angle_zero;///< ゼロとみなす角度
				///<Angle that should be regarded as zero angle.
	double m_line_zero;				///<Two curves distance that should be regarded as coincidence.
				///<２曲線を等しいとみなす曲線間の距離
	double m_max_knot_ratio;///<Maximum ratio of neighboring two spans of a knot vector.
		///<When (t(i)-t(i-1))/(t(i+1)-t(i)) >= m_max_knot_ratio, t(i+1) and
		///<t(i) are regarded same(multiple) knots.
		///<When (t(i)-t(i-1))/(t(i+1)-t(i)) <= 1/m_max_knot_ratio, t(i) and
		///<t(i-1) are regarded same(multiple) knots.
		///< 隣り合うKnotの比の最大値

///Tolerance Stack.    
	int m_count;				   ///< スタックカウンタ
	double m_mach_zero_stack[MG_MAX_TOL_STACK_SIZE];	///< m_mach_zero スタック
	double m_wc_zero_stack[MG_MAX_TOL_STACK_SIZE]; ///< m_wc_zero スタック
	double m_rc_zero_stack[MG_MAX_TOL_STACK_SIZE]; ///< m_rc_zero スタック
	double m_angle_zero_stack[MG_MAX_TOL_STACK_SIZE]; ///< m_angle_zero スタック
	double m_line_zero_stack[MG_MAX_TOL_STACK_SIZE]; ///< m_line_zero スタック
	double m_max_knot_ratio_stack[MG_MAX_TOL_STACK_SIZE]; ///< m_max_knot_ratio スタック

};

/// Class to set the relative coordinate zero tolerance.

/// On destruction, the saved tolerance is restored if restore was not invoked.
class MG_DLL_DECLR mgTolSetRCZero{
private:
	double m_RCzeroSave;//Save the rczero so far. When minus value is set,
				//restore-action is not performed on dtor.
public:
	mgTolSetRCZero(double rczero){m_RCzeroSave=MGTolerance::set_rc_zero(rczero);};
	~mgTolSetRCZero() { restore(); };
	void restore();
	void update(double errorNew);
};

/// Class to set the world coordinate zero tolerance.

/// On destruction, the saved tolerance is restored if restore was not invoked.
class MG_DLL_DECLR mgTolSetWCZero {
private:
	double m_WCzeroSave;//Save the wczero so far. When minus value is set,
				//restore-action is not performed on dtor.
public:
	mgTolSetWCZero(double wczero){m_WCzeroSave=MGTolerance::set_wc_zero(wczero);};
	~mgTolSetWCZero(){restore();};

	///Restore saved tolerance when ctor invoked.
	void restore();

	///Update the tolerance.
	void update(double errorNew);
};

/// Class to set the line zero tolerance.

/// On destruction, the saved tolerance is restored if restore was not invoked.
class MG_DLL_DECLR mgTolSetLineZero {
private:
	double m_LinezeroSave;//Save the line zero so far. When minus value is set,
				//restore-action is not performed on dtor.
public:
	mgTolSetLineZero(double linezero){m_LinezeroSave=MGTolerance::set_line_zero(linezero);};
	~mgTolSetLineZero() { restore(); };

	///Restore saved tolerance when ctor invoked.
	void restore();

	///Update the tolerance.
	void update(double errorNew);
};

/// Class to set the world coordinate and line zero tolerances.

/// On destruction, the save tolerances are restored if restore was not invoked.
class MG_DLL_DECLR mgTolSetWCLineZero {
private:
	double m_WCzeroSave, m_LineZeroSave;
		//Save the wczero and line zero so far. When minus values are set,
		//restore-action is not performed on dtor.
public:
	mgTolSetWCLineZero(double wczero, double linezero);
	~mgTolSetWCLineZero() { restore(); };

	///Restore saved tolerance when ctor invoked.
	void restore();
};

/// Class to set the relative coordinate zero tolerance.

/// On destruction, the saved tolerance is restored if restore was not invoked.
class MG_DLL_DECLR mgTolSetAngleZero {
private:
	double m_angleZeroSave;//Save the rczero so far. When minus value is set,
				//restore-action is not performed on dtor.
public:
	mgTolSetAngleZero(double angError){m_angleZeroSave=MGTolerance::set_angle_zero(angError);};
	~mgTolSetAngleZero() { restore(); };

	///Restore saved tolerance when ctor invoked.
	void restore();

	///Update the tolerance.
	void update(double errorNew);
};


///String stream function.
MG_DLL_DECLR std::ostream& operator << (std::ostream&, const MGTolerance&);

///Test if data is machine zero. Machine Zero version
MG_DLL_DECLR bool MGMZero(double data);

///Test if two double is equal in world coordinate.
MG_DLL_DECLR bool MGAEqual(double data1, double data2);

///Test if data is world coordinate zero.
MG_DLL_DECLR bool MGAZero(double data);

///Test if difference of two data is less than MGTolerance::rc_zero(). 
MG_DLL_DECLR bool MGREqual(double data1, double data2);

///Test if difference of two data is less than MGTolerance::rc_zero().

///Test is performed after changing data1 and data2 proportionally for the bigger one to be 1.
MG_DLL_DECLR bool MGREqual2(double data1, double data2);

///Test if difference of two data is equal.

///Comparison is:
///test if abs(data1-data2)/base_length is less than MGTolerance::rc_zero().
MG_DLL_DECLR bool MGREqual_base(double data1, double data2, double base_length);
MG_DLL_DECLR bool MGREqual_baseEReal(MGEReal data1, MGEReal data2, const MGEReal& base_length);

///Test if data is less or equal to rc_zero().

///トレランスを考慮して与えられた値が０か調べる1-----Relative Version
MG_DLL_DECLR bool MGRZero(double data);

///Test if data is less or equal to rc_zero() compared to base_length.

///Comparison is done after data and base_length are so changed
///that base_length is 1.
///If base_length is zero, MGRZero2 returns always false.
MG_DLL_DECLR bool MGRZero2(double data, double base_length);
MG_DLL_DECLR bool MGRZero2(double data, const MGEReal& base_length);

///Test if angle is right from the value cos(angle).
MG_DLL_DECLR bool MGRight_angle(double cos_data);

///Test if data is zero angle taking tolerance into account.

/// If data is small enough, data is almost equal to sin(data).
/// This fact says we can use sine value instead of radian as data input.
MG_DLL_DECLR bool MGZero_angle(double data);

/** @} */ // end of BASE group
#endif
