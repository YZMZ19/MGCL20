/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

/** @file perimTP.h
 *  Class mgPeimTP is internal for class mgTPmaker.
 *
 *  Declaration for class %mgPerimTP and its family.
 */
#if !defined(__mgPerimTP_H__)
#define __mgPerimTP_H__
#include <iostream>
#include <vector>
#include "mg/LBRep.h"
#include "mg/TPOnFaceBoundary.h"

// Forward declaration.
class mgTPmaker;

/** @addtogroup GEORelated
 *  @{
 */

//! @brief %mgPerimTP represents a perimeter context of quadrangular surface.
//!
//! Declaration of class mgPerimTP.
//! mgPerimTP represents a whole perimeter of MGSBRep to buid(buildByBlendWithTP()).
//! (mgTPOnFaceBoundary constitutes mgPerimTP).
//!
//! The mgTPmaker contains mgPerimTP[4] for each perimeter of MGSBRep to build. 
//! @sa mgTPOnFaceBoundary.
class mgPerimTP{
	friend std::ostream& operator<<(std::ostream&, const mgPerimTP&);

public:

	///Element & vector name definition.
	using MYELM=std::unique_ptr<mgTPOnFaceBoundary>;
	using MYVEC=std::vector<MYELM>;
	
private:	//////// Member data ///////
	mgTPmaker* m_tpmaker;//mgTPmaker's pointer that contais this.
	int m_perimNum;//Perimeter number of this.

	const MGLBRep*  m_perimeter;// A perimeter curve to build MGSBRep.
	MYVEC m_tpVec;/// The material of TP.

public:

	////////Special member functions/////////
	mgPerimTP();///Default Constructor.
	~mgPerimTP()=default;
	mgPerimTP(const mgPerimTP&)=delete;///Copy constructor.
	mgPerimTP& operator= (const mgPerimTP&)=delete;///Copy assignment.
	mgPerimTP(mgPerimTP&&)=default;		///Move constructor.
	mgPerimTP& operator= (mgPerimTP&&)=default;///Move assignment.

	//! Creates an object of class mgPerimTP with a perimeter curve of MGSBRep.
	//! and initializes all its pointer menber objects to null pointers.
	mgPerimTP(
		mgTPmaker* tpmaker,//the parent mgTPmaker that holds this as a perimeter tp data.
		int periNum,//perimeter number of the target MGSBRep to build.
		const MGLBRep* perim //perimeter periNum curve of the target MGSBRep.
	) :m_tpmaker(tpmaker), m_perimNum(periNum),m_perimeter(perim) {;};

// Member functions
public:

	///construct mgPerimTP.
	void build(mgTPmaker* qsmain,int periNum,const MGLBRep* perim);

	//! @brief Computes the TP curve.
	UniqueLBRep create_tp();

	///Test if this has no boudary conditions.
	bool is_0B_perim() const;

	/// Returns the number of objects of class mgTPOnFaceBoundary this object holds.
	size_t number_of_subtp() const;

	/// Outputs the context of an object of class %mgPerimTP.
	void toString(std::ostream& os = std::cout) const;

	//! @brief Creates segments of the tp curve.
	//! @param faces Normal data.
	//! @return Number of objects of class %mgTPOnFaceBoundary this perimeter has.
	size_t buidTpVec(const std::vector<const MGFSurface*>& faces);

	/// Returns the perimeter curve of an object of class %mgPerimTP.
	const MGLBRep& perimeter() const{ return *m_perimeter;}

	//! @brief 
	const MYVEC& subtp_container() const{ return m_tpVec; };

	//! Access to the firstt mgTPOnFaceBoundary object in the perimeter.
	const MYELM& subtp_first() const;
	MYELM& subtp_first();

	//! Access to the last mgTPOnFaceBoundary object in the perimeter.
	const MYELM& subtp_last() const;

	/// Returns the normal vector at the endpoint of this perimeter.
	const MGVector& unit_normal_start() const;
	const MGVector& unit_normal_end() const;

protected:
	//determine if face normal to necessary to negate.
	bool need_negate(double cmnargs, double cmnarge) const;	
};
/** @} */ // end of GEORelated group

#endif // __mgPerimTP_H__
