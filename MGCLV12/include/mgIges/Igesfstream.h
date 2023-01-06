/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#if !defined( __IGESFSTREAM_H__)
#define __IGESFSTREAM_H__

#include "mg/MGCL.h"
#include "mgIges/Iges.h"
#include "mgIges/IgesGSec.h"
#include "mgIges/IgesDirectoryEntry.h"
#include "mgIges/IgesParamLine.h"
class MGBox;

/** @addtogroup FileInputOutput
 *  @{
 */

///MGIgesFstream is a super class for MGIgesIfstream and MGIgesOfstream.

///MGIgesFstream holds the data for IGES, and provides a common functions
///to write and read IGES data.
class MG_DLL_DECLR MGIgesFstream{
public:
	using UniqueDE = std::unique_ptr<MGIgesDirectoryEntry>;
	using UniqueDEVec = std::vector<UniqueDE>;
	
////////Special member functions/////////
	MGIgesFstream(){ ; };
	virtual ~MGIgesFstream(){ ; };//Destructor.
	MGIgesFstream(const MGIgesFstream& rhs)=delete;//Copy constructor.
	MGIgesFstream& operator=(const MGIgesFstream& gel2)=delete;//Copy assignment.
	MGIgesFstream(MGIgesFstream&& rhs)=default;//Move constructor.
	MGIgesFstream& operator=(MGIgesFstream&& rhs)=default;//Move assignment.

	///Initialize all the member data to the state of no_value_holding.
	virtual void initialize(const TCHAR* filename=0);

///Function's return value is the directory entry pointer pushed back.
	int push_back_DE(MGIgesDirectoryEntry* de);

	///Return directory entry point of DEid.
	UniqueDE& directoryEntry(int DEid){
		return m_DirectoryEntries[DEid];
	};
	const UniqueDE& directoryEntry(int DEid)const{
		return m_DirectoryEntries[DEid];
	};

	void clearStartSection(){m_StartSection=std::string();};
	void clearGSection(){m_GSection=MGIgesGSec();};
	void clearDirectoryEntries(){m_DirectoryEntries.clear();};
	void clear(){clearStartSection(); clearDirectoryEntries(); clearDirectoryEntries();};
	void set_GSec_max_coordinate_value(const MGBox* bx=0);
	void set_initial_StartSection();

	const MGIgesGSec& GSection()const{return m_GSection;};
	MGIgesGSec& GSection(){return m_GSection;};

	///get the output line number of Start Section.
	int get_line_number_of_SS()const{return (int)(m_StartSection.size()/72+1);};

	///get the output line number of Global Sections.
	int get_line_number_of_GS()const{return m_nlineGSec;};

	///get the output line number of Directory Entries.
	int get_line_number_of_DE()const{return ((int)(m_DirectoryEntries.size()-1)*2);};

protected:

// Data members.

	std::string m_StartSection;///<Start section string data.
	MGIgesGSec m_GSection;///<Global section data.
	int m_nlineGSec;	///<Number of the lines of Global section.
	UniqueDEVec m_DirectoryEntries;///<Directry entry data vector.
				///<One pair of directory entry lines are stored in m_DirectryEntry[i].
				///<m_DirectoryEntries[0] is a dummy entry and has no meaning.
};

/** @} */ // end of FileInputOutput group
#endif // __IGESFSTREAM_H__
