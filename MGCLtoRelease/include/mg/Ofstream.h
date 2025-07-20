/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#pragma once

#include "mg/MGCL.h"

#include <stdio.h>
#include <fstream>
#include <map>

//forward declerations
class MGObject;
class MGGel;
class MGPCell;
class MGBCell;
class MGGroup;

/** @addtogroup FileInputOutput
 *  @{
 */

///MGOfstream is a class to serialize all of the subclasses of MGGel.

///Generally MGGroup is a class to hold all of the subclasses of MGGel,
///and the followig is the standard serialize sequence:
///(1) Generate MGGroup that includes necessary MGGel.
///(2) construct MGOfstream object like: MGOfstream ofstrm(file_name);
///(3) then, invoke operator<< like: ofstrm<<MGGroup;
class MG_DLL_DECLR MGOfstream{

typedef std::map<const MGGel*, long> MGOutPtrMap;
typedef std::map<const MGGel*, long>::iterator mapitr;

private:
	// Object's pointer map.
	MGOutPtrMap m_mapGelPid;

	int m_position;///<Output buffer's current position in MGOfstream by byte.

public:

/// Target input data stream. m_file is a public member.
/// All of the CFile's member functions can be used.
	std::ofstream* m_file;
////////////////////////////////////////////////

////////////constructor & destructor////////////

///Default constructor.
MGOfstream();

///Ordinal constructor. File name file is used to open the file.
MGOfstream(const TCHAR* file);

///Destructor.
virtual ~MGOfstream();

///////////////operator overload////////////////

///Write out an object to ofs. Written objects by this fucction are able
///to read by the following global function:
///  MGIfstream& MGIfstream::operator>>(MGGel*& gel);
MGOfstream& operator<< (const MGGel& gel);

	/// 基本型のファイル出力関数
MGOfstream& operator<<(char ch){write1Byte(&ch);return *this;};
MGOfstream& operator<<(unsigned char uch){write1Byte(&uch);return *this;};
MGOfstream& operator<<(signed char sch){write1Byte(&sch);return *this;};
MGOfstream& operator<<(short s){write2Byte(&s);return *this;};
MGOfstream& operator<<(unsigned short us){write2Byte(&us);return *this;};
MGOfstream& operator<<(int n){write4Byte(&n);return *this;};
MGOfstream& operator<<(unsigned int un){write4Byte(&un);return *this;};
MGOfstream& operator<<(long l){write4Byte(&l);return *this;};
MGOfstream& operator<<(unsigned long ul){write4Byte(&ul);return *this;};
MGOfstream& operator<<(float f){write4Byte(&f);return *this;};
MGOfstream& operator<<(double d){write8Byte(&d);return *this;};

////////////////member function////////////////

///Close the file. This can be used even open() was not used.
///Users need not use this close() if need not specify the file close
///before the destruction of the MGOfstream.
void close();

///Find the input prt's map address.
///If found, <position, 'true'> will be returned.
///If not found, <0, 'false'> will be returned.
///Here position means std::stremap of the m_file file where the ptr' is stored.
long find(const MGGel* ptr);

///Clear the map area m_mapGelPid.
void mapClear();	///Clear map

///Insert the ptr into the map. Function's return value is:
///True: if ptr did not exist in the map and insertion succeeded.
///False: if ptr did exist in the map and insertion failed.
bool insert(const MGGel* ptr, long pid);

/// fileのオープン
///Open the file. This is valid only when default constructor MGOfstream()
///is used.
///=0: open succeeded.
///=1: file not found.
int open(const TCHAR* file);

///Tell the current optput position of m_file(same as m_file.tellp).
int tellp();

///n bytes character data version.
MGOfstream& writenChar(const char* ps, int n);

/// Pointer base のオブジェクトファイル出力関数
/// 戻り値はオブジェクトのPID(識別ID)。
///This is an internal program. Ordinary users should not use this function.
///operator<< should be used instead.
void WritePointer(const MGGel* obj);

private:

///UNIXスタイルのbinary write関数

///Write out n bytes date to the buffer ps.
///This write data is row data, and the sequence will not be changed
///like write nByte.
void write(const void* ps, int n);

///1Byte版
MGOfstream& write1Byte(const void* ps2);

///2Byte版
MGOfstream& write2Byte(const void* ps2);

///4Byte版
MGOfstream& write4Byte(const void* ps4);

///8Byte版
MGOfstream& write8Byte(const void* ps8);


friend class MGGroup;
};

/** @} */ // end of FileInputOutput group
