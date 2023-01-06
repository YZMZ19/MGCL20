/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
//////////////////////////////////////////////////////////////////////
// MGOfstream.cpp: MGOfstream クラスのインプリメンテーション
//////////////////////////////////////////////////////////////////////
#include "StdAfx.h"
#include "mg/Ofstream.h"
#include "mg/Object.h"
#include "mg/Group.h"
#include "topo/BCell.h"
#include "topo/PCell.h"
#include "mgGL/VBO.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using std::ofstream;
using std::ios;

//////////////////////////////////////////////////////////////////////
// Constructor and destructor.
//////////////////////////////////////////////////////////////////////

MGOfstream::MGOfstream()
:m_file(0),m_position(0){}

MGOfstream::MGOfstream(const TCHAR *file)
:m_file(0),m_position(0){
	open(file);
}

MGOfstream::~MGOfstream(){
	close();
}

///////////////Operator overload//////////////////

//Write out an object to ofs. Written objects by this function are able
//to read using MGIfstream.
MGOfstream& MGOfstream::operator<< (const MGGel& gel){
	assert(abstractGelId(gel.identify_type())!=-1);

	long tid = gel.identify_type();
	if(tid==MGPLANEIMAGE_TID)
		return *this;

	MGOfstream& ostrm = *this;
	ostrm<<0xffffffffL;//Flag that indicates an object is followed.
	ostrm<<tid;
	insert(&gel, tellp());
	gel.WriteMembers(ostrm);
	return ostrm;
}

//Write out the pointer into ostrm.
//This is an internal program. Ordinary users should not use this function.
//operator<< should be used instead.
void MGOfstream::WritePointer(
	const MGGel* gel//input pointer to store into stream.
){
	MGOfstream&  ostrm = *this;//target stream.
	if(gel){
		//自分自身が既にmapに登録されているかどうか調べる
		long pid = ostrm.find(gel);
		if(pid){//When gel is found in the map.
			ostrm<<pid;
		} else{
			ostrm<<0xffffffffL;//Flag that indicates an object is followed.
			long tid = gel->identify_type();//input type id of pointer class.
			ostrm<<tid;
			pid = ostrm.tellp();  // PIDを決める(streamのPosition)
			ostrm.insert(gel, pid);	 //mapに自分自身のPIDを登録
			gel->WriteMembers(ostrm);//データメンバの書き出し
		}
	}else{//When null pointer.
		ostrm<<0x00000000L;
	}
}

///////////////Member function//////////////////

void MGOfstream::close(){
	if(m_file) delete m_file;
	m_file=0;
}

//Function's return value is:
//=0: open succeeded.
//=1: file not found, or could not be opened.
//=2: file found, but, the format is not MGCL format.
int MGOfstream::open(const TCHAR* file){
	int error=0;
	std::ofstream* file_new=0;
	file_new = new ofstream(file, ios::out | ios::binary);

	if(file_new && file_new->is_open()){
		if(m_file){
			m_file->close();
			delete m_file;
		}
		m_file=file_new;
		write(MGCL::Version(),8);
		write(MGCL::File_validity(),24);
		mapClear();
	}else{
		delete file_new; file_new=0;
		error=1;
	}
	return error;
}

//Write out n bytes date in th buffer ps.
//This read data is row data, and the sequence will not be changed
//like read nByte.
void MGOfstream::write(const void* ps, int n){
	m_file->write((char*)ps,n);
	m_position+=n;
}
