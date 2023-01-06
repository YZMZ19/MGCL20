/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
//////////////////////////////////////////////////////////////////////
// MGIfstream.cpp: MGIfstream クラスのインプリメンテーション
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "mg/Tolerance.h"
#include "mg/Ifstream.h"
#include "mg/Gel.h"
#include "mg/Group.h"
#include "mgGL/Context.h"
#include "topo/BCell.h"
#include "topo/PCell.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using std::ifstream;
using std::ios;

//////////////////////////////////////////////////////////////////////
// Constructor and destructor.
//////////////////////////////////////////////////////////////////////

MGIfstream::MGIfstream()
:m_file(NULL),m_position(0){;}

MGIfstream::MGIfstream(const TCHAR *file)
:m_file(0),m_position(0){
	open(file);
}

MGIfstream::~MGIfstream(){
	close();
}

///////////////Operator overload//////////////////

///Read in an object into gel. Returned gel is:
///(1) null if end of file or illegal file.
///(2) non null newed object pointer. User must delete it.
///If you know gel's concrete class, use its function ReadMembers()(which
///is written out by WriteMembers()).
MGIfstream& MGIfstream::operator>>(MGGel*& gel){
	gel = nullptr;
	long gelHeader=0;
	(*this)>>gelHeader; assert(gelHeader == 0xffffffffL);
	//if(gelHeader != 0xffffffffL) return *this;

	long tid; (*this) >> tid;
	gel=MGNullGel(tid);	assert(gel);
	long pid = tellg();
	insert(pid, gel);	//mapに自分自身のPIDを登録
	gel->ReadMembers(*this);	
	return *this;
}

//Read in header part and make MGGel*.
//The header part has 3 kinds:
//(1) 0x00000000L: This indicates nullptr.
//(2) 0xffffffffL: This indicates type id and member data follow.
//(3) otherwise:   This indicates this is already read object pointer(id or name of the object)
//                 and the id and object pointer is already registerd in m_mapPidGel.
//In cases (1) and (3), read members is unnecessary.
MGGel* MGIfstream::ReadPointer(){
	MGIfstream& istrm = *this;
	long header;
	istrm>>header;
	MGGel* gel = nullptr;
	if(header==0xffffffffL){
		long typeId;
		istrm>>typeId;
		gel = MGNullGel(typeId);
		assert(gel);
		istrm.insert(istrm.tellg(), gel);//mapに自分自身のPIDを登録
		gel->ReadMembers(istrm);
	}else if(header){
		//In this case, the header is pid, an object pointer.
		gel = istrm.findGel(header);
		assert(gel);//Assert that object found.
	}
	return gel;//means a null object pointer was stored.
}

//Check if this is the right MGCL file.
//Function's return value is:
//true: if this is the right MGCL file.
//false: if this is not the right MGCL file.
bool MGIfstream::MGCLHeader(std::ifstream& ar){
	ar.read(m_version,8);
	m_position=8;
	const char* large_version="MGCL1200";
	int comp1=strncmp(large_version, m_version, 4);
	if(comp1!=0)
		return false;
	int comp2=strncmp(large_version+4, m_version+4, 4);
	if(comp2>0)
		return false;

	char fielValid[24];
	ar.read(fielValid,24);
	m_position+=24;
	int comp3=strncmp(MGCL::File_validity(), fielValid, 24);
	if(comp3!=0)
		return false;
	return true;
}

void MGIfstream::close(){
	if(m_file) delete m_file;
	m_file=0;
}

//Function's return value is:
//=0: open succeeded.
//=1: file not found, or could not be opened.
//=2: file found, but, the format is not MGCL format.
int MGIfstream::open(const TCHAR* file){
	int error=0;
	std::ifstream* file_new=0;
	file_new = new ifstream(file, ios::in | ios::binary);

	if(file_new && file_new->is_open()){
		if(MGCLHeader(*file_new)){
			error=0;
			if(m_file) delete m_file;
			m_file=file_new;
			mapClear();
		}else{
			error=2;
			delete file_new;
		}
	}else{
		delete file_new; file_new=0;
		error=1;
	}
	return error;
}

//Read in n bytes date in th buffer ps.
//This read data is row data, and the sequence will not be changed
//like read nByte.
void MGIfstream::read(void* ps, int n){
	m_file->read((char*)ps,n);
	m_position+=n;
}
