/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
//////////////////////////////////////////////////////////////////////
// MGIfstream.cpp: MGIfstream クラスのインプリメンテーション
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "mg/Ifstream.h"
#include "mgGL/VBO.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

//Find the input pid's map address.
//If found, MGGel* will be returned.
//If not found, null pointer(0) will be returned.
MGGel* MGIfstream::findGel(long pid){
	MGGel* ptr = NULL;
	mapitr mitr = m_mapPidGel.find(pid);
	if (mitr != m_mapPidGel.end()){
		ptr = (*mitr).second;
	}
	return ptr;
}

//Insert the ptr into the map. Function's return value is:
//True: if ptr did not exist in the map and insertion succeeded.
//False: if ptr did exist in the map and insertion failed.
bool MGIfstream::insert(long pid, MGGel* ptr){
	std::pair<mapitr, bool> ret = m_mapPidGel.insert(std::make_pair(pid,ptr));
	return ret.second;
}

void MGIfstream::insertSharedBCell(SharedBCell * sbcel){
	MGBCell* bcelP = sbcel->get();
	m_SharedBCellMap.insert(std::make_pair(bcelP, sbcel));
}

SharedBCell * MGIfstream::findSharedBCell(MGBCell * bcel){
	SharedBCell* ptr = nullptr;
	sharedBCellMapItr mitr = m_SharedBCellMap.find(bcel);
	if(mitr != m_SharedBCellMap.end()){
		ptr = (*mitr).second;
	}
	return ptr;
}
//Clear maps
void MGIfstream::mapClear(){
	m_mapPidGel.clear();
	m_SharedBCellMap.clear();
}

//バイト列の並び順をUNIXスタイルに変更して出力する関数群
MGIfstream& MGIfstream::read1Byte(void *ps){
	read(ps,1);
	return *this;
}

//バイト列の並び順をUNIXスタイルに変更して出力する関数群
MGIfstream& MGIfstream::read2Byte(void *ps2){
	char buf[2];
	read(buf,2);
	char* tmp = (char*)ps2;
	tmp[0]=buf[1]; tmp[1]=buf[0];
	return *this;
}

MGIfstream& MGIfstream::read4Byte(void *ps4){
	char buf[4];
	read(buf,4);
	char* tmp = (char*)ps4;
	tmp[0]=buf[3]; tmp[1]=buf[2]; tmp[2]=buf[1]; tmp[3]=buf[0];
	return *this;
}

MGIfstream& MGIfstream::read8Byte(void *ps8){
	char buf[8];
	read(buf,8);
	char* tmp = (char*)ps8;
	tmp[0]=buf[7]; tmp[1]=buf[6]; tmp[2]=buf[5]; tmp[3]=buf[4];
	tmp[4]=buf[3]; tmp[5]=buf[2]; tmp[6]=buf[1]; tmp[7]=buf[0];
	return *this;
}

MGIfstream& MGIfstream::readnChar(char* ps, int n){
	read(ps,n);
	return *this;
}

long MGIfstream::tellg(){
	return m_position;
}

