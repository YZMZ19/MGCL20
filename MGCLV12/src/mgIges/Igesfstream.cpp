/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno             */
/* All rights reserved.                                             */
/********************************************************************/

//! @file
//!	@brief  Declaration for class MGIgesFstream.
//!	@author System fugen

#include "StdAfx.h"
#include "mg/Box.h"
#include "mgIges/Igesfstream.h"

//Initialize all the member data to the state of no_value_holding.
void MGIgesFstream::initialize(const TCHAR* filename){
	m_StartSection=std::string();
	m_GSection=MGIgesGSec(filename);
	m_DirectoryEntries.clear();
	MGIgesDirectoryEntry* de=new MGIgesDirectoryEntry;
	m_DirectoryEntries.emplace_back(de);//Set the dummy record.
}

//Function's return value is the directory entry pointer pushed back.
int MGIgesFstream::push_back_DE(MGIgesDirectoryEntry* de){
	int deNum=(int)m_DirectoryEntries.size();
	m_DirectoryEntries.emplace_back(de);
	return deNum;
}

void MGIgesFstream::set_GSec_max_coordinate_value(const MGBox* bx){
	double maxCValue=10000.;
	if(bx){
		const MGBox& box=*bx;
		int sd=box.sdim();
		for(int i=0; i<sd; i++){
			const MGInterval& rngi=box[i];
			double maxi=rngi.high_point();
			if(i){
				if(maxi>maxCValue)
					maxCValue=maxi;
			}else
				maxCValue=maxi;
		}
	}
	m_GSection.m_max_coordinate_value=maxCValue;
}
