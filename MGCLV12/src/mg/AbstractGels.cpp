/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/

#include "StdAfx.h"
#include "mg/AbstractGels.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGAbstractGels Class.


///Construct MGAbstractGels of a MGAbstractGel. This is a conversion contructor.
MGAbstractGels::MGAbstractGels(const MGAbstractGel& agell):m_vec(1,agell){;};

///////////////operator overloaded//////////////

//////////Member Function//////////

// Output virtual function.
std::ostream& operator<<(std::ostream& ostrm, const MGAbstractGels& agells){
	ostrm<<"MGAbstractGels::number of agells="<<agells.size()<<std::endl;
	MGAbstractGels::const_iterator i=agells.begin(), ie=agells.end();	
	for(int j=0; i!=ie; i++, j++){
		ostrm<<"agel-"<<j<<":"<<(*i)<<std::endl;
	}
	return ostrm;
}
