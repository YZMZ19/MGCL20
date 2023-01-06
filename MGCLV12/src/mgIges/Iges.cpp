/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
//! @file
//!	@brief  Implementaion for class MGIgesIfstream.
//!	@author System fugen

#include "StdAfx.h"
#include "mgiges/Iges.h"
#include "mgGL/Color.h"
#include "mgiges/IgesIfstream.h"
using namespace std;

//Convert an MGIgesFstream's DE pointer to the line number to store in IGES file.
//line_number=2*DEpointer-1;
int MGIges::DEpointer_to_lnumber(int DEpointer){
	if(DEpointer){
		DEpointer=DEpointer<<1;//Multiply by 2 since one DirectoryEntry contains a pair of lines.
		if(DEpointer>0)
			DEpointer-=1;
		else
			DEpointer+=1;
	}
	return DEpointer;
}

//Convert a line number stored in IGES file to the MGIgesFstream's DE pointer.
//DEpointer=(line_number+1)/2
int MGIges::lnumber_to_DEpointer(int line_number){
	if(line_number>0){
		line_number++;
		return (line_number>>1);//Divided by 2 since one DirectoryEntry contains a pair of lines.
	}else if(line_number<0){
		line_number--;
		return line_number/2;
	}else
		return line_number;
}

//convert the line id into int(sequence), inputting one line.
void MGIges::get_ID_sequence(
	const std::string& line,//Input whole line data(1-80)
	char& sectionID_letter,	//section identification letter of the line.
	int& sequence			//ascending sequence number of the line.
){
	sectionID_letter=line[72];
	stringstream seqstream(line.substr(73,7));
	seqstream>>sequence;
}

#define MAX_BUF_LEN 80	//Maximum buffer length able to store in local variable.
//Read in Hollerith_string into strngData.
//Function's return value is
//  true: when value specified, strngData.size() be >0.
//  false:when value not specified, strngData.size() be 0.
bool MGIges::get_Hollerith_string(
	char pDelimeter,	//parameter delimeter
	std::istringstream& istrm,	//Input string stream that contains Hollerith data.
		//The stream pointer will be advanced to the start position of the next item.
	std::string& strngData//output string data that is converted from the istrm's Hollerith data.
){
	char buffer[MAX_BUF_LEN];
	char* bufP=buffer;
	string dummy;
	bool specified;

	int nchar=0;
	istrm>>nchar;
	if(nchar==0){
		istrm.clear();
		specified=false;
	}else{
		getline(istrm,dummy,'H');
		if(nchar>=MAX_BUF_LEN)
			bufP=new char[nchar+1];
		istrm.read(bufP,nchar);
		specified=true;
	}
	bufP[nchar]=0;
	// TODO Unicodeでエラーになるので暫定的にコメントにする
	//strngData.assign(bufP);
	if(nchar>=MAX_BUF_LEN)
		delete[] bufP;
	getline(istrm,dummy,pDelimeter);//(ch==pDelimeter || ch==record_delimeter);
	//int nc2=istrm.gcount();
	return specified;
}

//Read in DE pointer into DEpointer.
//Line number in the istrm is converted to DE pointer.
//Function's return value is
//  true: when value specified.
//  false:when value not specified, DEpointer be 0.
bool MGIges::get_DEpointer(
	char pDelimeter,	//parameter delimeter
	std::istringstream& istrm,	//Input string stream that contains integer data.
		//The stream pointer will be advanced to the start position of the next item.
	int& DEpointer	//output integer data that is converted from the istrm data.
){
	int lnum;
	bool specified=MGCL::get_integer(pDelimeter,istrm,lnum);
	DEpointer=lnumber_to_DEpointer(lnum);
	return specified;
}

//Put integer data into plines, converting into string.
//Except for string data, one integer or double data is output
//into one line, not striding over more than one lines.
void MGIges::put_integer(
	int idata,  //integer to output.
	const MGIgesGSec& gsec,
	std::vector<std::string>& plines, //output plines.
		//lines will be added to input plines.
		//When plines.size()==0, a string whose size()<=63 is appended. 
		//When plines.size()>0, let endpline=plines.back(). When (*endpline)
		//can hold the input data, the data will be appended onto (*endpline).
		//If one string was not enough to store data(greater than 63
		//characters), new lines of string will be created to append after
		//plines.
	int line_len//line length to output, =64(for Parameter data section) or 72.
){
	char pDelimeter=gsec.paramDelimeter(); //parameter delimeter
	stringstream sstr; sstr<<idata;//Convert integer to string.
	const string& istrng=sstr.str();
	if(plines.size()==0){//case of brand new plines.
		plines.push_back(std::move(istrng));
	}else{				//case of plines already has string data.
		string& plend=plines.back();
		int nplendchar=(int)plend.size();
		if(nplendchar>=line_len){
			string newline(1,pDelimeter);
			newline+=istrng;
			plines.push_back(std::move(newline));
		}else{
			plend+=pDelimeter;
			if(nplendchar+istrng.size() < size_t(line_len))
				//1 character space is necessary for the delimeter.
				plend+= istrng;
			else
				plines.push_back(std::move(istrng));
		}
	}
}

//Put real data into plines, converting into string.
//Except for string data, one integer or double data is output
//into one line, not striding over more than one lines.
void MGIges::put_real(
	double rdata,  //double to output.
	const MGIgesGSec& gsec,
	std::vector<std::string>& plines, //output plines.
		//lines will be added to input plines.
		//When plines.size()==0, a new-ed string whose size()<=63
		//will be appended. 
		//When plines.size()>0, let endpline=plines.back(). When (*endpline)
		//can hold the input data, the data will be appended onto (*endpline).
		//If one string was not enough to store data(greater than 63
		//characters), new lines of string will be created to append after
		//plines.
	int line_len//line length to output, =64(for Parameter data section) or 72.
){
	char pDelimeter=gsec.paramDelimeter(); //parameter delimeter
	stringstream sstr;
	sstr<<setprecision(gsec.m_significance_double_precision)
		<<uppercase<<showpoint<<rdata;//Convert double to string.
	const string& rstrng=sstr.str();
	if(plines.size()==0){//case of brand new plines.
		plines.push_back(std::move(rstrng));
	}else{				//case of plines already has string data.
		string& plend=plines.back();
		int nplendchar=(int)plend.size();
		if(nplendchar>=line_len){
			string newline(1,pDelimeter);
			newline+= rstrng;
			plines.push_back(std::move(newline));
		}else{
			plend+=pDelimeter;
			if(nplendchar+rstrng.size() < size_t(line_len))
				//1 character space is necessary for the delimeter.
				plend+= rstrng;
			else
				plines.push_back(std::move(rstrng));
		}
	}
}

//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
void MGIges::put_Hollerith_string(
	const std::string& strngData,//string data to output,
								//will be converted to Hollerith data.
	const MGIgesGSec& gsec,
	std::vector<std::string>& plines, //output plines.
		//lines will be added to input plines.
		//When plines.size()==0, a new-ed string whose size()<=63
		//will be appended. 
		//When plines.size()>0, let endpline=plines.back(). When (*endpline)
		//can hold the input data, the data will be appended onto (*endpline).
		//If one string was not enough to store data(greater than 63
		//characters), new lines of string will be created to append after
		//plines.
	int line_len//line length to output, =64(for Parameter data section) or 72.
){
	// convert to Hollerith data
	int strSize = (int)strngData.length();
	stringstream convStream;
	if(plines.size()){
		convStream << gsec.paramDelimeter();
	}
	if(strSize)
		convStream << strSize << 'H' << strngData;
	// create converted string
	const string& convStr=convStream.str();
	int convSize = (int)convStr.size();

	int pos = 0;	// output start position
	while(pos < convSize){
		string outStr;	// output string
		int appendLength = line_len;
		if(pos == 0){
			// first output
			if(plines.size() == 0){
				if(convSize < line_len){
					appendLength = convSize;
				}
				// case of brand new plines.
				outStr.append(convStr, pos, appendLength);
				plines.push_back(std::move(outStr));
			}else{
				// case of plines already has string data.
				string& plend=plines.back();
				appendLength = (int)(line_len - plend.length());
				if(appendLength == 0){
					// last line is having max length string. so create new line
					if(convSize < line_len){
						appendLength = convSize;
					}else{
						appendLength = line_len;
					}
					outStr.append(convStr, pos, appendLength);
					plines.push_back(std::move(outStr));
				}else{
					// add to last line
					if(convSize < appendLength){
						appendLength = convSize;
					}
					outStr.append(convStr, pos, appendLength);
					plend += std::move(outStr);
				}
			}
		}else{
			// case of plines already has some lines.
			if(convSize - pos < line_len){
				appendLength = convSize - pos;
			}
			outStr.append(convStr, pos, appendLength);
			plines.push_back(std::move(outStr));
		}
		pos += appendLength;
	}
}

void MGIges::put_DEpointer(
	int DEpointer,  //DE pointer to output.
	const MGIgesGSec& gsec,
	std::vector<std::string>& plines //output plines.
		//lines will be added to input plines.
		//When plines.size()==0, a new-ed string whose size()<=63
		//will be appended. 
		//When plines.size()>0, let endpline=plines.back(). When (*endpline)
		//can hold the input data, the data will be appended onto (*endpline).
		//If one string was not enough to store data(greater than 63
		//characters), new lines of string will be created to append after
		//plines.
){
	put_integer(MGIges::DEpointer_to_lnumber(DEpointer),gsec,plines,64);
}

//append record delimeter to plines.
void MGIges::append_record_delimeter(
	char record_del,
	std::vector<std::string>& plines,
	int line_len
){
	string& lastline=plines.back();
	if((int)lastline.length()<line_len)
		lastline+=record_del;
	else{
		plines.push_back(string(1,record_del));
	}
}
