/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/tolerance.h"
#include "mg/OfStream.h"
#include "mg/IfStream.h"
#include "mg/MGStl.h"
#include "mgGL/Image.h"

// �t�@�C���̎�ނ𒲂ׂ�
// in[in]�F�t�@�C���X�g���[��
// �߂�l�F�o�C�i���t�@�C��(true)�AAscii�t�@�C��(false)
// ���O�����F���łɃI�[�v������Ă���t�@�C���X�g���[����n��
bool IsBinaryFile(
	 std::ifstream& in // [in/out]:���ɃI�[�v�����Ă���t�@�C���X�g���[��
){
	char buff[85];
	in.read(buff, 85);
	size_t readByte = (size_t)in.gcount();
	for(size_t i = 80; i < readByte; i++){
		if(buff[i] == 0){
			// �V�[�N�ʒu���t�@�C���̐擪�ɖ߂�
			in.seekg(0, std::ios::beg);
			return true;
		}
	}
	// �V�[�N�ʒu���t�@�C���̐擪�ɖ߂�
	in.seekg(0, std::ios::beg);
	return false;
}

//constructor from triangle data, index+vertices 
MGStl::MGStl(
	int nTriang, // �O�p�`�̐�
	const int* triang, // ���_���W�̃C���f�b�N�X�̔z��
	const double* verts // ���_�̍��W
){
	// Keep vertex position and index of vertex
	triangleMap VertexMap;
	// index of vertex
	int i0, i1, i2;
	for(int j = 0; j < nTriang; j++){
		int j3 = 3 * j;
		// Get index of cordinate of vertex
		i0 = triang[j3++] - 1;
		i1 = triang[j3++] - 1;
		i2 = triang[j3] - 1;
		// Get cordinate of vertex
		const MGPosition pos1(3, verts+3*i0); m_box.expand(pos1);
		const MGPosition pos2(3,verts+3*i1); m_box.expand(pos2);
		const MGPosition pos3(3,verts+3*i2); m_box.expand(pos3);
		push_back_triangle(pos1, pos2, pos3, VertexMap);
	}
}

///Constructor from triangle data, index+vertices.
///This constructor uses the wc_zero to identify
///different two positions as the same input position.
MGStl::MGStl(
	const std::vector<MGPosition>& vertices,/// trang[id0], [id1], [id2] make
		//a triangle for id0=indices[i*3], id1=indices[i*3+1], id2=indices[i*3+2]
		//for i=0,...,indices.size()/3.
	const std::vector<int>& indices /// ���_��id.
		//(0,1,2),...(i*3, i*3+1, i*3+2) for i=0,...,indices.size()/3.
){
	// Keep vertex position and index of vertex
	triangleMap VertexMap;

	size_t nTriang=indices.size()/3;
	// index of vertex
	size_t i0, i1, i2;
	for(size_t j=0; j<nTriang; j++){
		size_t j3 = j*3;

		// Get index of cordinate of vertex
		i0 = indices[j3++] - 1;
		i1 = indices[j3++] - 1;
		i2 = indices[j3] - 1;
		// Get cordinate of vertex
		const MGPosition& pos1=vertices[i0]; m_box.expand(pos1);
		const MGPosition& pos2=vertices[i1]; m_box.expand(pos2);
		const MGPosition& pos3=vertices[i2]; m_box.expand(pos3);
		push_back_triangle(pos1, pos2, pos3, VertexMap);
	}
}

// �S�Ă̒��_�ƃ{�b�N�X�̍��W�Ɏw�肵��MGVector�̒l�����Z
MGStl& MGStl::operator+=(const MGVector& v){
	size_t nPos = m_vecPos.size();
	for(size_t i = 0; i < nPos; i++){
		m_vecPos[i] += v;
	}
	m_box += v;
	return *this;
}

// �S�Ă̒��_�ƃ{�b�N�X�̍��W����w�肵��MGVector�̒l������
MGStl& MGStl::operator-=(const MGVector& v){
	size_t nPos = m_vecPos.size();
	for(size_t i = 0; i < nPos; i++){
		m_vecPos[i] -= v;
	}
	m_box -= v;
	return *this;
}

// �S�Ă̒��_�ƃ{�b�N�X�̍��W�Ɏw�肵���l���|����
MGStl& MGStl::operator*=(double scale){
	size_t nPos = m_vecPos.size();
	for(size_t i = 0; i < nPos; i++){
		m_vecPos[i] *= scale; 
	}
	m_box *= scale;
	return *this;
}

// �^����ꂽ�ϊ����s��
MGStl& MGStl::operator*=(const MGMatrix& mat){
	size_t nVecPos = m_vecPos.size();
	for(size_t i = 0; i < nVecPos; i++){
		m_vecPos[i] *= mat;
	}
	size_t nVecTriang = m_vecNormlTriang.size();
	for(size_t j = 0; j < nVecTriang; j++){
		m_vecNormlTriang[j] *= mat;
	}
	m_box *= mat;
	return *this;
}

// �^����ꂽ�ϊ��ɂ��g�����X�t�H�[�����s��
MGStl& MGStl::operator*=(const MGTransf& tr){
	size_t nPos = m_vecPos.size();
	for(size_t i = 0; i < nPos; i++)	{
		m_vecPos[i] *= tr;
	}
	size_t nVecTriang = m_vecNormlTriang.size();
	for(size_t j = 0; j < nVecTriang; j++)	{
		m_vecNormlTriang[j] *= tr.affine();
	}
	m_box *= tr;
	return *this;
}

// ���MGStl�I�u�W�F�N�g������������r����
bool MGStl::operator==(const MGStl& stl){
	// �{�b�N�X�g�̔�r
	if(this->m_box != stl.m_box)
		return false;

	// ���W�̔�r
	if(m_vecPos != stl.m_vecPos){
			return false;
	}

	// �ʂ̖@���x�N�g���̔�r
	if(m_vecNormlTriang != stl.m_vecNormlTriang){
		return false;
	}

	if(m_indices != stl.m_indices){
		return false;
	}

	// ��v���Ă���ꍇ��true��ԋp����
	return true;
}

//Construct new object by copying to newed area.
//User must delete this copied object by "delete".
MGStl* MGStl::clone()const{
	return new MGStl(*this);
}

// �w�肵���O�p�`�̒��_���W��������Ă���z��̃C���f�b�N�X���擾����
// ���O����:0����O�p�`�̖���-1�𒴂����͈͂̒l���P�����Ɏw�肵�Ȃ�����
// �������:pos[3]�Ɏw�肵���O�p�`�̊e���_���W�̔z��̓Y�������i�[�����
void MGStl::GetVertIndices(
	 int i,// [in]�F�O�p�`�̃C���f�b�N�X(0 <= i < GetTriangleCount)
	 int pos[3]// [out]�F�w�肵���O�p�`�̊e���_���W�̔z��̓Y����[i, j, k]
)const{
	// �w�肵���O�p�`�̒��_�ԍ��̔z��̃C���f�b�N�X�����߂Ă���
	int i3 = i*3;
	assert(size_t(i3+2) < m_indices.size());

	// ���_���W�̃C���f�b�N�X����
	pos[0] = m_indices[i3++];
	pos[1] = m_indices[i3++];
	pos[2] = m_indices[i3];
}

///Compute box of the geometry.
///Compute the box from the scratch.
void MGStl::compute_box(MGBox& bx) const{
	bx.set_null();
	for(auto& P: m_vecPos)
		bx.expand(P);
}

// �����Ŏw�肵���p�X��STL�t�@�C����ǂݍ��݃����o�ɒl��ݒ肷��
// �܂�MGTolerance::wc_zero�ɒl��ݒ肷��
// �߂�l: =0�t�@�C���̓ǂݍ�������
//		  !=0 �ǂݍ��܂Ȃ������B�܂��͎��s����(std::ifstream�̃G���[�R�[�h�j
// �������:m_vecPos, m_vecNorml, m_indices, m_box�ɐ}�`�̏�񂪊i�[�����
//			MGTolerance::wc_zero�ɒl���ݒ肳���
int MGStl::LoadFile(
	const TCHAR* strFilePath // [in]:�ǂݍ���STL�t�@�C���ւ̃p�X
){
	// �X�g���[�����I�[�v��
	std::ifstream in(strFilePath, std::ios_base::binary);
	if(!in){
		int state = in.rdstate();
		// �I�[�v�������s�����ꍇ�A�G���[�R�[�h��ԋp
		return state;
	}

	Initialize();
	int error; // �t�@�C���ǂݍ��݂̌��ʂ�ێ�����
	std::vector<MGPosition> vecPos;
	if(IsBinaryFile(in)){
		// Load Binary File
		error = LoadBinary(in, vecPos);
	}else{
		// Load Ascii File
		error = LoadAscii(in, vecPos);
	}
	
	if(error){ // 0:����ɓǂݍ��݂����� ����ȊO�̒l:�ǂݍ��݂����s
		// �ǂݍ��݂����s�����ꍇ�A�G���[�R�[�h��ԋp
		return error;
	}

	// �g�������X��ݒ�
	double wc_zero_old = MGTolerance::set_wc_zero(m_box.length()*MGTolerance::rc_zero());
	// �ǂݍ��񂾍��W�l����m_indices, m_vecNormlTriang��ݒ�
	set_mesh_data(vecPos);
	MGTolerance::set_wc_zero(wc_zero_old);

	// �G���[�R�[�h��ԋp(����ɓǂݍ��݂�����)
	return error;
}

// �t�@�C������ǂݍ��񂾑S�Ă̒��_�̍��W�l���i�[����Ă���
// vecPos����e�O�p�`�̖@�����v�Z���Am_vecNormlTriang�ɒǉ�
// vecPos�̊e�v�f�̍��W�l�̏d�����g�������X�����ɔ��f��
// �d������菜���Am_vecPos�ɒǉ�
// vecPos�̊e�v�f�̃C���f�b�N�X��m_indices�ɒǉ�
// �������:m_vecNormlTriang, m_indices�ɒl���ݒ肳���
void MGStl::set_mesh_data(
	const std::vector<MGPosition>& vecPos // [in]:�t�@�C������ǂݍ��񂾍��W�l�̔z��
){
	m_box.set_null();

	// �t�@�C������ǂ񂾒��_������O�p�`�̐����擾
	size_t nTriang = vecPos.size()/3;
	// ���W�l�ƒ��_�̃C���f�b�N�X���ꎞ�I�ɕێ����Ă����}�b�v
	triangleMap VertexMap;
	for(size_t i = 0; i < nTriang; i++){
		// �O�p�`�̒��_�ԍ��̔z��̃C���f�b�N�X�����߂Ă���
		size_t i3=i*3;
		const MGPosition& p0 = vecPos[i3++]; m_box.expand(p0);
		const MGPosition& p1=vecPos[i3++]; m_box.expand(p1);
		const MGPosition& p2=vecPos[i3]; m_box.expand(p2);
		// �܂��A���͂��ꂽ3�_����ʂ̖@�������߁Am_vecNormlTriang��push_back����
		// ���ɁA�w�肵�����_�̍��W�����ɕۑ�����Ă��邩����
		// �ۑ�����Ă���ꍇ�́A�Y�����钸�_�̃C���f�b�N�X���󂯎��
		// �����m_indices��push_back����
		// �ۑ�����Ă��Ȃ��ꍇ�́A���_��m_vecPos��push_back��
		// ���ɃC���f�b�N�X��+1���A������󂯎��m_indices��push_back����
		push_back_triangle(p0, p1, p2, VertexMap);
	}
}

// Ascii�`���̃t�@�C����ǂݍ��ݑS�Ă̍��W�l���擾����
// �܂��}�`�̃{�b�N�X�g��ݒ肷��
// �߂�l: =0�t�@�C���̓ǂݍ�������
//		  !=0 �ǂݍ��݂Ɏ��s����(std::ifstream�̃G���[�R�[�h)
// ���O����:���ɃI�[�v�����ꂽ�t�@�C���X�g���[����n��
// �������:vecPos�ɑS�Ă̍��W�l���i�[����A box�Ƀ{�b�N�X�g�̍��W�l���ݒ肳���
// �t�@�C���X�g���[�����i��
int MGStl::LoadAscii(
	std::ifstream& in, // [in/out]:���ɃI�[�v�����ꂽ�t�@�C���X�g���[��
	std::vector<MGPosition>& vecPos // [out]:�t�@�C������ǂݍ��񂾍��W�l�̔z��
){
	// Keep one line data read from file
	std::string strLine;
	// Keep cordinate of vertex
	MGPosition position1(3), position2(3), position3(3);

	if(in.eof()){ // �t�@�C������̏ꍇ
		return std::ios_base::eofbit;
	}

	m_box.set_null();

	// �t�@�C������facet normal�����邩�������ϐ�
	bool bStl = false;
	// �G���[�`�F�b�N�̌��ʂ��i�[����ϐ�
	int readState;
	// Read file and set date to member
	while(getline(in, strLine)){
		// �ǂݍ��݃G���[�`�F�b�N
		readState=in.rdstate();
		if(readState == std::ios::failbit || readState == std::ios::badbit){
			return readState; // �G���[�R�[�h��ԋp(�ǂݍ��݂����s)
		}

		if(strLine.find("facet normal") == -1){// �s����facet normal�����݂��邩
			continue;
		}
		// �t�@�C������facet normal�𔭌�
		bStl = true;
		// Triangle data
		std::string dummy;
		// outer loop�s��ǂݔ�΂�
		in.ignore(10000, '\n');
		// Get cordinate of vertex
		in >> dummy >> position1(0) >> position1(1) >> position1(2);
		in >> dummy >> position2(0) >> position2(1) >> position2(2);
		in >> dummy >> position3(0) >> position3(1) >> position3(2);

		// �ǂݍ��݃G���[�`�F�b�N
		readState = in.rdstate();
		if(readState){
			return readState; // �G���[�R�[�h��ԋp(�ǂݍ��݂����s)
		}

		// �t�@�C������ǂݍ��񂾍��W�l��������vector��push_back����
		vecPos.push_back(position1);
		vecPos.push_back(position2);
		vecPos.push_back(position3);
		// �{�b�N�X�g��3�����Ƃ��Aexpand
		m_box.expand(position1);
		m_box.expand(position2);
		m_box.expand(position3);
		// ���s����
		in.seekg(1, std::ios::cur);
		// endloop endfacet�s��ǂݔ�΂�
		in.ignore(10000, '\n');
		in.ignore(10000, '\n');
	}

	if(!bStl){
		// facet normal����x���ǂ܂��ɖ��s�܂œǂݍ��񂾏ꍇ
		// �܂�A�ʂ̎�ނ̃t�@�C����ǂ񂾏ꍇ�̏���
		return std::ios_base::eofbit;
	}

	// ����I���������l��ԋp
	return 0;
}

// Binary�`���̃t�@�C����ǂݍ��ݑS�Ă̍��W�l���擾����
// �܂��}�`�̃{�b�N�X�g��ݒ肷��
// �߂�l: =0�t�@�C���̓ǂݍ�������
//		  !=0 �ǂݍ��݂Ɏ��s����(std::ifstream�̃G���[�R�[�h)
// ���O����:���ɃI�[�v�����ꂽ�t�@�C���X�g���[����n��
// �������:vecPos�ɑS�Ă̍��W�l���i�[����A m_box�Ƀ{�b�N�X�g�̍��W�l���ݒ肳���
//			�t�@�C���X�g���[�����i��
int MGStl::LoadBinary(
	std::ifstream& in, // [in/out]:���ɃI�[�v�����ꂽ�t�@�C���X�g���[��
	std::vector<MGPosition>& vecPos //[out]:�t�@�C������ǂݍ��񂾍��W�l�̔z��
){
	// Keep the point of vertex of triangle
	MGPosition position1(3), position2(3), position3(3);	

	// Skip the top of 80 byte
	in.seekg(80, std::ios::beg);
	// Get count of triangle
	int nTriang;
	in.read((char*)&nTriang, 4);

	// �������݃G���[�`�F�b�N
	int readState = in.rdstate();
	if(readState){
		return readState; // �G���[�R�[�h��ԋp(�ǂݍ��݂����s)
	}

	// Keep vertex and normal of triangle temporary 
	float fVertex1[3], fVertex2[3], fVertex3[3];

	m_box.set_null();
	// Read and set cordinate of vertex and normal of triangle
	for(int i = 0; i < nTriang; i++){
		// normal���i�[����Ă��镔����ǂݔ�΂�
		in.seekg(12, std::ios::cur);
		// Read cordinate of vertex from file
		for(int j = 0; j < 3; j++){
			in.read((char*)&fVertex1[j], 4);
			position1(j) = fVertex1[j];
		}
		for(int j = 0; j < 3; j++){
			in.read((char*)&fVertex2[j], 4);
			position2(j) = fVertex2[j];
		}
		for(int j = 0; j < 3; j++){
			in.read((char*)&fVertex3[j], 4);
			position3(j) = fVertex3[j];
		}

		// �t�@�C������ǂݍ��񂾍��W�l��������vector��push_back����
		vecPos.push_back(position1);
		vecPos.push_back(position2);
		vecPos.push_back(position3);
		// �{�b�N�X�g��expand
		m_box.expand(position1);
		m_box.expand(position2);
		m_box.expand(position3);

		if(i == nTriang - 1){
			// �Ō�̎O�p�`�̏ꍇ�͓ǂݍ��ݏI��
			break;
		}
		// ���̎O�p�`������ꍇ�A2�o�C�g�X�L�b�v����
		in.seekg(2, SEEK_CUR);

		// �ǂݍ��݃G���[�`�F�b�N
		readState = in.rdstate();
		if(readState){
			return readState; // �G���[�R�[�h��ԋp(�ǂݍ��݂����s)
		}
	}

	// �������݃G���[�`�F�b�N
	readState = in.rdstate();
	if(readState == std::ios::failbit || readState == std::ios::badbit){
		return readState; // �G���[�R�[�h��ԋp(�ǂݍ��݂����s)
	}

	// ����I���������l��ԋp
	return 0;
}

// �w�肳�ꂽ�p�X��Ascii�`����STL�t�@�C����ۑ�����
// �߂�l: =0�t�@�C���̏������݂�����
//		  !=0 �������܂Ȃ������B�܂��͎��s����(std::ofstream�̃G���[�R�[�h)
// ���O����:m_vecPos��m_vecNorml��m_indices�ɐ�������񂪓�����Ă���
// �������:rSTLFilePath�Ɏw�肵���p�X��Ascii�`����STL�t�@�C�����ۑ������
//�d�l�ύX�F�����ɃX�g���[����n�����悤�ɕύX
int MGStl::SaveAscii(
	std::ofstream& fout
	//const char* rSTLFilePath // [in]:�t�@�C���̕ۑ���̃p�X
)const{
	//std::ofstream fout(rSTLFilePath);

	// �������݃G���[�`�F�b�N
	int writeState = fout.rdstate();
	if(writeState){
		return writeState; // �G���[�R�[�h��ԋp(�ǂݍ��݂����s)
	}

	fout.setf(std::ios::scientific);
	// solid
	fout << "solid ascii" << std::endl;
	// �O�p�`�̐����擾(�O�p�`�̐� == �ʖ@���x�N�g���̌�)
	size_t nTriang = m_vecNormlTriang.size();
	// output normal of triangle and cordinate of vertex to file
	for(size_t i = 0; i < nTriang; i++){
		// Get triangle normal
		const MGUnit_vector& normal = m_vecNormlTriang[i];
		// facet normal
		fout << "    facet normal  " << normal(0) << " " << normal(1)<< " " << normal(2) << " " << std::endl;
		// outerloop
		fout << "        outer loop" << std::endl;
		
		// ���_�̔ԍ��̔z��̓Y�������v�Z���Ă���
		size_t index=i*3;
		// output cordinate of vertex to file
		const MGPosition& position1 = m_vecPos[m_indices[index++]];
		fout << "        vertex  " << position1(0) << " " << position1(1) << " " << position1(2) << " " << std::endl;
		const MGPosition& position2 = m_vecPos[m_indices[index++]];
		fout << "        vertex  " << position2(0) << " " << position2(1) << " " << position2(2) << " " << std::endl;
		const MGPosition& position3 = m_vecPos[m_indices[index]];
		fout << "        vertex  " << position3(0) << " " << position3(1) << " " << position3(2) << " " << std::endl;
		// outerloop
		fout << "        endloop" << std::endl;
		fout << "    endfacet" << std::endl;

		// �������݃G���[�`�F�b�N
		writeState = fout.rdstate();
		if(writeState){
			return writeState; // �G���[�R�[�h��ԋp(�ǂݍ��݂����s)
		}
	}
	// endsolid
	fout << "endsolid" << std::endl;
	fout.unsetf(std::ios::scientific);

	writeState = fout.rdstate();
	if(writeState == std::ios::failbit || writeState == std::ios::badbit){
		return writeState; // �G���[�R�[�h��ԋp(�������݂����s)
	}

	//fout.close();
	return 0;
}

// �w�肳�ꂽ�p�X��Binary�`����STL�t�@�C����ۑ�����
// �߂�l: =0�t�@�C���̏������݂�����
//		  !=0 �������܂Ȃ������B�܂��͎��s����(std::ofstream�̃G���[�R�[�h)
// ���O����:m_vecPos��m_vecNorml��m_indices�ɐ�������񂪓�����Ă���
// �������:rSTLFilePath�Ɏw�肵���p�X��Binary�`����STL�t�@�C�����ۑ������
int MGStl::SaveBinary(
	const TCHAR* rSTLFilePath  // [in]:�t�@�C���̕ۑ���̃p�X
)const{
	// File open
	std::ofstream fout(rSTLFilePath, std::ios::binary);
	// �X�g���[���̃G���[�`�F�b�N
	int writeState = fout.rdstate();
	if(writeState)
		return writeState;

	// �O�p�`�̐����擾(�O�p�`�̐� == �ʖ@���x�N�g���̌�)
	int nTriang = (int)m_vecNormlTriang.size();
	// �t�@�C���擪80�o�C�g���X�L�b�v
	fout.seekp(80, std::ios::beg);
	// �t�@�C���ɎO�p�`�̖�������������
	fout.write((char*)&nTriang, 4);

	// �������݃G���[�`�F�b�N
	writeState = fout.rdstate();
	if(writeState){
		return writeState; // �G���[�R�[�h��ԋp(�������݂Ɏ��s)
	}

	// Keep vertex and normal of triangle temporary 
	float fNormal, fVertex;
	// �e�O�p�`�̏�����������
	for(int i = 0; i < nTriang; i++){
		// �O�p�`�̖ʖ@���x�N�g������������
		const MGUnit_vector& normal = m_vecNormlTriang[i];
		for(int j = 0; j < 3; j++){
			fNormal = float(normal[j]);
			fout.write((char*)&fNormal, 4);
		}

		// �O�p�`�̒��_�ԍ��̔z��̃C���f�b�N�X�����߂Ă���
		int i3 = i*3;

		// �x�N�^�[���璸�_�̃C���f�b�N�X���w���Ă�����W�l�����o��
		// �O�p�`�̒��_�̍��W����������
		const MGPosition& position1 = m_vecPos[m_indices[i3]];
		for(int j = 0; j < 3; j++){
			fVertex = float(position1[j]);
			fout.write((char*)&fVertex, 4);
		}
		
		const MGPosition& position2 = m_vecPos[m_indices[i3+1]];
		for(int j = 0; j < 3; j++)		{
			fVertex = float(position2[j]);
			fout.write((char*)&fVertex, 4);
		}
		
		const MGPosition& position3 = m_vecPos[m_indices[i3+2]];
		for(int j = 0; j < 3; j++){
			fVertex = float(position3[j]);
			fout.write((char*)&fVertex, 4);
		}

		// �Ō�̎O�p�`�̏ꍇ�͏������݂��I��
		if(i == nTriang - 1){
			break;
		}

		// �O�p�`���Ƃ�2�o�C�g�Ԋu���󂯂�
		fout.seekp(2, std::ios::cur);

		// �������݃G���[�`�F�b�N
		writeState = fout.rdstate();
		if(writeState){
			return writeState; // �G���[�R�[�h��ԋp(�������݂Ɏ��s)
		}
	}
	// �������݃G���[�`�F�b�N
	writeState = fout.rdstate();
	if(writeState == std::ios::failbit || writeState == std::ios::badbit){
		return writeState;// �G���[�R�[�h��ԋp(�������݂����s)
	}

	fout.close();
	return 0;
}

// ���͂ł���position��VertexMap�ɂ��łɓo�^����Ă���΂���m_vecPos�Y������Ԃ�
// ���o�^�̏ꍇ�A�V����m_vecPos�Ɋi�[���Aposition��m_vecPos�Y������map��VertexMap��
// �o�^����
int MGStl::IdentifyPosition(
	const MGPosition& position, // ���_�̍��W
	triangleMap& VertexMap // �d���̂Ȃ����_�̍��W�A���_�̃C���f�b�N�X��ێ�����
){
	m_box.expand(position);

	int nVertexIndex = (int)m_vecPos.size();
	//iterator��position�̍��W�l�ȏ�ƂȂ�v�f���͂��߂Č����ꏊ���w���悤�ɂ���
	triangleMap::iterator itAft = VertexMap.lower_bound(position), itPre;
	if(itAft != VertexMap.end()){	//iterator���I�����w���Ă��Ȃ��ꍇ
		if(itAft->first==position)
			return itAft->second;	//�w���Ă���v�f�ƈ�v�����ꍇ�C����index��Ԃ�
	}
	if(itAft != VertexMap.begin()){	//iterator���n�߂��w���Ă��Ȃ��ꍇ
		itPre=itAft;itPre--;	//��O�̗v�f��position���r����
		if(itPre->first==position){
			return itPre->second;	//��O�̗v�f�ƈ�v�����ꍇ�C����index��Ԃ�
		}
	}

	VertexMap.insert(itAft, std::make_pair(position, nVertexIndex));
	m_vecPos.push_back(position);
	return nVertexIndex;
}

// �����o�f�[�^���������ފ֐�
void MGStl::WriteMembers(MGOfstream& buf)const{
	MGObject::WriteMembers(buf);
	// �O�p�`�̖@���x�N�g������������
	int nTriangle = (int)m_vecNormlTriang.size();
	buf << nTriangle;
	for(int i = 0; i < nTriangle; i++){
		m_vecNormlTriang[i].dump(buf);
	}

	// �O�p�`�̏d�����Ȃ����W�̒��_����������
	int nPos = (int)m_vecPos.size();
	buf << nPos;
	for(int k = 0; k < nPos; k++){
		m_vecPos[k].dump(buf);
	}

	// �O�p�`�̊e���_�̃C���f�b�N�X����������
	int nIndex = (int)m_indices.size();
	buf << nIndex;
	for(int l = 0; l < nIndex; l++){
		buf << m_indices[l];
	}
}

// �����o�f�[�^��ǂݍ��ފ֐�
void MGStl::ReadMembers(MGIfstream& buf){
	MGObject::ReadMembers(buf);
	// �O�p�`�̖@���x�N�g����ǂݍ���
	int nTriangle;
	buf >> nTriangle;
	m_vecNormlTriang.resize(nTriangle);
	for(int i = 0; i < nTriangle; i++){	
		m_vecNormlTriang[i].restore(buf);
	}
	

	// �O�p�`�̊e���_�̍��W���d�������ǂݍ���
	m_box.set_null();
	int nPos;
	buf >> nPos;
	m_vecPos.resize(nPos);
	for(int k = 0; k < nPos; k++){
		MGPosition& Pk=m_vecPos[k];
		Pk.restore(buf);
		m_box.expand(Pk);
	}

	// �O�p�`�̊e���_�̃C���f�b�N�X��ǂݍ���
	int nIndex;
	buf >> nIndex;
	m_indices.resize(nIndex);
	for(int l = 0; l < nIndex; l++){
		buf >> m_indices[l];
	}
}

// ��x�t�@�C����ǂ݁A�����o�ɒl���ݒ肳�ꂽ��ԂōĂуt�@�C����ǂނ�
// �����o�̒l���㏑�������B���������邽�߁A�����o�̒l��S�ă��Z�b�g����
void MGStl::Initialize(){
	m_vecPos.clear();
	m_vecNormlTriang.clear();
	m_indices.clear();
}

// 3�_����쐬�����O�p�`�̏��������o�ϐ��ɒl��ݒ肷��
// ���O����:���͂���pos1, pos2, pos3�͔����v���ɂȂ��Ă��邱��
// �������:m_vecPos, m_vecNormlTriang, m_indices�ɎO�p�`�̏�񂪐ݒ肳���
// �܂��A���͂��ꂽ3�_����ʂ̖@�������߁Am_vecNormlTriang��push_back����
// ���ɁA�w�肵�����_�̍��W�����ɕۑ�����Ă��邩����
// �ۑ�����Ă���ꍇ�́A�Y�����钸�_�̃C���f�b�N�X���󂯎��
// �����m_indices��push_back����
// �ۑ�����Ă��Ȃ��ꍇ�́A���_��m_vecPos��push_back��
// ���ɃC���f�b�N�X��+1���A������󂯎��m_indices��push_back����
void MGStl::push_back_triangle(
	const MGPosition& pos1,
	const MGPosition& pos2,
	const MGPosition& pos3,
	triangleMap& VertexMap
){
	// 3�_����ʂ̖@�������߁A�����o�ϐ���push_back����
	m_vecNormlTriang.push_back(UnitNormal(pos1, pos2, pos3));

	// m_indices�ɒ��_�̃C���f�b�N�X��push_back��
	// m_vecPos�ɏd���Ȃ�MGPosition��push_back����
	m_indices.push_back(IdentifyPosition(pos1, VertexMap)); m_box.expand(pos1);
	m_indices.push_back(IdentifyPosition(pos2, VertexMap)); m_box.expand(pos2);
	m_indices.push_back(IdentifyPosition(pos3, VertexMap)); m_box.expand(pos3);
}

// �����o�[�f�[�^�𒼐ڃZ�b�g����
void MGStl::set_all_data(
	const std::vector<MGPosition>& vertices,
	const std::vector<int>& indices
){
	assert(indices.size() % 3 == 0);

	m_vecPos = vertices;
	m_indices = indices;

	// update normals and box
	compute_box(m_box);
	update_normals();
}

// Update the all normals of the triangles.
void MGStl::update_normals(){
	size_t nTri = m_indices.size() / 3;
	std::vector<MGUnit_vector> work(nTri);

	for(size_t i = 0; i < nTri; ++i){
		int vid[3];
		GetVertIndices((int)i, vid);
		assert(size_t(vid[0]) < m_vecPos.size());
		assert(size_t(vid[1]) < m_vecPos.size());
		assert(size_t(vid[2]) < m_vecPos.size());

		work[i] = UnitNormal(
			m_vecPos[vid[0]],
			m_vecPos[vid[1]],
			m_vecPos[vid[2]]);
	}

	work.swap(m_vecNormlTriang);
}

// �ǂݍ���stl�f�[�^��obj�t�H�[�}�b�g�ɕϊ����A�o�͂���
// �߂�l: =0 �t�@�C���̏������݂�����
//        !=0 �������܂Ȃ������B�܂��͎��s����(std::ofstream�̃G���[�R�[�h)
// �������:�w�肵���p�X��obj�t�@�C�����ۑ������.
int MGStl::SaveObjFormatFromStl(
	std::ofstream& fout
)const{

	// �������݃G���[�`�F�b�N
	int writeState = fout.rdstate();
	if(writeState){
		return writeState; // �G���[�R�[�h��ԋp(�ǂݍ��݂����s)
	}

	// �����t���O��ݒ�
	// TODO ����̖ړI�̓z�C�[���p�̃f�[�^(�P��:mm)���o�͂���̂ŁA
	// �����t���O��fixed�ɂ��Ă���
	fout.setf(std::ios::fixed);

	// stl�f�[�^�œǂݍ��񂾎O�p�`�̌����擾
	int nTriang = GetTriangleCount();

	// ���_���W�̏o��
	for(int i0 = 0; i0 < nTriang; i0++){
		// ���_�̔ԍ��̔z��̓Y�������v�Z���Ă���
		int index = i0 * 3;

		// (�O�p�`��)���_���W���擾���A�o�͂���
		for(int i1 = 0; i1 < 3; i1++){
			const MGPosition& position = m_vecPos[m_indices[index + i1]];
			fout << "v " << position(0) << " " << position(1) << " " << position(2) << std::endl;
		}
	}

	// �@���x�N�g���̏o��
	for(int j = 0; j < nTriang; j++){
		// �@���x�N�g���̎擾
		const MGUnit_vector& normal = m_vecNormlTriang[j];
		fout << "vn " << normal(0) << " " << normal(1) << " " << normal(2) << std::endl;
	}

	// �ʏ��̏o��
	for(int k = 0; k < nTriang; k++){
		// stl�`���̃f�[�^��ǂݍ���ł��邽�߁A�ʃf�[�^�͑S�ĎO�p�`�݂̂ɂȂ��Ă���
		// �o�͂����f�[�^�̌`���́A
		// v1//n1 v2//n1 v3//n1 �݂̂ł���
		int n1 = k + 1; // �@���x�N�g��
		fout << "f " << k * 3 + 1 << "//" << n1 << " " 
			         << k * 3 + 2 << "//" << n1 << " " 
					 << k * 3 + 3 << "//" << n1 << std::endl; 
	}

	// �����t���O������
	fout.unsetf(std::ios::fixed);

	// �������݃G���[�`�F�b�N
	writeState = fout.rdstate();
	if(writeState == std::ios::failbit || writeState == std::ios::badbit){
		return writeState; // �G���[�R�[�h��ԋp(�������݂����s)
	}

	return 0; // ����I��
}
