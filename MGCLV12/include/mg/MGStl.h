/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#ifndef _MGStl_HH_
#define _MGStl_HH_

#include <map>
#include "mg/drawParam.h"
#include "mg/object.h"
#include "mg/Curve.h"
#include "mg/Surface.h"
#include "mg/FSurface.h"
#include "topo/Face.h"
#include "topo/Shell.h"

class mgTL2Triangles;
class mgVBO;

/** @addtogroup MGObjectRelated
 *  @{
 */

//Define MGStl Class.

///MGStl is a concrete class which represents an STL information.

///���̃N���X�͕����̎O�p�`��z��Ƃ��ĕێ�����.
///�S�Ă̎O�p�`�̒��_�̍��W�̓t�@�C������ǂݍ��񂾏��ŁA�d������菜���Ĕz��Ɋi�[�����.
///m_indices�ɂ̓t�@�C������ǂݍ��񂾏��ԂŊe�O�p�`�̒��_�̕��т��i�[����Ă���
///���̊e�v�f�ɂ͊Y�����钸�_���W�̔z��̓Y�������i�[�����.
///�Ⴆ�΁Ai�Ԗڂ̎O�p�`�̊e���_�̍��W�̓Y������indices[i*3]�Aindices[i*3+1]�Aindices[i*3+2]
///�Ƃ����ȏ�̏����Ŋi�[����Ă���.
class MG_DLL_DECLR MGStl : public MGObject{

/// map�ɂ�1�̃L�[�ɑ΂��ĕ����̒l�͎��߂��Ȃ�����
/// �O�p�`���\�����钸�_��ID�����߂��\���̂𗘗p����
class vertId{
public:
	int id1;
	int id2;
	int id3;
};

/// �֐��I�u�W�F�N�g���i�[����\����
class positionComp{
public:
	bool operator()(const MGPosition& p1, const MGPosition& p2) const{
		/// lexicographical_compare style
		if(p1(0) < p2(0)){
			return true;
		}else if(p2(0) < p1(0)){
			return false;
		}
		if(p1(1) < p2(1)){
			return true;
		}else if(p2(1) < p1(1)){
			return false;
		}
		if(p1(2) < p2(2)){
			return true;
		}else if(p2(2) < p1(2)){
			return false;
		}
		/// p1 == p2
		return false;
	}		
};

/// �O�p�`�̒��_�̍��W�A���_�̔ԍ����i�[����map�̕ʖ����`.
typedef std::map<MGPosition, int, positionComp> triangleMap;

/// ���_�̃C���f�b�N�X�A���_�̍��W���i�[����map�̕ʖ����`.
typedef std::map<int, MGPosition> IndexPosMap;

/// �O�p�`�̃C���f�b�N�X�A3�̒��_�̃C���f�b�N�X���i�[����map�̕ʖ����`.
typedef std::map<int, vertId> TriVertMap;

public:

////////Special member functions/////////
MGStl()=default;
~MGStl()=default;
MGStl(const MGStl&)=default;///Copy constructor.
MGStl& operator= (const MGStl&)=default;///Copy assignment.
MGStl(MGStl&&)=default;		///Move constructor.
MGStl& operator= (MGStl&&)=default;///Move assignment.


///conversion constructor from tessellation data.
MGStl(
	const mgTL2Triangles& tris///<mgTL2Triangles whose data depend on tris.get_kind();
);
MGStl(
	double error,	///<Error to regard two points are the same.
	const mgTL2Triangles& tris///<mgTL2Triangles whose data depend on tris.get_kind();
);

///conversion constructor from tessellation data.
MGStl(const std::vector<mgTL2Triangles>& tlDataVector);

///Constructor from triangle data, index+vertices.
///This constructor uses the wc_zero to identify
///different two positions as the same input position.
MGStl(
	int nTriang, ///< �O�p�`�̐�
	const int* triang,///< ���_��id. The id's start from 1.
		///<(i*3, i*3+1, i*3+2) for i=0,...,indices.size()/3.
	const double* verts ///< ���_�̍��W. verts[id0], [id1], [id2] make
		///< a triangle for id0=triang[i*3], id1=triang[i*3+1], id2=triang[i*3+2]
		///<  for i=0,...,nTriang.
);

///Constructor from triangle data, index+vertices.
///This constructor uses the wc_zero to identify
///different two positions as the same input position.
MGStl(
	const std::vector<MGPosition>& vertices,///< vertices[id0], [id1], [id2] make
		///<a triangle for id0=indices[i*3], id1=indices[i*3+1], id2=indices[i*3+2]
		///<for i=0,...,indices.size()/3.
	const std::vector<int>& indices ///< ���_��id. The id's start from 1.
		///<  (0,1,2),...(i*3, i*3+1, i*3+2) for i=0,...,indices.size()/3.
);

/// �S�Ă̒��_�ƃ{�b�N�X�̍��W�Ɏw�肵��MGVector�̒l�����Z.
MGStl& operator+=(const MGVector& v);

/// �S�Ă̒��_�ƃ{�b�N�X�̍��W����w�肵��MGVector�̒l������.
MGStl& operator-=(const MGVector& v);

/// �S�Ă̒��_�ƃ{�b�N�X�̍��W�Ɏw�肵���l���|����.
MGStl& operator*=(double scale);

/// �^����ꂽ�ϊ����s��.
MGStl& operator*=(const MGMatrix& mat);

/// �^����ꂽ�ϊ��ɂ��g�����X�t�H�[�����s��.
MGStl& operator*=(const MGTransf& tr);

/// 2��MGStl�I�u�W�F�N�g�������������肷��.
bool operator==(const MGStl& stl);

std::ostream& toString(std::ostream& ostrm)const;

///Construct new object by copying to newed area.
///User must delete this copied object by "delete".
MGStl* clone()const;

///Compute box of the geometry.
///Compute the box from the scratch.
void compute_box(MGBox& bx) const;

/// Return This object's typeID.
long identify_type()const;

///Draw 3D curve in world coordinates.
///The object is converted to curve(s) and is drawn.
void drawWire(
	mgVBO& vbo,///<The target graphic object.
	int line_density=1	///<line density to draw a surface in wire mode.
)const;

///Shade the object in world coordinates.
void shade(
	mgVBO& vbo,///<The target graphic object.
	const MGDrawParam& para,///<Parameters to draw.
	MGCL::DRAW_TARGET target= MGCL::SHADING///<The target vbo element.
)const;

///Get manifold dimension.
int manifold_dimension()const{return 2;};

///Write all member data.
void WriteMembers(MGOfstream& buf)const;

///Read all member data.
void ReadMembers(MGIfstream& buf);

///return the ref to m_vecPos.
const std::vector<MGPosition>& positions()const{return m_vecPos;};
std::vector<MGPosition>& positions(){return m_vecPos;};

///return the ref to m_vecPos.
const std::vector<MGUnit_vector>& normals()const{return m_vecNormlTriang;};
std::vector<MGUnit_vector>& normals(){return m_vecNormlTriang;};

/// �O�p�`�̌����擾����.
/// �߂�l�F�O�p�`�̌�.
int GetTriangleCount()const{return int(m_vecNormlTriang.size());};

/// �w�肵���O�p�`�̒��_���W��������Ă���z��̃C���f�b�N�X���擾����.
void GetVertIndices(
	int i, ///< [in]�F�O�p�`�̃C���f�b�N�X(0 <= i < GetTriangleCount)
	int pos[3] ///< [out]�F�w�肵���O�p�`�̊e���_���W�̔z��̓Y����[i, j, k]
)const;

///Intersection of a STL and an object, not supported.
MGisects isect(const MGObject& obj2)const override{
	return MGisects();
};

/// �����Ŏw�肵���p�X��STL�t�@�C����ǂݍ��݃����o�ɒl��ݒ肷��.

/// �܂�MGTolerance::wc_zero�ɒl��ݒ肷��.
/// �߂�l: =0�t�@�C���̓ǂݍ�������
///		  !=0 �ǂݍ��܂Ȃ������B�܂��͎��s����(std::ifstream�̃G���[�R�[�h�j
/// �������:m_vecPos, m_vecNorml, m_indices, m_box�ɐ}�`�̏�񂪊i�[�����.
///			MGTolerance::wc_zero�ɒl���ݒ肳���.
int LoadFile(
	const TCHAR* strFilePath ///< [in]:�ǂݍ���STL�t�@�C���ւ̃p�X
);

/// �w�肳�ꂽ�p�X��Ascii�`����STL�t�@�C����ۑ�����.

/// �߂�l: =0�t�@�C���̏������݂�����
///		  !=0 �������܂Ȃ������B�܂��͎��s����(std::ofstream�̃G���[�R�[�h)
/// �������:rSTLFilePath�Ɏw�肵���p�X��Ascii�`����STL�t�@�C�����ۑ������.
int SaveAscii(
	std::ofstream& fout
)const;

/// �w�肳�ꂽ�p�X��Binary�`����STL�t�@�C����ۑ�����.

/// �߂�l: =0�t�@�C���̏������݂�����
///		  !=0 �������܂Ȃ������B�܂��͎��s����(std::ofstream�̃G���[�R�[�h)
/// �������:rSTLFilePath�Ɏw�肵���p�X��Binary�`����STL�t�@�C�����ۑ������.
///�d�l�ύX�F�����ɃX�g���[����n�����悤�ɕύX09/09/17.
int SaveBinary(
	const TCHAR* rSTLFilePath ///< [in]:�t�@�C���̕ۑ���̃p�X
)const;

/// 3�_����쐬�����O�p�`�̏��������o�ϐ��ɒl��ݒ肷��

/// ���O����:���͂���pos1, pos2, pos3�͔����v���ɂȂ��Ă��邱��
/// �������:m_vecPos, m_vecNormlTriang, m_indices�ɎO�p�`�̏�񂪐ݒ肳���
/// �܂��A���͂��ꂽ3�_����ʂ̖@�������߁Am_vecNormlTriang��push_back����
/// ���ɁA�w�肵�����_�̍��W�����ɕۑ�����Ă��邩����
/// �ۑ�����Ă���ꍇ�́A�Y�����钸�_�̃C���f�b�N�X���󂯎��
/// �����m_indices��push_back����
/// �ۑ�����Ă��Ȃ��ꍇ�́A���_��m_vecPos��push_back��
/// ���ɃC���f�b�N�X��+1���A������󂯎��m_indices��push_back����
void push_back_triangle(
	const MGPosition& pos1,
	const MGPosition& pos2,
	const MGPosition& pos3,
	triangleMap& VertexMap
);

/// �O�p�`���Ƃɖ@����\������.
void display_arrows(mgSysGL& sgl)const;

/// �����o�[�f�[�^�𒼐ڃZ�b�g����.
void set_all_data(
	const std::vector<MGPosition>& vertices,
	const std::vector<int>& indices);

/// Update the all normals of the triangles.
void update_normals();

///Get the name of the class.
std::string whoami()const{return "Stl";};

/// �ǂݍ���stl�f�[�^��obj�t�H�[�}�b�g�ɕϊ����A�o�͂���
/// �߂�l: =0 �t�@�C���̏������݂�����
///        !=0 �������܂Ȃ������B�܂��͎��s����(std::ofstream�̃G���[�R�[�h)
/// �������:�w�肵���p�X��obj�t�@�C�����ۑ������.
int SaveObjFormatFromStl(std::ofstream& fout)const;

private:

/// ��x�t�@�C����ǂ݁A�����o�ɒl���ݒ肳�ꂽ��ԂōĂуt�@�C����ǂނ�.
/// �����o�̒l���㏑�������B���������邽�߁A�����o�̒l��S�ă��Z�b�g����.
void Initialize();

/// Ascii�`���̃t�@�C����ǂݍ��ݑS�Ă̍��W�l���擾����.
/// �܂��}�`�̃{�b�N�X�g��ݒ肷��.
/// �߂�l: =0�t�@�C���̓ǂݍ�������
///		  !=0 �ǂݍ��݂Ɏ��s����(std::ifstream�̃G���[�R�[�h)
/// ���O����:���ɃI�[�v�����ꂽ�t�@�C���X�g���[����n��.
/// �������:vecPos�ɑS�Ă̍��W�l���i�[����A box�Ƀ{�b�N�X�g�̍��W�l���ݒ肳���.
/// �t�@�C���X�g���[�����i��.
int LoadAscii(
	std::ifstream& in, ///< [in/out]:���ɃI�[�v�����ꂽ�t�@�C���X�g���[�� .
	std::vector<MGPosition>& vecPos ///<[out]:�t�@�C������ǂݍ��񂾍��W�l�̔z��.
	);

/// Binary�`���̃t�@�C����ǂݍ��ݑS�Ă̍��W�l���擾����.
/// �܂��}�`�̃{�b�N�X�g��ݒ肷��.
/// �߂�l: =0�t�@�C���̓ǂݍ�������
///		  !=0 �ǂݍ��݂Ɏ��s����(std::ifstream�̃G���[�R�[�h)
/// ���O����:���ɃI�[�v�����ꂽ�t�@�C���X�g���[����n��.
/// �������:vecPos�ɑS�Ă̍��W�l���i�[����A m_box�Ƀ{�b�N�X�g�̍��W�l���ݒ肳���.
///			�t�@�C���X�g���[�����i��
int LoadBinary(
	std::ifstream& in, ///< [in/out]:���ɃI�[�v�����ꂽ�t�@�C���X�g���[��.
	std::vector<MGPosition>& vecPos ///<[out]:�t�@�C������ǂݍ��񂾍��W�l�̔z��.
	);

/// �t�@�C������ǂݍ��񂾑S�Ă̒��_�̍��W�l���i�[����Ă���.
/// vecPos����e�O�p�`�̖@�����v�Z���Am_vecNormlTriang�ɒǉ�.
/// vecPos�̊e�v�f�̍��W�l�̏d�����g�������X�����ɔ��f��
/// �d������菜���Am_vecPos�ɒǉ�.
/// vecPos�̊e�v�f�̃C���f�b�N�X��m_indices�ɒǉ�.
/// �������:m_vecNormlTriang, m_indices�ɒl���ݒ肳���.
void set_mesh_data(
	const std::vector<MGPosition>& vecPos ///< [in]:�t�@�C������ǂݍ��񂾍��W�l�̔z��.
);

/// ���͂ł���position��VertexMap�ɂ��łɓo�^����Ă��邩�`�F�b�N����.
/// �o�^����Ă���΂���m_vecPos�Y������Ԃ�,
/// ���o�^�̏ꍇ�A�V����m_vecPos�Ɋi�[���Aposition��m_vecPos�Y������map��VertexMap��
/// �o�^����.
int IdentifyPosition(
	const MGPosition& position, ///< [in]:���_�̍��W.
	triangleMap& VertexMap ///< [in/out]:���_�̍��W�A���_�̃C���f�b�N�X��ێ�����.
);

///mgTL2Triangles��MGStl�ɒǉ�����
void AddTL2Data(
	const mgTL2Triangles& tris,///<mgTL2Triangles whose data depend on tris.get_kind();
	triangleMap& VertexMap
);

///mgTL2Triangles��MGStl�ɒǉ�����
void AddTL2Data(
	double error,	///<Error to regard two points are the same.
	const mgTL2Triangles& tris,///<mgTL2Triangles whose data depend on tris.get_kind();
	triangleMap& VertexMap
);

/// STL�t�@�C���̐}�`���\������e�O�p�`�́A���W�̏d�����Ȃ����_���W�̔z��.
/// �t�@�C������ǂݍ��񂾏��ŁA���W�̏d������菜���A���W�l���i�[����Ă���.
std::vector<MGPosition> m_vecPos;

/// STL�t�@�C���̐}�`���\������e�O�p�`�̖@���x�N�g���̔z��.
/// �t�@�C������ǂݍ��܂ꂽ���ŎO�p�`�̖@���x�N�g�����i�[����Ă���.
///m_vecNormlTriang.size()*3=m_indices.size().
///m_vecNormlTriang[i] is the normal of the triangle  m_indices[i], [i+1], [i+2] for
///i=0, ..., m_indices.size()/3.
std::vector<MGUnit_vector> m_vecNormlTriang;

/// STL�t�@�C���̐}�`���\������e�O�p�`�̊e���_�ɑΉ�����
/// ���W�̔z��̃C���f�b�N�X���i�[����z��.
/// �t�@�C������ǂݍ��񂾏��ԂŊe�O�p�`�̒��_�̕��т��i�[����Ă���.
/// �e�v�f�ɂ͊Y�����钸�_���W�̔z��̓Y�������i�[����Ă���.
/// ��Fi�Ԗڂ̎O�p�`�̊e���_�̍��W�̔z��̓Y������
/// m_indices[i*3]�A[i*3+1]�A[i*3+2]�Ƃ����ȏ�̏����Ŏ擾�ł���.
std::vector<int> m_indices;

};

/** @} */ // end of MGObjectRelated group

#endif
