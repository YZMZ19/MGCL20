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

///このクラスは複数の三角形を配列として保持する.
///全ての三角形の頂点の座標はファイルから読み込んだ順で、重複を取り除いて配列に格納される.
///m_indicesにはファイルから読み込んだ順番で各三角形の頂点の並びが格納されており
///その各要素には該当する頂点座標の配列の添え字が格納される.
///例えば、i番目の三角形の各頂点の座標の添え字はindices[i*3]、indices[i*3+1]、indices[i*3+2]
///という以上の順序で格納されている.
class MG_DLL_DECLR MGStl : public MGObject{

/// mapには1つのキーに対して複数の値は収められないため
/// 三角形を構成する頂点のIDを収めた構造体を利用する
class vertId{
public:
	int id1;
	int id2;
	int id3;
};

/// 関数オブジェクトを格納する構造体
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

/// 三角形の頂点の座標、頂点の番号を格納するmapの別名を定義.
typedef std::map<MGPosition, int, positionComp> triangleMap;

/// 頂点のインデックス、頂点の座標を格納するmapの別名を定義.
typedef std::map<int, MGPosition> IndexPosMap;

/// 三角形のインデックス、3つの頂点のインデックスを格納するmapの別名を定義.
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
	int nTriang, ///< 三角形の数
	const int* triang,///< 頂点のid. The id's start from 1.
		///<(i*3, i*3+1, i*3+2) for i=0,...,indices.size()/3.
	const double* verts ///< 頂点の座標. verts[id0], [id1], [id2] make
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
	const std::vector<int>& indices ///< 頂点のid. The id's start from 1.
		///<  (0,1,2),...(i*3, i*3+1, i*3+2) for i=0,...,indices.size()/3.
);

/// 全ての頂点とボックスの座標に指定したMGVectorの値を加算.
MGStl& operator+=(const MGVector& v);

/// 全ての頂点とボックスの座標から指定したMGVectorの値を引く.
MGStl& operator-=(const MGVector& v);

/// 全ての頂点とボックスの座標に指定した値を掛ける.
MGStl& operator*=(double scale);

/// 与えられた変換を行う.
MGStl& operator*=(const MGMatrix& mat);

/// 与えられた変換によるトランスフォームを行う.
MGStl& operator*=(const MGTransf& tr);

/// 2つのMGStlオブジェクトが等しいか判定する.
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

/// 三角形の個数を取得する.
/// 戻り値：三角形の個数.
int GetTriangleCount()const{return int(m_vecNormlTriang.size());};

/// 指定した三角形の頂点座標が入れられている配列のインデックスを取得する.
void GetVertIndices(
	int i, ///< [in]：三角形のインデックス(0 <= i < GetTriangleCount)
	int pos[3] ///< [out]：指定した三角形の各頂点座標の配列の添え字[i, j, k]
)const;

///Intersection of a STL and an object, not supported.
MGisects isect(const MGObject& obj2)const override{
	return MGisects();
};

/// 引数で指定したパスのSTLファイルを読み込みメンバに値を設定する.

/// またMGTolerance::wc_zeroに値を設定する.
/// 戻り値: =0ファイルの読み込が成功
///		  !=0 読み込まなかった。または失敗した(std::ifstreamのエラーコード）
/// 事後条件:m_vecPos, m_vecNorml, m_indices, m_boxに図形の情報が格納される.
///			MGTolerance::wc_zeroに値が設定される.
int LoadFile(
	const TCHAR* strFilePath ///< [in]:読み込むSTLファイルへのパス
);

/// 指定されたパスにAscii形式のSTLファイルを保存する.

/// 戻り値: =0ファイルの書き込みが成功
///		  !=0 書き込まなかった。または失敗した(std::ofstreamのエラーコード)
/// 事後条件:rSTLFilePathに指定したパスにAscii形式のSTLファイルが保存される.
int SaveAscii(
	std::ofstream& fout
)const;

/// 指定されたパスにBinary形式のSTLファイルを保存する.

/// 戻り値: =0ファイルの書き込みが成功
///		  !=0 書き込まなかった。または失敗した(std::ofstreamのエラーコード)
/// 事後条件:rSTLFilePathに指定したパスにBinary形式のSTLファイルが保存される.
///仕様変更：引数にストリームを渡されるように変更09/09/17.
int SaveBinary(
	const TCHAR* rSTLFilePath ///< [in]:ファイルの保存先のパス
)const;

/// 3点から作成した三角形の情報をメンバ変数に値を設定する

/// 事前条件:入力するpos1, pos2, pos3は反時計回りになっていること
/// 事後条件:m_vecPos, m_vecNormlTriang, m_indicesに三角形の情報が設定される
/// まず、入力された3点から面の法線を求め、m_vecNormlTriangへpush_backする
/// 次に、指定した頂点の座標が既に保存されているか調べ
/// 保存されている場合は、該当する頂点のインデックスを受け取る
/// これをm_indicesへpush_backする
/// 保存されていない場合は、頂点をm_vecPosへpush_backし
/// 次にインデックスを+1し、それを受け取りm_indicesにpush_backする
void push_back_triangle(
	const MGPosition& pos1,
	const MGPosition& pos2,
	const MGPosition& pos3,
	triangleMap& VertexMap
);

/// 三角形ごとに法線を表示する.
void display_arrows(mgSysGL& sgl)const;

/// メンバーデータを直接セットする.
void set_all_data(
	const std::vector<MGPosition>& vertices,
	const std::vector<int>& indices);

/// Update the all normals of the triangles.
void update_normals();

///Get the name of the class.
std::string whoami()const{return "Stl";};

/// 読み込んだstlデータをobjフォーマットに変換し、出力する
/// 戻り値: =0 ファイルの書き込みが成功
///        !=0 書き込まなかった。または失敗した(std::ofstreamのエラーコード)
/// 事後条件:指定したパスにobjファイルが保存される.
int SaveObjFormatFromStl(std::ofstream& fout)const;

private:

/// 一度ファイルを読み、メンバに値が設定された状態で再びファイルを読むと.
/// メンバの値が上書きされる。これを避けるため、メンバの値を全てリセットする.
void Initialize();

/// Ascii形式のファイルを読み込み全ての座標値を取得する.
/// また図形のボックス枠を設定する.
/// 戻り値: =0ファイルの読み込が成功
///		  !=0 読み込みに失敗した(std::ifstreamのエラーコード)
/// 事前条件:既にオープンされたファイルストリームを渡す.
/// 事後条件:vecPosに全ての座標値が格納され、 boxにボックス枠の座標値が設定される.
/// ファイルストリームが進む.
int LoadAscii(
	std::ifstream& in, ///< [in/out]:既にオープンされたファイルストリーム .
	std::vector<MGPosition>& vecPos ///<[out]:ファイルから読み込んだ座標値の配列.
	);

/// Binary形式のファイルを読み込み全ての座標値を取得する.
/// また図形のボックス枠を設定する.
/// 戻り値: =0ファイルの読み込が成功
///		  !=0 読み込みに失敗した(std::ifstreamのエラーコード)
/// 事前条件:既にオープンされたファイルストリームを渡す.
/// 事後条件:vecPosに全ての座標値が格納され、 m_boxにボックス枠の座標値が設定される.
///			ファイルストリームが進む
int LoadBinary(
	std::ifstream& in, ///< [in/out]:既にオープンされたファイルストリーム.
	std::vector<MGPosition>& vecPos ///<[out]:ファイルから読み込んだ座標値の配列.
	);

/// ファイルから読み込んだ全ての頂点の座標値が格納されている.
/// vecPosから各三角形の法線を計算し、m_vecNormlTriangに追加.
/// vecPosの各要素の座標値の重複をトレランスを元に判断し
/// 重複を取り除き、m_vecPosに追加.
/// vecPosの各要素のインデックスをm_indicesに追加.
/// 事後条件:m_vecNormlTriang, m_indicesに値が設定される.
void set_mesh_data(
	const std::vector<MGPosition>& vecPos ///< [in]:ファイルから読み込んだ座標値の配列.
);

/// 入力であるpositionがVertexMapにすでに登録されているかチェックする.
/// 登録されていればそのm_vecPos添え字を返し,
/// 未登録の場合、新たにm_vecPosに格納し、positionとm_vecPos添え字のmapをVertexMapに
/// 登録する.
int IdentifyPosition(
	const MGPosition& position, ///< [in]:頂点の座標.
	triangleMap& VertexMap ///< [in/out]:頂点の座標、頂点のインデックスを保持する.
);

///mgTL2TrianglesをMGStlに追加する
void AddTL2Data(
	const mgTL2Triangles& tris,///<mgTL2Triangles whose data depend on tris.get_kind();
	triangleMap& VertexMap
);

///mgTL2TrianglesをMGStlに追加する
void AddTL2Data(
	double error,	///<Error to regard two points are the same.
	const mgTL2Triangles& tris,///<mgTL2Triangles whose data depend on tris.get_kind();
	triangleMap& VertexMap
);

/// STLファイルの図形を構成する各三角形の、座標の重複がない頂点座標の配列.
/// ファイルから読み込んだ順で、座標の重複を取り除き、座標値が格納されている.
std::vector<MGPosition> m_vecPos;

/// STLファイルの図形を構成する各三角形の法線ベクトルの配列.
/// ファイルから読み込まれた順で三角形の法線ベクトルが格納されている.
///m_vecNormlTriang.size()*3=m_indices.size().
///m_vecNormlTriang[i] is the normal of the triangle  m_indices[i], [i+1], [i+2] for
///i=0, ..., m_indices.size()/3.
std::vector<MGUnit_vector> m_vecNormlTriang;

/// STLファイルの図形を構成する各三角形の各頂点に対応する
/// 座標の配列のインデックスを格納する配列.
/// ファイルから読み込んだ順番で各三角形の頂点の並びが格納されている.
/// 各要素には該当する頂点座標の配列の添え字が格納されている.
/// 例：i番目の三角形の各頂点の座標の配列の添え字は
/// m_indices[i*3]、[i*3+1]、[i*3+2]という以上の順序で取得できる.
std::vector<int> m_indices;

};

/** @} */ // end of MGObjectRelated group

#endif
