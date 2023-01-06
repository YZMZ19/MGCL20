/********************************************************************/
/* Copyright (c) 2019 System fugen G.K. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/********************************************************************/
#include "StdAfx.h"
#include "mg/GelFactory.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

namespace
{
	bool g_valid = false; // クラスオブジェクトが有効かどうか
}

// コンストラクター
MGGelFactoryRegistry::MGGelFactoryRegistry()
{
	g_valid = true;
}

// デストラクター
MGGelFactoryRegistry::~MGGelFactoryRegistry()
{
	g_valid = false;
}

// シングルトンアクセス。
MGGelFactoryRegistry* MGGelFactoryRegistry::get_instance()
{
	static MGGelFactoryRegistry registry;
	return &registry;
}

// クラスオブジェクトが有効かどうかを返す。
bool MGGelFactoryRegistry::is_valid()
{
	return g_valid;
}

// 名前からコマンドオブジェクトを作成する。
MGGel* MGGelFactoryRegistry::create_gel(const KeyType& name) const
{
	TypeMap::const_iterator it = m_map.find(name);
	if(it == m_map.end()){
		//throw UnknownGelError(name);
		return 0;
	}

	return it->second->create_gel();
}

// ファクトリーオブジェクトを登録する。
void MGGelFactoryRegistry::register_factory(
	const KeyType& name,
	MGGelFactoryBase* factory)
{
	m_map.insert(std::make_pair(name, factory));
}

// ファクトリーオブジェクトを登録から削除する。
void MGGelFactoryRegistry::unregister_factory(const KeyType& name)
{
	m_map.erase(name);
}

// ファクトリーオブジェクトを登録から削除する。
void MGGelFactoryRegistry::unregister_factory(MGGelFactoryBase* factory)
{
	TypeMap::iterator first = m_map.begin(), last = m_map.end();
	for(; first != last;){
		if(first->second == factory){
			first = m_map.erase(first);
		}
		else{
			++first;
		}
	}
}
