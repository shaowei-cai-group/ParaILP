/*=====================================================================================

    Filename:     LocalCon.cpp

    Description:
        Version:  1.0

    Author:       Peng Lin, penglincs@outlook.com

    Organization: Shaowei Cai Group,
                  State Key Laboratory of Computer Science,
                  Institute of Software, Chinese Academy of Sciences,
                  Beijing, China

=====================================================================================*/
#include "ModelCon.h"

ModelCon::ModelCon(
    const string &_name,
    size_t _idx)
    : name(_name),
      idx(_idx),
      isEqual(false),
      isLess(false),
      rhs(0),
      inferSAT(false),
      termNum(-1)

{
  // coeffSet.reserve(55615);
  // varIdxSet.reserve(55615);
}

ModelCon::~ModelCon()
{
  coeffSet.clear();
  varIdxSet.clear();
}

ModelConUtil::ModelConUtil()
    : conNum(-1)
{
  // name2idx.reserve(17602870);
}

ModelConUtil::~ModelConUtil()
{
  conSet.clear();
  name2idx.clear();
}

size_t ModelConUtil::MakeCon(
    const string &name)
{
  auto iter = name2idx.find(name);
  if (iter != name2idx.end())
    return iter->second;
  int con_idx = conSet.size();
  conSet.emplace_back(name, con_idx);
  name2idx[name] = con_idx;
  return con_idx;
}

size_t ModelConUtil::GetConIdx(
    const string &name)
{
  if (name == objName)
    return 0;
  auto iter = name2idx.find(name);
  return iter->second;
}

const ModelCon &ModelConUtil::GetCon(
    size_t idx) const
{
  return conSet[idx];
}

ModelCon &ModelConUtil::GetCon(
    size_t idx)
{
  return conSet[idx];
}

ModelCon &ModelConUtil::GetCon(
    const string &name)
{
  if (name == objName)
    return conSet[0];
  return conSet[name2idx[name]];
}
