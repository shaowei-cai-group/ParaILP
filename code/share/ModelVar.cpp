/*=====================================================================================

    Filename:     ModelVar.cpp

    Description:
        Version:  1.0

    Author:       Peng Lin, penglincs@outlook.com

    Organization: Shaowei Cai Group,
                  State Key Laboratory of Computer Science,
                  Institute of Software, Chinese Academy of Sciences,
                  Beijing, China

=====================================================================================*/
#include "ModelVar.h"

ModelVar::ModelVar(const string &_name, size_t _idx)
    : name(_name),
      idx(_idx),
      upperBound(DefaultUpperBound),
      lowerBound(DefaultLowerBound),
      type(VarType::Binary),
      termNum(-1)
{
  // con_idxs.reserve(5000);
  // pos_in_con.reserve(5000);
}

ModelVar::~ModelVar()
{
  conIdxs.clear();
  posInCon.clear();
}

void ModelVar::SetType(VarType varType)
{
  type = varType;
}

bool ModelVar::inbound(Integer v) const
{
  return lowerBound <= v && v <= upperBound;
}

ModelVarUtil::ModelVarUtil()
    : objBias(0),
      varNum(-1),
      integerNum(0),
      binaryNum(0),
      fixedNum(0),
      isBin(true)
{
  // name2idx.reserve(37709950);
}
ModelVarUtil::~ModelVarUtil()
{
  varIdx2ObjIdx.clear();
  name2idx.clear();
  varSet.clear();
}

size_t ModelVarUtil::MakeVar(const string &name)
{
  auto iter = name2idx.find(name);
  if (iter != name2idx.end())
    return iter->second;
  size_t var_idx = varSet.size();
  varSet.emplace_back(name, var_idx);
  name2idx[name] = var_idx;
  return var_idx;
}

const ModelVar &ModelVarUtil::GetVar(size_t idx) const
{
  return varSet[idx];
}

ModelVar &ModelVarUtil::GetVar(size_t idx)
{
  return varSet[idx];
}

ModelVar &ModelVarUtil::GetVar(const string &name)
{
  auto iter = name2idx.find(name);
  return varSet[name2idx[name]];
}