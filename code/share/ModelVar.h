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
#pragma once
#include "../utils/header.h"

class ModelVar
{
public:
  string name;
  size_t idx;
  Integer upperBound; // [lower_bound, upper_bound]
  Integer lowerBound;
  size_t termNum;
  VarType type;
  vector<size_t> conIdxs;
  vector<size_t> posInCon; // position in cons[con_idxs[i]]

  ModelVar(
      const string &_name,
      size_t _idx);
  ~ModelVar();
  void SetType(
      VarType varType);
  bool inbound(
      Integer value) const;
};

class ModelVarUtil
{
public:
  unordered_map<string, size_t> name2idx;
  vector<ModelVar> varSet;
  vector<size_t> varIdx2ObjIdx;
  Integer objBias;
  size_t varNum;
  size_t integerNum;
  size_t binaryNum;
  size_t fixedNum;
  bool isBin;

  ModelVarUtil();
  ~ModelVarUtil();
  size_t MakeVar(
      const string &name);
  const ModelVar &GetVar(
      size_t idx) const;
  ModelVar &GetVar(
      size_t idx);
  ModelVar &GetVar(
      const string &name);
};