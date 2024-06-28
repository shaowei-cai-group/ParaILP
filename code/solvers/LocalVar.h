/*=====================================================================================

    Filename:     LocalVar.h

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

class LocalVar
{
public:
  Integer nowValue;
  Integer bestValue;
  long allowIncStep;
  long allowDecStep;
  long lastIncStep;
  long lastDecStep;

  LocalVar();
  ~LocalVar();
};

class LocalVarUtil
{
public:
  vector<LocalVar> varSet;
  vector<Integer> lowerDeltaInLiftMove;
  vector<Integer> upperDeltaInLifiMove;
  vector<Integer> tempDeltas;
  vector<size_t> tempVarIdxs;
  vector<bool> scoreTable;
  unordered_set<size_t> affectedVar;
  vector<size_t> tempTwoVarIdxs_1;
  vector<Integer> tempTwoDeltas_1;
  vector<size_t> tempTwoVarIdxs_2;
  vector<Integer> tempTwoDeltas_2;

  LocalVarUtil();
  ~LocalVarUtil();
  void Allocate(
      size_t varNum,
      size_t varNumInObj);
  LocalVar &GetVar(
      size_t idx);
};