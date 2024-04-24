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
#include "LocalVar.h"

LocalVar::LocalVar()
    : allowIncStep(0),
      allowDecStep(0),
      lastIncStep(0),
      lastDecStep(0)
{
}

LocalVar::~LocalVar()
{
}

LocalVarUtil::LocalVarUtil()
{
}

void LocalVarUtil::Allocate(
    size_t var_number,
    size_t var_number_of_obj)
{
  tempDeltas.reserve(var_number);
  tempVarIdxs.reserve(var_number);
  tempTwoVarIdxs_1.reserve(var_number);
  tempTwoDeltas_1.reserve(var_number);
  tempTwoVarIdxs_2.reserve(var_number);
  tempTwoDeltas_2.reserve(var_number);
  affectedVar.reserve(var_number);
  varSet.resize(var_number);
  scoreTable.resize(var_number, false);
  lowerDeltaInLiftMove.resize(var_number_of_obj);
  upperDeltaInLifiMove.resize(var_number_of_obj);
}

LocalVarUtil::~LocalVarUtil()
{
  lowerDeltaInLiftMove.clear();
  upperDeltaInLifiMove.clear();
  scoreTable.clear();
  affectedVar.clear();
  tempTwoVarIdxs_1.clear();
  tempTwoDeltas_1.clear();
  tempTwoVarIdxs_2.clear();
  tempTwoDeltas_2.clear();
  varSet.clear();
  tempDeltas.clear();
  tempVarIdxs.clear();
}

LocalVar &LocalVarUtil::GetVar(
    size_t idx)
{
  return varSet[idx];
}