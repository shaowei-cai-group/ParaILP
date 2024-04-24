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

// \sigma ax <= key
class LocalCon
{
public:
  long weight;
  Integer gap; // gap = \sigma a_i \dot x_i - key
  size_t posInUnsatConIdxs;
  Integer rhs;

  LocalCon();
  ~LocalCon();
};

class LocalConUtil
{
public:
  vector<LocalCon> conSet;
  vector<size_t> unsatConIdxs;
  vector<size_t> tempUnsatConIdxs;
  vector<size_t> tempSatConIdxs;
  vector<Integer> newConGap;
  vector<bool> isNewConGap;

  LocalConUtil();
  ~LocalConUtil();
  void Allocate(
      size_t conNum);
  LocalCon &GetCon(
      size_t idx);
  void insertUnsat(
      size_t conIdx);
  void RemoveUnsat(
      size_t conIdx);
};
