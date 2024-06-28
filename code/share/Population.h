/*=====================================================================================

    Filename:     Population.h

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
#include "solvers/LocalVar.h"

class Config
{
public:
  size_t sampleUnsat;
  size_t bmsUnsat;
  size_t sampleSat;
  size_t bmsSat;
  size_t bmsRandom;
  size_t bmsPair;
  size_t bmsRan;
  size_t weight_upb;
  size_t obj_weight_upb;
  bool greadyScore;
  bool satTightMove;
  bool pairTightMove;
  double rvd;
  Config(const Config &_config)
  {
    sampleUnsat = _config.sampleUnsat;
    bmsUnsat = _config.bmsUnsat;
    sampleSat = _config.sampleSat;
    bmsSat = _config.bmsSat;
    bmsRandom = _config.bmsRandom;
    bmsPair = _config.bmsPair;
    bmsRan = _config.bmsRan;
    weight_upb = _config.weight_upb;
    obj_weight_upb = _config.obj_weight_upb;
    greadyScore = _config.greadyScore;
    satTightMove = _config.satTightMove;
    pairTightMove = _config.pairTightMove;
    rvd = _config.rvd;
  }
  Config(size_t _sample_unsat, size_t _bms_unsat, size_t _sample_sat, size_t _bms_sat,
         size_t _sample_pair, size_t _bms_pair, size_t _bms_ran, size_t _weight_upb, size_t _obj_weight_upb,
         bool _gready_score, bool _sat_tight_move, bool _pair_tight_move, double _rvd)
      : sampleUnsat(_sample_unsat), bmsUnsat(_bms_unsat), sampleSat(_sample_sat),
        bmsSat(_bms_sat), bmsRandom(_sample_pair), bmsPair(_bms_pair), bmsRan(_bms_ran),
        weight_upb(_weight_upb), obj_weight_upb(_obj_weight_upb), greadyScore(_gready_score),
        satTightMove(_sat_tight_move), pairTightMove(_pair_tight_move), rvd(_rvd) {}
  const string tostring() const
  {
    string s;
    s += to_string(sampleUnsat) + "\t";
    s += to_string(bmsUnsat) + "\t";
    s += to_string(sampleSat) + "\t";
    s += to_string(bmsSat) + "\t";
    s += to_string(bmsRandom) + "\t";
    s += to_string(bmsPair) + "\t";
    s += to_string(bmsRan) + "\t";
    s += to_string(weight_upb) + "\t";
    s += to_string(obj_weight_upb) + "\t";
    s += to_string(greadyScore) + "\t";
    s += to_string(satTightMove) + "\t";
    s += to_string(pairTightMove) + "\t";
    s += to_string(rvd) + "\t";
    return s;
  }
  bool operator==(const Config &_config) const
  {
    return sampleUnsat == _config.sampleUnsat &&
           bmsUnsat == _config.bmsUnsat &&
           sampleSat == _config.sampleSat &&
           bmsSat == _config.bmsSat &&
           bmsRandom == _config.bmsRandom &&
           bmsPair == _config.bmsPair &&
           bmsRan == _config.bmsRan &&
           weight_upb == _config.weight_upb &&
           obj_weight_upb == _config.obj_weight_upb &&
           greadyScore == _config.greadyScore &&
           satTightMove == _config.satTightMove &&
           pairTightMove == _config.pairTightMove &&
           rvd == _config.rvd;
  }
};

class Solution
{
public:
  vector<Integer> sol;
  Integer ori_obj;
  Integer obj;
  Integer diff;
  int score;
  Config config;

  Solution(const vector<LocalVar> &var, const Integer _obj, const Config &_config);

  Integer &operator[](int index);

  bool operator==(const Solution &other) const;
};

class Population
{
public:
  vector<Solution> solutions;
  vector<Integer> temp_diff;
  mutable boost::mutex mtx_pool;
  mutable boost::mutex mtx_sol;
  mutable boost::mutex mtx_submit;
  Population();
  // mutable boost::shared_mutex mtx_pool;
  Solution &operator[](int index);

  int size();
};