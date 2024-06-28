/*=====================================================================================

    Filename:     Worker.h

    Description:
        Version:  1.0

    Author:       Peng Lin, penglincs@outlook.com

    Organization: Shaowei Cai Group,
                  State Key Laboratory of Computer Science,
                  Institute of Software, Chinese Academy of Sciences,
                  Beijing, China

=====================================================================================*/
#pragma once
#include <chrono>
#include "Master.hpp"
#include "solvers/local_ilp.h"
#include "utils/header.h"

class Master;
class LocalVarUtil;
class ModelVarUtil;
class LocalVar;
class LocalCon;
class modelVar;
class modelCon;

class Worker
{
public:
  int tid;
  unsigned int seed;
  Integer best_obj;
  double bestTime;
  Master *master;
  vector<Integer> bestSol;
  bool isFeasible;
  Local_ILP *solver;

  void Diversity();
  int Solve();
  void Terminate();
  void InitSol();
  void ShuffleInitSol();
  void PolarityInitSol();
  void PrintResult();
  Worker(int _tid, Master *_master);
  ~Worker();
  friend void cbkUpdateSol(
      void *_worker,
      const Integer obj);
  friend void cbkUpdatePopulation(
      void *_worker,
      const vector<LocalVar> &var,
      const Integer obj,
      const Config &config);
  friend bool cbkGetSol(
      void *_worker);
};