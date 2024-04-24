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
#include "Worker.h"

void cbkUpdateSol(void *_worker, const Integer obj)
{
  Worker *worker = (Worker *)_worker;
  Master *master = worker->master;
  if (worker->isFeasible == false || worker->best_obj > obj)
  {
    worker->best_obj = obj;
    worker->isFeasible = true;
    auto now = chrono::high_resolution_clock::now();
    worker->bestTime = chrono::duration_cast<chrono::milliseconds>(now - master->clk_st_global).count() / 1000.0;
    master->BestSolUpdate(obj);
  }
}

void cbkUpdatePopulation(void *_worker, const vector<LocalVar> &var, const Integer obj, const Config &config)
{
  Worker *worker = (Worker *)_worker;
  Master *master = worker->master;
  if (master->EnterPopulation(obj))
    master->PopulationUpdate(var, obj, config);
}

bool cbkGetSol(void *_worker)
{
  Worker *worker = (Worker *)_worker;
  Master *master = worker->master;
  return master->PopulationSharing(worker->solver);
}

Worker::Worker(int _tid, Master *_master) : tid(_tid), master(_master), best_obj(MaxValue)
{
  solver = new Local_ILP(master->modelConUtil, master->modelVarUtil);
  solver->worker = this;
  solver->cbkUpdateSol = cbkUpdateSol;
  solver->cbkUpdatePool = cbkUpdatePopulation;
  solver->cbkGetSol_Config = cbkGetSol;
  solver->terminated = 0;
  solver->isBin = master->modelVarUtil.isBin;
  solver->optValue = master->optValue;
  solver->mt.seed(tid == 0 ? 2832 : tid + OPT(nThreads) * OPT(seed));
  solver->sharingPeriod = master->sharingPeriod;
  solver->tid = tid;
  this->isFeasible = false;
  this->seed = tid + OPT(nThreads) * OPT(seed);
}

void Worker::Diversity()
{
  size_t con_num = solver->modelConUtil.conNum;
  if (tid == 0)
  {
    solver->sharingPeriod = -2;
    if (solver->weightUpperBound < con_num)
      solver->weightUpperBound = con_num;
    solver->objWeightUpperBound = solver->weightUpperBound / 10;
    return;
  }
  solver->sampleUnsat = 2 + rand_r(&seed) % 4;
  solver->bmsUnsat = 1000 + rand_r(&seed) % 2000;
  solver->sampleSat = 15 + rand_r(&seed) % 30;
  solver->bmsSat = 175 + rand_r(&seed) % 350;
  solver->bmsRandom = 30 + rand_r(&seed) % 60;
  solver->bmsPair = 38 + rand_r(&seed) % 77;
  solver->bmsRan = 75 + rand_r(&seed) % 150;
  solver->restartStep = 750000 + rand_r(&seed) % 1500000;
  // solver->weight_upb = 500 + rand_r(&seed) % 10000;
  if (solver->weightUpperBound < con_num)
    solver->weightUpperBound = con_num;
  solver->objWeightUpperBound = solver->weightUpperBound / ((rand_r(&seed) % 20) + 1);
  solver->greadyScore = (tid + 1) % 3 == 0 ? false : true;
  solver->satTightMove = (tid + 1) % 2 == 0 ? false : true;
  solver->pairTightMove = (tid + 1) % 5 == 0 ? false : true;
  solver->rvd = tid % 2 ? 1 : 0.5;
  return;
}

void Worker::InitSol()
{
  solver->Allocate();
  if (tid == 0)
  {
    LocalVarUtil &var_util = solver->localVarUtil;
    const ModelVarUtil &modelVarUtil = solver->modelVarUtil;
    for (int var_idx = 0; var_idx < var_util.varSet.size(); var_idx++)
    {
      LocalVar &var = var_util.GetVar(var_idx);
      const ModelVar &modelVar = modelVarUtil.GetVar(var_idx);
      if (modelVar.lowerBound > 0)
        var.nowValue = modelVar.lowerBound;
      else if (modelVar.upperBound < 0)
        var.nowValue = modelVar.upperBound;
      else
        var.nowValue = 0;
      if (!modelVar.inbound(var.nowValue))
        var.nowValue = (Integer)(modelVar.lowerBound / 2.0 + modelVar.upperBound / 2.0);
    }
    return;
  }
  PolarityInitSol();
}

void Worker::PolarityInitSol()
{
  LocalVarUtil &localVarUtil = solver->localVarUtil;
  const ModelVarUtil &modelVarUtil = solver->modelVarUtil;
  for (size_t varIdx = 0; varIdx < modelVarUtil.varNum; varIdx++)
  {
    auto &modelVar = modelVarUtil.GetVar(varIdx);
    LocalVar &localVar = localVarUtil.GetVar(varIdx);
    if ((rand_r(&seed) % 100 >= 50) && modelVar.inbound(0))
    {
      localVar.nowValue = 0;
      continue;
    }
    if (master->Polarity[varIdx] < 0)
      localVar.nowValue = modelVar.upperBound;
    else if (master->Polarity[varIdx] > 0)
      localVar.nowValue = modelVar.lowerBound;
    if (master->Polarity[varIdx] == 0 ||
        localVar.nowValue == InfiniteUpperBound ||
        localVar.nowValue == InfiniteLowerBound)
    {
      if (modelVar.lowerBound > 0)
        localVar.nowValue = modelVar.lowerBound;
      else if (modelVar.upperBound < 0)
        localVar.nowValue = modelVar.upperBound;
      else
        localVar.nowValue = 0;
      if (!modelVar.inbound(localVar.nowValue))
        localVar.nowValue = (Integer)(modelVar.lowerBound / 2.0 + modelVar.upperBound / 2.0);
    }
  }
}

int Worker::Solve()
{
  return solver->LocalSearch();
}

void Worker::Terminate()
{
  solver->terminated = 1;
}

void Worker::PrintResult()
{
  solver->PrintResult();
}

Worker::~Worker()
{
  delete solver;
  master = nullptr;
  bestSol.clear();
}
