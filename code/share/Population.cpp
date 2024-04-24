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
#include "Population.h"
#include "../Master.hpp"

Solution::Solution(const vector<LocalVar> &var, const Integer _obj, const Config &_config)
    : obj(_obj), config(_config)
{
  ori_obj = _obj;
  diff = 0;
  sol.reserve(var.size());
  for (auto &v : var)
    sol.push_back(v.nowValue);
}

Integer &Solution::operator[](int index)
{
  return sol[index];
}

bool Solution::operator==(const Solution &other) const
{
  if (other.ori_obj != ori_obj || sol.size() != other.sol.size())
    return false;
  for (int i = 0; i < sol.size(); i++)
    if (sol[i] != other.sol[i])
      return false;
  return true;
}

Population::Population()
{
  temp_diff.resize(PoolSize + 1, 0);
}

Solution &Population::operator[](int index)
{
  return solutions[index];
}

int Population::size()
{
  return solutions.size();
}

void Master::BestSolUpdate(const Integer obj)
{
  boost::mutex::scoped_lock lock(population.mtx_sol);
  if (bestObj <= obj)
    return;
  bestObj = obj;
  auto now = chrono::high_resolution_clock::now();
  double bestTime = chrono::duration_cast<chrono::milliseconds>(now - clk_st_global).count() / 1000.0;
  printf("n %-40s\t%-20f\n", itos(bestObj + modelVarUtil.objBias).c_str(), bestTime);
}

bool Master::EnterPopulation(const Integer obj)
{
  boost::mutex::scoped_lock lock(population.mtx_submit);
  if (submitObj.find(obj) != submitObj.end())
    return false;
  submitObj.insert(obj);
  // printf("c obj: %s\n", itos(obj).c_str());
  return true;
}

void Master::PopulationUpdate(const vector<LocalVar> &var, const Integer obj, const Config &config)
{
  boost::mutex::scoped_lock lock(population.mtx_pool);
  if (population.size() == 0)
  {
    population.solutions.emplace_back(var, obj, config);
    return;
  }
  population.solutions.emplace_back(var, obj, config);
  int pre_add = population.size() - 1;
  int new_idx = population.size() - 1;
  population.temp_diff[new_idx] = 0;
  for (int i = 0; i < pre_add; i++)
  {
    population.temp_diff[i] = population[i].diff;
    Integer new_diff = 0;
    for (int j = 0; j < var.size(); j++)
      new_diff += abs(long(population[i][j] - var[j].nowValue));
    population.temp_diff[i] += new_diff;
    population.temp_diff[new_idx] += new_diff;
  }
  if (population.size() <= PoolSize)
  {
    for (int i = 0; i < population.size(); i++)
      population[i].diff = population.temp_diff[i];
    return;
  }
  int worst_score, worst_index = -1;
  for (int i = 0; i < population.size(); ++i)
  {
    int obj_rank = 1, diff_rank = 1;
    for (int j = 0; j < population.size(); ++j)
    {
      if (i == j)
        continue;
      obj_rank += (population[i].obj < population[j].obj);
      diff_rank += (population.temp_diff[i] > population.temp_diff[j]);
    }
    population[i].score = diff_rank * 20 + obj_rank * 80; // max is better
    if (worst_index == -1 || population[i].score < worst_score)
    {
      worst_index = i;
      worst_score = population[i].score;
    }
  }
  if (worst_index != new_idx)
  {
    for (int i = 0; i < population.size(); i++)
    {
      if (i == worst_index)
        continue;
      Integer del_diff = 0;
      for (int j = 0; j < var.size(); j++)
        del_diff += abs(long(population[i][j] - population[worst_index][j]));
      population.temp_diff[i] -= del_diff;
      population[i].diff = population.temp_diff[i];
    }
    population[worst_index].ori_obj = obj;
    population[worst_index].obj = obj;
    population[worst_index].config = config;
    population[worst_index].diff = population[new_idx].diff;
    for (size_t idx = 0; idx < var.size(); idx++)
      population[worst_index][idx] = var[idx].nowValue;
  }
  population.solutions.pop_back();
  return;
}

bool Master::PopulationSharing(Local_ILP *solver)
{
  boost::mutex::scoped_lock lock(population.mtx_pool);
  // boost::shared_lock<boost::shared_mutex> lock(pool.mtx_pool);
  if (population.size() == 0)
    return false;
  size_t sol_idx_1 = solver->mt() % population.size();
  size_t config_idx = solver->mt() % population.size();
  const auto &sol_1 = population[sol_idx_1].sol;
  double punish_coeff = 1 + Punish;
  if (population[sol_idx_1].obj < 0)
    punish_coeff = 1 - Punish;
  population[sol_idx_1].obj *= punish_coeff;
  if (solver->mt() % 100 < 50)
  {
    for (size_t idx = 0; idx < sol_1.size(); idx++)
      solver->localVarUtil.varSet[idx].nowValue = sol_1[idx];
    if (config_idx == sol_idx_1)
      if (config_idx + 1 < population.size())
        config_idx++;
      else if (config_idx > 0)
        config_idx--;
  }
  else
  {
    size_t sol_idx_2 = solver->mt() % population.size();
    if (sol_idx_2 == sol_idx_1)
      if (sol_idx_2 + 1 < population.size())
        sol_idx_2++;
      else if (sol_idx_2 > 0)
        sol_idx_2--;
    const auto &sol_2 = population[sol_idx_2].sol;
    if (sol_idx_2 == sol_idx_1)
      for (size_t idx = 0; idx < sol_1.size(); idx++)
      {
        auto &var = solver->localVarUtil.GetVar(idx);
        const auto &common_var = modelVarUtil.GetVar(idx);
        if (solver->mt() % 100 < 25)
          var.nowValue = common_var.lowerBound + solver->mt() % (common_var.upperBound + 1 - common_var.lowerBound);
      }
    else
    {
      for (size_t idx = 0; idx < sol_1.size(); idx++)
      {
        auto &var = solver->localVarUtil.GetVar(idx);
        const auto &common_var = modelVarUtil.GetVar(idx);
        var.nowValue = sol_1[idx];
        if (var.nowValue != sol_2[idx])
        {
          if (solver->mt() % 100 < 50)
            var.nowValue = sol_2[idx];
          if (solver->mt() % 100 < 25)
            var.nowValue = common_var.lowerBound + solver->mt() % (common_var.upperBound + 1 - common_var.lowerBound);
        }
      }
      population[sol_idx_2].obj *= punish_coeff;
    }
  }
  const auto &_config = population[config_idx].config;
  if (solver->mt() % 100 < 75)
    solver->sampleUnsat = _config.sampleUnsat;
  else
    solver->sampleUnsat = 2 + solver->mt() % 4;
  if (solver->mt() % 100 < 75)
    solver->bmsUnsat = _config.bmsUnsat;
  else
    solver->bmsUnsat = 1000 + solver->mt() % 2000;
  if (solver->mt() % 100 < 75)
    solver->sampleSat = _config.sampleSat;
  else
    solver->sampleSat = 15 + solver->mt() % 30;
  if (solver->mt() % 100 < 75)
    solver->bmsSat = _config.bmsSat;
  else
    solver->bmsSat = 175 + solver->mt() % 350;
  if (solver->mt() % 100 < 75)
    solver->bmsRandom = _config.bmsRandom;
  else
    solver->bmsRandom = 30 + solver->mt() % 60;
  if (solver->mt() % 100 < 75)
    solver->bmsPair = _config.bmsPair;
  else
    solver->bmsPair = 38 + solver->mt() % 77;
  if (solver->mt() % 100 < 75)
    solver->bmsRan = _config.bmsRan;
  else
    solver->bmsRan = 75 + solver->mt() % 150;
  if (solver->mt() % 100 < 75)
    solver->weightUpperBound = _config.weight_upb;
  if (solver->mt() % 100 < 75)
    solver->objWeightUpperBound = _config.obj_weight_upb;
  else
    solver->objWeightUpperBound = solver->weightUpperBound / ((solver->mt() % 20) + 1);
  if (solver->mt() % 100 < 75)
    solver->greadyScore = _config.greadyScore;
  else
    solver->greadyScore = solver->mt() % 2;
  if (solver->mt() % 100 < 75)
    solver->satTightMove = _config.satTightMove;
  else
    solver->satTightMove = solver->mt() % 2;
  if (solver->mt() % 100 < 75)
    solver->pairTightMove = _config.pairTightMove;
  else
    solver->pairTightMove = solver->mt() % 2;
  if (solver->mt() % 100 < 75)
    solver->rvd = _config.rvd;
  else
    solver->rvd = solver->mt() % 2 ? 1 : 0.5;
  return true;
}
