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
#include "local_ilp.h"

// obj min. \sigma ax
// \sigma ax - key <= 0, where key is cur best obj + epsilon

int Local_ILP::LocalSearch()
{
  InitState();
  auto &obj = localConUtil.conSet[0];
  curStep = 0;
  while (true)
  {
    // printf("c %-10d unsat: %-10d\n", curStep, localConUtil.unsatConIdxs.size());
    if (localConUtil.unsatConIdxs.empty())
    {
      if (!isFoundFeasible || obj.gap <= 0)
      {
        UpdateBestSolution();
        isFoundFeasible = true;
        cbkUpdateSol(worker, obj.rhs + 1);
      }
      LiftMove();
      if (GetObjValue() == optValue)
        return 1;
      ++curStep;
      if (terminated)
        break;
      continue;
    }
    if (terminated)
      break;
    if (sharingPeriod == -2)
    {
      if (curStep - lastImproveStep > restartStep)
      {
        Restart();
        continue;
      }
    }
    else if (curStep - lastImproveStep > sharingPeriod)
    {
      Sharing();
      continue;
    }

    if (!UnsatTightMove())
    {
      if (mt() % 10000 > smoothProbability)
        UpdateWeight();
      else
        SmoothWeight();
      RandomTightMove();
    }
    ++curStep;
  }
  return 0;
}

Integer Local_ILP::GetObjValue()
{
  return localConUtil.conSet[0].rhs + 1 + modelVarUtil.objBias;
}

void Local_ILP::InitState()
{
  for (size_t con_idx = 1; con_idx < localConUtil.conSet.size(); ++con_idx)
  {
    auto &con = localConUtil.conSet[con_idx];
    auto &common_con = modelConUtil.conSet[con_idx];
    // calculate gap
    con.gap = 0;
    for (size_t i = 0; i < common_con.termNum; ++i)
    {
      con.gap += common_con.coeffSet[i] * localVarUtil.GetVar(common_con.varIdxSet[i]).nowValue;
    }
    con.gap -= con.rhs;
    if (con.gap > 0)
      localConUtil.insertUnsat(con_idx);
  }
  // Obj
  auto &obj = localConUtil.conSet[0];
  auto &modelObj = modelConUtil.conSet[0];
  isFoundFeasible = false;
  obj.rhs = 1e38;
  obj.gap = 0;
  for (size_t i = 0; i < modelObj.termNum; ++i)
    obj.gap += modelObj.coeffSet[i] * localVarUtil.GetVar(modelObj.varIdxSet[i]).nowValue;
  obj.gap -= obj.rhs;
}

void Local_ILP::UpdateBestSolution()
{
  lastImproveStep = curStep;
  for (auto &var : localVarUtil.varSet)
  {
    var.bestValue = var.nowValue;
  }
  auto &obj = localConUtil.conSet[0];
  auto &modelObj = modelConUtil.conSet[0];
  obj.rhs = 0;
  for (size_t i = 0; i < modelObj.termNum; ++i)
  {
    size_t var_idx = modelObj.varIdxSet[i];
    auto &var = localVarUtil.GetVar(var_idx);
    Integer coeff = modelObj.coeffSet[i];
    obj.rhs += var.nowValue * coeff;
  }

  obj.rhs -= 1;
  obj.gap = 1;
}

void Local_ILP::makePair(vector<size_t> &var_idxs, vector<Integer> &deltas,
                         vector<size_t> &var_idxs_1, vector<Integer> &deltas_1,
                         vector<size_t> &var_idxs_2, vector<Integer> &deltas_2)
{
  size_t bms_list_size = var_idxs.size();
  if (var_idxs.size() > bmsRandom)
  {
    bms_list_size = bmsRandom;
    for (size_t i = 0; i < bmsRandom; ++i)
    {
      size_t ran_idx = mt() % (var_idxs.size() - i);
      size_t var_idx = var_idxs[ran_idx + i];
      Integer delta = deltas[ran_idx + i];
      var_idxs[ran_idx + i] = var_idxs[i];
      deltas[ran_idx + i] = deltas[i];
      var_idxs[i] = var_idx;
      deltas[i] = delta;
    }
  }
  for (int q = 1; q <= 2; ++q)
  {
    if (q == 2 && var_idxs_1.size() > 0)
      break;
    for (size_t i = 0; i < bms_list_size; ++i)
    {
      size_t var_idx = var_idxs[i];
      Integer delta = deltas[i];
      auto &var = localVarUtil.GetVar(var_idx);
      auto &common_var = modelVarUtil.GetVar(var_idx);
      for (size_t j = 0; j < common_var.conIdxs.size(); ++j)
      {
        size_t con_idx = common_var.conIdxs[j];
        size_t pos = common_var.posInCon[j];
        auto &con = localConUtil.conSet[con_idx];
        auto &common_con = modelConUtil.conSet[con_idx];
        if (con_idx == 0)
          continue;
        Integer new_gap = con.gap + common_con.coeffSet[pos] * delta;
        if (q == 1 && new_gap > 0 && con.gap == 0 || q == 2 && new_gap > 0 && con.gap <= 0)
          for (size_t idx = 0; idx < common_con.termNum; ++idx)
          {
            size_t rel_var_idx = common_con.varIdxSet[idx];
            if (rel_var_idx == var_idx)
              continue;
            double temp_delta = -((double)new_gap / common_con.coeffSet[idx]);
            Integer rel_delta;
            if (common_con.coeffSet[idx] > 0)
            {
              rel_delta = floor(temp_delta);
            }
            else
            {
              rel_delta = ceil(temp_delta);
            }
            auto &rel_var = localVarUtil.GetVar(rel_var_idx);
            auto &common_rel_var = modelVarUtil.GetVar(rel_var_idx);
            if (rel_delta < 0 && curStep < rel_var.allowDecStep ||
                rel_delta > 0 && curStep < rel_var.allowIncStep)
            { // tabu
              continue;
            }
            if (!common_rel_var.inbound(rel_var.nowValue + rel_delta))
            {
              if (rel_delta < 0)
                rel_delta = common_rel_var.lowerBound - var.nowValue;
              else
                rel_delta = common_rel_var.upperBound - var.nowValue;
            }
            if (rel_delta != 0 && common_rel_var.inbound(rel_var.nowValue + rel_delta))
            {
              var_idxs_1.push_back(var_idx);
              deltas_1.push_back(delta);
              var_idxs_2.push_back(rel_var_idx);
              deltas_2.push_back(rel_delta);
            }
          }
      }
    }
  }
}
double Local_ILP::towScore(const ModelVar &var_1, Integer delta_1, const ModelVar &var_2, Integer delta_2)
{
  double score = 0;
  auto &new_gap_arr = localConUtil.newConGap;
  auto &is_new_gap_arr = localConUtil.isNewConGap;
  vector<size_t> changed_idx;
  // auto &obj = con_util.cons[0];
  for (size_t i = 0; i < var_1.conIdxs.size(); ++i)
  {
    size_t con_idx = var_1.conIdxs[i];
    auto &con = localConUtil.conSet[con_idx];
    auto &common_con = modelConUtil.conSet[con_idx];
    size_t pos = var_1.posInCon[i];

    if (con_idx == 0)
    { // obj
      if (isFoundFeasible)
      {
        Integer new_gap = con.gap + common_con.coeffSet[pos] * delta_1;
        if (1 > new_gap)
        { // batter
          score += con.weight;
        }
        else if (greadyScore || 1 < new_gap)
        { // worse
          score -= con.weight;
        }
        new_gap_arr[con_idx] = new_gap;
        is_new_gap_arr[con_idx] = true;
        changed_idx.push_back(con_idx);
      }
      continue;
    }

    // TODO break ties by obj
    // TODO obj score

    Integer new_gap = con.gap + common_con.coeffSet[pos] * delta_1;

    bool is_pre_sat = con.gap <= 0;
    bool is_now_sat = new_gap <= 0;

    // TODO distance
    if (!is_pre_sat && is_now_sat)
    {
      score += con.weight;
    }
    else if (is_pre_sat && !is_now_sat)
    {
      score -= con.weight;
    }
    else if (!is_pre_sat && !is_now_sat)
    {
      if (con.gap > new_gap)
      {                            // better
        score += con.weight * rvd; //*(1-((double)new_gap/(double)con.gap));
      }
      else if (greadyScore || con.gap < new_gap)
      {                            // worse
        score -= con.weight * rvd; //*(1-((double)con.gap/(double)new_gap));
      }
    }
    new_gap_arr[con_idx] = new_gap;
    is_new_gap_arr[con_idx] = true;
    changed_idx.push_back(con_idx);
  }

  for (size_t i = 0; i < var_2.conIdxs.size(); ++i)
  {
    size_t con_idx = var_2.conIdxs[i];
    auto &con = localConUtil.conSet[con_idx];
    auto &common_con = modelConUtil.conSet[con_idx];
    size_t pos = var_2.posInCon[i];
    if (is_new_gap_arr[con_idx])
    {
      if (con_idx == 0)
      { // obj
        if (isFoundFeasible)
        {
          auto gap_1 = new_gap_arr[con_idx];
          Integer gap_2 = gap_1 + common_con.coeffSet[pos] * delta_2;
          if (gap_1 < 1)
            score -= con.weight;
          else if (greadyScore || 1 < gap_1)
            score += con.weight;
          if (gap_2 < 1)
            score += con.weight;
          else if (greadyScore || 1 < gap_2)
            score -= con.weight;
          // if (1 <= gap_1 && 1 > gap_2)
          // { // batter
          //   score += con.weight << 1;
          // }
          // else if (1 > gap_1 && 1 <= gap_2)
          // {
          //   score -= con.weight << 1;
          // }
        }
        continue;
      }

      // TODO break ties by obj
      // TODO obj score

      auto gap_1 = new_gap_arr[con_idx];
      Integer gap_2 = gap_1 + common_con.coeffSet[pos] * delta_2;

      bool is_pre_sat = con.gap <= 0;
      bool is_now_sat_1 = gap_1 <= 0;
      bool is_now_sat_2 = gap_2 <= 0;

      // TODO distance
      if (is_pre_sat && is_now_sat_1 && is_now_sat_2)
      {
        continue;
      }
      else if (is_pre_sat && is_now_sat_1 && !is_now_sat_2)
      {
        score -= con.weight;
      }
      else if (is_pre_sat && !is_now_sat_1 && is_now_sat_2)
      {
        score += con.weight;
      }
      else if (is_pre_sat && !is_now_sat_1 && !is_now_sat_2)
      {
        continue;
      }
      else if (!is_pre_sat && is_now_sat_1 && is_now_sat_2)
      {
        continue;
      }
      else if (!is_pre_sat && is_now_sat_1 && !is_now_sat_2)
      {
        score -= con.weight;
        if (gap_2 < con.gap)
          score += con.weight * rvd;
        else if (greadyScore || con.gap < gap_2)
          score -= con.weight * rvd;
      }
      else if (!is_pre_sat && !is_now_sat_1 && is_now_sat_2)
      {
        if (gap_1 < con.gap)
          score -= con.weight * rvd;
        else if (greadyScore || con.gap < gap_1)
          score += con.weight * rvd;
        score += con.weight;
      }
      else if (!is_pre_sat && !is_now_sat_1 && !is_now_sat_2)
      {
        if (gap_1 < con.gap)
          score -= con.weight * rvd;
        else if (greadyScore || con.gap < gap_1)
          score += con.weight * rvd;
        if (gap_2 < con.gap)
          score += con.weight * rvd;
        else if (greadyScore || con.gap < gap_2)
          score -= con.weight * rvd;
      }
    }
    else
    {
      if (con_idx == 0)
      { // obj
        if (isFoundFeasible)
        {
          Integer new_gap = con.gap + common_con.coeffSet[pos] * delta_2;
          if (1 > new_gap)
          { // batter
            score += con.weight;
          }
          else if (greadyScore || 1 < new_gap)
          { // worse
            score -= con.weight;
          }
        }
        continue;
      }

      // TODO break ties by obj
      // TODO obj score

      Integer new_gap = con.gap + common_con.coeffSet[pos] * delta_2;

      bool is_pre_sat = con.gap <= 0;
      bool is_now_sat = new_gap <= 0;

      // TODO distance
      if (!is_pre_sat && is_now_sat)
      {
        score += con.weight;
      }
      else if (is_pre_sat && !is_now_sat)
      {
        score -= con.weight;
      }
      else if (!is_pre_sat && !is_now_sat)
      {
        if (con.gap > new_gap)
        {                            // better
          score += con.weight * rvd; //*(1-((double)new_gap/(double)con.gap));
        }
        else if (greadyScore || con.gap < new_gap)
        {                            // worse
          score -= con.weight * rvd; //*(1-((double)con.gap/(double)new_gap));
        }
      }
    }
  }
  for (auto idx : changed_idx)
    is_new_gap_arr[idx] = false;
  return score;
}
bool Local_ILP::towMove(vector<size_t> &var_idxs, vector<Integer> &deltas)
{
  double best_score = 0;
  long best_last_move_step = std::numeric_limits<long>::max();
  size_t best_var_idx_1 = -1; // maximum size_t value
  Integer best_delta_1 = 0;
  size_t best_var_idx_2 = -1; // maximum size_t value
  Integer best_delta_2 = 0;
  vector<size_t> &var_idxs_1 = localVarUtil.tempTwoVarIdxs_1;
  vector<Integer> &deltas_1 = localVarUtil.tempTwoDeltas_1;
  vector<size_t> &var_idxs_2 = localVarUtil.tempTwoVarIdxs_2;
  vector<Integer> &deltas_2 = localVarUtil.tempTwoDeltas_2;
  var_idxs_1.clear();
  deltas_1.clear();
  var_idxs_2.clear();
  deltas_2.clear();
  makePair(var_idxs, deltas, var_idxs_1, deltas_1, var_idxs_2, deltas_2);
  size_t bms_list_size = var_idxs_1.size();
  if (var_idxs_1.size() > bmsPair)
  {
    bms_list_size = bmsPair;
    for (size_t i = 0; i < bmsPair; ++i)
    {
      size_t ran_idx = mt() % (var_idxs_1.size() - i);
      size_t var_idx_1 = var_idxs_1[ran_idx + i];
      Integer delta_1 = deltas_1[ran_idx + i];
      size_t var_idx_2 = var_idxs_2[ran_idx + i];
      Integer delta_2 = deltas_2[ran_idx + i];
      var_idxs_1[ran_idx + i] = var_idxs_1[i];
      deltas_1[ran_idx + i] = deltas_1[i];
      var_idxs_2[ran_idx + i] = var_idxs_2[i];
      deltas_2[ran_idx + i] = deltas_2[i];
      var_idxs_1[i] = var_idx_1;
      deltas_1[i] = delta_1;
      var_idxs_2[i] = var_idx_2;
      deltas_2[i] = delta_2;
    }
  }

  for (size_t i = 0; i < bms_list_size; ++i)
  {
    size_t var_idx_1 = var_idxs_1[i];
    Integer delta_1 = deltas_1[i];
    size_t var_idx_2 = var_idxs_2[i];
    Integer delta_2 = deltas_2[i];
    auto &var_1 = localVarUtil.GetVar(var_idx_1);
    auto &common_var_1 = modelVarUtil.GetVar(var_idx_1);
    auto &var_2 = localVarUtil.GetVar(var_idx_2);
    auto &common_var_2 = modelVarUtil.GetVar(var_idx_2);
    double score = towScore(common_var_1, delta_1, common_var_2, delta_2);
    long last_move_step_1 = delta_1 < 0 ? var_1.lastDecStep : var_1.lastIncStep;
    long last_move_step_2 = delta_2 < 0 ? var_2.lastDecStep : var_2.lastIncStep;
    long last_move_step = last_move_step_1 + last_move_step_2;
    if (best_score < score ||
        best_score == score && last_move_step < best_last_move_step)
    {
      best_score = score;
      best_var_idx_1 = var_idx_1;
      best_delta_1 = delta_1;
      best_var_idx_2 = var_idx_2;
      best_delta_2 = delta_2;
      best_last_move_step = last_move_step;
    }
  }

  if (best_score > 0)
  {
    ++tightStep;
    ApplyMove(best_var_idx_1, best_delta_1);
    ApplyMove(best_var_idx_2, best_delta_2);
    return true;
  }
  return false;
}

bool Local_ILP::SatTightMove(vector<bool> &score_table,
                             vector<size_t> &score_idx)
{
  long bestScore = 0;
  size_t bestLastMoveStep = std::numeric_limits<size_t>::max();
  size_t bestVarIdx = -1;
  Integer bestDelta = 0;
  vector<size_t> &neighborVarIdxs = localVarUtil.tempVarIdxs;
  vector<Integer> &neighborDeltas = localVarUtil.tempDeltas;
  neighborVarIdxs.clear();
  neighborDeltas.clear();
  auto &neighborConIdxs = localConUtil.tempSatConIdxs;
  neighborConIdxs.clear();
  for (size_t conIdx = 1; conIdx < modelConUtil.conNum; ++conIdx)
  {
    if (localConUtil.conSet[conIdx].gap <= 0 && !modelConUtil.conSet[conIdx].inferSAT)
      neighborConIdxs.push_back(conIdx);
  }
  size_t neighborSize = neighborConIdxs.size();
  if (sampleSat < neighborSize)
  {
    neighborSize = sampleSat;
    for (size_t sampleIdx = 0; sampleIdx < sampleSat; ++sampleIdx)
    {
      size_t randomIdx = mt() % (neighborConIdxs.size() - sampleIdx);
      size_t temp = neighborConIdxs[sampleIdx];
      neighborConIdxs[sampleIdx] = neighborConIdxs[randomIdx + sampleIdx];
      neighborConIdxs[randomIdx + sampleIdx] = temp;
    }
  }
  for (size_t neighborIdx = 0; neighborIdx < neighborSize; ++neighborIdx)
  {
    auto &localCon = localConUtil.conSet[neighborConIdxs[neighborIdx]];
    auto &modelCon = modelConUtil.conSet[neighborConIdxs[neighborIdx]];
    for (size_t termIdx = 0; termIdx < modelCon.termNum; ++termIdx)
    {
      size_t varIdx = modelCon.varIdxSet[termIdx];
      auto &localVar = localVarUtil.GetVar(varIdx);
      auto &modelVar = modelVarUtil.GetVar(varIdx);
      Integer delta;
      if (!TightDelta(localCon, modelCon, termIdx, delta))
        if (modelCon.coeffSet[termIdx] > 0)
          delta = modelVar.upperBound - localVar.nowValue;
        else
          delta = modelVar.lowerBound - localVar.nowValue;
      if (delta < 0 && curStep < localVar.allowDecStep ||
          delta > 0 && curStep < localVar.allowIncStep)
        continue;
      if (delta == 0)
        continue;
      neighborVarIdxs.push_back(varIdx);
      neighborDeltas.push_back(delta);
    }
  }
  size_t scoreSize = neighborVarIdxs.size();
  if (neighborVarIdxs.size() > bmsSat)
  {
    scoreSize = bmsSat;
    for (size_t bmsIdx = 0; bmsIdx < bmsSat; ++bmsIdx)
    {
      size_t randomIdx = (mt() % (neighborVarIdxs.size() - bmsIdx)) + bmsIdx;
      size_t varIdx = neighborVarIdxs[randomIdx];
      Integer delta = neighborDeltas[randomIdx];
      neighborVarIdxs[randomIdx] = neighborVarIdxs[bmsIdx];
      neighborDeltas[randomIdx] = neighborDeltas[bmsIdx];
      neighborVarIdxs[bmsIdx] = varIdx;
      neighborDeltas[bmsIdx] = delta;
    }
  }
  for (size_t idx = 0; idx < scoreSize; ++idx)
  {
    size_t varIdx = neighborVarIdxs[idx];
    Integer delta = neighborDeltas[idx];
    auto &localVar = localVarUtil.GetVar(varIdx);
    auto &modelVar = modelVarUtil.GetVar(varIdx);
    if (isBin)
    {
      if (score_table[varIdx])
        continue;
      else
      {
        score_table[varIdx] = true;
        score_idx.push_back(varIdx);
      }
    }
    long score = TightScore(modelVar, delta);
    size_t lastMoveStep =
        delta < 0 ? localVar.lastDecStep : localVar.lastIncStep;
    if (bestScore < score ||
        bestScore == score && lastMoveStep < bestLastMoveStep)
    {
      bestScore = score;
      bestVarIdx = varIdx;
      bestDelta = delta;
      bestLastMoveStep = lastMoveStep;
    }
  }
  if (isBin)
    for (auto idx : score_idx)
      score_table[idx] = false;

  if (bestScore > 0)
  {
    ++tightStep;
    ApplyMove(bestVarIdx, bestDelta);
    return true;
  }
  return false;
}

bool Local_ILP::UnsatTightMove()
{
  vector<bool> &scoreTable = localVarUtil.scoreTable;
  vector<size_t> scoreIdxs;
  long bestScore = 0;
  size_t bestLastMoveStep = std::numeric_limits<size_t>::max();
  size_t bestVarIdx = -1;
  Integer bestDelta = 0;
  vector<size_t> &neighborVarIdxs = localVarUtil.tempVarIdxs;
  vector<Integer> &neighborDeltas = localVarUtil.tempDeltas;
  neighborVarIdxs.clear();
  neighborDeltas.clear();
  size_t neighborSize = localConUtil.unsatConIdxs.size();
  vector<size_t> *neighborConIdxs = &localConUtil.unsatConIdxs;
  if (sampleUnsat < neighborSize)
  {
    neighborSize = sampleUnsat;
    neighborConIdxs = &localConUtil.tempUnsatConIdxs;
    neighborConIdxs->clear();
    neighborConIdxs->assign(localConUtil.unsatConIdxs.begin(), localConUtil.unsatConIdxs.end());
    for (size_t sampleIdx = 0; sampleIdx < sampleUnsat; ++sampleIdx)
    {
      size_t randomIdx = mt() % (neighborConIdxs->size() - sampleIdx);
      size_t temp = neighborConIdxs->at(sampleIdx);
      neighborConIdxs->at(sampleIdx) = neighborConIdxs->at(randomIdx + sampleIdx);
      neighborConIdxs->at(randomIdx + sampleIdx) = temp;
    }
  }
  for (size_t neighborIdx = 0; neighborIdx < neighborSize; ++neighborIdx)
  {
    auto &localCon = localConUtil.conSet[neighborConIdxs->at(neighborIdx)];
    auto &modelCon = modelConUtil.conSet[neighborConIdxs->at(neighborIdx)];
    for (size_t termIdx = 0; termIdx < modelCon.termNum; ++termIdx)
    {
      size_t varIdx = modelCon.varIdxSet[termIdx];
      auto &localVar = localVarUtil.GetVar(varIdx);
      auto &modelVar = modelVarUtil.GetVar(varIdx);
      Integer delta;
      if (!TightDelta(localCon, modelCon, termIdx, delta))
        if (modelCon.coeffSet[termIdx] > 0)
          delta = modelVar.lowerBound - localVar.nowValue;
        else
          delta = modelVar.upperBound - localVar.nowValue;
      if (delta < 0 && curStep < localVar.allowDecStep ||
          delta > 0 && curStep < localVar.allowIncStep)
        continue;
      if (delta == 0)
        continue;
      neighborVarIdxs.push_back(varIdx);
      neighborDeltas.push_back(delta);
    }
  }
  size_t scoreSize = neighborVarIdxs.size();
  if (!isFoundFeasible && scoreSize > bmsUnsat)
  {
    scoreSize = bmsUnsat;
    for (size_t bmsIdx = 0; bmsIdx < bmsUnsat; ++bmsIdx)
    {
      size_t randomIdx = (mt() % (neighborVarIdxs.size() - bmsIdx)) + bmsIdx;
      size_t varIdx = neighborVarIdxs[randomIdx];
      Integer delta = neighborDeltas[randomIdx];
      neighborVarIdxs[randomIdx] = neighborVarIdxs[bmsIdx];
      neighborDeltas[randomIdx] = neighborDeltas[bmsIdx];
      neighborVarIdxs[bmsIdx] = varIdx;
      neighborDeltas[bmsIdx] = delta;
    }
  }
  for (size_t idx = 0; idx < scoreSize; ++idx)
  {
    size_t varIdx = neighborVarIdxs[idx];
    Integer delta = neighborDeltas[idx];
    auto &localVar = localVarUtil.GetVar(varIdx);
    auto &modelVar = modelVarUtil.GetVar(varIdx);
    if (isBin)
    {
      if (scoreTable[varIdx])
        continue;
      else
      {
        scoreTable[varIdx] = true;
        scoreIdxs.push_back(varIdx);
      }
    }
    long score = TightScore(modelVar, delta);
    size_t lastMoveStep =
        delta < 0 ? localVar.lastDecStep : localVar.lastIncStep;
    if (bestScore < score ||
        bestScore == score && lastMoveStep < bestLastMoveStep)
    {
      bestScore = score;
      bestVarIdx = varIdx;
      bestDelta = delta;
      bestLastMoveStep = lastMoveStep;
    }
  }
  if (bestScore > 0)
  {
    ++tightStep;
    ApplyMove(bestVarIdx, bestDelta);
    if (isBin)
      for (auto idx : scoreIdxs)
        scoreTable[idx] = false;
    return true;
  }
  else
  {
    if (isFoundFeasible && satTightMove)
      return SatTightMove(scoreTable, scoreIdxs);
    else if (!isFoundFeasible && pairTightMove)
      return towMove(neighborVarIdxs, neighborDeltas);
  }
  if (isBin)
    for (auto idx : scoreIdxs)
      scoreTable[idx] = false;
  return false;
}

void Local_ILP::LiftMove()
{
  auto &localObj = localConUtil.conSet[0];
  auto &modelObj = modelConUtil.conSet[0];
  vector<Integer> &lowerDelta = localVarUtil.lowerDeltaInLiftMove;
  vector<Integer> &upperDelta = localVarUtil.upperDeltaInLifiMove;
  if (!isKeepFeas)
  {
    for (size_t termIdx = 0; termIdx < modelObj.termNum; ++termIdx)
    {
      size_t varIdx = modelObj.varIdxSet[termIdx];
      auto &localVar = localVarUtil.GetVar(varIdx);
      auto &modelVar = modelVarUtil.GetVar(varIdx);
      lowerDelta[termIdx] = modelVar.lowerBound - localVar.nowValue;
      upperDelta[termIdx] = modelVar.upperBound - localVar.nowValue;
      for (size_t j = 0; j < modelVar.termNum; ++j)
      {
        size_t conIdx = modelVar.conIdxs[j];
        auto &localCon = localConUtil.conSet[conIdx];
        auto &modelCon = modelConUtil.conSet[conIdx];
        size_t posInCon = modelVar.posInCon[j];
        Integer coeff = modelCon.coeffSet[posInCon];
        if (conIdx == 0)
          continue;
        Integer delta;
        if (!TightDelta(localCon, modelCon, posInCon, delta))
          continue;
        else
        {
          if (coeff > 0)
          {
            if (delta < upperDelta[termIdx])
              upperDelta[termIdx] = delta;
          }
          else if (coeff < 0)
          {
            if (delta > lowerDelta[termIdx])
              lowerDelta[termIdx] = delta;
          }
        }
        if (lowerDelta[termIdx] >= upperDelta[termIdx])
          break;
      }
    }
  }
  Integer bestObjDelta = 0;
  size_t bestVarIdx = -1;
  Integer bestVarDelta = 0;
  Integer varDelta;
  Integer objDelta;
  Integer objDelta_l;
  Integer objDelta_u;
  for (size_t termIdx = 0; termIdx < modelObj.termNum; ++termIdx)
  {
    size_t varIdx = modelObj.varIdxSet[termIdx];
    Integer coeff = modelObj.coeffSet[termIdx];
    Integer l_d = lowerDelta[termIdx];
    Integer u_d = upperDelta[termIdx];
    if (l_d == u_d)
      continue;
    auto &localVar = localVarUtil.GetVar(varIdx);
    auto &modelVar = modelVarUtil.GetVar(varIdx);
    if (!modelVar.inbound(localVar.nowValue + l_d))
      lowerDelta[termIdx] = l_d = 0;
    if (!modelVar.inbound(localVar.nowValue + u_d))
      upperDelta[termIdx] = u_d = 0;
    // assert(modelVar.InBound(localVar.nowValue + l_d) &&
    //        modelVar.InBound(localVar.nowValue + u_d));
    objDelta_l = coeff * l_d;
    objDelta_u = coeff * u_d;
    if (objDelta_l < objDelta_u)
    {
      objDelta = objDelta_l;
      varDelta = l_d;
    }
    else
    {
      objDelta = objDelta_u;
      varDelta = u_d;
    }
    if (objDelta < bestObjDelta)
    {
      bestObjDelta = objDelta;
      bestVarIdx = varIdx;
      bestVarDelta = varDelta;
    }
  }

  if (bestVarIdx != -1 && bestVarDelta != 0)
  {
    ++liftStep;
    ApplyMove(bestVarIdx, bestVarDelta);
    isKeepFeas = true;
    unordered_set<size_t> &affectedVar = localVarUtil.affectedVar;
    affectedVar.clear();
    auto &bestLocalVar = localVarUtil.GetVar(bestVarIdx);
    auto &bestModelVar = modelVarUtil.GetVar(bestVarIdx);
    for (auto conIdx : bestModelVar.conIdxs)
    {
      if (conIdx == 0)
        continue;
      auto &localCon = localConUtil.GetCon(conIdx);
      auto &modelCon = modelConUtil.GetCon(conIdx);
      for (auto varIdx : modelCon.varIdxSet)
        affectedVar.insert(varIdx);
    }
    for (auto varIdx : affectedVar)
    {
      size_t idxInObj = modelVarUtil.varIdx2ObjIdx[varIdx];
      if (idxInObj == -1)
        continue;
      auto &localVar = localVarUtil.GetVar(varIdx);
      auto &modelVar = modelVarUtil.GetVar(varIdx);
      lowerDelta[idxInObj] = modelVar.lowerBound - localVar.nowValue;
      upperDelta[idxInObj] = modelVar.upperBound - localVar.nowValue;
      for (size_t termIdx = 0; termIdx < modelVar.termNum; ++termIdx)
      {
        size_t conIdx = modelVar.conIdxs[termIdx];
        auto &localCon = localConUtil.conSet[conIdx];
        auto &modelCon = modelConUtil.conSet[conIdx];
        size_t posInCon = modelVar.posInCon[termIdx];
        Integer coeff = modelCon.coeffSet[posInCon];
        if (conIdx == 0)
          continue;
        Integer delta;
        if (!TightDelta(localCon, modelCon, posInCon, delta))
          continue;
        else
        {
          if (coeff > 0)
          {
            if (delta < upperDelta[idxInObj])
              upperDelta[idxInObj] = delta;
          }
          else if (coeff < 0)
          {
            if (delta > lowerDelta[idxInObj])
              lowerDelta[idxInObj] = delta;
          }
        }
        if (lowerDelta[idxInObj] >= upperDelta[idxInObj])
          break;
      }
    }
  }
  else
  {
    isKeepFeas = false;
    if (true)
    {
      Config c(sampleUnsat, bmsUnsat, sampleSat, bmsSat, bmsRandom, bmsPair, bmsRan,
               weightUpperBound, objWeightUpperBound, greadyScore, satTightMove, pairTightMove, rvd);
      cbkUpdatePool(worker, localVarUtil.varSet, localObj.rhs + localObj.gap, c);
      // printf("c obj: %s\n", itos(obj.key + obj.gap).c_str());
    }
    for (size_t termIdx = 0; termIdx < modelObj.termNum; ++termIdx)
    {
      size_t randomIdx = mt() % (modelObj.termNum);
      size_t varIdx = modelObj.varIdxSet[randomIdx];
      Integer coeff = modelObj.coeffSet[randomIdx];
      auto &localVar = localVarUtil.GetVar(varIdx);
      auto &modelVar = modelVarUtil.GetVar(varIdx);
      Integer varDelta;
      if (coeff > 0)
        varDelta = -1;
      else
        varDelta = 1;
      if (!modelVar.inbound(varDelta + localVar.nowValue))
        continue;
      else
      {
        ++liftStep;
        ApplyMove(varIdx, varDelta);
        break;
      }
    }
  }
}

double Local_ILP::TightScore(
    const ModelVar &modelVar,
    Integer delta)
{
  long score = 0;
  size_t conIdx;
  size_t posInCon;
  Integer newGap;
  bool isPreSat;
  bool isNowSat;
  for (size_t termIdx = 0; termIdx < modelVar.termNum; ++termIdx)
  {
    conIdx = modelVar.conIdxs[termIdx];
    posInCon = modelVar.posInCon[termIdx];
    auto &localCon = localConUtil.conSet[conIdx];
    auto &modelCon = modelConUtil.conSet[conIdx];
    if (conIdx == 0)
    {
      if (isFoundFeasible)
      {
        newGap =
            localCon.gap + modelCon.coeffSet[posInCon] * delta;
        if (1 > newGap)
          score += localCon.weight;
        else if (greadyScore || 1 < newGap)
          score -= localCon.weight;
      }
    }
    else
    {
      newGap =
          localCon.gap + modelCon.coeffSet[posInCon] * delta;
      isPreSat = localCon.gap <= 0;
      isNowSat = newGap <= 0;
      if (!isPreSat && isNowSat)
        score += localCon.weight;
      else if (isPreSat && !isNowSat)
        score -= localCon.weight;
      else if (!isPreSat && !isNowSat)
        if (localCon.gap > newGap)
          score += localCon.weight * rvd;
        else if (greadyScore || localCon.gap < newGap)
          score -= localCon.weight * rvd;
    }
  }
  return score;
}

// return delta_x
// a * delta_x + gap <= 0
bool Local_ILP::TightDelta(
    LocalCon &con,
    const ModelCon &modelCon,
    size_t termIdx,
    Integer &res)
{
  auto varIdx = modelCon.varIdxSet[termIdx];
  auto &localVar = localVarUtil.GetVar(varIdx);
  auto &modelVar = modelVarUtil.GetVar(varIdx);
  double delta =
      -((double)con.gap / modelCon.coeffSet[termIdx]);
  if (modelCon.coeffSet[termIdx] > 0)
    res = floor(delta);
  else
    res = ceil(delta);
  if (modelVar.inbound(localVar.nowValue + res))
    return true;
  else
    return false;
}

void Local_ILP::ApplyMove(size_t var_idx, Integer delta)
{
  // printf("c %-10d varIdx: %-10d delta: %-10d\n", curStep, var_idx, delta);
  auto &var = localVarUtil.GetVar(var_idx);
  auto &common_var = modelVarUtil.GetVar(var_idx);
  var.nowValue += delta;

  for (size_t i = 0; i < common_var.conIdxs.size(); ++i)
  {
    size_t con_idx = common_var.conIdxs[i];
    size_t pos = common_var.posInCon[i];
    auto &con = localConUtil.conSet[con_idx];
    auto &common_con = modelConUtil.conSet[con_idx];
    Integer new_gap = con.gap + common_con.coeffSet[pos] * delta;

    if (con_idx == 0)
    { // obj
      con.gap = new_gap;
      continue;
    }

    bool is_pre_sat = con.gap <= 0;
    bool is_now_sat = new_gap <= 0;

    // update unsat_con_idxs
    if (is_pre_sat && !is_now_sat)
    {
      localConUtil.insertUnsat(con_idx);
    }
    else if (!is_pre_sat && is_now_sat)
    {
      localConUtil.RemoveUnsat(con_idx);
    }

    // update gaps
    con.gap = new_gap;
  }

  // update tabu and last_move_step
  if (delta > 0)
  {
    var.lastIncStep = curStep;
    var.allowDecStep = curStep + tabuBase + mt() % tabuVariation;
  }
  else
  {
    // assert(delta < 0);
    var.lastDecStep = curStep;
    var.allowIncStep = curStep + tabuBase + mt() % tabuVariation;
  }
}

void Local_ILP::RandomTightMove_2()
{
  auto max_weight_idx = localConUtil.unsatConIdxs[0];
  auto max_weight = localConUtil.conSet[max_weight_idx].weight;
  for (size_t i = 1; i < localConUtil.unsatConIdxs.size(); ++i)
  {
    auto weight_idx = localConUtil.unsatConIdxs[i];
    auto weight = localConUtil.conSet[weight_idx].weight;
    if (weight > max_weight)
    {
      max_weight = weight;
      max_weight_idx = weight_idx;
    }
    if (max_weight >= weightUpperBound)
      break;
  }
  auto &con = localConUtil.conSet[max_weight_idx];
  auto &common_con = modelConUtil.conSet[max_weight_idx];
  size_t ran_idx = mt() % common_con.termNum;
  auto var_idx = common_con.varIdxSet[ran_idx];
  auto coffe = common_con.coeffSet[ran_idx];
  Integer delta = 0;
  if (coffe > 0)
    delta = modelVarUtil.varSet[var_idx].lowerBound - localVarUtil.varSet[var_idx].nowValue;
  else
    delta = modelVarUtil.varSet[var_idx].upperBound - localVarUtil.varSet[var_idx].nowValue;

  if (delta != 0)
  {
    ++randomStep;
    ApplyMove(var_idx, delta);
    return;
  }
}

// random select one unsat constraint
void Local_ILP::RandomTightMove()
{
  long bestScore = -100000000000;
  size_t bestLastMoveStep = std::numeric_limits<size_t>::max();
  size_t bestVarIdx = -1;
  Integer bestDelta = 0;
  size_t conIdx = localConUtil.unsatConIdxs[mt() % localConUtil.unsatConIdxs.size()];
  auto &localCon = localConUtil.conSet[conIdx];
  auto &modelCon = modelConUtil.conSet[conIdx];
  vector<size_t> &neighborVarIdxs = localVarUtil.tempVarIdxs;
  vector<Integer> &neighborDeltas = localVarUtil.tempDeltas;
  neighborVarIdxs.clear();
  neighborDeltas.clear();
  for (size_t termIdx = 0; termIdx < modelCon.termNum; ++termIdx)
  {
    size_t varIdx = modelCon.varIdxSet[termIdx];
    auto &localVar = localVarUtil.GetVar(varIdx);
    auto &modelVar = modelVarUtil.GetVar(varIdx);
    Integer delta;
    if (!TightDelta(localCon, modelCon, termIdx, delta))
    {
      if (modelCon.coeffSet[termIdx] > 0)
        delta = modelVar.lowerBound - localVar.nowValue;
      else
        delta = modelVar.upperBound - localVar.nowValue;
    }
    if (
        (delta < 0 && curStep == localVar.lastIncStep + 1 ||
         delta > 0 && curStep == localVar.lastDecStep + 1))
      continue;
    if (delta == 0)
      continue;
    neighborVarIdxs.push_back(varIdx);
    neighborDeltas.push_back(delta);
  }
  size_t scoreSize = neighborVarIdxs.size();
  if (neighborVarIdxs.size() > bmsRandom)
  {
    scoreSize = bmsRandom;
    for (size_t bmsIdx = 0; bmsIdx < bmsRandom; ++bmsIdx)
    {
      size_t randomIdx = (mt() % (neighborVarIdxs.size() - bmsIdx)) + bmsIdx;
      size_t varIdx = neighborVarIdxs[randomIdx];
      Integer delta = neighborDeltas[randomIdx];
      neighborVarIdxs[randomIdx] = neighborVarIdxs[bmsIdx];
      neighborDeltas[randomIdx] = neighborDeltas[bmsIdx];
      neighborVarIdxs[bmsIdx] = varIdx;
      neighborDeltas[bmsIdx] = delta;
    }
  }
  for (size_t idx = 0; idx < scoreSize; ++idx)
  {
    size_t varIdx = neighborVarIdxs[idx];
    Integer delta = neighborDeltas[idx];
    auto &localVar = localVarUtil.GetVar(varIdx);
    auto &modelVar = modelVarUtil.GetVar(varIdx);
    long score = TightScore(modelVar, delta);
    size_t lastMoveStep =
        delta < 0 ? localVar.lastDecStep : localVar.lastIncStep;
    if (bestScore < score ||
        bestScore == score && lastMoveStep < bestLastMoveStep)
    {
      bestScore = score;
      bestVarIdx = varIdx;
      bestDelta = delta;
      bestLastMoveStep = lastMoveStep;
    }
  }

  if (bestVarIdx != -1 && bestDelta != 0)
  {
    ++randomStep;
    ApplyMove(bestVarIdx, bestDelta);
    return;
  }
  // else
  //   RandomTightMove_2();
  // TODO no availiable one-step critical_move
}

void Local_ILP::Sharing()
{
  lastImproveStep = curStep;
  ++sharingTimes;
  bool res = cbkGetSol_Config(worker);
  if (res == false)
    return;
  for (auto unsat_idx : localConUtil.unsatConIdxs)
    localConUtil.RemoveUnsat(unsat_idx);
  for (int var_idx = 0; var_idx < modelVarUtil.varNum; var_idx++)
  {
    auto &var = localVarUtil.GetVar(var_idx);
    var.lastDecStep = curStep;
    var.allowIncStep = 0;
    var.lastIncStep = curStep;
    var.allowDecStep = 0;
  }
  for (size_t con_idx = 1; con_idx < localConUtil.conSet.size(); ++con_idx)
  {
    auto &con = localConUtil.conSet[con_idx];
    auto &common_con = modelConUtil.conSet[con_idx];
    // calculate gap
    con.gap = 0;
    for (size_t i = 0; i < common_con.termNum; ++i)
      con.gap += common_con.coeffSet[i] * localVarUtil.GetVar(common_con.varIdxSet[i]).nowValue;

    con.gap -= con.rhs;

    // setup unsat_con_idxs and total_weight
    if (con.gap > 0)
      localConUtil.insertUnsat(con_idx);
    if (tid > 16)
      con.weight = 1;
  }

  // Obj
  auto &obj = localConUtil.conSet[0];
  auto &modelObj = modelConUtil.conSet[0];
  obj.gap = 0;
  if (tid > 16)
    obj.weight = 1;
  for (size_t i = 0; i < modelObj.termNum; ++i)
    obj.gap += modelObj.coeffSet[i] * localVarUtil.GetVar(modelObj.varIdxSet[i]).nowValue;
  obj.gap -= obj.rhs;
}

void Local_ILP::Restart()
{
  lastImproveStep = curStep;
  ++restartTimes;
  for (auto unsat_idx : localConUtil.unsatConIdxs)
  {
    localConUtil.RemoveUnsat(unsat_idx);
  }
  // for (auto &var : var_util.vars)
  for (int var_idx = 0; var_idx < modelVarUtil.varNum; var_idx++)
  {
    auto &var = localVarUtil.GetVar(var_idx);
    auto &common_var = modelVarUtil.GetVar(var_idx);
    if (common_var.upperBound + 1 == common_var.lowerBound)
      var.nowValue = (Integer)(common_var.lowerBound / 2.0 + common_var.upperBound / 2.0); // avoid overflow for free var
    else
      var.nowValue = common_var.lowerBound + mt() % (common_var.upperBound + 1 - common_var.lowerBound);
    if (!common_var.inbound(var.nowValue))
    {
      var.nowValue = (Integer)(common_var.lowerBound / 2.0 + common_var.upperBound / 2.0); // avoid overflow for free var
    }
    if (isFoundFeasible && mt() % 100 > 50)
      var.nowValue = var.bestValue;
    var.lastDecStep = curStep;
    var.allowIncStep = 0;
    var.lastIncStep = curStep;
    var.allowDecStep = 0;
  }
  for (size_t con_idx = 1; con_idx < localConUtil.conSet.size(); ++con_idx)
  {
    auto &con = localConUtil.conSet[con_idx];
    auto &common_con = modelConUtil.conSet[con_idx];
    // calculate gap
    con.gap = 0;
    for (size_t i = 0; i < common_con.termNum; ++i)
    {
      con.gap += common_con.coeffSet[i] *
                 localVarUtil.GetVar(common_con.varIdxSet[i]).nowValue;
    }

    con.gap -= con.rhs;

    // setup unsat_con_idxs and total_weight
    if (con.gap > 0)
    {
      localConUtil.insertUnsat(con_idx);
    }
    con.weight = 1;
  }

  // Obj
  auto &obj = localConUtil.conSet[0];
  auto &modelObj = modelConUtil.conSet[0];
  obj.gap = 0;
  obj.weight = 1;
  for (size_t i = 0; i < modelObj.termNum; ++i)
  {
    obj.gap += modelObj.coeffSet[i] *
               localVarUtil.GetVar(modelObj.varIdxSet[i]).nowValue;
  }
  obj.gap -= obj.rhs;
}

bool Local_ILP::VerifySolution()
{
  for (size_t var_idx = 0; var_idx < modelVarUtil.varNum; var_idx++)
  {
    auto &var = localVarUtil.GetVar(var_idx);
    auto &common_var = modelVarUtil.GetVar(var_idx);
    if (!common_var.inbound(var.bestValue))
      return false;
  }

  for (size_t con_idx = 1; con_idx < localConUtil.conSet.size(); ++con_idx)
  {
    auto &con = localConUtil.conSet[con_idx];
    auto &common_con = modelConUtil.conSet[con_idx];
    // calculate gap
    Integer con_gap = 0;
    for (size_t i = 0; i < common_con.termNum; ++i)
      con_gap += common_con.coeffSet[i] * localVarUtil.GetVar(common_con.varIdxSet[i]).bestValue;
    // setup unsat_con_idxs and total_weight
    if (con_gap > common_con.rhs)
      cout << "con_gap:\t" << (double)con_gap << endl
           << "common_con.ori_key:\t" << (double)common_con.rhs << endl;
    if (con_gap > con.rhs)
      return false;
  }

  // Obj
  auto &obj = localConUtil.conSet[0];
  auto &modelObj = modelConUtil.conSet[0];
  Integer obj_gap = 0;
  for (size_t i = 0; i < modelObj.termNum; ++i)
    obj_gap += modelObj.coeffSet[i] * localVarUtil.GetVar(modelObj.varIdxSet[i]).bestValue;
  obj_gap -= obj.rhs;
  return obj_gap == 1;
}

void Local_ILP::PrintResult()
{
  if (isFoundFeasible && VerifySolution())
  {
    if (OPT(PrintSol))
      PrintSol();
  }
  else
    cout << "solution verify failed." << endl;
}

void Local_ILP::PrintSol()
{
  for (int i = 0; i < modelVarUtil.varNum; i++)
  {
    const auto &var = localVarUtil.GetVar(i);
    const auto &common_var = modelVarUtil.GetVar(i);
    if (var.bestValue)
      printf("%24.24s        %lld\n", common_var.name.c_str(), (long long)var.bestValue);
  }
}
// TODO obj case
void Local_ILP::UpdateWeight()
{
  for (size_t idx : localConUtil.unsatConIdxs)
  {
    auto &con = localConUtil.conSet[idx];
    if (con.weight < weightUpperBound)
      ++con.weight;
  }
  //  obj case
  auto &obj = localConUtil.conSet[0];
  if (isFoundFeasible && obj.gap > 0 && obj.weight <= objWeightUpperBound)
    ++obj.weight;
}

void Local_ILP::SmoothWeight()
{
  for (auto &con : localConUtil.conSet)
    if (con.gap <= 0 && con.weight > 0)
      --con.weight;
}

void Local_ILP::Allocate()
{
  localVarUtil.Allocate(modelVarUtil.varNum, modelConUtil.conSet[0].termNum);
  localConUtil.Allocate(modelConUtil.conNum);
  for (size_t con_idx = 1; con_idx < modelConUtil.conNum; con_idx++)
    localConUtil.conSet[con_idx].rhs = modelConUtil.conSet[con_idx].rhs;
}

Local_ILP::Local_ILP(const ModelConUtil &_modelConUtil, const ModelVarUtil &_modelVarUtil)
    : modelConUtil(_modelConUtil), modelVarUtil(_modelVarUtil)
{
  smoothProbability = 3;
  tabuBase = 3;
  tabuVariation = 10;
  isFoundFeasible = false;
  liftStep = 0;
  tightStep = 0;
  randomStep = 0;
  weightUpperBound = 1000;
  objWeightUpperBound = 100;
  lastImproveStep = 0;
  restartTimes = 0;
  sharingTimes = 0;
  isBin = true;
  isKeepFeas = false;
  sampleUnsat = 3;
  bmsUnsat = 2000;
  sampleSat = 30;
  bmsSat = 350;
  bmsRandom = 60;
  bmsPair = 77;
  bmsRan = 150;
  restartStep = 1500000;
  greadyScore = true;
  rvd = 0.5;
  satTightMove = true;
  pairTightMove = false;
}
Local_ILP::~Local_ILP()
{
}
