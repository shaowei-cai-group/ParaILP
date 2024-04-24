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
#include "Master.hpp"

std::atomic<int> terminated;

void *WorkerSolve(void *arg)
{
  Worker *worker = (Worker *)arg;
  Master *prs = worker->master;
  worker->Diversity();
  worker->InitSol();
  int res = worker->Solve();
  if (res == 1)
  {
    terminated = 1;
    prs->TerminateWorker();
  }
  return NULL;
}

void Master::InitWorkerSet()
{
  terminated = 0;
  for (int tid = 0; tid < OPT(nThreads); tid++)
  {
    Worker *solver = new Worker(tid, this);
    workerSet.push_back(solver);
  }
}

void Master::TerminateWorker()
{
  for (int tid = 0; tid < OPT(nThreads); tid++)
    workerSet[tid]->Terminate();
}

void Master::Solve()
{
  printf("c -----------------solve start----------------------\n");
  printf("c %-40s\t%-20s\n", "Obj Value", "Found Time");
  pthread_t *ptr = new pthread_t[OPT(nThreads)];
  for (int tid = 0; tid < OPT(nThreads); tid++)
    pthread_create(&ptr[tid], NULL, WorkerSolve, workerSet[tid]);

  while (!terminated)
  {
    usleep(1000000);
    auto clk_now = chrono::high_resolution_clock::now();
    int solve_time = chrono::duration_cast<chrono::seconds>(clk_now - clk_st_global).count();
    if (solve_time >= OPT(cutoff))
    {
      terminated = 1;
      TerminateWorker();
    }
  }
  printf("c -----------------ending solve----------------------\n");

  for (int i = 0; i < OPT(nThreads); i++)
    pthread_join(ptr[i], NULL);
  printf("c -----------------ending join----------------------\n");

  int winner_tid = -1;
  for (int tid = 0; tid < OPT(nThreads); tid++)
  {
    if (workerSet[tid]->isFeasible == false)
      continue;
    if (winner_tid == -1)
      winner_tid = tid;
    else if (workerSet[tid]->best_obj < workerSet[winner_tid]->best_obj ||
             workerSet[tid]->best_obj == workerSet[winner_tid]->best_obj &&
                 workerSet[tid]->bestTime < workerSet[winner_tid]->bestTime)
      winner_tid = tid;
  }
  if (winner_tid != -1)
  {
    printf("o winner thread is %d\no found time is %.3f seconds.\no best found objective value is %s\n",
           winner_tid, workerSet[winner_tid]->bestTime, itos(workerSet[winner_tid]->best_obj).c_str());
    workerSet[winner_tid]->PrintResult();
  }
  else
    printf("o no feasible solution found\n");
  delete[] ptr;
}

void Master::Run()
{
  ReadMPS();
  InitWorkerSet();
  Solve();
}

Master::Master()
    : bestObj(MaxValue)
{
  clk_st_global = chrono::high_resolution_clock::now();
  fileName = (char *)OPT(instance).c_str();
  optValue = __global_paras.identify_opt(fileName);
  if (optValue > -1e37)
    printf("c optValue: %s\n", itos(optValue).c_str());
}
Master::~Master()
{
  for (int tid = 0; tid < workerSet.size(); tid++)
    delete (workerSet[tid]);
}

void Master::TightenBound()
{
  for (size_t conIdx = 1; conIdx < modelConUtil.conNum; ++conIdx)
  {
    ModelCon &modelCon = modelConUtil.conSet[conIdx];
    if (modelCon.varIdxSet.size() == 1)
    {
      TightenBoundVar(modelCon);
    }
  }
}

void Master::TightenBoundVar(ModelCon &modelCon)
{
  Integer coeff = modelCon.coeffSet[0];
  ModelVar &modelVar = modelVarUtil.GetVar(modelCon.varIdxSet[0]);
  double bound = (double)modelCon.rhs / coeff;
  if (coeff > 0 && bound <= modelVar.upperBound) // x <= bound
    modelVar.upperBound = floor(bound);
  else if (coeff < 0 && modelVar.lowerBound <= bound) // x >= bound
    modelVar.lowerBound = ceil(bound);
}

bool Master::TightBoundGlobally()
{
  vector<size_t> fixedIdxs;
  for (auto &modelVar : modelVarUtil.varSet)
    if (modelVar.lowerBound == modelVar.upperBound)
    {
      modelVar.SetType(VarType::Fixed);
      fixedIdxs.push_back(modelVar.idx);
    }
  while (fixedIdxs.size() > 0)
  {
    size_t removeVarIdx = fixedIdxs.back();
    fixedIdxs.pop_back();
    deleteVarNum++;
    auto &removeVar = modelVarUtil.GetVar(removeVarIdx);
    Integer removeVarValue = removeVar.lowerBound;
    for (int termIdx = 0; termIdx < removeVar.conIdxs.size(); termIdx++)
    {
      size_t conIdx = removeVar.conIdxs[termIdx];
      size_t posInCon = removeVar.posInCon[termIdx];
      auto &modelCon = modelConUtil.GetCon(conIdx);
      auto coeff = modelCon.coeffSet[posInCon];
      size_t movedVarIdx = modelCon.varIdxSet.back();
      auto movedCoeff = modelCon.coeffSet.back();
      size_t movedPosInVar = modelCon.posInVar.back();
      modelCon.varIdxSet[posInCon] = movedVarIdx;
      modelCon.coeffSet[posInCon] = movedCoeff;
      modelCon.posInVar[posInCon] = movedPosInVar;
      auto &movedVar = modelVarUtil.GetVar(movedVarIdx);
      assert(movedVar.conIdxs[movedPosInVar] == conIdx);
      movedVar.posInCon[movedPosInVar] = posInCon;
      modelCon.varIdxSet.pop_back();
      modelCon.coeffSet.pop_back();
      modelCon.posInVar.pop_back();
      if (conIdx == 0)
        modelVarUtil.objBias += coeff * removeVarValue;
      else
      {
        modelCon.rhs -= coeff * removeVarValue;
        if (modelCon.varIdxSet.size() == 1)
        {
          TightenBoundVar(modelCon);
          auto &relatedVar = modelVarUtil.GetVar(modelCon.varIdxSet[0]);
          if (relatedVar.type != VarType::Fixed &&
              relatedVar.lowerBound == relatedVar.upperBound)
          {
            relatedVar.SetType(VarType::Fixed);
            fixedIdxs.push_back(relatedVar.idx);
            inferVarNum++;
          }
        }
        else if (modelCon.varIdxSet.size() == 0)
        {
          assert(modelCon.varIdxSet.size() == 0 && modelCon.posInVar.size() == 0);
          if (modelCon.rhs >= 0)
          {
            modelCon.inferSAT = true;
            deleteConNum++;
          }
          else
          {
            printf("c con.rhs %d\n", modelCon.rhs);
            return false;
          }
        }
      }
    }
  }
  return true;
}

bool Master::SetVarType()
{
  for (size_t varIdx = 0; varIdx < modelVarUtil.varNum; varIdx++)
  {
    auto &modelVar = modelVarUtil.GetVar(varIdx);
    modelVar.termNum = modelVar.conIdxs.size();
    if (modelVar.lowerBound == modelVar.upperBound)
    {
      modelVarUtil.fixedNum++;
      modelVar.SetType(VarType::Fixed);
    }
    else if (
        modelVar.lowerBound == 0 &&
        modelVar.upperBound == 1)
    {
      modelVarUtil.binaryNum++;
      modelVar.SetType(VarType::Binary);
    }
    else
    {
      modelVar.SetType(VarType::Integer);
      modelVarUtil.integerNum++;
      modelVarUtil.isBin = false;
    }
  }
  for (size_t conIdx = 0; conIdx < modelConUtil.conNum; conIdx++)
  {
    auto &modelCon = modelConUtil.GetCon(conIdx);
    modelCon.termNum = modelCon.varIdxSet.size();
    if (modelCon.inferSAT)
      assert(modelCon.termNum == 0);
  }
  return true;
}

void Master::PushCoeffVarIdx(
    const size_t conIdx,
    const Integer coeff,
    const string &varName)
{
  ModelCon &modelCon = modelConUtil.conSet[conIdx];
  size_t varIdx = modelVarUtil.MakeVar(varName);
  ModelVar &modelVar = modelVarUtil.GetVar(varIdx);

  modelVar.conIdxs.push_back(conIdx);
  modelVar.posInCon.push_back(modelCon.varIdxSet.size());

  modelCon.coeffSet.push_back(coeff);
  modelCon.varIdxSet.push_back(varIdx);
  modelCon.posInVar.push_back(modelVar.conIdxs.size() - 1);
}

inline void Master::IssSetup(istringstream &iss, const string &line)
{
  iss.clear();
  iss.str(line);
  iss.seekg(0, ios::beg);
}

void Master::ReadMPS()
{
  istringstream iss;
  string readLine;
  ifstream inFile(fileName);
  string tmpStr;
  char conType;
  string conName;
  size_t conIdx;
  string inverseConName;
  size_t inverseConIdx;
  string varName;
  Integer coefficient;
  double tempVal;
  Integer rhs;
  string varType;
  double inputBound;

  if (!inFile)
  {
    printf("The input filename %s is invalid\n", fileName);
    exit(-1);
  }

  while (getline(inFile, readLine)) // NAME section
  {
    if (readLine[0] == '*' || readLine.length() < 1)
      continue;
    if (readLine[0] == 'R')
      break;
    IssSetup(iss, readLine);
    if (!(iss >> tmpStr >> modelName))
      continue;
    if (modelName.length() < 1)
      continue;
    printf("c model name: %s\n", modelName.c_str());
  }

  modelConUtil.conSet.emplace_back("", 0); // con[0] is objective function
  while (getline(inFile, readLine))        // ROWS section
  {
    if (readLine[0] == '*' || readLine.length() < 1)
      continue;
    if (readLine[0] == 'C')
      break;
    IssSetup(iss, readLine);
    if (!(iss >> conType >> conName))
      continue;
    if (conType == 'L')
    {
      conIdx = modelConUtil.MakeCon(conName);
      modelConUtil.conSet[conIdx].isEqual = false;
      modelConUtil.conSet[conIdx].isLess = true;
    }
    else if (conType == 'E')
    {
      conIdx = modelConUtil.MakeCon(conName);
      modelConUtil.conSet[conIdx].isEqual = true;
      modelConUtil.conSet[conIdx].isLess = false;
      inverseConName = conName + "!";
      inverseConIdx = modelConUtil.MakeCon(inverseConName);
      modelConUtil.conSet[inverseConIdx].isEqual = true;
      modelConUtil.conSet[inverseConIdx].isLess = false;
    }
    else if (conType == 'G')
    {
      conIdx = modelConUtil.MakeCon(conName);
      modelConUtil.conSet[conIdx].isEqual = false;
      modelConUtil.conSet[conIdx].isLess = false;
    }
    else
    {
      assert(conType == 'N'); // type=='N',this con is obj
      modelConUtil.objName = conName;
      modelConUtil.conSet[0].isEqual = false;
      modelConUtil.conSet[0].isLess = false;
    }
  }

  while (getline(inFile, readLine)) // COLUMNS section
  {
    if (readLine[0] == '*' ||
        readLine.length() < 1)
      continue;
    if (readLine[0] == 'R')
      break;
    IssSetup(iss, readLine);
    if (!(iss >> varName >> conName))
      continue;
    if (conName == "\'MARKER\'")
    {
      iss >> tmpStr;
      if (tmpStr == "\'INTORG\'" || tmpStr == "\'INTEND\'")
        continue;
      else
      {
        printf("c error %s\n", readLine.c_str());
        exit(-1);
      }
    }
    iss >> tempVal;
    coefficient = tempVal * ZoomTimes;
    conIdx = modelConUtil.GetConIdx(conName);
    PushCoeffVarIdx(conIdx, coefficient, varName);
    if (modelConUtil.conSet[conIdx].isEqual)
      PushCoeffVarIdx(conIdx + 1, -coefficient, varName);
    if (iss >> conName)
    {
      iss >> tempVal;
      coefficient = tempVal * ZoomTimes;
      conIdx = modelConUtil.GetConIdx(conName);
      PushCoeffVarIdx(conIdx, coefficient, varName);
      if (modelConUtil.conSet[conIdx].isEqual)
        PushCoeffVarIdx(conIdx + 1, -coefficient, varName);
    }
  }

  while (getline(inFile, readLine)) // RHS  section
  {
    if (readLine[0] == '*' ||
        readLine.length() < 1)
      continue;
    if (readLine[0] == 'B' ||
        readLine[0] == 'E')
      break;
    // assert(line[0] != 'R'); // do not handle RANGS and SOS
    IssSetup(iss, readLine);
    if (!(iss >> tmpStr >> conName >> tempVal))
      continue;
    if (conName.length() < 1)
      continue;
    rhs = tempVal * ZoomTimes;
    conIdx = modelConUtil.GetConIdx(conName);
    modelConUtil.conSet[conIdx].rhs = rhs;
    if (modelConUtil.conSet[conIdx].isEqual)
      modelConUtil.conSet[conIdx + 1].rhs = -rhs;

    if (iss >> conName)
    {
      iss >> tempVal;
      rhs = tempVal * ZoomTimes;
      conIdx = modelConUtil.GetConIdx(conName);
      modelConUtil.conSet[conIdx].rhs = rhs;
      if (modelConUtil.conSet[conIdx].isEqual)
        modelConUtil.conSet[conIdx + 1].rhs = -rhs;
    }
  }

  while (getline(inFile, readLine)) // BOUNDS section
  {
    if (readLine[0] == '*' ||
        readLine.length() < 1)
      continue;
    if (readLine[0] == 'E')
      break;
    assert(readLine[0] != 'I'); // do not handle INDICATORS
    IssSetup(iss, readLine);
    if (!(iss >> varType >> tmpStr >> varName))
      continue;
    iss >> inputBound;
    ModelVar &var = modelVarUtil.GetVar(varName);
    if (var.type == VarType::Binary)
    {
      var.SetType(VarType::Integer);
      var.upperBound = InfiniteUpperBound;
    }
    if (varType == "UP")
      var.upperBound = floor(inputBound);
    else if (varType == "LO")
      var.lowerBound = ceil(inputBound);
    else if (varType == "BV")
    {
      var.upperBound = 1;
      var.lowerBound = 0;
    }
    else if (varType == "LI")
      var.lowerBound = ceil(inputBound);
    else if (varType == "UI")
      var.upperBound = floor(inputBound);
    else if (varType == "FX")
    {
      var.lowerBound = inputBound;
      var.upperBound = inputBound;
    }
    else if (varType == "FR")
    {
      var.upperBound = InfiniteUpperBound;
      var.lowerBound = InfiniteLowerBound;
    }
    else if (varType == "MI")
      var.lowerBound = InfiniteLowerBound;
    else if (varType == "PL")
      var.upperBound = InfiniteUpperBound;
  }
  inFile.close();
  modelConUtil.conNum = modelConUtil.conSet.size();
  modelVarUtil.varNum = modelVarUtil.varSet.size();
  for (size_t conIdx = 1; conIdx < modelConUtil.conNum; ++conIdx) // deal con_type =='G'
  {
    ModelCon &con = modelConUtil.conSet[conIdx];
    if (con.isLess == false && con.isEqual == false)
    {
      con.isLess = true;
      for (Integer &inverseCoefficient : con.coeffSet)
        inverseCoefficient = -inverseCoefficient;
      con.rhs = -con.rhs;
    }
  }
  auto now = chrono::high_resolution_clock::now();
  double readTime = chrono::duration_cast<chrono::milliseconds>(now - clk_st_global).count() / 1000.0;
  TightenBound();
  if (!TightBoundGlobally())
  {
    printf("c model is infeasible.\n");
    exit(-1);
  }
  SetVarType();
  modelVarUtil.varIdx2ObjIdx.resize(modelVarUtil.varNum, -1);
  ModelCon &modelObj = modelConUtil.conSet[0];
  for (size_t idx = 0; idx < modelObj.termNum; ++idx)
    modelVarUtil.varIdx2ObjIdx[modelObj.varIdxSet[idx]] = idx;
  sharingPeriod = SharingPeriod;
  lhsAbsSum.resize(modelConUtil.conNum, 0);
  for (size_t conIdx = 0; conIdx < modelConUtil.conNum; conIdx++)
    for (size_t idx = 0; idx < modelConUtil.conSet[conIdx].termNum; idx++)
    {
      Integer coeff = modelConUtil.conSet[conIdx].coeffSet[idx];
      size_t varIdx = modelConUtil.conSet[conIdx].varIdxSet[idx];
      ModelVar &modelVar = modelVarUtil.GetVar(varIdx);
      Integer varRange = modelVar.upperBound - modelVar.lowerBound;
      double lhs_i = ((double)coeff / (double)ZoomTimes) * varRange;
      lhsAbsSum[conIdx] += lhs_i > 0 ? lhs_i : -lhs_i;
    }
  Polarity.resize(modelVarUtil.varNum, 0);
  for (size_t varIdx = 0; varIdx < modelVarUtil.varNum; varIdx++)
  {
    auto &modelVar = modelVarUtil.GetVar(varIdx);
    double IC = 0;
    double IO = 0;
    Integer varRange = modelVar.upperBound - modelVar.lowerBound;
    if (varRange == 0)
      continue;
    for (size_t idx = 0; idx < modelVar.conIdxs.size(); ++idx)
    {
      size_t conIdx = modelVar.conIdxs[idx];
      size_t pos = modelVar.posInCon[idx];
      Integer coeff = modelConUtil.conSet[conIdx].coeffSet[pos];
      double lhs_i = (coeff / (double)ZoomTimes) * varRange;
      if (conIdx == 0)
        IO = lhs_i / (lhsAbsSum[conIdx]);
      else
        IC += lhs_i / (lhsAbsSum[conIdx]);
    }
    if (modelVar.conIdxs.size() > 1)
      IC /= (modelVar.conIdxs.size() - 1);
    Polarity[varIdx] = IC + IO;
  }
}
