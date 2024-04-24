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
#include "LocalVar.h"
#include "LocalCon.h"
#include "share/ModelVar.h"
#include "share/ModelCon.h"
#include "share/Population.h"
#include "../utils/paras.h"
class ModelVarUtil;
class ModelVar;
class ModelConUtil;
class ModelCon;
class Config;

class Local_ILP
{
public:
    /*++++++++++++++++++++++for parrale control++++++++++++*/

    void *worker;
    int terminated;
    Integer optValue;
    long long sharingPeriod;
    int tid;

    void (*cbkUpdateSol)(void *_worker, const Integer obj);
    void (*cbkUpdatePool)(void *_worker, const vector<LocalVar> &var, const Integer obj, const Config &config);
    bool (*cbkGetSol_Config)(void *_worker);

    /*----------------------original----------------------*/
    LocalVarUtil localVarUtil;
    LocalConUtil localConUtil;
    const ModelVarUtil &modelVarUtil;
    const ModelConUtil &modelConUtil;
    long curStep;
    std::mt19937 mt;
    int smoothProbability;
    int tabuBase;
    int tabuVariation;
    bool isFoundFeasible;
    size_t liftStep;
    size_t tightStep;
    size_t randomStep;
    size_t weightUpperBound;
    size_t objWeightUpperBound;
    size_t lastImproveStep;
    size_t restartTimes;
    size_t sharingTimes;
    bool isBin;
    bool isKeepFeas;
    size_t sampleUnsat;
    size_t bmsUnsat;
    size_t sampleSat;
    size_t bmsSat;
    size_t bmsRandom;
    size_t bmsPair;
    size_t bmsRan;
    size_t restartStep;
    bool greadyScore;
    double rvd;
    bool satTightMove;
    bool pairTightMove;

public:
    Local_ILP(
        const ModelConUtil &_modelConUtil,
        const ModelVarUtil &_modelVarUtil);
    ~Local_ILP();
    int LocalSearch();
    void PrintResult();
    void PrintSol();
    void Allocate();

private:
    bool VerifySolution();
    Integer GetObjValue();
    void InitState();
    void UpdateBestSolution();
    void Restart();
    void Sharing();
    bool UnsatTightMove();
    void RandomTightMove();
    void RandomTightMove_2();
    void LiftMove();
    bool SatTightMove(
        vector<bool> &score_table,
        vector<size_t> &score_idx);
    void makePair(
        vector<size_t> &var_idxs, vector<Integer> &deltas,
        vector<size_t> &var_idxs_1, vector<Integer> &deltas_1,
        vector<size_t> &var_idxs_2, vector<Integer> &deltas_2);
    double towScore(
        const ModelVar &var_1,
        Integer delta_1,
        const ModelVar &var_2,
        Integer delta_2);
    bool towMove(
        vector<size_t> &var_idxs,
        vector<Integer> &deltas);
    void UpdateWeight();
    void SmoothWeight();
    void ApplyMove(
        size_t var_idx,
        Integer delta);
    double TightScore(
        const ModelVar &var,
        Integer delta);
    bool TightDelta(
        LocalCon &con,
        const ModelCon &common_con,
        size_t i,
        Integer &res);
};