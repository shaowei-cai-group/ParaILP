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
#include "utils/header.h"
#include "utils/paras.h"
#include "workers/Worker.h"
#include "share/Population.h"
#include "share/ModelVar.h"
#include "share/ModelCon.h"
#include "solvers/LocalVar.h"
#include "solvers/local_ilp.h"
class Worker;
class Population;
class Solution;
class ModelVarUtil;
class ModelVar;
class ModelConUtil;
class ModelCon;
class LocalVar;
class Local_ILP;

class Master
{
public:
	char *fileName; //*.mps file
	string modelName;
	vector<Worker *> workerSet;
	Population population;
	size_t sharingPeriod;
	Integer bestObj;
	Integer optValue; // hard instance opt value
	ModelVarUtil modelVarUtil;
	ModelConUtil modelConUtil;
	vector<double> lhsAbsSum;
	vector<double> Polarity;
	std::chrono::_V2::system_clock::time_point clk_st_global;
	unordered_set<Integer> submitObj;

	Master();
	~Master();
	void InitWorkerSet();
	void Run();
	void Solve();
	void TerminateWorker();
	void ReadMPS();

	void PopulationUpdate(
			const vector<LocalVar> &var,
			const Integer obj,
			const Config &config);
	bool PopulationSharing(
			Local_ILP *solver);
	void BestSolUpdate(const Integer obj);
	bool EnterPopulation(const Integer obj);

private:
	inline void IssSetup(
			istringstream &iss,
			const string &line);
	void PushCoeffVarIdx(
			const size_t con_idx,
			const Integer coeff,
			const string &var_name);
	void TightenBound();
	void TightenBoundVar(
			ModelCon &modelCon);
	bool TightBoundGlobally();
	bool SetVarType();
	size_t deleteVarNum = 0;
	size_t inferVarNum = 0;
	size_t deleteConNum = 0;
};
