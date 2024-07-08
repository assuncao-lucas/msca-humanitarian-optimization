#pragma once

#include <unordered_map>
#include <vector>
#include <mutex>
#include "src/instance.h"
#include "src/timer.h"
#include "src/heuristic_solution.h"

class ALNS
{
public:
	ALNS();
	~ALNS();
	void Init(Instance &instance, std::string algo, std::string folder, std::string file_name);
	void Init(Instance &instance, HeuristicSolution *initial_sol);
	void Reset();
	void Run();
	void RunOneThread(int num_thread, int num_iterations);
	void PrintPool();
	ALNSHeuristicSolution *best_solution();
	ALNSHeuristicSolution *worst_solution();

private:
	Instance *curr_instance_;
	std::vector<ALNSHeuristicSolution *> pool_;
	int num_elements_in_pool_;
	int pos_best_sol_;
	int pos_worst_sol_;
	int last_improve_iteration_;
	int non_improve_iterations_counter_;
	double time_spent_generating_initial_solution_ = 0.0;
	std::mutex mutex_;

	bool AddSolutionToPool(ALNSHeuristicSolution *sol, int iter);
	ALNSHeuristicSolution *CopyRandomSolutionFromPool();

	ALNSHeuristicSolution *DoPathRelinking(ALNSHeuristicSolution *sol);
	ALNSHeuristicSolution *DoPathRelinkingIter(ALNSHeuristicSolution *guiding_sol, ALNSHeuristicSolution *new_sol);

	bool CheckPoolIntegrity();
	ALNSHeuristicSolution *BuildBestSolutionFromPool();
	// ALNSHeuristicSolution *BuildBestSolutionFromGraphInducedByPool();
};