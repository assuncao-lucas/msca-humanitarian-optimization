#pragma once

#include <unordered_map>
#include <vector>
#include <mutex>
#include "src/instance.h"
#include "src/heuristic_solution.h"

class ALNS
{
public:
	explicit ALNS() = delete;
	explicit ALNS(Instance &instance, std::string algo, std::string folder, std::string file_name, int max_pool_size);
	explicit ALNS(Instance &instance, HeuristicSolution *initial_sol, int max_pool_size);
	~ALNS();
	// void Reset();
	void Run(int num_iterations, bool multithreading);
	void PrintPool();
	MetaHeuristicSolution *best_solution();
	MetaHeuristicSolution *worst_solution();

private:
	void RunOneThread(int num_thread, int num_iterations);
	Instance *curr_instance_;
	std::vector<MetaHeuristicSolution *> pool_;
	int num_elements_in_pool_ = 0;
	int pos_best_sol_ = -1;
	int pos_worst_sol_ = -1;
	int last_improve_iteration_ = 0;
	int non_improve_iterations_counter_ = 0;
	int max_pool_size_ = 0;
	int total_iter_ = 0;
	double time_spent_generating_initial_solution_ = 0.0;

	std::mutex mutex_;

	bool AddSolutionToPool(MetaHeuristicSolution *sol);
	MetaHeuristicSolution *CopyRandomSolutionFromPool();

	MetaHeuristicSolution *DoPathRelinking(MetaHeuristicSolution *sol);
	MetaHeuristicSolution *DoPathRelinkingIter(MetaHeuristicSolution *guiding_sol, MetaHeuristicSolution *new_sol);

	bool CheckPoolIntegrity();
	MetaHeuristicSolution *BuildBestSolutionFromPool();
	// MetaHeuristicSolution *BuildBestSolutionFromGraphInducedByPool();
};