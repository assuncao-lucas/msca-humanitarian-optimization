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
	std::vector<std::pair<int, int>> ordered_profits_;
	std::vector<ALNSHeuristicSolution *> pool_;
	int num_elements_in_pool_;
	int pos_best_sol_;
	int pos_worst_sol_;
	int last_improve_iteration_;
	int non_improve_iterations_counter_;
	std::mutex mutex_;

	void RemoveVerticesFromSolution(ALNSHeuristicSolution *solution, double percentage);
	bool TryToInsertUnvisitedVertices(ALNSHeuristicSolution *solution, int type = -1);
	bool TryToInsertUnvisitedVerticesOneRoute(ALNSHeuristicSolution *solution, int curr_route);
	bool DoLocalSearchImprovements(ALNSHeuristicSolution *solution);

	bool ShiftingAndInsertion(ALNSHeuristicSolution *solution);
	bool ShiftingOneVertex(ALNSHeuristicSolution *sol, int vertex, double &global_variation);

	bool DoAllInterRoutesImprovements(ALNSHeuristicSolution *solution);

	bool Do_1_1_Improvement(ALNSHeuristicSolution *solution, int i1, int i2);
	bool Do_1_0_Improvement(ALNSHeuristicSolution *solution, int i1, int i2);
	bool Do_2_1_Improvement(ALNSHeuristicSolution *solution, int i1, int i2);

	bool Do_1_1_Unrouted_Improvements(ALNSHeuristicSolution *solution);
	bool Do_2_1_Unrouted_Improvements(ALNSHeuristicSolution *solution);

	bool DoIntraRouteImprovementOneRoute(ALNSHeuristicSolution *sol, int route);

	bool DoAllIntraRouteImprovements(ALNSHeuristicSolution *solution);

	bool DoReplacementImprovements(ALNSHeuristicSolution *sol);

	bool AddSolutionToPool(ALNSHeuristicSolution *sol, int iter);
	ALNSHeuristicSolution *CopyRandomSolutionFromPool();

	ALNSHeuristicSolution *DoPathRelinking(ALNSHeuristicSolution *sol);
	ALNSHeuristicSolution *DoPathRelinkingIter(ALNSHeuristicSolution *guiding_sol, ALNSHeuristicSolution *new_sol);

	bool CheckPoolIntegrity();
	ALNSHeuristicSolution *BuildBestSolutionFromPool();
	// ALNSHeuristicSolution *BuildBestSolutionFromGraphInducedByPool();
};