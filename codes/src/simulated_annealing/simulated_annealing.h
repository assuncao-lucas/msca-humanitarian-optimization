#pragma once

#include <mutex>
#include "src/instance.h"
#include "src/heuristic_solution.h"

class SimulatedAnnealing
{
public:
    explicit SimulatedAnnealing() = delete;
    explicit SimulatedAnnealing(Instance &instance, std::string algo, std::string folder, std::string file_name, double initial_temperature, double temperature_decrease_rate);
    explicit SimulatedAnnealing(Instance &instance, HeuristicSolution *initial_sol, double initial_temperature, double temperature_decrease_rate);
    ~SimulatedAnnealing() = default;
    // void Reset();
    void Run(bool multithreading);
    ALNSHeuristicSolution *best_solution();

private:
    void CheckUpdateBestSolution(ALNSHeuristicSolution *current_solution);
    void RunOneThread(int num_thread, ALNSHeuristicSolution *initial_sol);
    Instance *curr_instance_;
    int last_improve_iteration_ = 0;
    int total_iter_ = 0;
    double time_spent_generating_initial_solution_ = 0.0;
    double temperature_decrease_rate_ = 0.0;
    double current_temparature_ = 0.0;
    ALNSHeuristicSolution *best_solution_ = nullptr;

    std::mutex mutex_;
};