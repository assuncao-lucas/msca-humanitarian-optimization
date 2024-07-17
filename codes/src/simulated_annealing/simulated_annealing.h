#pragma once

#include <mutex>
#include "src/instance.h"
#include "src/heuristic_solution.h"

class SimulatedAnnealing
{
public:
    explicit SimulatedAnnealing() = delete;
    explicit SimulatedAnnealing(Instance &instance, std::string algo, std::string folder, std::string file_name);
    explicit SimulatedAnnealing(Instance &instance, HeuristicSolution *initial_sol);
    ~SimulatedAnnealing();
    // void Reset();
    void Run(double temperature_decrease_rate, bool multithreading);
    MetaHeuristicSolution *best_solution();

private:
    void CheckUpdateBestSolution(MetaHeuristicSolution *current_solution);
    void RunOneThread(int num_thread, MetaHeuristicSolution *initial_sol);
    MetaHeuristicSolution *RunOneStep(MetaHeuristicSolution *current_solution);
    void ComputeAndSetInitialTemperature(int sampling_size, double target_acceptance_probability);
    Instance *curr_instance_;
    int last_improve_iteration_ = 0;
    int total_iter_ = 0;
    double time_spent_generating_initial_solution_ = 0.0;
    double temperature_decrease_rate_ = 0.0;
    double current_temperature_ = 0.0;
    MetaHeuristicSolution *best_solution_ = nullptr;

    std::mutex mutex_;
};