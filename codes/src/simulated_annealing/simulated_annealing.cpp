#include <thread>
#include "src/simulated_annealing/simulated_annealing.h"
#include "src/general.h"
#include "src/timer.h"
#include "src/local_searches/local_searches.h"

SimulatedAnnealing::SimulatedAnnealing(Instance &instance, std::string algo, std::string folder, std::string file_name, double initial_temperature, double temperature_decrease_rate)
{
    current_temparature_ = initial_temperature;
    temperature_decrease_rate_ = temperature_decrease_rate;
    curr_instance_ = &instance;
    const Graph *graph = instance.graph();
    int num_vertices = graph->num_vertices(), num_arcs = graph->num_arcs();
    int num_routes = instance.num_vehicles();

    // std::cout << "read from file " << file_name << std::endl;
    ALNSHeuristicSolution *sol = new ALNSHeuristicSolution(num_vertices, num_arcs, num_routes);
    sol->ReadFromFile(instance, algo, folder, file_name);
    best_solution_ = sol;
    best_solution_->BuildBitset(instance);
    time_spent_generating_initial_solution_ = sol->total_time_spent_;
    // std::cout << *sol << std::endl;
}

SimulatedAnnealing::SimulatedAnnealing(Instance &instance, HeuristicSolution *initial_sol, double initial_temperature, double temperature_decrease_rate)
{
    current_temparature_ = initial_temperature;
    temperature_decrease_rate_ = temperature_decrease_rate;
    curr_instance_ = &instance;
    // std::cout << "read from file " << file_name << std::endl;
    best_solution_ = new ALNSHeuristicSolution(initial_sol);
    best_solution_->BuildBitset(instance);
    time_spent_generating_initial_solution_ = initial_sol->total_time_spent_;
    // std::cout << *sol << std::endl;
}

ALNSHeuristicSolution *SimulatedAnnealing::best_solution()
{
    return best_solution_;
}

void SimulatedAnnealing::RunOneThread(int num_thread, ALNSHeuristicSolution *initial_solution)
{
    bool converged = false, allow_worse_solution = false;
    ALNSHeuristicSolution *possibly_improved_solution = nullptr, *current_solution = new ALNSHeuristicSolution(initial_solution);

    while (!converged)
    {
        possibly_improved_solution = new ALNSHeuristicSolution(current_solution);
        LocalSearches::RemoveVerticesFromSolution(*curr_instance_, possibly_improved_solution, K_ALNS_PERTURBATION_PERCENTAGE);

        do
        {

            bool continue_search = false;
            do
            {
                continue_search = false;
                if (LocalSearches::DoLocalSearchImprovements(*curr_instance_, possibly_improved_solution))
                {
                    continue_search = true;
                    LocalSearches::TryToInsertUnvisitedVertices(*curr_instance_, possibly_improved_solution);
                }

                if (LocalSearches::DoReplacementImprovements(*curr_instance_, possibly_improved_solution))
                {
                    continue_search = true;
                    LocalSearches::TryToInsertUnvisitedVertices(*curr_instance_, possibly_improved_solution);
                }
            } while (continue_search);

            // IMPORTANT: ShiftingAndInsertion might lead to a worse solution, so we should check BEFORE if the current solution is better than the best found so far as to not lose track of it!
            CheckUpdateBestSolution(possibly_improved_solution);
        } while (LocalSearches::ShiftingAndInsertion(*curr_instance_, possibly_improved_solution));

        (mutex_).lock();
        // std::cout << "temp: " << current_temparature_ << std::endl;
        // std::cout << "curr profits: " << current_solution->profits_sum_ << std::endl;
        // std::cout << "new profits: " << possibly_improved_solution->profits_sum_ << std::endl;
        auto curr_percentage = 100.0 * exp((possibly_improved_solution->profits_sum_ - current_solution->profits_sum_) / current_temparature_);
        auto sorted = 1.0 * (rand() % 101);
        // std::cout << "percentage: " << curr_percentage << std::endl;
        // std::cout << "sorted: " << sorted << std::endl;
        allow_worse_solution = double_greater(current_temparature_, 0.0) ? !double_greater(sorted, curr_percentage) : false;
        (mutex_).unlock();

        if ((allow_worse_solution) || double_greater(possibly_improved_solution->profits_sum_, current_solution->profits_sum_))
        {
            delete current_solution;
            current_solution = possibly_improved_solution;
            CheckUpdateBestSolution(current_solution);
            // getchar();
            // getchar();
        }
        else
        {
            delete possibly_improved_solution;
            possibly_improved_solution = nullptr;
        }

        (mutex_).lock();
        ++total_iter_;
        current_temparature_ *= temperature_decrease_rate_;
        if (double_equals(current_temparature_, 0.0))
            converged = true;
        (mutex_).unlock();
    }

    // std::cout << "* " << best_solution()->profits_sum_ << std::endl;

    // delete tread's solution if not the best found.
    delete current_solution;
    current_solution = nullptr;
}

void SimulatedAnnealing::CheckUpdateBestSolution(ALNSHeuristicSolution *current_solution)
{
    (mutex_).lock();
    if (double_greater(current_solution->profits_sum_, (best_solution())->profits_sum_))
    {
        delete best_solution_;
        best_solution_ = new ALNSHeuristicSolution(current_solution);
        best_solution_->BuildBitset(*curr_instance_);
        last_improve_iteration_ = total_iter_;
    }

    (mutex_).unlock();
}

void SimulatedAnnealing::Run(bool multithreading)
{
    Timestamp *ti = NewTimestamp();
    Timer *timer = GetTimer();
    timer->Clock(ti);
    int initial_solution_profits_sum = 0;
    total_iter_ = 0;
    ALNSHeuristicSolution *curr_sol = nullptr;
    curr_sol = best_solution();

    if (!(curr_sol->is_infeasible_) && (curr_sol->is_feasible_))
    {
        // std::cout << *curr_sol << std::endl;
        // std::cout << curr_sol->profits_sum_ << std::endl;
        LocalSearches::TryToInsertUnvisitedVertices(*curr_instance_, curr_sol, 0);
        // std::cout << curr_sol->profits_sum_ << std::endl;
        initial_solution_profits_sum = curr_sol->profits_sum_;
        // std::cout << *curr_sol << std::endl;

        bool continue_search = false;
        do
        {
            continue_search = false;
            // std::cout << *curr_sol << std::endl;
            // std::cout << "tentou local search" << std::endl;
            if (LocalSearches::DoLocalSearchImprovements(*curr_instance_, curr_sol))
            {
                // std::cout << "fez local search" << std::endl;
                // std::cout << *curr_sol << std::endl;
                continue_search = true;
                LocalSearches::TryToInsertUnvisitedVertices(*curr_instance_, curr_sol, 0);
                // std::cout << *curr_sol << std::endl;
                // std::cout << "after inserts 1" << std::endl;
                // curr_sol->BuildBitset(*(curr_instance_));
                // curr_sol->CheckCorrectness(*curr_instance_);
            }

            // std::cout << "after local searches" << std::endl;
            // // // std::cout << *curr_sol << std::endl;
            // curr_sol->BuildBitset(*(curr_instance_));
            // curr_sol->CheckCorrectness(*curr_instance_);

            // std::cout << "tentou replacements" << std::endl;

            if (LocalSearches::DoReplacementImprovements(*curr_instance_, curr_sol))
            {
                // std::cout << "fez replacements" << std::endl;
                continue_search = true;
                LocalSearches::TryToInsertUnvisitedVertices(*curr_instance_, curr_sol, 0);
            }

            // std::cout << "after replacements" << std::endl;
            // curr_sol->BuildBitset(*(curr_instance_));
            // curr_sol->CheckCorrectness(*curr_instance_);
        } while (continue_search);
        // std::cout << "saiu do loop" << std::endl;
        // curr_sol->BuildBitset(*(curr_instance_));
        // curr_sol->CheckCorrectness(*curr_instance_);

        // std::cout << "FIM" << std::endl;
        // getchar(); getchar();

        // std::cout << curr_sol->profits_sum_ << std::endl;
        // update bitset.
        curr_sol->BuildBitset(*(curr_instance_));
        // std::cout << curr_sol->profits_sum_ << " ";

        // if(curr_sol->CheckCorrectness(*(curr_instance_)) == false){ std::cout << "deu merda" << std::endl; getchar(); getchar();}

        // std::cout << *curr_sol << std::endl;

        int num_cores = (int)std::thread::hardware_concurrency();
        // std::cout << "cores: " << num_cores << std::endl;

        if ((multithreading) && (num_cores > 1))
        {
            std::vector<std::thread> v(num_cores - 1);

            for (size_t i = 0; i < v.size(); ++i)
            {
                v[i] = std::thread(&SimulatedAnnealing::RunOneThread, this, i, curr_sol);
            }

            RunOneThread(num_cores - 1, curr_sol);

            for (size_t i = 0; i < v.size(); ++i)
            {
                (v[i]).join();
            }
        }
        else
            RunOneThread(0, curr_sol);

        curr_sol = best_solution();
        /*std::cout << curr_sol->profits_sum_ << std::endl;
        std::cout << "comecou matheuristica" << std::endl;
        curr_sol = BuildBestSolutionFromGraphInducedByPool();
        std::cout << "* " << curr_sol->profits_sum_ << std::endl;*/
        curr_sol->initial_solution_profits_sum_ = initial_solution_profits_sum;
        curr_sol->last_improve_iteration_ = last_improve_iteration_;
        curr_sol->num_iterations_ = total_iter_;
        curr_sol->total_time_spent_ = time_spent_generating_initial_solution_ + timer->CurrentElapsedTime(ti);
        // std::cout << *curr_sol << std::endl;
        // if(curr_sol->CheckCorrectness(*(curr_instance_)) == false){ std::cout << "deu merda" << std::endl; getchar(); getchar();}
        // std::cout << curr_sol->profits_sum_ << " " << std::endl;
        // std::cout << "last iter improvement: " << last_improve_iteration_ << std::endl;
        //(pool_)[pos_best_sol_]->WriteSolutionFile(curr_instance_->curr_demand(),"teste_improved");
        // std::cout << " ** TERMINOU!" << std::endl;
        // getchar();getchar();

        // delete curr_sol;
        // curr_sol = nullptr;
    }
    delete ti;
    ti = nullptr;
}