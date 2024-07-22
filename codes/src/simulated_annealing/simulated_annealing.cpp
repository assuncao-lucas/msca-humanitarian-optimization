#include <thread>
#include "src/simulated_annealing/simulated_annealing.h"
#include "src/general.h"
#include "src/timer.h"
#include "src/local_searches/local_searches.h"

SimulatedAnnealing::SimulatedAnnealing(Instance &instance, std::string algo, std::string folder, std::string file_name)
{
    curr_instance_ = &instance;
    const Graph *graph = instance.graph();
    int num_vertices = graph->num_vertices(), num_arcs = graph->num_arcs();
    int num_routes = instance.num_vehicles();

    // std::cout << "read from file " << file_name << std::endl;
    MetaHeuristicSolution *sol = new MetaHeuristicSolution(num_vertices, num_arcs, num_routes);
    sol->ReadFromFile(instance, algo, folder, file_name);
    best_solution_ = sol;
    best_solution_->BuildBitset(instance);
    time_spent_generating_initial_solution_ = sol->total_time_spent_;
    // std::cout << *sol << std::endl;
}

SimulatedAnnealing::SimulatedAnnealing(Instance &instance, HeuristicSolution *initial_sol)
{
    curr_instance_ = &instance;
    // std::cout << "read from file " << file_name << std::endl;
    best_solution_ = new MetaHeuristicSolution(initial_sol);
    best_solution_->BuildBitset(instance);
    time_spent_generating_initial_solution_ = initial_sol->total_time_spent_;
    // std::cout << *sol << std::endl;
}

SimulatedAnnealing::~SimulatedAnnealing()
{
    if (best_solution_)
    {
        delete best_solution_;
        best_solution_ = nullptr;
    }
}

MetaHeuristicSolution *SimulatedAnnealing::best_solution()
{
    return best_solution_;
}

MetaHeuristicSolution *SimulatedAnnealing::RunOneStep(MetaHeuristicSolution *current_solution)
{
    MetaHeuristicSolution *possibly_improved_solution = new MetaHeuristicSolution(current_solution);
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

    return possibly_improved_solution;
}

void SimulatedAnnealing::ComputeAndSetInitialTemperature(size_t sampling_size, double target_acceptance_probability)
{
    std::vector<std::pair<double, double>> degrading_sampling;
    MetaHeuristicSolution *initial_solution = new MetaHeuristicSolution(best_solution_), *possibly_degraded_solution = nullptr;

    // Do initial sampling.
    bool converged_sampling = false;
    size_t count = 0, sampling_iter = 0;
    double previous_profit_sum = 0.0;
    do
    {
        previous_profit_sum = initial_solution->profits_sum_;
        possibly_degraded_solution = RunOneStep(initial_solution);

        // if step degrades solution quality/cost, add it to sampling.
        if (double_less(possibly_degraded_solution->profits_sum_, previous_profit_sum))
        {
            degrading_sampling.push_back({possibly_degraded_solution->profits_sum_, previous_profit_sum});
            ++count;
        }
        else
        {
            // use the sampling procedure to also explore improvement of the solution.
            CheckUpdateBestSolution(possibly_degraded_solution);
        }

        ++sampling_iter;
        // std::cout << sampling_iter << " " << count << std::endl;
        ++total_iter_; // already counts as an iteration of the SA.
        if (sampling_iter >= sampling_size)
            converged_sampling = true;

        delete initial_solution;
        initial_solution = possibly_degraded_solution;
    } while (!converged_sampling);

    delete possibly_degraded_solution;
    possibly_degraded_solution = nullptr;

    // compute input temperature.
    double input_temp = K_SA_DEFAULT_INITIAL_TEMPERATURE;
    if (count >= K_SA_MIN_SAMPLING)
    {
        double sum_deltas = 0.0;
        for (const auto &[min, max] : degrading_sampling)
            sum_deltas += (max - min);
        input_temp = -sum_deltas / (degrading_sampling.size() * log(target_acceptance_probability));
    }

    // std::cout << "input temp: " << input_temp << std::endl;

    bool converged = false;
    double convergence_tolerance = 0.001, current_acceptance_probability = 0.0;
    double sum_min = 0.0, sum_max = 0.0;
    double current_temp = input_temp, previous_temp = 0.0;
    double temp_delta = 0.0, previous_temp_delta = 0.0;
    double p_factor = 1.0;

    // compute the ideal initial temperature based on sampling and input temperature.
    do
    {
        // std::cout << "current_temp: " << current_temp << std::endl;
        sum_min = 0.0;
        sum_max = 0.0;
        for (const auto &[min, max] : degrading_sampling)
        {
            sum_min += exp(-min / current_temp);
            sum_max += exp(-max / current_temp);
        }
        current_acceptance_probability = sum_max / sum_min;
        // std::cout << "current_acceptance_probability: " << current_acceptance_probability << std::endl;

        if (!double_greater(fabs(current_acceptance_probability - target_acceptance_probability), convergence_tolerance))
        {
            converged = true;
        }
        else
        {
            previous_temp = current_temp;
            current_temp *= pow(log(current_acceptance_probability) / log(target_acceptance_probability), 1.0 / p_factor);
            previous_temp_delta = temp_delta;
            temp_delta = current_temp - previous_temp;

            // check oscilation in temperature increase and adjust parameter p to gaarantee convergence.
            if (double_less(temp_delta * previous_temp_delta, 0.0))
            {
                p_factor *= 2.0;
                // std::cout << "fez ajuste!" << std::endl;
            }
        }

    } while (!converged);

    // std::cout << "initial temp: " << current_temp << std::endl;
    current_temperature_ = current_temp;
}

void SimulatedAnnealing::RunOneThread(int num_thread, MetaHeuristicSolution *initial_solution)
{
    bool converged = false, allow_worse_solution = false;
    MetaHeuristicSolution *possibly_improved_solution = nullptr, *current_solution = new MetaHeuristicSolution(initial_solution);
    double curr_percentage = 0.0, sorted = 0.0;
    while (!converged)
    {
        possibly_improved_solution = RunOneStep(current_solution);

        (mutex_).lock();
        // std::cout << "temp: " << current_temperature_ << std::endl;
        // std::cout << "curr profits: " << current_solution->profits_sum_ << std::endl;
        // std::cout << "new profits: " << possibly_improved_solution->profits_sum_ << std::endl;
        curr_percentage = 100.0 * exp((possibly_improved_solution->profits_sum_ - current_solution->profits_sum_) / current_temperature_);
        sorted = 1.0 * (rand() % 101);
        // std::cout << "percentage: " << curr_percentage << std::endl;
        // std::cout << "sorted: " << sorted << std::endl;
        allow_worse_solution = double_greater(current_temperature_, 0.0) ? !double_greater(sorted, curr_percentage) : false;
        (mutex_).unlock();

        if ((allow_worse_solution) || double_greater(possibly_improved_solution->profits_sum_, current_solution->profits_sum_))
        {
            delete current_solution;
            current_solution = possibly_improved_solution;
            CheckUpdateBestSolution(current_solution);
        }
        else
        {
            delete possibly_improved_solution;
            possibly_improved_solution = nullptr;
        }

        (mutex_).lock();
        ++total_iter_;
        // std::cout << total_iter_ << std::endl;
        // std::cout << current_temperature_ << std::endl;
        current_temperature_ *= temperature_decrease_rate_;
        if (double_equals(current_temperature_, 0.0))
            converged = true;
        (mutex_).unlock();
    }

    // std::cout << "* " << best_solution()->profits_sum_ << std::endl;

    // delete tread's solution if not the best found.
    delete current_solution;
    current_solution = nullptr;
}

void SimulatedAnnealing::CheckUpdateBestSolution(MetaHeuristicSolution *current_solution)
{
    (mutex_).lock();
    if (double_greater(current_solution->profits_sum_, (best_solution())->profits_sum_))
    {
        delete best_solution_;
        best_solution_ = new MetaHeuristicSolution(current_solution);
        best_solution_->BuildBitset(*curr_instance_);
        last_improve_iteration_ = total_iter_;
    }

    (mutex_).unlock();
}

void SimulatedAnnealing::Run(double temperature_decrease_rate, size_t sampling_size, double target_acceptance_probability, bool multithreading)
{
    temperature_decrease_rate_ = temperature_decrease_rate;
    Timestamp *ti = NewTimestamp();
    Timer *timer = GetTimer();
    timer->Clock(ti);
    int initial_solution_profits_sum = 0;
    total_iter_ = 0;
    MetaHeuristicSolution *curr_sol = best_solution_;

    if (!(curr_sol->is_infeasible_) && (curr_sol->is_feasible_))
    {
        LocalSearches::TryToInsertUnvisitedVertices(*curr_instance_, curr_sol, 0);
        initial_solution_profits_sum = curr_sol->profits_sum_;

        bool continue_search = false;
        do
        {
            continue_search = false;
            if (LocalSearches::DoLocalSearchImprovements(*curr_instance_, curr_sol))
            {
                continue_search = true;
                LocalSearches::TryToInsertUnvisitedVertices(*curr_instance_, curr_sol, 0);
            }

            if (LocalSearches::DoReplacementImprovements(*curr_instance_, curr_sol))
            {
                continue_search = true;
                LocalSearches::TryToInsertUnvisitedVertices(*curr_instance_, curr_sol, 0);
            }
        } while (continue_search);

        curr_sol->BuildBitset(*(curr_instance_));

        ComputeAndSetInitialTemperature(sampling_size, target_acceptance_probability);
        // NOTE: at this point, the curr_sol might have been deleted, because ComputeAndSetInitialTemperature might update best_solution!!
        // this is why, either in single or multithread, the best solution is retrieved once again!
        // current_temperature_ = 100;

        int num_cores = (int)std::thread::hardware_concurrency();
        // std::cout << "cores: " << num_cores << std::endl;

        if ((multithreading) && (num_cores > 1))
        {
            curr_sol = new MetaHeuristicSolution(best_solution()); // need to create copy because, in multithreading, when the best solution is
                                                                   // updated before the init of another thread, the reference to this initial best solution would already be lost.
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
            delete curr_sol;
        }
        else
            RunOneThread(0, best_solution());

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