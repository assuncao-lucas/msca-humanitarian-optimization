#include "src/initial_solution/initial_solution.h"
#include "src/kernel_search/kernel_search.h"
#include "src/general.h"
#include <boost/dynamic_bitset.hpp>

namespace InitalSolutionGenerator
{
    HeuristicSolution *GenerateInitialSolution(Instance &inst, bool multithreading)
    {
        Timestamp *ti = NewTimestamp();
        Timer *timer = GetTimer();
        timer->Clock(ti);

        auto formulation = Formulation::single_commodity;
        KernelSearch ks(inst);
        auto graph = inst.graph();
        auto num_vertices = graph->num_vertices();
        auto num_mandatory = inst.num_mandatory();
        KSHeuristicSolution *solution = new KSHeuristicSolution(num_vertices, graph->num_arcs(), inst.num_vehicles());

        if (num_mandatory == 0) // if no mandatory vertex, return solution with empty routes.
        {
            solution->is_feasible_ = solution->found_x_integer_ = true;
            ks.curr_best_solution_value_ = 0.0;
            ks.BuildHeuristicSolution(solution);
        }
        else
        {
            ks.InitCplex(multithreading);
            ks.BuildModel(formulation, false, true, false);

            // enable just origin and mandatory vertices!
            auto active_vertices = boost::dynamic_bitset<>(num_vertices, 0);
            for (int i = 0; i <= num_mandatory; ++i)
                active_vertices[i] = 1;

            auto empty_bitset = boost::dynamic_bitset<>(num_vertices, 0);
            ks.UpdateModelVarBounds(active_vertices, empty_bitset, active_vertices);

            if (ks.cplex_->solve())
            {
                ks.RetrieveSolutionArcVertexValues();
                solution->is_feasible_ = solution->found_x_integer_ = true;
                ks.curr_best_solution_value_ = ks.cplex_->getObjValue();
                ks.BuildHeuristicSolution(solution);
            }
            else if ((ks.cplex_->getCplexStatus() == IloCplex::Infeasible) || (ks.cplex_->getCplexStatus() == IloCplex::InfOrUnbd))
            {
                solution->is_infeasible_ = true;
            }
        }

        solution->total_time_spent_ = timer->CurrentElapsedTime(ti);
        delete ti;
        ti = nullptr;

        return solution;
    }
}