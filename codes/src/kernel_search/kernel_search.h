#pragma once

#include "src/heuristic_solution.h"
#include "src/benders_generic_callback.h"
#include <vector>
#include <list>
#include <boost/dynamic_bitset.hpp>

class Instance;

class KernelSearch
{
public:
    explicit KernelSearch(Instance &instance);
    virtual ~KernelSearch();
    void Run();

private:
    IloEnv env_;     // Cplex environment
    IloCplex cplex_; // Cplex solver
    IloModel model_; // Cplex model

    IloNumVarArray f_;
    IloNumVar slack_;
    MasterVariables master_vars_;

    IloNumArray curr_x_values_;
    IloNumArray curr_y_values_;
    IloNumVarArray curr_mip_start_vars_;
    IloNumArray curr_mip_start_vals_;

    boost::dynamic_bitset<> curr_int_y_;
    boost::dynamic_bitset<> curr_int_x_;
    double curr_best_solution_value_ = -1;

    double *R0_ = nullptr;
    double *Rn_ = nullptr;

    const Instance &instance_;
    FPHeuristicSolution solution_;

    void CreateKernelAndBuckets();
    void InitCplex();
    void ResetCplex();
    void BuildModel(bool linearly_relaxed, bool disable_all_binary_vars, bool export_model);
    void RetrieveSolutionArcVertexValues();
    void BuildHeuristicSolution();
    void PrintKernelAndBuckets();
    void UpdateModelVarBounds(boost::dynamic_bitset<> &vars_entering_kernel, boost::dynamic_bitset<> &vars_leaving_kernel, boost::dynamic_bitset<> &curr_reference_kernel);

    boost::dynamic_bitset<> curr_kernel_bitset_;
    std::vector<boost::dynamic_bitset<>> buckets_bitsets_;

    bool found_int_x_ = false;
};