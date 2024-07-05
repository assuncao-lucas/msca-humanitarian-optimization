#pragma once

#include "src/heuristic_solution.h"
#include "src/instance.h"

namespace InitalSolutionGenerator
{
    HeuristicSolution *GenerateInitialSolution(Instance &inst);
}