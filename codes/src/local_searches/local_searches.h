#pragma once

#include "src/heuristic_solution.h"
#include "src/instance.h"

namespace LocalSearches
{
    void RemoveVerticesFromSolution(Instance &instance, HeuristicSolution *solution, double percentage);
    bool TryToInsertUnvisitedVertices(Instance &instance, HeuristicSolution *solution, int type = -1);
    bool TryToInsertUnvisitedVerticesOneRoute(Instance &instance, HeuristicSolution *solution, int curr_route);
    bool DoLocalSearchImprovements(Instance &instance, HeuristicSolution *solution);

    bool ShiftingAndInsertion(Instance &instance, HeuristicSolution *solution);
    bool ShiftingOneVertex(Instance &instance, HeuristicSolution *sol, int vertex, double &global_variation);

    bool DoAllInterRoutesImprovements(Instance &instance, HeuristicSolution *solution);

    bool Do_1_1_Improvement(Instance &instance, HeuristicSolution *solution, int i1, int i2);
    bool Do_1_0_Improvement(Instance &instance, HeuristicSolution *solution, int i1, int i2);
    bool Do_2_1_Improvement(Instance &instance, HeuristicSolution *solution, int i1, int i2);

    bool Do_1_1_Unrouted_Improvements(Instance &instance, HeuristicSolution *solution);
    bool Do_2_1_Unrouted_Improvements(Instance &instance, HeuristicSolution *solution);

    bool DoIntraRouteImprovementOneRoute(Instance &instance, HeuristicSolution *sol, int route);

    bool DoAllIntraRouteImprovements(Instance &instance, HeuristicSolution *solution);

    bool DoReplacementImprovements(Instance &instance, HeuristicSolution *sol);
}