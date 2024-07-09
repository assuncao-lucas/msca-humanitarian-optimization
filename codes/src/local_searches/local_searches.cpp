
#include "src/local_searches/local_searches.h"

namespace LocalSearches
{

    bool DoReplacementImprovements(Instance &instasnce, HeuristicSolution *sol)
    {
        bool global_improvement = false;
        bool improved = false;

        do
        {
            improved = (Do_1_1_Unrouted_Improvements(instasnce, sol) || Do_2_1_Unrouted_Improvements(instasnce, sol));
            global_improvement = (global_improvement || improved);
        } while (improved);

        return global_improvement;
    }

    bool DoLocalSearchImprovements(Instance &instance, HeuristicSolution *sol)
    {
        bool global_improvement = false;
        bool improved = false;

        do
        {
            improved = DoAllInterRoutesImprovements(instance, sol);
            global_improvement = (global_improvement || improved);

            // std::cout << "after inter route" << std::endl;
            // sol->BuildBitset(instance);
            // sol->CheckCorrectness(instance);

            improved = DoAllIntraRouteImprovements(instance, sol);
            global_improvement = (global_improvement || improved);
        } while (improved);

        return global_improvement;
    }

    bool DoIntraRouteImprovementOneRoute(Instance &instance, HeuristicSolution *sol, int route)
    {
        bool iteration_improvement = false;
        bool improved = false;

        do
        {
            // std::cout << "start 2 opt" << std::endl;
            improved = sol->Do2OptImprovement(instance, route);
            iteration_improvement = (iteration_improvement || improved);
            // std::cout << "end 2 opt" << std::endl;
            // sol->BuildBitset(instance);
            // sol->CheckCorrectness(instance);
        } while (improved);

        // do
        // {
        //   improved = sol->Do3OptImprovement(graph, route);
        //   iteration_improvement = (iteration_improvement || improved);
        // } while (improved);

        return iteration_improvement;
    }

    bool DoAllIntraRouteImprovements(Instance &instance, HeuristicSolution *solution)
    {
        bool global_improvement = false, improved = false;
        int num_vehicles = instance.num_vehicles();

        for (int i = 0; i < num_vehicles; ++i)
        {
            improved = DoIntraRouteImprovementOneRoute(instance, solution, i);
            global_improvement = (global_improvement || improved);
        }

        return global_improvement;
    }

    bool Do_1_1_Unrouted_Improvements(Instance &instance, HeuristicSolution *solution)
    {
        Route *curr_route = nullptr;
        std::list<int>::iterator it_i1, it_i2, it_f1, prox_it;
        double profit_variation1 = 0.0;
        double time_variation1 = 0.0;
        int pos_i1 = 0, pos_f1 = 0;
        int max_pos1 = 0;
        int num_routes = instance.num_vehicles();
        int curr_vertex = 0;
        VertexStatus *curr_status = nullptr;
        int num_vertices = (instance.graph())->num_vertices(), num_mandatory = instance.num_mandatory();
        std::list<int>::iterator next_vertex_to_be_replaced_it, first_vertex_to_be_replaced_it;

        const auto &ordered_profits = instance.ordered_profits();

        int type = rand() % 2;
        for (size_t i = 0; i < ordered_profits.size(); ++i)
        {
            if (type == 0)
                curr_vertex = (ordered_profits[num_vertices - num_mandatory - 2 - i]).first;
            if (type == 1)
                curr_vertex = (ordered_profits[i]).first;
            curr_status = &((solution->vertex_status_vec_)[curr_vertex]);

            if (!(curr_status->selected_))
            {
                // for(std::list<int>::iterator it_i2 = (solution->unvisited_vertices_).begin(); it_i2 != (solution->unvisited_vertices_).end();)
                //{
                // prox_it = it_i2;
                //++prox_it;
                it_i2 = curr_status->pos_;
                for (int r1 = 0; r1 < num_routes; ++r1)
                {
                    curr_route = &((solution->routes_vec_)[r1]);
                    max_pos1 = (int)((curr_route->vertices_).size()) - 1;
                    next_vertex_to_be_replaced_it = (curr_route->vertices_).begin();
                    for (int i = 0; i <= max_pos1; ++i)
                    {
                        first_vertex_to_be_replaced_it = next_vertex_to_be_replaced_it;
                        ++next_vertex_to_be_replaced_it;
                        pos_i1 = pos_f1 = i;
                        if (*first_vertex_to_be_replaced_it > num_mandatory)
                        {
                            if (solution->PreviewInterRouteSwapUnrouted(instance, r1, pos_i1, pos_f1, it_i1, it_f1, profit_variation1, time_variation1, it_i2))
                            {
                                if (double_greater(profit_variation1, 0.0) || (double_equals(profit_variation1, 0.0) && double_less(time_variation1, 0.0)))
                                {
                                    // std::cout << "before inter route swap unrouted" << std::endl;
                                    // solution->BuildBitset(instance);
                                    // solution->CheckCorrectness(instance);
                                    solution->InterRouteSwapUnrouted(r1, it_i1, it_f1, profit_variation1, time_variation1, it_i2);
                                    // std::cout << "after inter route swap unrouted" << std::endl;
                                    // solution->BuildBitset(instance);
                                    // solution->CheckCorrectness(instance);
                                    return true;
                                }
                            }
                            // std::cout << "after preview inter route unrouted FAILEDA" << std::endl;
                            // solution->BuildBitset(instance);
                            // solution->CheckCorrectness(instance);
                        }
                        // first_vertex_to_be_replaced_it = next_vertex_to_be_replaced_it;
                    }
                }

                // it_i2 = prox_it;
            }
        }

        return false;
    }

    bool Do_2_1_Unrouted_Improvements(Instance &instance, HeuristicSolution *solution)
    {
        Route *curr_route = nullptr;
        std::list<int>::iterator it_i1, it_i2, it_f1, prox_it;
        double profit_variation1 = 0.0;
        double time_variation1 = 0.0;
        int pos_i1 = 0, pos_f1 = 0;
        int max_pos1 = 0;
        int num_routes = instance.num_vehicles();
        int curr_vertex = 0;
        VertexStatus *curr_status = nullptr;
        int num_vertices = (instance.graph())->num_vertices(), num_mandatory = instance.num_mandatory();
        ;
        int type = rand() % 2;
        std::list<int>::iterator next_vertex_to_be_replaced_it, first_vertex_to_be_replaced_it;

        const auto &ordered_profits = instance.ordered_profits();

        for (size_t i = 0; i < ordered_profits.size(); ++i)
        {
            if (type == 0)
                curr_vertex = (ordered_profits[num_vertices - num_mandatory - 2 - i]).first;
            if (type == 1)
                curr_vertex = (ordered_profits[i]).first;
            curr_status = &((solution->vertex_status_vec_)[curr_vertex]);

            if (!(curr_status->selected_))
            {
                // for(std::list<int>::iterator it_i2 = (solution->unvisited_vertices_).begin(); it_i2 != (solution->unvisited_vertices_).end();)
                //{
                // prox_it = it_i2;
                //++prox_it;
                it_i2 = curr_status->pos_;
                for (int r1 = 0; r1 < num_routes; ++r1)
                {
                    curr_route = &((solution->routes_vec_)[r1]);
                    max_pos1 = (int)((curr_route->vertices_).size()) - 1;
                    next_vertex_to_be_replaced_it = (curr_route->vertices_).begin();
                    for (int i = 0; i <= max_pos1 - 1; ++i)
                    {
                        first_vertex_to_be_replaced_it = next_vertex_to_be_replaced_it;
                        ++next_vertex_to_be_replaced_it; // it's never the .end() iterator, since loop only goes until max_pos1-1

                        pos_i1 = i;
                        pos_f1 = i + 1;

                        if ((*first_vertex_to_be_replaced_it > num_mandatory) && (*next_vertex_to_be_replaced_it > num_mandatory))
                        {
                            if (solution->PreviewInterRouteSwapUnrouted(instance, r1, pos_i1, pos_f1, it_i1, it_f1, profit_variation1, time_variation1, it_i2))
                            {
                                if (double_greater(profit_variation1, 0.0) || (double_equals(profit_variation1, 0.0) && double_less(time_variation1, 0.0)))
                                {
                                    // std::cout << "before inter route swap unrouted" << std::endl;
                                    // solution->BuildBitset(instance);
                                    // solution->CheckCorrectness(instance);
                                    solution->InterRouteSwapUnrouted(r1, it_i1, it_f1, profit_variation1, time_variation1, it_i2);
                                    // std::cout << "after inter route swap unrouted" << std::endl;
                                    // solution->BuildBitset(instance);
                                    // solution->CheckCorrectness(instance);
                                    return true;
                                }
                            }
                            // std::cout << "after preview inter route unrouted FAILEDA" << std::endl;
                            // solution->BuildBitset(instance);
                            // solution->CheckCorrectness(instance);
                        }
                    }
                }
                // it_i2 = prox_it;
            }
        }

        return false;
    }

    bool Do_1_1_Improvement(Instance &instance, HeuristicSolution *solution, int r1, int r2)
    {
        Route *route1 = &((solution->routes_vec_)[r1]), *route2 = &((solution->routes_vec_)[r2]);

        if (((route1->vertices_).empty()) || ((route2->vertices_).empty()))
            return false;

        std::list<int>::iterator it_i1, it_i2, it_f1, it_f2;
        double profit_variation1 = 0.0, profit_variation2 = 0.0;
        double time_variation1 = 0.0, time_variation2 = 0.0;
        double total_time_variation = 0.0, total_profit_variation = 0.0;
        int pos_i1 = 0, pos_f1 = 0, pos_i2 = 0, pos_f2 = 0;
        int max_pos1 = (int)((route1->vertices_).size()) - 1, max_pos2 = (int)((route2->vertices_).size()) - 1;

        for (int i = 0; i <= max_pos1; ++i)
        {
            pos_i1 = pos_f1 = i;
            for (int j = 0; j <= max_pos2; ++j)
            {
                pos_i2 = pos_f2 = j;
                // std::cout << *solution << std::endl;
                // std::cout << "preview swap " << pos_i1 << " to " << pos_f1 << " from " << r1 << " with " << pos_i2 << " to " << pos_f2 << " of " << r2 << std::endl;
                if (solution->PreviewInterRouteSwap(instance, r1, pos_i1, pos_f1, r2, pos_i2, pos_f2, it_i1, it_f1, profit_variation1, time_variation1, it_i2, it_f2, profit_variation2, time_variation2))
                {
                    total_time_variation = time_variation1 + time_variation2;
                    total_profit_variation = profit_variation1 + profit_variation2;
                    if (double_greater(total_profit_variation, 0.0) || (double_equals(total_profit_variation, 0.0) && double_less(total_time_variation, 0.0)))
                    {
                        // std::cout << "before inter route swap vertex" << std::endl;
                        // solution->BuildBitset(instance);
                        // solution->CheckCorrectness(instance);
                        solution->InterRouteSwap(r1, r2, it_i1, it_f1, profit_variation1, time_variation1, it_i2, it_f2, profit_variation2, time_variation2);
                        // std::cout << "after inter route swap vertex" << std::endl;
                        // solution->BuildBitset(instance);
                        // solution->CheckCorrectness(instance);
                        return true;
                    }
                }
                // std::cout << "after preview inter route FAILED" << std::endl;
                // // std::cout << *solution << std::endl;
                // solution->BuildBitset(instance);
                // solution->CheckCorrectness(instance);
            }
        }

        return false;
    }

    bool Do_1_0_Improvement(Instance &instance, HeuristicSolution *solution, int r1, int r2)
    {
        Route *route1 = &((solution->routes_vec_)[r1]), *route2 = &((solution->routes_vec_)[r2]);

        if ((route1->vertices_).empty())
            return false;

        std::list<int>::iterator it;
        double profit_variation1 = 0.0, profit_variation2 = 0.0;
        double time_variation1 = 0.0, time_variation2 = 0.0;
        double total_time_variation = 0.0, total_profit_variation = 0.0;
        int max_pos2 = (int)((route2->vertices_).size());
        int vertex = 0;

        auto it_vertex = (route1->vertices_).begin();
        while (it_vertex != (route1->vertices_).end())
        {
            vertex = *it_vertex;
            ++it_vertex; // already increment it as to avoid losing track of the original iterator reference. This is necessary because, at each PreviewInterRoute call, the current vertex is removed from the route and added back within a new iterator.
            for (int pos = 0; pos <= max_pos2; ++pos)
            {
                // std::cout << *solution << std::endl;
                // std::cout << "preview move vertex " << vertex << " to route " << r2 << " pos " << pos << std::endl;
                if (solution->PreviewInterRouteMoveVertex(instance, vertex, r2, pos, it, profit_variation1, time_variation1, profit_variation2, time_variation2))
                {
                    total_time_variation = time_variation1 + time_variation2;
                    total_profit_variation = profit_variation1 + profit_variation2;
                    if (double_greater(total_profit_variation, 0.0) || (double_equals(total_profit_variation, 0.0) && double_less(total_time_variation, 0.0)))
                    {
                        // std::cout << "before inter route move vertex" << std::endl;
                        // solution->BuildBitset(instance);
                        // solution->CheckCorrectness(instance);
                        solution->InterRouteMoveVertex(vertex, r2, it, profit_variation1, time_variation1, profit_variation2, time_variation2);
                        // std::cout << "after inter route move vertex" << std::endl;
                        // solution->BuildBitset(instance);
                        // solution->CheckCorrectness(instance);
                        return true;
                    }
                }
                // std::cout << "after preview inter route FAILED" << std::endl;
                // // std::cout << *solution << std::endl;
                // solution->BuildBitset(instance);
                // solution->CheckCorrectness(instance);
            }
        }

        return false;
    }

    bool Do_2_1_Improvement(Instance &instance, HeuristicSolution *solution, int r1, int r2)
    {
        Route *route1 = &((solution->routes_vec_)[r1]), *route2 = &((solution->routes_vec_)[r2]);

        if (((route1->vertices_).size() < 2) || ((route2->vertices_).empty()))
            return false;

        std::list<int>::iterator it_i1, it_i2, it_f1, it_f2;
        double profit_variation1 = 0.0, profit_variation2 = 0.0;
        double time_variation1 = 0.0, time_variation2 = 0.0;
        double total_time_variation = 0.0, total_profit_variation = 0.0;
        int pos_i1 = 0, pos_f1 = 0, pos_i2 = 0, pos_f2 = 0;
        int max_pos1 = (int)((route1->vertices_).size()) - 1, max_pos2 = (int)((route2->vertices_).size()) - 1;

        for (int i = 0; i <= max_pos1 - 1; ++i)
        {
            pos_i1 = i;
            pos_f1 = i + 1;
            for (int j = 0; j <= max_pos2; ++j)
            {
                pos_i2 = pos_f2 = j;
                // std::cout << "preview swap " << pos_i1 << " to  " << pos_f1 << " from " << r1 << " with " << pos_i2 << " to " << pos_f2 << " of " << r2 << std::endl;
                if (solution->PreviewInterRouteSwap(instance, r1, pos_i1, pos_f1, r2, pos_i2, pos_f2, it_i1, it_f1, profit_variation1, time_variation1, it_i2, it_f2, profit_variation2, time_variation2))
                {
                    total_time_variation = time_variation1 + time_variation2;
                    total_profit_variation = profit_variation1 + profit_variation2;
                    if (double_greater(total_profit_variation, 0.0) || (double_equals(total_profit_variation, 0.0) && double_less(total_time_variation, 0.0)))
                    {
                        // std::cout << "before inter route swap vertex" << std::endl;
                        // solution->BuildBitset(instance);
                        // solution->CheckCorrectness(instance);
                        solution->InterRouteSwap(r1, r2, it_i1, it_f1, profit_variation1, time_variation1, it_i2, it_f2, profit_variation2, time_variation2);
                        // std::cout << "after inter route swap vertex" << std::endl;
                        // solution->BuildBitset(instance);
                        // solution->CheckCorrectness(instance);
                        return true;
                    }
                }
                // std::cout << "after preview inter route FAILED" << std::endl;
                // solution->BuildBitset(instance);
                // solution->CheckCorrectness(instance);
            }
        }

        return false;
    }

    bool DoAllInterRoutesImprovements(Instance &instance, HeuristicSolution *solution)
    {
        bool improved = false;
        bool global_improvement = false;
        int num_routes = instance.num_vehicles();

        do
        {
            improved = false;
            // para todas as combinações entre duas rotas
            for (int i1 = 0; i1 < num_routes; ++i1)
            {
                for (int i2 = 0; i2 < num_routes; ++i2)
                {
                    // Faz inter improvement apenas se as rotas são diferentes
                    if (i1 != i2)
                    {
                        // std::cout << "start 1 1 improvemente" << std::endl;
                        improved = (Do_1_1_Improvement(instance, solution, i1, i2) || improved);
                        global_improvement = (global_improvement || improved);
                        // std::cout << "end 1 1 improvemente" << std::endl;

                        // std::cout << "start 1 0 improvemente" << std::endl;
                        improved = (Do_1_0_Improvement(instance, solution, i1, i2) || improved);
                        global_improvement = (global_improvement || improved);
                        // std::cout << "end 1 0 improvemente" << std::endl;

                        // std::cout << "star 2 1 improvemente" << std::endl;
                        improved = (Do_2_1_Improvement(instance, solution, i1, i2) || improved);
                        global_improvement = (global_improvement || improved);
                        // std::cout << "end 2 1 improvemente" << std::endl;
                    }
                }
            }
        } while (improved);

        return global_improvement;
    }

    bool ShiftingOneVertex(Instance &instance, HeuristicSolution *sol, int vertex, double &global_variation)
    {
        double profit_variation1 = 0.0, profit_variation2 = 0.0;
        double time_variation1 = 0.0, time_variation2 = std::numeric_limits<double>::infinity();
        double curr_profit_variation = 0.0;
        double curr_time_variation = 0.0;
        std::list<int>::iterator curr_it, it;
        bool can_add = false;
        VertexStatus *status = &((sol->vertex_status_vec_)[vertex]);
        int route = 0, base_route = status->route_;
        int num_routes = instance.num_vehicles();

        if (!(sol->PreviewRemoveVertex(instance, vertex, profit_variation1, time_variation1)))
            return false;

        status->selected_ = false; // suppose vertex is removed!
        for (int curr_route = 0; curr_route < num_routes; ++curr_route)
        {
            if (curr_route != base_route)
            {
                // std::cout << "before preview add vertex minimum" << std::endl;
                // std::cout << *sol << std::endl;
                // sol->BuildBitset(instance);
                // status->selected_ = true;
                // sol->CheckCorrectness(instance);
                // status->selected_ = false;
                if (sol->PreviewAddVertexToRouteWithinMaximumProfitIncrease(instance, vertex, curr_route, curr_it, curr_profit_variation, curr_time_variation))
                {
                    if (double_greater(curr_profit_variation, profit_variation2) || (double_equals(curr_profit_variation, profit_variation2) && double_less(curr_time_variation, time_variation2)))
                    {
                        can_add = true;
                        route = curr_route;
                        it = curr_it;
                        profit_variation2 = curr_profit_variation;
                        time_variation2 = curr_time_variation;
                    }
                }
                // std::cout << "after preview add vertex minimum" << std::endl;
                // sol->BuildBitset(instance);
                // status->selected_ = true;
                // sol->CheckCorrectness(instance);
                // status->selected_ = false;
            }
        }
        status->selected_ = true;

        if (can_add)
        {
            // std::cout << "before inter route move vertex" << std::endl;
            // sol->BuildBitset(instance);
            // sol->CheckCorrectness(instance);
            sol->InterRouteMoveVertex(vertex, route, it, profit_variation1, time_variation1, profit_variation2, time_variation2);
            // std::cout << "after inter route move vertex" << std::endl;
            // sol->BuildBitset(instance);
            // sol->CheckCorrectness(instance);
            global_variation += (time_variation1 + time_variation2);
            return true;
        }

        return false;
    }

    bool TryToInsertUnvisitedVertices(Instance &instance, HeuristicSolution *solution, int type)
    {
        int curr_vertex = 0, curr_route = 0;
        const Graph *graph = instance.graph();
        int num_vertices = graph->num_vertices();
        int num_mandatory = instance.num_mandatory();

        double time_variation = 0.0, profit_variation = 0.0;
        std::list<int>::iterator vertex_it;
        bool added_vertex = false;
        VertexStatus *curr_status = nullptr;

        // Route * route = nullptr;

        if (type == -1)
            type = rand() % 2;

        const auto &ordered_profits = instance.ordered_profits();
        for (size_t i = 0; i < ordered_profits.size(); ++i)
        {
            if (type == 0)
                curr_vertex = (ordered_profits[num_vertices - num_mandatory - 2 - i]).first;
            if (type == 1)
                curr_vertex = (ordered_profits[i]).first;
            curr_status = &((solution->vertex_status_vec_)[curr_vertex]);

            if (!(curr_status->selected_))
            {
                if (solution->PreviewAddVertexWithinMaximumProfitIncrease(instance, curr_vertex, curr_route, vertex_it, profit_variation, time_variation))
                {
                    // std::cout << "before add vertex" << std::endl;
                    // solution->BuildBitset(instance);
                    // solution->CheckCorrectness(instance);
                    solution->AddVertex(curr_vertex, curr_route, vertex_it, profit_variation, time_variation);
                    // std::cout << "after add vertex" << std::endl;
                    // solution->BuildBitset(instance);
                    // solution->CheckCorrectness(instance);
                    added_vertex = true;
                }
                // std::cout << "after preview insert minimum FAILED" << std::endl;
                // solution->BuildBitset(instance);
                // solution->CheckCorrectness(instance);
            }
        }

        return added_vertex;
    }

    bool TryToInsertUnvisitedVerticesOneRoute(Instance &instance, HeuristicSolution *solution, int curr_route)
    {
        int curr_vertex = 0;
        int num_vertices = instance.graph()->num_vertices();
        int num_mandatory = instance.num_mandatory();
        double time_variation = 0.0, profit_variation = 0;
        bool added_vertex = false;
        std::list<int>::iterator vertex_it;
        VertexStatus *curr_status = nullptr;

        int type = rand() % 2;

        const auto &ordered_profits = instance.ordered_profits();
        for (int i = 0; i < (int)(ordered_profits.size()); ++i)
        {
            if (type == 0)
                curr_vertex = (ordered_profits[num_vertices - num_mandatory - 2 - i]).first;
            if (type == 1)
                curr_vertex = (ordered_profits[i]).first;
            curr_status = &((solution->vertex_status_vec_)[curr_vertex]);

            if (!(curr_status->selected_))
            {
                if (solution->PreviewAddVertexToRouteWithinMaximumProfitIncrease(instance, curr_vertex, curr_route, vertex_it, profit_variation, time_variation))
                {
                    // route = &((solution->routes_vec_)[curr_route]);
                    // std::cout << "before add vertex" << std::endl;
                    // solution->BuildBitset(instance);
                    // solution->CheckCorrectness(instance);
                    solution->AddVertex(curr_vertex, curr_route, vertex_it, profit_variation, time_variation);
                    // std::cout << "after add vertex" << std::endl;
                    // solution->BuildBitset(instance);
                    // solution->CheckCorrectness(instance);
                    added_vertex = true;
                }
                // std::cout << "after preview insert minimum FAILED" << std::endl;
                // solution->BuildBitset(instance);
                // solution->CheckCorrectness(instance);
            }
        }

        return added_vertex;
    }

    bool ShiftingAndInsertion(Instance &instance, HeuristicSolution *solution)
    {
        bool shifting = false;
        double global_variation = 0.0;
        int num_vehicles = instance.num_vehicles();
        Route *curr_route = nullptr;
        double initial_profits_sum = solution->profits_sum_;

        for (int i = 0; i < num_vehicles; ++i)
        {
            curr_route = &((solution->routes_vec_)[i]);
            std::list<int> temp = curr_route->vertices_;
            shifting = false;
            for (std::list<int>::iterator it = temp.begin(); it != temp.end(); ++it)
            {
                shifting = ((ShiftingOneVertex(instance, solution, *it, global_variation)) || shifting);
            }

            if (shifting)
            {
                TryToInsertUnvisitedVerticesOneRoute(instance, solution, i);
            }
        }

        // return true if improve total profits or decrease total time without decreasing total profits.
        if (double_greater(solution->profits_sum_, initial_profits_sum) || (double_equals(solution->profits_sum_, initial_profits_sum) && double_less(global_variation, 0.0)))
            return true;
        return false;
    }

    void RemoveVerticesFromSolution(Instance &instance, HeuristicSolution *sol, double percentage)
    {
        Route *curr_route = nullptr;
        VertexStatus *curr_status = nullptr;
        const Graph *graph = instance.graph();
        int num_vertices = graph->num_vertices(), route = 0;
        double time_variation = 0.0, profit_variation = 0.0;
        int num_mandatory = instance.num_mandatory();
        int num_visited_vertices = num_vertices - 2 - (int)((sol->unvisited_vertices_).size());
        int num_visited_profitable_vertices = num_visited_vertices - num_mandatory;
        int max_removed = 0;
        int num_vertices_to_remove = 0;
        int cont = 0, sub_cont = 0, type = 0;
        int num_routes = sol->num_routes_;
        int vertex_position = 0, curr_vertex = 0;
        std::list<int>::iterator vertex_it;

        type = rand() % 3;
        // std::cout << "Remover " << num_removed_poles << " de " << sol->num_visited_poles_ << std::endl;
        if (type == 0)
            max_removed = (int)(std::floor(percentage * num_visited_vertices));
        else
            max_removed = (int)(std::floor(percentage * num_visited_profitable_vertices));

        if (max_removed <= 0)
            return;
        else
            num_vertices_to_remove = 1 + rand() % max_removed;
        // std::cout << "num_visited_vertices:" << num_visited_vertices << std::endl;
        // std::cout << "max removed:" << max_removed << std::endl;
        // std::cout << "num_vertices_to_remove:" << num_vertices_to_remove << std::endl;
        // std::cout << "type: " << type << std::endl;

        const auto &ordered_profits = instance.ordered_profits();
        do
        {
            switch (type)
            {
            // remove randomly
            case 0:
            {
                route = rand() % num_routes;
                curr_route = &((sol->routes_vec_)[route]);
                while ((curr_route->vertices_).empty())
                {
                    route = (route + 1) % num_routes;
                    curr_route = &((sol->routes_vec_)[route]);
                }

                vertex_position = rand() % ((curr_route->vertices_).size());
                vertex_it = (curr_route->vertices_).begin();
                if (vertex_position > 0)
                    std::advance(vertex_it, vertex_position);
                curr_vertex = *vertex_it;

                break;
            }
            // menor prioridade para maior prioridade
            case 1:
            {
                do
                {
                    curr_vertex = (ordered_profits[sub_cont]).first;
                    curr_status = &((sol->vertex_status_vec_)[curr_vertex]);

                    ++sub_cont;
                } while (!(curr_status->selected_));

                break;
            }

            // maior prioridade para menor prioridade
            case 2:
            {
                do
                {
                    curr_vertex = (ordered_profits[num_vertices - num_mandatory - 2 - sub_cont]).first;
                    curr_status = &((sol->vertex_status_vec_)[curr_vertex]);

                    ++sub_cont;
                } while (!(curr_status->selected_));

                break;
            }
            }

            // std::cout << "try to remove " << curr_vertex << std::endl;
            if ((curr_vertex > num_mandatory) && (sol->PreviewRemoveVertex(instance, curr_vertex, profit_variation, time_variation)))
            {
                // std::cout << "before remove vertex" << std::endl;
                // sol->BuildBitset(instance);
                // sol->CheckCorrectness(instance);
                sol->RemoveVertex(curr_vertex, profit_variation, time_variation);
                // std::cout << "after remove vertex" << std::endl;
                // sol->BuildBitset(instance);
                // sol->CheckCorrectness(instance);
            }

            // std::cout << *sol << std::endl;
            // getchar(); getchar();

            ++cont;
        } while (cont < num_vertices_to_remove);
        // std::cout << "Restaram " << sol->num_visited_poles_ << std::endl;
    }
}