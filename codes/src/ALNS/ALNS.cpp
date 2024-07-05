#include <thread>
#include "ALNS.h"
#include "src/general.h"
#include "src/graph_algorithms.h"
#include "src/formulations.h"

ALNS::ALNS()
{
  num_elements_in_pool_ = 0;
  pool_ = std::vector<ALNSHeuristicSolution *>(K_SIZE_OF_ALNS_POOL, nullptr);
  Reset();
}

ALNS::~ALNS()
{
  for (int i = 0; i < num_elements_in_pool_; ++i)
  {
    delete ((pool_)[i]);
    (pool_)[i] = nullptr;
  }
}

void ALNS::Reset()
{
  for (int i = 0; i < num_elements_in_pool_; ++i)
  {
    delete ((pool_)[i]);
    (pool_)[i] = nullptr;
  }

  num_elements_in_pool_ = 0;
  pos_best_sol_ = -1;
  pos_worst_sol_ = -1;
  last_improve_iteration_ = 0;
  non_improve_iterations_counter_ = 0;
}

bool ALNS::CheckPoolIntegrity()
{
  //(mutex_).lock();
  ALNSHeuristicSolution *curr_sol = nullptr;
  double max_profits = -1.0, min_profits = std::numeric_limits<double>::max();
  int max_profits_index = -1, min_profits_index = -1;
  if (num_elements_in_pool_ > K_SIZE_OF_ALNS_POOL)
    throw "Inconsistent pool 1";
  // if(num_elements_in_pool_ != (int)((pool_).size())) throw "Inconsistent pool 2";

  for (int i = 0; i < num_elements_in_pool_; ++i)
  {
    curr_sol = (pool_)[i];
    if (!(curr_sol->CheckCorrectness(*(curr_instance_))))
      throw "Inconsistent pool 3";

    if (double_greater(curr_sol->profits_sum_, max_profits))
    {
      max_profits = curr_sol->profits_sum_;
      max_profits_index = i;
    }

    if (curr_sol->profits_sum_ < min_profits)
    {
      min_profits = curr_sol->profits_sum_;
      min_profits_index = i;
    }
  }

  if (((pool_)[max_profits_index])->profits_sum_ != ((pool_)[pos_best_sol_])->profits_sum_)
    throw "Inconsistent pool 4";
  if (((pool_)[min_profits_index])->profits_sum_ != ((pool_)[pos_worst_sol_])->profits_sum_)
    throw "Inconsistent pool 5";

  //(mutex_).unlock();
  return true;
}

void ALNS::Init(Instance &instance, std::string algo, std::string folder, std::string file_name)
{
  curr_instance_ = &instance;
  const Graph *graph = instance.graph();
  int num_vertices = graph->num_vertices(), num_arcs = graph->num_arcs(), num_mandatory = instance.num_mandatory();
  int num_routes = instance.num_vehicles();

  struct
  {
    bool operator()(const std::pair<int, int> &v1, const std::pair<int, int> &v2) const
    {
      return v1.second < v2.second;
    }
  } order_less;

  ordered_profits_ = std::vector<std::pair<int, int>>(num_vertices - num_mandatory - 1); // -1 stands for the origin vertex. Last element of vector is: num_vertices - num_mandatory - 2!
  const auto vertices_info = graph->vertices_info();
  for (int i = num_mandatory + 1; i < num_vertices; ++i)
    (ordered_profits_)[i - num_mandatory - 1] = std::pair<int, int>(i, vertices_info[i].profit_);
  std::sort((ordered_profits_).begin(), (ordered_profits_).end(), order_less);

  /*for(size_t i = 0; i < (ordered_profits_).size(); ++i)
  {
    std::cout << "(" << ((ordered_profits_)[i]).first << "," << ((ordered_profits_)[i]).second << ")[" << (graph->profits())[((ordered_profits_)[i]).first] << "]" << std::endl;
  }*/

  // std::cout << "read from file " << file_name << std::endl;
  ALNSHeuristicSolution *sol = new ALNSHeuristicSolution(num_vertices, num_arcs, num_routes);
  sol->ReadFromFile(instance, algo, folder, file_name);
  AddSolutionToPool(sol, 0);
  // std::cout << *sol << std::endl;
}

void ALNS::Init(Instance &instance, HeuristicSolution *initial_sol)
{
  curr_instance_ = &instance;
  const Graph *graph = instance.graph();
  int num_vertices = graph->num_vertices(), num_arcs = graph->num_arcs(), num_mandatory = instance.num_mandatory();
  int num_routes = instance.num_vehicles();

  struct
  {
    bool operator()(const std::pair<int, int> &v1, const std::pair<int, int> &v2) const
    {
      return v1.second < v2.second;
    }
  } order_less;

  ordered_profits_ = std::vector<std::pair<int, int>>(num_vertices - num_mandatory - 1); // -1 stands for the origin vertex. Last element of vector is: num_vertices - num_mandatory - 2!
  const auto vertices_info = graph->vertices_info();
  for (int i = num_mandatory + 1; i < num_vertices; ++i)
    (ordered_profits_)[i - num_mandatory - 1] = std::pair<int, int>(i, vertices_info[i].profit_);
  std::sort((ordered_profits_).begin(), (ordered_profits_).end(), order_less);

  /*for(size_t i = 0; i < (ordered_profits_).size(); ++i)
  {
    std::cout << "(" << ((ordered_profits_)[i]).first << "," << ((ordered_profits_)[i]).second << ")[" << (graph->profits())[((ordered_profits_)[i]).first] << "]" << std::endl;
  }*/

  // std::cout << "read from file " << file_name << std::endl;
  ALNSHeuristicSolution *sol = new ALNSHeuristicSolution(initial_sol);
  AddSolutionToPool(sol, 0);
  // std::cout << *sol << std::endl;
}

void ALNS::PrintPool()
{
  (mutex_).lock();
  std::cout << "pool size: " << num_elements_in_pool_ << std::endl;
  std::cout << "pos best: " << pos_best_sol_ << " pos worst: " << pos_worst_sol_ << std::endl;
  for (size_t i = 0; i < num_elements_in_pool_; ++i)
  {
    std::cout << (pool_)[i]->profits_sum_ << " ";
  }
  std::cout << std::endl;
  (mutex_).unlock();
}

ALNSHeuristicSolution *ALNS::best_solution()
{
  if (pos_best_sol_ != -1)
    return (pool_)[pos_best_sol_];
  return nullptr;
}

ALNSHeuristicSolution *ALNS::worst_solution()
{
  if (pos_worst_sol_ != -1)
    return (pool_)[pos_worst_sol_];
  return nullptr;
}

bool ALNS::AddSolutionToPool(ALNSHeuristicSolution *sol, int iter)
{
  ALNSHeuristicSolution *curr_sol = nullptr;
  if (K_SIZE_OF_ALNS_POOL <= 0)
    return false;

  int pos_second_worst_sol = -1;
  bool sol_is_in_pool = false;

  sol->BuildBitset(*(curr_instance_));
  // sol->CheckCorrectness(*curr_instance_);

  // only makes sense if profits do not decay over time!
  // if ((sol->unvisited_vertices_).empty())
  //   sol->is_optimal_ = true;

  if (num_elements_in_pool_ == 0)
  {
    (pool_)[num_elements_in_pool_] = sol;
    (num_elements_in_pool_)++;
    pos_best_sol_ = 0;
    pos_worst_sol_ = 0;
    last_improve_iteration_ = iter;
  }
  else
  {
    // check if sol is already in pool
    sol_is_in_pool = false;
    pos_second_worst_sol = -1;
    if (num_elements_in_pool_ > 1)
      pos_second_worst_sol = (pos_worst_sol_ + 1) % (num_elements_in_pool_);

    for (int i = 0; i < num_elements_in_pool_; ++i)
    {
      curr_sol = (pool_)[i];
      if ((*sol) == (*curr_sol))
      {
        sol_is_in_pool = true;
        break;
      }

      if ((num_elements_in_pool_ > 1) && (i != pos_worst_sol_) && (curr_sol->profits_sum_ < ((pool_)[pos_second_worst_sol])->profits_sum_))
        pos_second_worst_sol = i;
    }

    if (!sol_is_in_pool)
    {
      if (num_elements_in_pool_ < K_SIZE_OF_ALNS_POOL)
      {
        (pool_)[num_elements_in_pool_] = sol;
        (num_elements_in_pool_)++;
        if (double_greater(sol->profits_sum_, ((pool_)[pos_best_sol_])->profits_sum_))
        {
          pos_best_sol_ = num_elements_in_pool_ - 1;
          last_improve_iteration_ = iter;
          // std::cout << "*Melhorou!!!" << std::endl;
        }

        if (sol->profits_sum_ < ((pool_)[pos_worst_sol_])->profits_sum_)
        {
          pos_worst_sol_ = num_elements_in_pool_ - 1;
        }
      }
      else if (num_elements_in_pool_ == K_SIZE_OF_ALNS_POOL)
      {
        // se nova solução é melhor que a pior, substitui a pior solução pela nova solução e atualiza pior solução e melhor solução!
        if (double_greater(sol->profits_sum_, ((pool_)[pos_worst_sol_])->profits_sum_))
        {
          if (double_greater(sol->profits_sum_, ((pool_)[pos_best_sol_])->profits_sum_))
          {
            pos_best_sol_ = pos_worst_sol_;
            last_improve_iteration_ = iter;
            // std::cout << "*Melhorou!!!" << std::endl;
          }

          delete ((pool_)[pos_worst_sol_]);
          (pool_)[pos_worst_sol_] = sol;

          // updates pos of worst solution
          if (double_greater(sol->profits_sum_, ((pool_)[pos_second_worst_sol])->profits_sum_))
            pos_worst_sol_ = pos_second_worst_sol;
        }
        else
        {
          delete sol;
          // CheckPoolIntegrity();
          return false;
        }
      }
    }
    else
    {
      delete sol;
      // CheckPoolIntegrity();
      return false;
    }
  }
  // CheckPoolIntegrity();
  return true;
}

ALNSHeuristicSolution *ALNS::DoPathRelinking(ALNSHeuristicSolution *sol)
{
  ALNSHeuristicSolution *curr_sol = nullptr, *curr_pool_sol = nullptr, *best_sol = sol;
  double curr_similarity = 0.0;

  sol->BuildBitset(*(curr_instance_));

  // creates a copy of the pool of solutions
  (mutex_).lock();
  std::vector<ALNSHeuristicSolution *> pool_copy(num_elements_in_pool_, nullptr);

  for (int i = 0; i < num_elements_in_pool_; ++i)
    pool_copy[i] = new ALNSHeuristicSolution((pool_)[i]);

  // for(int i = 0; i < num_elements_in_pool_; ++i)
  // std::cout << ((pool_)[i])->profits_sum_ << " == " << (pool_copy[i])->profits_sum_ << std::endl;

  // getchar();getchar();

  (mutex_).unlock();

  // for(int i = 0; i < K_NUM_POOL_SOLUTIONS_SELECTED_FOR_PATH_RELINKING+2; ++i)
  for (size_t i = 0; i < pool_copy.size(); ++i)
  {
    // std::cout << " * " << i << std::endl;
    /*(mutex_).lock();
    if(i == 0) curr_pool_sol = new ALNSHeuristicSolution(best_solution());
    else if(i == 1) curr_pool_sol = new ALNSHeuristicSolution(worst_solution());
    else curr_pool_sol = CopyRandomSolutionFromPool();
    (mutex_).unlock();*/

    curr_pool_sol = pool_copy[i];

    // std::cout << curr_pool_sol->bitset_vertices_ << " " << curr_pool_sol->bitset_vertices_.count() << std::endl;
    // std::cout << sol->bitset_vertices_ << " " << sol->bitset_vertices_.count() << std::endl;
    // std::cout << (sol->bitset_vertices_ & curr_pool_sol->bitset_vertices_) << " " << (sol->bitset_vertices_ & curr_pool_sol->bitset_vertices_ ).count() << std::endl;
    curr_similarity = (2.0 * (((curr_pool_sol->bitset_vertices_) & (sol->bitset_vertices_)).count())) / (1.0 * ((curr_pool_sol->bitset_vertices_).count() + (sol->bitset_vertices_).count()));
    // std::cout << curr_similarity << std::endl;
    // getchar(); getchar();
    if (double_less(curr_similarity, K_SIMILARITY_THRESHOLD))
    {
      curr_sol = DoPathRelinkingIter(curr_pool_sol, sol);
      if (double_greater(curr_sol->profits_sum_, best_sol->profits_sum_))
      {
        if (best_sol != sol)
          delete best_sol;
        best_sol = curr_sol;
      }
      else
      {
        delete curr_sol;
        curr_sol = nullptr;
      }

      curr_sol = DoPathRelinkingIter(sol, curr_pool_sol);
      if (double_greater(curr_sol->profits_sum_, best_sol->profits_sum_))
      {
        if (best_sol != sol)
          delete best_sol;
        best_sol = curr_sol;
      }
      else
      {
        delete curr_sol;
        curr_sol = nullptr;
      }
    }

    delete curr_pool_sol;
    curr_pool_sol = nullptr;
  }

  if (best_sol != sol)
    return best_sol;
  return nullptr;
}

ALNSHeuristicSolution *ALNS::DoPathRelinkingIter(ALNSHeuristicSolution *guiding_sol, ALNSHeuristicSolution *new_sol)
{
  const Graph *graph = (curr_instance_)->graph();
  size_t num_routes = (size_t)((curr_instance_)->num_vehicles());
  int num_vertices = graph->num_vertices();
  int num_mandatory = curr_instance_->num_mandatory();
  std::list<int> infeasible_routes;
  int curr_vertex = -1;
  double route_limit = (curr_instance_)->limit();
  ALNSHeuristicSolution *best_sol = new ALNSHeuristicSolution(new_sol), *curr_sol = new ALNSHeuristicSolution(new_sol);

  double time_variation = std::numeric_limits<double>::infinity();
  double curr_time_variation = std::numeric_limits<double>::infinity();
  double profit_variation = 0.0;
  double curr_profit_variation = 0.0;

  int route = -1;
  std::list<int>::iterator it, curr_it;

  std::vector<std::pair<int, int>> ordered_candidate_vertices;

  // std::cout << "before" << std::endl;
  // std::cout << *curr_sol << std::endl;

  // std::cout << "candidates: ";

  auto vertices_info = graph->vertices_info();
  for (int i = 1; i < num_vertices - 1; ++i)
  {
    if (((guiding_sol->bitset_vertices_)[i] == 1) && ((new_sol->bitset_vertices_)[i] == 0))
    {
      ordered_candidate_vertices.push_back(std::pair<int, int>(i, vertices_info[i].profit_));
      // std::cout << i << " ";
    }
  }
  // std::cout << std::endl;
  // getchar(); getchar();

  struct
  {
    bool operator()(const std::pair<int, int> &v1, const std::pair<int, int> &v2) const
    {
      return v1.second < v2.second;
    }
  } order_less;

  std::sort(ordered_candidate_vertices.begin(), ordered_candidate_vertices.end(), order_less);

  VertexStatus *curr_status = nullptr;
  int type = rand() % 2;
  // std::cout << "type: " << type << std::endl;
  int i = 0;
  int size = (int)(ordered_candidate_vertices.size());
  bool can_add = false;

  while (i < size)
  {
    // add vertices (may make routes infeasible)
    while ((infeasible_routes.size() < num_routes) && (i < size))
    {
      if (type == 0)
        curr_vertex = ((ordered_candidate_vertices)[size - 1 - i]).first;
      if (type == 1)
        curr_vertex = ((ordered_candidate_vertices)[i]).first;
      curr_status = &((curr_sol->vertex_status_vec_)[curr_vertex]);

      // std::cout << i << " < " << size << ": " << curr_vertex << " " << (graph->profits())[curr_vertex] << std::endl;
      //  if vertex already added to solution, skip it
      //  this might happen when a candidate vertex is added by the local search in a previous iteration of the outer while loop!
      if (!(curr_status->selected_))
      {
        // try to add to any feasible route with minimum time increase
        can_add = false;
        time_variation = std::numeric_limits<double>::infinity();
        profit_variation = 0;

        for (int curr_route = 0; curr_route < num_routes; ++curr_route)
        {
          // check if route is still feasible
          if (!double_greater(((curr_sol->routes_vec_)[curr_route]).time_, route_limit))
          {
            // std::cout << "try to add " << curr_vertex << " to route " << curr_route << std::endl;
            if (curr_sol->PreviewAddVertexToRouteWithinMaximumProfitIncrease(*(curr_instance_), curr_vertex, curr_route, curr_it, curr_profit_variation, curr_time_variation, true))
            {
              // std::cout << "deu!" << std::endl;
              if (double_less(curr_time_variation, time_variation))
              {
                can_add = true;
                route = curr_route;
                it = curr_it;
                profit_variation = curr_profit_variation;
                time_variation = curr_time_variation;
              }
            } // else std::cout << "nao deu!" << std::endl;
          }
        }

        if (can_add)
        {
          curr_sol->AddVertex(curr_vertex, route, it, profit_variation, time_variation);
          // if route becomes infeasible after add
          if (double_greater(((curr_sol->routes_vec_)[route]).time_, route_limit))
            infeasible_routes.push_back(route);
        }
      }
      ++i;
    }

    // std::cout << "after" << std::endl;
    // std::cout << *curr_sol << std::endl;
    //!(curr_sol->CheckCorrectness(*(curr_instance_)));
    // std::cout << "passou" << std::endl;

    // getchar(); getchar();
    // return new ALNSHeuristicSolution(new_sol);

    // remove vertices to restore feassearchesibility
    Route *curr_route = nullptr;
    for (auto it2 = infeasible_routes.begin(); it2 != infeasible_routes.end(); ++it2)
    {
      curr_route = &((curr_sol->routes_vec_)[*it2]);
      /*std::cout << "before" << std::endl;
      std::cout << *curr_route << std::endl;*/
      std::list<std::pair<int, int>> ordered_vertices_route;

      for (auto it3 = (curr_route->vertices_).begin(); it3 != (curr_route->vertices_).end(); ++it3)
      {
        if (*it3 > num_mandatory)
          ordered_vertices_route.push_back(std::pair<int, int>(*it3, vertices_info[*it3].profit_));
      }

      ordered_vertices_route.sort(order_less);

      // for(auto it3 = ordered_vertices_route.begin(); it3 != ordered_vertices_route.end(); ++it3) std::cout << (*it3).second << " | ";
      // std::cout << std::endl;
      // getchar(); getchar();

      // remove vertices from infeasible route
      std::list<std::pair<int, int>>::iterator it_route_vertices = ordered_vertices_route.begin();
      curr_vertex = (*it_route_vertices).first;
      do
      {
        if (curr_sol->PreviewRemoveVertex(*(curr_instance_), curr_vertex, curr_profit_variation, curr_time_variation, true))
        {
          // std::cout << " * removeu " << curr_vertex << std::endl;
          curr_sol->RemoveVertex(curr_vertex, curr_profit_variation, curr_time_variation);
        }

        ++it_route_vertices;
        curr_vertex = (*it_route_vertices).first;
      } while (double_greater(curr_route->time_, route_limit));

      /*std::cout << "after" << std::endl;
      std::cout << *curr_route << std::endl;
      getchar();
      getchar();*/
    }

    infeasible_routes.clear();
    // curr_sol->CheckCorrectness(*(curr_instance_));

    // Do local searches
    bool continue_search = false;
    do
    {
      continue_search = false;
      do
      {
        continue_search = false;
        if (DoLocalSearchImprovements(curr_sol))
        {
          continue_search = TryToInsertUnvisitedVertices(curr_sol);
        }

        if (DoReplacementImprovements(curr_sol))
        {
          continue_search = true;
          TryToInsertUnvisitedVertices(curr_sol);
        }
      } while (continue_search);

      if (double_greater(curr_sol->profits_sum_, best_sol->profits_sum_))
      {
        delete best_sol;
        best_sol = new ALNSHeuristicSolution(curr_sol);
        // AddSolutionToPool(new ALNSHeuristicSolution(curr_sol),ALNS_iter);
      }
    } while (ShiftingAndInsertion(curr_sol));

    // check is current solution is the best so far at the path relinking
    if (double_greater(curr_sol->profits_sum_, best_sol->profits_sum_))
    {
      delete best_sol;
      best_sol = new ALNSHeuristicSolution(curr_sol);
    }

    // curr_sol->CheckCorrectness(*(curr_instance_));
  }

  // std::cout << guiding_sol->profits_sum_ << " " << new_sol->profits_sum_ << " " << best_sol->profits_sum_ << std::endl;
  // getchar(); getchar();
  delete curr_sol;
  curr_sol = nullptr;
  return best_sol;
}

ALNSHeuristicSolution *ALNS::CopyRandomSolutionFromPool()
{
  if (num_elements_in_pool_ == 0)
    return nullptr;
  return new ALNSHeuristicSolution((pool_)[rand() % (num_elements_in_pool_)]);
}

void ALNS::RunOneThread(int num_thread, int num_iterations)
{
  int iter = 0;
  ALNSHeuristicSolution *curr_sol1 = nullptr, *curr_sol2 = nullptr, *path_relinking_sol = nullptr;
  bool continue_inner_loop = false, converged = false;

  while (!converged)
  {
    // std::cout << num_thread << std::endl;
    // std::cout << iter << std::endl;
    ++iter;
    (mutex_).lock();
    curr_sol1 = CopyRandomSolutionFromPool();
    // std::cout << "after remove from pool" << std::endl;
    // curr_sol1->BuildBitset(*(curr_instance_));
    // curr_sol1->CheckCorrectness(*curr_instance_);
    (mutex_).unlock();
    RemoveVerticesFromSolution(curr_sol1, K_ALNS_PERTURBATION_PERCENTAGE);

    do
    {
      /*DoLocalSearchImprovements(curr_sol1);
      TryToInsertUnvisitedVertices(curr_sol1);
      //std::cout << "1" << std::endl;
      //std::cout << *curr_sol1 << std::endl;
      //std::cout << "- " << curr_sol1->profits_sum_ << std::endl;
      if(DoReplacementImprovements(curr_sol1))
      {
        //std::cout << "trocou" << std::endl;
        //std::cout << *curr_sol1 << std::endl;
        //getchar(); getchar();
        TryToInsertUnvisitedVertices(curr_sol1);
      }*/
      // std::cout << "2" << std::endl;

      bool continue_search = false;
      do
      {
        continue_search = false;
        if (DoLocalSearchImprovements(curr_sol1))
        {
          // std::cout << "after local search" << std::endl;
          // curr_sol1->BuildBitset(*(curr_instance_));
          // curr_sol1->CheckCorrectness(*curr_instance_);
          // std::cout << *curr_sol << std::endl;
          continue_search = true;
          TryToInsertUnvisitedVertices(curr_sol1);
        }

        if (DoReplacementImprovements(curr_sol1))
        {
          // std::cout << "after replacements" << std::endl;
          // curr_sol1->BuildBitset(*(curr_instance_));
          // curr_sol1->CheckCorrectness(*curr_instance_);
          continue_search = true;
          TryToInsertUnvisitedVertices(curr_sol1);
        }
      } while (continue_search);

      (mutex_).lock();
      if (double_greater(curr_sol1->profits_sum_, (best_solution())->profits_sum_))
      {
        // std::cout << *curr_sol1 << std::endl;
        // getchar();
        // getchar();
        AddSolutionToPool(new ALNSHeuristicSolution(curr_sol1), iter);
      }
      (mutex_).unlock();

      // std::cout << *curr_sol1 << std::endl;
    } while (ShiftingAndInsertion(curr_sol1));

    if (K_PATH_RELINKING)
    {
      if (non_improve_iterations_counter_ >= K_NON_IMPROVE_ITERATIONS_LIMIT)
      {

        if (K_STOP_AT_FIRST_PR_DETECTION)
          converged = true;
        else
        {
          // std::cout << " RODOU PR" << std::endl;
          non_improve_iterations_counter_ = 0;
          path_relinking_sol = DoPathRelinking(curr_sol1);

          if (path_relinking_sol != nullptr)
          {
            (mutex_).lock();
            AddSolutionToPool(path_relinking_sol, iter);
            (mutex_).unlock();
          }
        }
      }
    }

    if (iter >= num_iterations)
      converged = true;

    (mutex_).lock();
    AddSolutionToPool(curr_sol1, iter);

    if (last_improve_iteration_ != iter)
      ++non_improve_iterations_counter_;
    else
      non_improve_iterations_counter_ = 0;

    (mutex_).unlock();
  }

  // std::cout << "* " << best_solution()->profits_sum_ << std::endl;
}

void ALNS::Run()
{
  Timestamp *ti = NewTimestamp();
  Timer *timer = GetTimer();
  timer->Clock(ti);
  int iter = 0, initial_solution_profits_sum = 0;
  ALNSHeuristicSolution *curr_sol = nullptr;

  // at this point, must have one solution at pool (generated from feasibility pump)
  curr_sol = best_solution();

  if (!(curr_sol->is_infeasible_) && (curr_sol->is_feasible_))
  {
    // std::cout << *curr_sol << std::endl;
    // std::cout << curr_sol->profits_sum_ << std::endl;
    TryToInsertUnvisitedVertices(curr_sol, 0);
    // std::cout << curr_sol->profits_sum_ << std::endl;
    initial_solution_profits_sum = curr_sol->profits_sum_;
    // std::cout << *curr_sol << std::endl;

    bool continue_search = false;
    do
    {
      continue_search = false;
      // std::cout << *curr_sol << std::endl;
      // std::cout << "tentou local search" << std::endl;
      if (DoLocalSearchImprovements(curr_sol))
      {
        // std::cout << "fez local search" << std::endl;
        // std::cout << *curr_sol << std::endl;
        continue_search = true;
        TryToInsertUnvisitedVertices(curr_sol, 0);
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

      if (DoReplacementImprovements(curr_sol))
      {
        // std::cout << "fez replacements" << std::endl;
        continue_search = true;
        TryToInsertUnvisitedVertices(curr_sol, 0);
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
    curr_sol->BuildBitset(*(curr_instance_));
    // std::cout << curr_sol->profits_sum_ << " ";

    // if(curr_sol->CheckCorrectness(*(curr_instance_)) == false){ std::cout << "deu merda" << std::endl; getchar(); getchar();}

    // std::cout << *curr_sol << std::endl;

    int num_cores = (int)std::thread::hardware_concurrency();
    int iterations_per_thread = (int)std::ceil((1.0 * K_NUM_ITERATIONS_ALNS) / num_cores);

    if ((K_ALNS_MULTI_THREAD) && (num_cores > 1))
    {
      std::vector<std::thread> v(num_cores - 1);

      for (size_t i = 0; i < v.size(); ++i)
      {
        v[i] = std::thread(&ALNS::RunOneThread, this, i, iterations_per_thread);
      }

      RunOneThread(num_cores - 1, iterations_per_thread);

      for (size_t i = 0; i < v.size(); ++i)
      {
        (v[i]).join();
      }
    }
    else
      RunOneThread(0, K_NUM_ITERATIONS_ALNS);

    curr_sol = best_solution();
    /*std::cout << curr_sol->profits_sum_ << std::endl;
    std::cout << "comecou matheuristica" << std::endl;
    curr_sol = BuildBestSolutionFromGraphInducedByPool();
    std::cout << "* " << curr_sol->profits_sum_ << std::endl;*/
    curr_sol->initial_solution_profits_sum_ = initial_solution_profits_sum;
    curr_sol->last_improve_iteration_ = last_improve_iteration_;
    curr_sol->num_iterations_ = iter;
    curr_sol->total_time_spent_ = timer->CurrentElapsedTime(ti);
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

ALNSHeuristicSolution *ALNS::BuildBestSolutionFromPool()
{
  IloEnv env;
  IloModel model(env);
  IloCplex cplex(env);
  cplex.setOut(env.getNullStream());

  const Graph *graph = curr_instance_->graph();
  int num_vehicles = (curr_instance_)->num_vehicles();
  int num_vertices = graph->num_vertices();
  int num_arcs = graph->num_arcs();
  int num_mandatory = curr_instance_->num_mandatory();

  Route *curr_route = nullptr;
  int curr_route_index = -1;
  ALNSHeuristicSolution *new_solution = new ALNSHeuristicSolution(num_vertices, num_arcs, num_vehicles);

  IloNumVarArray h(env, (num_elements_in_pool_)*num_vehicles, 0.0, 1.0, ILOINT);
  std::vector<std::list<int>> vertex_to_route_map(num_vertices);

  IloExpr obj(env);
  IloExpr exp(env);

  for (int i = 0; i < num_elements_in_pool_; ++i)
  {
    // std::cout << i << std::endl;
    for (int j = 0; j < num_vehicles; ++j)
    {
      // std::cout << j << std::endl;
      curr_route = &((((pool_)[i])->routes_vec_)[j]);
      curr_route_index = i * num_vehicles + j;
      // std::cout << *curr_route << std::endl;

      for (auto it = (curr_route->vertices_).begin(); it != (curr_route->vertices_).end(); ++it)
      {
        (vertex_to_route_map[*it]).push_back(curr_route_index);
      }

      char strnum[26];
      sprintf(strnum, "h(%d)(%d)", i, j);
      h[curr_route_index].setName(strnum);
      exp += h[curr_route_index];

      obj += operator*(curr_route->sum_profits_, h[curr_route_index]);
    }
  }

  // std::cout << "* 1" << std::endl;

  model.add(exp <= num_vehicles);
  exp.end();

  model.add(IloMaximize(env, obj));
  obj.end();

  for (int i = 1; i < num_vertices - 1; ++i)
  {
    // std::cout << i << std::endl;
    IloExpr exp2(env);
    for (auto it = (vertex_to_route_map[i]).begin(); it != (vertex_to_route_map[i]).end(); ++it)
      exp2 += h[*it];

    if (i <= num_mandatory)
      model.add(exp2 == 1);
    else
      model.add(exp2 <= 1);
    exp2.end();
  }

  // std::cout << "* 2" << std::endl;
  cplex.extract(model);
  /*cplex.exportModel("Mat_ALNS.lp");
  getchar();getchar();*/

  cplex.solve();

  IloNumArray h_values(env);
  int route_count = 0;

  cplex.getValues(h_values, h);

  for (int i = 0; i < num_elements_in_pool_; ++i)
  {
    for (int j = 0; j < num_vehicles; ++j)
    {
      curr_route = &((((pool_)[i])->routes_vec_)[j]);
      curr_route_index = i * num_vehicles + j;

      if (double_equals(h_values[curr_route_index], 1.0))
      {
        (new_solution->routes_vec_)[route_count] = *curr_route;
        (new_solution->profits_sum_) += (curr_route->sum_profits_);
        ++route_count;
      }
    }
  }

  cplex.end();
  env.end();

  return new_solution;
}

// ALNSHeuristicSolution *ALNS::BuildBestSolutionFromGraphInducedByPool()
// {

//   const Graph *graph = curr_instance_->graph();
//   int num_vehicles = curr_instance_->num_vehicles();
//   int num_vertices = graph->num_vertices();
//   // int num_arcs = graph->num_arcs();
//   int num_mandatory = curr_instance_->num_mandatory();
//   int v1 = 0, v2 = 0;

//   int *profits = new int[num_vertices];

//   auto vertices_info = graph->vertices_info();
//   for (int i = 0; i < num_vertices; ++i)
//     profits[i] = vertices_info[i].profit_;

//   Graph *new_graph = new Graph(num_vertices, num_mandatory, profits);

//   Route *curr_route = nullptr;

//   for (int i = 0; i < num_elements_in_pool_; ++i)
//   {
//     for (int j = 0; j < num_vehicles; ++j)
//     {
//       curr_route = &((((pool_)[i])->routes_vec_)[j]);
//       v1 = v2 = 0;
//       // std::cout << *curr_route << std::endl;

//       for (auto it = (curr_route->vertices_).begin(); it != (curr_route->vertices_).end(); ++it)
//       {
//         v1 = v2;
//         v2 = *it;
//         if ((*new_graph)[v1][v2] == nullptr)
//         {
//           new_graph->AddArc(v1, v2, ((*graph)[v1][v2])->distance());
//         }
//       }
//       v1 = v2;
//       v2 = num_vertices - 1;

//       if ((*new_graph)[v1][v2] == nullptr)
//       {
//         new_graph->AddArc(v1, v2, ((*graph)[v1][v2])->distance());
//       }
//     }
//   }

//   Solution<int> *sol = new Solution<int>(graph->num_vertices());

//   double *R = Dijkstra(new_graph, false, false);
//   double *Rn = Dijkstra(new_graph, true, false);

//   (curr_instance_)->set_graph(new_graph);

//   CapacitatedSingleCommodity(*(curr_instance_), R, Rn, -1, sol, false, false, false, nullptr, nullptr);

//   int num_arcs = new_graph->num_arcs();
//   ALNSHeuristicSolution *new_solution = new ALNSHeuristicSolution(num_vertices, num_arcs, num_vehicles);
//   new_solution->profits_sum_ = sol->lb_;

//   delete sol;
//   sol = nullptr;

//   delete[] R;
//   R = nullptr;

//   delete[] Rn;
//   Rn = nullptr;

//   (curr_instance_)->set_graph(graph);

//   delete new_graph;
//   new_graph = nullptr;

//   return new_solution;
// }

bool ALNS::DoReplacementImprovements(ALNSHeuristicSolution *sol)
{
  bool global_improvement = false;
  bool improved = false;

  do
  {
    improved = (Do_1_1_Unrouted_Improvements(sol) || Do_2_1_Unrouted_Improvements(sol));
    global_improvement = (global_improvement || improved);
  } while (improved);

  return global_improvement;
}

bool ALNS::DoLocalSearchImprovements(ALNSHeuristicSolution *sol)
{
  bool global_improvement = false;
  bool improved = false;

  do
  {
    improved = DoAllInterRoutesImprovements(sol);
    global_improvement = (global_improvement || improved);

    // std::cout << "after inter route" << std::endl;
    // sol->BuildBitset(*(curr_instance_));
    // sol->CheckCorrectness(*curr_instance_);

    improved = DoAllIntraRouteImprovements(sol);
    global_improvement = (global_improvement || improved);
  } while (improved);

  return global_improvement;
}

bool ALNS::DoIntraRouteImprovementOneRoute(ALNSHeuristicSolution *sol, int route)
{
  bool iteration_improvement = false;
  bool improved = false;

  do
  {
    // std::cout << "start 2 opt" << std::endl;
    improved = sol->Do2OptImprovement(*curr_instance_, route);
    iteration_improvement = (iteration_improvement || improved);
    // std::cout << "end 2 opt" << std::endl;
    // sol->BuildBitset(*(curr_instance_));
    // sol->CheckCorrectness(*curr_instance_);
  } while (improved);

  // do
  // {
  //   improved = sol->Do3OptImprovement(graph, route);
  //   iteration_improvement = (iteration_improvement || improved);
  // } while (improved);

  return iteration_improvement;
}

bool ALNS::DoAllIntraRouteImprovements(ALNSHeuristicSolution *solution)
{
  bool global_improvement = false, improved = false;
  int num_vehicles = (curr_instance_)->num_vehicles();

  for (int i = 0; i < num_vehicles; ++i)
  {
    improved = DoIntraRouteImprovementOneRoute(solution, i);
    global_improvement = (global_improvement || improved);
  }

  return global_improvement;
}

bool ALNS::Do_1_1_Unrouted_Improvements(ALNSHeuristicSolution *solution)
{
  Route *curr_route = nullptr;
  std::list<int>::iterator it_i1, it_i2, it_f1, prox_it;
  double profit_variation1 = 0.0;
  double time_variation1 = 0.0;
  int pos_i1 = 0, pos_f1 = 0;
  int max_pos1 = 0;
  int num_routes = (curr_instance_)->num_vehicles();
  int curr_vertex = 0;
  VertexStatus *curr_status = nullptr;
  int num_vertices = ((curr_instance_)->graph())->num_vertices(), num_mandatory = curr_instance_->num_mandatory();
  std::list<int>::iterator next_vertex_to_be_replaced_it, first_vertex_to_be_replaced_it;

  int type = rand() % 2;
  for (size_t i = 0; i < (ordered_profits_).size(); ++i)
  {
    if (type == 0)
      curr_vertex = ((ordered_profits_)[num_vertices - num_mandatory - 2 - i]).first;
    if (type == 1)
      curr_vertex = ((ordered_profits_)[i]).first;
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
            if (solution->PreviewInterRouteSwapUnrouted(*(curr_instance_), r1, pos_i1, pos_f1, it_i1, it_f1, profit_variation1, time_variation1, it_i2))
            {
              if (double_greater(profit_variation1, 0.0) || (double_equals(profit_variation1, 0.0) && double_less(time_variation1, 0.0)))
              {
                // std::cout << "before inter route swap unrouted" << std::endl;
                // solution->BuildBitset(*(curr_instance_));
                // solution->CheckCorrectness(*curr_instance_);
                solution->InterRouteSwapUnrouted(r1, it_i1, it_f1, profit_variation1, time_variation1, it_i2);
                // std::cout << "after inter route swap unrouted" << std::endl;
                // solution->BuildBitset(*(curr_instance_));
                // solution->CheckCorrectness(*curr_instance_);
                return true;
              }
            }
            // std::cout << "after preview inter route unrouted FAILEDA" << std::endl;
            // solution->BuildBitset(*(curr_instance_));
            // solution->CheckCorrectness(*curr_instance_);
          }
          // first_vertex_to_be_replaced_it = next_vertex_to_be_replaced_it;
        }
      }

      // it_i2 = prox_it;
    }
  }

  return false;
}

bool ALNS::Do_2_1_Unrouted_Improvements(ALNSHeuristicSolution *solution)
{
  Route *curr_route = nullptr;
  std::list<int>::iterator it_i1, it_i2, it_f1, prox_it;
  double profit_variation1 = 0.0;
  double time_variation1 = 0.0;
  int pos_i1 = 0, pos_f1 = 0;
  int max_pos1 = 0;
  int num_routes = (curr_instance_)->num_vehicles();
  int curr_vertex = 0;
  VertexStatus *curr_status = nullptr;
  int num_vertices = ((curr_instance_)->graph())->num_vertices(), num_mandatory = curr_instance_->num_mandatory();
  ;
  int type = rand() % 2;
  std::list<int>::iterator next_vertex_to_be_replaced_it, first_vertex_to_be_replaced_it;

  for (size_t i = 0; i < (ordered_profits_).size(); ++i)
  {
    if (type == 0)
      curr_vertex = ((ordered_profits_)[num_vertices - num_mandatory - 2 - i]).first;
    if (type == 1)
      curr_vertex = ((ordered_profits_)[i]).first;
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
            if (solution->PreviewInterRouteSwapUnrouted(*(curr_instance_), r1, pos_i1, pos_f1, it_i1, it_f1, profit_variation1, time_variation1, it_i2))
            {
              if (double_greater(profit_variation1, 0.0) || (double_equals(profit_variation1, 0.0) && double_less(time_variation1, 0.0)))
              {
                // std::cout << "before inter route swap unrouted" << std::endl;
                // solution->BuildBitset(*(curr_instance_));
                // solution->CheckCorrectness(*curr_instance_);
                solution->InterRouteSwapUnrouted(r1, it_i1, it_f1, profit_variation1, time_variation1, it_i2);
                // std::cout << "after inter route swap unrouted" << std::endl;
                // solution->BuildBitset(*(curr_instance_));
                // solution->CheckCorrectness(*curr_instance_);
                return true;
              }
            }
            // std::cout << "after preview inter route unrouted FAILEDA" << std::endl;
            // solution->BuildBitset(*(curr_instance_));
            // solution->CheckCorrectness(*curr_instance_);
          }
        }
      }
      // it_i2 = prox_it;
    }
  }

  return false;
}

bool ALNS::Do_1_1_Improvement(ALNSHeuristicSolution *solution, int r1, int r2)
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
      if (solution->PreviewInterRouteSwap(*(curr_instance_), r1, pos_i1, pos_f1, r2, pos_i2, pos_f2, it_i1, it_f1, profit_variation1, time_variation1, it_i2, it_f2, profit_variation2, time_variation2))
      {
        total_time_variation = time_variation1 + time_variation2;
        total_profit_variation = profit_variation1 + profit_variation2;
        if (double_greater(total_profit_variation, 0.0) || (double_equals(total_profit_variation, 0.0) && double_less(total_time_variation, 0.0)))
        {
          // std::cout << "before inter route swap vertex" << std::endl;
          // solution->BuildBitset(*(curr_instance_));
          // solution->CheckCorrectness(*curr_instance_);
          solution->InterRouteSwap(r1, r2, it_i1, it_f1, profit_variation1, time_variation1, it_i2, it_f2, profit_variation2, time_variation2);
          // std::cout << "after inter route swap vertex" << std::endl;
          // solution->BuildBitset(*(curr_instance_));
          // solution->CheckCorrectness(*curr_instance_);
          return true;
        }
      }
      // std::cout << "after preview inter route FAILED" << std::endl;
      // std::cout << *solution << std::endl;
      // solution->BuildBitset(*(curr_instance_));
      // solution->CheckCorrectness(*curr_instance_);
    }
  }

  return false;
}

bool ALNS::Do_1_0_Improvement(ALNSHeuristicSolution *solution, int r1, int r2)
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
      if (solution->PreviewInterRouteMoveVertex(*(curr_instance_), vertex, r2, pos, it, profit_variation1, time_variation1, profit_variation2, time_variation2))
      {
        total_time_variation = time_variation1 + time_variation2;
        total_profit_variation = profit_variation1 + profit_variation2;
        if (double_greater(total_profit_variation, 0.0) || (double_equals(total_profit_variation, 0.0) && double_less(total_time_variation, 0.0)))
        {
          // std::cout << "before inter route move vertex" << std::endl;
          // solution->BuildBitset(*(curr_instance_));
          // solution->CheckCorrectness(*curr_instance_);
          solution->InterRouteMoveVertex(vertex, r2, it, profit_variation1, time_variation1, profit_variation2, time_variation2);
          // std::cout << "after inter route move vertex" << std::endl;
          // solution->BuildBitset(*(curr_instance_));
          // solution->CheckCorrectness(*curr_instance_);
          return true;
        }
      }
      // std::cout << "after preview inter route FAILED" << std::endl;
      // // std::cout << *solution << std::endl;
      // solution->BuildBitset(*(curr_instance_));
      // solution->CheckCorrectness(*curr_instance_);
    }
  }

  return false;
}

bool ALNS::Do_2_1_Improvement(ALNSHeuristicSolution *solution, int r1, int r2)
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
      if (solution->PreviewInterRouteSwap(*(curr_instance_), r1, pos_i1, pos_f1, r2, pos_i2, pos_f2, it_i1, it_f1, profit_variation1, time_variation1, it_i2, it_f2, profit_variation2, time_variation2))
      {
        total_time_variation = time_variation1 + time_variation2;
        total_profit_variation = profit_variation1 + profit_variation2;
        if (double_greater(total_profit_variation, 0.0) || (double_equals(total_profit_variation, 0.0) && double_less(total_time_variation, 0.0)))
        {
          // std::cout << "before inter route swap vertex" << std::endl;
          // solution->BuildBitset(*(curr_instance_));
          // solution->CheckCorrectness(*curr_instance_);
          solution->InterRouteSwap(r1, r2, it_i1, it_f1, profit_variation1, time_variation1, it_i2, it_f2, profit_variation2, time_variation2);
          // std::cout << "after inter route swap vertex" << std::endl;
          // solution->BuildBitset(*(curr_instance_));
          // solution->CheckCorrectness(*curr_instance_);
          return true;
        }
      }
      // std::cout << "after preview inter route FAILED" << std::endl;
      // solution->BuildBitset(*(curr_instance_));
      // solution->CheckCorrectness(*curr_instance_);
    }
  }

  return false;
}

bool ALNS::DoAllInterRoutesImprovements(ALNSHeuristicSolution *solution)
{
  bool improved = false;
  bool global_improvement = false;
  int num_routes = (curr_instance_)->num_vehicles();

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
          improved = (Do_1_1_Improvement(solution, i1, i2) || improved);
          global_improvement = (global_improvement || improved);
          // std::cout << "end 1 1 improvemente" << std::endl;

          // std::cout << "start 1 0 improvemente" << std::endl;
          improved = (Do_1_0_Improvement(solution, i1, i2) || improved);
          global_improvement = (global_improvement || improved);
          // std::cout << "end 1 0 improvemente" << std::endl;

          // std::cout << "star 2 1 improvemente" << std::endl;
          improved = (Do_2_1_Improvement(solution, i1, i2) || improved);
          global_improvement = (global_improvement || improved);
          // std::cout << "end 2 1 improvemente" << std::endl;
        }
      }
    }
  } while (improved);

  return global_improvement;
}

bool ALNS::ShiftingOneVertex(ALNSHeuristicSolution *sol, int vertex, double &global_variation)
{
  double profit_variation1 = 0.0, profit_variation2 = 0.0;
  double time_variation1 = 0.0, time_variation2 = std::numeric_limits<double>::infinity();
  double curr_profit_variation = 0.0;
  double curr_time_variation = 0.0;
  std::list<int>::iterator curr_it, it;
  bool can_add = false;
  VertexStatus *status = &((sol->vertex_status_vec_)[vertex]);
  int route = 0, base_route = status->route_;
  int num_routes = (curr_instance_)->num_vehicles();

  if (!(sol->PreviewRemoveVertex(*(curr_instance_), vertex, profit_variation1, time_variation1)))
    return false;

  status->selected_ = false; // suppose vertex is removed!
  for (int curr_route = 0; curr_route < num_routes; ++curr_route)
  {
    if (curr_route != base_route)
    {
      // std::cout << "before preview add vertex minimum" << std::endl;
      // std::cout << *sol << std::endl;
      // sol->BuildBitset(*(curr_instance_));
      // sol->CheckCorrectness(*curr_instance_);
      if (sol->PreviewAddVertexToRouteWithinMaximumProfitIncrease(*(curr_instance_), vertex, curr_route, curr_it, curr_profit_variation, curr_time_variation))
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
      // sol->BuildBitset(*(curr_instance_));
      // sol->CheckCorrectness(*curr_instance_);
    }
  }
  status->selected_ = true;

  if (can_add)
  {
    // std::cout << "before inter route move vertex" << std::endl;
    // sol->BuildBitset(*(curr_instance_));
    // sol->CheckCorrectness(*curr_instance_);
    sol->InterRouteMoveVertex(vertex, route, it, profit_variation1, time_variation1, profit_variation2, time_variation2);
    // std::cout << "after inter route move vertex" << std::endl;
    // sol->BuildBitset(*(curr_instance_));
    // sol->CheckCorrectness(*curr_instance_);
    global_variation += (time_variation1 + time_variation2);
    return true;
  }

  return false;
}

bool ALNS::TryToInsertUnvisitedVertices(ALNSHeuristicSolution *solution, int type)
{
  int curr_vertex = 0, curr_route = 0;
  int num_vertices = ((curr_instance_)->graph())->num_vertices();
  int num_mandatory = curr_instance_->num_mandatory();

  double time_variation = 0.0, profit_variation = 0.0;
  std::list<int>::iterator vertex_it;
  const Graph *graph = (curr_instance_)->graph();
  bool added_vertex = false;
  VertexStatus *curr_status = nullptr;

  // Route * route = nullptr;

  if (type == -1)
    type = rand() % 2;

  for (size_t i = 0; i < (ordered_profits_).size(); ++i)
  {
    if (type == 0)
      curr_vertex = ((ordered_profits_)[num_vertices - num_mandatory - 2 - i]).first;
    if (type == 1)
      curr_vertex = ((ordered_profits_)[i]).first;
    curr_status = &((solution->vertex_status_vec_)[curr_vertex]);

    if (!(curr_status->selected_))
    {
      if (solution->PreviewAddVertexWithinMaximumProfitIncrease(*(curr_instance_), curr_vertex, curr_route, vertex_it, profit_variation, time_variation))
      {
        if (double_greater(profit_variation, 0.0) || (double_equals(profit_variation, 0.0) && double_less(time_variation, 0.0)))
        {
          // std::cout << "before add vertex" << std::endl;
          // solution->BuildBitset(*(curr_instance_));
          // solution->CheckCorrectness(*curr_instance_);
          solution->AddVertex(curr_vertex, curr_route, vertex_it, profit_variation, time_variation);
          // std::cout << "after add vertex" << std::endl;
          // solution->BuildBitset(*(curr_instance_));
          // solution->CheckCorrectness(*curr_instance_);
          added_vertex = true;
        }
      }
      // std::cout << "after preview insert minimum FAILED" << std::endl;
      // solution->BuildBitset(*(curr_instance_));
      // solution->CheckCorrectness(*curr_instance_);
    }
  }

  return added_vertex;
}

bool ALNS::TryToInsertUnvisitedVerticesOneRoute(ALNSHeuristicSolution *solution, int curr_route)
{
  int curr_vertex = 0;
  int num_vertices = ((curr_instance_)->graph())->num_vertices();
  int num_mandatory = curr_instance_->num_mandatory();
  double time_variation = 0.0, profit_variation = 0;
  bool added_vertex = false;
  std::list<int>::iterator vertex_it;
  VertexStatus *curr_status = nullptr;

  int type = rand() % 2;

  for (int i = 0; i < (int)((ordered_profits_).size()); ++i)
  {
    if (type == 0)
      curr_vertex = ((ordered_profits_)[num_vertices - num_mandatory - 2 - i]).first;
    if (type == 1)
      curr_vertex = ((ordered_profits_)[i]).first;
    curr_status = &((solution->vertex_status_vec_)[curr_vertex]);

    if (!(curr_status->selected_))
    {
      if (solution->PreviewAddVertexToRouteWithinMaximumProfitIncrease(*(curr_instance_), curr_vertex, curr_route, vertex_it, profit_variation, time_variation))
      {
        // route = &((solution->routes_vec_)[curr_route]);
        // std::cout << "before add vertex" << std::endl;
        // solution->BuildBitset(*(curr_instance_));
        // solution->CheckCorrectness(*curr_instance_);
        solution->AddVertex(curr_vertex, curr_route, vertex_it, profit_variation, time_variation);
        // std::cout << "after add vertex" << std::endl;
        // solution->BuildBitset(*(curr_instance_));
        // solution->CheckCorrectness(*curr_instance_);
        added_vertex = true;
      }
      // std::cout << "after preview insert minimum FAILED" << std::endl;
      // solution->BuildBitset(*(curr_instance_));
      // solution->CheckCorrectness(*curr_instance_);
    }
  }

  return added_vertex;
}

bool ALNS::ShiftingAndInsertion(ALNSHeuristicSolution *solution)
{
  bool shifting = false;
  double global_variation = 0.0;
  int num_vehicles = (curr_instance_)->num_vehicles();
  Route *curr_route = nullptr;
  double initial_profits_sum = solution->profits_sum_;

  for (int i = 0; i < num_vehicles; ++i)
  {
    curr_route = &((solution->routes_vec_)[i]);
    std::list<int> temp = curr_route->vertices_;
    shifting = false;
    for (std::list<int>::iterator it = temp.begin(); it != temp.end(); ++it)
    {
      shifting = ((ShiftingOneVertex(solution, *it, global_variation)) || shifting);
    }

    if (shifting)
    {
      TryToInsertUnvisitedVerticesOneRoute(solution, i);
    }
  }

  // return true if improve total profits or decrease total time without decreasing total profits.
  if (double_greater(solution->profits_sum_, initial_profits_sum) || (double_equals(solution->profits_sum_, initial_profits_sum) && double_less(global_variation, 0.0)))
    return true;
  return false;
}

void ALNS::RemoveVerticesFromSolution(ALNSHeuristicSolution *sol, double percentage)
{
  Route *curr_route = nullptr;
  VertexStatus *curr_status = nullptr;
  const Graph *graph = curr_instance_->graph();
  int num_vertices = graph->num_vertices(), route = 0;
  double time_variation = 0.0, profit_variation = 0.0;
  int num_mandatory = curr_instance_->num_mandatory();
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
        curr_vertex = ((ordered_profits_)[sub_cont]).first;
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
        curr_vertex = ((ordered_profits_)[num_vertices - num_mandatory - 2 - sub_cont]).first;
        curr_status = &((sol->vertex_status_vec_)[curr_vertex]);

        ++sub_cont;
      } while (!(curr_status->selected_));

      break;
    }
    }

    // std::cout << "try to remove " << curr_vertex << std::endl;
    if ((curr_vertex > num_mandatory) && (sol->PreviewRemoveVertex(*(curr_instance_), curr_vertex, profit_variation, time_variation)))
    {
      // std::cout << "before remove vertex" << std::endl;
      // sol->BuildBitset(*(curr_instance_));
      // sol->CheckCorrectness(*curr_instance_);
      sol->RemoveVertex(curr_vertex, profit_variation, time_variation);
      // std::cout << "after remove vertex" << std::endl;
      // sol->BuildBitset(*(curr_instance_));
      // sol->CheckCorrectness(*curr_instance_);
    }

    // std::cout << *sol << std::endl;
    // getchar(); getchar();

    ++cont;
  } while (cont < num_vertices_to_remove);
  // std::cout << "Restaram " << sol->num_visited_poles_ << std::endl;
}
