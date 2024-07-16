#include <thread>
#include "src/ALNS/ALNS.h"
#include "src/general.h"
#include "src/timer.h"
#include "src/graph_algorithms.h"
#include "src/formulations.h"
#include "src/local_searches/local_searches.h"

ALNS::~ALNS()
{
  for (int i = 0; i < num_elements_in_pool_; ++i)
  {
    delete ((pool_)[i]);
    (pool_)[i] = nullptr;
  }
}

// void ALNS::Reset()
// {
//   for (int i = 0; i < num_elements_in_pool_; ++i)
//   {
//     delete ((pool_)[i]);
//     (pool_)[i] = nullptr;
//   }

//   num_elements_in_pool_ = 0;
//   pos_best_sol_ = -1;
//   pos_worst_sol_ = -1;
//   last_improve_iteration_ = 0;
//   non_improve_iterations_counter_ = 0;
// }

bool ALNS::CheckPoolIntegrity()
{
  //(mutex_).lock();
  ALNSHeuristicSolution *curr_sol = nullptr;
  double max_profits = -1.0, min_profits = std::numeric_limits<double>::max();
  int max_profits_index = -1, min_profits_index = -1;
  if (num_elements_in_pool_ > max_pool_size_)
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

ALNS::ALNS(Instance &instance, std::string algo, std::string folder, std::string file_name, int pool_size)
{
  max_pool_size_ = pool_size;
  pool_ = std::vector<ALNSHeuristicSolution *>(pool_size, nullptr);
  curr_instance_ = &instance;
  const Graph *graph = instance.graph();
  int num_vertices = graph->num_vertices(), num_arcs = graph->num_arcs();
  int num_routes = instance.num_vehicles();

  // std::cout << "read from file " << file_name << std::endl;
  ALNSHeuristicSolution *sol = new ALNSHeuristicSolution(num_vertices, num_arcs, num_routes);
  sol->ReadFromFile(instance, algo, folder, file_name);
  AddSolutionToPool(sol);
  time_spent_generating_initial_solution_ = sol->total_time_spent_;
  // std::cout << *sol << std::endl;
}

ALNS::ALNS(Instance &instance, HeuristicSolution *initial_sol, int pool_size)
{
  max_pool_size_ = pool_size;
  pool_ = std::vector<ALNSHeuristicSolution *>(pool_size, nullptr);
  curr_instance_ = &instance;
  // std::cout << "read from file " << file_name << std::endl;
  ALNSHeuristicSolution *sol = new ALNSHeuristicSolution(initial_sol);
  AddSolutionToPool(sol);
  time_spent_generating_initial_solution_ = initial_sol->total_time_spent_;
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

bool ALNS::AddSolutionToPool(ALNSHeuristicSolution *sol)
{
  ALNSHeuristicSolution *curr_sol = nullptr;
  if (max_pool_size_ <= 0)
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
    last_improve_iteration_ = total_iter_;
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
      if (num_elements_in_pool_ < max_pool_size_)
      {
        (pool_)[num_elements_in_pool_] = sol;
        (num_elements_in_pool_)++;
        if (double_greater(sol->profits_sum_, ((pool_)[pos_best_sol_])->profits_sum_))
        {
          pos_best_sol_ = num_elements_in_pool_ - 1;
          last_improve_iteration_ = total_iter_;
          // std::cout << "*Melhorou!!!" << std::endl;
        }

        if (sol->profits_sum_ < ((pool_)[pos_worst_sol_])->profits_sum_)
        {
          pos_worst_sol_ = num_elements_in_pool_ - 1;
        }
      }
      else if (num_elements_in_pool_ == max_pool_size_)
      {
        // se nova solução é melhor que a pior, substitui a pior solução pela nova solução e atualiza pior solução e melhor solução!
        if (double_greater(sol->profits_sum_, ((pool_)[pos_worst_sol_])->profits_sum_))
        {
          if (double_greater(sol->profits_sum_, ((pool_)[pos_best_sol_])->profits_sum_))
          {
            pos_best_sol_ = pos_worst_sol_;
            last_improve_iteration_ = total_iter_;
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
        if (LocalSearches::DoLocalSearchImprovements(*curr_instance_, curr_sol))
        {
          continue_search = LocalSearches::TryToInsertUnvisitedVertices(*curr_instance_, curr_sol);
        }

        if (LocalSearches::DoReplacementImprovements(*curr_instance_, curr_sol))
        {
          continue_search = true;
          LocalSearches::TryToInsertUnvisitedVertices(*curr_instance_, curr_sol);
        }
      } while (continue_search);

      if (double_greater(curr_sol->profits_sum_, best_sol->profits_sum_))
      {
        delete best_sol;
        best_sol = new ALNSHeuristicSolution(curr_sol);
        // AddSolutionToPool(new ALNSHeuristicSolution(curr_sol),ALNS_iter);
      }
    } while (LocalSearches::ShiftingAndInsertion(*curr_instance_, curr_sol));

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
  ALNSHeuristicSolution *curr_sol1 = nullptr, *path_relinking_sol = nullptr;
  bool converged = false;

  while (!converged)
  {
    // std::cout << num_thread << std::endl;
    // std::cout << iter << std::endl;
    ++iter;
    (mutex_).lock();
    ++total_iter_;
    curr_sol1 = CopyRandomSolutionFromPool();
    // std::cout << "after remove from pool" << std::endl;
    // curr_sol1->BuildBitset(*(curr_instance_));
    // curr_sol1->CheckCorrectness(*curr_instance_);
    (mutex_).unlock();
    LocalSearches::RemoveVerticesFromSolution(*curr_instance_, curr_sol1, K_ALNS_PERTURBATION_PERCENTAGE);

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
        if (LocalSearches::DoLocalSearchImprovements(*curr_instance_, curr_sol1))
        {
          // std::cout << "after local search" << std::endl;
          // curr_sol1->BuildBitset(*(curr_instance_));
          // curr_sol1->CheckCorrectness(*curr_instance_);
          // std::cout << *curr_sol << std::endl;
          continue_search = true;
          LocalSearches::TryToInsertUnvisitedVertices(*curr_instance_, curr_sol1);
        }

        if (LocalSearches::DoReplacementImprovements(*curr_instance_, curr_sol1))
        {
          // std::cout << "after replacements" << std::endl;
          // curr_sol1->BuildBitset(*(curr_instance_));
          // curr_sol1->CheckCorrectness(*curr_instance_);
          continue_search = true;
          LocalSearches::TryToInsertUnvisitedVertices(*curr_instance_, curr_sol1);
        }
      } while (continue_search);

      (mutex_).lock();
      // IMPORTANT: ShiftingAndInsertion might lead to a worse solution, so we should check BEFORE if the current solution is better than the best found so far as to not lose track of it!
      if (double_greater(curr_sol1->profits_sum_, (best_solution())->profits_sum_))
      {
        // std::cout << *curr_sol1 << std::endl;
        // getchar();
        // getchar();
        AddSolutionToPool(new ALNSHeuristicSolution(curr_sol1));
      }
      (mutex_).unlock();

      // std::cout << *curr_sol1 << std::endl;
    } while (LocalSearches::ShiftingAndInsertion(*curr_instance_, curr_sol1));

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
            AddSolutionToPool(path_relinking_sol);
            (mutex_).unlock();
          }
        }
      }
    }

    if (iter >= num_iterations)
      converged = true;

    (mutex_).lock();
    AddSolutionToPool(curr_sol1);

    if (last_improve_iteration_ != total_iter_)
      ++non_improve_iterations_counter_;
    else
      non_improve_iterations_counter_ = 0;

    (mutex_).unlock();
  }

  // std::cout << "* " << best_solution()->profits_sum_ << std::endl;
}

void ALNS::Run(int num_iterations, bool multithreading)
{
  Timestamp *ti = NewTimestamp();
  Timer *timer = GetTimer();
  timer->Clock(ti);
  int initial_solution_profits_sum = 0;
  total_iter_ = 0;
  ALNSHeuristicSolution *curr_sol = nullptr;

  // at this point, must have one solution at pool (generated from feasibility pump)
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
    curr_sol->BuildBitset(*(curr_instance_));
    // std::cout << curr_sol->profits_sum_ << " ";

    // if(curr_sol->CheckCorrectness(*(curr_instance_)) == false){ std::cout << "deu merda" << std::endl; getchar(); getchar();}

    // std::cout << *curr_sol << std::endl;

    int num_cores = (int)std::thread::hardware_concurrency();
    // std::cout << "cores: " << num_cores << std::endl;
    int iterations_per_thread = (int)std::floor((1.0 * num_iterations) / num_cores);
    int rest_of_iterations = num_iterations % num_cores;
    // std::cout << "iterations_per_thread: " << iterations_per_thread << std::endl;
    // std::cout << "rest_of_iterations: " << rest_of_iterations << std::endl;

    if ((multithreading) && (num_cores > 1))
    {
      std::vector<std::thread> v(num_cores - 1);

      for (size_t i = 0; i < v.size(); ++i)
      {
        v[i] = std::thread(&ALNS::RunOneThread, this, i, iterations_per_thread);
      }

      RunOneThread(num_cores - 1, iterations_per_thread + rest_of_iterations);

      for (size_t i = 0; i < v.size(); ++i)
      {
        (v[i]).join();
      }
    }
    else
      RunOneThread(0, num_iterations);

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