#include <iostream>
#include <algorithm>
#include <queue>
#include <dirent.h>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <unistd.h>
#include <getopt.h>
#include <cstdlib>
#include <boost/dynamic_bitset.hpp>
#include "src/graph.h"
#include "src/instance.h"
#include "src/timer.h"
#include "src/matrix.hpp"
#include "src/general.h"
#include "src/graph_algorithms.h"
#include "src/solution.hpp"
#include "src/feasibility_pump/feasibility_pump.h"
#include "src/local_branching.h"
#include "src/ALNS/ALNS.h"

void GenerateALNSLatexTable(std::vector<std::string> dirs, bool stop)
{
  // std::cout << "3" << std::endl; getchar(); getchar();
  std::string curr_file, algo_best_primal_bound = "cb_csc3";
  std::vector<std::string> algorithms;

  // algorithms.push_back("csc3_initial_primal_bound");
  // algorithms.push_back("cb_csc3");
  // algorithms.push_back("relax_cb_angle_0.03_clique_ccs_lcis_csc3");
  // algorithms.push_back("fp_alpha2_0_p_fixed_10");
  // algorithms.push_back("fp_alpha2_100_p_fixed_10");

  // algorithms.push_back("fp_alpha2_0_p_fixed_10_ALNS_1000");
  // algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_ALNS_1000");
  // algorithms.push_back("fp_alpha2_100_p_fixed_10_ALNS_1000");
  // algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_ALNS_1000");

  // algorithms.push_back("fp_alpha2_0_p_fixed_10_ALNS_2000");
  // algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_ALNS_2000");
  // algorithms.push_back("fp_alpha2_100_p_fixed_10_ALNS_2000");
  // algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_ALNS_2000");

  // algorithms.push_back("fp_alpha2_0_p_fixed_10_PR_ALNS_5000");
  // algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_5000");
  // algorithms.push_back("fp_alpha2_100_p_fixed_10_PR_ALNS_5000");
  // algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_5000");

  // algorithms.push_back("fp_alpha2_0_p_fixed_10_PR_ALNS_5000");
  // algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_5000");
  // algorithms.push_back("fp_alpha2_100_p_fixed_10_PR_ALNS_5000");
  // algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_5000");

  // algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_100000");
  // algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_100000");

  // algorithms.push_back("fp_alpha2_0_p_fixed_10_ALNS_CONV");
  // algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_ALNS_CONV");
  algorithms.push_back("fp_alpha2_100_p_fixed_10_ALNS_CONV");
  algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_ALNS_CONV");

  std::fstream output;
  std::string output_name;
  if (stop)
    output_name = ".//tables//latex//stop_table_ALNS.txt";
  else
    output_name = ".//tables//latex//table_ALNS.txt";
  output.open(output_name.c_str(), std::fstream::out);

  if (!output.is_open())
  {
    std::cout << "Could not open file " << output_name << std::endl;
    throw 1;
  }

  // std::cout << output_name << std::endl;

  output << std::setprecision(2) << std::fixed;

  std::vector<std::vector<double>> total_primal_gap_per_algo(algorithms.size(), std::vector<double>());
  std::vector<double> total_avg_primal_gap_per_algo(algorithms.size(), 0.0);

  std::vector<std::vector<double>> total_best_primal_gap_per_algo(algorithms.size(), std::vector<double>());
  std::vector<double> total_avg_best_primal_gap_per_algo(algorithms.size(), 0.0);

  std::vector<std::vector<double>> total_profits_per_algo(algorithms.size(), std::vector<double>());
  std::vector<double> total_avg_profits_per_algo(algorithms.size(), 0.0);

  std::vector<std::vector<double>> total_best_profits_per_algo(algorithms.size(), std::vector<double>());
  std::vector<double> total_avg_best_profits_per_algo(algorithms.size(), 0.0);

  std::vector<std::vector<double>> total_max_improve_iter_per_algo(algorithms.size(), std::vector<double>());
  std::vector<double> total_avg_max_improve_iter_per_algo(algorithms.size(), 0.0);

  std::vector<std::vector<double>> total_time_per_algo1(algorithms.size(), std::vector<double>());
  std::vector<double> total_avg_time_per_algo1(algorithms.size(), 0.0);

  for (size_t k = 0; k < dirs.size(); ++k)
  {
    std::vector<std::vector<double>> primal_gap_per_algo(algorithms.size(), std::vector<double>());
    std::vector<double> avg_primal_gap_per_algo(algorithms.size(), 0.0);
    std::vector<double> st_dev_primal_gap_per_algo(algorithms.size(), 0.0);

    std::vector<std::vector<double>> best_primal_gap_per_algo(algorithms.size(), std::vector<double>());
    std::vector<double> avg_best_primal_gap_per_algo(algorithms.size(), 0.0);
    std::vector<double> st_dev_best_primal_gap_per_algo(algorithms.size(), 0.0);

    std::vector<std::vector<double>> time_per_algo1(algorithms.size(), std::vector<double>());
    std::vector<double> avg_time_per_algo1(algorithms.size(), 0.0);

    std::vector<std::vector<double>> profits_per_algo(algorithms.size(), std::vector<double>());
    std::vector<double> avg_profits_per_algo(algorithms.size(), 0.0);

    std::vector<std::vector<double>> max_improve_iter_per_algo(algorithms.size(), std::vector<double>());
    std::vector<double> avg_max_improve_iter_per_algo(algorithms.size(), 0.0);

    std::vector<std::vector<double>> best_profits_per_algo(algorithms.size(), std::vector<double>());
    std::vector<double> avg_best_profits_per_algo(algorithms.size(), 0.0);

    std::vector<std::string> instances;
    std::string folder = dirs[k].substr(20);
    AddInstancesFromDirectory(dirs[k], instances, false);
    for (size_t i = 0; i < instances.size(); ++i)
    {
      bool feasible = true;
      double heuristic_lb = 0.0, best_lb = 0.0, lb = 0.0, best_primal_gap = 0.0, curr_primal_gap = 0.0, best_heuristic_lb = 0.0, curr_max_improve_iter = 0;

      std::fstream input;
      std::string status;
      std::string line;

      // retrieve best primal bound (2h of execution)
      curr_file = ".//solutions//";
      curr_file.append(folder);
      curr_file.append("s_");
      if (stop)
        curr_file.append("stop_");
      curr_file.append(algo_best_primal_bound);
      curr_file.append("_");
      curr_file.append(instances[i]);
      // std::cout << curr_file << std::endl;

      // std::fstream input;
      input.open(curr_file.c_str(), std::fstream::in);

      if (!input.is_open())
      {
        std::cout << "Could not open file " << curr_file << std::endl;
        throw 4;
      }

      std::stringstream s_lb;
      lb = 0.0;

      getline(input, line);
      size_t pos = line.find_first_of(":");
      status = line.substr(pos + 2);

      std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
      status.erase(end_pos, status.end());

      // std::cout << status << std::endl;
      getline(input, line);
      getline(input, line);
      pos = line.find_first_of(":");
      s_lb << line.substr(pos + 2);
      if (s_lb.str() == "-inf")
        lb = -1.0;
      else
        s_lb >> lb;
      // std::cout << s_lb.str() << std::endl;
      // std::cout << lb << std::endl;

      if (status == "INFEASIBLE")
        feasible = false;

      if ((double_equals(lb, -1)) || (status == "INFEASIBLE"))
      {
        best_lb = -1.0;
      }
      else
        best_lb = lb;

      // std::cout << "best_lb: " << best_lb << std::endl;
      input.close();

      for (size_t j = 0; j < algorithms.size(); j++)
      {
        int count = 0;
        double avg_lb = 0.0, avg_t1 = 0.0, avg_profits = 0.0, avg_max_improve_iter = 0.0;

        best_heuristic_lb = 0.0;

        for (int seed = 1; seed <= 10; ++seed)
        {
          lb = 0.0;
          double t1 = 0.0;
          // retrieve primal bound at root node
          curr_file = ".//solutions//";
          curr_file.append(folder);
          curr_file.append("s_");
          if (stop)
            curr_file.append("stop_");
          curr_file.append(algorithms[j]);
          curr_file.append("_seed_");
          curr_file.append(std::to_string(seed));
          curr_file.append("_");
          curr_file.append(instances[i]);
          // std::cout << curr_file << std::endl;

          // std::fstream input;
          input.open(curr_file.c_str(), std::fstream::in);

          if (!input.is_open())
          {
            std::cout << "Could not open file " << curr_file << std::endl;
            throw 4;
          }

          std::stringstream s_lb1, s_t1, s_max_improve_iter;

          getline(input, line);
          size_t pos = line.find_first_of(":");
          status = line.substr(pos + 2);

          std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
          status.erase(end_pos, status.end());

          // std::cout << status << std::endl;
          getline(input, line);
          // getline(input,line);
          pos = line.find_first_of(":");
          s_lb1 << line.substr(pos + 2);
          s_lb1 >> lb;

          if (!(double_equals(best_lb, -1)) && !(double_equals(lb, -1)) && (status != "INFEASIBLE") && (status != "POSSIBLYFEASIBLE") && (status != "FAILEDTOFINDAFEASIBLESOLUTION"))
          {
            avg_lb += lb;
            avg_profits += lb;

            if (double_greater(lb, best_heuristic_lb))
              best_heuristic_lb = lb;

            /*if(double_greater(lb,best_lb))
            {
              std::cout << instances[i] << " seed " << seed << std::endl;
              std::cout << lb << " > " << best_lb << std::endl;
            }*/
            count++;
          }

          getline(input, line);
          getline(input, line);
          getline(input, line);

          pos = line.find_first_of(":");
          s_t1 << line.substr(pos + 2);
          s_t1 >> t1;
          avg_t1 += t1;

          getline(input, line);
          pos = line.find_first_of(":");
          s_max_improve_iter << line.substr(pos + 2);
          s_max_improve_iter >> curr_max_improve_iter;
          // std::cout << "root_lb: " << root_lb << std::endl;

          avg_max_improve_iter += curr_max_improve_iter;
          input.close();
        }

        if (count > 0)
          heuristic_lb = avg_lb / (1.0 * count);
        else
          heuristic_lb = -1.0;

        avg_t1 /= 10.0;
        avg_max_improve_iter /= 10.0;
        // if(count > 0) avg_profits /= (1.0*count);
        avg_profits /= 10.0;

        if (double_greater(avg_profits, best_heuristic_lb))
          std::cout << avg_profits << " > " << best_heuristic_lb << std::endl;

        if (!feasible)
        {
          curr_primal_gap = 0.0;
          best_primal_gap = 0.0;
        }
        else if ((double_equals(best_lb, -1)) || (double_equals(heuristic_lb, -1)))
        {
          curr_primal_gap = 100.0;
          best_primal_gap = 100.0;
        }
        else
        {
          if (double_equals(best_lb, heuristic_lb))
            curr_primal_gap = 0.0;
          else
            curr_primal_gap = (100 * (best_lb - heuristic_lb)) / best_lb;
          if (double_equals(best_lb, best_heuristic_lb))
            best_primal_gap = 0.0;
          else
            best_primal_gap = (100 * (best_lb - best_heuristic_lb)) / best_lb;
        }

        // std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

        primal_gap_per_algo[j].push_back(curr_primal_gap);
        total_primal_gap_per_algo[j].push_back(curr_primal_gap);
        avg_primal_gap_per_algo[j] += curr_primal_gap;
        total_avg_primal_gap_per_algo[j] += curr_primal_gap;

        best_primal_gap_per_algo[j].push_back(best_primal_gap);
        total_best_primal_gap_per_algo[j].push_back(best_primal_gap);
        avg_best_primal_gap_per_algo[j] += best_primal_gap;
        total_avg_best_primal_gap_per_algo[j] += best_primal_gap;

        profits_per_algo[j].push_back(avg_profits);
        total_profits_per_algo[j].push_back(avg_profits);
        avg_profits_per_algo[j] += avg_profits;
        total_avg_profits_per_algo[j] += avg_profits;

        max_improve_iter_per_algo[j].push_back(avg_max_improve_iter);
        total_max_improve_iter_per_algo[j].push_back(avg_max_improve_iter);
        avg_max_improve_iter_per_algo[j] += avg_max_improve_iter;
        total_avg_max_improve_iter_per_algo[j] += avg_max_improve_iter;

        best_profits_per_algo[j].push_back(best_heuristic_lb);
        total_best_profits_per_algo[j].push_back(best_heuristic_lb);
        avg_best_profits_per_algo[j] += best_heuristic_lb;
        total_avg_best_profits_per_algo[j] += best_heuristic_lb;

        time_per_algo1[j].push_back(avg_t1);
        total_time_per_algo1[j].push_back(avg_t1);
        avg_time_per_algo1[j] += avg_t1;
        total_avg_time_per_algo1[j] += avg_t1;
      }

      // getchar();getchar();
    }

    output << k + 1;

    for (size_t h = 0; h < algorithms.size(); h++)
    {
      avg_primal_gap_per_algo[h] /= (1.0 * ((primal_gap_per_algo[h]).size()));
      st_dev_primal_gap_per_algo[h] = StDev(primal_gap_per_algo[h], avg_primal_gap_per_algo[h]);

      avg_best_primal_gap_per_algo[h] /= (1.0 * ((best_primal_gap_per_algo[h]).size()));
      st_dev_best_primal_gap_per_algo[h] = StDev(best_primal_gap_per_algo[h], avg_best_primal_gap_per_algo[h]);

      avg_profits_per_algo[h] /= (1.0 * ((profits_per_algo[h]).size()));

      avg_max_improve_iter_per_algo[h] /= (1.0 * ((max_improve_iter_per_algo[h]).size()));

      avg_best_profits_per_algo[h] /= (1.0 * ((best_profits_per_algo[h]).size()));

      avg_time_per_algo1[h] /= (1.0 * ((time_per_algo1[h]).size()));

      output << " && " << avg_primal_gap_per_algo[h] << " & " << st_dev_primal_gap_per_algo[h] << " & " << avg_profits_per_algo[h] << " & " << avg_best_profits_per_algo[h] << " & " << avg_time_per_algo1[h];
    }
    output << "\\\\" << std::endl;
  }

  output << "Total";
  for (size_t j = 0; j < algorithms.size(); j++)
  {
    // std::cout << total_improvement_per_algo[j].size() << std::endl;
    total_avg_primal_gap_per_algo[j] /= (1.0 * ((total_primal_gap_per_algo[j]).size()));
    total_avg_best_primal_gap_per_algo[j] /= (1.0 * ((total_best_primal_gap_per_algo[j]).size()));
    total_avg_profits_per_algo[j] /= (1.0 * ((total_profits_per_algo[j]).size()));
    total_avg_max_improve_iter_per_algo[j] /= (1.0 * ((total_max_improve_iter_per_algo[j]).size()));
    total_avg_best_profits_per_algo[j] /= (1.0 * ((total_best_profits_per_algo[j]).size()));
    total_avg_time_per_algo1[j] /= (1.0 * ((total_time_per_algo1[j]).size()));

    output << " && " << total_avg_primal_gap_per_algo[j] << " & " << StDev(total_primal_gap_per_algo[j], total_avg_primal_gap_per_algo[j]) << " & " << total_avg_profits_per_algo[j] << " & " << total_avg_best_profits_per_algo[j] << " & " << total_avg_time_per_algo1[j];
  }

  output << "\\\\" << std::endl;

  output.close();
}

void GenerateALNSLatexTable2(std::vector<std::string> dirs, bool stop)
{
  // std::cout << "3" << std::endl; getchar(); getchar();
  std::string curr_file, algo_best_primal_bound = "cb_csc3";
  std::vector<std::string> algorithms;

  // algorithms.push_back("csc3_initial_primal_bound");
  // algorithms.push_back("cb_csc3");
  // algorithms.push_back("relax_cb_angle_0.03_clique_ccs_lcis_csc3");
  // algorithms.push_back("fp_alpha2_0_p_fixed_10");
  // algorithms.push_back("fp_alpha2_100_p_fixed_10");

  // algorithms.push_back("fp_alpha2_0_p_fixed_10_ALNS_1000");
  // algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_ALNS_1000");
  // algorithms.push_back("fp_alpha2_100_p_fixed_10_ALNS_1000");
  // algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_ALNS_1000");

  // algorithms.push_back("fp_alpha2_0_p_fixed_10_ALNS_2000");
  // algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_ALNS_2000");
  // algorithms.push_back("fp_alpha2_100_p_fixed_10_ALNS_2000");
  // algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_ALNS_2000");

  // algorithms.push_back("fp_alpha2_0_p_fixed_10_PR_ALNS_5000");
  // algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_5000");
  // algorithms.push_back("fp_alpha2_100_p_fixed_10_PR_ALNS_5000");
  // algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_5000");

  // algorithms.push_back("fp_alpha2_0_p_fixed_10_PR_ALNS_5000");
  algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_5000");
  // algorithms.push_back("fp_alpha2_100_p_fixed_10_PR_ALNS_5000");
  algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_5000");

  std::fstream output;
  std::string output_name;
  if (stop)
    output_name = ".//tables//latex//stop_table_ALNS2.txt";
  else
    output_name = ".//tables//latex//table_ALNS2.txt";
  output.open(output_name.c_str(), std::fstream::out);

  if (!output.is_open())
  {
    std::cout << "Could not open file " << output_name << std::endl;
    throw 1;
  }

  // std::cout << output_name << std::endl;

  output << std::setprecision(2) << std::fixed;

  std::vector<int> total_num_reached_best_per_algo(algorithms.size(), 0);
  std::vector<int> total_num_improved_best_per_algo(algorithms.size(), 0);

  size_t total_num_instances = 0;
  for (size_t k = 0; k < dirs.size(); ++k)
  {
    std::vector<int> num_reached_best_per_algo(algorithms.size(), 0);
    std::vector<int> num_improved_best_per_algo(algorithms.size(), 0);

    std::vector<std::string> instances;
    std::string folder = dirs[k].substr(20);
    AddInstancesFromDirectory(dirs[k], instances, false);
    size_t num_instances = instances.size();
    total_num_instances += num_instances;
    for (size_t i = 0; i < num_instances; ++i)
    {
      bool feasible = true;
      double heuristic_lb = 0.0, best_lb = 0.0, lb = 0.0, best_heuristic_lb = 0.0;

      std::fstream input;
      std::string status;
      std::string line;

      // retrieve best primal bound (2h of execution)
      curr_file = ".//solutions//";
      curr_file.append(folder);
      curr_file.append("s_");
      if (stop)
        curr_file.append("stop_");
      curr_file.append(algo_best_primal_bound);
      curr_file.append("_");
      curr_file.append(instances[i]);
      // std::cout << curr_file << std::endl;

      // std::fstream input;
      input.open(curr_file.c_str(), std::fstream::in);

      if (!input.is_open())
      {
        std::cout << "Could not open file " << curr_file << std::endl;
        throw 4;
      }

      std::stringstream s_lb;
      lb = 0.0;

      getline(input, line);
      size_t pos = line.find_first_of(":");
      status = line.substr(pos + 2);

      std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
      status.erase(end_pos, status.end());

      // std::cout << status << std::endl;
      getline(input, line);
      getline(input, line);
      pos = line.find_first_of(":");
      s_lb << line.substr(pos + 2);
      if (s_lb.str() == "-inf")
        lb = -1.0;
      else
        s_lb >> lb;
      // std::cout << s_lb.str() << std::endl;
      // std::cout << lb << std::endl;

      if (status == "INFEASIBLE")
        feasible = false;

      if ((double_equals(lb, -1)) || (status == "INFEASIBLE"))
      {
        best_lb = -1.0;
      }
      else
        best_lb = lb;

      // std::cout << "best_lb: " << best_lb << std::endl;
      input.close();

      for (size_t j = 0; j < algorithms.size(); j++)
      {
        int count = 0;
        double avg_lb = 0.0;

        best_heuristic_lb = 0.0;

        for (int seed = 1; seed <= 10; ++seed)
        {
          lb = 0.0;
          double t1 = 0.0;
          // retrieve primal bound at root node
          curr_file = ".//solutions//";
          curr_file.append(folder);
          curr_file.append("s_");
          if (stop)
            curr_file.append("stop_");
          curr_file.append(algorithms[j]);
          curr_file.append("_seed_");
          curr_file.append(std::to_string(seed));
          curr_file.append("_");
          curr_file.append(instances[i]);
          // std::cout << curr_file << std::endl;

          // std::fstream input;
          input.open(curr_file.c_str(), std::fstream::in);

          if (!input.is_open())
          {
            std::cout << "Could not open file " << curr_file << std::endl;
            throw 4;
          }

          std::stringstream s_lb1, s_t1, s_max_improve_iter;

          getline(input, line);
          size_t pos = line.find_first_of(":");
          status = line.substr(pos + 2);

          std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
          status.erase(end_pos, status.end());

          // std::cout << status << std::endl;
          getline(input, line);
          // getline(input,line);
          pos = line.find_first_of(":");
          s_lb1 << line.substr(pos + 2);
          s_lb1 >> lb;

          if (!(double_equals(best_lb, -1)) && !(double_equals(lb, -1)) && (status != "INFEASIBLE") && (status != "POSSIBLYFEASIBLE") && (status != "FAILEDTOFINDAFEASIBLESOLUTION"))
          {
            avg_lb += lb;

            if (double_greater(lb, best_heuristic_lb))
              best_heuristic_lb = lb;

            /*if(double_greater(lb,best_lb))
            {
              std::cout << instances[i] << " seed " << seed << std::endl;
              std::cout << lb << " > " << best_lb << std::endl;
            }*/
            count++;
          }

          input.close();
        }

        if (count > 0)
          heuristic_lb = avg_lb / (1.0 * count);
        else
          heuristic_lb = -1.0;

        if (!feasible)
        {
          ++(num_reached_best_per_algo[j]);
          ++(total_num_reached_best_per_algo[j]);
        }
        else if ((!double_equals(best_lb, -1)) && (!double_equals(heuristic_lb, -1)))
        {
          if (!double_less(best_heuristic_lb, best_lb))
          {
            ++(num_reached_best_per_algo[j]);
            ++(total_num_reached_best_per_algo[j]);
          }

          if (double_greater(best_heuristic_lb, best_lb))
          {
            ++(num_improved_best_per_algo[j]);
            ++(total_num_improved_best_per_algo[j]);
          }
        }
      }

      // getchar();getchar();
    }

    output << k + 1;

    for (size_t h = 0; h < algorithms.size(); h++)
    {
      output << "& & " << num_reached_best_per_algo[h] << "/" << num_instances << " & " << num_improved_best_per_algo[h] << "/" << num_instances;
    }
    output << "\\\\" << std::endl;
  }

  output << "Total";
  for (size_t j = 0; j < algorithms.size(); j++)
  {
    output << " & &" << total_num_reached_best_per_algo[j] << "/" << total_num_instances << " & " << total_num_improved_best_per_algo[j] << "/" << total_num_instances;
  }

  output << "\\\\" << std::endl;

  output.close();
}

void GenerateALNSLatexAppendix(std::vector<std::string> dirs)
{
  bool stop = false;
  // std::cout << "3" << std::endl; getchar(); getchar();
  std::vector<std::string> algorithms;

  // algorithms.push_back("csc3_initial_primal_bound");
  // algorithms.push_back("cb_csc3");
  // algorithms.push_back("relax_cb_angle_0.03_clique_ccs_lcis_csc3");
  // algorithms.push_back("fp_alpha2_0_p_fixed_10");
  // algorithms.push_back("fp_alpha2_100_p_fixed_10");

  // algorithms.push_back("fp_alpha2_0_p_fixed_10_ALNS_1000");
  // algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_ALNS_1000");
  // algorithms.push_back("fp_alpha2_100_p_fixed_10_ALNS_1000");
  // algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_ALNS_1000");

  // algorithms.push_back("fp_alpha2_0_p_fixed_10_ALNS_2000");
  // algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_ALNS_2000");
  // algorithms.push_back("fp_alpha2_100_p_fixed_10_ALNS_2000");
  // algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_ALNS_2000");

  // algorithms.push_back("fp_alpha2_0_p_fixed_10_PR_ALNS_5000");
  // algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_5000");
  // algorithms.push_back("fp_alpha2_100_p_fixed_10_PR_ALNS_5000");
  // algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_5000");

  // algorithms.push_back("fp_alpha2_0_p_fixed_10_PR_ALNS_5000");
  algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_5000");
  // algorithms.push_back("fp_alpha2_100_p_fixed_10_PR_ALNS_5000");
  algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_5000");

  std::fstream output;
  std::string output_name;
  output_name = ".//tables//latex//table_appendix_ALNS.txt";
  output.open(output_name.c_str(), std::fstream::out);

  if (!output.is_open())
  {
    std::cout << "Could not open file " << output_name << std::endl;
    throw 1;
  }

  // std::cout << output_name << std::endl;

  output << std::setprecision(2) << std::fixed;

  for (size_t k = 0; k < dirs.size(); ++k)
  {
    std::vector<std::string> instances;
    std::string folder = dirs[k].substr(20);
    AddInstancesFromDirectory(dirs[k], instances, false);

    std::sort(instances.begin(), instances.end());
    size_t num_instances = instances.size();

    for (size_t i = 0; i < num_instances; ++i)
    {
      output << instances[i].substr(0, instances[i].size() - 4) << "(\\_5\\%)";

      double heuristic_lb = 0.0, lb = 0.0, best_heuristic_lb = 0.0;

      std::fstream input;
      std::string status;
      std::string line;

      for (size_t l = 0; l <= 1; ++l)
      {
        if (l == 0)
          stop = false;
        else
          stop = true;

        for (size_t j = 0; j < algorithms.size(); j++)
        {
          int count = 0;
          double avg_lb = 0.0;

          best_heuristic_lb = 0.0;

          for (int seed = 1; seed <= 10; ++seed)
          {
            lb = 0.0;
            double t1 = 0.0;
            // retrieve primal bound at root node
            std::string curr_file = ".//solutions//";
            curr_file.append(folder);
            curr_file.append("s_");
            if (stop)
              curr_file.append("stop_");
            curr_file.append(algorithms[j]);
            curr_file.append("_seed_");
            curr_file.append(std::to_string(seed));
            curr_file.append("_");
            curr_file.append(instances[i]);
            // std::cout << curr_file << std::endl;

            // std::fstream input;
            input.open(curr_file.c_str(), std::fstream::in);

            if (!input.is_open())
            {
              std::cout << "Could not open file " << curr_file << std::endl;
              throw 4;
            }

            std::stringstream s_lb1, s_t1, s_max_improve_iter;

            getline(input, line);
            size_t pos = line.find_first_of(":");
            status = line.substr(pos + 2);

            std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
            status.erase(end_pos, status.end());

            // std::cout << status << std::endl;
            getline(input, line);
            // getline(input,line);
            pos = line.find_first_of(":");
            s_lb1 << line.substr(pos + 2);
            s_lb1 >> lb;

            if (!(double_equals(lb, -1)) && (status != "INFEASIBLE") && (status != "POSSIBLYFEASIBLE") && (status != "FAILEDTOFINDAFEASIBLESOLUTION"))
            {
              avg_lb += lb;

              if (double_greater(lb, best_heuristic_lb))
                best_heuristic_lb = lb;

              /*if(double_greater(lb,best_lb))
              {
                std::cout << instances[i] << " seed " << seed << std::endl;
                std::cout << lb << " > " << best_lb << std::endl;
              }*/
              count++;
            }

            input.close();
          }

          if (count > 0)
            heuristic_lb = avg_lb / (1.0 * count);
          else
            heuristic_lb = -1.0;

          if (!((l == 0) && (j == 0)))
            output << " &";

          if (double_equals(heuristic_lb, -1.0))
            output << "& -- & -- ";
          else
            output << "& " << heuristic_lb << " & " << best_heuristic_lb;
        }
      }
      output << "\\\\" << std::endl;
    }
    if (k != dirs.size() - 1)
      output << "\\\\" << std::endl;
  }

  output.close();
}

void GenerateALNSCSVTableFridman(std::vector<std::string> dirs, bool stop, int type, int num_executions)
{
  // std::cout << "3" << std::endl; getchar(); getchar();
  std::string curr_file, algo_best_primal_bound = "cb_csc3";
  std::vector<std::string> algorithms;

  // algorithms.push_back("csc3_initial_primal_bound");
  // algorithms.push_back("cb_csc3");
  // algorithms.push_back("relax_cb_angle_0.03_clique_ccs_lcis_csc3");
  // algorithms.push_back("fp_alpha2_0_p_fixed_10");
  // algorithms.push_back("fp_alpha2_100_p_fixed_10");

  // algorithms.push_back("fp_alpha2_0_p_fixed_10_PR_ALNS_1000");
  // algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_1000");
  // algorithms.push_back("fp_alpha2_100_p_fixed_10_PR_ALNS_1000");
  // algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_1000");

  /*algorithms.push_back("fp_alpha2_0_p_fixed_10_PR_ALNS_2000");
  algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_2000");
  algorithms.push_back("fp_alpha2_100_p_fixed_10_PR_ALNS_2000");
  algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_2000");*/

  // algorithms.push_back("fp_alpha2_0_p_fixed_10_PR_ALNS_5000");
  // algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_5000");
  // algorithms.push_back("fp_alpha2_100_p_fixed_10_PR_ALNS_5000");
  // algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_5000");

  // algorithms.push_back("fp_alpha2_0_p_fixed_10_PR_ALNS_5000");
  // algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_10000");
  // algorithms.push_back("fp_alpha2_100_p_fixed_10_PR_ALNS_5000");
  // algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_10000");

  /*algorithms.push_back("fp_alpha2_0_p_fixed_10_PR_ALNS_1000");
  algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_1000");
  algorithms.push_back("fp_alpha2_100_p_fixed_10_PR_ALNS_1000");
  algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_1000");*/

  // algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_1000");
  // algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_1000");
  /*algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_2000");
  algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_2000");
  algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_5000");
  algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_5000");*/

  algorithms.push_back("fp_alpha2_0_p_fixed_10_ALNS_CONV");
  algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_ALNS_CONV");
  algorithms.push_back("fp_alpha2_100_p_fixed_10_ALNS_CONV");
  algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_ALNS_CONV");

  std::fstream output;
  std::string output_name;
  if (stop)
    output_name = ".//tables//CSV//stop_table_ALNS_friedman_" + std::string(1, type) + ".csv";
  else
    output_name = ".//tables//CSV//table_ALNS_friedman_" + std::string(1, type) + ".csv";
  output.open(output_name.c_str(), std::fstream::out);

  if (!output.is_open())
  {
    std::cout << "Could not open file " << output_name << std::endl;
    throw 1;
  }

  // std::cout << output_name << std::endl;

  output << std::setprecision(2) << std::fixed;

  for (size_t j = 0; j < algorithms.size(); j++)
  {
    output << algorithms[j];
    if (j < algorithms.size() - 1)
      output << ",";
  }

  output << std::endl;

  for (size_t k = 0; k < dirs.size(); ++k)
  {
    std::vector<std::string> instances;
    std::string folder = dirs[k].substr(20);
    AddInstancesFromDirectory(dirs[k], instances, false);
    for (size_t i = 0; i < instances.size(); ++i)
    {
      bool feasible = true;
      double heuristic_lb = 0.0, best_lb = 0.0, lb = 0.0, best_primal_gap = 0.0, curr_primal_gap = 0.0, best_heuristic_lb = 0.0, curr_max_improve_iter = 0;

      std::fstream input;
      std::string status;
      std::string line;

      // retrieve best primal bound (2h of execution)
      curr_file = ".//solutions//";
      curr_file.append(folder);
      curr_file.append("s_");
      if (stop)
        curr_file.append("stop_");
      curr_file.append(algo_best_primal_bound);
      curr_file.append("_");
      curr_file.append(instances[i]);
      // std::cout << curr_file << std::endl;

      // std::fstream input;
      input.open(curr_file.c_str(), std::fstream::in);

      if (!input.is_open())
      {
        std::cout << "Could not open file " << curr_file << std::endl;
        throw 4;
      }

      std::stringstream s_lb;
      lb = 0.0;

      getline(input, line);
      size_t pos = line.find_first_of(":");
      status = line.substr(pos + 2);

      std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
      status.erase(end_pos, status.end());

      // std::cout << status << std::endl;
      getline(input, line);
      getline(input, line);
      pos = line.find_first_of(":");
      s_lb << line.substr(pos + 2);
      if (s_lb.str() == "-inf")
        lb = -1.0;
      else
        s_lb >> lb;
      // std::cout << s_lb.str() << std::endl;
      // std::cout << lb << std::endl;

      if (status == "INFEASIBLE")
        feasible = false;

      if ((double_equals(lb, -1)) || (status == "INFEASIBLE"))
      {
        best_lb = -1.0;
      }
      else
        best_lb = lb;

      // std::cout << "best_lb: " << best_lb << std::endl;
      input.close();

      for (size_t j = 0; j < algorithms.size(); j++)
      {
        int count = 0;
        double avg_lb = 0.0, avg_t1 = 0.0, avg_profits = 0.0, avg_max_improve_iter = 0.0;

        best_heuristic_lb = 0.0;

        for (int seed = 1; seed <= 10; ++seed)
        {
          lb = 0.0;
          double t1 = 0.0;
          // retrieve primal bound at root node
          curr_file = ".//solutions//";
          curr_file.append(folder);
          curr_file.append("s_");
          if (stop)
            curr_file.append("stop_");
          curr_file.append(algorithms[j]);
          curr_file.append("_seed_");
          curr_file.append(std::to_string(seed));
          curr_file.append("_");
          curr_file.append(instances[i]);
          // std::cout << curr_file << std::endl;

          // std::fstream input;
          input.open(curr_file.c_str(), std::fstream::in);

          if (!input.is_open())
          {
            std::cout << "Could not open file " << curr_file << std::endl;
            throw 4;
          }

          std::stringstream s_lb1, s_t1, s_max_improve_iter;

          getline(input, line);
          size_t pos = line.find_first_of(":");
          status = line.substr(pos + 2);

          std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
          status.erase(end_pos, status.end());

          // std::cout << status << std::endl;
          getline(input, line);
          // getline(input,line);
          pos = line.find_first_of(":");
          s_lb1 << line.substr(pos + 2);
          s_lb1 >> lb;

          if (!(double_equals(best_lb, -1)) && !(double_equals(lb, -1)) && (status != "INFEASIBLE") && (status != "POSSIBLYFEASIBLE") && (status != "FAILEDTOFINDAFEASIBLESOLUTION"))
          {
            avg_lb += lb;
            avg_profits += lb;

            if (double_greater(lb, best_heuristic_lb))
              best_heuristic_lb = lb;

            /*if(double_greater(lb,best_lb))
            {
              std::cout << instances[i] << " seed " << seed << std::endl;
              std::cout << lb << " > " << best_lb << std::endl;
            }*/
            count++;
          }

          getline(input, line);
          getline(input, line);
          getline(input, line);

          pos = line.find_first_of(":");
          s_t1 << line.substr(pos + 2);
          s_t1 >> t1;
          avg_t1 += t1;

          getline(input, line);
          pos = line.find_first_of(":");
          s_max_improve_iter << line.substr(pos + 2);
          s_max_improve_iter >> curr_max_improve_iter;
          // std::cout << "root_lb: " << root_lb << std::endl;

          avg_max_improve_iter += curr_max_improve_iter;
          input.close();
        }

        if (count > 0)
          heuristic_lb = avg_lb / (1.0 * count);
        else
          heuristic_lb = -1.0;

        avg_t1 /= 10.0;
        avg_max_improve_iter /= (1.0 * num_executions);
        // if(count > 0) avg_profits /= (1.0*count);
        avg_profits /= (1.0 * num_executions);

        if (double_greater(avg_profits, best_heuristic_lb))
          std::cout << avg_profits << " > " << best_heuristic_lb << std::endl;

        if (!feasible)
        {
          curr_primal_gap = 0.0;
          best_primal_gap = 0.0;
        }
        else if ((double_equals(best_lb, -1)) || (double_equals(heuristic_lb, -1)))
        {
          curr_primal_gap = 100.0;
          best_primal_gap = 100.0;
        }
        else
        {
          if (double_equals(best_lb, heuristic_lb))
            curr_primal_gap = 0.0;
          else
            curr_primal_gap = (100 * (best_lb - heuristic_lb)) / best_lb;
          if (double_equals(best_lb, best_heuristic_lb))
            best_primal_gap = 0.0;
          else
            best_primal_gap = (100 * (best_lb - best_heuristic_lb)) / best_lb;
        }

        // std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

        switch (type)
        {
        case '1':
          output << (100.0 - curr_primal_gap);
          break;
        case '2':
          output << (100.0 - best_primal_gap);
          break;
        case '3':
          output << (5000.0 - avg_max_improve_iter);
          break;
        }
        if (j < algorithms.size() - 1)
          output << ",";
      }

      output << std::endl;
      // getchar();getchar();
    }
  }
}

void GenerateALNSCSVTableFridmanGaps(std::vector<std::string> dirs, bool stop, int num_executions = 10)
{
  GenerateALNSCSVTableFridman(dirs, stop, '1', num_executions);
}

void GenerateALNSCSVTableFridmanBestGaps(std::vector<std::string> dirs, bool stop, int num_executions = 10)
{
  GenerateALNSCSVTableFridman(dirs, stop, '2', num_executions);
}

void GenerateALNSCSVTableFridmanIterations(std::vector<std::string> dirs, bool stop, int num_executions = 10)
{
  GenerateALNSCSVTableFridman(dirs, stop, '3', num_executions);
}

void GenerateCSVTableFeasibilityPumpFriedmanGaps(std::vector<std::string> dirs, bool stop, int num_seeds = 10)
{
  // std::cout << "3" << std::endl; getchar(); getchar();
  std::string curr_file, algo_best_primal_bound = "cb_csc3";
  std::vector<std::string> algorithms;

  // algorithms.push_back("csc3_initial_primal_bound");
  // algorithms.push_back("cb_csc3");
  // algorithms.push_back("relax_cb_angle_0.03_clique_ccs_lcis_csc3");
  // algorithms.push_back("fp_alpha2_0_p_fixed_10");
  // algorithms.push_back("fp_alpha2_100_p_fixed_10");

  algorithms.push_back("fp_alpha2_0_p_fixed_10");
  algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10");

  algorithms.push_back("fp_alpha2_100_p_fixed_10");
  algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10");

  std::fstream output;
  std::string output_name;
  if (stop)
    output_name = ".//tables//CSV//stop_table_fp_fried_gaps.csv";
  else
    output_name = ".//tables//CSV//table_fp_fried_gaps.csv";
  output.open(output_name.c_str(), std::fstream::out);

  if (!output.is_open())
  {
    std::cout << "Could not open file " << output_name << std::endl;
    throw 1;
  }

  // std::cout << output_name << std::endl;

  output << std::setprecision(2) << std::fixed;

  for (size_t j = 0; j < algorithms.size(); j++)
  {
    output << algorithms[j];
    if (j < algorithms.size() - 1)
      output << ",";
  }

  output << std::endl;

  for (size_t k = 0; k < dirs.size(); ++k)
  {
    std::vector<std::string> instances;
    std::string folder = dirs[k].substr(20);
    AddInstancesFromDirectory(dirs[k], instances, false);
    for (size_t i = 0; i < instances.size(); ++i)
    {
      bool feasible = true;
      double heuristic_lb = 0.0, best_lb = 0.0, lb = 0.0, curr_primal_gap = 0.0, best_heuristic_lb = 0.0;

      std::fstream input;
      std::string status;
      std::string line;

      // retrieve best primal bound (2h of execution)
      curr_file = ".//solutions//";
      curr_file.append(folder);
      curr_file.append("s_");
      if (stop)
        curr_file.append("stop_");
      curr_file.append(algo_best_primal_bound);
      curr_file.append("_");
      curr_file.append(instances[i]);
      // std::cout << curr_file << std::endl;

      // std::fstream input;
      input.open(curr_file.c_str(), std::fstream::in);

      if (!input.is_open())
      {
        std::cout << "Could not open file " << curr_file << std::endl;
        throw 4;
      }

      std::stringstream s_lb;
      lb = 0.0;

      getline(input, line);
      size_t pos = line.find_first_of(":");
      status = line.substr(pos + 2);

      std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
      status.erase(end_pos, status.end());

      // std::cout << status << std::endl;
      getline(input, line);
      getline(input, line);
      pos = line.find_first_of(":");
      s_lb << line.substr(pos + 2);
      if (s_lb.str() == "-inf")
        lb = -1.0;
      else
        s_lb >> lb;
      // std::cout << s_lb.str() << std::endl;
      // std::cout << lb << std::endl;

      if (status == "INFEASIBLE")
        feasible = false;

      if ((double_equals(lb, -1)) || (status == "INFEASIBLE"))
      {
        best_lb = -1.0;
      }
      else
        best_lb = lb;

      // std::cout << "best_lb: " << best_lb << std::endl;
      input.close();

      for (size_t j = 0; j < algorithms.size(); j++)
      {
        int count = 0;
        double avg_lb = 0.0, avg_t1 = 0.0, avg_profits = 0.0;

        // best_heuristic_lb = 0.0;

        for (int seed = 1; seed <= num_seeds; ++seed)
        {
          lb = 0.0;
          double t1 = 0.0;
          // retrieve primal bound at root node
          curr_file = ".//solutions//";
          curr_file.append(folder);
          curr_file.append("s_");
          if (stop)
            curr_file.append("stop_");
          curr_file.append(algorithms[j]);
          curr_file.append("_seed_");
          curr_file.append(std::to_string(seed));
          curr_file.append("_");
          curr_file.append(instances[i]);
          // std::cout << curr_file << std::endl;

          // std::fstream input;
          input.open(curr_file.c_str(), std::fstream::in);

          if (!input.is_open())
          {
            std::cout << "Could not open file " << curr_file << std::endl;
            throw 4;
          }

          std::stringstream s_lb1, s_t1;

          getline(input, line);
          size_t pos = line.find_first_of(":");
          status = line.substr(pos + 2);

          std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
          status.erase(end_pos, status.end());

          // std::cout << status << std::endl;
          getline(input, line);
          // getline(input,line);
          pos = line.find_first_of(":");
          s_lb1 << line.substr(pos + 2);
          s_lb1 >> lb;

          if (!(double_equals(best_lb, -1)) && !(double_equals(lb, -1)) && (status != "INFEASIBLE") && (status != "POSSIBLYFEASIBLE") && (status != "FAILEDTOFINDAFEASIBLESOLUTION"))
          {
            avg_lb += lb;
            avg_profits += lb;

            count++;
          }

          input.close();
        }

        if (count > 0)
          heuristic_lb = avg_lb / (1.0 * count);
        else
          heuristic_lb = -1.0;

        // if(count > 0) avg_profits /= (1.0*count);
        avg_profits /= (10.0);

        // if(double_greater(avg_profits,best_heuristic_lb)) std::cout << avg_profits << " > " << best_heuristic_lb << std::endl;

        if (!feasible)
        {
          curr_primal_gap = 0.0;
        }
        else if ((double_equals(best_lb, -1)) || (double_equals(heuristic_lb, -1)))
        {
          curr_primal_gap = 100.0;
        }
        else
        {
          curr_primal_gap = (100 * (best_lb - heuristic_lb)) / best_lb;
        }

        // std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

        output << (100.0 - curr_primal_gap);
        if (j < algorithms.size() - 1)
          output << ",";
      }

      output << std::endl;

      // getchar();getchar();
    }
  }

  output.close();
}

void GenerateFeasibilityPumpPerInstanceResultsCSV(std::vector<std::string> dirs, bool stop, int num_seeds = 10)
{
  // std::cout << "3" << std::endl; getchar(); getchar();
  std::string curr_file;
  std::vector<std::string> algorithms;

  // algorithms.push_back("csc3_initial_primal_bound");
  // algorithms.push_back("cb_csc3");
  // algorithms.push_back("relax_cb_angle_0.03_clique_ccs_lcis_csc3");

  // algorithms.push_back("fp_alpha2_0_p_fixed_10");
  algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10");
  // algorithms.push_back("fp_alpha2_100_p_fixed_10");
  algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10");

  std::fstream output;
  std::string output_name;
  if (stop)
    output_name = ".//tables//CSV//stop_table_FP_detailed_per_instance.csv";
  else
    output_name = ".//tables//CSV//table_FP_detailed_per_instance.csv";
  output.open(output_name.c_str(), std::fstream::out);

  if (!output.is_open())
  {
    std::cout << "Could not open file " << output_name << std::endl;
    throw 1;
  }

  output << "Instance;";

  for (size_t i = 0; i < algorithms.size(); ++i)
  {
    for (size_t j = 1; j <= 4; ++j)
      output << algorithms[i] << ";";
  }

  output << std::endl;

  output << ";";

  for (size_t i = 0; i < algorithms.size(); ++i)
  {
    output << "avg lb;best lb;avg time(s);avg iterations;";
  }

  output << std::endl;
  // std::cout << output_name << std::endl;

  output << std::setprecision(2) << std::fixed;

  for (size_t k = 0; k < dirs.size(); ++k)
  {
    std::vector<std::string> instances;
    std::string folder = dirs[k].substr(20);
    AddInstancesFromDirectory(dirs[k], instances, false);

    std::sort(instances.begin(), instances.end());

    for (size_t i = 0; i < instances.size(); ++i)
    {
      bool feasible = true;
      double heuristic_lb = 0.0, lb = 0.0;

      if (stop)
        output << instances[i].substr(0, instances[i].size() - 4) << "_5%;";
      else
        output << instances[i].substr(0, instances[i].size() - 4) << ";";
      std::cout << instances[i] << std::endl;

      std::fstream input;
      std::string status;
      std::string line;

      for (size_t j = 0; j < algorithms.size(); j++)
      {
        int count = 0;
        double avg_lb = 0.0, avg_t1 = 0.0, avg_t2 = 0.0, avg_num_iter = 0.0, best_heuristic_lb = -1.0;

        for (int seed = 1; seed <= num_seeds; ++seed)
        {
          lb = 0.0;
          double t1 = 0.0, t2 = 0.0, num_iter = 0.0;
          // retrieve primal bound at root node
          curr_file = ".//solutions//";
          curr_file.append(folder);
          curr_file.append("s_");
          if (stop)
            curr_file.append("stop_");
          curr_file.append(algorithms[j]);
          curr_file.append("_seed_");
          curr_file.append(std::to_string(seed));
          curr_file.append("_");
          curr_file.append(instances[i]);
          // std::cout << curr_file << std::endl;

          // std::fstream input;
          input.open(curr_file.c_str(), std::fstream::in);

          if (!input.is_open())
          {
            std::cout << "Could not open file " << curr_file << std::endl;
            throw 4;
          }

          std::stringstream s_lb1, s_t1, s_t2, s_num_iter;

          getline(input, line);
          size_t pos = line.find_first_of(":");
          status = line.substr(pos + 2);

          std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
          status.erase(end_pos, status.end());

          // std::cout << status << std::endl;
          getline(input, line);
          // getline(input,line);
          pos = line.find_first_of(":");
          s_lb1 << line.substr(pos + 2);
          s_lb1 >> lb;

          if (!(double_equals(lb, -1)) && (status != "INFEASIBLE") && (status != "POSSIBLYFEASIBLE") && (status != "FAILEDTOFINDAFEASIBLESOLUTION"))
          {
            if (double_greater(lb, best_heuristic_lb))
              best_heuristic_lb = lb;
            avg_lb += lb;
            count++;
          }

          getline(input, line);
          getline(input, line);
          getline(input, line);
          getline(input, line);
          getline(input, line);

          pos = line.find_first_of(":");
          s_t1 << line.substr(pos + 2);
          s_t1 >> t1;

          getline(input, line);
          getline(input, line);
          pos = line.find_first_of(":");
          s_num_iter << line.substr(pos + 2);
          s_num_iter >> num_iter;

          getline(input, line);
          getline(input, line);
          getline(input, line);
          pos = line.find_first_of(":");
          s_t2 << line.substr(pos + 2);
          s_t2 >> t2;

          // std::cout << "root_lb: " << root_lb << std::endl;

          avg_t1 += t1;
          avg_t2 += t2;
          avg_num_iter += (num_iter + 1); // + 1 from stage 1
          input.close();
        }

        if (count > 0)
          heuristic_lb = avg_lb / (1.0 * count);
        else
          heuristic_lb = -1.0;

        avg_t1 /= 10.0;
        avg_t2 /= 10.0;
        avg_num_iter /= 10.0;

        if (double_equals(heuristic_lb, -1.0))
          output << " - ; - ;";
        else
          output << heuristic_lb << ";" << best_heuristic_lb << ";";
        output << avg_t1 + avg_t2 << ";" << avg_num_iter << ";";
      }

      // getchar();getchar();

      output << std::endl;
    }
  }

  output.close();
}

void GenerateCutsConfigurationsPerInstanceResultsCSV(std::vector<std::string> dirs, bool stop)
{
  std::string curr_file;
  std::vector<std::string> algorithms;

  algorithms.push_back("relax_FBs_cb_csc3");
  algorithms.push_back("relax_FBs_GCCs_cb_csc3");
  algorithms.push_back("relax_FBs_CCs_cb_csc3");
  algorithms.push_back("relax_FBs_CCCs_cb_csc3");
  algorithms.push_back("relax_FBs_LCIs_cb_csc3");
  algorithms.push_back("relax_FBs_AVICs_cb_csc3");

  algorithms.push_back("relax_FBs_GCCs_CCs_cb_csc3");
  algorithms.push_back("relax_FBs_GCCs_CCs_LCIs_cb_csc3");
  algorithms.push_back("relax_FBs_GCCs_CCs_LCIs_AVICs_cb_csc3");
  algorithms.push_back("relax_FBs_CCCs_LCIs_cb_csc3");
  algorithms.push_back("relax_FBs_CCCs_LCIs_AVICs_cb_csc3");

  std::fstream output;
  std::string output_name;
  if (stop)
    output_name = ".//tables//CSV//stop_table_LP_bounds.csv";
  else
    output_name = ".//tables//CSV//table_LP_bounds.csv";
  output.open(output_name.c_str(), std::fstream::out);

  if (!output.is_open())
  {
    std::cout << "Could not open file " << output_name << std::endl;
    throw 1;
  }

  // std::cout << output_name << std::endl;

  output << "Configuration;";

  for (size_t i = 0; i < algorithms.size(); ++i)
  {
    output << i << ";";
  }

  output << std::endl;

  output << "Instance;";

  for (size_t i = 0; i < algorithms.size(); ++i)
    output << "upper bound;";

  output << std::endl;
  output << std::setprecision(2) << std::fixed;

  for (size_t j = 0; j < dirs.size(); j++)
  {

    std::vector<std::string> instances;
    std::string folder = dirs[j].substr(20);
    AddInstancesFromDirectory(dirs[j], instances, false);

    std::sort(instances.begin(), instances.end());

    for (size_t i = 0; i < instances.size(); i++)
    {
      std::cout << instances[i] << std::endl;

      if (stop)
        output << instances[i].substr(0, instances[i].size() - 4) << "_5%;";
      else
        output << instances[i].substr(0, instances[i].size() - 4) << ";";

      for (size_t j = 0; j < algorithms.size(); j++)
      {
        double curr_improvement = 0.0;
        curr_file = ".//solutions//";
        curr_file.append(folder);
        curr_file.append("s_");
        if (stop)
          curr_file.append("stop_");
        curr_file.append(algorithms[j]);
        curr_file.append("_");
        curr_file.append(instances[i]);
        // std::cout << curr_file << std::endl;

        std::fstream input;
        input.open(curr_file.c_str(), std::fstream::in);

        if (!input.is_open())
        {
          std::cout << "Could not open file " << curr_file << std::endl;
          throw 4;
        }

        std::stringstream s_lp;
        std::string status;
        std::string line;

        getline(input, line);
        size_t pos = line.find_first_of(":");
        status = line.substr(pos + 2);

        std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
        status.erase(end_pos, status.end());

        getline(input, line);
        pos = line.find_first_of(":");
        s_lp << line.substr(pos + 2);

        if (status == "INFEASIBLE")
          output << " - ;";
        else
          output << s_lp.str() << ";";

        // output << s_lp.str() << ";";
        input.close();
      }

      // getchar();getchar();
      output << std::endl;
    }
  }

  output.close();
}

void GenerateLNSPerInstanceResultsCSV(std::vector<std::string> dirs, bool stop, int num_seeds = 10)
{
  // std::cout << "3" << std::endl; getchar(); getchar();
  std::string curr_file;
  std::vector<std::string> algorithms;

  // algorithms.push_back("fp_alpha2_0_p_fixed_10_ALNS_1000");
  algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_ALNS_1000");
  // algorithms.push_back("fp_alpha2_100_p_fixed_10_ALNS_1000");
  algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_ALNS_1000");

  // algorithms.push_back("fp_alpha2_0_p_fixed_10_ALNS_2000");
  algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_ALNS_2000");
  // algorithms.push_back("fp_alpha2_100_p_fixed_10_ALNS_2000");
  algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_ALNS_2000");

  // algorithms.push_back("fp_alpha2_0_p_fixed_10_PR_ALNS_5000");
  algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_5000");
  // algorithms.push_back("fp_alpha2_100_p_fixed_10_PR_ALNS_5000");
  algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_5000");

  // algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_100000");
  // algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_100000");

  algorithms.push_back("fp_alpha2_0_p_fixed_10_ALNS_CONV");
  algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_ALNS_CONV");
  algorithms.push_back("fp_alpha2_100_p_fixed_10_ALNS_CONV");
  algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_ALNS_CONV");

  std::fstream output;
  std::string output_name;
  if (stop)
    output_name = ".//tables//CSV//stop_table_LNS_detailed_per_instance.csv";
  else
    output_name = ".//tables//CSV//table_LNS_detailed_per_instance.csv";
  output.open(output_name.c_str(), std::fstream::out);

  if (!output.is_open())
  {
    std::cout << "Could not open file " << output_name << std::endl;
    throw 1;
  }

  output << "Instance;";

  for (size_t i = 0; i < algorithms.size(); ++i)
  {
    for (size_t j = 1; j <= 4; ++j)
      output << algorithms[i] << ";";
  }

  output << std::endl;

  output << ";";

  for (size_t i = 0; i < algorithms.size(); ++i)
  {
    output << "avg lb;best lb;avg time(s);avg max improve iter;";
  }

  output << std::endl;
  // std::cout << output_name << std::endl;

  output << std::setprecision(2) << std::fixed;

  for (size_t k = 0; k < dirs.size(); ++k)
  {
    std::vector<std::string> instances;
    std::string folder = dirs[k].substr(20);
    AddInstancesFromDirectory(dirs[k], instances, false);

    std::sort(instances.begin(), instances.end());

    for (size_t i = 0; i < instances.size(); ++i)
    {
      bool feasible = true;
      double heuristic_lb = 0.0, lb = 0.0;

      if (stop)
        output << instances[i].substr(0, instances[i].size() - 4) << "_5%;";
      else
        output << instances[i].substr(0, instances[i].size() - 4) << ";";

      std::cout << instances[i] << std::endl;

      std::fstream input;
      std::string status;
      std::string line;

      for (size_t j = 0; j < algorithms.size(); j++)
      {
        int count = 0;
        double avg_lb = 0.0, avg_t1 = 0.0, avg_profits = 0.0, avg_max_improve_iter = 0.0;

        double best_heuristic_lb = 0.0;

        for (int seed = 1; seed <= 10; ++seed)
        {
          lb = 0.0;
          double t1 = 0.0;
          double curr_max_improve_iter = 0.0;
          // retrieve primal bound at root node
          curr_file = ".//solutions//";
          curr_file.append(folder);
          curr_file.append("s_");
          if (stop)
            curr_file.append("stop_");
          curr_file.append(algorithms[j]);
          curr_file.append("_seed_");
          curr_file.append(std::to_string(seed));
          curr_file.append("_");
          curr_file.append(instances[i]);
          // std::cout << curr_file << std::endl;

          // std::fstream input;
          input.open(curr_file.c_str(), std::fstream::in);

          if (!input.is_open())
          {
            std::cout << "Could not open file " << curr_file << std::endl;
            throw 4;
          }

          std::stringstream s_lb1, s_t1, s_max_improve_iter;

          getline(input, line);
          size_t pos = line.find_first_of(":");
          status = line.substr(pos + 2);

          std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
          status.erase(end_pos, status.end());

          // std::cout << status << std::endl;
          getline(input, line);
          // getline(input,line);
          pos = line.find_first_of(":");
          s_lb1 << line.substr(pos + 2);
          s_lb1 >> lb;

          if (!(double_equals(lb, -1)) && (status != "INFEASIBLE") && (status != "POSSIBLYFEASIBLE") && (status != "FAILEDTOFINDAFEASIBLESOLUTION"))
          {
            avg_lb += lb;
            avg_profits += lb;

            if (double_greater(lb, best_heuristic_lb))
              best_heuristic_lb = lb;

            /*if(double_greater(lb,best_lb))
            {
              std::cout << instances[i] << " seed " << seed << std::endl;
              std::cout << lb << " > " << best_lb << std::endl;
            }*/
            count++;
          }

          getline(input, line);
          getline(input, line);
          getline(input, line);

          pos = line.find_first_of(":");
          s_t1 << line.substr(pos + 2);
          s_t1 >> t1;
          avg_t1 += t1;

          getline(input, line);
          pos = line.find_first_of(":");
          s_max_improve_iter << line.substr(pos + 2);
          s_max_improve_iter >> curr_max_improve_iter;
          // std::cout << "root_lb: " << root_lb << std::endl;

          avg_max_improve_iter += curr_max_improve_iter;
          input.close();
        }

        if (count > 0)
          heuristic_lb = avg_lb / (1.0 * count);
        else
          heuristic_lb = -1.0;

        avg_t1 /= 10.0;
        avg_max_improve_iter /= 10.0;
        // if(count > 0) avg_profits /= (1.0*count);
        avg_profits /= 10.0;

        if (double_equals(heuristic_lb, -1.0))
          output << " - ; - ;";
        else
          output << heuristic_lb << ";" << best_heuristic_lb << ";";
        output << avg_t1 << ";" << avg_max_improve_iter << ";";
      }
      // getchar();getchar();

      output << std::endl;
    }
  }

  output.close();
}

void GenerateFeasibilityPumpLatexTable(std::vector<std::string> dirs, bool stop, int num_seeds = 10)
{
  // std::cout << "3" << std::endl; getchar(); getchar();
  std::string curr_file, algo_best_primal_bound = "cb_csc3";
  std::vector<std::string> algorithms;

  // algorithms.push_back("csc3_initial_primal_bound");
  // algorithms.push_back("cb_csc3");
  // algorithms.push_back("relax_cb_angle_0.03_clique_ccs_lcis_csc3");

  // algorithms.push_back("fp_alpha2_0_p_fixed_10");
  algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10");
  // algorithms.push_back("fp_alpha2_100_p_fixed_10");
  algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10");

  std::fstream output;
  std::string output_name;
  if (stop)
    output_name = ".//tables//latex//stop_table_FP.txt";
  else
    output_name = ".//tables//latex//table_FP.txt";
  output.open(output_name.c_str(), std::fstream::out);

  if (!output.is_open())
  {
    std::cout << "Could not open file " << output_name << std::endl;
    throw 1;
  }

  // std::cout << output_name << std::endl;

  output << std::setprecision(2) << std::fixed;

  std::vector<std::vector<double>> total_primal_gap_per_algo(algorithms.size(), std::vector<double>());
  std::vector<double> total_avg_primal_gap_per_algo(algorithms.size(), 0.0);

  std::vector<std::vector<double>> total_iter_per_algo(algorithms.size(), std::vector<double>());
  std::vector<double> total_avg_iter_per_algo(algorithms.size(), 0.0);

  std::vector<std::vector<double>> total_time_per_algo1(algorithms.size(), std::vector<double>());
  std::vector<double> total_avg_time_per_algo1(algorithms.size(), 0.0);

  std::vector<std::vector<double>> total_time_per_algo2(algorithms.size(), std::vector<double>());
  std::vector<double> total_avg_time_per_algo2(algorithms.size(), 0.0);

  for (size_t k = 0; k < dirs.size(); ++k)
  {
    std::vector<std::vector<double>> primal_gap_per_algo(algorithms.size(), std::vector<double>());
    std::vector<double> avg_primal_gap_per_algo(algorithms.size(), 0.0);
    std::vector<double> st_dev_primal_gap_per_algo(algorithms.size(), 0.0);

    std::vector<std::vector<double>> time_per_algo1(algorithms.size(), std::vector<double>());
    std::vector<double> avg_time_per_algo1(algorithms.size(), 0.0);

    std::vector<std::vector<double>> time_per_algo2(algorithms.size(), std::vector<double>());
    std::vector<double> avg_time_per_algo2(algorithms.size(), 0.0);

    std::vector<std::vector<double>> iter_per_algo(algorithms.size(), std::vector<double>());
    std::vector<double> avg_iter_per_algo(algorithms.size(), 0.0);

    std::vector<std::string> instances;
    std::string folder = dirs[k].substr(20);
    AddInstancesFromDirectory(dirs[k], instances, false);
    for (size_t i = 0; i < instances.size(); ++i)
    {
      bool feasible = true;
      double heuristic_lb = 0.0, best_lb = 0.0, lb = 0.0, curr_primal_gap = 0.0;

      std::fstream input;
      std::string status;
      std::string line;

      // retrieve best primal bound (2h of execution)
      curr_file = ".//solutions//";
      curr_file.append(folder);
      curr_file.append("s_");
      if (stop)
        curr_file.append("stop_");
      curr_file.append(algo_best_primal_bound);
      curr_file.append("_");
      curr_file.append(instances[i]);
      // std::cout << curr_file << std::endl;

      // std::fstream input;
      input.open(curr_file.c_str(), std::fstream::in);

      if (!input.is_open())
      {
        std::cout << "Could not open file " << curr_file << std::endl;
        throw 4;
      }

      std::stringstream s_lb;
      lb = 0.0;

      getline(input, line);
      size_t pos = line.find_first_of(":");
      status = line.substr(pos + 2);

      std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
      status.erase(end_pos, status.end());

      // std::cout << status << std::endl;
      getline(input, line);
      getline(input, line);
      pos = line.find_first_of(":");
      s_lb << line.substr(pos + 2);
      if (s_lb.str() == "-inf")
        lb = -1.0;
      else
        s_lb >> lb;
      // std::cout << s_lb.str() << std::endl;
      // std::cout << lb << std::endl;

      if (status == "INFEASIBLE")
        feasible = false;

      if ((double_equals(lb, -1)) || (status == "INFEASIBLE"))
      {
        best_lb = -1.0;
      }
      else
        best_lb = lb;

      // std::cout << "best_lb: " << best_lb << std::endl;
      input.close();

      for (size_t j = 0; j < algorithms.size(); j++)
      {
        int count = 0;
        double avg_lb = 0.0, avg_t1 = 0.0, avg_t2 = 0.0, avg_num_iter = 0.0;

        for (int seed = 1; seed <= num_seeds; ++seed)
        {
          lb = 0.0;
          double t1 = 0.0, t2 = 0.0, num_iter = 0.0;
          // retrieve primal bound at root node
          curr_file = ".//solutions//";
          curr_file.append(folder);
          curr_file.append("s_");
          if (stop)
            curr_file.append("stop_");
          curr_file.append(algorithms[j]);
          curr_file.append("_seed_");
          curr_file.append(std::to_string(seed));
          curr_file.append("_");
          curr_file.append(instances[i]);
          // std::cout << curr_file << std::endl;

          // std::fstream input;
          input.open(curr_file.c_str(), std::fstream::in);

          if (!input.is_open())
          {
            std::cout << "Could not open file " << curr_file << std::endl;
            throw 4;
          }

          std::stringstream s_lb1, s_t1, s_t2, s_num_iter;

          getline(input, line);
          size_t pos = line.find_first_of(":");
          status = line.substr(pos + 2);

          std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
          status.erase(end_pos, status.end());

          // std::cout << status << std::endl;
          getline(input, line);
          // getline(input,line);
          pos = line.find_first_of(":");
          s_lb1 << line.substr(pos + 2);
          s_lb1 >> lb;

          if (!(double_equals(lb, -1)) && (status != "INFEASIBLE") && (status != "POSSIBLYFEASIBLE") && (status != "FAILEDTOFINDAFEASIBLESOLUTION"))
          {
            avg_lb += lb;
            count++;
          }

          getline(input, line);
          getline(input, line);
          getline(input, line);
          getline(input, line);
          getline(input, line);

          pos = line.find_first_of(":");
          s_t1 << line.substr(pos + 2);
          s_t1 >> t1;

          getline(input, line);
          getline(input, line);
          pos = line.find_first_of(":");
          s_num_iter << line.substr(pos + 2);
          s_num_iter >> num_iter;

          // if(num_iter > 1998){ std::cout << num_iter << std::endl; getchar();}

          getline(input, line);
          getline(input, line);
          getline(input, line);
          pos = line.find_first_of(":");
          s_t2 << line.substr(pos + 2);
          s_t2 >> t2;

          // std::cout << "root_lb: " << root_lb << std::endl;

          avg_t1 += t1;
          avg_t2 += t2;
          avg_num_iter += (num_iter + 1); // + 1 from stage 1
          input.close();
        }

        if (count > 0)
          heuristic_lb = avg_lb / (1.0 * count);
        else
          heuristic_lb = -1.0;

        avg_t1 /= 10.0;
        avg_t2 /= 10.0;
        avg_num_iter /= 10.0;

        if (!feasible)
        {
          curr_primal_gap = 0.0;
        }
        else if ((double_equals(best_lb, -1)) || (double_equals(heuristic_lb, -1)))
        {
          curr_primal_gap = 100.0;
        }
        else
        {
          if (double_equals(best_lb, heuristic_lb))
            curr_primal_gap = 0.0;
          else
            curr_primal_gap = (100 * (best_lb - heuristic_lb)) / best_lb;
        }

        // std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

        primal_gap_per_algo[j].push_back(curr_primal_gap);
        total_primal_gap_per_algo[j].push_back(curr_primal_gap);
        avg_primal_gap_per_algo[j] += curr_primal_gap;
        total_avg_primal_gap_per_algo[j] += curr_primal_gap;

        iter_per_algo[j].push_back(avg_num_iter);
        total_iter_per_algo[j].push_back(avg_num_iter);
        avg_iter_per_algo[j] += avg_num_iter;
        total_avg_iter_per_algo[j] += avg_num_iter;

        time_per_algo1[j].push_back(avg_t1);
        total_time_per_algo1[j].push_back(avg_t1);
        avg_time_per_algo1[j] += avg_t1;
        total_avg_time_per_algo1[j] += avg_t1;

        time_per_algo2[j].push_back(avg_t2);
        total_time_per_algo2[j].push_back(avg_t2);
        avg_time_per_algo2[j] += avg_t2;
        total_avg_time_per_algo2[j] += avg_t2;
      }

      // getchar();getchar();
    }

    output << k + 1;
    if (stop)
      output << "\\_5\\%";

    for (size_t h = 0; h < algorithms.size(); h++)
    {
      avg_primal_gap_per_algo[h] /= (1.0 * ((primal_gap_per_algo[h]).size()));
      st_dev_primal_gap_per_algo[h] = StDev(primal_gap_per_algo[h], avg_primal_gap_per_algo[h]);

      avg_iter_per_algo[h] /= (1.0 * ((iter_per_algo[h]).size()));

      avg_time_per_algo1[h] /= (1.0 * ((time_per_algo1[h]).size()));

      avg_time_per_algo2[h] /= (1.0 * ((time_per_algo2[h]).size()));

      output << " & & " << avg_primal_gap_per_algo[h] << " & " << st_dev_primal_gap_per_algo[h] << " & " << avg_iter_per_algo[h] << " & " << avg_time_per_algo1[h] + avg_time_per_algo2[h];
    }
    output << "\\\\" << std::endl;
  }

  output << "Total";
  for (size_t j = 0; j < algorithms.size(); j++)
  {
    // std::cout << total_improvement_per_algo[j].size() << std::endl;
    total_avg_primal_gap_per_algo[j] /= (1.0 * ((total_primal_gap_per_algo[j]).size()));
    total_avg_iter_per_algo[j] /= (1.0 * ((total_iter_per_algo[j]).size()));
    total_avg_time_per_algo1[j] /= (1.0 * ((total_time_per_algo1[j]).size()));
    total_avg_time_per_algo2[j] /= (1.0 * ((total_time_per_algo2[j]).size()));

    output << " & & " << total_avg_primal_gap_per_algo[j] << " & " << StDev(total_primal_gap_per_algo[j], total_avg_primal_gap_per_algo[j]) << " & " << total_avg_iter_per_algo[j] << " & " << total_avg_time_per_algo1[j] + total_avg_time_per_algo2[j];
  }

  output << "\\\\" << std::endl;

  output.close();
}

void GeneratePrimalAndDualGapsLatexTable(std::vector<std::string> dirs, bool stop)
{
  // std::cout << "3" << std::endl; getchar(); getchar();
  std::string curr_file, algo_initial_primal_bound = "fp_cuts_alpha2_0_p_fixed_10_ALNS_5000", algo_best_primal_bound = "cb_csc3";
  std::vector<std::string> algorithms;

  // algorithms.push_back("csc3_initial_primal_bound");
  // algorithms.push_back("cb_csc3");
  // algorithms.push_back("relax_cb_angle_0.03_clique_ccs_lcis_csc3");
  algorithms.push_back("relax_cb_avics_clique_ccs_lcis_csc3");

  std::fstream output;
  std::string output_name;
  if (stop)
    output_name = ".//tables//latex//stop_table_primal_dual_gaps_ALNS.txt";
  else
    output_name = ".//tables//latex//table_primal_dual_gaps_ALNS.txt";
  output.open(output_name.c_str(), std::fstream::out);

  if (!output.is_open())
  {
    std::cout << "Could not open file " << output_name << std::endl;
    throw 1;
  }

  // std::cout << output_name << std::endl;

  output << std::setprecision(2) << std::fixed;

  std::vector<std::vector<double>> total_primal_gap_per_algo(algorithms.size(), std::vector<double>());
  std::vector<std::vector<double>> total_dual_gap_per_algo(algorithms.size(), std::vector<double>());
  std::vector<double> total_avg_primal_gap(algorithms.size(), 0.0);
  std::vector<double> total_avg_dual_gap(algorithms.size(), 0.0);

  for (size_t j = 0; j < dirs.size(); j++)
  {
    std::vector<std::vector<double>> primal_gap_per_algo(algorithms.size(), std::vector<double>());
    std::vector<std::vector<double>> dual_gap_per_algo(algorithms.size(), std::vector<double>());
    std::vector<double> avg_primal_gap(algorithms.size(), 0.0);
    std::vector<double> avg_dual_gap(algorithms.size(), 0.0);
    std::vector<double> st_dev_primal_gap(algorithms.size(), 0.0);
    std::vector<double> st_dev_dual_gap(algorithms.size(), 0.0);

    std::vector<std::string> instances;
    std::string folder = dirs[j].substr(20);
    AddInstancesFromDirectory(dirs[j], instances, false);
    for (size_t i = 0; i < instances.size(); i++)
    {
      bool feasible = true;
      double root_lb = 0.0, root_ub = 0.0, best_lb = 0.0;
      double lb = 0.0, ub = 0.0, lp = 0.0;
      double avg_lb = 0.0;
      int count = 0;

      for (int seed = 1; seed <= 10; ++seed)
      {
        lb = 0.0;
        // retrieve primal bound at root node
        curr_file = ".//solutions//";
        curr_file.append(folder);
        curr_file.append("s_");
        if (stop)
          curr_file.append("stop_");
        curr_file.append(algo_initial_primal_bound);
        curr_file.append("_seed_");
        curr_file.append(std::to_string(seed));
        curr_file.append("_");
        curr_file.append(instances[i]);
        // std::cout << curr_file << std::endl;

        std::fstream input;
        input.open(curr_file.c_str(), std::fstream::in);

        if (!input.is_open())
        {
          std::cout << "Could not open file " << curr_file << std::endl;
          throw 4;
        }

        std::stringstream s_lb1;
        std::string status;
        std::string line;

        getline(input, line);
        size_t pos = line.find_first_of(":");
        status = line.substr(pos + 2);

        std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
        status.erase(end_pos, status.end());

        // std::cout << status << std::endl;
        getline(input, line);
        // getline(input,line);
        pos = line.find_first_of(":");
        s_lb1 << line.substr(pos + 2);
        if (s_lb1.str() == "-inf")
          lb = -1.0;
        else
          s_lb1 >> lb;

        if (!(double_equals(lb, -1)) && (status != "INFEASIBLE") && (status != "POSSIBLYFEASIBLE") && (status != "FAILEDTOFINDAFEASIBLESOLUTION"))
        {
          avg_lb += lb;
          count++;
        }

        // std::cout << "root_lb: " << root_lb << std::endl;
        input.close();
      }

      if (count > 0)
        root_lb = avg_lb / (1.0 * count);
      else
        root_lb = -1.0;

      std::fstream input;
      std::stringstream s_lb1;
      std::string status;
      std::string line;

      // retrieve best primal bound (2h of execution)
      curr_file = ".//solutions//";
      curr_file.append(folder);
      curr_file.append("s_");
      if (stop)
        curr_file.append("stop_");
      curr_file.append(algo_best_primal_bound);
      curr_file.append("_");
      curr_file.append(instances[i]);
      // std::cout << curr_file << std::endl;

      // std::fstream input;
      input.open(curr_file.c_str(), std::fstream::in);

      if (!input.is_open())
      {
        std::cout << "Could not open file " << curr_file << std::endl;
        throw 4;
      }

      std::stringstream s_lb;
      lb = 0.0;

      getline(input, line);
      size_t pos = line.find_first_of(":");
      status = line.substr(pos + 2);

      std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
      status.erase(end_pos, status.end());

      std::cout << status << std::endl;
      getline(input, line);
      getline(input, line);
      pos = line.find_first_of(":");
      s_lb << line.substr(pos + 2);
      if (s_lb.str() == "-inf")
        lb = -1.0;
      else
        s_lb >> lb;
      // std::cout << s_lb.str() << std::endl;
      // std::cout << lb << std::endl;

      if (status == "INFEASIBLE")
        feasible = false;

      if ((double_equals(lb, -1)) || (status == "INFEASIBLE"))
      {
        best_lb = -1.0;
      }
      else
        best_lb = lb;

      // std::cout << "best_lb: " << best_lb << std::endl;
      input.close();

      for (size_t j = 0; j < algorithms.size(); j++)
      {
        double curr_primal_gap = 0.0, curr_dual_gap = 0.0;
        ub = 0.0;
        root_ub = 0.0;
        curr_file = ".//solutions//";
        curr_file.append(folder);
        curr_file.append("s_");
        if (stop)
          curr_file.append("stop_");
        curr_file.append(algorithms[j]);
        curr_file.append("_");
        curr_file.append(instances[i]);
        // std::cout << curr_file << std::endl;

        // std::fstream input;
        input.open(curr_file.c_str(), std::fstream::in);

        if (!input.is_open())
        {
          std::cout << "Could not open file " << curr_file << std::endl;
          throw 4;
        }

        std::stringstream s_ub;
        getline(input, line);
        pos = line.find_first_of(":");
        status = line.substr(pos + 2);

        end_pos = std::remove(status.begin(), status.end(), ' ');
        status.erase(end_pos, status.end());

        // std::cout << status << std::endl;
        getline(input, line);
        pos = line.find_first_of(":");
        s_ub << line.substr(pos + 2);
        if (s_ub.str() == "inf")
          ub = -1;
        else
          s_ub >> ub;
        // std::cout << s_lp.str() << std::endl;
        // std::cout << lp << std::endl;
        // getchar(); getchar();

        if ((double_equals(ub, -1)) || (status == "INFEASIBLE"))
        {
          root_ub = -1.0;
        }
        else
          root_ub = ub;

        std::cout << "root_ub: " << root_ub << std::endl;

        if (!feasible)
        {
          curr_primal_gap = 0.0;
          curr_dual_gap = 0.0;
        }
        else if ((double_equals(root_ub, -1)) || (double_equals(best_lb, -1)) || (double_equals(root_lb, -1)))
        {
          curr_primal_gap = 100.0;
          curr_dual_gap = 100.0;
        }
        else
        {
          curr_primal_gap = (100 * (best_lb - root_lb)) / best_lb;
          curr_dual_gap = (100 * (root_ub - best_lb)) / best_lb;
        }

        // std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

        primal_gap_per_algo[j].push_back(curr_primal_gap);
        dual_gap_per_algo[j].push_back(curr_dual_gap);
        total_primal_gap_per_algo[j].push_back(curr_primal_gap);
        total_dual_gap_per_algo[j].push_back(curr_dual_gap);
        total_avg_primal_gap[j] += curr_primal_gap;
        total_avg_dual_gap[j] += curr_dual_gap;
        avg_primal_gap[j] += curr_primal_gap;
        avg_dual_gap[j] += curr_dual_gap;
        input.close();
      }

      // getchar();getchar();
    }

    output << j + 1;

    for (size_t k = 0; k < algorithms.size(); k++)
    {
      avg_primal_gap[k] /= (1.0 * (instances.size()));
      st_dev_primal_gap[k] = StDev(primal_gap_per_algo[k], avg_primal_gap[k]);
      avg_dual_gap[k] /= (1.0 * (instances.size()));
      st_dev_dual_gap[k] = StDev(dual_gap_per_algo[k], avg_dual_gap[k]);
      output << " & & " << avg_primal_gap[k] << " & " << st_dev_primal_gap[k] << " & & " << avg_dual_gap[k] << " & " << st_dev_dual_gap[k];
    }
    output << "\\\\" << std::endl;
  }

  output << "Total";
  for (size_t j = 0; j < algorithms.size(); j++)
  {
    // std::cout << total_improvement_per_algo[j].size() << std::endl;
    total_avg_primal_gap[j] /= (1.0 * ((total_primal_gap_per_algo[j]).size()));
    total_avg_dual_gap[j] /= (1.0 * ((total_dual_gap_per_algo[j]).size()));
    output << "& & " << total_avg_primal_gap[j] << " & " << StDev(total_primal_gap_per_algo[j], total_avg_primal_gap[j]) << " & & " << total_avg_dual_gap[j] << " & " << StDev(total_dual_gap_per_algo[j], total_avg_dual_gap[j]);
  }

  output << "\\\\" << std::endl;

  output.close();
}

void GenerateHeuristicSolutionImprovementLatexTable(std::vector<std::string> dirs, bool stop)
{
  // std::cout << "3" << std::endl; getchar(); getchar();
  std::string curr_file, algo_initial_primal_bound = "cb_avics_clique_ccs_lcis_csc3_primal_bound_nodelim_1";
  std::vector<std::string> algorithms;

  // algorithms.push_back("csc3_initial_primal_bound");
  // algorithms.push_back("cb_csc3");
  // algorithms.push_back("relax_cb_angle_0.03_clique_ccs_lcis_csc3");
  /*algorithms.push_back("fp_alpha2_0_p_fixed_10");
  algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10");
  algorithms.push_back("fp_alpha2_100_p_fixed_10");*/
  algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10");
  algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_ALNS_2000");
  algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_MT_ALNS_2000");
  algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_ALNS_5000");
  algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_MT_ALNS_5000");

  std::fstream output;
  std::string output_name;
  if (stop)
    output_name = ".//tables//latex//stop_table_primal_heuristic_improvements.txt";
  else
    output_name = ".//tables//latex//table_primal_heuristic_improvements.txt";
  output.open(output_name.c_str(), std::fstream::out);

  if (!output.is_open())
  {
    std::cout << "Could not open file " << output_name << std::endl;
    throw 1;
  }

  // std::cout << output_name << std::endl;

  output << std::setprecision(2) << std::fixed;

  std::vector<std::vector<double>> total_improvement_per_algo(algorithms.size(), std::vector<double>());
  std::vector<double> total_avg_improvement(algorithms.size(), 0.0);

  for (size_t j = 0; j < dirs.size(); j++)
  {
    std::vector<std::vector<double>> improvement_per_algo(algorithms.size(), std::vector<double>());
    std::vector<double> avg_improvement(algorithms.size(), 0.0);
    std::vector<double> st_dev(algorithms.size(), 0.0);

    std::vector<std::string> instances;
    std::string folder = dirs[j].substr(20);
    AddInstancesFromDirectory(dirs[j], instances, false);
    for (size_t i = 0; i < instances.size(); i++)
    {
      std::cout << "* " << instances[i] << std::endl;
      bool feasible = true;
      double root_lb = 0.0, lb = 0.0;

      // retrieve primal bound at root node
      curr_file = ".//solutions//";
      curr_file.append(folder);
      curr_file.append("s_");
      if (stop)
        curr_file.append("stop_");
      curr_file.append(algo_initial_primal_bound);
      curr_file.append("_");
      curr_file.append(instances[i]);
      // std::cout << curr_file << std::endl;

      std::fstream input;
      input.open(curr_file.c_str(), std::fstream::in);

      if (!input.is_open())
      {
        std::cout << "Could not open file " << curr_file << std::endl;
        throw 4;
      }

      std::stringstream s_lb1;
      std::string status;
      std::string line;

      getline(input, line);
      size_t pos = line.find_first_of(":");
      status = line.substr(pos + 2);

      std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
      status.erase(end_pos, status.end());

      // std::cout << status << std::endl;
      getline(input, line);
      getline(input, line);
      pos = line.find_first_of(":");
      s_lb1 << line.substr(pos + 2);
      s_lb1 >> lb;

      if ((double_equals(lb, -1)) || (status == "INFEASIBLE"))
      {
        root_lb = -1.0;
        feasible = false;
      }
      else
        root_lb = lb;

      // std::cout << "root_lb: " << root_lb << std::endl;
      input.close();

      if (feasible)
      {
        for (size_t j = 0; j < algorithms.size(); j++)
        {
          double curr_avg_improvement = 0.0;
          for (int seed = 1; seed <= 10; ++seed)
          {
            double curr_improvement = 0.0;

            curr_file = ".//solutions//";
            curr_file.append(folder);
            curr_file.append("s_");
            if (stop)
              curr_file.append("stop_");
            curr_file.append(algorithms[j]);
            curr_file.append("_seed_");
            curr_file += (std::to_string(seed) + "_");
            curr_file.append(instances[i]);

            std::fstream input;
            input.open(curr_file.c_str(), std::fstream::in);

            if (!input.is_open())
            {
              std::cout << "Could not open file " << curr_file << std::endl;
              throw 4;
            }

            std::stringstream s_bound;
            std::string status;
            double bound = 0.0;
            std::string line;

            getline(input, line);
            size_t pos = line.find_first_of(":");
            status = line.substr(pos + 2);

            std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
            status.erase(end_pos, status.end());

            getline(input, line);
            pos = line.find_first_of(":");
            s_bound << line.substr(pos + 2);
            s_bound >> bound;

            if ((double_equals(bound, -1)) || (status == "INFEASIBLE") || (double_equals(root_lb, -1)))
            {
              curr_improvement = -1;
              feasible = false;
            }
            else
              curr_improvement = (100 * (bound - root_lb)) / root_lb;

            // std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

            if (feasible)
            {
              // std::cout << bound << std::endl;
              curr_avg_improvement += curr_improvement;
            }
          }

          if (feasible)
          {
            curr_avg_improvement /= 10.0;

            // std::cout << "curr_avg_improv: " << curr_avg_improvement << std::endl;

            improvement_per_algo[j].push_back(curr_avg_improvement);
            total_improvement_per_algo[j].push_back(curr_avg_improvement);
            total_avg_improvement[j] += curr_avg_improvement;
            avg_improvement[j] += curr_avg_improvement;
          }
        }
      }
    }

    output << j + 1;

    for (size_t j = 0; j < algorithms.size(); j++)
    {
      avg_improvement[j] /= (1.0 * ((improvement_per_algo[j]).size()));
      st_dev[j] = StDev(improvement_per_algo[j], avg_improvement[j]);
      output << " & & " << avg_improvement[j] << " & " << st_dev[j];
    }
    output << "\\\\" << std::endl;
  }

  output << "Total";
  for (size_t j = 0; j < algorithms.size(); j++)
  {
    // std::cout << total_improvement_per_algo[j].size() << std::endl;
    total_avg_improvement[j] /= (1.0 * ((total_improvement_per_algo[j]).size()));
    output << "& & " << total_avg_improvement[j] << " & " << StDev(total_improvement_per_algo[j], total_avg_improvement[j]);
  }

  output << "\\\\" << std::endl;

  output.close();
}

void GenerateCSVTableALNS(std::vector<std::string> dirs, bool stop, int num_executions = 10)
{
  std::string curr_file;
  std::vector<std::string> algorithms;

  /*algorithms.push_back("fp_alpha2_0_p_25");
  algorithms.push_back("fp_alpha2_100_p_25");

  algorithms.push_back("fp_alpha2_0_p_15");
  algorithms.push_back("fp_alpha2_100_p_15");

  algorithms.push_back("fp_alpha2_0_p_10");
  algorithms.push_back("fp_alpha2_100_p_10");*/

  // algorithms.push_back("fp_alpha2_0_p_fixed_10_ALNS_5000");
  // algorithms.push_back("fp_alpha2_100_p_fixed_10_ALNS_5000");

  // algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_ALNS_5000");
  // algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_ALNS_2000");
  // algorithms.push_back("fp_alpha2_0_p_fixed_10_ALNS_5000");
  // algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_ALNS_5000");
  // algorithms.push_back("new_fp_cuts_alpha2_0_p_fixed_10_ALNS_5000");
  // algorithms.push_back("new2_fp_cuts_alpha2_0_p_fixed_10_ALNS_5000");

  /*algorithms.push_back("fp_alpha1_0_alpha2_0_p_25");
  algorithms.push_back("fp_alpha1_100_alpha2_0_p_25");
  algorithms.push_back("fp_alpha1_0_alpha2_100_p_25");
  algorithms.push_back("fp_alpha1_100_alpha2_100_p_25");*/

  algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_5000");
  algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_5000");
  // algorithms.push_back("fp_alpha2_0_p_fixed_10_MT_PR_ALNS_5000_BW");

  std::fstream output;
  std::string output_name;
  if (stop)
    output_name = ".//tables//CSV//stop_table_ALNS_correct.csv";
  else
    output_name = ".//tables//CSV//table_ALNS_correct.csv";
  output.open(output_name.c_str(), std::fstream::out);

  if (!output.is_open())
  {
    std::cout << "Could not open file " << output_name << std::endl;
    throw 1;
  }

  // std::cout << output_name << std::endl;

  output << std::setprecision(2) << std::fixed;

  double curr_num_iter = 0.0, curr_time = 0.0, curr_profits_sum = 0.0, curr_improve_iter = 0.0, curr_max_improve_iter = 0.0;

  output << "algo;";
  for (size_t j = 0; j < algorithms.size(); j++)
  {
    for (int i = 1; i <= 4; i++)
      output << algorithms[j] << ";";
  }
  output << std::endl;

  output << "instance;";
  for (size_t j = 0; j < algorithms.size(); j++)
  {
    output << "profits_sum;avg_improve_iter;max_improve_iter;time;";
  }
  output << std::endl;

  for (size_t j = 0; j < dirs.size(); j++)
  {
    std::vector<std::string> instances;
    std::string folder = dirs[j].substr(20);
    AddInstancesFromDirectory(dirs[j], instances, false);

    for (size_t i = 0; i < instances.size(); i++)
    {
      output << instances[i] << ";";
      double original_lp = 0.0;
      for (size_t j = 0; j < algorithms.size(); j++)
      {
        curr_num_iter = curr_time = curr_profits_sum = curr_improve_iter = curr_max_improve_iter = 0.0;
        for (size_t seed = 1; seed <= 10; seed++)
        {
          curr_file = ".//solutions//";
          curr_file.append(folder);
          curr_file.append("s_");
          if (stop)
            curr_file.append("stop_");
          curr_file.append(algorithms[j]);
          curr_file.append("_seed_");
          curr_file += (std::to_string(seed) + "_");
          curr_file.append(instances[i]);

          std::cout << curr_file << std::endl;
          // getchar();getchar();
          std::fstream input;
          input.open(curr_file.c_str(), std::fstream::in);

          if (!input.is_open())
          {
            std::cout << "Could not open file " << curr_file << std::endl;
            throw 4;
          }

          std::stringstream s_iter, s_t1, s_profits_sum, s_improve_iter;
          std::string status;
          double iter = 0.0, t1 = 0.0, profits_sum = 0.0, improve_iter = 0.0;
          std::string line;

          getline(input, line);
          size_t pos = line.find_first_of(":");
          status = line.substr(pos + 2);

          std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
          status.erase(end_pos, status.end());

          // std::cout << status << std::endl;

          getline(input, line);
          pos = line.find_first_of(":");
          s_profits_sum << line.substr(pos + 2);
          s_profits_sum >> profits_sum;

          // std::cout << "sum: " << profits_sum << std::endl;

          getline(input, line);
          getline(input, line);
          pos = line.find_first_of(":");
          s_iter << line.substr(pos + 2);
          s_iter >> iter;

          // std::cout << "iter1: " << iter1 << std::endl;

          getline(input, line);
          pos = line.find_first_of(":");
          s_t1 << line.substr(pos + 2);
          s_t1 >> t1;

          // std::cout << "t1: " << t1 << std::endl;

          getline(input, line);
          pos = line.find_first_of(":");
          s_improve_iter << line.substr(pos + 2);
          s_improve_iter >> improve_iter;

          // std::cout << "t2: " << t2 << std::endl;
          input.close();

          curr_profits_sum += profits_sum;
          curr_num_iter += iter;
          curr_time += t1;
          curr_improve_iter += improve_iter;

          if (double_greater(improve_iter, curr_max_improve_iter))
            curr_max_improve_iter = improve_iter;

          // getchar(); getchar();
        }
        output << curr_profits_sum / num_executions << ";" << curr_improve_iter / num_executions << ";" << curr_max_improve_iter << ";" << curr_time / num_executions << ";";
      }
      output << std::endl;
    }
  }

  output.close();
}

void GenerateCSVTableFeasibilityPumpFriedmanIterations(std::vector<std::string> dirs, bool stop, int num_executions = 10)
{
  std::string curr_file;
  std::vector<std::string> algorithms;

  /*algorithms.push_back("fp_alpha2_0_p_25");
  algorithms.push_back("fp_alpha2_100_p_25");

  algorithms.push_back("fp_alpha2_0_p_15");
  algorithms.push_back("fp_alpha2_100_p_15");

  algorithms.push_back("fp_alpha2_0_p_10");
  algorithms.push_back("fp_alpha2_100_p_10");*/

  algorithms.push_back("fp_alpha2_0_p_fixed_10");
  algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10");

  algorithms.push_back("fp_alpha2_100_p_fixed_10");
  algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10");

  /*algorithms.push_back("fp_alpha1_0_alpha2_0_p_25");
  algorithms.push_back("fp_alpha1_100_alpha2_0_p_25");
  algorithms.push_back("fp_alpha1_0_alpha2_100_p_25");
  algorithms.push_back("fp_alpha1_100_alpha2_100_p_25");*/

  std::fstream output;
  std::string output_name;
  if (stop)
    output_name = ".//tables//CSV//stop_table_fp_fried_pumps.csv";
  else
    output_name = ".//tables//CSV//table_fp_fried_pumps.csv";
  output.open(output_name.c_str(), std::fstream::out);

  if (!output.is_open())
  {
    std::cout << "Could not open file " << output_name << std::endl;
    throw 1;
  }

  // std::cout << output_name << std::endl;

  output << std::setprecision(2) << std::fixed;

  double curr_num_iter = 0.0, curr_time = 0.0, curr_profits_sum = 0.0, curr_num_perturb = 0.0, curr_num_stalls = 0.0;

  // output << "algo;";
  for (size_t j = 0; j < algorithms.size(); j++)
  {
    output << algorithms[j];
    if (j < algorithms.size() - 1)
      output << ",";
  }
  output << std::endl;

  // output << "instance;";
  for (size_t j = 0; j < algorithms.size(); j++)
  {
    // output << "max_iter-minus-num_iter;";
  }
  // output << std::endl;

  for (size_t j = 0; j < dirs.size(); j++)
  {
    std::vector<std::string> instances;
    std::string folder = dirs[j].substr(20);
    AddInstancesFromDirectory(dirs[j], instances, false);

    for (size_t i = 0; i < instances.size(); i++)
    {
      // output << instances[i] << ";";
      double original_lp = 0.0;
      for (size_t j = 0; j < algorithms.size(); j++)
      {
        curr_num_iter = curr_time = curr_profits_sum = curr_num_perturb = curr_num_stalls = 0.0;
        for (size_t seed = 1; seed <= 10; seed++)
        {
          curr_file = ".//solutions//";
          curr_file.append(folder);
          curr_file.append("s_");
          if (stop)
            curr_file.append("stop_");
          curr_file.append(algorithms[j]);
          curr_file.append("_seed_");
          curr_file += (std::to_string(seed) + "_");
          curr_file.append(instances[i]);

          // std::cout << curr_file << std::endl;
          // getchar();getchar();
          std::fstream input;
          input.open(curr_file.c_str(), std::fstream::in);

          if (!input.is_open())
          {
            std::cout << "Could not open file " << curr_file << std::endl;
            throw 4;
          }

          std::stringstream s_iter1, s_iter2, s_t1, s_t2, s_profits_sum, s_num_perturb1, s_num_perturb2, s_num_stalls1, s_num_stalls2;
          std::string status;
          double iter1 = 0.0, iter2 = 0.0, t1 = 0.0, t2 = 0.0, profits_sum = 0.0, num_perturb1 = 0.0, num_perturb2 = 0.0, num_stalls1 = 0.0, num_stalls2 = 0.0;
          std::string line;

          getline(input, line);
          size_t pos = line.find_first_of(":");
          status = line.substr(pos + 2);

          std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
          status.erase(end_pos, status.end());

          // std::cout << status << std::endl;

          getline(input, line);
          pos = line.find_first_of(":");
          s_profits_sum << line.substr(pos + 2);
          s_profits_sum >> profits_sum;

          // std::cout << "sum: " << profits_sum << std::endl;

          getline(input, line);
          getline(input, line);
          pos = line.find_first_of(":");
          s_iter1 << line.substr(pos + 2);
          s_iter1 >> iter1;

          // std::cout << "iter1: " << iter1 << std::endl;

          getline(input, line);
          pos = line.find_first_of(":");
          s_num_perturb1 << line.substr(pos + 2);
          s_num_perturb1 >> num_perturb1;

          getline(input, line);
          pos = line.find_first_of(":");
          s_num_stalls1 << line.substr(pos + 2);
          s_num_stalls1 >> num_stalls1;

          getline(input, line);
          pos = line.find_first_of(":");
          s_t1 << line.substr(pos + 2);
          s_t1 >> t1;

          // std::cout << "t1: " << t1 << std::endl;

          getline(input, line);
          getline(input, line);
          pos = line.find_first_of(":");
          s_iter2 << line.substr(pos + 2);
          s_iter2 >> iter2;

          // std::cout << "iter2: " << iter2 << std::endl;

          getline(input, line);

          pos = line.find_first_of(":");
          s_num_perturb2 << line.substr(pos + 2);
          s_num_perturb2 >> num_perturb2;

          getline(input, line);
          pos = line.find_first_of(":");
          s_num_stalls2 << line.substr(pos + 2);
          s_num_stalls2 >> num_stalls2;

          getline(input, line);
          pos = line.find_first_of(":");
          s_t2 << line.substr(pos + 2);
          s_t2 >> t2;

          // std::cout << "t2: " << t2 << std::endl;
          input.close();

          curr_profits_sum += profits_sum;
          curr_num_iter += (iter2);
          curr_time += (0 + t2);
          curr_num_perturb += (num_perturb1 + num_perturb2);
          curr_num_stalls += (num_stalls1 + num_stalls2);

          // getchar(); getchar();
        }
        output << (1.0 * max_iter_stage2 - curr_num_iter / num_executions);
        if (j < algorithms.size() - 1)
          output << ",";
        // output << curr_profits_sum/num_executions << ";" << curr_num_perturb/num_executions << ";" << curr_num_stalls/num_executions << ";" << curr_num_iter/num_executions << ";" << curr_time/num_executions << ";";
      }
      output << std::endl;
    }
  }

  output.close();
}

void GenerateCSVTableFeasibilityPump(std::vector<std::string> dirs, bool stop, int num_executions = 10)
{
  std::string curr_file;
  std::vector<std::string> algorithms;

  /*algorithms.push_back("fp_alpha2_0_p_25");
  algorithms.push_back("fp_alpha2_100_p_25");

  algorithms.push_back("fp_alpha2_0_p_15");
  algorithms.push_back("fp_alpha2_100_p_15");

  algorithms.push_back("fp_alpha2_0_p_10");
  algorithms.push_back("fp_alpha2_100_p_10");*/

  algorithms.push_back("fp_alpha2_0_p_fixed_10");
  algorithms.push_back("fp_cuts_alpha2_0_p_fixed_10");

  algorithms.push_back("fp_alpha2_100_p_fixed_10");
  algorithms.push_back("fp_cuts_alpha2_100_p_fixed_10");

  /*algorithms.push_back("fp_alpha1_0_alpha2_0_p_25");
  algorithms.push_back("fp_alpha1_100_alpha2_0_p_25");
  algorithms.push_back("fp_alpha1_0_alpha2_100_p_25");
  algorithms.push_back("fp_alpha1_100_alpha2_100_p_25");*/

  std::fstream output;
  std::string output_name;
  if (stop)
    output_name = ".//tables//CSV//stop_table_fp.csv";
  else
    output_name = ".//tables//CSV//table_fp.csv";
  output.open(output_name.c_str(), std::fstream::out);

  if (!output.is_open())
  {
    std::cout << "Could not open file " << output_name << std::endl;
    throw 1;
  }

  // std::cout << output_name << std::endl;

  output << std::setprecision(2) << std::fixed;

  double curr_num_iter = 0.0, curr_time = 0.0, curr_profits_sum = 0.0, curr_num_perturb = 0.0, curr_num_stalls = 0.0;

  output << "algo;";
  for (size_t j = 0; j < algorithms.size(); j++)
  {
    for (int i = 1; i <= 5; i++)
      output << algorithms[j] << ";";
  }
  output << std::endl;

  output << "instance;";
  for (size_t j = 0; j < algorithms.size(); j++)
  {
    output << "profits_sum;num_perturb;num_restarts;num_iter;time;";
  }
  output << std::endl;

  for (size_t j = 0; j < dirs.size(); j++)
  {
    std::vector<std::string> instances;
    std::string folder = dirs[j].substr(20);
    AddInstancesFromDirectory(dirs[j], instances, false);

    for (size_t i = 0; i < instances.size(); i++)
    {
      output << instances[i] << ";";
      double original_lp = 0.0;
      for (size_t j = 0; j < algorithms.size(); j++)
      {
        curr_num_iter = curr_time = curr_profits_sum = curr_num_perturb = curr_num_stalls = 0.0;
        for (size_t seed = 1; seed <= 10; seed++)
        {
          curr_file = ".//solutions//";
          curr_file.append(folder);
          curr_file.append("s_");
          if (stop)
            curr_file.append("stop_");
          curr_file.append(algorithms[j]);
          curr_file.append("_seed_");
          curr_file += (std::to_string(seed) + "_");
          curr_file.append(instances[i]);

          std::cout << curr_file << std::endl;
          // getchar();getchar();
          std::fstream input;
          input.open(curr_file.c_str(), std::fstream::in);

          if (!input.is_open())
          {
            std::cout << "Could not open file " << curr_file << std::endl;
            throw 4;
          }

          std::stringstream s_iter1, s_iter2, s_t1, s_t2, s_profits_sum, s_num_perturb1, s_num_perturb2, s_num_stalls1, s_num_stalls2;
          std::string status;
          double iter1 = 0.0, iter2 = 0.0, t1 = 0.0, t2 = 0.0, profits_sum = 0.0, num_perturb1 = 0.0, num_perturb2 = 0.0, num_stalls1 = 0.0, num_stalls2 = 0.0;
          std::string line;

          getline(input, line);
          size_t pos = line.find_first_of(":");
          status = line.substr(pos + 2);

          std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
          status.erase(end_pos, status.end());

          // std::cout << status << std::endl;

          getline(input, line);
          pos = line.find_first_of(":");
          s_profits_sum << line.substr(pos + 2);
          s_profits_sum >> profits_sum;

          // std::cout << "sum: " << profits_sum << std::endl;

          getline(input, line);
          getline(input, line);
          pos = line.find_first_of(":");
          s_iter1 << line.substr(pos + 2);
          s_iter1 >> iter1;

          // std::cout << "iter1: " << iter1 << std::endl;

          getline(input, line);
          pos = line.find_first_of(":");
          s_num_perturb1 << line.substr(pos + 2);
          s_num_perturb1 >> num_perturb1;

          getline(input, line);
          pos = line.find_first_of(":");
          s_num_stalls1 << line.substr(pos + 2);
          s_num_stalls1 >> num_stalls1;

          getline(input, line);
          pos = line.find_first_of(":");
          s_t1 << line.substr(pos + 2);
          s_t1 >> t1;

          // std::cout << "t1: " << t1 << std::endl;

          getline(input, line);
          getline(input, line);
          pos = line.find_first_of(":");
          s_iter2 << line.substr(pos + 2);
          s_iter2 >> iter2;

          // std::cout << "iter2: " << iter2 << std::endl;

          getline(input, line);

          pos = line.find_first_of(":");
          s_num_perturb2 << line.substr(pos + 2);
          s_num_perturb2 >> num_perturb2;

          getline(input, line);
          pos = line.find_first_of(":");
          s_num_stalls2 << line.substr(pos + 2);
          s_num_stalls2 >> num_stalls2;

          getline(input, line);
          pos = line.find_first_of(":");
          s_t2 << line.substr(pos + 2);
          s_t2 >> t2;

          // std::cout << "t2: " << t2 << std::endl;
          input.close();

          curr_profits_sum += profits_sum;
          curr_num_iter += (iter1 + iter2);
          curr_time += (0 + t2);
          curr_num_perturb += (num_perturb1 + num_perturb2);
          curr_num_stalls += (num_stalls1 + num_stalls2);

          // getchar(); getchar();
        }
        output << curr_profits_sum / num_executions << ";" << curr_num_perturb / num_executions << ";" << curr_num_stalls / num_executions << ";" << curr_num_iter / num_executions << ";" << curr_time / num_executions << ";";
      }
      output << std::endl;
    }
  }

  output.close();
}

void SolveFeasibilityPump(std::vector<std::string> &instances)
{
  double time_limit = -1;
  int option = -1;
  Graph *graph = NULL;
  Solution<int> *sol = NULL;
  double *r = NULL, *R = NULL, *Rn = NULL;
  int seed = 1;

  std::string folder;
  std::string file_name;

  std::cout << "COMPACT FORMULATIONS:" << std::endl;

  do
  {
    std::cout << "*******************************************" << std::endl
              << " * Time limit in seconds (-1 to be unlimited): " << std::endl;
    std::cin >> time_limit;
  } while (double_less(time_limit, 0.0) && !(double_equals(time_limit, -1)));

  do
  {
    std::cout << "*******************************************" << std::endl
              << " 0 - Run ALL" << std::endl
              << " 1 - Run Butt&Ryan single commodity" << std::endl
              << " 2 - Run Butt&Ryan multi commodity" << std::endl
              << " 3 - Run Butt&Ryan two commodity" << std::endl
              << " 4 - Run Butt&Ryan multi two commodity" << std::endl
              << " 5 - Run Capacitated single commodity" << std::endl
              << " 6 - Run Capacitated multi commodity" << std::endl
              << " 7 - Run Capacitated two commodity" << std::endl
              << " 8 - Run Capacitated multi two commodity" << std::endl
              << " 9 - Run Kulkarni and Bhave" << std::endl
              << " 10 - Run Bianchessi et al." << std::endl
              << " 11 - Run Capacitated single commodity Extended" << std::endl
              << "Option: ";
    std::cin >> option;
    switch (option)
    {
    case 0:
    {
      for (size_t i = 0; i < instances.size(); i++)
      {
        try
        {
          int profit_variation = 0, profit_variation2 = 0;
          double time_variation = 0, time_variation2 = 0, global_variation = 0;
          std::list<int>::iterator it;

          split_file_path(instances[i], folder, file_name);
          std::cout << "* " << file_name << std::endl;
          Instance inst(instances[i]);
          // std::cout << inst.graph()->num_mandatory() << std::endl;
          /*ALNS y;
          y.Init(inst,"3opt",folder,file_name);

          std::list<int>::iterator vertex_it;

          if( (y.best_solution())->CheckCorrectness(inst)) std::cout << *(y.best_solution()) << std::endl;
          if(y.Do_1_0_Improvement((y.best_solution()),0,2)) std::cout << "moveu" << std::endl;
          else std::cout << "nao moveu" << std::endl;
          if( (y.best_solution())->CheckCorrectness(inst)) std::cout << *(y.best_solution()) << std::endl;
          if(y.Do_1_0_Improvement((y.best_solution()),1,2)) std::cout << "moveu" << std::endl;
          else std::cout << "nao moveu" << std::endl;
          if( (y.best_solution())->CheckCorrectness(inst)) std::cout << *(y.best_solution()) << std::endl;*/
          /*if((y.best_solution())->PreviewAddVertexToRouteWithinMinimumDistanceIncrease(inst,3,0,vertex_it,profit_variation,time_variation))
          {
            std::cout << "pode add before ";
            if(vertex_it != (y.best_solution())->routes_vec_[0].vertices_.end()) std::cout << *vertex_it;
            else std::cout << inst.graph()->num_vertices() -1;

            std::cout << " with time variation of " << time_variation << std::endl;
          }else std::cout << "nao pode add" << std::endl;*/

          /*std::cout << " * Before remove" << std::endl;
          std::cout << *(y.best_solution()) << std::endl;
          getchar(); getchar();
          y.RemoveVerticesFromSolution(y.best_solution(), 0.75);

          std::cout << " * After remove" << std::endl;
          std::cout << *(y.best_solution()) << std::endl;
          getchar(); getchar();
          y.TryToInsertUnvisitedVertices(y.best_solution());

          std::cout << " * After add" << std::endl;
          std::cout << *(y.best_solution()) << std::endl;
          getchar(); getchar();*/

          /*HeuristicSolution * curr_sol = NULL;
          curr_sol = new HeuristicSolution(y.best_solution());
          (curr_sol->bitset_arcs_)[0] = 1;
          curr_sol->profits_sum_ = 2;
          y.AddSolutionToPool(curr_sol,0);

          curr_sol = new HeuristicSolution(y.best_solution());
          (curr_sol->bitset_arcs_)[1] = 1;
          curr_sol->profits_sum_ = 1;

          y.AddSolutionToPool(curr_sol,0);

          curr_sol = new HeuristicSolution(y.best_solution());
          (curr_sol->bitset_arcs_)[2] = 1;
          curr_sol->profits_sum_ = 6;

          y.AddSolutionToPool(curr_sol,0);

          curr_sol = new HeuristicSolution(y.best_solution());
          (curr_sol->bitset_arcs_)[3] = 1;
          curr_sol->profits_sum_ = 4;

          y.AddSolutionToPool(curr_sol,0);

          y.PrintPool();

          curr_sol = new HeuristicSolution(y.best_solution());
          (curr_sol->bitset_arcs_)[4] = 1;
          curr_sol->profits_sum_ = 10;

          y.AddSolutionToPool(curr_sol,0);

          curr_sol = new HeuristicSolution(y.best_solution());
          (curr_sol->bitset_arcs_)[5] = 1;
          curr_sol->profits_sum_ = 2;

          y.AddSolutionToPool(curr_sol,0);

          y.PrintPool();

          getchar();getchar();
          continue;*/

          /*FeasibilityPump x;
          x.Init(inst);

          (x.solution_).ReadFromFile(inst,std::string("3opt"), folder, file_name);
          (x.solution_).BuildBitset(inst);
          std::cout << x.solution_.bitset_arcs_ << std::endl;
          std::cout << x.solution_ << std::endl;
          getchar();getchar();

          if(!((x.solution_).CheckCorrectness(inst))){ std::cout << "deu merda" << std::endl; getchar();getchar();}

          std::cout << x.solution_ << std::endl;

          std::cout << "* preview add to route 0" << std::endl;
          if((x.solution_).PreviewAddVertexToRouteWithinMinimumDistanceIncrease(inst.graph(),5,0,it,profit_variation,time_variation))
          {
            std::cout << "profit var: " << profit_variation << " time var: " << time_variation << std::endl;
            if(it != ((((x.solution_).routes_vec_)[0]).vertices_).end()) std::cout << "route 0 before " << *it << std::endl;
            else std::cout << "route 0 before " << inst.graph()->num_vertices() - 1 << std::endl;
            //(x.solution_).AddVertex(5,0,it,profit_variation,time_variation);
          }
          else std::cout << "cannot add" << std::endl;
          std::cout << x.solution_ << std::endl;

          std::cout << "* preview add to route 0" << std::endl;
          int curr_route = 0;
          if((x.solution_).PreviewAddVertexWithinMinimumDistanceIncrease(inst.graph(),5,curr_route,it,profit_variation,time_variation))
          {
            std::cout << "profit var: " << profit_variation << " time var: " << time_variation << std::endl;
            if(it != ((((x.solution_).routes_vec_)[curr_route]).vertices_).end()) std::cout << "route " << curr_route << " before " << *it << std::endl;
            else std::cout << "route " << curr_route << " before " << inst.graph()->num_vertices() - 1 << std::endl;
            (x.solution_).AddVertex(5,curr_route,it,profit_variation,time_variation);
          }
          else std::cout << "cannot add" << std::endl;
          std::cout << x.solution_ << std::endl;

          getchar(); getchar();*/

          /*(x.solution_).Do3OptImprovement(inst.graph(),0);

          std::cout << x.solution_ << std::endl;
          if(!((x.solution_).CheckCorrectness(inst))){ std::cout << "deu merda" << std::endl; getchar();getchar();}

          getchar(); getchar();*/

          /*std::cout << "* preview add 1" << std::endl;
          (x.solution_).PreviewAddVertex(inst.graph(),1,1,-1, profit_variation, time_variation);
          std::cout << "profit var: " << profit_variation << " time var: " << time_variation << std::endl;
          getchar(); getchar();

          std::cout << "* preview add 2" << std::endl;
          (x.solution_).PreviewAddVertex(inst.graph(),2,1,-1, profit_variation, time_variation);
          std::cout << "profit var: " << profit_variation << " time var: " << time_variation << std::endl;
          getchar(); getchar();*/

          /*if((x.solution_).PreviewAddVertex(inst.graph(),2,1,-1,it,profit_variation,time_variation)) (x.solution_).AddVertex(2,1,it,profit_variation,time_variation);
          std::cout << x.solution_;
          if(!((x.solution_).CheckCorrectness(inst))){ std::cout << "deu merda" << std::endl; getchar();getchar();}
          //getchar(); getchar();

          if((x.solution_).PreviewAddVertex(inst.graph(),3,1,-1,it,profit_variation,time_variation)) (x.solution_).AddVertex(3,1,it,profit_variation,time_variation);
          std::cout << x.solution_;
          //getchar(); getchar();
          if(!((x.solution_).CheckCorrectness(inst))){ std::cout << "deu merda" << std::endl; getchar();getchar();}

          if((x.solution_).PreviewAddVertex(inst.graph(),1,1,0,it,profit_variation,time_variation)) (x.solution_).AddVertex(1,1,it,profit_variation,time_variation);
          std::cout << x.solution_;
          if(!((x.solution_).CheckCorrectness(inst))){ std::cout << "deu merda" << std::endl; getchar();getchar();}*/
          // getchar(); getchar();
          // getchar(); getchar();
          // getchar(); getchar();
          //(x.solution_).RemoveVertex(2,0,0);
          // std::cout << x.solution_;
          // getchar(); getchar();
          //(x.solution_).AddVertex(2,1,1,0,0);
          // std::cout << x.solution_;
          // getchar(); getchar();

          //(x.solution_).PreviewRemoveVertex(inst.graph(), int vertex, profit_variation, time_variation);
          /*std::cout << "* preview remove 2" << std::endl;
          (x.solution_).PreviewRemoveVertex(inst.graph(),2, profit_variation, time_variation);
          std::cout << "profit var: " << profit_variation << " time var: " << time_variation << std::endl;
          getchar(); getchar();
          std::cout << "* preview remove 1" << std::endl;
          (x.solution_).PreviewRemoveVertex(inst.graph(),1, profit_variation, time_variation);
          std::cout << "profit var: " << profit_variation << " time var: " << time_variation << std::endl;
          getchar(); getchar();
          std::cout << "* preview remove 3" << std::endl;
          (x.solution_).PreviewRemoveVertex(inst.graph(),3, profit_variation, time_variation);
          std::cout << "profit var: " << profit_variation << " time var: " << time_variation << std::endl;
          getchar(); getchar();*/
          //(x.solution_).PreviewAddVertex(inst.graph(), int vertex, int route, int pos, int& profit_variation,time_variation);

          // bool AddVertex(int vertex, int route, int pos, double profit_variation, double time_variation);
          // bool RemoveVertex(int vertex, double profit_variation, double time_variation);
          // MoveVertex(int vertex, int r1, int r2, int pos, double profit_variation1, double time_variation1, double profit_variation2, double time_variation2);

          /*std::cout << "* preview move 2 to route 0 pos -1" << std::endl;
          bool can_do = (x.solution_).PreviewInterRouteMoveVertex(inst.graph(),2,0,-1,it,profit_variation,time_variation,profit_variation2,time_variation2);
          std::cout << "profit var: " << profit_variation << " time var: " << time_variation << std::endl;
          std::cout << "profit var2: " << profit_variation2 << " time var2: " << time_variation2 << std::endl;
          if(can_do) (x.solution_).InterRouteMoveVertex(2,0,it,profit_variation,time_variation,profit_variation2,time_variation2);
          std::cout << x.solution_;
          if(!((x.solution_).CheckCorrectness(inst))){ std::cout << "deu merda" << std::endl; getchar();getchar();}
          //getchar(); getchar();

          std::cout << "* preview move 1 to route 0 pos 0" << std::endl;
          can_do = (x.solution_).PreviewInterRouteMoveVertex(inst.graph(),1,0,0,it,profit_variation,time_variation,profit_variation2,time_variation2);
          std::cout << "profit var: " << profit_variation << " time var: " << time_variation << std::endl;
          std::cout << "profit var2: " << profit_variation2 << " time var2: " << time_variation2 << std::endl;
          if(can_do) (x.solution_).InterRouteMoveVertex(1,0,it,profit_variation,time_variation,profit_variation2,time_variation2);
          std::cout << x.solution_;
          if(!((x.solution_).CheckCorrectness(inst))){ std::cout << "deu merda" << std::endl; getchar();getchar();}
          //getchar(); getchar();

          std::list<int>::iterator it_i1, it_f1, it_i2, it_f2;
          std::cout << "* preview swap route 0 [0,0] and route 1 [0,0]" << std::endl;
          can_do = (x.solution_).PreviewInterRouteSwap(inst.graph(),0,0,0,1,0,0,it_i1,it_f1,profit_variation,time_variation,it_i2,it_f2,profit_variation2,time_variation2);
          std::cout << "profit var: " << profit_variation << " time var: " << time_variation << std::endl;
          std::cout << "profit var2: " << profit_variation2 << " time var2: " << time_variation2 << std::endl;
          if(!can_do) std::cout << "INVALIDO" << std::endl;
          if(can_do) (x.solution_).InterRouteSwap(0,1,it_i1,it_f1,profit_variation,time_variation,it_i2,it_f2,profit_variation2,time_variation2);
          std::cout << x.solution_;
          if(!((x.solution_).CheckCorrectness(inst))){ std::cout << "deu merda" << std::endl; getchar();getchar();}

          std::cout << "* preview swap route 0 [1,1] and route 1 [0,0]" << std::endl;
          can_do = (x.solution_).PreviewInterRouteSwap(inst.graph(),0,1,1,1,0,0,it_i1,it_f1,profit_variation,time_variation,it_i2,it_f2,profit_variation2,time_variation2);
          std::cout << "profit var: " << profit_variation << " time var: " << time_variation << std::endl;
          std::cout << "profit var2: " << profit_variation2 << " time var2: " << time_variation2 << std::endl;
          if(!can_do) std::cout << "INVALIDO" << std::endl;
          if(can_do) (x.solution_).InterRouteSwap(0,1,it_i1,it_f1,profit_variation,time_variation,it_i2,it_f2,profit_variation2,time_variation2);
          std::cout << x.solution_;
          if(!((x.solution_).CheckCorrectness(inst))){ std::cout << "deu merda" << std::endl; getchar();getchar();}*/
          // getchar();getchar();
          // if(can_do) (x.solution_).InterRouteMoveVertex(1,0,it,profit_variation,time_variation,profit_variation2,time_variation2);
          // std::cout << x.solution_;
          /*std::cout << "* preview swap route 0 [0,1] and route 1 [0,0]" << std::endl;
          can_do = (x.solution_).PreviewInterRouteSwap(inst.graph(),0,0,1,1,0,0,it_i1,it_f1,profit_variation,time_variation,it_i2,it_f2,profit_variation2,time_variation2);
          std::cout << "profit var: " << profit_variation << " time var: " << time_variation << std::endl;
          std::cout << "profit var2: " << profit_variation2 << " time var2: " << time_variation2 << std::endl;
          if(!can_do) std::cout << "INVALIDO" << std::endl;
          if(can_do) (x.solution_).InterRouteSwap(0,1,it_i1,it_f1,profit_variation,time_variation,it_i2,it_f2,profit_variation2,time_variation2);
          std::cout << x.solution_;*/

          /*for(int i = 1; i < (inst.graph())->num_vertices() -1; i++)
          {
            std::cout << i << ":" << std::endl;
            std::cout << "route: " << (((x.solution_).vertex_status_vec_)[i]).route_ << std::endl;
            std::cout << "pos: " << *((((x.solution_).vertex_status_vec_)[i]).pos_) << std::endl;
          }

          getchar();getchar();*/

          /*std::cout << "* preview remove 1" << std::endl;
          if((x.solution_).PreviewRemoveVertex(inst.graph(),1, profit_variation, time_variation))
            (x.solution_).RemoveVertex(1, profit_variation, time_variation);

          std::cout << "profit var: " << profit_variation << " time var: " << time_variation << std::endl;
          std::cout << x.solution_;
          if(!((x.solution_).CheckCorrectness(inst))){ std::cout << "deu merda" << std::endl; getchar();getchar();}
          //getchar(); getchar();

          std::cout << "* preview remove 2" << std::endl;
          if((x.solution_).PreviewRemoveVertex(inst.graph(),2, profit_variation, time_variation))
            (x.solution_).RemoveVertex(2, profit_variation, time_variation);

          std::cout << "profit var: " << profit_variation << " time var: " << time_variation << std::endl;
          std::cout << x.solution_;
          if(!((x.solution_).CheckCorrectness(inst))){ std::cout << "deu merda" << std::endl; getchar();getchar();}
          //getchar(); getchar();

          std::cout << "* preview remove 3" << std::endl;
          if((x.solution_).PreviewRemoveVertex(inst.graph(),3, profit_variation, time_variation))
            (x.solution_).RemoveVertex(3, profit_variation, time_variation);

          std::cout << "profit var: " << profit_variation << " time var: " << time_variation << std::endl;
          std::cout << x.solution_;
          if(!((x.solution_).CheckCorrectness(inst))){ std::cout << "deu merda" << std::endl; getchar();getchar();}
          //getchar(); getchar();


          x.Run();

          (x.solution_).WriteToFile((x.solution_).GenerateFileName() + "_seed_" + std::to_string(seed), folder, file_name);

          (x.solution_).ReadFromFile(inst,(x.solution_).GenerateFileName() + "_seed_" + std::to_string(seed), folder, file_name);

          if(!((x.solution_).CheckCorrectness(inst))){ std::cout << "deu merda" << std::endl; getchar();getchar();}
          x.Reset();

          continue;*/

          graph = inst.graph();

          sol = new Solution<int>(graph->num_vertices());

          r = Dijkstra(graph, false, true);
          R = Dijkstra(graph, false, false);
          Rn = Dijkstra(graph, true, false);

          sol->write_to_file("csc3", folder, file_name);

          delete[] R;
          R = NULL;

          delete[] Rn;
          Rn = NULL;

          delete[] r;
          r = NULL;

          delete sol;
          sol = NULL;
        }
        catch (const std::bad_alloc &ex)
        {
          sol->out_of_memory_ = true;
          delete[] R;
          R = NULL;

          delete[] Rn;
          Rn = NULL;

          delete[] r;
          r = NULL;

          delete sol;
          sol = NULL;
          std::cout << "Out of memory exception" << std::endl;
          continue;
        }
      }
      break;
    }
    default:
      option = -1;
      std::cout << "Invalid option!" << std::endl;
      break;
    }
  } while (option == -1);
}

void AddInstances(std::vector<std::string> &instances)
{
  // AddInstancesFromDirectory(".//instances//STOP//unsolved//",instances);
  AddInstancesFromDirectory(".//instances//STOP//current//", instances);
  /*AddInstancesFromDirectory(".//instances//STOP//Set_21_234//Set_0.05//",instances);
  AddInstancesFromDirectory(".//instances//STOP//Set_32_234//Set_0.05//",instances);
  AddInstancesFromDirectory(".//instances//STOP//Set_33_234//Set_0.05//",instances);
  AddInstancesFromDirectory(".//instances//STOP//Set_64_234//Set_0.05//",instances);
  AddInstancesFromDirectory(".//instances//STOP//Set_66_234//Set_0.05//",instances);
  AddInstancesFromDirectory(".//instances//STOP//Set_100_234//Set_0.05//",instances);
  AddInstancesFromDirectory(".//instances//STOP//Set_102_234//Set_0.05//",instances);*/
}

static const struct option longOpts[] = {
    {"feasibility-pump", no_argument, NULL, '0'},
    {"instance", required_argument, NULL, '1'},
    {"seed", required_argument, NULL, '2'},
    {"LNS", no_argument, NULL, '3'},
    {"simulated-annealing", no_argument, NULL, '4'},
    {NULL, no_argument, NULL, 0}};

void ParseArgumentsAndRun(int argc, char *argv[])
{
  std::string instance, folder, file_name;
  int c = 0, seed = 0;
  bool solve_feasibility_pump = false;
  bool solve_ALNS = false;
  bool solve_simulated_annealing = false;

  while ((c = getopt_long(argc, argv, "01:2:3", longOpts, NULL)) != -1)
  {
    switch (c)
    {
    case '0':
      solve_feasibility_pump = true;
      break;
    case '1':
      if (optarg)
        instance = std::string(optarg);
      break;
    case '2':
      if (optarg)
        seed = std::atoi(optarg);
      break;
    case '3':
      solve_ALNS = true;
      break;
    case '4':
      solve_simulated_annealing = true;
      break;
    }
  }

  srand(seed);
  split_file_path(instance, folder, file_name);
  std::cout << "* " << file_name << std::endl;

  Instance inst(instance);
  int num_vertices = (inst.graph())->num_vertices();

  // add to graph a slack arc to allow empty routes!
  GArc *slack = (*(inst.graph()))[0][num_vertices - 1];
  if (slack != NULL)
  {
    slack->set_dist(0.0);
  }
  else
    (inst.graph())->AddArc(0, num_vertices - 1, 0.0);

  if (solve_feasibility_pump)
  {
    FeasibilityPump fp;
    fp.Init(inst);
    std::cout << (fp.solution_).GenerateFileName() + "_seed_" + std::to_string(seed) << std::endl;
    fp.Run();
    (fp.solution_).WriteToFile((fp.solution_).GenerateFileName() + "_seed_" + std::to_string(seed), folder, file_name);
  }
  // std::cout << (fp.solution_) << std::endl;
  // if(!((fp.solution_).CheckCorrectness(inst))){ std::cout << "deu merda" << std::endl; getchar();getchar();}
  // fp.Reset();

  if (solve_ALNS)
  {
    // std::cout << ALNSHeuristicSolution::GenerateFileName() << "_seed_" << seed << std::endl;
    ALNS alns;
    // std::cout << "antes" << std::endl;
    alns.Init(inst, FPHeuristicSolution::GenerateFileName() + "_seed_" + std::to_string(seed), folder, file_name);
    // std::cout << "depois" << std::endl;
    // if(!((alns.best_solution())->CheckCorrectness(inst))){ std::cout << "deu merda" << std::endl; getchar();getchar();}
    // std::cout << "passou" << std::endl;

    alns.Run();
    std::cout << ALNSHeuristicSolution::GenerateFileName() << "_seed_" << seed << std::endl;

    (alns.best_solution())->WriteToFile(ALNSHeuristicSolution::GenerateFileName() + "_seed_" + std::to_string(seed), folder, file_name);
  }

  if (solve_local_branching)
  {
    LocalBranching lb;

    lb.Init(inst);

    HeuristicSolution initial_sol;
    std::cout << "K = " << K_LOCAL_BRANCHING_K << std::endl;

    if (K_LOCAL_BRANCHING_SOLVE_ALNS)
      initial_sol.ReadFromFile(inst, ALNSHeuristicSolution::GenerateFileName() + "_seed_" + std::to_string(seed), folder, file_name);
    else
      initial_sol.ReadFromFile(inst, FPHeuristicSolution::GenerateFileName() + "_seed_" + std::to_string(seed), folder, file_name);

    // std::cout << initial_sol << std::endl;
    // getchar(); getchar();
    lb.Run(180, &initial_sol);

    //(alns.best_solution())->WriteToFile(ALNSHeuristicSolution::GenerateFileName() + "_seed_" + std::to_string(seed), folder, file_name);
  }
}

int main(int argc, char *argv[])
{
  /*std::list<std::list<int>> conflicts_list;
  std::vector<std::list<int>> conflict_graph(3,std::list<int>());
  conflict_graph[0].push_back(1);
  conflict_graph[1].push_back(0);

    FindAllMaximalConflictCliques2(&conflict_graph,conflicts_list);
  return 0;*/

  double time_limit = 7200;
  try
  {

    if (K_GETOPT)
    {
      ParseArgumentsAndRun(argc, argv);
      DeleteTimer();
      DeleteCallbackSelection();
      return 0;
    }
    std::vector<std::string> instances;

    int option = -1;
    double mandatory_percentage = 0.05;

    std::cout << "#INSTANCES SELECTED HARD CODED." << std::endl;
    do
    {
      std::cout << "*******************************************" << std::endl
                << " 0 - Generate result tables" << std::endl
                << "*******************************************" << std::endl
                << " 1 - Solve Feasibility Pump" << std::endl
                << "Option: ";
      std::cin >> option;
      switch (option)
      {
      case 0:
      {
        std::vector<std::string> dirs;
        // dirs.push_back(".//instances//STOP//current//");
        // dirs.push_back(".//instances//STOP//Set_32_234//Set_0.1//");
        // dirs.push_back(".//instances//STOP//Set_21_234//Set_0.1//");
        // dirs.push_back(".//instances//STOP//Set_33_234//Set_0.1//");
        // dirs.push_back(".//instances//STOP//Set_100_234//Set_0.1//");
        // dirs.push_back(".//instances//STOP//Set_66_234//Set_0.1//");
        // dirs.push_back(".//instances//STOP//Set_64_234//Set_0.1//");
        // dirs.push_back(".//instances//STOP//Set_102_234//Set_0.1//");

        dirs.clear();
        dirs.push_back(".//instances//STOP//Set_32_234//Set_0.05//");
        dirs.push_back(".//instances//STOP//Set_21_234//Set_0.05//");
        dirs.push_back(".//instances//STOP//Set_33_234//Set_0.05//");
        dirs.push_back(".//instances//STOP//Set_100_234//Set_0.05//");
        dirs.push_back(".//instances//STOP//Set_66_234//Set_0.05//");
        dirs.push_back(".//instances//STOP//Set_64_234//Set_0.05//");
        dirs.push_back(".//instances//STOP//Set_102_234//Set_0.05//");
        // dirs.push_back(".//instances//STOP//unsolved//");

        // GenerateHeuristicSolutionImprovementLatexTable(dirs, false);
        // GeneratePrimalAndDualGapsLatexTable(dirs,false);
        // GenerateFeasibilityPumpLatexTable(dirs,true,10);
        // GenerateFeasibilityPumpPerInstanceResultsCSV(dirs,true,10);
        // GenerateFeasibilityPumpPerInstanceResultsCSV(dirs,false,10);
        GenerateCutsConfigurationsPerInstanceResultsCSV(dirs, true);
        GenerateCutsConfigurationsPerInstanceResultsCSV(dirs, false);
        // GenerateLNSPerInstanceResultsCSV(dirs,false,10);
        /*GenerateALNSLatexTable(dirs,true);
        GenerateALNSLatexTable2(dirs,true);*/

        // GenerateFeasibilityPumpLatexTable(dirs,false,10);
        // GenerateALNSLatexTable(dirs,true);
        // GenerateALNSLatexTable2(dirs,true);
        // GenerateALNSCSVTableFridmanGaps(dirs,true,10);
        // GenerateALNSCSVTableFridmanBestGaps(dirs,true,10);
        // GenerateALNSCSVTableFridmanIterations(dirs,true,10);
        // GenerateCSVTableFeasibilityPump(dirs,false,10);
        /*GenerateCSVTableFeasibilityPumpFriedmanIterations(dirs,true,10);
        GenerateCSVTableFeasibilityPumpFriedmanGaps(dirs,true,10);
        //GenerateCSVTableALNS(dirs,false,10);
        //GenerateALNSLatexAppendix(dirs);*/
        break;
      }
      case 1:
      {
        AddInstances(instances);
        SolveFeasibilityPump(instances);
        break;
      }
      default:
        option = -1;
        std::cout << "Invalid option!" << std::endl;
        break;
      }
    } while (option == -1);

    DeleteTimer();
    DeleteCallbackSelection();
  }
  catch (const std::runtime_error &re)
  {
    std::cout << "Runtime error: " << re.what() << std::endl;
  }
  catch (const std::exception &ex)
  {
    std::cout << "Error occurred: " << ex.what() << std::endl;
  }
  catch (const int &error)
  {
    std::cout << "Error occurred: " << error << std::endl;
  }
  catch (IloException &e)
  {
    std::cout << "Concert Exception: " << e << std::endl;
  }
  catch (const char *e)
  {
    std::cout << e << std::endl;
  }
  catch (const std::string &e)
  {
    std::cout << e << std::endl;
  }
  catch (...)
  {
    std::cout << "Unknown failure occurred. Possible memory corruption" << std::endl;
  }

  return 0;
}
