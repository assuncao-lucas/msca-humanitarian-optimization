#include <iostream>
#include <string>
#include <dirent.h>
#include <sstream>
#include <fstream>
#include "src/instance.h"
#include "src/graph.h"
#include "src/general.h"
#include "src/formulations.h"
#include "src/timer.h"
#include "src/graph_algorithms.h"

void GenerateAlgorithmsLatexTable(double total_time_limit)
{
	std::string curr_file;
	std::vector<std::string> algorithms;
	std::vector<std::string> instance_sizes;
	std::vector<std::string> instance_types;
	std::vector<std::string> instance_limit_quantiles;

	std::vector<std::string> num_vehicles_vec{"2","4"};
	std::vector<std::string> service_time_deviation_vec{"0.50"};
	std::vector<std::string> uncertainty_budget_vec{"1","5"};

	std::vector<std::string> instances_terminations;

	for(auto num_vehicles: num_vehicles_vec)
		for(auto service_time_deviation: service_time_deviation_vec)
			for (auto uncertainty_budget: uncertainty_budget_vec)
				instances_terminations.push_back("_v" + num_vehicles + "_d" + service_time_deviation + "_b" + uncertainty_budget + ".txt");

	algorithms.push_back("baseline");
	//algorithms.push_back("cb_baseline");
	//algorithms.push_back("csc");
	//algorithms.push_back("cb_csc");

	instance_types.push_back("C");
	instance_types.push_back("R");
	instance_types.push_back("RC");

	instance_sizes.push_back("25");
	instance_sizes.push_back("50");
	//instance_sizes.push_back("100");

	instance_limit_quantiles.push_back("0.8");
	instance_limit_quantiles.push_back("1");
	//instance_limit_quantiles.push_back("2");
	instance_limit_quantiles.push_back("3");
	//instance_limit_quantiles.push_back("4");

	std::fstream output;
	std::string output_name;
	output_name = "..//tables//latex//table_algorithms.txt";
	output.open(output_name.c_str(),std::fstream::out);

	if(!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	//std::cout << output_name << std::endl;

	output << std::setprecision(2) << std::fixed;

	std::vector<std::vector<double>> total_time_per_algo(algorithms.size(),std::vector<double>());
	std::vector<std::vector<double>> total_gap_per_algo(algorithms.size(),std::vector<double>());
	std::vector<double> total_avg_time(algorithms.size(),0.0);
	std::vector<double> total_avg_gap(algorithms.size(),0.0);
	std::vector<int> total_num_optimal(algorithms.size(),0);
	size_t num_inst_per_vertex_size_type_quantile = num_vehicles_vec.size()*service_time_deviation_vec.size()*uncertainty_budget_vec.size();
	size_t num_inst_per_vertex_size_quantile = instance_types.size() * num_inst_per_vertex_size_type_quantile;
	size_t num_inst_per_vertex_size = num_inst_per_vertex_size_quantile * instance_limit_quantiles.size();
	size_t total_num_instances = num_inst_per_vertex_size * instance_sizes.size();

	for(auto instance_size : instance_sizes)
	{
		std::vector<std::vector<double>> time_per_algo_inst_size(algorithms.size(),std::vector<double>());
		std::vector<std::vector<double>> gap_per_algo_inst_size(algorithms.size(),std::vector<double>());
		std::vector<double> avg_time_inst_size(algorithms.size(),0.0);
		std::vector<double> avg_gap_inst_size(algorithms.size(),0.0);
		std::vector<double> st_dev_inst_size(algorithms.size(),0.0);

		std::vector<int> num_optimal_inst_size(algorithms.size(),0);

		for(auto instance_type : instance_types)
		{
			for(auto instance_limit_quantile : instance_limit_quantiles)
			{
				std::vector<std::vector<double>> time_per_algo_quantile(algorithms.size(),std::vector<double>());
				std::vector<std::vector<double>> gap_per_algo_quantile(algorithms.size(),std::vector<double>());
				std::vector<double> avg_time_quantile(algorithms.size(),0.0);
				std::vector<double> avg_gap_quantile(algorithms.size(),0.0);
				std::vector<double> st_dev_quantile(algorithms.size(),0.0);

				std::vector<int> num_optimal_quantile(algorithms.size(),0);

				std::string instance_prefix = instance_type + instance_size + "_" + instance_limit_quantile;

				for(auto instances_termination: instances_terminations)
				{
					std::string instance = instance_prefix + instances_termination;
					//std::cout << instance << std::endl;
					double original_lp = 0.0;
					for(size_t algo = 0; algo < algorithms.size(); ++algo)
					{
						curr_file = "..//solutions//";
						curr_file.append("s_");
						curr_file.append(algorithms[algo]);
						curr_file.append("_");
						curr_file.append(instance);

						//std::cout << curr_file << std::endl;

						std::fstream input;
						input.open(curr_file.c_str(),std::fstream::in);

						if(!input.is_open())
						{
							std::cout << "Could not open file " << curr_file << std::endl;
							throw 4;
						}

						std::stringstream s_lb, s_ub, s_time;
						std::string status;
						double lb = 0.0, ub = 0.0, time = 0.0, gap = 0.0;
						std::string line;

						getline(input,line);
						size_t pos = line.find_first_of(":");
						status = line.substr(pos + 2);

						std::string::iterator end_pos = std::remove(status.begin(),status.end(),' ');
						status.erase(end_pos, status.end());

						getline(input,line);
						getline(input,line);
						pos = line.find_first_of(":");
						s_lb << line.substr(pos + 2);
						if(s_lb.str() == "-inf") lb = -1;
						else s_lb >> lb;

						getline(input,line);
						pos = line.find_first_of(":");
						s_ub << line.substr(pos + 2);
						if(s_ub.str() == "inf") ub = -1;
						else s_ub >> ub;

						getline(input,line);
						getline(input,line);
						pos = line.find_first_of(":");
						s_time << line.substr(pos + 2);
						s_time >> time;

						if((status == "OPTIMAL")||(status == "INFEASIBLE"))
							std::cout << instance << std::endl;

						if((status != "OPTIMAL") && (status != "INFEASIBLE"))
						{
							if((!(double_equals(lb,-1))) && (!(double_equals(ub,-1))))
							{
								if((!double_greater(ub,lb)))
								//if(double_less(time,total_time_limit))
								{
									++(num_optimal_inst_size[algo]);
									++(num_optimal_quantile[algo]);
									++(total_num_optimal[algo]);
									time_per_algo_inst_size[algo].push_back(time);
									time_per_algo_quantile[algo].push_back(time);
									total_time_per_algo[algo].push_back(time);

									total_avg_time[algo] += time;
									avg_time_inst_size[algo] += time;
									avg_time_quantile[algo] += time;
								}else gap = (100.0*(ub-lb))/ub;
							}else gap = 100.0;
						}else
						{
							++(num_optimal_inst_size[algo]);
							++(num_optimal_quantile[algo]);
							++(total_num_optimal[algo]);
							time_per_algo_inst_size[algo].push_back(time);
							time_per_algo_quantile[algo].push_back(time);
							total_time_per_algo[algo].push_back(time);

							total_avg_time[algo] += time;
							avg_time_inst_size[algo] += time;
							avg_time_quantile[algo] += time;
						}

						//std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

						if(!double_equals(gap,0.0))
						{
							gap_per_algo_inst_size[algo].push_back(gap);
							gap_per_algo_quantile[algo].push_back(gap);
							total_gap_per_algo[algo].push_back(gap);

							total_avg_gap[algo] += gap;
							avg_gap_quantile[algo] += gap;
							avg_gap_inst_size[algo] += gap;
						}
						input.close();
					}

					//getchar();getchar();
				}

				output << instance_type + instance_size + "Q" + instance_limit_quantile;

				for(size_t algo = 0; algo < algorithms.size(); ++algo)
				{
					if((time_per_algo_quantile[algo]).size() > 0) avg_time_quantile[algo]/=(1.0*((time_per_algo_quantile[algo]).size()));
					else avg_time_quantile[algo] = -1;
					if((gap_per_algo_quantile[algo]).size() > 0)
					{
						avg_gap_quantile[algo]/=(1.0*((gap_per_algo_quantile[algo]).size()));
						st_dev_quantile[algo] = StDev(gap_per_algo_quantile[algo],avg_gap_quantile[algo]);
					}else avg_gap_quantile[algo] = st_dev_quantile[algo] = -1;
					output << " & & " << num_optimal_quantile[algo] << "/" << num_inst_per_vertex_size_type_quantile << " & " << avg_time_quantile[algo] << " & " << avg_gap_quantile[algo] << " & " << st_dev_quantile[algo];
				}
				output << "\\\\" << std::endl;
			}
		}
		output << "Sub-total";

		for(size_t algo = 0; algo < algorithms.size(); ++algo)
		{
			if((time_per_algo_inst_size[algo]).size() > 0) avg_time_inst_size[algo]/=(1.0*((time_per_algo_inst_size[algo]).size()));
			else avg_time_inst_size[algo] = -1;
			if((gap_per_algo_inst_size[algo]).size() > 0)
			{
				avg_gap_inst_size[algo]/=(1.0*((gap_per_algo_inst_size[algo]).size()));
				st_dev_inst_size[algo] = StDev(gap_per_algo_inst_size[algo],avg_gap_inst_size[algo]);
			}else avg_gap_inst_size[algo] = st_dev_inst_size[algo] = -1;
			output << " & & " << num_optimal_inst_size[algo] << "/" << num_inst_per_vertex_size << " & " << avg_time_inst_size[algo] << " & " << avg_gap_inst_size[algo] << " & " << st_dev_inst_size[algo];
		}
		output << "\\\\" << std::endl;
	}

	output << "Total";
	for(size_t algo = 0; algo < algorithms.size(); ++algo)
	{
		//std::cout << total_improvement_per_algo[j].size() << std::endl;
		if((total_time_per_algo[algo]).size() > 0) total_avg_time[algo]/=(1.0*((total_time_per_algo[algo]).size()));
		else total_avg_time[algo] = -1;
		if((total_gap_per_algo[algo]).size() > 0) total_avg_gap[algo]/=(1.0*((total_gap_per_algo[algo]).size()));
		else total_avg_gap[algo] = -1;
		output << "& & " << total_num_optimal[algo] << "/" << total_num_instances << " & " << total_avg_time[algo] << " & " << total_avg_gap[algo] << " & " << StDev(total_gap_per_algo[algo],total_avg_gap[algo]);
	}

	output << "\\\\" << std::endl;

	output.close();
}

void GenerateLPImprovementsLatexTable()
{
	std::string curr_file;
	std::vector<std::string> algorithms;
	std::vector<std::string> instance_sizes;
	std::vector<std::string> instance_types;
	std::vector<std::string> instance_limit_quantiles;

	std::vector<std::string> num_vehicles_vec{"2","3","4","5"};
	std::vector<std::string> service_time_deviation_vec{"0.10","0.25","0.50"};
	std::vector<std::string> uncertainty_budget_vec{"0","1","5"};

	std::vector<std::string> instances_terminations;

	for(auto num_vehicles: num_vehicles_vec)
		for(auto service_time_deviation: service_time_deviation_vec)
			for (auto uncertainty_budget: uncertainty_budget_vec)
				instances_terminations.push_back("_v" + num_vehicles + "_d" + service_time_deviation + "_b" + uncertainty_budget + ".txt");

	algorithms.push_back("relax_baseline");
	algorithms.push_back("relax_cb_baseline");
	algorithms.push_back("relax_csc");
	algorithms.push_back("relax_cb_csc");

	instance_types.push_back("C");
	instance_types.push_back("R");
	instance_types.push_back("RC");

	instance_sizes.push_back("25");
	instance_sizes.push_back("50");
	//instance_sizes.push_back("100");

	instance_limit_quantiles.push_back("0.8");

	// instance_limit_quantiles.push_back("1");
	// instance_limit_quantiles.push_back("2");
	// instance_limit_quantiles.push_back("3");
	// instance_limit_quantiles.push_back("4");

	std::fstream output;
	std::string output_name = "..//tables//latex//table_LP_improvements.txt";
	output.open(output_name.c_str(),std::fstream::out);

	if(!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	//std::cout << output_name << std::endl;

	output << std::setprecision(2) << std::fixed;

	std::vector<std::vector<double>> total_improvement_per_algo(algorithms.size(),std::vector<double>());
	std::vector<double> total_avg_improvement(algorithms.size(),0.0);

	for(auto instance_size : instance_sizes)
	{
		std::vector<std::vector<double>> improvement_per_algo_inst_size(algorithms.size(),std::vector<double>());
		std::vector<double> avg_improvement_per_inst_size(algorithms.size(),0.0);
		std::vector<double> st_dev_per_inst_size(algorithms.size(),0.0);
		for(auto instance_type : instance_types)
		{
			for(auto instance_limit_quantile : instance_limit_quantiles)
			{
				std::vector<std::vector<double>> improvement_per_algo(algorithms.size(),std::vector<double>());
				std::vector<double> avg_improvement(algorithms.size(),0.0);
				std::vector<double> st_dev(algorithms.size(),0.0);

				std::string instance_prefix = instance_type + instance_size + "_" + instance_limit_quantile;

				for(auto instances_termination: instances_terminations)
				{
					std::string instance = instance_prefix + instances_termination;
					std::cout << instance << std::endl;
					double original_lp = 0.0;
					for(size_t algo = 0; algo < algorithms.size(); ++algo)
					{
						double curr_improvement = 0.0;
						curr_file = "..//solutions//";
						curr_file.append("s_");
						curr_file.append(algorithms[algo]);
						curr_file.append("_");
						curr_file.append(instance);
						//std::cout << curr_file << std::endl;

						std::fstream input;
						input.open(curr_file.c_str(),std::fstream::in);

						if(!input.is_open())
						{
							std::cout << "Could not open file " << curr_file << std::endl;
							continue;
						}

						std::stringstream s_lp;
						std::string status;
						double lp = 0.0;
						std::string line;

						getline(input,line);
						size_t pos = line.find_first_of(":");
						status = line.substr(pos + 2);

						std::string::iterator end_pos = std::remove(status.begin(),status.end(),' ');
						status.erase(end_pos, status.end());

						getline(input,line);
						pos = line.find_first_of(":");
						s_lp << line.substr(pos + 2);

						if(s_lp.str() == "inf") lp = -1;
						else s_lp >> lp;

						if(algo == 0)
						{
							if((double_equals(lp,-1)) || (status == "INFEASIBLE"))
							{
								original_lp = -1.0;
							}else original_lp = lp;
							curr_improvement = 0.0;
						}else
						{
							if((double_equals(lp,-1)) || (status == "INFEASIBLE") || (double_equals(original_lp,-1)))
							{
								curr_improvement = -1;
							}else
							{
								if(double_equals(original_lp,0.0)) curr_improvement = 0.0;
								else curr_improvement = (100*(original_lp-lp))/original_lp;

								if(double_less(curr_improvement,0)) std::cout << original_lp << " - " << lp << std::endl;
								if(!double_greater(original_lp,lp)) std::cout << original_lp << " - " << lp << std::endl;
							}
						}

						//std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

						if(!double_equals(curr_improvement,-1))
						{
							improvement_per_algo[algo].push_back(curr_improvement);
							improvement_per_algo_inst_size[algo].push_back(curr_improvement);
							total_improvement_per_algo[algo].push_back(curr_improvement);
							total_avg_improvement[algo] += curr_improvement;
							avg_improvement[algo] += curr_improvement;
							avg_improvement_per_inst_size[algo] += curr_improvement;
						}
						input.close();
					}

					//getchar();getchar();
				}

				output << instance_type + instance_size + "Q" + instance_limit_quantile;

				for(size_t algo = 1; algo < algorithms.size(); ++algo)
				{
					if(improvement_per_algo[algo].size() > 0)
						avg_improvement[algo]/=(1.0*((improvement_per_algo[algo]).size()));
					else avg_improvement[algo] = -1;
					st_dev[algo] = StDev(improvement_per_algo[algo],avg_improvement[algo]);
					output << " & & "<< avg_improvement[algo] << " & " << st_dev[algo];
				}
				output << "\\\\" << std::endl;
			}
		}
		output << "Subtotal";

		for(size_t algo = 1; algo < algorithms.size(); ++algo)
		{
			if ((improvement_per_algo_inst_size[algo]).size() > 0)
				avg_improvement_per_inst_size[algo]/=(1.0*((improvement_per_algo_inst_size[algo]).size()));
			else avg_improvement_per_inst_size[algo] = -1;
			st_dev_per_inst_size[algo] = StDev(improvement_per_algo_inst_size[algo],avg_improvement_per_inst_size[algo]);
			output << " & & "<< avg_improvement_per_inst_size[algo] << " & " << st_dev_per_inst_size[algo];
		}
		output << "\\\\" << std::endl;
	}

	output << "Total";
	for(size_t algo = 1; algo < algorithms.size(); algo++)
	{
		//std::cout << total_improvement_per_algo[j].size() << std::endl;
		if((total_improvement_per_algo[algo]).size() >0)
			total_avg_improvement[algo]/=(1.0*((total_improvement_per_algo[algo]).size()));
		else total_avg_improvement[algo] = -1;
		output << "& & "<< total_avg_improvement[algo] << " & " << StDev(total_improvement_per_algo[algo],total_avg_improvement[algo]);
	}

	output << "\\\\" << std::endl;

	output.close();
}

int main()
{
	// GenerateLPImprovementsLatexTable();
	//GenerateAlgorithmsLatexTable(3600);
	//return 0;
	int time_limit = -1;
	Instance inst("/home/lucas/Documentos/Research/msca-humanitarian-optimization/instances/R-STOP-DP/","test.txt",1,0.25,1,false);

	std::vector<bool> * CALLBACKS_SELECTION = GetCallbackSelection();
	(*CALLBACKS_SELECTION)[K_TYPE_CLIQUE_CONFLICT_CUT] = true;
	(*CALLBACKS_SELECTION)[K_TYPE_INITIAL_ARC_VERTEX_INFERENCE_CUT] = false;

	//std::cout << inst << std::endl;
	const Graph * graph = inst.graph();

	//std::cout << *graph << std::endl;

	double * R0 = Dijkstra(graph,false,false);
	double * Rn = Dijkstra(graph,true,false);

	// for (int i =0; i < graph->num_vertices(); ++i)
	// {
	// 	std::cout << "R0," << i << " = " << R0[i] << std::endl;
	// }

	//inst.WriteToFile("../","teste.txt");
	bool force_use_all_vehicles = false;
	bool export_model = true;
	bool solve_relaxed = false;
	bool use_valid_inequalities = false;
	bool combine_feas_op_cuts = true;
	bool apply_benders_generic_callback = false;
	auto root_cuts = new std::list<UserCutGeneral*>();
	Solution<double> solution(graph->num_vertices());
	BendersFormulation formulation = BendersFormulation::single_commodity;

	// CompactBaseline(inst,R0,Rn,-1,true,false,false,nullptr,nullptr,force_use_all_vehicles,export_model, nullptr,solution);
	// std::cout <<  " LP: " << solution.lp_ << std::endl;
	// std::cout << "num cuts: " << solution.num_cuts_found_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << solution.num_cuts_added_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;

	// solution.reset();
	// CompactBaseline(inst,R0,Rn,-1,true,use_valid_inequalities,true,nullptr,nullptr,force_use_all_vehicles,export_model, root_cuts,solution);
	// std::cout <<  " LP: " << solution.lp_ << std::endl;
	// std::cout << "num cuts: " << solution.num_cuts_found_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << solution.num_cuts_added_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;

	// solution.reset();
	// CompactBaseline(inst,R0,Rn,-1,false,use_valid_inequalities,false,root_cuts,nullptr,force_use_all_vehicles,export_model, nullptr,solution);
	
	// if(!solution.is_feasible_)
	// 	std::cout << "Infeasible" << std::endl;
	// solution.is_optimal_? std::cout <<  " optimal: " << solution.lb_ << std::endl
	// : std::cout <<  " non optimal: [" << solution.lb_ << ", " << solution.ub_ << "]" << std::endl;
	// std::cout << "num cuts: " << solution.num_cuts_found_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << solution.num_cuts_added_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;
	//solution.reset();
	//DeleteCuts(root_cuts);

	CompactBaseline(inst,R0,Rn,time_limit,solve_relaxed,use_valid_inequalities,false,nullptr,nullptr,force_use_all_vehicles,export_model, nullptr,solution);

	if (solve_relaxed)
	{
		std::cout <<  " LP: " << solution.lp_ << std::endl;
	}
	else
	{
		solution.is_optimal_? std::cout <<  " optimal: " << solution.lb_ << std::endl
		: std::cout <<  " non optimal: [" << solution.lb_ << ", " << solution.ub_ << "]" << std::endl;
	}
	std::cout << "num cuts: " << solution.num_cuts_found_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << solution.num_cuts_added_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;

	solution.reset();
	// CompactBaseline(inst,R0,Rn,-1,solve_relaxed,use_valid_inequalities,false,nullptr,nullptr,force_use_all_vehicles,export_model, root_cuts,solution);

	// if (solve_relaxed)
	// {
	// 	std::cout <<  " LP: " << solution.lp_ << std::endl;
	// }
	// else
	// {
	// 	solution.is_optimal_? std::cout <<  " optimal: " << solution.lb_ << std::endl
	// 	: std::cout <<  " non optimal: [" << solution.lb_ << ", " << solution.ub_ << "]" << std::endl;
	// }
	// std::cout << "num cuts: " << solution.num_cuts_found_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << solution.num_cuts_added_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;

	// solution.reset();
	// auto conflicts = inst.conflicts_list();
	// std::cout << "conflicts: " << std::endl;
	// for (auto &conflict: conflicts)
	// {
	// 	for( auto& vertex: conflict)
	// 		std::cout << vertex << " ";
		
	// 	std::cout << std::endl;
	// }

	// CompactSingleCommodity(inst,R0,Rn,-1,false,use_valid_inequalities,false,root_cuts,nullptr,force_use_all_vehicles,export_model,nullptr,solution);

	// if(!solution.is_feasible_)
	// 	std::cout << "Infeasible" << std::endl;
	// solution.is_optimal_? std::cout <<  " optimal: " << solution.lb_ << std::endl
	// : std::cout <<  " non optimal: [" << solution.lb_ << ", " << solution.ub_ << "]" << std::endl;
	// std::cout << "num cuts: " << solution.num_cuts_found_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << solution.num_cuts_added_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;

	solution.reset();
	DeleteCuts(root_cuts);

	Benders(inst,formulation,R0,Rn,time_limit,apply_benders_generic_callback,false,solve_relaxed,use_valid_inequalities,false,nullptr,nullptr,force_use_all_vehicles,export_model,nullptr,solution);
	if (solve_relaxed)
	{
		std::cout <<  " LP: " << solution.lp_ << std::endl;
	}
	else
	{
		solution.is_optimal_? std::cout <<  " optimal: " << solution.lb_ << std::endl
		: std::cout <<  " non optimal: [" << solution.lb_ << ", " << solution.ub_ << "]" << std::endl;
	}
	std::cout << "# opt cuts: " << solution.num_benders_opt_cuts_ << std::endl;
	std::cout << "# feas cuts: " << solution.num_benders_feas_cuts_ << std::endl;
	std::cout << "num cuts: " << solution.num_cuts_found_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << solution.num_cuts_added_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;

	solution.reset();

	Benders(inst,formulation,R0,Rn,time_limit,apply_benders_generic_callback,true,solve_relaxed,use_valid_inequalities,false,nullptr,nullptr,force_use_all_vehicles,export_model,nullptr,solution);
	if (solve_relaxed)
	{
		std::cout <<  " LP: " << solution.lp_ << std::endl;
	}
	else
	{
		solution.is_optimal_? std::cout <<  " optimal: " << solution.lb_ << std::endl
		: std::cout <<  " non optimal: [" << solution.lb_ << ", " << solution.ub_ << "]" << std::endl;
	}
	std::cout << "# opt cuts: " << solution.num_benders_opt_cuts_ << std::endl;
	std::cout << "# feas cuts: " << solution.num_benders_feas_cuts_ << std::endl;
	std::cout << "num cuts: " << solution.num_cuts_found_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << solution.num_cuts_added_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;
	solution.reset();




	// IloEnv env;
	// IloNumArray y_values(env,inst.graph()->num_vertices());
	// // for(int i = 0; i < inst.graph()->num_vertices(); ++i)
	// // {
	// //   y_values[i] = 0;
	// // }
	// y_values[0] = 1;
	// y_values[1] = 1;

	// IloNumArray x_values(env,inst.graph()->num_arcs());
	
	// for(int i  = 0; i < inst.graph()->num_vertices(); ++i)
	// {
	//   for (const auto & j: inst.graph()->AdjVerticesOut(i))
	//   {
	//     x_values[inst.graph()->pos(i,j)] = 0;
	//   }
	// }

	// x_values[inst.graph()->pos(0,1)] = 1;
	// x_values[inst.graph()->pos(1,0)] = 1;

	// // PrimalSubproblemCompactBaseline(inst,x_values,y_values,R0,Rn,-1,export_model,solution);
	// // std::cout << solution.lp_ << std::endl;
	// // solution.reset();
	// // DualSubproblemCompactBaseline(inst,x_values,y_values,R0,Rn,-1,export_model,solution);
	// // std::cout << solution.lp_ << std::endl;
	// // solution.reset();

	// PrimalSubproblemCompactSingleCommodity(inst,x_values,y_values,R0,Rn,-1,export_model,solution);
	// std::cout << solution.lp_ << std::endl;
	// solution.reset();
	// DualSubproblemCompactSingleCommodity(inst,x_values,y_values,R0,Rn,-1,export_model,solution);
	// std::cout << solution.lp_ << std::endl;

	// env.end();

	delete [] R0; 
	R0 = nullptr;
	delete [] Rn;
	Rn = nullptr;

	DeleteTimer();
	DeleteCallbackSelection();

	return 0;
}
