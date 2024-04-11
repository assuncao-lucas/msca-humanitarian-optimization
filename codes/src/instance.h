#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <iterator>
#include <tuple>
#include <unordered_map>
#include <boost/dynamic_bitset.hpp>
#include "src/matrix.hpp"
#include "src/graph.h"

class Route; // forward declaration.
class Instance
{
private:
	Graph *graph_ = nullptr;
	int num_vehicles_ = 0;
	int uncertainty_budget_ = 0;
	double service_time_deviation_ = 0.0;
	double limit_ = 0.0;
	std::string raw_file_name_;
	int num_mandatory_ = 0;
	bool found_maximal_cliques_ = false;
	double time_spent_in_preprocessing_ = 0.0;
	std::vector<std::list<int>> *conflict_graph_ = nullptr;
	Matrix<bool> conflict_matrix_;
	std::list<int> mandatory_list_;
	std::vector<std::list<int>> conflicts_list_; // list of cliques

	std::vector<std::list<int>> map_vertices_to_cliques_;
	boost::dynamic_bitset<> active_conflict_cliques_;
	std::vector<int> map_reordered_vertices_to_original_positions_;

	void GeneratePreProcessedGraph();
	void ReorderMandatoryVertices();
	struct PairHash
	{
		template <class T1, class T2>
		std::size_t operator()(const std::pair<T1, T2> &p) const
		{
			auto h1 = std::hash<T1>{}(p.first);
			auto h2 = std::hash<T2>{}(p.second);
			return h1 ^ h2;
		}
	};

	struct PairEqual
	{
		template <class T1, class T2>
		bool operator()(const std::pair<T1, T2> &p1, const std::pair<T1, T2> &p2) const
		{
			return p1.first == p2.first && p1.second == p2.second;
		}
	};

	typedef std::unordered_map<std::pair<int, int>, std::tuple<double, double>, PairHash, PairEqual> VertexBudgetHash;

public:
	explicit Instance() = default;
	explicit Instance(std::string dir_path, std::string file_name, int num_vehicles, double service_time_deviation, int uncertainty_budged, bool pre_process_graph);
	virtual ~Instance();
	std::string GetInstanceName() const;
	void FillInstanceFromFile(std::string dir_path, std::string file_name, double service_time_deviation);
	void AddMandatoryVertices(double mandatory_percentage);
	void WriteToFile(std::string dir_path, std::string curr_file);
	void BuildVerticesToCliquesMapping();
	void FindAllMaximalConflictCliquesBronKerosch();
	void FindAllMaximalConflictCliquesTomita();
	void SelectMaximumCliquesPerVertex();
	void ComputeConflictGraph();
	void ResetConflictsCliques();
	int getOriginalVertexPosition(int new_pos) { return map_reordered_vertices_to_original_positions_[new_pos]; } // position before reordering vertices.
	int num_vehicles() const { return num_vehicles_; }
	int num_mandatory() const { return num_mandatory_; }
	bool found_maximal_cliques() { return found_maximal_cliques_; }
	bool FoundConflictGraph() { return conflict_graph_ != nullptr; }
	double limit() const { return limit_; }
	double time_spent_in_preprocessing() const { return time_spent_in_preprocessing_; }
	int uncertainty_budget() const { return uncertainty_budget_; }
	double service_time_deviation() const { return service_time_deviation_; }
	boost::dynamic_bitset<> active_conflict_cliques() { return active_conflict_cliques_; }
	const std::vector<std::list<int>> &map_vertices_to_cliques() const { return map_vertices_to_cliques_; }
	const std::vector<std::list<int>> &conflicts_list() const { return conflicts_list_; }
	void set_graph(Graph *graph);
	const Graph *graph() const { return graph_; }
	std::tuple<double, double> ComputeRouteCostsRecIter(Route &route, const std::list<int>::reverse_iterator &it_vertex, int budget, VertexBudgetHash *cache) const;
	std::tuple<double, double> ComputeRouteCostsRec(Route &route, bool memoization) const;
	std::tuple<double, double> ComputeRouteCosts(Route &route) const;
	friend std::ostream &operator<<(std::ostream &out, const Instance &instance);
};
