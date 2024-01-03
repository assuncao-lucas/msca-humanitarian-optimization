#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include "instance.h"
#include "general.h"
#include "graph_algorithms.h"
#include "timer.h"

void Instance::FillInstanceFromFile(std::string dir_path, std::string file_name, double service_time_deviation)
{
	std::string curr_file = dir_path;
	curr_file.append(file_name);
	Graph * graph = nullptr;

	//std::cout << dir_path << " " << file_name << std::endl;

	std::fstream file;
	file.open(curr_file.c_str(), std::fstream::in);
	if(!file.is_open()){
		std::cout << "Could not open file" << std::endl;
		throw 1;
		return;
	}

	int num_vertices = 0, mandatory_iter = 0;

    file >> num_vertices >> num_mandatory_;

    //std::cout << num_vertices << " " << num_mandatory << std::endl;

    for(int i= 0; i < num_mandatory_; ++i)
	{
		file >> mandatory_iter;
		mandatory_list_.push_back(mandatory_iter);
	}

	std::vector<std::pair<double,double>> coordinates(num_vertices);
	Graph::VertexInfo * vertices_info = new Graph::VertexInfo[num_vertices];

    // fill Vertices info.
    double garbage = 0.0;
	for(int i = 0; i < num_vertices; ++i)
    {
		file >> vertices_info[i].coordinates_.first >> vertices_info[i].coordinates_.second >> vertices_info[i].profit_ >> vertices_info[i].decay_ratio_ >> vertices_info[i].nominal_service_time_ >> garbage;
        vertices_info[i].dev_service_time_ = round_decimals(vertices_info[i].nominal_service_time_ * service_time_deviation,2);
    }
    assert(vertices_info[0].profit_ == 0);
    assert(double_equals(vertices_info[0].decay_ratio_,0.0));
    assert(double_equals(vertices_info[0].nominal_service_time_,0.0));
    assert(double_equals(vertices_info[0].dev_service_time_,0.0));

    // read route time limit.
    file >> limit_;

	file.close();

	graph_ = new Graph(num_vertices,vertices_info);

    for(int i = 0; i < num_vertices; ++i)
        for(int j = i+1; j < num_vertices; ++j)
            graph_->AddEdge(i,j,round_decimals(euclidian_distance(vertices_info[i].coordinates_, vertices_info[j].coordinates_),2));
}

void Instance::AddMandatoryVertices(double mandatory_percentage)
{
    int num_vertices = graph_->num_vertices();
    Graph::VertexInfo * vertices_info = new Graph::VertexInfo[num_vertices];

    int num_mandatory_ = (int)ceil(mandatory_percentage * num_vertices);
    std::vector<bool> selected_vertices(num_vertices, false);
    int iter_mandatory = 0, cont_mandatory = 0;
    while(cont_mandatory < num_mandatory_)
    {
        iter_mandatory = rand()%(num_vertices-2) + 1;
        if(!selected_vertices[iter_mandatory])
        {
            selected_vertices[iter_mandatory] = true;
            ++cont_mandatory;
            mandatory_list_.push_back(iter_mandatory);
            vertices_info[iter_mandatory].profit_ = 0;
            vertices_info[iter_mandatory].decay_ratio_ = 0;
        }
    }
}

void Instance::ReorderMandatoryVertices()
{
    int num_vertices = graph_->num_vertices();
    const Graph::VertexInfo * vertices_info = graph_->vertices_info();
    Graph::VertexInfo * reordered_vertices_info = new Graph::VertexInfo[num_vertices];
    Graph * new_graph = nullptr;
    std::vector<bool> reordered(num_vertices,false);
    map_reordered_vertices_to_original_positions_ = std::vector<int>(num_vertices,-1);
    std::vector<int> new_pos(num_vertices,-1);
    reordered[0] = true;
    new_pos[0] = 0;
    map_reordered_vertices_to_original_positions_[0] = 0; // origin remains at the same position.
    reordered_vertices_info[0] = vertices_info[0];
    int cont = 1;
    for(auto mandatory_vertex: mandatory_list_)
    {
        reordered[mandatory_vertex] = true;
        map_reordered_vertices_to_original_positions_[cont] = mandatory_vertex;
        new_pos[mandatory_vertex] = cont;
        reordered_vertices_info[cont] = vertices_info[mandatory_vertex];
        // as they are mandatory now, set as zero.
        reordered_vertices_info[cont].profit_ = 0;
        reordered_vertices_info[cont].decay_ratio_ = 0.0;
        ++cont;
    }

    // skip origin (zero), as already added.
    for(int i = 1; i < num_vertices; ++i)
    {
        if(!(reordered[i]))
        {
            reordered[i] = true;
            map_reordered_vertices_to_original_positions_[cont] = i;
            new_pos[i] = cont;
            reordered_vertices_info[cont] = vertices_info[i];
            ++cont;
        }
    }

    new_graph = new Graph(num_vertices,reordered_vertices_info);

    // add updated arcs to new graph.
    for(int i = 0; i < num_vertices; ++i)
    {
        for(auto adj_vertex: graph_->AdjVerticesOut(i))
        {
            GArc * curr_arc = (*graph_)[i][adj_vertex];
            new_graph->AddArc(new_pos[i],new_pos[adj_vertex],curr_arc->distance());
        }
    }

    // std::cout << "OLD" << std::endl;
    // std::cout << *graph_ << std::endl;
    // std::cout << "NEW" << std::endl;
    // std::cout << *new_graph << std::endl;

    delete graph_;
    graph_ = new_graph;
}

Instance::Instance(std::string dir_path, std::string file_name, int num_vehicles, double service_time_deviation, int uncertainty_budged, bool pre_process_graph):
num_vehicles_(num_vehicles), service_time_deviation_(service_time_deviation), uncertainty_budget_(uncertainty_budged), raw_file_name_(file_name)
{
    FillInstanceFromFile(dir_path,file_name,service_time_deviation);

    ReorderMandatoryVertices();

	if(pre_process_graph) GeneratePreProcessedGraph();
}

Instance::~Instance()
{
    if(graph_)
    {
        delete graph_;
        graph_ = nullptr;
    }

    if(conflict_graph_)
    {
        delete conflict_graph_;
        conflict_graph_ = nullptr;
    }
}

std::string Instance::GetInstanceName() const
{
    std::ostringstream stream;
    stream << std::fixed;
    stream << std::setprecision(2);
    stream << service_time_deviation_;

    size_t pos = raw_file_name_.find(".txt");

    std::string raw_file_name_without_txt = raw_file_name_.substr(0,raw_file_name_.length() - pos+1);

    return raw_file_name_without_txt + "_v" + std::to_string(num_vehicles_) + "_d" + stream.str() + "_b" + std::to_string(uncertainty_budget_) + ".txt";
}

void Instance::set_graph(Graph * graph)
{
  graph_ = graph;
}

void Instance::GeneratePreProcessedGraph()
{
    Timer * timer = GetTimer();
    Timestamp * ti = NewTimestamp();
    timer->Clock(ti);
    Matrix<double> * min_dists = nullptr;
    Graph * graph = graph_;
    int v1 = 0, v2 = 0;
    int num_vertices = 0;
    int cont = 0;
    GArc * curr_arc =  nullptr;
    const auto * vertices_info = graph->vertices_info();
    double route_length = 0.0;

    if(graph)
    {
        num_vertices = graph->num_vertices();
        min_dists = FloydWarshall(graph);

        std::list<int>::iterator it_to_remove;

        for(v1 = 0; v1 < graph->num_vertices(); ++v1)
        {
            for(std::list<int>::iterator it = (graph->AdjVerticesOut(v1)).begin(); it != (graph->AdjVerticesOut(v1)).end();)
            {
                v2 = *it;
                curr_arc = (*graph)[v1][v2];
                route_length = (*min_dists)[0][v1] + vertices_info[v1].nominal_service_time_ + curr_arc->distance() + vertices_info[v2].nominal_service_time_ + (*min_dists)[v2][0];
                if(double_greater(route_length,limit_))
                {
                    // iterates before removing element from list!
                    it_to_remove = it;
                    ++it;
                    (graph->AdjVerticesOut(v1)).erase(it_to_remove);
                    (graph->AdjVerticesIn(v2)).remove(v1);
                    // removes arc from graph, since it cannot be in any feasible solution
                    graph->DeleteArc(v1,v2,false);
                    //std::cout << "removeu arco (" << v1 << "," << v2 << "): " << route_length << std::endl;
                    ++cont;
                }else ++it;
            }
        }

        delete min_dists;
        min_dists = nullptr;
        std::cout << "# arcs deleted: " << cont << std::endl;
        cont = 0;
        graph->set_indexes(new Matrix<int>(num_vertices,num_vertices,-1));

        //updates positions of arcs
        for(v1 = 0; v1 < graph->num_vertices(); ++v1)
        {
            auto& adj_vertices = graph->AdjVerticesOut(v1);
            for(std::list<int>::iterator it = adj_vertices.begin(); it != adj_vertices.end(); ++it)
            {
                v2 = *it;
                graph->SetArcIndex(v1,v2,cont);
                ++cont;
            }
        }
    }

    time_spent_in_preprocessing_ = timer->CurrentElapsedTime(ti);
    delete ti;
    ti = NULL;
}

void Instance::BuildVerticesToCliquesMapping()
{
	// sort cliques in non-increasing order of size
	if(!K_ONLY_MAXIMUM_CLIQUES_PER_VERTEX)
	{
		std::sort((conflicts_list_).begin(), (conflicts_list_).end(), [] (const std::list<int>& lhs, const std::list<int>& rhs)
		{
	    		return lhs.size() > rhs.size();
		});
	}

    	int num_vertices = (graph_)->num_vertices(), cont = 0;
	map_vertices_to_cliques_ = std::vector<std::list<int>>(num_vertices,std::list<int>());

	for(std::vector<std::list<int>>::iterator it = (conflicts_list_).begin(); it != (conflicts_list_).end(); ++it)
	{
		//std::cout << (*it).size() << std::endl;
		for(std::list<int>::iterator it2 = (*it).begin(); it2 != (*it).end(); ++it2)
		{
			((map_vertices_to_cliques_)[*it2]).push_back(cont);
		}
		++cont;
	}
}

void Instance::SelectMaximumCliquesPerVertex()
{
	// sort cliques in non-increasing order of size
	std::sort((conflicts_list_).begin(), (conflicts_list_).end(), [] (const std::list<int>& lhs, const std::list<int>& rhs)
	{
    		return lhs.size() > rhs.size();
	});
	//std::cout << (graph_)->num_vertices() << std::endl;
	//std::cout << "before: " << conflicts_list_.size() << std::endl;

	std::vector<std::list<int>> new_conflicts_list;
	boost::dynamic_bitset<> vertices_set_to_a_maximum_clique((graph_)->num_vertices(),0);

	for(auto it = (conflicts_list_).begin(); it != (conflicts_list_).end(); ++it)
	{
		bool add_clique = false;
		for(auto it2 = it->begin(); it2 != it->end(); ++it2)
		{
			if(vertices_set_to_a_maximum_clique[*it2] == 0)
			{
				add_clique = true;
				vertices_set_to_a_maximum_clique[*it2] = 1;
			}
		}
		if(add_clique) new_conflicts_list.push_back(*it);
	}

	conflicts_list_ = new_conflicts_list;

	//std::cout << "after: " << conflicts_list_.size() << std::endl;
	//getchar();getchar();
}

void Instance::FindAllMaximalConflictCliquesBronKerosch()
{
	conflicts_list_.clear();
	size_t num_vertices = conflict_graph_->size();
	std::unordered_set<int> P,X;
	std::list<int> R;

	for(int i = 0; i < num_vertices; i++) P.insert(i);

	BronKerosch(conflict_graph_,R,P,X,conflicts_list_);
	if(K_ONLY_MAXIMUM_CLIQUES_PER_VERTEX) SelectMaximumCliquesPerVertex();
	found_maximal_cliques_ = true;
}

void Instance::FindAllMaximalConflictCliquesTomita()
{
	conflicts_list_.clear();
	size_t num_vertices = conflict_graph_->size();
	boost::dynamic_bitset<> P(num_vertices), X(num_vertices);
	P.set();
	std::list<int> R;

	// build conflic graph based on bitsets
	std::vector<boost::dynamic_bitset<>> bitsets_conflict_graph(num_vertices,boost::dynamic_bitset<>(num_vertices));
	for(size_t i = 0; i < num_vertices; ++i)
	{
		for(std::list<int>::iterator it = ((*conflict_graph_)[i]).begin(); it != ((*conflict_graph_)[i]).end(); ++it)
			bitsets_conflict_graph[i][*it] = 1;

		//std::cout << bitsets_conflict_graph[i] << std::endl;
	}

	Tomita(&bitsets_conflict_graph,R,P,X,conflicts_list_);

	if(K_ONLY_MAXIMUM_CLIQUES_PER_VERTEX) SelectMaximumCliquesPerVertex();
	found_maximal_cliques_ = true;
}

void Instance::ComputeConflictGraph()
{
    Graph * graph = graph_;
    int num_vertices = graph->num_vertices();
    Timer * timer = GetTimer();
    Timestamp * ti = NewTimestamp();
    timer->Clock(ti);
    const auto &vertices_info = graph->vertices_info();

    conflict_matrix_ = Matrix<bool>(num_vertices,num_vertices,false);

    if(conflict_graph_)
    {
        delete conflict_graph_;
        conflict_graph_ = nullptr;
    }
    conflict_graph_ = new std::vector<std::list<int>>(num_vertices,std::list<int>());

    Matrix<double> * min_paths_dist = nullptr;
    double dist1 = 0.0, dist2 = 0.0, dist_j_from_i = 0.0, dist_i_from_j = 0.0, deadline_i = 0.0, deadline_j = 0.0;
    int cont = 0;

    min_paths_dist = FloydWarshall(graph);

    for(int i = 1; i < num_vertices-1; i++)
    {
        if(!(graph->AdjVerticesOut(i).empty())) // if vertice is reachable
        {
            for(int j = i+1; j < num_vertices-1; j++)
            {
                if(!( graph->AdjVerticesOut(j).empty())) // if vertice is reachable
                {
                    dist1 = (*min_paths_dist)[0][i] + vertices_info[i].nominal_service_time_ + vertices_info[j].nominal_service_time_ + (*min_paths_dist)[i][j] + (*min_paths_dist)[j][0];
                    dist2 = (*min_paths_dist)[0][j] + vertices_info[i].nominal_service_time_ + vertices_info[j].nominal_service_time_ + (*min_paths_dist)[j][i] + (*min_paths_dist)[i][0];

                    // if vertices i and j do not belong to a path from s to t that satisfies the length limit.
                    if(double_greater(dist1,limit_) && double_greater(dist2,limit_))
                    {
                        ++cont;
                        ((*(conflict_graph_))[i]).push_back(j);
                        ((*(conflict_graph_))[j]).push_back(i);
                        (conflict_matrix_)[i][j] = (conflict_matrix_)[j][i] =  true;
                        continue;
                    }

                    dist_j_from_i = (*min_paths_dist)[0][i] + vertices_info[i].nominal_service_time_ + (*min_paths_dist)[i][j];
                    dist_i_from_j = (*min_paths_dist)[0][j] + vertices_info[j].nominal_service_time_ + (*min_paths_dist)[j][i];
                    
                    deadline_i = round_decimals(vertices_info[i].profit_/vertices_info[i].decay_ratio_,2);
                    deadline_j = round_decimals(vertices_info[j].profit_/vertices_info[j].decay_ratio_,2);
                    if( double_greater(dist_j_from_i,deadline_j) && double_greater(dist_i_from_j,deadline_i) )
                    {
                        ++cont;
                        std::cout << " COLOCOU A MAIS!" << std::endl;
                        ((*(conflict_graph_))[i]).push_back(j);
                        ((*(conflict_graph_))[j]).push_back(i);
                        (conflict_matrix_)[i][j] = (conflict_matrix_)[j][i] =  true;
                    }
                }
            }
        }
    }

    //std::cout << "num conflicts: " << cont << std::endl;

    delete min_paths_dist;
    min_paths_dist = NULL;
    (time_spent_in_preprocessing_)+= (timer->CurrentElapsedTime(ti));
    delete ti;
}

void Instance::ResetConflictsCliques()
{
	if(!((conflicts_list_).empty()))
	{
		active_conflict_cliques_ = boost::dynamic_bitset<>((conflicts_list_).size(),0);
		(active_conflict_cliques_).set();
	}
}

std::ostream& operator<< (std::ostream &out, const Instance &instance)
{
    out << instance.num_vehicles() << " " << instance.limit() << " "
        << instance.service_time_deviation() << " " << instance.uncertainty_budget() << std::endl;
    
    out << *(instance.graph());
                
    return out;
}

void Instance::WriteToFile(std::string folder, std::string curr_file)
{
    std::fstream output;

    Graph * graph = graph_;
    int num_vertices = graph->num_vertices();
    double limit = limit_;
    const Graph::VertexInfo * vertices_info = graph->vertices_info();

    std::string file_path = folder + curr_file;
	output.open(file_path.c_str(),std::fstream::out);

    //std::cout << folder << " " << curr_file << std::endl;

    if(!output.is_open())
    {
        std::cout << "Could not open file" << std::endl;
        throw 1;
        return;
    }

    output << num_vertices << "\t" << num_mandatory_ << std::endl;

    std::list<int>::iterator last_element = mandatory_list_.end();
    --last_element;
    for(std::list<int>::iterator it = mandatory_list_.begin(); it != mandatory_list_.end(); ++it)
        it != last_element? output << map_reordered_vertices_to_original_positions_[*it] << "\t" : output << map_reordered_vertices_to_original_positions_[*it] << std::endl;

    //output << std::setprecision(2) << std::fixed;

    for(int i = 0; i < num_vertices; ++i)
        output << vertices_info[i].coordinates_.first << "\t"
            << vertices_info[i].coordinates_.second << "\t"
            << vertices_info[i].profit_ << "\t"
            << vertices_info[i].decay_ratio_ << "\t"
            << vertices_info[i].nominal_service_time_ << "\t"
            << vertices_info[i].nominal_service_time_ << std::endl;

    output << limit_;
    output.close();
}
