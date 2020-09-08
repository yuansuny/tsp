#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <sstream>

namespace TSP {
  class Graph {
    // Number of nodes in the graph.
    // Nodes go from 0 to n_nodes - 1.
    std::uint32_t n_nodes;
    const std::string file_name;
    const std::string input_dir;
    std::vector<std::vector<double>> costs;
    std::vector<double> min_costs;
    std::vector<double> max_costs;
    std::vector<double> mean_costs;
    std::vector<std::vector<bool>> optimal_solution;
    std::string name;
    std::string type;
    std::string comment;
    std::string edgeWeightType;
    std::string edgeWeightFormat;
    std::uint64_t opt_obj;
    std::vector<std::uint32_t> opt_tour;

    void read_graph();
    void read_optimal_solution();

    bool checkKeyword(std::string keyword, std::string value);

    std::string trim(std::string);

    
  public:
    
    // Created a new (random) graph with n_nodes nodes
    explicit Graph(std::string file_name, std::string input_dir);

    // Size of the graph (i.e. number of nodes).
    std::uint32_t size() const { return n_nodes; }
    
    // Cost of arc (i,j).
    double cost(std::uint32_t i, std::uint32_t j) const { return costs[i][j]; }

    double get_min_cost(std::uint32_t i) const { return min_costs[i]; }

    double get_max_cost(std::uint32_t i) const { return max_costs[i]; }

    double get_mean_cost(std::uint32_t i) const { return mean_costs[i]; }

    // Optimal solution of edge x[i][j]
    bool get_optimal_value(std::uint32_t i, std::uint32_t j) const { return optimal_solution[i][j]; }

    std::string get_file_name() const { return file_name; }
    
    // Prints out the cost matrix of the graph.
    friend std::ostream& operator<<(std::ostream& out, const Graph& g);
  };
}

#endif
