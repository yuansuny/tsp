#ifndef SOLVER_H
#define SOLVER_H

#include "graph.h"
#include "reduce.h"
#include <iostream>
#include <vector>
#include <omp.h>

// Magic tricks to have CPLEX behave well:
#ifndef IL_STD
#define IL_STD
#endif
#include <cstring>
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN
// End magic tricks

namespace TSP {
  class Solver {
    // The graph on which we are solving the TSP.
    const Graph& g;
    const Reduce& r;
    const double cutoff;
    double objVal;
    std::vector<std::uint32_t> opt_sol;
    
    // Prints the integer solution obtained by CPLEX to stdout.
    void print_solution(const IloCplex& cplex, const IloArray<IloNumVarArray>& x);
    
  public:
    
    // Builds a solver for graph g.
    explicit Solver(const Graph& g, const Reduce& r, const float cutoff) : g{g}, r{r}, cutoff{cutoff} {}
    
    // Solves the TSP with CPLEX and prints the result.
    void solve_tsp_cplex();

//    void solve_tsp_concorde();

    double get_objective_value() const { return objVal; }

    std::vector<std::uint32_t> get_optimal_solution() const { return opt_sol; }

  };
}

#endif
