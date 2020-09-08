#include "solver.h"
#include <cmath>
#include <limits>

//extern "C" {
//#include <concorde.h>
//}

namespace TSP {
//    //solve TSP instances using Concorde
//    void Solver::solve_tsp_concorde() {
//        int rval = 0;
//        int n = g.size();
//        int ecount = r.get_reduced_problem_size()/2;
////        std::uint64_t ecount = (n * (n - 1)) / 2;
//
////        int *elist = CC_SAFE_MALLOC(ecount*2, int);
////        int *elen = CC_SAFE_MALLOC(ecount, int);
//        int *elist = new int[ecount*2];
//        int *elen = new int[ecount];
//        int *in_tour = (int *) NULL; //Gives a starting tour in node node node format (it can be NULL)
////        int *in_tour = new int[n];
////        std::vector<std::uint32_t> sol = r.get_best_sol_sampling();
////        for (int i = 0; i < n; ++i){
////            in_tour[i] = sol[i];
////        }
////        double in_val = r.get_objective_value_sampling();
//
//        int *out_tour = (int *) NULL; //Optimal tour (it can be NULL, if it is not NULL then it should point to an array of length at least ncount).
//        out_tour = CC_SAFE_MALLOC(n, int);
//        double *in_val = (double *) NULL; //Can be used to specify an initial upperbound (it can be NULL)
//        double optval; //Value of the optimal tour
//        int success; //1 if the run finished normally, and set to 0 if the search was terminated early (by hitting some predefined limit)
//        int foundtour; //1 if a tour has been found (if success is 0, then it may not be the optimal tour)
////        char *filename = (char *) NULL; //Specifes a char string that will be used to name various files that are written during the branch and bound search
//        char *filename = new char[g.get_file_name().size()];
//        strcpy(filename, g.get_file_name().c_str());
//        double timebound = cutoff; //Run time limit
////        double *timebound = (double *) NULL; //Run time limit
//        int hit_timebound = 0; //1 if timebound was reached
//        static int silent = 1; //Suppress most output if set to a nonzero value
//
//        CCrandstate rstate;
//        int seed = std::rand();
//        CCutil_sprand(seed, &rstate); //Initialize the portable random number generator
//
//        std::uint64_t count = 0u;
//        for (int i = 0; i < n; i++){
//            for (int j = i+1; j < n; j++){
//                if (r.get_predicted_value(i,j) == 1){
//                    elist[2*count] = i;
//                    elist[2*count+1] = j;
//                    elen[count] = (int) g.cost(i, j);
//                    count++;
//                }
////                else{
////                    elist[2*count] = i;
////                    elist[2*count+1] = j;
////                    elen[count] = 99999999;
////                    count++;
////                }
//            }
//        }
//
//        rval = CCtsp_solve_sparse(n, ecount, elist, elen, in_tour, out_tour, in_val, &optval, &success, \
//            &foundtour, filename, &timebound, &hit_timebound, silent, &rstate);
//        std::cout << "rval is " << rval << std::endl;
//
//        //Print solution
//        objVal = optval;
//        std::cout << "optimal objective value found is " << optval << std::endl;
//
//        opt_sol = std::vector<std::uint32_t>(n);
//        for (int i = 0; i < n; ++i){
//            opt_sol[i] = out_tour[i];
//        }
//
////        szeit = CCutil_zeit();
//
//        //free memory
//        CC_IFFREE (elist, int);
//        CC_IFFREE (elen, int);
//        CC_IFFREE (out_tour, int);
//        CC_IFFREE (filename, char);
//    }

    void Solver::solve_tsp_cplex() {
        const auto n = g.size();

        // CPLEX environment. Takes care of everything, including memory management for CPLEX objects.
        IloEnv env;

        // CPLEX model. We put variables and constraints in it!
        IloModel model(env);

        // Model:
        //
        // BINARY VARIABLE x[i][j]    For all i,j = 0, ..., n - 1
        //    x[i][j] == 1            If arc (i,j) is selected
        //    x[i][j] == 0            Otherwise
        //
        // INTEGER VARIABLE t[i]      For all i = 0, ..., n - 1
        //    t[i] == k               Iff node i is the k-th node in the tour
        //    t[0] == 1
        //    t[i] in [2, ..., n]     For all i = 1, ... n - 1
        //
        // OBJECTIVE FUNCTION
        //    MIN sum((i,j), c[i][j] * x[i][j])
        //
        // CONSTRAINTS
        //    1) sum(j, x[j][i]) == 1                    For all i
        //    2) sum(j, x[i][j]) == 1                    For all i
        //    3) t[i] - t[j] + 1 <= n * (1 - x[i][j])    For all i,j = 1, ..., n - 1
        //       Can be written as:
        //       t[i] - t[j] + n * x[i][j] <= n - 1

        // Variables
        IloArray<IloNumVarArray> x(env, n);
        IloNumVarArray t(env, n);

        // Constraints
        IloRangeArray inbound_arcs(env, n);  // Constraints 1)
        IloRangeArray outbound_arcs(env, n); // Constraints 2)
        IloArray<IloRangeArray> mtz(env, n); // Constraints 3)
        //    IloRangeArray removed_arcs(env, 2); // Constraints 4)

        // We use this stringstream to create variable and constraint names
        std::stringstream name;

        // Create variable t[0] and fix it to value 1
        // This breaks symmetry, because it fixes node 0 as the starting node of the tour
        t[0] = IloNumVar(env, 1, 1, IloNumVar::Int, "t_0");

        // Create variables t[1], ..., t[n]
        for(auto i = 1u; i < n; ++i) {
            name << "t_" << i;
            t[i] = IloNumVar(env, 2, n, IloNumVar::Int, name.str().c_str());
            name.str(""); // Clean name
        }

        // Create variables x
        for(auto i = 0u; i < n; ++i) {
            x[i] = IloNumVarArray(env, n);
            for(auto j = 0u; j < n; ++j) {
            name << "x_" << i << "_" << j;
            x[i][j] = IloNumVar(env, 0, 1, IloNumVar::Bool, name.str().c_str());
            name.str(""); // Clean name
            }
        }

        IloExpr expr(env);

        // removing variables by setting upper bound = lower bound = 0
        for (auto i = 0u; i < n; ++i){
            for (auto j = i+1; j < n; ++j){
                if (r.get_predicted_value(i,j) == 0){
                    x[i][j].setUB(0);
                    x[j][i].setUB(0);
                }
            }
        }


        // Create constraints 1)
        for(auto i = 0u; i < n; ++i) {
            for(auto j = 0u; j < n; ++j) {
                expr += x[j][i];
            }
            name << "inbound_" << i;
            inbound_arcs[i] = IloRange(env, 1, expr, 1, name.str().c_str());
            name.str(""); // Clean name
            expr.clear(); // Clean expr
        }

        // Add constraints 1) to the model
        model.add(inbound_arcs);

        // Create constraints 2)
        for(auto i = 0u; i < n; ++i) {
            for(auto j = 0u; j < n; ++j) {
                expr += x[i][j];
            }
            name << "outbound_" << i;
            outbound_arcs[i] = IloRange(env, 1, expr, 1, name.str().c_str());
            name.str(""); // Clean name
            expr.clear(); // Clean expr
        }
        // Add constraints 2) to the model
        model.add(outbound_arcs);

        // Create constraints 3)
        // The constraint is for i = 1,...,n and therefore we add empty constraints for i == 0
        mtz[0] = IloRangeArray(env);
        // We then continue normally for all other i > 0
        for(auto i = 1u; i < n; ++i) {
            mtz[i] = IloRangeArray(env, n);
            for(auto j = 1u; j < n; ++j) {
                expr = t[i] - t[j] + static_cast<int>(n) * x[i][j];

                name << "mtz_" << i << "_" << j;
                mtz[i][j] = IloRange(env, -IloInfinity, expr, n - 1, name.str().c_str());
                name.str(""); // Clean name
                expr.clear(); // Clean expr
            }
            // Add constraints 3)[i] to the model
            model.add(mtz[i]);
        }

        // Create objective function
        for(auto i = 0u; i < n; ++i) {
            for(auto j = 0u; j < n; ++j) {
                expr += g.cost(i, j) * x[i][j];
            }
        }
        IloObjective obj(env, expr, IloObjective::Minimize);

        // Add the objective function to the model
        model.add(obj);

        // Free the memory used by expr
        expr.end();

        // Create the solver object
        IloCplex cplex(model);

        // Export model to file (useful for debugging!)
        cplex.exportModel("model.lp");

        //Set cutoff time
        cplex.setParam(IloCplex::TiLim, cutoff);

        // Set number of threads
        cplex.setParam(IloCplex::Threads, omp_get_max_threads());

        //Set random seed
        srand (time(NULL));
        cplex.setParam(IloCplex::RandomSeed, rand()%10000);


        // warm start
        IloNumVarArray startVar(env);
        IloNumArray startVal(env);
        std::vector<std::uint32_t> sol = r.get_best_sol_sampling();
        std::vector<std::vector<std::uint32_t>> start(n, std::vector<std::uint32_t>(n, 0));
        for (auto i = 0u; i < n - 1; ++i){
            start[sol[i]][sol[i+1]] = 1;
        }
        start[sol[n-1]][sol[0]] = 1;
        for(auto i = 0u; i < n; ++i){
            for (auto j = 0u; j < n; ++j) {
                startVar.add(x[i][j]);
                startVal.add(start[i][j]);
            }
        }
        cplex.addMIPStart(startVar, startVal);
        startVal.end();
        startVar.end();


        bool solved = false;
        try {
            // Try to solve with CPLEX (and hope it does not raise an exception!)
            solved = cplex.solve();
        } catch(const IloException& e) {
            std::cerr << "\n\nCPLEX Raised an exception:\n";
            std::cerr << e << "\n";
            env.end();
            throw;
        }

        if(solved) {
            // If CPLEX successfully solved the model, print the results
            std::cout << "\n\nCplex success!\n";
            std::cout << "\tStatus: " << cplex.getStatus() << "\n";
            objVal = cplex.getObjValue();
            std::cout << "\tObjective value: " << cplex.getObjValue() << "\n";
            print_solution(cplex, x);
        } else {
            std::cerr << "\n\nCplex error!\n";
            std::cerr << "\tStatus: " << cplex.getStatus() << "\n";
            std::cerr << "\tSolver status: " << cplex.getCplexStatus() << "\n";
        }
        env.end();
    }

    void Solver::print_solution(const IloCplex& cplex, const IloArray<IloNumVarArray>& x) {
        const auto n = g.size();
        assert(x.getSize() == n);
        opt_sol = std::vector<std::uint32_t>(n);
        std::cout << "\n\nTour: ";
        const auto starting_vertex = 0u;
        auto current_vertex = starting_vertex;
        auto k = 0u;
        do {
            opt_sol[k++] = current_vertex;
            std::cout << current_vertex << " ";
            for(auto i = 0u; i < n; ++i) {
                if(cplex.getValue(x[current_vertex][i]) > .5) {
                    current_vertex = i;
                    break;
                }
            }
        } while(current_vertex != starting_vertex);
        std::cout << "\n";
    }
}
