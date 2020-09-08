#ifndef REDUCE_H
#define REDUCE_H
#include "graph.h"

extern "C" {
#include "svm_predict_model.h"
#include "linear_svm_predict_model.h"
}

#include <random>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <numeric>      // std::iota
#include <vector>
#include <cstring>
#include <string>
#include <time.h>
#include <sys/time.h>
#include <iomanip>
#include <omp.h>

namespace TSP{
    class Reduce{
        // The graph on which we are solving the TSP.
        const Graph& g;
        const float threshold_r = 0.5;
        const float threshold_c = -0.00;
        const std::string test_data_dir = "../test_data/";
        const std::string training_data_dir = "../train_data/";
        const std::string training_model_name = "train_model";
        const std::string test_data_name = "_test_data_";
        const std::string output_file_name = "_predicted_value_";
        const std::uint32_t sample_factor = 100;
        std::uint32_t sample_size = 0;
        int probability = 0;
        int num_thread = 0;
        std::vector<std::vector<std::uint32_t>> samples;
        std::vector<double> objs;
        std::vector<std::vector<double>> mean_x;
        std::vector<std::vector<double>> S1;
        std::vector<double> diff_y;
        double best_obj_sampling;
        std::vector<std::uint32_t> best_sol_sampling;
        std::vector<std::vector<bool>> predicted_value;
        std::vector<std::vector<float>> ranking_scores;
        std::vector<std::vector<float>> corr_xy;
        float min_cbm, max_rbm;
        std::vector<std::uint32_t> rank;
        std::uint64_t num_edge_left;
        void multi_thread_sampling();
        void constructing_test_data();
        void loading_output_data();

//        // Prints the integer solution obtained by CPLEX to stdout.
//        void print_solution(const IloCplex& cplex, const IloArray<IloNumVarArray>& x) const;
        public:
            // Builds a solver for graph g.
            explicit Reduce(const Graph& g);
            void removing_variables_rbm();
            void removing_variables_cbm();
            void removing_variables_svm();
            void removing_variables_svm_linear();
            void repair();
            double get_objective_value_sampling() const { return best_obj_sampling; }
            bool get_predicted_value(std::uint32_t i, std::uint32_t j) const { return predicted_value[i][j]; }
            float get_cbm_value(std::uint32_t i, std::uint32_t j) const { return corr_xy[i][j]; }
            float get_rbm_value(std::uint32_t i, std::uint32_t j) const { return ranking_scores[i][j]; }
            float get_min_cbm() const { return min_cbm; }
            float get_max_rbm() const { return max_rbm; }
            void compute_reduced_problem_size();
            void compute_correlation_based_measure();
            void compute_ranking_based_measure();
            std::uint64_t get_reduced_problem_size() const { return num_edge_left; }
            std::vector<std::uint32_t> get_best_sol_sampling() const { return best_sol_sampling; }
    };
}

#endif
