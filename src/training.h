#ifndef TRAINING_H
#define TRAINING_H
#include "graph.h"
#include "reduce.h"

extern "C" {
#include "svm_train_model.h"
#include "linear_svm_train_model.h"
}


#include <random>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <numeric>      // std::iota
#include <vector>
#include <cstring>
#include <time.h>
#include <sys/time.h>
#include<cmath>
#include<limits>
#include<iomanip>

namespace TSP{
    class Training{
        // The graph on which we are building a training model.
        std::vector<std::string> training_files;
        std::string input_dir;
        std::string train_data_dir = "../train_data/";
        std::string train_model_name = "train_model";
        std::string train_file_name = "train_data.csv";
        int alpha = 10; //penalty of missing classifying positive training instances
        int kernel_type = 2;
        int prob = 0;
        double weight = 1;

        public:
            // Builds a solver for graph g.
            explicit Training(std::vector<std::string> training_files, std::string input_dir);
            void construct_training_set();
            void generate_training_model_svm();
            void generate_training_model_svm_linear();
    };
}

#endif
