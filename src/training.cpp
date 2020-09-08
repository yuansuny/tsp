#include "training.h"

namespace TSP {
    Training::Training(std::vector<std::string> training_files, std::string input_dir) : training_files{training_files}, input_dir{input_dir}{
        std::cout << "number of training graphs is: " << training_files.size() << "\n";
        construct_training_set();
    }

    void Training::construct_training_set(){
        std::string train_s = train_data_dir  + train_file_name;
        char train_data[train_s.size()+1];
        strcpy(train_data, train_s.c_str());

        std::ofstream train_file(train_data, std::ios::trunc);
        if (! train_file.is_open()){
            std::cout << "Cannot open the output file " <<  train_data << "\n";
            return;
        }
        train_file.close();

        std::uint64_t num0 = 0;
        std::uint64_t num1 = 0;
        for (auto d = 0u; d < training_files.size(); ++d){
            const auto graph = Graph(training_files[d], input_dir);
            auto reduce = Reduce(graph);
            reduce.compute_correlation_based_measure();
            reduce.compute_ranking_based_measure();
            train_file.open(train_data, std::ios::app);
            for (auto i = 0u; i < graph.size(); ++i){
                for (auto j = i+1; j < graph.size(); ++j){
                    if (graph.get_optimal_value(i,j) == 1)
                        num1++;
                    else
                        num0++;
                    train_file << graph.get_optimal_value(i,j) << " ";
                    train_file << "1:" << std::fixed << std::setprecision(6) << reduce.get_cbm_value(i,j) / reduce.get_min_cbm() << " ";
                    train_file << "2:" << std::fixed << std::setprecision(6) << reduce.get_rbm_value(i,j) / reduce.get_max_rbm() << " ";
                    train_file << "3:" << std::fixed << std::setprecision(6) << (graph.cost(i,j) - graph.get_min_cost(i)) / (graph.get_max_cost(i) - graph.get_min_cost(i)) << " ";
                    train_file << "4:" << std::fixed << std::setprecision(6) << (graph.cost(i,j) - graph.get_min_cost(j)) / (graph.get_max_cost(j) - graph.get_min_cost(j)) << " ";
                    train_file << "5:" << std::fixed << std::setprecision(6) << (graph.cost(i,j) - graph.get_mean_cost(i)) / (graph.get_max_cost(i) - graph.get_min_cost(i)) << " ";
                    train_file << "6:" << std::fixed << std::setprecision(6) << (graph.cost(i,j) - graph.get_mean_cost(j)) / (graph.get_max_cost(j) - graph.get_min_cost(j)) << " " <<"\n";
                }
            }
            train_file.close();
        }
        std::cout << "num0 is " << num0 << "; " << "num1 is " << num1 <<  std::endl;
        weight = (double)alpha * num0/num1;
    }

    void Training::generate_training_model_svm(){
        std::string train_s = train_data_dir  + train_file_name;
        std::string model_s = train_data_dir + train_model_name;
        char train_data[train_s.size()+1];
        char model_file[model_s.size()+1];
        strcpy(train_data, train_s.c_str());
        strcpy(model_file, model_s.c_str());
        std::cout << train_data << std::endl;
        std::cout << model_file << std::endl;
        std::cout << weight << std::endl;
        std::cout << kernel_type << std::endl;
        std::cout << prob << std::endl;
        svm_train_model(train_data, model_file, weight, kernel_type, prob);

        const int rem_result = remove(train_data);
        if(rem_result == 0){
            std::cout << "Successfully remove training data file" << std::endl;
        } else {
            std::cout << "No such training data file " << std::endl;
        }
    }

    void Training::generate_training_model_svm_linear(){
        std::string train_s = train_data_dir  + train_file_name;
        std::string model_s = train_data_dir + train_model_name;
        char train_data[train_s.size()+1];
        char model_file[model_s.size()+1];
        strcpy(train_data, train_s.c_str());
        strcpy(model_file, model_s.c_str());
        linear_svm_train_model(train_data, model_file, weight);

        const int rem_result = remove(train_data);
        if(rem_result == 0){
            std::cout << "Successfully remove training data file" << std::endl;
        } else {
            std::cout << "No such training data file " << std::endl;
        }
    }
}
