#include "graph.h"
#include <cmath>
#include <random>
#include <iterator>
#include <algorithm>

namespace TSP {
    Graph::Graph(std::string file_name, std::string input_dir) :  file_name{file_name}, input_dir{input_dir} {
        std::cout << file_name << "\n";

        read_graph();
        read_optimal_solution();
    }

    void Graph::read_optimal_solution(){
        std::string opt_file_name = input_dir + file_name + ".sol";
        std::ifstream opt_file(opt_file_name);
        if (!opt_file){
            std::cout << "optimal solution is not provided \n";
            return;
        }
        std::string line;
        bool READ_Point = 0;
        std::uint32_t node;
        while(!opt_file.eof()) {
            getline(opt_file, line);
//            std::cout << line << "\n";
            if (line == "EOF" || line == "-1") {
                break;
            }
            if (READ_Point){
                std::stringstream stream(line);
                while (stream >> node) {
                    opt_tour.push_back(node - 1);
//                    std::cout << node << ", ";
                }
//                std::cout << "\n";
            }
            if (line == "TOUR_SECTION"){
                READ_Point = 1;
            }
        }
        opt_obj = 0;
        optimal_solution = std::vector<std::vector<bool>>(n_nodes, std::vector<bool>(n_nodes, 0));
        for (auto i = 0u; i < n_nodes - 1; ++i){
            opt_obj += costs[opt_tour[i]][opt_tour[i+1]];
            optimal_solution[opt_tour[i]][opt_tour[i+1]] = 1;
            optimal_solution[opt_tour[i+1]][opt_tour[i]] = 1;
        }
        opt_obj += costs[opt_tour[n_nodes-1]][opt_tour[0]];
        optimal_solution[opt_tour[n_nodes-1]][opt_tour[0]] = 1;
        optimal_solution[opt_tour[0]][opt_tour[n_nodes-1]] = 1;
        std::cout << "optimal objective value is " << opt_obj << "\n";
    }

    void Graph::read_graph(){
        std::string input_file = input_dir  + file_name + ".tsp";
        std::ifstream file(input_file);
        std::string line;
        bool COORD = 0, isMatrix = 0, successRead = 1;
        char delimiter = ':';
        std::vector<double> x_coord, y_coord;
        std::uint32_t node, count = 0;
        double xc, yc, dis;
        std::vector<double> inp_dis;

        while(!file.eof()) {
            getline(file, line);
//            std::cout << line << "\n";
            if (line == "EOF" || line == "DISPLAY_DATA_SECTION") {
                break;
            }
            if (line.find(delimiter) != std::string::npos) {
                std::string keyword = line.substr(0, line.find(delimiter));
                std::string value = line.substr(line.find(delimiter) + 1, line.npos);
//                std::cout << keyword << value << "\n";
                if (!checkKeyword(trim(keyword), trim(value))) {
                    successRead = 0;
                    break;
                }
            }
            if (COORD){
                std::stringstream stream(line);
                stream >> node >> xc >> yc;
                x_coord.push_back(xc);
                y_coord.push_back(yc);
//                std::cout << node << ", " << x_coord[count] << ", " << y_coord[count] << "\n";
                count++;
            }
            if (isMatrix){
                std::stringstream stream(line);
                while (stream >> dis) {
                    inp_dis.push_back(dis);
//                    std::cout << dis << ", ";
                }
//                std::cout << "\n";
            }
            if (line == "NODE_COORD_SECTION"){
                COORD = 1;
            }
            if (line == "EDGE_WEIGHT_SECTION") {
                isMatrix = 1;
            }
        }

        costs = std::vector<std::vector<double>>(n_nodes, std::vector<double>(n_nodes, 999999999));
        if (COORD){
            if(edgeWeightType == "EUC_2D"){
                for(auto i = 0u; i < n_nodes; ++i) {
                    for(auto j = i+1; j < n_nodes; ++j) {
                        costs[i][j] = round(std::sqrt(std::pow(x_coord[i] - x_coord[j], 2.0) + std::pow(y_coord[i] - y_coord[j], 2.0)));
                        costs[j][i] = costs[i][j];
                    }
                }
            }
//            else if (edgeWeightType == "GEO"){
//                for (auto i = 0u; i < n_nodes; ++i){
//                    double PI = 3.141592;
//                    double deg = round(x_coord[i]);
//                    double min = x_coord[i] - deg;
//                    x_coord[i] = PI * ( deg + 5.0 * min / 3.0 ) / 180.0;
//                    deg = round(y_coord[i]);
//                    min = y_coord[i] - deg;
//                    y_coord[i] = PI * ( deg + 5.0 * min / 3.0 ) / 180.0;
////                    x_coord[i] = 3.141592*(round(x_coord[i]) + 5.0*(x_coord[i]-round(x_coord[i]))/3.0) / 180.0;
////                    y_coord[i] = 3.141592*(round(y_coord[i]) + 5.0*(y_coord[i]-round(y_coord[i]))/3.0) / 180.0;
//                }
//                for(auto i = 0u; i < n_nodes; ++i) {
//                    for(auto j = i+1; j < n_nodes; ++j) {
//                        double RRR = 6378.388;
//                        double q1 = cos( y_coord[i] - y_coord[j] );
//                        double q2 = cos( x_coord[i] - x_coord[j] );
//                        double q3 = cos( x_coord[i] + x_coord[j] );
//                        costs[i][j] = (int)( RRR * acos( 0.5 * ( ( 1.0 + q1 ) * q2 - (1.0 - q1) * q3 ) ) + 1.0 );
//                        costs[j][i] = costs[i][j];
//                    }
//                }
//            }
            else if(edgeWeightType == "GEO"){
                for (auto i = 0u; i < n_nodes; ++i){
                    x_coord[i] = 3.141592*((int)x_coord[i] + 5.0*(x_coord[i]-(int)x_coord[i])/3.0) / 180.0;
                    y_coord[i] = 3.141592*((int)y_coord[i] + 5.0*(y_coord[i]-(int)y_coord[i])/3.0) / 180.0;
                }
                for(auto i = 0u; i < n_nodes; ++i) {
                    for(auto j = i+1; j < n_nodes; ++j) {
                        costs[i][j] = (int)(6378.388 * acos(0.5*((1.0+cos(y_coord[i]-y_coord[j]))*cos(x_coord[i]-x_coord[j])-(1.0-cos(y_coord[i]-y_coord[j]))*cos(x_coord[i]+x_coord[j]))) + 1.0);
//                        costs[i][j] = ceil(6378.388 * acos(0.5*((1.0+cos(y_coord[i]-y_coord[j]))*cos(x_coord[i]-x_coord[j])-(1.0-cos(y_coord[i]-y_coord[j]))*cos(x_coord[i]+x_coord[j]))));
                        costs[j][i] = costs[i][j];
                    }
                }
            }
            else if(edgeWeightType == "MAX_2D"){
                for(auto i = 0u; i < n_nodes; ++i) {
                    for(auto j = i+1; j < n_nodes; ++j) {
                        costs[i][j] = round(std::max(fabs(x_coord[i] - x_coord[j]), fabs(y_coord[i] - y_coord[j])));
                        costs[j][i] = costs[i][j];
                    }
                }
            }
            else if(edgeWeightType == "ATT"){
                for(auto i = 0u; i < n_nodes; ++i) {
                    for(auto j = i+1; j < n_nodes; ++j) {
                        costs[i][j] = ceil(sqrt((std::pow(x_coord[i] - x_coord[j], 2.0) + std::pow(y_coord[i] - y_coord[j], 2.0)) / 10.0));
                        costs[j][i] = costs[i][j];
                    }
                }
            }
            else if(edgeWeightType == "CEIL_2D"){
                for(auto i = 0u; i < n_nodes; ++i) {
                    for(auto j = i+1; j < n_nodes; ++j) {
                        costs[i][j] = ceil(std::sqrt(std::pow(x_coord[i] - x_coord[j], 2.0) + std::pow(y_coord[i] - y_coord[j], 2.0)));
                        costs[j][i] = costs[i][j];
                    }
                }
            }
            else {
                std::cout << "Error : Unknown EdgeWeightType " << edgeWeightType << std::endl;
                return;
            }
        }
        if (isMatrix){
            if (edgeWeightFormat == "FULL_MATRIX") {
                for (auto i = 0u; i < n_nodes; ++i) {
                    for (auto j = 0u; j < n_nodes; ++j) {
                        if (i != j) {
                            costs[i][j] = inp_dis[i * n_nodes + j];
                        }
                    }
                }
            }
            if ((edgeWeightFormat == "UPPER_ROW") or (edgeWeightFormat == "LOWER_COL")) {
                count = 0;
                for (auto i = 0u; i < n_nodes - 1; ++i) {
                    for (auto j = i + 1; j < n_nodes; ++j) {
                        costs[i][j] = inp_dis[count];
                        costs[j][i] = inp_dis[count];
                        count++;
                    }
                }
            }
            if ((edgeWeightFormat == "LOWER_ROW") or (edgeWeightFormat == "UPPER_COL")) {
                count = 0;
                for (auto i = 1u; i < n_nodes; ++i) {
                    for (auto j = 0u; j < i; ++j) {
                        costs[i][j] = inp_dis[count];
                        costs[j][i] = inp_dis[count];
                        count++;
                    }
                }
            }
            if ((edgeWeightFormat == "UPPER_DIAG_ROW") or (edgeWeightFormat == "LOWER_DIAG_COL")) {
                count = 0;
                for (auto i = 0u; i < n_nodes; ++i) {
                    for (auto j = i; j < n_nodes; ++j) {
                        if (i != j) {
                            costs[i][j] = inp_dis[count];
                            costs[j][i] = inp_dis[count];
                        }
                        count++;
                    }
                }
            }
            if ((edgeWeightFormat == "LOWER_DIAG_ROW") or (edgeWeightFormat == "UPPER_DIAG_COL")) {
                count = 0;
                for (auto i = 0u; i < n_nodes; ++i) {
                    for (auto j = 0u; j < i + 1; ++j) {
                        if (j != i) {
                            costs[i][j] = inp_dis[count];
                            costs[j][i] = inp_dis[count];
                        }
                        count++;
                    }
                }
            }
        }
        if (!successRead) {
            std::cout << "No success on read input file. Exiting..." << std::endl;
            exit(0);
        }

        min_costs = std::vector<double>(n_nodes, 999999999);
        max_costs = std::vector<double>(n_nodes, 0);
        mean_costs = std::vector<double>(n_nodes, 0);
        for (auto i = 0u; i < n_nodes; ++i) {
            for (auto j = 0u; j < n_nodes; ++j) {
                if (i != j){
                    if (costs[i][j] < min_costs[i]){
                        min_costs[i] = costs[i][j];
                    }
                    if (costs[i][j] > max_costs[i]){
                        max_costs[i] = costs[i][j];
                    }
                    mean_costs[i] += costs[i][j];
                }
            }
            mean_costs[i] /= (n_nodes - 1);
        }
    }

    // Check and store the value of each keyword
    bool Graph::checkKeyword(std::string keyword, std::string value) {
        if (keyword == "NAME") {
            name = value;
        }
        else if (keyword == "TYPE") {
            if ((value == "TSP") || (value == "ATSP")) {
                type = value;
            }
            else {
                std::cout << keyword << " not supported" << std::endl;
                return 0;
            }
        }
        else if (keyword == "COMMENT") {
            comment = value;
        }
        else if (keyword == "DIMENSION") {
            n_nodes = std::stoi(value);
        }
        else if (keyword == "EDGE_WEIGHT_TYPE") {
//            if (value == "EXPLICIT") {
                edgeWeightType = value;
//            }
//            else {
//                std::cout << keyword << " not supported" << std::endl;
//                return 0;
//            }
        }
        else if (keyword == "EDGE_WEIGHT_FORMAT") {
            if ((value == "FULL_MATRIX")  ||
                    (value == "UPPER_ROW") ||
                    (value == "LOWER_ROW") ||
                    (value == "UPPER_DIAG_ROW") ||
                    (value == "LOWER_DIAG_ROW") ||
                    (value == "UPPER_COL") ||
                    (value == "LOWER_COL") ||
                    (value == "UPPER_DIAG_COL") ||
                    (value == "LOWER_DIAG_COL")) {
                edgeWeightFormat = value;
            }
            else {
                std::cout << keyword << " not supported" << std::endl;
                return 0;
            }
        }
        else if (keyword == "DISPLAY_DATA_TYPE") {
        }
        else if (keyword == "NODE_COORD_TYPE"){
        }
        else {
            std::cout << "Error : Unknown keyword: " << keyword << std::endl;
            return false;
        }
        return true;
    }

    std::ostream& operator<<(std::ostream& out, const Graph& g) {
        for(const auto& row : g.costs) {
            std::copy(row.begin(), row.end(), std::ostream_iterator<double>(out, "\t"));
            out << "\n";
        }
        return out;
    }

    std::string Graph::trim(std::string s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), [](int c){return !std::isspace(c);}).base(), s.end());
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int c){return !std::isspace(c);}));
        return s;
    }
}