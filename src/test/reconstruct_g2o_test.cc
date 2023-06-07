#include "core/tt_test.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

struct VariableData
{
    int id;
    double data[3];
};

struct FactorData
{
    int variable_ids[2];
    double data[9];
};

void ReadG2O(const std::string& data_file, std::vector<VariableData>& variables, std::vector<FactorData>& factors)
{
    std::ifstream fin(data_file);
    std::string line;

    while(std::getline(fin, line)){
        std::stringstream ss(line);
        std::string data_type;
        ss >> data_type;
        if(data_type == "VERTEX_SE2"){
            VariableData variable_data;
            ss >> variable_data.id >> variable_data.data[0] >> variable_data.data[1] >> variable_data.data[2];
            variables.push_back(variable_data);
        }else if(data_type == "EDGE_SE2"){
            FactorData factor_data;
            ss >> factor_data.variable_ids[0] >> factor_data.variable_ids[1] >> factor_data.data[0] >> factor_data.data[1] >> factor_data.data[2] 
                >> factor_data.data[3] >> factor_data.data[4] >> factor_data.data[5] >> factor_data.data[6] >> factor_data.data[7] >> factor_data.data[8];
            factors.push_back(factor_data);
        }else{
            std::cout << "Unhandled type: " << data_type << std::endl;
        }
    }
}

// JUST_RUN_TEST(reconstruct_g2o, test)
TEST(reconstruct_g2o, test)
{
    std::string data_file = "../data/input_INTEL_g2o.g2o";
    std::string reconstruct_g2o_file = "../data/reconstruct_g2o.g2o";

    std::vector<VariableData> variables;
    std::vector<FactorData> factors;
    ReadG2O(data_file, variables, factors);

    std::ofstream fout(reconstruct_g2o_file);

    for(auto& variable : variables){
        fout << "VERTEX_SE2 " << variable.id << " " << variable.data[0] << " " << variable.data[1] << " " << variable.data[2] << std::endl;
        for(auto& factor : factors){
            if(factor.variable_ids[1] == variable.id){
                fout << "EDGE_SE2 " << factor.variable_ids[0] << " " << factor.variable_ids[1] << " " << factor.data[0] << " " << factor.data[1] << " " << factor.data[2] << " " << factor.data[3] << " " << factor.data[4] << " " << factor.data[5] << " " << factor.data[6] << " " << factor.data[7] << " " << factor.data[8] << std::endl;
            }
        }
    }

    fout.close();
}