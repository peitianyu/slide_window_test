#include "core/tt_test.h"

#include "factor_graph/core/factor_graph.h"
#include "factor_graph/types/pose_between_factor.h"


#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>

class TestSlam2dG2O
{
public:
    TestSlam2dG2O() = default;

    void Run(const std::string &data_path)
    {
        std::vector<Pose2d *> vars;
        std::vector<PoseBetweenFactor *> factors;
        LoadG2O(data_path, vars, factors);

        for(Pose2d *var : vars) m_factor_graph.AddVariable(var);
        
        uint cnt = 0;
        for(PoseBetweenFactor *factor : factors){
            if(cnt++ > 100) break;
            m_factor_graph.AddFactor(factor);
        }
        
        // 只优化前5个因子
        // m_factor_graph.OptimizeGN(std::vector<Factor *>(factors.begin() + 3, factors.begin() + 7));

        m_factor_graph.OptimizeGN();
    }
private:
    void LoadG2O(const std::string &data_path, std::vector<Pose2d *>& vars, std::vector<PoseBetweenFactor *>& factors)
    {
        std::ifstream fin(data_path);
        if(!fin.is_open())
        {
            std::cerr << "Failed to open file: " << data_path << std::endl;
            return;
        }

        std::string line;
        while(std::getline(fin, line))
        {
            std::stringstream ss(line);
            std::string type;
            ss >> type;
            if(type == "VERTEX_SE2")
            {
                int id;
                double x, y, theta;
                ss >> id >> x >> y >> theta;
                Pose2d *pose = new Pose2d(x, y, theta);
                vars.push_back(pose);
            }
            else if(type == "EDGE_SE2")
            {
                int id1, id2;
                double x, y, theta, i11, i12, i13, i22, i23, i33;
                ss >> id1 >> id2 >> x >> y >> theta >> i11 >> i12 >> i13 >> i22 >> i23 >> i33;
               Eigen::Matrix3d info_matrix;
                info_matrix << i11, i12, i13,
                                i12, i22, i23,
                                i13, i23, i33;
                PoseBetweenFactor *factor = new PoseBetweenFactor(vars[id1], vars[id2], x, y, theta, info_matrix);
                factors.push_back(factor);
            }
        }

        fin.close();
    }
private:
    FactorGraph m_factor_graph;
};

JUST_RUN_TEST(factor_graph_construct, test)
TEST(factor_graph_construct, test)
{
    std::string data_path = "../data/input_INTEL_g2o.g2o";
    TestSlam2dG2O test;
    test.Run(data_path);
}