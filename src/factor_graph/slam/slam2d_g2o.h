#ifndef __SLAM2D_G2O_H__
#define __SLAM2D_G2O_H__

#include <iostream>
#include <Eigen/Core>
#include<fstream>

#include"factor_graph/core/graph_optimize.h"
#include"read_data.h"

class Slam2dG2o
{
public:
    void Run(const std::string &g2o_path, const std::string &log_file)
    {
        std::vector<Variable*> variables; 
        std::vector<Factor*> factors;
        if (!LoadG2O(g2o_path, variables, factors))
            std::cerr<<"Cant load g2o file"<<std::endl;
        
        m_graph.AddVariables(variables);
        m_graph.AddFactors(factors);

        // Fix the first variable.
        m_graph.GetVariables()[0]->fixed = true;

        Optimize(log_file);
    }
private:
    void Optimize(const std::string &log_file)
    {
        DumpPoses(log_file, m_graph);
        GraphOptimize::Option option = GraphOptimize::Option();
        GraphOptimize graph_optimize = GraphOptimize(option);
        graph_optimize.OptimizeGN(&m_graph);
        DumpPoses(log_file, m_graph);
    }

    void DumpPoses(const std::string &filename, const FactorGraph &graph)
    {
        std::ofstream file(filename);
        if (!file.is_open())
        {
            GRAPH_LOG("Failed to open file: %s", filename.c_str());
            return;
        }

        const std::vector<Variable *> &variables = graph.GetVariables();
        for (int i = 0, count = variables.size(); i < count; ++i)
        {
            Pose2d *p = static_cast<Pose2d *>(variables[i]);
            file << i << " " << p->x() << " " << p->y() << " " << p->yaw_rad() << std::endl;
        }
    }

private:
    FactorGraph m_graph;
};

#endif // __SLAM2D_G2O_H__