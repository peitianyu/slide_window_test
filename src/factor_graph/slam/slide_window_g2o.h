#ifndef __SLIDE_WINDOW_G2O_H__
#define __SLIDE_WINDOW_G2O_H__

#include <iostream>
#include <Eigen/Core>
#include<fstream>

#include"factor_graph/core/graph_optimize.h"
#include"read_data.h"

class SlideWindowG2O
{
public:
    SlideWindowG2O(uint max_window_size): m_max_window_size(max_window_size)
    {}

    void Run(const std::string &g2o_path, const std::string &log_file)
    {
        std::map<int, Variable*> variable_map; 
        std::vector<Factor*> factors;
        if (!LoadG2O(g2o_path, variable_map, factors))
            std::cerr<<"Cant load g2o file"<<std::endl;
        
        for(uint i = 0; i < factors.size(); i++)
        {
            if(m_graph.GetVariables().size() < m_max_window_size){
                m_graph.AddVariable(variable_map[factors[i]->GetVariableId(0)]);
                m_graph.AddVariable(variable_map[factors[i]->GetVariableId(1)]);
                m_graph.AddFactor(factors[i]);
            }
            else{
                m_graph.AddVariable(variable_map[factors[i]->GetVariableId(1)]);
                m_graph.AddFactor(factors[i]);
                m_graph.RemoveVariable(m_graph.GetVariables()[0]);
            }
        }
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

        Pose2d *p = static_cast<Pose2d *>(graph.GetVariables().back());
        file << i << " " << p->x() << " " << p->y() << " " << p->yaw_rad() << std::endl;
    }
private:
    uint m_max_window_size;
    FactorGraph m_graph;
};

#endif // __SLAM2D_G2O_H__