#include"read_data.h"

#include <map>

static double NormalizeAngle(double theta_rad)
{
    // Normalize the angle to the range [-pi, pi).
    constexpr double kPI = 3.14159265358979323846;
    constexpr double k2PI = 2.0 * kPI;
    return (theta_rad - k2PI * std::floor((theta_rad + kPI) / k2PI));
}

bool LoadG2O(const std::string &filename, std::map<int, Variable*>& variables, std::vector<Factor*>& factors) 
{
    variables.clear();
    factors.clear();

    std::ifstream file(filename);
    if (!file.is_open())
    {
        GRAPH_LOG("Failed to open file: %s", filename.c_str());
        return false;
    }

    std::string line;
    std::map<int, Pose2d *> id_to_pose;
    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        std::string data_type;
        ss >> data_type;
        if (data_type == "VERTEX_SE2")
        {
            int id;
            double x, y, th;
            ss >> id >> x >> y >> th;
            Pose2d *p = new Pose2d(x, y, NormalizeAngle(th));
            variables[id] = p;
            id_to_pose[id] = p;
        }
        else if (data_type == "EDGE_SE2")
        {
            int id_a, id_b;
            double dx, dy, d_yaw, i11, i12, i13, i22, i23, i33;
            ss >> id_a >> id_b >> dx >> dy >> d_yaw >> i11 >> i12 >> i13 >> i22 >> i23 >> i33;
            Eigen::Matrix3d info_mtrx = (Eigen::Matrix3d() << i11, i12, i13,
                                                                i12, i22, i23,
                                                                i13, i23, i33)
                                                                .finished();
            GRAPH_ASSERT(id_to_pose.count(id_a) != 0);
            GRAPH_ASSERT(id_to_pose.count(id_b) != 0);
            factors.push_back(new Pose2PoseFactor(id_to_pose[id_a], id_to_pose[id_b], dx, dy, d_yaw, info_mtrx
                                                               ));
        }
        else
        {
            GRAPH_LOG("Unhandled type: %s", data_type.c_str());
            return false;
        }
    }
    
    return true;
}




