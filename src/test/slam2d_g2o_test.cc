#include "core/tt_test.h"

#include "factor_graph/slam/slam2d_g2o.h"

// JUST_RUN_TEST(slam2d_g2o, test)
TEST(slam2d_g2o, test)
{
    std::string data_file = "../data/input_INTEL_g2o.g2o";
    std::string log_file = "../log/result.log";

    Slam2dG2o slam2d_g2o;
    slam2d_g2o.Run(data_file, log_file);
}