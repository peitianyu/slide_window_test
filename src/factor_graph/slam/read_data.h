#ifndef __READ_DATA_H__
#define __READ_DATA_H__

#include <iostream>
#include <Eigen/Core>
#include<fstream>

#include"factor_graph/types/univariate_factor.h"
#include"factor_graph/types/bivariate_factor.h"
#include"factor_graph/types/variable_type.h"
#include"factor_graph/core/utils.h"

bool LoadG2O(const std::string &filename, std::vector<Variable*>& variables, std::vector<Factor*>& factors);

#endif // __READ_DATA_H__