#include "core/tt_test.h"

#include <Eigen/Dense>
#include <iostream>


void TestSlideWindow()
{
    // 初始H与b
    Eigen::Matrix<double, 6, 6> H;
    H << 8, 1, 4, 1, 0, 0,
         1, 1, 1, 1, 0, 0,
         4, 1, 8, 1, 4, 1,
         1, 1, 1, 1, 1, 1,
         0, 0, 4, 1, 8, 1,
         0, 0, 1, 1, 1, 1;
    Eigen::Matrix<double, 6, 1> b;
    b << 1, 2, 3, 4, 5, 6;

    std::cout << "H = \n" << H << std::endl;
    std::cout << "b = \n" << b.transpose() << std::endl;

    // 1. 使用和要边缘化掉的量无关的因子，构建剩余变量对应的Hessian矩阵和b
    Eigen::Matrix<double, 4, 4> Hrr_a = H.block<4, 4>(2, 2);
    Eigen::Matrix<double, 4, 1> br_a = b.tail(4);

    std::cout << "-----------------------------" << std::endl;
    std::cout << "Hrr_a = \n" << Hrr_a << std::endl;
    std::cout << "br_a = \n" << br_a.transpose() << std::endl;

    // 2. 挑出和要边缘化掉的量相关的因子，构建待边缘化的Hessian矩阵，并使用舒尔补做边缘化, 生成先验因子
    Eigen::Matrix<double, 2, 2> Hrr_b = Hrr_a.block<2, 2>(0, 0);
    Eigen::Matrix<double, 2, 2> Hrm_b = H.block<2, 2>(2, 0);
    Eigen::Matrix<double, 2, 2> Hmm_b = H.block<2, 2>(0, 0);
    Eigen::Matrix<double, 2, 2> Hmr_b = H.block<2, 2>(0, 2);
    Eigen::Matrix<double, 2, 1> br_b = br_a.head(2);
    std::cout << "-----------------------------" << std::endl;
    std::cout << "H_rr_b = \n" << Hrr_b << std::endl;
    std::cout << "H_rm_b = \n" << Hrm_b << std::endl;
    std::cout << "H_mm_b = \n" << Hmm_b << std::endl;
    std::cout << "H_mr_b = \n" << Hmr_b << std::endl;
    std::cout << "b_r_b = \n" << br_b.transpose() << std::endl;

    Eigen::Matrix<double, 2, 2> H_prior = Hrr_b - Hrm_b * Hmm_b.inverse() * Hmr_b; 
    Eigen::Matrix<double, 2, 1> b_prior = br_b - Hrm_b * Hmm_b.inverse() * b.head(2);

    std::cout << "-----------------------------" << std::endl;
    std::cout << "H_prior = \n" << H_prior << std::endl;
    std::cout << "b_prior = \n" << b_prior.transpose() << std::endl;
    
    // 3. 结合先验因子，构建待边缘化的Hessian矩阵，并使用舒尔补做边缘化, 生成先验因子
    Eigen::Matrix<double, 4, 4> Hrr = Hrr_a;
    Hrr.block<2, 2>(0, 0) += H_prior;

    Eigen::Matrix<double, 4, 1> br = br_a;
    br.head(2) += b_prior;

    std::cout << "-----------------------------" << std::endl;
    std::cout << "Hrr = \n" << Hrr << std::endl;
    std::cout << "br = \n" << br.transpose() << std::endl;


    // 结果对比
    std::cout << "-----------------------------" << std::endl;
    Eigen::Matrix<double, 6, 1> x = H.inverse() * b;
    std::cout << "x= \n " << x.transpose() << std::endl;
}

void TestSlideWindow1()
{
    // 初始H与b
    Eigen::Matrix<double, 6, 6> H;
    H << 13.7143, 1, 4, 1, 0, 0,
         1, 1, 1, 1, 0, 0,
         4, 1, 8, 1, 4, 1,
         1, 1, 1, 1, 1, 1,
         0, 0, 4, 1, 8, 1,
         0, 0, 1, 1, 1, 1;
    Eigen::Matrix<double, 6, 1> b;
    b << 4.42857, 6, 5, 6, 2, 3;

    std::cout << "H = \n" << H << std::endl;
    std::cout << "b = \n" << b.transpose() << std::endl;

    // 1. 使用和要边缘化掉的量无关的因子，构建剩余变量对应的Hessian矩阵和b
    Eigen::Matrix<double, 4, 4> Hrr_a = H.block<4, 4>(2, 2);
    Eigen::Matrix<double, 4, 1> br_a = b.tail(4);

    std::cout << "-----------------------------" << std::endl;
    std::cout << "Hrr_a = \n" << Hrr_a << std::endl;
    std::cout << "br_a = \n" << br_a.transpose() << std::endl;

    // 2. 挑出和要边缘化掉的量相关的因子，构建待边缘化的Hessian矩阵，并使用舒尔补做边缘化, 生成先验因子
    Eigen::Matrix<double, 2, 2> Hrr_b = Hrr_a.block<2, 2>(0, 0);
    Eigen::Matrix<double, 2, 2> Hrm_b = H.block<2, 2>(2, 0);
    Eigen::Matrix<double, 2, 2> Hmm_b = H.block<2, 2>(0, 0);
    Eigen::Matrix<double, 2, 2> Hmr_b = H.block<2, 2>(0, 2);
    Eigen::Matrix<double, 2, 1> br_b = br_a.head(2);
    std::cout << "-----------------------------" << std::endl;
    std::cout << "H_rr_b = \n" << Hrr_b << std::endl;
    std::cout << "H_rm_b = \n" << Hrm_b << std::endl;
    std::cout << "H_mm_b = \n" << Hmm_b << std::endl;
    std::cout << "H_mr_b = \n" << Hmr_b << std::endl;
    std::cout << "b_r_b = \n" << br_b.transpose() << std::endl;

    Eigen::Matrix<double, 2, 2> H_prior = Hrr_b - Hrm_b * Hmm_b.inverse() * Hmr_b; 
    Eigen::Matrix<double, 2, 1> b_prior = br_b - Hrm_b * Hmm_b.inverse() * b.head(2);

    std::cout << "-----------------------------" << std::endl;
    std::cout << "H_prior = \n" << H_prior << std::endl;
    std::cout << "b_prior = \n" << b_prior.transpose() << std::endl;
    
    // 3. 结合先验因子，构建待边缘化的Hessian矩阵，并使用舒尔补做边缘化, 生成先验因子
    Eigen::Matrix<double, 4, 4> Hrr = Hrr_a;
    Hrr.block<2, 2>(0, 0) += H_prior;

    Eigen::Matrix<double, 4, 1> br = br_a;
    br.head(2) += b_prior;

    std::cout << "-----------------------------" << std::endl;
    std::cout << "Hrr = \n" << Hrr << std::endl;
    std::cout << "br = \n" << br.transpose() << std::endl;


    // 结果对比
    std::cout << "-----------------------------" << std::endl;
    Eigen::Matrix<double, 6, 1> x = H.inverse() * b;
    std::cout << "x= \n " << x.transpose() << std::endl;
}

JUST_RUN_TEST(slide_window, test)
TEST(slide_window, test)
{
    TestSlideWindow();
    TestSlideWindow1();
}
