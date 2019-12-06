// Copyright 2019 Konnov Sergey
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <mpi.h>
#include <vector>
#include <functional>
#include "./monte_carlo_multiple_integration.h"

#define abs_error 0.5  // ?

TEST(monteCarloIntegraion, One_Dimensional_Function_On_Small_Interval) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::function<double(const std::vector<double>&)>f = 
                    [](const std::vector<double>& v) {return v[0]*v[0];};
    int count_of_dots = 10000;
    double lower_limit = 0., upper_limit = 5.;
    double ans = monteCarloMultipleIntegraion({lower_limit}, {upper_limit}, count_of_dots, f);
    if(rank == 0) {
        std::function<double(double)> f = [](double x) {return x*x*x/3;};
        double ans_check = f(upper_limit) - f(lower_limit);
        ASSERT_NEAR(ans, ans_check, abs_error);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
