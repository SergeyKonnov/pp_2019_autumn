// Copyright 2019 Konnov Sergey
#include "../../../modules/task_3/konnov_s_monte_carlo_multiple_integration/monte_carlo_multiple_integration.h"
#include <mpi.h>
#include <random>
#include <ctime>
#include <iostream>
#include <vector>

double monteCarloMultipleIntegraion(const std::vector<double>& lower_limits,
                                    const std::vector<double>& upper_limits,
                                    int count_of_dots,
                                    const std::function<double(const std::vector<double>&)>& f,
                                    int seed = -1) {
    if (lower_limits.empty() || upper_limits.empty())
        throw "count of limits must be postive";
    if (lower_limits.size() != upper_limits.size())
        throw "count of lower and upper limits must be equal";
    if (count_of_dots <= 0)
        throw "count of dots must be positive";

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double delta = (upper_limits[0] - lower_limits[0])/(static_cast<double>(size));
    double count_of_dots_proc = count_of_dots/size + (rank < count_of_dots%size?1:0);
    double lower_limit = lower_limits[0] + static_cast<double>(rank)*delta;
    double upper_limit = lower_limit + delta;

    std::mt19937 mt;
    if (seed == -1) {
        mt = std::mt19937(time(0));
    } else {
        mt = std::mt19937(seed);
    }

    std::vector<std::uniform_real_distribution<double>> rand(lower_limits.size());
    rand[0] = std::uniform_real_distribution<double>(lower_limit, upper_limit);
    for (int i = 1; i < lower_limits.size(); i++)
        rand[i] = std::uniform_real_distribution<double>(lower_limits[i], upper_limits[i]);

    double ans = 0.;
    for (int i = 0; i < count_of_dots_proc; i++) {
        std::vector<double> tmp(lower_limits.size());
        for (int i = 0; i < lower_limits.size(); i++)
            tmp[i] = rand[i](mt);
        ans += f(tmp);
    }
    double global_ans = 0;
    ans *= (upper_limit-lower_limit);
    for (int i = 1; i < lower_limits.size(); i++) {
        ans *= (upper_limits[i]-lower_limits[i]);
    }
    ans /= count_of_dots_proc;
    MPI_Reduce(&ans, &global_ans, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    return global_ans;
}

double monteCarloMultipleIntegraionSequentional(const std::vector<double>& lower_limits,
                                                const std::vector<double>& upper_limits,
                                                int count_of_dots,
                                                const std::function<double(const std::vector<double>&)>& f,
                                                int seed = -1) {
    if (lower_limits.empty() || upper_limits.empty())
        throw "count of limits must be postive";
    if (lower_limits.size() != upper_limits.size())
        throw "count of lower and upper limits must be equal";
    if (count_of_dots <= 0)
        throw "count of dots must be positive";
    std::mt19937 mt;
    if (seed == -1) {
        mt = std::mt19937(time(0));
    } else {
        mt = std::mt19937(seed);
    }
    std::vector<std::uniform_real_distribution<double>> rand(lower_limits.size());
    double ans = 0.;
    for (int i = 0; i < lower_limits.size(); i++)
        rand[i] = std::uniform_real_distribution<double>(lower_limits[i], upper_limits[i]);
    for (int i = 0; i < count_of_dots; i++) {
        std::vector<double> tmp(lower_limits.size());
        for (int i = 0; i < lower_limits.size(); i++)
            tmp[i] = rand[i](mt);
        ans += f(tmp);
    }
    for (int i = 0; i < lower_limits.size(); i++) {
        ans *= (upper_limits[i]-lower_limits[i]);
    }
    ans /= count_of_dots;
    return ans;
}
