// Copyright 2019 Konnov Sergey
#include "../../../modules/task_1/konnov_s_trapez_integr/trapez_integr.h"
#include <mpi.h>
#include <functional>
#include <iostream>
double GetTrapezIntegrParallel(int l, int r, int n, const std::function<double(double)>& f) {
    const double step = static_cast<double>(r - l) / static_cast<double>(n);
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (size == 1)
        return GetTrapezIntegrSequential(l, r, step, f);

    const int delta = (r - l) / size;
    int left, right;
    if (rank == 0) {
        left = l, right = l + delta - 1;
        for (int i = 1; i < size - 1; i++) {
            int mes[] = { left + delta * i, left + delta * (i + 1) - 1 };
            MPI_Send(&mes, 2, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
        int mes[] = { left + delta * (size - 1), r };
        MPI_Send(&mes, 2, MPI_INT, size - 1, 0, MPI_COMM_WORLD);
    } else {
        int mes[2];
        MPI_Status status;
        MPI_Recv(&mes, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        left = mes[0], right = mes[1];
    }

    double local_integral;
    if (l > r) {
        local_integral = 0;
    } else {
        local_integral = GetTrapezIntegrSequential(left, right, step, f);
        if (l != left)
            local_integral += f(left)/2;
        if (r != right)
            local_integral += f(right)/2;
    }

    double global_integral;
    MPI_Reduce(&local_integral, &global_integral, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    return global_integral;
}

double GetTrapezIntegrSequential(int l, int r, double step, const std::function<double(double)>& f) {
    double integral = f(static_cast<double>(l)) / 2 + f(static_cast<double>(r)) / 2;
    for (double i = static_cast<double>(l) + step; i < static_cast<double>(r); i += step) {
        integral += f(i);
    }
    integral*=step;
    return integral;
}
