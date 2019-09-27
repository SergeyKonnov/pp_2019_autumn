// Copyright 2019 Konnov Sergey
#include <mpi.h>
#include "../../../modules/task_1/konnov_s_trapez_integr/trapez_integr.h"
#include <functional>
#include <iostream>
double GetTrapezIntegrParallel(int l, int r, int n, std::function<double(double)>& f) {
    
    const double step = double(r - l) / double(n);
    int size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const int delta = (r - l) / rank;
    int left, right;
    
    if (rank == 0) {
        left = l, right = l + delta;
        for (int i = 1; i < size - 1; i++) {
            int mes[] = { left + delta * i, left + delta * (i + 1) };
            MPI_Send(&mes, 2, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
        int mes[] = { left + delta * (size - 1), r };
        MPI_Send(&mes, 2, MPI_INT, size - 1, 0, MPI_COMM_WORLD);
    }
    else {
        int mes[2];
        MPI_Status status;
        MPI_Recv(&mes, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        left = mes[0], right = mes[1];
    }
    
    double local_integral = GetTrapezIntegrSequential(left, right, step, f);
    if (rank != 0)
        local_integral += f(left)/2;
    if (rank != size - 1)
        local_integral += f(right)/2;
    double global_integral = MPI_Reduce(&local_integral, &global_integral, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    return global_integral;
}

double GetTrapezIntegrSequential(int l, int r, double step, std::function<double(double)>& f) {
    double integral = f(double(l)) / 2 + f(double(r)) / 2;
    double i = double(l) + step;
    for (; i < double(r); i+=step) {
        integral += f(i);
    }
    integral*=step;
    return integral;
}