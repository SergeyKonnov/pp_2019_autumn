// Copyright 2019 Konnov Sergey

#ifndef MODULES_TEST_TASKS_TRAPEZ_INTEGR_MPI_H_
#define MODULES_TEST_TASKS_TRAPEZ_INTEGR_MPI_H_
#include <functional>

double GetTrapezIntegrParallel(int l, int r, int n, std::function<double(double)>& f);
double GetTrapezIntegrSequential(int l, int r, double step, std::function<double(double)>& f);

#endif  // MODULES_TEST_TASKS_TRAPEZ_INTEGR_MPI_H_