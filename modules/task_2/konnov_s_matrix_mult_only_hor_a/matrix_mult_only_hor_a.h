// Copyright 2019 Konnov Sergey
#ifndef MODULES_TASK_1_KONNOV_S_MATRIX_MULT_ONLY_HOR_A_H_
#define MODULES_TASK_1_KONNOV_S_MATRIX_MULT_ONLY_HOR_A_H_
#include <vector>

std::vector<int> matrix_mult_parallel(std::vector<int>& a, std::vector<int>&b, int msize);
std::vector<int> matrix_mult_sequential(const std::vector<int>& a, const std::vector<int>&b, int msize);
std::vector<int>generate_matrix(int size, int time);

#endif  // MODULES_TASK_1_KONNOV_S_MATRIX_MULT_ONLY_HOR_A_H_
