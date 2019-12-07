// Copyright 2019 Pinaev Danil
#include <mpi.h>
#include <math.h>
#include<cmath>
#include<stack>
#include <vector>
#include <random>
//#include <ctime>
//#include <algorithm>
#include "../../../modules/task_3/pinaev_d_graham_alg/graham.h"

const double PI = 3.1415;

int Lowest_point(std::vector<point>& points){
    int l = 0;
    for (size_t i = 0; i < points.size(); ++i) {
        if (points[i].y < points[l].y) {
            l = i;
        }
    }
    return l;
}

void Get_Angles(std::vector<double>& ang, std::vector<point>& points) {
    ang.resize(points.size());
    int j = Lowest_point(points);
    for (size_t i = 0; i < points.size(); i++)
        ang[i] = atan2(points[i].y - points[j].y, points[i].x - points[j].x) * 180 / PI;
}

void Sort(std::vector<double> &a, std::vector<point>& p)
{
    size_t n = p.size() < a.size() ? p.size() : a.size(); 
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            if (a[j] > a[i]) {
                std::swap(a[i], a[j]);
                std::swap(p[i], p[j]);
            }
        }
    }
}

double result(const point &p1, const point &p2, const point &p3)
{
    return (p2.x - p1.x)*(p3.y - p1.y) - (p3.x - p1.x)*(p2.y - p1.y);
}

bool leftTurn(const point &p1, const point &p2, const point &p3)
{
    return result(p1, p2, p3) > 0;
}

void Stack(std::vector<double> &ang, std::vector<point>& points)
{
    std::stack<int>S;
    S.push(0);
    S.push(1);

    for (size_t c = 2; c < ang.size(); c++) {
        int a, b;
        do {
            b = S.top();
            S.pop();
            a = S.top();
        } while (!leftTurn(points[a], points[b], points[c]));

        S.push(b);
        S.push(c);
    }
}


//int main()
//{
//    ifstream inFile;
//    inFile.open("point.txt");
// 
//    if (!inFile.is_open())
//    {
//        cout << "Could not open the file! ";
//        system("pause");
//        exit(EXIT_FAILURE);
//    }
// 
//    int N;
//    inFile >> N;
// 
//    point *ps = new point[N];
// 
//    int i = 0;
//    while (i < N && inFile.good())
//    {
//        inFile >> ps[i].x >> ps[i].y;
//        ++i;
//    }
//    vector<double>ps2;
//
//    //core functianality
//    Angle(ps2, ps, N);
//    Sort(ps2, ps, N);
//    Stack(ps2, ps);
// 
//
//    delete[] ps;
//    system("pause");
//    return 0;
//}