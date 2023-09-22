#include <cmath>
#include <vector>
#include <memory>
#include <cassert>
#include <fstream>
#include <cstring>
#include <iostream>
#include "cubic_hermite_spline.h"
#include <random>
#include <math.h>



int main(int argc, char* argv[])
{
    // using T = double;
    // std::vector<T> x = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
    // std::vector<T> y(x.size(), 0);
    // for(int i=0; i<x.size(); ++i) y[i] = std::sin(x[i]);
    // MonotoneCubicInterpolation<T> tester(x.data(), y.data(), x.size());
    // std::ofstream file("tester.out");
    // const int nPoints = 100;
    // for(int i=0; i<nPoints; ++i) {
    //     auto x = i / static_cast<T>(nPoints);
    //     file << x << " " << tester(x) << "\n";
    // }
    const int nrand=10; 
    int t0 = 0, dt = 10;
    double x0 = 10, y0 = 30, v0 = 5, f0 = 60, ep = 0.03;
    std::default_random_engine generator;
    std::normal_distribution<double> distx0(x0,ep);
    std::normal_distribution<double> disty0(y0,ep);

    double X0[nrand]={};
    double Y0[nrand]={};

    for (int i=0; i<nrand; ++i) {
        double number1 = distx0(generator);
        double number2 = disty0(generator);
        X0[i] = number1;
        Y0[i] = number2;
    }
    for (int i = 0; i < nrand; i++)
    {
        std::cout<< "x0="<< X0[i] <<" " <<"y0="<< Y0[i]<< std::endl;
    }
    // Calculate the mean
    double sumx0 = 0.0;
    double sumy0 = 0.0;
    for (int i = 0; i < nrand; ++i) {
        sumx0 += X0[i];
        sumy0 += Y0[i];
    }
    double meanx0 = sumx0 / nrand;
    double meany0 = sumy0 / nrand;

    // Calculate the variance
    double sumSquaredDiffx0 = 0.0;
    double sumSquaredDiffy0 = 0.0;
    for (int i = 0; i < nrand; ++i) {
        double diffx0 = X0[i] - meanx0;
        sumSquaredDiffx0 += diffx0 * diffx0;
        double diffy0 = Y0[i] - meany0;
        sumSquaredDiffy0 += diffy0 * diffy0;
    }
    double variancex0 = sumSquaredDiffx0 / (nrand -1);
    double variancey0 = sumSquaredDiffy0 / (nrand -1);

    std::cout<< "mean x0=" << meanx0 <<" " << "var x0=" << variancex0 << std::endl;
    std::cout<< "mean y0=" << meany0 <<" " << "var y0=" << variancey0 << std::endl;

    int t1= t0 + dt;
    double s0 = dt * v0;
    double frad0 = f0 * M_PI / 180;
    std::cout<< "t1=" << t1 << std::endl;
    std::cout<< "s0=" << s0 << std::endl;
    std::cout<< "frad0=" << frad0 << std::endl;
    

    double x1 = meanx0 + s0 * cos(frad0);
    double y1 = meany0 + s0 * sin(frad0);

    std::normal_distribution<double> distx1(x1,ep);
    std::normal_distribution<double> disty1(y1,ep);

    double X1[nrand]={};
    double Y1[nrand]={};

    for (int i=0; i<nrand; ++i) {
        double number1 = distx1(generator);
        double number2 = disty1(generator);
        X1[i] = number1;
        Y1[i] = number2;
    }
    for (int i = 0; i < nrand; i++)
    {
        std::cout<< "x1="<< X1[i] <<" " <<"y1="<< Y1[i]<< std::endl;
    }
    // Calculate the mean
    double sumx1 = 0.0;
    double sumy1 = 0.0;
    for (int i = 0; i < nrand; ++i) {
        sumx1 += X1[i];
        sumy1 += Y1[i];
    }
    double meanx1 = sumx1 / nrand;
    double meany1 = sumy1 / nrand;

    // Calculate the variance
    double sumSquaredDiffx1 = 0.0;
    double sumSquaredDiffy1 = 0.0;
    for (int i = 0; i < nrand; ++i) {
        double diffx1 = X1[i] - meanx1;
        sumSquaredDiffx1 += diffx1 * diffx1;
        double diffy1 = Y1[i] - meany1;
        sumSquaredDiffy1 += diffy1 * diffy1;
    }
    double variancex1 = sumSquaredDiffx1 / (nrand -1);
    double variancey1 = sumSquaredDiffy1 / (nrand -1);

    std::cout<< "mean x1=" << meanx1 <<" " << "var x1=" << variancex1 << std::endl;
    std::cout<< "mean y1=" << meany1 <<" " << "var y1=" << variancey1 << std::endl;


    double S1[nrand]={};
    for (int i=0; i<nrand; ++i) {
        double number = sqrt((X1[i] - X0[i]) * (X1[i] - X0[i]) + (Y1[i] - Y0[i])*(Y1[i] - Y0[i]));
        S1[i] = number;
    }
    for (int i = 0; i < nrand; i++)
    {
        std::cout<< "s1="<< S1[i] << std::endl;
    }
    // Calculate the mean
    double sums1 = 0.0;
    for (int i = 0; i < nrand; ++i) {
        sums1 += S1[i];
    }
    double means1 = sums1 / nrand;

    // Calculate the variance
    double sumSquaredDiffs1 = 0.0;
    for (int i = 0; i < nrand; ++i) {
        double diffs1 = S1[i] - means1;
        sumSquaredDiffs1 += diffs1 * diffs1;
    }
    double variances1 = sumSquaredDiffs1 / (nrand -1);

    std::cout<< "mean s1=" << means1 <<" " << "var s1=" << variances1 << std::endl;

    double V1[nrand]={};
    for (int i=0; i<nrand; ++i) {
        double number = S1[i] / dt;
        V1[i] = number;
    }
    for (int i = 0; i < nrand; i++)
    {
        std::cout<< "v1="<< V1[i] << std::endl;
    }
    // Calculate the mean
    double sumv1 = 0.0;
    for (int i = 0; i < nrand; ++i) {
        sumv1 += V1[i];
    }
    double meanv1 = sumv1 / nrand;

    // Calculate the variance
    double sumSquaredDiffv1 = 0.0;
    for (int i = 0; i < nrand; ++i) {
        double diffv1 = V1[i] - meanv1;
        sumSquaredDiffv1 += diffv1 * diffv1;
    }
    double variancev1 = sumSquaredDiffv1 / (nrand -1);

    std::cout<< "mean v1=" << meanv1 <<" " << "var v1=" << variancev1 << std::endl;

    std::ofstream file1("tester1.out");
    file1<< "-----x0 y0-----" << "\n";
    for(int i=0; i<nrand; ++i) {
        file1 << "x0="<< X0[i] <<" " <<"y0="<< Y0[i] << "\n";
    }
    file1 << "mean x0=" << meanx0 <<" " << "var x0=" << variancex0 << "\n";
    file1 << "mean y0=" << meany0 <<" " << "var y0=" << variancey0 << "\n";
    file1<< "t1=" << t1 << "\n";
    file1<< "s0=" << s0 << "\n";
    file1<< "frad0=" << frad0 << "\n";
    file1<< "-----x1 y1-----" << "\n";
    for(int i=0; i<nrand; ++i) {
        file1 << "x1="<< X1[i] <<" " <<"y1="<< Y1[i] << "\n";
    }
    file1 << "mean x1=" << meanx1 <<" " << "var x1=" << variancex1 << "\n";
    file1 << "mean y1=" << meany1 <<" " << "var y1=" << variancey1 << "\n";
    file1<< "-----s1-----" << "\n";
    for(int i=0; i<nrand; ++i) {
        file1 << "s1="<< S1[i] << "\n";
    }
    file1 << "mean s1=" << means1 <<" " << "var s1=" << variances1 << "\n";
    file1<< "-----v1-----" << "\n";
    for(int i=0; i<nrand; ++i) {
        file1 << "v1="<< V1[i] << "\n";
    }
    file1 << "mean v1=" << meanv1 <<" " << "var v1=" << variancev1 << "\n";

    return 0;    
}


