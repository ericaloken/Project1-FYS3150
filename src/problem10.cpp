/*
Problem 10
Program that measures the time used by the general algorithm and the 
special algorithm. The code from problem 7 and problem 9 will be pasted
into this file and I will use chrono to measure the time. 
*/

#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include <vector>
#include <tuple>
#include <chrono>

// ---------------------- code from problem 7 ----------------------

std::vector<double> general_algorithm(const std::vector<double>& a,
                                      const std::vector<double>& b,
                                      const std::vector<double>& c,
                                      const std::vector<double>& g)
{
    // empty vectors for btilde, gtilde and v of length n (number of unknowns):
    const int n = a.size();
    std::vector<double> bt(n), gt(n), v(n);

    // initializing
    bt[0] = b[0];
    gt[0] = g[0];

    // forward substitution
    for (int i = 1; i < n; i++)
    {
        double m = a[i] / bt[i - 1];
        bt[i] = b[i] - m * c[i - 1];
        gt[i] = g[i] - m * gt[i - 1];
    }

    // back substitution
    v[n - 1] = gt[n - 1] / bt[n - 1];       // last element v_n
    for (int i = n - 2; i >= 0; i--)
    {
        v[i] = (gt[i] - c[i] * v[i + 1]) / bt[i];
    }

    return v;
}


std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>,std::vector<double>>
define_matrix_eq_general(int n_steps)
{
    // interior unknowns (exluding endpoints)
    const int n_int = n_steps - 1;

    // step size
    double h = 1.0 / n_steps;

    // defining sub, main and super diagonal in A
    std::vector<double> a(n_int, -1.0);
    std::vector<double> b(n_int,  2.0);
    std::vector<double> c(n_int, -1.0);

    // setting a_1 = 0 and c_n = 0 (these are not in the matrix)
    a[0]         = 0.0;
    c[n_int - 1] = 0.0;

    // empty vector g
    std::vector<double> g(n_int);

    // building RHS g
    // g_i = h^2 * f(x_i), f(x) = 100e^-10x (only interior points)
    for (int i = 0; i < n_int; i++)
    {
        double x_i = (i + 1) * h;
        g[i] = h * h * 100.0 * std::exp(-10.0 * x_i);
    }

    // making vector x including endpoints x(0) = 0, x(n) = 1 
    std::vector<double> x_complete(n_steps + 1);
    for (int i = 0; i <= n_steps; i++)
    {
        x_complete[i] = i * h;
    }

    return {a, b, c, g, x_complete};
}




// ---------------------- code from problem 9 ----------------------

std::vector<double> special_algorithm(const std::vector<double>& g)
{
    // empty vectors for btilde, gtilde and v of length n (number of unknowns):
    const int n = g.size();
    std::vector<double> bt(n), gt(n), v(n);

    // initializing
    bt[0] = 2.0;
    gt[0] = g[0];

    // forward substitution
    for (int i = 1; i < n; i++)
    {
        bt[i] = 2.0 - 1.0 / bt[i - 1];
        gt[i] = g[i] + gt[i - 1] / bt[i - 1];
    }

    // back substitution
    v[n - 1] = gt[n - 1] / bt[n - 1];       // last element v_n
    for (int i = n - 2; i >= 0; i--)
    {
        v[i] = (gt[i] + v[i + 1]) / bt[i];
    }

    return v;
}

std::tuple< std::vector<double>, std::vector<double>>
define_matrix_eq_special(int n_steps)
{
    // interior unknowns (exluding endpoints)
    const int n_int = n_steps - 1;

    // step size
    double h = 1.0 / n_steps;
    
    // building RHS g 
    std::vector<double> g(n_int);

    for (int i = 0; i < n_int; i++)
    {
        double x_i = (i + 1) * h;
        g[i] = h * h * 100.0 * std::exp(-10.0 * x_i);
    }

    // making vector x including endpoints x(0) = 0, x(n) = 1 
    std::vector<double> x_complete(n_steps + 1);
    for (int i = 0; i <= n_steps; i++)
    {
        x_complete[i] = i * h;
    }

    return {g, x_complete};
}


// building and solving equation for n_steps with general and special algorithms
void problem7(int n_steps)
{
    auto [a, b, c, g, x_complete] = define_matrix_eq_general(n_steps);
    std::vector<double> v_int_general = general_algorithm(a, b, c, g);
}

void problem9(int n_steps)
{
    auto [g, x_complete] = define_matrix_eq_special(n_steps);
    std::vector<double> v_int_special = special_algorithm(g);
}


int main()
{
    // n-values up to 1e6
    std::vector<int> n_list = {10, 100, 1000, 100000, 1000000};

    // number of timings for each n
    int runs = 100;

    // timing general algorithm
    for (int n_steps : n_list)
    {
        // creating and opening a file 
        std::ofstream ofile("../results/output_problem10_general_n" + std::to_string(n_steps) + ".txt");
        ofile.setf(std::ios::scientific);
        ofile << std::setprecision(10);

        for (int r = 0; r < runs; r++)
        {
            // Start measuring time
            auto t1 = std::chrono::high_resolution_clock::now();

            // perform algorithm
            problem7(n_steps);

            // Stop measuring time
            auto t2 = std::chrono::high_resolution_clock::now();

            // Calculate the elapsed time
            double duration_seconds = std::chrono::duration<double>(t2 - t1).count();
            ofile << duration_seconds << "\n";
        }
        ofile.close();
    }

    // timing special algorithm
    for (int n_steps : n_list)
    {
        // creating and opening a file 
        std::ofstream ofile("../results/output_problem10_special_n" + std::to_string(n_steps) + ".txt");
        ofile.setf(std::ios::scientific);
        ofile << std::setprecision(10);

        for (int r = 0; r < runs; r++)
        {
            // Start measuring time
            auto t1 = std::chrono::high_resolution_clock::now();

            // cperform algorithm
            problem9(n_steps);

            // Stop measuring time
            auto t2 = std::chrono::high_resolution_clock::now();

            // Calculate the elapsed time
            double duration_seconds = std::chrono::duration<double>(t2 - t1).count();
            ofile << duration_seconds << "\n";
        }
        ofile.close();
    }
    return 0;
}