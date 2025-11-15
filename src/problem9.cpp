/*
Problem 9c
Program that uses the special algorithm to solve the matrix equation
Av=g, where A is the tridiagnoal matrix with signature (-1,2,-1).
The solution v and corresponding x is written to a file.
*/

#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include <vector>
#include <tuple>

// --------------- special algorithm ---------------------

std::vector<double> special_algorithm(const std::vector<double>& g)
{
    // empty vectors for btilde, gtilde and v of length n (number of unknowns):
    const int n = g.size();
    std::vector<double> bt(n), gt(n), v(n);

    // initializing first loop
    bt[0] = 2.0;                // btilde_1 = b_1 = 2
    gt[0] = g[0];               // gtilde_1 = g_1

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


// ------- function building the matrix equation (g and x) -------
std::tuple< std::vector<double>, std::vector<double>>
define_matrix_eq_special(int n_steps)
{
    // interior unknowns (exluding endpoints)
    const int n_int = n_steps - 1;

    // step size
    double h = 1.0 / n_steps;

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

    return {g, x_complete};
}


int main()
{
    // defining n-values
    std::vector<int> n_list = {10, 100, 1000, 100000, 10000000};

    // looping and solving the equation for each n
    for (int n_steps : n_list)
    {
        // building RHS g and x_complete
        auto [g, x_complete] = define_matrix_eq_special(n_steps);

        // solving A v = g using special algorithm (without endpoints)
        std::vector<double> v_int = special_algorithm(g);

        // making complete solution v_complete by adding boundary values
        std::vector<double> v_complete(n_steps + 1);
        v_complete[0]       = 0.0;
        v_complete[n_steps] = 0.0;

        for (int i = 0; i < n_steps - 1; i++)
        {
            v_complete[i + 1] = v_int[i];
        }

        // creating and opening a file 
        std::ofstream ofile("../results/output_problem9_n" + std::to_string(n_steps) + ".txt");
        ofile.setf(std::ios::scientific);
        ofile << std::setprecision(10);

        // writing x and v to file (including endpoints)
        for (int i = 0; i <= n_steps; i++)
        {
            ofile << x_complete[i] << " " << v_complete[i] << "\n";
        }

        // closing file
        ofile.close();
    }

    return 0;
}