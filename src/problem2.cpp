/* 
Problem 2
This program evaluates u(x) at some points x defined in a vector and saves
the data in a file with two coloumns.
*/

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>

// defining u(x)
double u_exact(double x)
{  
    return 1.0 - (1.0 - std::exp(-10.0)) * x - std::exp(-10.0 * x);
}


int main()
{
    // creating and opening a file
    std::string filename = "../results/output_problem2.txt";
    std::ofstream ofile;
    ofile.open(filename);

    // creating vectors of length n
    const int n = 100000;
    std::vector<double> uvals(n);
    std::vector<double> xvals(n);

    // step size (x1-x0 )/ (n-1)
    double h = 1.0 / (n-1);

    // computing u(x) for each x-value and writing to file
    ofile << std::scientific << std::setprecision(10);
    for (int i = 0; i < n; i++)
    {
        xvals[i] = i * h;
        uvals[i] = u_exact(xvals[i]);
        ofile << xvals[i] << "  " << uvals[i]  << "\n";
    }

    // closing file
    ofile.close();

    return 0;
}