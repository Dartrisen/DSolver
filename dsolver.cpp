//---------------------------------------------------------
//  Solver v 2.1
//---------------------------------------------------------
//  code produced by UID_0
//  begin at 05/24/2015
//  modified at 01/09/2025
//  all rights reserved
//---------------------------------------------------------
//  function f is rhs of the equation y''= z'= f(x,y,z)
//  function g is rhs of y'= z = g(x,y,z)
//  y'' + 4y = cos(3x)
//---------------------------------------------------------
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <vector>
#include <chrono>

class DSolver {
    public:
        // initialize solver with chosen amount of steps
        DSolver(size_t steps = 100) : num_steps(steps) {
            x.resize(num_steps);
            y.resize(num_steps);
            z.resize(num_steps);
        }

        static double y_analytic(double x) {
            return cos(2 * x) + sin(2 * x) - 0.2 * cos(3 * x);
        }

        // getters
        const std::vector<double>& get_x() const { return x; }
        const std::vector<double>& get_y() const { return y; }
        const std::vector<double>& get_z() const { return z; }

        void solve(double x0, double y0, double z0, double x_end) {
            if (x0 >= x_end) {
                throw std::invalid_argument("Start point must be less than end point");
            }
            const double h = (x_end - x0) / (num_steps - 1);
            x[0] = x0;
            y[0] = y0;
            z[0] = z0;

            // Corrected x array population
            for (size_t i = 1; i < num_steps; ++i) {
                x[i] = x[i - 1] + h;
            }

            // Fixed RK4 loop
            for (size_t i = 0; i < num_steps - 1; ++i) {
                step_rk4(i, h);
            }
        }

        void write(const std::string& filename) const {
            std::ofstream outfile(filename);
            if (!outfile) {
                throw std::runtime_error("Could not open output file");
            }

            outfile.precision(6);
            for (size_t i = 0; i < num_steps; ++i) {
                outfile << std::scientific
                        << x[i] << " "
                        << y[i] << " "
                        << z[i] << "\n";
            }
        }
    private:
        size_t num_steps;
        std::vector<double> x, y, z;

        static double g(double, double, double z) {
            return z;
        }

        static double f(double x, double y, double) {
            return std::cos(3 * x) - 4 * y;
        }

        void step_rk4(size_t i, double h) {
            const double xi = x[i];
            const double yi = y[i];
            const double zi = z[i];

            const double k1 = h * f(xi, yi, zi);
            const double q1 = h * g(xi, yi, zi);
            const double k2 = h * f(xi + h/2, yi + q1/2, zi + k1/2);
            const double q2 = h * g(xi + h/2, yi + q1/2, zi + k1/2);
            const double k3 = h * f(xi + h/2, yi + q2/2, zi + k2/2);
            const double q3 = h * g(xi + h/2, yi + q2/2, zi + k2/2);
            const double k4 = h * f(xi + h, yi + q3, zi + k3);
            const double q4 = h * g(xi + h, yi + q3, zi + k3);

            z[i + 1] = zi + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
            y[i + 1] = yi + (q1 + 2.0 * q2 + 2.0 * q3 + q4) / 6.0;
        }
};

int main() {
    try {
        size_t num = 100;
        DSolver solver(num);

        auto start_time = std::chrono::high_resolution_clock::now(); // Start timing
        solver.solve(0.0, 0.8, 2.0, 1.0);
        auto end_time = std::chrono::high_resolution_clock::now(); // End timing

        std::chrono::duration<double> duration = end_time - start_time;
        std::cout << "Time taken for computation: " << duration.count() << " seconds\n";
        solver.write("solution.dat");

        const auto& x_num = solver.get_x();
        const auto& y_num = solver.get_y();

        double error_sum = 0.0;
        for (size_t i = 0; i < num; ++i) {
            double error = std::fabs(y_num[i] - solver.y_analytic(x_num[i]));
            error_sum += error;
            std::cout << "Error at step " << i << ": " << error << std::endl;
        }
        double avg_error = error_sum / num;
        std::cout << "Average Error: " << avg_error << std::endl;

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
