#include <iostream>
#include <vector>
#include <cmath>
#include <functional>

class EulerMaruyama {
private:
    double dt;  // Time step
    double t;   // Current time
    std::vector<double> state; // {x1, x2, x3}

public:
    // Constructor initializes the time step, start time, and initial conditions
    EulerMaruyama(double dt, double t0, const std::vector<double>& initialState)
        : dt(dt), t(t0), state(initialState) {}

    // Function to advance one step using Euler-Maruyama method
    void step(std::function<double(double, double, double)> f1,
              std::function<double(double, double, double)> f2,
              std::function<double(double, double, double)> f3) {
        
        double x1 = state[0];
        double x2 = state[1];
        double x3 = state[2];

        // Update based on the ODE functions
        state[0] += f1(x1, x2, x3) * dt;  // Euler update for r
        state[1] += f2(x1, x2, x3) * dt;  // Euler update for s
        state[2] += f3(x1, x2, x3) * dt;  // Euler update for c

        // Increment time
        t += dt;
    }

    // Function to get the current state
    std::vector<double> getState() const {
        return state;
    }

    // Function to get the current time
    double getTime() const {
        return t;
    }

    // Run the solver for a number of steps
    void run(int steps,
             std::function<double(double, double, double)> f1,
             std::function<double(double, double, double)> f2,
             std::function<double(double, double, double)> f3) {
        
        for (int i = 0; i < steps; ++i) {
            step(f1, f2, f3);
        }
    }
};

// Example usage
int main() {
    // Initial conditions x1, x2, x3 at time t0 = 0
    std::vector<double> initialState = {1.0, 0.0, 0.0};
    double dt = 0.01;  // Time step
    double t0 = 0.0;   // Initial time

    EulerMaruyama solver(dt, t0, initialState);

    // Define the coupled ODEs
    auto f1 = [](double x1, double x2, double x3) {
        return -x1 + x2 * x3;  // Example ODE for r
    };
    auto f2 = [](double x1, double x2, double x3) {
        return -x2 + x1 * x3;  // Example ODE for s
    };
    auto f3 = [](double x1, double x2, double x3) {
        return -x3 + x1 * x2;  // Example ODE for c
    };

    // Run the solver for 1000 steps
    solver.run(1000, f1, f2, f3);

    // Output the final state
    std::vector<double> finalState = solver.getState();
    std::cout << "Final state after 1000 steps:\n";
    std::cout << "x1 = " << finalState[0] << "\n";
    std::cout << "x2 = " << finalState[1] << "\n";
    std::cout << "x3 = " << finalState[2] << "\n";

    return 0;
}
