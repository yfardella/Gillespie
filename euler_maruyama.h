#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <random>

class EulerMaruyama {
private:
    private:
    int r, s, c, a, b;  // Populations of species A, B, and C
    double j_off, g_off, k_off, k_u, k_b, alpha;  // Reaction rate constants
    double time, dt;  // Current simulation time
    double J_r, J_rs, J_s;
    std::random_device rd;
    std::mt19937 gen;
    std::string string_r = "", string_s = "", string_c = "", string_time = "";

    void savestring(double s_prev, double r_prev, double c_prev, double time_prev){
        string_r += std::to_string(r_prev) + " ";
        string_s += std::to_string(s_prev) + " ";
        string_c += std::to_string(c_prev) + " ";
        string_time += std::to_string(time_prev) + " ";
    }


    // Print current state of the system
    void printState() const {
        std::cout << "Time: " << time << " | R: " << r << " | S: " << s << " | C: " << c << std::endl;
    }

    std::vector<double> rd_vector(int size){
        std::vector<double> rd_vec(size,0.);
        std::normal_distribution dis{0., 1.0};
        
        for(int k=0; k < size; k++){
            rd_vec[k] = dis(gen);
        }
        return rd_vec;
    }

public:
    // Constructor initializes the time step, start time, and initial conditions
    EulerMaruyama(double s_eq, int r_init, int s_init, int c_init, double J_on, double j_off, double G_on,
                       double g_off, double K_on, double k_off, double k_u, double k_b, int a, int b)
        : r(r_init), s(s_init), c(c_init), J_r(J_on * s_eq), j_off(j_off), J_rs(G_on * s_eq), g_off(g_off),
        J_s(K_on * s_eq), k_off(k_off), k_u(k_u), k_b(k_b / s_eq), a(a), b(b), alpha(double(a)/double(b)), time(0.0) {
        gen.seed(rd());
    }


    // Function to advance one step using Euler-Maruyama method
    void step() {
        // Sampling noise from the normal distribution
        std::normal_distribution noise{0., std::sqrt(dt)};
        double r_prev = r, s_prev = s, c_prev = c;
        // Update based on the ODE functions
        r += (J_r + J_rs + k_u * c - j_off * r - g_off * r - alpha * (k_b * r * s)) * dt + noise(gen);  // Euler update for r
        s += (J_rs + J_s + k_u * c - g_off * r - k_off * s - (k_b * r * s)) * dt + noise(gen);  // Euler update for s
        c += ((k_b * r * s) - k_u * c) * dt;  // Euler update for c

        // Increment time
        time += dt;
    }

    // Run the solver for a number of steps
    void run(double end_time, double dt) {
        savestring(s, r, c, time);
        while (time < end_time) {
            step();
            savestring(s, r, c, time);
            printState();
        }
    }

    std::string receptors(){
        return string_r;
    }
    std::string scaffolds(){
        return string_s;
    }
    std::string complex(){
        return string_c;
    }
    std::string time_str(){
        return string_time;
    }
};