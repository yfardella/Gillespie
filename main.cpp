#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>


class GillespieSimulator {
public:
    GillespieSimulator(double s_eq, int r_init, int s_init, int c_init, double J_on, double j_off, double G_on,
                       double g_off, double K_on, double k_off, double k_u, double k_b, int a, int b)
        : r(r_init), s(s_init), c(c_init), J_r(J_on * s_eq), j_off(j_off), J_rs(G_on * s_eq), g_off(g_off),
        J_s(K_on * s_eq), k_off(k_off), k_u(k_u), k_b(k_b / s_eq), a(a), b(b), time(0.0) {
        gen.seed(rd());
    }

    // Perform one step of the Gillespie algorithm
    void step() {
        // Calculate propensities (reaction rates)
        r1 = j_off * r;  // r -> r - 1
        r2 = g_off * r;  // r -> r - 1 & s -> s - 1
        r3 = k_u * c / double(b);  // r -> r + a & s -> s + b & c -> c - b
        r4 = (k_b * r * s) / double(b);  // r -> r - a & s -> s - b & c -> c + b
        r5 = k_off * s;  // s -> s - 1

        // Total propensity
        double total_rate = r1 + r2 + r3 + r4 + r5 + J_r + J_rs + J_s;

        // If no reactions are possible, stop the simulation
        if (total_rate == 0.0) return;

        // Generate two random numbers
        std::uniform_real_distribution<> dis(0.0, 1.0);
        double r1_rand = dis(gen);
        double r2_rand = dis(gen);

        // Time until next reaction (exponentially distributed)
        double tau = -std::log(r1_rand) / total_rate;
        time += tau;

        // Determine which reaction occurs
        double threshold = r2_rand * total_rate;
        if (threshold < r1) { // r -> r - 1
            r--;
        } else if (threshold < r1 + r2) { // r -> r - 1 & s -> s - 1
            r--;
            s--;
        } else if (threshold < r1 + r2 + r3) { // r -> r + a & s -> s + b & c -> c - b
            r += a;
            s += b;
            c -= b;
        } else if (threshold < r1 + r2 + r3 + r4) { // r -> r - a & s -> s - b & c -> c + b
            r -= a;
            s -= b;
            c += b;
        } else if (threshold < r1 + r2 + r3 + r4 + r5) { // s -> s - 1
            s--;
        } else if (threshold < r1 + r2 + r3 + r4 + r5 + J_r) { // r -> r + 1
            r++;
        } else if (threshold < r1 + r2 + r3 + r4 + r5 + J_r + J_rs) { // r -> r + 1 & s -> s + 1
            r++;
            s++;
        } else { // s -> s + 1
            s++;
        }
    }

    // Run simulation for a specified amount of time
    void run(double end_time, double dt) {
        double discretized_time = 5., s_prev, r_prev, c_prev, time_prev;
        int n=0;
        while (time < end_time) {
            s_prev = s; r_prev = r; c_prev = c; time_prev = time;
            step();
            
            if((discretized_time <= time) && (time < dt + discretized_time) && (time_prev < discretized_time)){
                savestring(s_prev, r_prev, c_prev, discretized_time);
                discretized_time += dt;
                n++;
            }
            else if((time >= discretized_time + dt) && (discretized_time <= time) && (time_prev < discretized_time)){
                int m = int((time - discretized_time)/dt);
                std::cout << time << " "<< m<< std::endl;
                while((time >= discretized_time) && (discretized_time < end_time)){
                    savestring(s_prev, r_prev, c_prev, discretized_time);
                    discretized_time += dt;
                    n++;
                }
            }
        }
        std::cout<<n;
        if(n != 1001){
            abort();
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

private:
    int r, s, c, a, b;  // Populations of species A, B, and C
    double j_off, g_off, k_off, k_u, k_b;  // Reaction rate constants
    double time;  // Current simulation time
    double r1, r2, r3, r4, r5, J_r, J_rs, J_s;
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
};

int main() {
    int N_sim = 100000;
    // Initial populations of r, s, c
    double s_eq = 100., dt = 0.01;
    int r_init = 100;
    int s_init = 100;
    int c_init = 0;

    // Reaction rate constants
    double J_on, J_on_hat = 3.0167, j_off = 3.4817, G_on = 0.1, g_off = 5.91, K_on = 0.1, k_off = 3.9019; // k_off, k_u, k_b, j_off, J_on_hat, g_off = [3.9018766517549692, 1.9799379861044861, 6.615504603306438,3.481669158554892, 3.016688368288114, 5.910017828840013]
    double k_u = 1.9799, k_b = 6.6155;
    int  a = 3, b = 2;
    double alpha = double(a)/double(b);

    G_on = J_on_hat/(1. +j_off/g_off);
    J_on = J_on_hat - G_on;

    K_on = k_off/(1.+(k_b/k_u)*J_on_hat/(j_off+g_off)) - G_on + g_off*J_on_hat/(j_off+g_off);

    // Theoretical steady state values
    double r_star = (J_on + G_on) / (j_off + g_off), s_star = (K_on + G_on -g_off*r_star) / k_off, c_star = k_b/k_u * r_star * s_star;


    std::ofstream file_r("./r.txt");
    std::ofstream file_s("./s.txt");
    std::ofstream file_c("./c.txt");
    std::ofstream file_time("./time.txt");

    for(size_t k=0; k<N_sim; k++){

        // Create simulator
        GillespieSimulator simulator(s_eq, r_init, s_init, c_init, J_on, j_off, G_on, g_off, K_on, k_off, k_u, k_b, a, b);

        // Run simulation for 100 units of time
        simulator.run(15.0, dt);

        // Save to file
        file_r << simulator.receptors() << std::endl;
        file_s << simulator.scaffolds() << std::endl;
        file_c << simulator.complex() << std::endl;
        file_time << simulator.time_str() << std::endl;
        std::cout << k << std::endl;
    }

    file_time.close(); file_r.close(); file_s.close(); file_c.close();


    std::cout << "R equilibrium: " << r_star * s_eq << " " << s_star * s_eq << " " << c_star * s_eq << std::endl;
    std::cout  << (alpha * k_b * s_star + j_off + g_off - k_u)/2. - std::sqrt((alpha * k_b * s_star + j_off + g_off - k_u) 
    * (alpha * k_b * s_star + j_off + g_off - k_u)/4. + alpha *k_u * k_b * s_star) << " " << (alpha * k_b * s_star + j_off + g_off - k_u)/2. + std::sqrt((alpha * k_b * s_star + j_off + g_off - k_u) 
    * (alpha * k_b * s_star + j_off + g_off - k_u)/4. + alpha *k_u * k_b * s_star);

    return 0;
}

