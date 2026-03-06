#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "particle.hpp"
#include "verlet.hpp"
#include "observables.hpp"

int main() {
    constexpr double PI = 3.14159265358979323846;
    int N = 10;   // Número de partículas deben ser 10
    double dt = 0.001; // debería ser 0.001
    int steps = 30000; // deben ser 30 mil
    std::vector<Particle2D> particles(N);

    // condiciones iniciales, probablemente no lo mas optimizado posiblemente, pero ahí camella 
    std::vector<double> pos_init_x(N);
    std::vector<double> pos_init_y(N);
    std::vector<double> theta(N);
    for(int i = 0; i < N; i++){
        theta[i]= 2.0 * PI * i / (1.0*N);
    }
    for (int k = 0; k < N; k++){
        pos_init_x[k] = 4.0 * cos(theta[k]);
        pos_init_y[k] = 4.0 * sin(theta[k]);
    }
    for(int i = 0; i < N; i++){
        particles[i].x =pos_init_x[i]; 
        particles[i].y =pos_init_y[i];
        particles[i].vx = 0.0; 
        particles[i].vy = 0.0;
    }


    //guardando los observables
    std::ofstream traj("trayectoria2Dring.xyz"); // Lo mismo que ya teniamos
    std::ofstream obs("observables2Dring.dat"); // Cambio para obtener y pintar los observables

    // encabezado para el archivo de observables
    obs << "# t K U E\n";

    for(int i = 0; i < steps; i++)
    {
        double t = i * dt;

        velocityVerlet2D(particles, dt);

        // Trayectoria partículas
        traj << particles.size() <<"\n";
        traj << "#t = " << t << "\n";

        for(size_t j = 0; j < particles.size(); j++){
             traj << j << " "
                  << particles[j].x << " "
                  << particles[j].y<< " "
                  << 0.0 << "\n";
        }

        // archivo observables
        double K = kineticEnergy2D(particles);
        double U = potentialEnergy2D(particles);
        double E = K + U;

        obs << t << " "
            << K << " "
            << U << " "
            << E << "\n";
    }

    traj.close();
    obs.close();
    return 0;
}