#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <random>
#include "verlet.hpp"
#include "observables.hpp"


enum class Dimension {D1, D2, D3};      // clase para las dimensiones posibles del experimento (1D, 2D, 3D)
enum class Placement {Sphere, Uniform}; // Clase para las posiciones posibles de las particulas (Esfera al rededor del orignes, deistribucion uniforme)

int main() {
    constexpr double PI = 3.14159265358979323846;
    int N = 750;
    double rho = 0.05;                     // número de partículas
    double v_initial = 1.0;         // velocidad inicial
    
    // === definicion condiciones basicas del sistema ===
    Dimension sim_dim = Dimension::D3;  //escogemos la dimension del sistema poniendo D1, D2, D3
    Placement Placement = Placement::Uniform; //escogemos si las particulas 

    // Switches para elegir condiones del sistema
    bool use_rotation = false;       // rotación 2D o 3D
    bool perturbation = true;        // perturbación en posiciones
    bool periodicB =true;
    bool reflectiveB= false;

    // === calculo dinamico de la caja (L) dependiendo de la dimension, N y rho
    double L = 0.0;
    switch (sim_dim){
        case Dimension::D1:
            L = N / rho;               //densidad lineal: rho = N / L
            break;
        case Dimension::D2:
            L = std::sqrt(N/rho);      //densidad superficial: rho = N / l^2
            break;
        case Dimension::D3;
            L = std::pow(N / rho, 1.0/3.0)  //densidad volumetrica: = N / L^3
            break;
        }

    //asiganmos dimensiones a la caja y determinamos que la esfera ocupa el 80% de ella
    double radius = 0.4 * L;
    double Lx = L, Ly = L, Lz = L;

    //en casos de 1D o 2D reducimos las dimensiones de la caja por coordenada
    if(sim_dim == Dimension::D1){Ly = 1.0; Lz = 1.0;}
    if(sim_dim == Dimension::D2){ Lz = 1.0;}

    //caja centrada en el origen
    double x_min = -Lx/2.0, x_max = Lx/2.0;
    double y_min = -Ly/2.0, y_max = Ly/2.0;
    double z_min = -Lz/2.0, z_max = Lz/2.0;

    std::cout << "Simulando N = " << N << " con densidad rho = " << rho << std::endl;
    std::cout << "Tamano de caja calculado L = " << L << std::endl;

    // Paso temporal
    double dt = 0.001;
    int steps = 30000;

    // Acoplamientos
    std::vector<double> k_harmonic = {0.0, 0.0, 0.0};  // z=0

    // Vector de partículas
    std::vector<Particle3D> particles(N);

    // Perturbación aleatoria
    double noise_amp = 0.2;
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist(-noise_amp, noise_amp);

    // Variables auxiliares 1D
    std::vector<double> pos_init(N);
    double start = -3.5;
    double step  = 1.0;

    for(int i=0; i<N; i++){
        //sintaxis de switches (mas organizado)
        switch(sim_dim){
            case Dimension::D1:{
            // === Inicialización 1D ===
            if(i==0) pos_init[i] = -8.0;
            else pos_init[i] = start + (i-1)*step;

            particles[i].x = pos_init[i];
            particles[i].y = 0.0;
            particles[i].z = 0.0;

            particles[i].vx = 0.0;
            particles[i].vy = 0.0;
            particles[i].vz = 0.0;
            break;
            }

            case Dimension::D2: {
                // === Inicialización 2D ===
                double base_x, base_y;

                do {
                    if(N <= 10){
                        double phi = 2.0 * PI * i / N;
                        base_x = radius * cos(phi);
                        base_y = radius * sin(phi);
                    }
                    else{
                        std::uniform_real_distribution<double> distx(-Lx/2.0, Lx/2.0);
                        std::uniform_real_distribution<double> disty(-Ly/2.0, Ly/2.0);
                        base_x = distx(gen);
                        base_y = disty(gen);
                    }

                    // Aplicamos la perturbación solo en el plano XY
                    if(perturbation){
                        base_x += dist(gen);
                        base_y += dist(gen);
                        // z se mantiene en 0.0 explícitamente por su ausencia
                    }

                // Usamos tooClose para asegurar que la perturbación no superponga partículas
                } while(tooClose(particles, base_x, base_y, 0.0, i, 1.0, periodicB, Lx, Ly, Lz));

                particles[i].x = base_x;
                particles[i].y = base_y;
                particles[i].z = 0.0;

                // Inicialización de velocidades
                if(use_rotation){
                    double r_planar = std::sqrt(base_x*base_x + base_y*base_y);
                    if(r_planar > 1e-12){
                        particles[i].vx = -v_initial * base_y / r_planar;
                        particles[i].vy =  v_initial * base_x / r_planar;
                    }
                    else{
                        particles[i].vx = particles[i].vy = 0.0;
                    }
                    particles[i].vz = 0.0; // Mantenemos vz en 0
                }
                else{
                    double mag = std::sqrt(base_x*base_x + base_y*base_y);
                    if(mag > 1e-12){
                        particles[i].vx = v_initial * base_x / mag;
                        particles[i].vy = v_initial * base_y / mag;
                    }
                    else{
                        particles[i].vx = particles[i].vy = 0.0;
                    }
                    particles[i].vz = 0.0; // Mantenemos vz en 0
                }
                break;
            }       

            case Dimension::D3:{
                // === Inicialización 3D ===    
                double base_x, base_y, base_z;

                do{
                
                    double phi = std::acos(1.0 - 2.0 * (i + 0.5) / (double)N);
                    double theta = std::sqrt(N * PI) * phi;
                
                    base_x = radius * std::sin(phi) * std::cos(theta);
                    base_y = radius * std::sin(phi) * std::sin(theta);
                    base_z = radius * std::cos(phi);
                
                    if(perturbation){
                        base_x += dist(gen);
                        base_y += dist(gen);
                        base_z += dist(gen);
                    }
                
                }while(tooClose(particles,
                                base_x, base_y, base_z,
                                i, 1.0,
                                periodicB, Lx, Ly, Lz));
                
                particles[i].x = base_x;
                particles[i].y = base_y;
                particles[i].z = base_z;
                
                if(use_rotation){
                
                    double r_planar = std::sqrt(base_x*base_x + base_y*base_y);
                
                    if(r_planar > 1e-12){
                        particles[i].vx = -v_initial * base_y / r_planar;
                        particles[i].vy =  v_initial * base_x / r_planar;
                        particles[i].vz = 0.0;
                    }
                    else{
                        particles[i].vx = particles[i].vy = particles[i].vz = 0.0;
                    }
                }
                else{
                
                    double mag = std::sqrt(base_x*base_x + base_y*base_y + base_z*base_z);
                
                    if(mag > 1e-12){
                        particles[i].vx = v_initial * base_x / mag;
                        particles[i].vy = v_initial * base_y / mag;
                        particles[i].vz = v_initial * base_z / mag;
                    }
                    else{
                        particles[i].vx = particles[i].vy = particles[i].vz = 0.0;
                    }
                }
                break;
            }       
        }       
    }

    // === generar nombres automáticamente para los archivos de datos ===
    std::stringstream ss;

    // Dimensionalidad
    switch (sim_dim){
    case Dimension::D1: ss << "1D_"; break;
    case Dimension::D2: ss << "2D_"; break;
    case Dimension::D3: ss << "3D_"; break;
    }

    // Tipo de movimiento
    ss << (use_rotation ? "ROT" : "RAD");
    // Velocidad inicial con 1 decimal
    ss << "_v" << std::fixed << std::setprecision(1) << v_initial;
    // Perturbación
    ss << (perturbation ? "_pert" : "_clean");
    // Fronteras
    ss << (periodicB ? "_periodic" : "_box");
    // Construir nombres de archivo finales
    std::string suffix = ss.str();
    std::string traj_filename = "trayectoria_" + suffix + ".dat";
    std::string obs_filename  = "observables_" + suffix + ".dat";
    // === Ahora traj_filename y obs_filename reflejan correctamente 1D, 2D o 3D ===

    //guardando los observables
    std::ofstream traj(traj_filename); // Lo mismo que ya teniamos
    std::ofstream obs(obs_filename); // Cambio para obtener y pintar los observables

    // encabezado para el archivo de observables
    obs << "# t K U E\n";

    for(int i = 0; i < steps; i++)
    {
        double t = i * dt;

        velocityVerlet3D(particles, dt, k_harmonic, x_min, x_max, y_min, y_max, z_min, z_max, reflectiveB, periodicB, Lx,Ly,Lz);

        // Trayectoria partículas
        traj << particles.size() <<"\n";
        traj << "#t = " << t << "\n";

        for(size_t j = 0; j < particles.size(); j++){
             traj << j << " "
                  << particles[j].x << " "
                  << particles[j].y << " "
                  << particles[j].z << "\n";
        }

        // archivo observables
        double K = kineticEnergy3D(particles);
        double U = potentialEnergy3D(particles, k_harmonic, periodicB,Lx,Ly,Lz);
        double E = K + U;

        obs << t << " " << K << " " << U << " " << E << "\n";
    }

    traj.close();
    obs.close();
    return 0;
}