#include <iostream>
#include <fstream>
#include "LatticeModel.h"
#include "Simulation.h"
#include "OutputFile.h"
#include <mpi.h>
#include "Subdomain.h"
#include "Subdomain.h"
#include "Input.h"
#include "ComplexGeometry.h"
#include "disclaimer.h"
#include <iostream>
#include <chrono>
#include <thread>
#include <iomanip>
#include <omp.h>
#include "Utility.h"
#include "GeometryReader.h"
#include "Input.h"

int main(int argc, char* argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    int myrank, num_procs;
    Utility Utility;
    Reader Reader;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    double dx = 0.025; // mm = 0.5e-4 m
    double cs = 1.0 / sqrt(3.0); // Speed of sound in lattice units
    double dt = dx * 1e-3 / cs; // Time step in seconds
    double inlet_velocity_physical = 0.01; // mm/s = 0.01 m/s
    double radius = 1.0; // mm = 1e-3 m

    // Calculate inlet velocity in lattice units
    double inlet_velocity_lattice = inlet_velocity_physical * dt / dx;

    // Physical viscosity
    double physical_viscosity = 1e-6;  // physical viscosity in m^2/s

    // Calculate Reynolds number
    double length_scale = 2 * radius * 1e-3; // 2 * radius in meters
    double Re = (inlet_velocity_physical * length_scale) / physical_viscosity;

    // Calculate relaxation time tau
    //double tau = 0.5 + (3 * physical_viscosity * dt) / (dx * dx * 1e-6);
    double tau = 0.67;

    // Calculate lattice viscosity and collision frequency omega
    double viscosity_lattice = (tau - 0.5) * dx * dx / (3 * dt);
    double omega = 1.0 / tau;

    if (myrank == 0) {
        Print_disclaimer();




        std::cout << "dx: " << dx << " mm" << std::endl;
        std::cout << "cs: " << cs << std::endl;
        std::cout << "dt: " << dt << " s" << std::endl;
        std::cout << "Reynolds number (Re): " << Re << std::endl;
        std::cout << "tau: " << tau << std::endl;
        std::cout << "inlet_velocity_lattice: " << inlet_velocity_lattice << std::endl;
        std::cout << "omega: " << omega << std::endl;

        omp_set_num_threads(omp_get_max_threads());
        int numThreads = omp_get_max_threads();
        omp_set_num_threads(numThreads);
        bool stopChrono = false;
        std::thread chronoThread([&]() {
            while (!stopChrono) {
                auto now = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = now - start;

                //std::cout << "\rElapsed time: " << std::fixed << std::setprecision(2)
                 //   << std::setw(8) << elapsed.count() << " seconds" << std::flush;

                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            }
            });
        chronoThread.detach();
    }

    std::shared_ptr<LatticeModel> velSet = std::make_shared<LatticeModel>();
    Geometry Geometry("Bulge.stl", velSet, dx, mdepth);
    SubDomain SubDomain(Geometry);
    Simulation(SubDomain, tau, pressure,dt);

    MPI_Finalize();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    if (myrank == 0) {
        //std::cout << "\rElapsed time: " << std::fixed << std::setprecision(2)
          //  << std::setw(8) << elapsed.count() << " seconds." << std::endl;
    }

    return 0;
}


