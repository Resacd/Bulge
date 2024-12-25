#include "Simulation.h"

Simulation::Simulation(SubDomain& SubDomain,double tau,double pressure,double del_t)
	: SubDomain_(SubDomain)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs_);
    SetRelaxation(tau);
    pressure_ = pressure;
    m_dt = del_t;
    m_dx = SubDomain_.Geometry_.spacing;
    m_omega = 1.0 / m_relaxation;


    Iterate();
}

Simulation::~Simulation()
{
}


void Simulation::Collision()
{

    for (int f = 0; f < SubDomain_.Involved_lat_.size(); f++) {
        if (SubDomain_.Involved_lat_[f]->NodeType_ != NodeType::Halo) {

            int i = SubDomain_.Involved_lat_[f]->x_position;
            int j = SubDomain_.Involved_lat_[f]->y_position;
            int k = SubDomain_.Involved_lat_[f]->z_position;

            SubDomain_.m_lattice[i][j][k]->density = 0.0;

            for (int dir = 0; dir < _nLatNodes; dir++) {

                SubDomain_.m_lattice[i][j][k]->density += SubDomain_.m_lattice[i][j][k]->m_distributions[dir];

            }


#ifdef D3Q19


            bool hasNullNeighbor = false;
            for (int dir = 0; dir < _nLatNodes; ++dir) {
                if (SubDomain_.Involved_lat_[f]->m_neighbours[dir] == nullptr) {
                    hasNullNeighbor = true;
                    break;
                }
            }

            if (SubDomain_.Involved_lat_[f] != nullptr &&
                SubDomain_.Involved_lat_[f]->BoundaryType_ != BoundaryType::Inlet &&
                SubDomain_.Involved_lat_[f]->BoundaryType_ != BoundaryType::Outlet && hasNullNeighbor) {

                SubDomain_.Involved_lat_[f]->velocity[0] = 0.0;
                SubDomain_.Involved_lat_[f]->velocity[1] = 0.0;
                SubDomain_.Involved_lat_[f]->velocity[2] = 0.0;


            }

            else {


                SubDomain_.m_lattice[i][j][k]->velocity[0] =
                    ((SubDomain_.m_lattice[i][j][k]->m_distributions[0] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[6] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[7] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[8] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[9]) -
                        (SubDomain_.m_lattice[i][j][k]->m_distributions[1] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[10] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[11] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[12] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[13])) / SubDomain_.m_lattice[i][j][k]->density;


                SubDomain_.m_lattice[i][j][k]->velocity[1] =
                    ((SubDomain_.m_lattice[i][j][k]->m_distributions[2] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[6] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[10] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[14] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[15]) -
                        (SubDomain_.m_lattice[i][j][k]->m_distributions[3] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[7] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[11] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[16] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[17])) / SubDomain_.m_lattice[i][j][k]->density;


                SubDomain_.m_lattice[i][j][k]->velocity[2] =
                    ((SubDomain_.m_lattice[i][j][k]->m_distributions[4] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[8] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[12] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[14] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[16]) -
                        (SubDomain_.m_lattice[i][j][k]->m_distributions[5] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[9] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[13] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[15] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[17])) / SubDomain_.m_lattice[i][j][k]->density;
#endif


#ifdef D3Q27


                SubDomain_.m_lattice[i][j][k]->velocity[0] =
                    ((SubDomain_.m_lattice[i][j][k]->m_distributions[1] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[7] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[8] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[11] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[13] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[19] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[21] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[23] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[25]) -
                        (SubDomain_.m_lattice[i][j][k]->m_distributions[2] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[9] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[10] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[12] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[14] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[20] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[22] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[24] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[26])) / SubDomain_.m_lattice[i][j][k]->density;




                SubDomain_.m_lattice[i][j][k]->velocity[1] =
                    ((SubDomain_.m_lattice[i][j][k]->m_distributions[3] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[7] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[9] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[15] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[16] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[19] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[21] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[24] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[26]) -
                        (SubDomain_.m_lattice[i][j][k]->m_distributions[4] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[8] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[10] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[17] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[18] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[20] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[22] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[23] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[25])) / SubDomain_.m_lattice[i][j][k]->density;



                SubDomain_.m_lattice[i][j][k]->velocity[1] =
                    ((SubDomain_.m_lattice[i][j][k]->m_distributions[5] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[11] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[12] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[15] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[17] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[19] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[22] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[23] +
                        SubDomain_.m_lattice[i][j][k]->m_distributions[26]) -
                        (SubDomain_.m_lattice[i][j][k]->m_distributions[6] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[13] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[14] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[16] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[18] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[20] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[21] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[24] +
                            SubDomain_.m_lattice[i][j][k]->m_distributions[25])) / SubDomain_.m_lattice[i][j][k]->density;



#endif 
            }

            for (int dir = 0; dir < _nLatNodes; dir++) {

                double du = SubDomain_.m_lattice[i][j][k]->velocity[0] * SubDomain_.m_velSet->GetDirection(dir)[0] +
                    SubDomain_.m_lattice[i][j][k]->velocity[1] * SubDomain_.m_velSet->GetDirection(dir)[1] +
                    SubDomain_.m_lattice[i][j][k]->velocity[2] * SubDomain_.m_velSet->GetDirection(dir)[2];

                double u_sqr = SubDomain_.m_lattice[i][j][k]->velocity[0] * SubDomain_.m_lattice[i][j][k]->velocity[0] +
                    SubDomain_.m_lattice[i][j][k]->velocity[1] * SubDomain_.m_lattice[i][j][k]->velocity[1] +
                    SubDomain_.m_lattice[i][j][k]->velocity[2] * SubDomain_.m_lattice[i][j][k]->velocity[2];

                SubDomain_.m_lattice[i][j][k]->m_eqdistributions[dir] = SubDomain_.m_velSet->GetWeight(dir) *

                    SubDomain_.m_lattice[i][j][k]->density * (1.0 + 3.0 * du + 4.5 * du * du - 1.5 * u_sqr);
            }


            for (int dir = 0; dir < _nLatNodes; dir++) {

                    SubDomain_.m_lattice[i][j][k]->m_newDistributions[dir] = SubDomain_.m_lattice[i][j][k]->m_distributions[dir] +
                        (1.0 / m_relaxation) * (SubDomain_.m_lattice[i][j][k]->m_eqdistributions[dir] - SubDomain_.m_lattice[i][j][k]->m_distributions[dir]);
            }
        }
    }
}

void Simulation::Streaming()
{
    for (int f = 0; f < SubDomain_.Involved_lat_.size(); f++) {

        int i = SubDomain_.Involved_lat_[f]->x_position;
        int j = SubDomain_.Involved_lat_[f]->y_position;
        int k = SubDomain_.Involved_lat_[f]->z_position;
        for (int dir = 0; dir < _nLatNodes; dir++) {

            if (SubDomain_.m_lattice[i][j][k]->m_neighbours[dir] != nullptr   && SubDomain_.m_lattice[i][j][k]->m_distributions[dir] != 0.0) {
                SubDomain_.m_lattice[i][j][k]->Stream(dir);
            }
        }
    }
}




void Simulation::FaceCommunicate(int iter) {

    int tag_xp = 0;
    int tag_xm = 1;
    int tag_yp = 2;
    int tag_ym = 3;
    int tag_zp = 4;
    int tag_zm = 5;

    int mpi_error;
    MPI_Status status;


    if (SubDomain_.my_left_ != MPI_PROC_NULL  && SubDomain_.need_communication_[0] == 1) {
        std::vector<double> send_data;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int j = 0; j < SubDomain_.my_Ny_; j++) {
                if (SubDomain_.m_lattice[1][j][k] == nullptr) {
                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        send_data.push_back(0.0);
                    }
                }
                else {
                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        send_data.push_back(SubDomain_.m_lattice[1][j][k]->m_distributions[dir]);
                    }
                }
            }
        }

        std::vector<double> recv_data(send_data.size());
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.my_left_, tag_xm,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.my_left_, tag_xp,
            MPI_COMM_WORLD, &status);

        if (mpi_error != MPI_SUCCESS) {
            std::cerr << "MPI_Sendrecv failed with error code: " << mpi_error << std::endl;
        }

        int index = 0;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int j = 0; j < SubDomain_.my_Ny_; j++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    SubDomain_.m_lattice[0][j][k]->m_distributions[dir] = recv_data[index];
                    index++;
                }
            }
        }
    }


    if (SubDomain_.my_right_ != MPI_PROC_NULL && SubDomain_.need_communication_[1] == 1) {
        std::vector<double> send_data;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int j = 0; j < SubDomain_.my_Ny_; j++) {
                if (SubDomain_.m_lattice[SubDomain_.my_Nx_ - 2][j][k] == nullptr) {
                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        send_data.push_back(0.0);
                    }
                }
                else {
                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        send_data.push_back(SubDomain_.m_lattice[SubDomain_.my_Nx_ - 2][j][k]->m_distributions[dir]);
                    }
                }
            }
        }


        std::vector<double> recv_data(send_data.size());
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.my_right_, tag_xp,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.my_right_, tag_xm,
            MPI_COMM_WORLD, &status);

        if (mpi_error != MPI_SUCCESS) {
            std::cerr << "MPI_Sendrecv failed with error code: " << mpi_error << std::endl;
        }

        int index = 0;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int j = 0; j < SubDomain_.my_Ny_; j++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    SubDomain_.m_lattice[SubDomain_.my_Nx_ - 1][j][k]->m_distributions[dir] = recv_data[index];
                    index++;
                }
            }
        }

    }


    if (SubDomain_.my_front_ != MPI_PROC_NULL && SubDomain_.need_communication_[2]==1) {
        std::vector<double> send_data;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                if (SubDomain_.m_lattice[i][1][k] == nullptr) {
                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        send_data.push_back(0.0);
                    }
                }
                else {
                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        send_data.push_back(SubDomain_.m_lattice[i][1][k]->m_distributions[dir]);
                    }
                }
            }
        }

        std::vector<double> recv_data(send_data.size(), 0.0);
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.my_front_, tag_ym,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.my_front_, tag_yp,
            MPI_COMM_WORLD, &status);

        if (mpi_error != MPI_SUCCESS) {
            std::cerr << "MPI_Sendrecv failed with error code: " << mpi_error << std::endl;
        }

        int index = 0;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    SubDomain_.m_lattice[i][0][k]->m_distributions[dir] = recv_data[index];
                    index++;
                }
            }
        }
    }


    if (SubDomain_.my_back_ != MPI_PROC_NULL && SubDomain_.need_communication_[3] == 1) {
        std::vector<double> send_data;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {

                if (SubDomain_.m_lattice[i][SubDomain_.my_Ny_ - 2][k] == nullptr) {
                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        send_data.push_back(0.0);
                    }
                }
                else {
                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        send_data.push_back(SubDomain_.m_lattice[i][SubDomain_.my_Ny_ - 2][k]->m_distributions[dir]);
                    }
                }
            }
        }

        std::vector<double> recv_data(send_data.size(), 0.0);
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.my_back_, tag_yp,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.my_back_, tag_ym,
            MPI_COMM_WORLD, &status);
        int index = 0;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    SubDomain_.m_lattice[i][SubDomain_.my_Ny_ - 1][k]->m_distributions[dir] = recv_data[index];
                    index++;
                }
            }
        }
    }


    if (SubDomain_.my_up_ != MPI_PROC_NULL && SubDomain_.need_communication_[4] == 1) {
        std::vector<double> send_data;
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {

                if (SubDomain_.m_lattice[i][j][SubDomain_.my_Nz_ - 2] == nullptr) {
                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        send_data.push_back(0.0);
                    }
                }
                else {
                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        send_data.push_back(SubDomain_.m_lattice[i][j][SubDomain_.my_Nz_ - 2]->m_distributions[dir]);
                    }
                }
            }
        }

        std::vector<double> recv_data(send_data.size(), 0.0);
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.my_up_, tag_zp,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.my_up_, tag_zm,
            MPI_COMM_WORLD, &status);
        if (mpi_error != MPI_SUCCESS) {
            std::cerr << "MPI_Sendrecv failed with error code: " << mpi_error << std::endl;
        }
        int index = 0;
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    SubDomain_.m_lattice[i][j][SubDomain_.my_Nz_ - 1]->m_distributions[dir] = recv_data[index];
                    index++;
                }
            }
        }
    }


    if (SubDomain_.my_down_ != MPI_PROC_NULL && SubDomain_.need_communication_[5] == 1) {
        std::vector<double> send_data;
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {

                if (SubDomain_.m_lattice[i][j][1] == nullptr) {
                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        send_data.push_back(0.0);
                    }
                }

                else {
                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        send_data.push_back(SubDomain_.m_lattice[i][j][1]->m_distributions[dir]);
                    }
                }
            }
        }

        std::vector<double> recv_data(send_data.size(), 0.0);
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.my_down_, tag_zm,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.my_down_, tag_zp,
            MPI_COMM_WORLD, &status);

        if (mpi_error != MPI_SUCCESS) {
            std::cerr << "MPI_Sendrecv failed with error code: " << mpi_error << std::endl;
        }
        int index = 0;
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int i = 0; i < SubDomain_.my_Nx_; i++) {
                for (int dir = 0; dir < _nLatNodes; dir++) {
                    SubDomain_.m_lattice[i][j][0]->m_distributions[dir] = recv_data[index];
                    index++;
                }
            }
        }
    }
}


void Simulation::EdgeCommunicate() {

    int tag_rup = 6;   int tag_ldwn = 7;  int tag_rdwn = 8;  int tag_lup = 9;   int tag_rfrt = 10; int tag_lbck = 11;
    int tag_rbck = 12; int tag_lfrt = 13; int tag_ufrt = 14; int tag_dbck = 15; int tag_ubck = 16; int tag_dfrt = 17;

    int mpi_error;
    MPI_Status status;


    if (SubDomain_.right_up_ != -1 && SubDomain_.need_communication_[6] == 1) {
        std::vector<double> send_data;
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int dir = 0; dir < _nLatNodes; dir++) {

                if (SubDomain_.m_lattice[SubDomain_.my_Nx_ - 2][j][SubDomain_.my_Nz_ - 2] == nullptr) {
                    send_data.push_back(0.0);
                }
                else {

                    send_data.push_back(SubDomain_.m_lattice[SubDomain_.my_Nx_ - 2][j][SubDomain_.my_Nz_ - 2]->m_distributions[dir]);
                }

            }
        }

        std::vector<double> recv_data(send_data.size());
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.right_up_, tag_rup,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.right_up_, tag_ldwn,
            MPI_COMM_WORLD, &status);

        if (mpi_error != MPI_SUCCESS) {
            std::cerr << "MPI_Sendrecv failed with error code: " << mpi_error << std::endl;
        }

        int index = 0;
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int dir = 0; dir < _nLatNodes; dir++) {
                SubDomain_.m_lattice[SubDomain_.my_Nx_ - 1][j][SubDomain_.my_Nz_ - 1]->m_distributions[dir] = recv_data[index];
                index++;
            }
        }
    }

    if (SubDomain_.right_down_ != -1 && SubDomain_.need_communication_[7] == 1) {

        std::vector<double> send_data;
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int dir = 0; dir < _nLatNodes; dir++) {

                if (SubDomain_.m_lattice[SubDomain_.my_Nx_ - 2][j][1] == nullptr) {

                    send_data.push_back(0.0);
                }

                else {

                    send_data.push_back(SubDomain_.m_lattice[SubDomain_.my_Nx_ - 2][j][1]->m_distributions[dir]);
                }
            }
        }

        std::vector<double> recv_data(send_data.size(), 0.0);
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.right_down_, tag_rdwn,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.right_down_, tag_lup,
            MPI_COMM_WORLD, &status);

        if (mpi_error != MPI_SUCCESS) {
            std::cerr << "MPI_Sendrecv failed with error code: " << mpi_error << std::endl;
        }
        int index = 0;
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int dir = 0; dir < _nLatNodes; dir++) {
                SubDomain_.m_lattice[SubDomain_.my_Nx_ - 1][j][0]->m_distributions[dir] = recv_data[index];
                index++;
            }
        }
    }

    if (SubDomain_.left_up_ != -1 && SubDomain_.need_communication_[8]==1) {

        std::vector<double> send_data;
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int dir = 0; dir < _nLatNodes; dir++) {

                if (SubDomain_.m_lattice[1][j][SubDomain_.my_Nz_ - 2] == nullptr) {

                    send_data.push_back(0.0);
                }

                else {

                    send_data.push_back(SubDomain_.m_lattice[1][j][SubDomain_.my_Nz_ - 2]->m_distributions[dir]);
                }
            }
        }

        std::vector<double> recv_data(send_data.size(), 0.0);
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.left_up_, tag_lup,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.left_up_, tag_rdwn,
            MPI_COMM_WORLD, &status);

        if (mpi_error != MPI_SUCCESS) {
            std::cerr << "MPI_Sendrecv failed with error code: " << mpi_error << std::endl;
        }
        int index = 0;
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int dir = 0; dir < _nLatNodes; dir++) {
                SubDomain_.m_lattice[0][j][SubDomain_.my_Nz_ - 1]->m_distributions[dir] = recv_data[index];
                index++;
            }
        }
    }

    if (SubDomain_.left_down_ != -1 && SubDomain_.need_communication_[9]==1) {
        std::vector<double> send_data;
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int dir = 0; dir < _nLatNodes; dir++) {

                if (SubDomain_.m_lattice[1][j][1] == nullptr) {
                    send_data.push_back(0.0);
                }
                else {

                    send_data.push_back(SubDomain_.m_lattice[1][j][1]->m_distributions[dir]);

                }
            }
        }

        std::vector<double> recv_data(send_data.size());
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.left_down_, tag_ldwn,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.left_down_, tag_rup,
            MPI_COMM_WORLD, &status);

        if (mpi_error != MPI_SUCCESS) {
            std::cerr << "MPI_Sendrecv failed with error code: " << mpi_error << std::endl;
        }

        int index = 0;
        for (int j = 0; j < SubDomain_.my_Ny_; j++) {
            for (int dir = 0; dir < _nLatNodes; dir++) {
                SubDomain_.m_lattice[0][j][1]->m_distributions[dir] = recv_data[index];
                index++;
            }
        }
    }

    if (SubDomain_.up_front_ != -1 && SubDomain_.need_communication_[10]==1) {

        std::vector<double> send_data;
        for (int i = 0; i < SubDomain_.my_Nx_; i++) {
            for (int dir = 0; dir < _nLatNodes; dir++) {

                if (SubDomain_.m_lattice[i][1][SubDomain_.my_Nz_ - 2] == nullptr) {

                    send_data.push_back(0.0);
                }

                else {

                    send_data.push_back(SubDomain_.m_lattice[i][1][SubDomain_.my_Nz_ - 2]->m_distributions[dir]);
                }
            }
        }

        std::vector<double> recv_data(send_data.size(), 0.0);
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.up_front_, tag_ufrt,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.up_front_, tag_dbck,
            MPI_COMM_WORLD, &status);

        if (mpi_error != MPI_SUCCESS) {
            std::cerr << "MPI_Sendrecv failed with error code: " << mpi_error << std::endl;
        }
        int index = 0;
        for (int i = 0; i < SubDomain_.my_Nx_; i++) {
            for (int dir = 0; dir < _nLatNodes; dir++) {
                SubDomain_.m_lattice[i][0][SubDomain_.my_Nz_ - 1]->m_distributions[dir] = recv_data[index];
                index++;
            }
        }
    }

    if (SubDomain_.up_back_ != -1 && SubDomain_.need_communication_[11]==1) {

        std::vector<double> send_data;
        for (int i = 0; i < SubDomain_.my_Nx_; i++) {
            for (int dir = 0; dir < _nLatNodes; dir++) {

                if (SubDomain_.m_lattice[i][SubDomain_.my_Ny_ - 2][SubDomain_.my_Nz_ - 2] == nullptr) {

                    send_data.push_back(0.0);
                }

                else {

                    send_data.push_back(SubDomain_.m_lattice[i][SubDomain_.my_Ny_ - 2][SubDomain_.my_Nz_ - 2]->m_distributions[dir]);
                }
            }
        }

        std::vector<double> recv_data(send_data.size(), 0.0);
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.up_back_, tag_ubck,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.up_back_, tag_dfrt,
            MPI_COMM_WORLD, &status);

        if (mpi_error != MPI_SUCCESS) {
            std::cerr << "MPI_Sendrecv failed with error code: " << mpi_error << std::endl;
        }
        int index = 0;
        for (int i = 0; i < SubDomain_.my_Nx_; i++) {
            for (int dir = 0; dir < _nLatNodes; dir++) {
                SubDomain_.m_lattice[i][SubDomain_.my_Ny_ - 2][SubDomain_.my_Nz_ - 2]->m_distributions[dir] = recv_data[index];
                index++;
            }
        }
    }

    if (SubDomain_.down_front_ != -1 && SubDomain_.need_communication_[12]==1) {

        std::vector<double> send_data;
        for (int i = 0; i < SubDomain_.my_Nx_; i++) {
            for (int dir = 0; dir < _nLatNodes; dir++) {

                if (SubDomain_.m_lattice[i][1][1] == nullptr) {

                    send_data.push_back(0.0);
                }

                else {

                    send_data.push_back(SubDomain_.m_lattice[i][1][1]->m_distributions[dir]);
                }
            }
        }

        std::vector<double> recv_data(send_data.size(), 0.0);
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.down_front_, tag_dfrt,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.down_front_, tag_ubck,
            MPI_COMM_WORLD, &status);

        if (mpi_error != MPI_SUCCESS) {
            std::cerr << "MPI_Sendrecv failed with error code: " << mpi_error << std::endl;
        }
        int index = 0;
        for (int i = 0; i < SubDomain_.my_Nx_; i++) {
            for (int dir = 0; dir < _nLatNodes; dir++) {
                SubDomain_.m_lattice[i][0][0]->m_distributions[dir] = recv_data[index];
                index++;
            }
        }

    }


    if (SubDomain_.down_back_ != -1 && SubDomain_.need_communication_[13]==1) {

        std::vector<double> send_data;
        for (int i = 0; i < SubDomain_.my_Nx_; i++) {
            for (int dir = 0; dir < _nLatNodes; dir++) {

                if (SubDomain_.m_lattice[i][SubDomain_.my_Ny_ - 2][1] == nullptr) {

                    send_data.push_back(0.0);
                }

                else {

                    send_data.push_back(SubDomain_.m_lattice[i][SubDomain_.my_Ny_ - 2][1]->m_distributions[dir]);
                }
            }
        }

        std::vector<double> recv_data(send_data.size(), 0.0);
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.down_back_, tag_dbck,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.down_back_, tag_ufrt,
            MPI_COMM_WORLD, &status);

        if (mpi_error != MPI_SUCCESS) {
            std::cerr << "MPI_Sendrecv failed with error code: " << mpi_error << std::endl;
        }
        int index = 0;
        for (int i = 0; i < SubDomain_.my_Nx_; i++) {
            for (int dir = 0; dir < _nLatNodes; dir++) {
                SubDomain_.m_lattice[i][SubDomain_.my_Ny_ - 1][0]->m_distributions[dir] = recv_data[index];
                index++;
            }
        }

    }
    if (SubDomain_.right_front_ != -1 && SubDomain_.need_communication_[14]==1) {

        std::vector<double> send_data;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int dir = 0; dir < _nLatNodes; dir++) {

                if (SubDomain_.m_lattice[SubDomain_.my_Nx_ - 2][1][k] == nullptr) {

                    send_data.push_back(0.0);
                }
                else {

                    send_data.push_back(SubDomain_.m_lattice[SubDomain_.my_Nx_ - 2][1][k]->m_distributions[dir]);
                }
            }
        }

        std::vector<double> recv_data(send_data.size(), 0.0);
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.right_front_, tag_rfrt,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.right_front_, tag_lbck,
            MPI_COMM_WORLD, &status);

        if (mpi_error != MPI_SUCCESS) {
            std::cerr << "MPI_Sendrecv failed with error code: " << mpi_error << std::endl;
        }

        int index = 0;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int dir = 0; dir < _nLatNodes; dir++) {
                SubDomain_.m_lattice[SubDomain_.my_Nx_ - 1][0][k]->m_distributions[dir] = recv_data[index];
                index++;
            }
        }
    }

    if (SubDomain_.right_back_ != -1 && SubDomain_.need_communication_[15]==1) {

        std::vector<double> send_data;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int dir = 0; dir < _nLatNodes; dir++) {

                if (SubDomain_.m_lattice[SubDomain_.my_Nx_ - 2][SubDomain_.my_Ny_ - 2][k] == nullptr) {

                    send_data.push_back(0.0);
                }

                else {

                    send_data.push_back(SubDomain_.m_lattice[SubDomain_.my_Nx_ - 2][SubDomain_.my_Ny_ - 2][k]->m_distributions[dir]);
                }
            }
        }

        std::vector<double> recv_data(send_data.size(), 0.0);
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.right_back_, tag_rbck,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.right_back_, tag_lfrt,
            MPI_COMM_WORLD, &status);

        if (mpi_error != MPI_SUCCESS) {
            std::cerr << "MPI_Sendrecv failed with error code: " << mpi_error << std::endl;
        }
        int index = 0;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int dir = 0; dir < _nLatNodes; dir++) {
                SubDomain_.m_lattice[SubDomain_.my_Nx_ - 1][SubDomain_.my_Ny_ - 1][k]->m_distributions[dir] = recv_data[index];
                index++;
            }
        }
    }

    if (SubDomain_.left_front_ != -1 && SubDomain_.need_communication_[16]==1) {

        std::vector<double> send_data;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int dir = 0; dir < _nLatNodes; dir++) {

                if (SubDomain_.m_lattice[1][1][k] == nullptr) {

                    send_data.push_back(0.0);
                }

                else {

                    send_data.push_back(SubDomain_.m_lattice[1][1][k]->m_distributions[dir]);
                }
            }
        }

        std::vector<double> recv_data(send_data.size(), 0.0);
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.left_front_, tag_lfrt,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.left_front_, tag_rbck,
            MPI_COMM_WORLD, &status);

        if (mpi_error != MPI_SUCCESS) {
            std::cerr << "MPI_Sendrecv failed with error code: " << mpi_error << std::endl;
        }
        int index = 0;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int dir = 0; dir < _nLatNodes; dir++) {
                SubDomain_.m_lattice[0][0][k]->m_distributions[dir] = recv_data[index];
                index++;
            }
        }

    }


    if (SubDomain_.left_back_ != -1 && SubDomain_.need_communication_[17] == 1) {

        std::vector<double> send_data;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int dir = 0; dir < _nLatNodes; dir++) {

                if (SubDomain_.m_lattice[1][SubDomain_.my_Ny_ - 2][k] == nullptr) {

                    send_data.push_back(0.0);
                }
                else {

                    send_data.push_back(SubDomain_.m_lattice[1][SubDomain_.my_Ny_ - 2][k]->m_distributions[dir]);
                }
            }
        }

        std::vector<double> recv_data(send_data.size(), 0.0);
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, SubDomain_.left_back_, tag_lbck,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, SubDomain_.left_back_, tag_rfrt,
            MPI_COMM_WORLD, &status);

        if (mpi_error != MPI_SUCCESS) {
            std::cerr << "MPI_Sendrecv failed with error code: " << mpi_error << std::endl;
        }

        int index = 0;
        for (int k = 0; k < SubDomain_.my_Nz_; k++) {
            for (int dir = 0; dir < _nLatNodes; dir++) {
                SubDomain_.m_lattice[0][SubDomain_.my_Ny_ - 1][k]->m_distributions[dir] = recv_data[index];
                index++;
            }
        }
    }
}


double Simulation::getInletVelocity(double x, double period) {
    double vel = 0.0;
    double normalizedX = fmod(x, period);



    if (normalizedX >= 0 && normalizedX < (0.1 * period)) {
        vel = Pulsatile_.InletVelequation1(normalizedX);
    }
    else if (normalizedX >= (0.1 * period) && normalizedX < (0.3 * period)) {
        vel = Pulsatile_.InletVelequation2(normalizedX);
    }
    else if (normalizedX >= (0.3 * period) && normalizedX < (0.9 * period)) {
        vel = Pulsatile_.InletVelequation3(normalizedX);
    }

    return vel;
}



void Simulation::Iterate()
{
    double t = 0.0;
    double dt = 1e-6;
    int iter = 0;
    while (t <= 2.7) { 

        Collision();
        FaceCommunicate(iter);
        Streaming();
        //InletVelocity_ = getInletVelocity(t, 0.9);
        InletVelocity_ = 0.1;
        ApplyZouHeBc();
        ApplyWallBc();
        ApplyOutletBc();
       // ComputeVelocityGradientAtWall();
       // ComputeWallShearStress();

        if ((iter/* + 1*/) % 500 == 0) {
            WriteOutput(iter);
        }

        if (myrank_ == 0) {
            std::cout << std::endl;
            std::cout << "Iteration " << " " << iter << " " << "has been performed successfully " << std::endl;
        }

        if (myrank_ == 0) {
            std::cout << std::endl;
            //std::cout << "Time " << " " << t << " " << "has been performed successfully " << std::endl;
        }

        t += dt;
        iter++;
    }
}



void Simulation::ApplyZouHeBc() {


    for (int lat = 0; lat < SubDomain_.InletBoundaryNode_.size(); lat++) {

        SubDomain_.InletBoundaryNode_[lat]->ApplyBc(InletVelocity_);

    }
}



void Simulation::ApplyOutletBc() {
    double rho = 1.01;
    double vx = 0.0;
    double vy = 0.0;
    double vz = 0.0;
    double Nxy = 0.0;
    double Nxz = 0.0;

    for (int lat = 0; lat < SubDomain_.OutletBoundaryNode_.size(); lat++) {
        int x_pos = SubDomain_.OutletBoundaryNode_[lat]->x_position;
        int y_pos = SubDomain_.OutletBoundaryNode_[lat]->y_position;
        int z_pos = SubDomain_.OutletBoundaryNode_[lat]->z_position;

        vx = -1 + (1.0 / rho) * (
            SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[2] +
            SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[3] +
            SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[4] +
            SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[5] +
            SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[14] +
            SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[15] +
            SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[16] +
            SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[17] +
            SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[18] +
            2.0 * (
                SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[0] +
                SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[6] +
                SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[7] +
                SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[8] +
                SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[9]
            )
        );

        Nxy = 0.5 * (
            SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[2] +
            SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[14] +
            SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[15] -
            (
                SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[3] +
                SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[16] +
                SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[17]
            )
        ) - (1.0 / 3.0) * rho * vy;

        Nxz = 0.5 * (
            SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[4] +
            SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[10] +
            SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[14] -
            (
                SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[5] +
                SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[15] +
                SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[17]
            )
        ) - (1.0 / 3.0) * rho * vz;

        SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[1] = 
            SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[0] - 
            (1.0 / 3.0) * rho * vx;

        SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[10] = 
            SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[7] +
            (1.0 / 6.0) * rho * (-vx + vy) - Nxy;

        SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[11] = 
            SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[6] +
            (1.0 / 6.0) * rho * (-vx - vy) + Nxy;

        SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[13] = 
            SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[8] +
            (1.0 / 6.0) * rho * (-vx - vz) + Nxz;

        SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[12] = 
            SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_distributions[9] +
            (1.0 / 6.0) * rho * (-vx + vz) - Nxz;
    }
}





void Simulation::ApplyWallBc() {


    for (int lat = 0; lat < SubDomain_.Involved_lat_.size(); lat++) {

        if (SubDomain_.Involved_lat_[lat] != nullptr &&
            SubDomain_.Involved_lat_[lat]->BoundaryType_ != BoundaryType::Inlet &&
            SubDomain_.Involved_lat_[lat]->BoundaryType_ != BoundaryType::Outlet) {

            for (int dir = 0; dir < _nLatNodes; dir++) {
                if (SubDomain_.Involved_lat_[lat]->m_neighbours[dir] == nullptr) {

                    int opp_dir = SubDomain_.GetVelSet()->OppositeDirection(dir);

                    SubDomain_.Involved_lat_[lat]->m_distributions[opp_dir] =
                        SubDomain_.Involved_lat_[lat]->m_distributions[dir];

                }
            }
        }
    }
}

void Simulation::ComputeVelocityGradientAtWall() {
    for (int lat = 0; lat < SubDomain_.WallBoundaryNode_.size(); lat++) {
        int x_pos = SubDomain_.WallBoundaryNode_[lat]->x_position;
        int y_pos = SubDomain_.WallBoundaryNode_[lat]->y_position;
        int z_pos = SubDomain_.WallBoundaryNode_[lat]->z_position;

        std::vector<double> velocity_gradient = { 0.0, 0.0, 0.0 };
        int neighbor_count = 0;

        for (int dir = 0; dir < _nLatNodes; dir++) {
            if (SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_neighbours[dir] != nullptr) {
                int nx = SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_neighbours[dir]->x_position;
                int ny = SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_neighbours[dir]->y_position;
                int nz = SubDomain_.m_lattice[x_pos][y_pos][z_pos]->m_neighbours[dir]->z_position;

                velocity_gradient[0] += ((SubDomain_.m_lattice[nx][ny][nz]->velocity[0] - SubDomain_.m_lattice[x_pos][y_pos][z_pos]->velocity[0]*1e-3));
                velocity_gradient[1] += ((SubDomain_.m_lattice[nx][ny][nz]->velocity[1] - SubDomain_.m_lattice[x_pos][y_pos][z_pos]->velocity[1]*1e-3));
                velocity_gradient[2] += ((SubDomain_.m_lattice[nx][ny][nz]->velocity[2] - SubDomain_.m_lattice[x_pos][y_pos][z_pos]->velocity[2]*1e-3));
                neighbor_count++;
            }
        }

        if (neighbor_count > 0) {
            velocity_gradient[0] /= neighbor_count;
            velocity_gradient[1] /= neighbor_count;
            velocity_gradient[2] /= neighbor_count;
        }

        SubDomain_.WallBoundaryNode_[lat]->velocity_gradient = velocity_gradient;
    }
}

void Simulation::ComputeWallShearStress() {
    double mu = 1.0 / 6.0 * (2.0 * m_relaxation - 1.0);

    for (int lat = 0; lat < SubDomain_.WallBoundaryNode_.size(); lat++) {
        auto velocity_gradient = SubDomain_.WallBoundaryNode_[lat]->velocity_gradient;
        double shear_stress[3];

        shear_stress[0] = mu * velocity_gradient[0];
        shear_stress[1] = mu * velocity_gradient[1];
        shear_stress[2] = mu * velocity_gradient[2];

        SubDomain_.WallBoundaryNode_[lat]->wall_shear_stress = sqrt(shear_stress[0] * shear_stress[0] +
            shear_stress[1] * shear_stress[1] +
            shear_stress[2] * shear_stress[2]);
    }
}


void Simulation::SetRelaxation(double tau)
{
	m_relaxation = tau;
}

void Simulation::WriteOutput(int iter)
{
	OutputFile file(iter);
    file.WriteVTK_LocalData(SubDomain_);
    file.WriteWallShearStressVTK(SubDomain_, iter);
}
