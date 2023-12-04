#include <cstdlib>
#include <thread>

#include "../../lib/taskflow/taskflow.hpp"
#include "rhs.h"

using namespace std;
using namespace nlsm;

namespace tfdendro {

char *executorThreads = std::getenv("DENDRO_TF_EXECUTOR_THREADS");
int threadsUse = executorThreads == NULL ? std::thread::hardware_concurrency()
                                         : atoi(executorThreads);

tf::Executor executor(threadsUse);

tf::Taskflow nlsmTasks;

dendrotf_rhs_data data;
}  // namespace tfdendro

void make_deriv_computation_taskflow(tf::Taskflow &tf) {
    auto chi_xx_deriv = tf.placeholder().name("chi xx deriv");
    auto chi_yy_deriv = tf.placeholder().name("chi yy deriv");
    auto chi_zz_deriv = tf.placeholder().name("chi zz deriv");

    chi_xx_deriv.work([chi_xx_deriv]() {
        std::cout << "in chi_xx_deriv, data pointer is " << chi_xx_deriv.data() << std::endl;
        auto d =
            *static_cast<tfdendro::dendrotf_rhs_data *>(chi_xx_deriv.data());

        unsigned int sz[3] = {d.nx, d.ny, d.nz};
        
        deriv_xx(d.grad2_0_0_chi, d.chi, d.hx, sz, d.bflag);
        std::cout << "finished chi_xx_deriv" << std::endl;
    });
    chi_yy_deriv.work([chi_yy_deriv]() {
        std::cout << "in chi_yy_deriv" << std::endl;
        auto d =
            *static_cast<tfdendro::dendrotf_rhs_data *>(chi_yy_deriv.data());

        unsigned int sz[3] = {d.nx, d.ny, d.nz};
        deriv_yy(d.grad2_1_1_chi, d.chi, d.hy, sz, d.bflag);
        std::cout << "finished chi_yy_deriv" << std::endl;
    });
    chi_zz_deriv.work([chi_zz_deriv]() {
        std::cout << "in chi_zz_deriv" << std::endl;
        auto d =
            *static_cast<tfdendro::dendrotf_rhs_data *>(chi_zz_deriv.data());

        unsigned int sz[3] = {d.nx, d.ny, d.nz};
        deriv_zz(d.grad2_2_2_chi, d.chi, d.hz, sz, d.bflag);
        std::cout << "finished chi_zz_deriv" << std::endl;
    });
}

void make_rhs_computation_taskflow(tf::Taskflow &tf) {
    auto phi_rhs_task = tf.placeholder().name("Phi RHS");
    phi_rhs_task.work([phi_rhs_task]() {
        auto d =
            *static_cast<tfdendro::dendrotf_rhs_data *>(phi_rhs_task.data());
        DendroRegister double x;
        DendroRegister double y;
        DendroRegister double z;
        DendroRegister unsigned int pp;
        double r;

        for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
            z = d.pmin[2] + k * d.hz;
            for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
                y = d.pmin[1] + j * d.hy;
                for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                    x = d.pmin[0] + i * d.hx;
                    pp = i + d.nx * (j + d.ny * k);
#ifdef NLSM_NONLINEAR
                    r = sqrt(x * x + y * y + z * z);
                    f(r > 1.0e-6) {
                        d.phi_rhs[pp] =
                            NLSM_WAVE_SPEED_X * d.grad2_0_0_chi[pp] +
                            NLSM_WAVE_SPEED_Y * d.grad2_1_1_chi[pp] +
                            NLSM_WAVE_SPEED_Z * d.grad2_2_2_chi[pp] -
                            sin(2 * d.chi[pp]) / pow(r, 2);
                    }
                    else {
                        d.chi_rhs[pp] = 0.0;
                    }
#else
                    d.phi_rhs[pp] = NLSM_WAVE_SPEED_X * d.grad2_0_0_chi[pp] +
                                    NLSM_WAVE_SPEED_Y * d.grad2_1_1_chi[pp] +
                                    NLSM_WAVE_SPEED_Z * d.grad2_2_2_chi[pp];
#endif
                }
            }
        }
    });

    auto chi_rhs_task = tf.placeholder().name("Chi RHS");
    chi_rhs_task.data(&tfdendro::data).work([chi_rhs_task]() {
        auto d =
            *static_cast<tfdendro::dendrotf_rhs_data *>(chi_rhs_task.data());
        DendroRegister double x;
        DendroRegister double y;
        DendroRegister double z;
        DendroRegister unsigned int pp;
        double r;

        for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
            z = d.pmin[2] + k * d.hz;
            for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
                y = d.pmin[1] + j * d.hy;
                for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                    x = d.pmin[0] + i * d.hx;
                    pp = i + d.nx * (j + d.ny * k);

#ifdef NLSM_NONLINEAR
                    r = sqrt(x * x + y * y + z * z);
                    if (r > 1.0e-6) {
                        d.chi_rhs[pp] = d.phi[pp];
                    } else {
                        d.chi_rhs[pp] = 0.0;
                    }
#else
                    d.chi_rhs[pp] = d.phi[pp];
#endif
                }
            }
        }
    });
}

void make_kodiss_computation_taskflow(tf::Taskflow &tf){

}

void buildNLSMRHSGraph(int mpiRank) {
    auto chi_xx_deriv = tfdendro::nlsmTasks.placeholder().name("chi xx deriv");
    auto chi_yy_deriv = tfdendro::nlsmTasks.placeholder().name("chi yy deriv");
    auto chi_zz_deriv = tfdendro::nlsmTasks.placeholder().name("chi zz deriv");

    chi_xx_deriv.data(&tfdendro::data).work([chi_xx_deriv]() {
        std::cout << "in chi_xx_deriv" << std::endl;
        auto d =
            *static_cast<tfdendro::dendrotf_rhs_data *>(chi_xx_deriv.data());
        std::cout << d.grad2_0_0_chi << std::endl;
        deriv_xx(d.grad2_0_0_chi, d.chi, d.hx, d.sz, d.bflag);
        std::cout << "finished chi_xx_deriv" << std::endl;
    });
    chi_yy_deriv.data(&tfdendro::data).work([chi_yy_deriv]() {
        std::cout << "in chi_yy_deriv" << std::endl;
        auto d =
            *static_cast<tfdendro::dendrotf_rhs_data *>(chi_yy_deriv.data());
        deriv_yy(d.grad2_1_1_chi, d.chi, d.hy, d.sz, d.bflag);
        std::cout << "finished chi_yy_deriv" << std::endl;
    });
    chi_zz_deriv.data(&tfdendro::data).work([chi_zz_deriv]() {
        std::cout << "in chi_zz_deriv" << std::endl;
        auto d =
            *static_cast<tfdendro::dendrotf_rhs_data *>(chi_zz_deriv.data());
        deriv_zz(d.grad2_2_2_chi, d.chi, d.hz, d.sz, d.bflag);
        std::cout << "finished chi_zz_deriv" << std::endl;
    });

    auto phi_rhs_task = tfdendro::nlsmTasks.placeholder().name("Phi RHS");
    phi_rhs_task.data(&tfdendro::data).work([phi_rhs_task]() {
        auto d =
            *static_cast<tfdendro::dendrotf_rhs_data *>(phi_rhs_task.data());
        DendroRegister double x;
        DendroRegister double y;
        DendroRegister double z;
        DendroRegister unsigned int pp;
        double r;

        for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
            z = d.pmin[2] + k * d.hz;
            for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
                y = d.pmin[1] + j * d.hy;
                for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                    x = d.pmin[0] + i * d.hx;
                    pp = i + d.nx * (j + d.ny * k);
#ifdef NLSM_NONLINEAR
                    r = sqrt(x * x + y * y + z * z);
                    f(r > 1.0e-6) {
                        d.phi_rhs[pp] =
                            NLSM_WAVE_SPEED_X * d.grad2_0_0_chi[pp] +
                            NLSM_WAVE_SPEED_Y * d.grad2_1_1_chi[pp] +
                            NLSM_WAVE_SPEED_Z * d.grad2_2_2_chi[pp] -
                            sin(2 * d.chi[pp]) / pow(r, 2);
                    }
                    else {
                        d.chi_rhs[pp] = 0.0;
                    }
#else
                    d.phi_rhs[pp] = NLSM_WAVE_SPEED_X * d.grad2_0_0_chi[pp] +
                                    NLSM_WAVE_SPEED_Y * d.grad2_1_1_chi[pp] +
                                    NLSM_WAVE_SPEED_Z * d.grad2_2_2_chi[pp];
#endif
                }
            }
        }
    });

    auto chi_rhs_task = tfdendro::nlsmTasks.placeholder().name("Chi RHS");
    chi_rhs_task.data(&tfdendro::data).work([chi_rhs_task]() {
        auto d =
            *static_cast<tfdendro::dendrotf_rhs_data *>(chi_rhs_task.data());
        DendroRegister double x;
        DendroRegister double y;
        DendroRegister double z;
        DendroRegister unsigned int pp;
        double r;

        for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
            z = d.pmin[2] + k * d.hz;
            for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
                y = d.pmin[1] + j * d.hy;
                for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                    x = d.pmin[0] + i * d.hx;
                    pp = i + d.nx * (j + d.ny * k);

#ifdef NLSM_NONLINEAR
                    r = sqrt(x * x + y * y + z * z);
                    if (r > 1.0e-6) {
                        d.chi_rhs[pp] = d.phi[pp];
                    } else {
                        d.chi_rhs[pp] = 0.0;
                    }
#else
                    d.chi_rhs[pp] = d.phi[pp];
#endif
                }
            }
        }
    });
    // set up the dependencies
    chi_xx_deriv.precede(phi_rhs_task);
    chi_yy_deriv.precede(phi_rhs_task);
    chi_zz_deriv.precede(phi_rhs_task);

    // derivs
    auto chi_x_deriv = tfdendro::nlsmTasks.placeholder().name("chi x deriv");
    auto chi_y_deriv = tfdendro::nlsmTasks.placeholder().name("chi y deriv");
    auto chi_z_deriv = tfdendro::nlsmTasks.placeholder().name("chi z deriv");

    chi_x_deriv.data(&tfdendro::data).work([chi_x_deriv]() {
        auto d =
            *static_cast<tfdendro::dendrotf_rhs_data *>(chi_x_deriv.data());
        if (!d.bflag) deriv_x(d.grad_0_chi, d.chi, d.hx, d.sz, d.bflag);
    });
    chi_y_deriv.data(&tfdendro::data).work([chi_y_deriv]() {
        auto d =
            *static_cast<tfdendro::dendrotf_rhs_data *>(chi_y_deriv.data());
        if (!d.bflag) deriv_y(d.grad_1_chi, d.chi, d.hy, d.sz, d.bflag);
    });
    chi_y_deriv.data(&tfdendro::data).work([chi_y_deriv]() {
        auto d =
            *static_cast<tfdendro::dendrotf_rhs_data *>(chi_y_deriv.data());
        if (!d.bflag) deriv_z(d.grad_2_chi, d.chi, d.hz, d.sz, d.bflag);
    });

    auto phi_x_deriv = tfdendro::nlsmTasks.placeholder().name("phi x deriv");
    auto phi_y_deriv = tfdendro::nlsmTasks.placeholder().name("phi y deriv");
    auto phi_z_deriv = tfdendro::nlsmTasks.placeholder().name("phi z deriv");

    phi_x_deriv.data(&tfdendro::data).work([phi_x_deriv]() {
        auto d =
            *static_cast<tfdendro::dendrotf_rhs_data *>(phi_x_deriv.data());
        if (!d.bflag) deriv_x(d.grad_0_phi, d.phi, d.hx, d.sz, d.bflag);
    });
    phi_y_deriv.data(&tfdendro::data).work([phi_y_deriv]() {
        auto d =
            *static_cast<tfdendro::dendrotf_rhs_data *>(phi_y_deriv.data());
        if (!d.bflag) deriv_y(d.grad_1_phi, d.phi, d.hy, d.sz, d.bflag);
    });
    phi_y_deriv.data(&tfdendro::data).work([phi_y_deriv]() {
        auto d =
            *static_cast<tfdendro::dendrotf_rhs_data *>(phi_y_deriv.data());
        if (!d.bflag) deriv_z(d.grad_2_phi, d.phi, d.hz, d.sz, d.bflag);
    });

    auto phi_rhs_boundary =
        tfdendro::nlsmTasks.placeholder().name("phi_rhs_boundary");
    auto chi_rhs_boundary =
        tfdendro::nlsmTasks.placeholder().name("chi_rhs_boundary");

    phi_rhs_boundary.data(&tfdendro::data).work([phi_rhs_boundary]() {
        auto d = *static_cast<tfdendro::dendrotf_rhs_data *>(
            phi_rhs_boundary.data());
        if (!d.bflag)
            nlsm_bcs(d.phi_rhs, d.phi, d.grad_0_phi, d.grad_1_phi, d.grad_2_phi,
                     d.pmin, d.pmax, 1.0, 0.0, d.sz, d.bflag);
    });
    chi_rhs_boundary.data(&tfdendro::data).work([chi_rhs_boundary]() {
        auto d = *static_cast<tfdendro::dendrotf_rhs_data *>(
            chi_rhs_boundary.data());
        if (!d.bflag)
            nlsm_bcs(d.chi_rhs, d.chi, d.grad_0_chi, d.grad_1_chi, d.grad_2_chi,
                     d.pmin, d.pmax, 1.0, 0.0, d.sz, d.bflag);
    });

    // gradient to boundary condition dependencies
    chi_x_deriv.precede(chi_rhs_boundary);
    chi_y_deriv.precede(chi_rhs_boundary);
    chi_z_deriv.precede(chi_rhs_boundary);
    phi_x_deriv.precede(phi_rhs_boundary);
    phi_y_deriv.precede(phi_rhs_boundary);
    phi_z_deriv.precede(phi_rhs_boundary);

    // now the rhs_boundary conditions must come after what they
    // computed
    phi_rhs_task.precede(phi_rhs_boundary);
    chi_rhs_task.precede(chi_rhs_boundary);

    auto chi_ko_x_deriv =
        tfdendro::nlsmTasks.placeholder().name("chi_ko_x_deriv");
    auto chi_ko_y_deriv =
        tfdendro::nlsmTasks.placeholder().name("chi_ko_y_deriv");
    auto chi_ko_z_deriv =
        tfdendro::nlsmTasks.placeholder().name("chi_ko_z_deriv");
    chi_ko_x_deriv.data(&tfdendro::data).work([chi_ko_x_deriv]() {
        std::cout << "in chi_ko_x_deriv" << std::endl;
        auto d =
            *static_cast<tfdendro::dendrotf_rhs_data *>(chi_ko_x_deriv.data());
        ko_deriv_x(d.grad_0_chi, d.chi, d.hx, d.sz, d.bflag);
        std::cout << "finished chi_ko_x_deriv" << std::endl;
    });
    chi_ko_y_deriv.data(&tfdendro::data).work([chi_ko_y_deriv]() {
        std::cout << "in chi_ko_y_deriv" << std::endl;
        auto d =
            *static_cast<tfdendro::dendrotf_rhs_data *>(chi_ko_y_deriv.data());
        ko_deriv_y(d.grad_1_chi, d.chi, d.hy, d.sz, d.bflag);
        std::cout << "finished chi_ko_y_deriv" << std::endl;
    });
    chi_ko_z_deriv.data(&tfdendro::data).work([chi_ko_z_deriv]() {
        std::cout << "in chi_ko_z_deriv" << std::endl;
        auto d =
            *static_cast<tfdendro::dendrotf_rhs_data *>(chi_ko_z_deriv.data());
        ko_deriv_z(d.grad_2_chi, d.chi, d.hz, d.sz, d.bflag);
        std::cout << "finished chi_ko_z_deriv" << std::endl;
    });
    chi_rhs_boundary.precede(chi_ko_x_deriv);
    chi_rhs_boundary.precede(chi_ko_y_deriv);
    chi_rhs_boundary.precede(chi_ko_z_deriv);

    auto phi_ko_x_deriv =
        tfdendro::nlsmTasks.placeholder().name("phi_ko_x_deriv");
    auto phi_ko_y_deriv =
        tfdendro::nlsmTasks.placeholder().name("phi_ko_y_deriv");
    auto phi_ko_z_deriv =
        tfdendro::nlsmTasks.placeholder().name("phi_ko_z_deriv");
    phi_ko_x_deriv.data(&tfdendro::data).work([phi_ko_x_deriv]() {
        auto d =
            *static_cast<tfdendro::dendrotf_rhs_data *>(phi_ko_x_deriv.data());
        ko_deriv_x(d.grad_0_phi, d.phi, d.hx, d.sz, d.bflag);
    });
    phi_ko_y_deriv.data(&tfdendro::data).work([phi_ko_y_deriv]() {
        auto d =
            *static_cast<tfdendro::dendrotf_rhs_data *>(phi_ko_y_deriv.data());
        ko_deriv_y(d.grad_1_phi, d.phi, d.hy, d.sz, d.bflag);
    });
    phi_ko_z_deriv.data(&tfdendro::data).work([phi_ko_z_deriv]() {
        auto d =
            *static_cast<tfdendro::dendrotf_rhs_data *>(phi_ko_z_deriv.data());
        ko_deriv_z(d.grad_2_phi, d.phi, d.hz, d.sz, d.bflag);
    });
    phi_rhs_boundary.precede(phi_ko_x_deriv);
    phi_rhs_boundary.precede(phi_ko_y_deriv);
    phi_rhs_boundary.precede(phi_ko_z_deriv);

    // then the KO calculation itself
    auto chi_ko_diss = tfdendro::nlsmTasks.placeholder().name("chi_ko_diss");
    chi_ko_diss.data(&tfdendro::data).work([chi_ko_diss]() {
        auto d =
            *static_cast<tfdendro::dendrotf_rhs_data *>(chi_ko_diss.data());
        DendroRegister unsigned int pp;
        for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
            for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
                for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                    pp = i + d.nx * (j + d.ny * k);

                    d.chi_rhs[pp] +=
                        d.sigma * (d.grad_0_chi[pp] + d.grad_1_chi[pp] +
                                   d.grad_2_chi[pp]);
                }
            }
        }
    });
    auto phi_ko_diss = tfdendro::nlsmTasks.placeholder().name("phi_ko_diss");
    phi_ko_diss.data(&tfdendro::data).work([phi_ko_diss]() {
        auto d =
            *static_cast<tfdendro::dendrotf_rhs_data *>(phi_ko_diss.data());
        DendroRegister unsigned int pp;
        for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
            for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
                for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                    pp = i + d.nx * (j + d.ny * k);

                    d.phi_rhs[pp] +=
                        d.sigma * (d.grad_0_phi[pp] + d.grad_1_phi[pp] +
                                   d.grad_2_phi[pp]);
                }
            }
        }
    });

    chi_ko_x_deriv.precede(chi_ko_diss);
    chi_ko_y_deriv.precede(chi_ko_diss);
    chi_ko_z_deriv.precede(chi_ko_diss);
    phi_ko_x_deriv.precede(phi_ko_diss);
    phi_ko_y_deriv.precede(phi_ko_diss);
    phi_ko_z_deriv.precede(phi_ko_diss);

    // and that should be good?

    if (!mpiRank) tfdendro::nlsmTasks.dump(std::cout);
}

void nlsmRHSTaskflowFromPrebuilt(double **unzipVarsRHS, const double **uZipVars,
                                 const unsigned int &offset, const double *pmin,
                                 const double *pmax, const unsigned int *sz,
                                 const unsigned int &bflag) {
    tfdendro::data.chi = (double *)&uZipVars[VAR::U_CHI][offset];
    tfdendro::data.phi = (double *)&uZipVars[VAR::U_PHI][offset];
    tfdendro::data.chi_rhs = (double *)&unzipVarsRHS[VAR::U_CHI][offset];
    tfdendro::data.phi_rhs = (double *)&unzipVarsRHS[VAR::U_PHI][offset];

    tfdendro::data.nx = sz[0];
    tfdendro::data.ny = sz[1];
    tfdendro::data.nz = sz[2];
    tfdendro::data.pmin = (double *)pmin;
    tfdendro::data.pmax = (double *)pmax;

    tfdendro::data.hx = (pmax[0] - pmin[0]) / (tfdendro::data.nx - 1);
    tfdendro::data.hy = (pmax[1] - pmin[1]) / (tfdendro::data.ny - 1);
    tfdendro::data.hz = (pmax[2] - pmin[2]) / (tfdendro::data.nz - 1);

    tfdendro::data.n = sz[0] * sz[1] * sz[2];
    tfdendro::data.PW = nlsm::NLSM_PADDING_WIDTH;

    tfdendro::data.BLK_SZ = tfdendro::data.n;

    tfdendro::data.sigma = KO_DISS_SIGMA;
    tfdendro::data.sz = (unsigned int *)sz;
    tfdendro::data.bflag = bflag;

    std::cout << tfdendro::data.n << std::endl;
    std::cout << tfdendro::data.chi << " " << &uZipVars[VAR::U_CHI][offset]
              << std::endl;

    tfdendro::data.grad_0_chi = new double[tfdendro::data.n];
    tfdendro::data.grad_1_chi = new double[tfdendro::data.n];
    tfdendro::data.grad_2_chi = new double[tfdendro::data.n];
    tfdendro::data.grad_0_phi = new double[tfdendro::data.n];
    tfdendro::data.grad_1_phi = new double[tfdendro::data.n];
    tfdendro::data.grad_2_phi = new double[tfdendro::data.n];
    tfdendro::data.grad2_0_0_chi = new double[tfdendro::data.n];
    tfdendro::data.grad2_1_1_chi = new double[tfdendro::data.n];
    tfdendro::data.grad2_2_2_chi = new double[tfdendro::data.n];

    std::cout << tfdendro::data.grad2_0_0_chi << std::endl;

    // run execution
    nlsm::timer::t_rhs.start();
    tfdendro::executor.run(tfdendro::nlsmTasks).wait();
    nlsm::timer::t_rhs.stop();

    delete[] tfdendro::data.grad_0_chi;
    delete[] tfdendro::data.grad_1_chi;
    delete[] tfdendro::data.grad_2_chi;
    delete[] tfdendro::data.grad_0_phi;
    delete[] tfdendro::data.grad_1_phi;
    delete[] tfdendro::data.grad_2_phi;
    delete[] tfdendro::data.grad2_0_0_chi;
    delete[] tfdendro::data.grad2_1_1_chi;
    delete[] tfdendro::data.grad2_2_2_chi;
}

void nlsmRHSTaskflowEntryBUILDONENTRY(double **unzipVarsRHS,
                                      const double **uZipVars,
                                      const unsigned int &offset,
                                      const double *pmin, const double *pmax,
                                      const unsigned int *sz,
                                      const unsigned int &bflag) {
    // extract the vars
    const double *chi = &uZipVars[VAR::U_CHI][offset];
    const double *phi = &uZipVars[VAR::U_PHI][offset];

    double *chi_rhs = &unzipVarsRHS[VAR::U_CHI][offset];
    double *phi_rhs = &unzipVarsRHS[VAR::U_PHI][offset];

    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];

    double hx = (pmax[0] - pmin[0]) / (nx - 1);
    double hy = (pmax[1] - pmin[1]) / (ny - 1);
    double hz = (pmax[2] - pmin[2]) / (nz - 1);

    int idx[3];

    unsigned int n = sz[0] * sz[1] * sz[2];
    const unsigned int PW = nlsm::NLSM_PADDING_WIDTH;

    nlsm::timer::t_deriv.start();

    double *grad_0_chi = new double[n];
    double *grad_1_chi = new double[n];
    double *grad_2_chi = new double[n];

    double *grad_0_phi = new double[n];
    double *grad_1_phi = new double[n];
    double *grad_2_phi = new double[n];

    double *grad2_0_0_chi = new double[n];
    double *grad2_1_1_chi = new double[n];
    double *grad2_2_2_chi = new double[n];

    const double sigma = KO_DISS_SIGMA;

    // the a timer is timing how long it takes to initialize the
    // graph when we enter the function
    nlsm::timer::t_rhs_a.stop();
    tf::Taskflow taskflow;

    auto chi_xx_deriv =
        taskflow.emplace([&]() { deriv_xx(grad2_0_0_chi, chi, hx, sz, bflag); })
            .name("chi_xx_deriv");
    auto chi_yy_deriv =
        taskflow.emplace([&]() { deriv_yy(grad2_1_1_chi, chi, hy, sz, bflag); })
            .name("chi_yy_deriv");
    auto chi_zz_deriv =
        taskflow.emplace([&]() { deriv_zz(grad2_2_2_chi, chi, hz, sz, bflag); })
            .name("chi_zz_deriv");

    // then set up the different variable tasks
    auto phi_rhs_task =
        taskflow
            .emplace([&]() {
                DendroRegister double x;
                DendroRegister double y;
                DendroRegister double z;
                DendroRegister unsigned int pp;
                double r;

                for (unsigned int k = PW; k < nz - PW; k++) {
                    z = pmin[2] + k * hz;
                    for (unsigned int j = PW; j < ny - PW; j++) {
                        y = pmin[1] + j * hy;
                        for (unsigned int i = PW; i < nx - PW; i++) {
                            x = pmin[0] + i * hx;
                            pp = i + nx * (j + ny * k);

#ifdef NLSM_NONLINEAR
                            r = sqrt(x * x + y * y + z * z);
                            if (r > 1.0e-6) {
                                phi_rhs[pp] =
                                    NLSM_WAVE_SPEED_X * grad2_0_0_chi[pp] +
                                    NLSM_WAVE_SPEED_Y * grad2_1_1_chi[pp] +
                                    NLSM_WAVE_SPEED_Z * grad2_2_2_chi[pp] -
                                    sin(2 * chi[pp]) / pow(r, 2);
                            } else {
                                chi_rhs[pp] = 0.0;
                            }
#else
                            phi_rhs[pp] =
                                NLSM_WAVE_SPEED_X * grad2_0_0_chi[pp] +
                                NLSM_WAVE_SPEED_Y * grad2_1_1_chi[pp] +
                                NLSM_WAVE_SPEED_Z * grad2_2_2_chi[pp];
#endif
                        }
                    }
                }
            })
            .name("phi_rhs_task");

    auto chi_rhs_task =
        taskflow
            .emplace([&]() {
                DendroRegister double x;
                DendroRegister double y;
                DendroRegister double z;
                DendroRegister unsigned int pp;
                double r;

                for (unsigned int k = PW; k < nz - PW; k++) {
                    z = pmin[2] + k * hz;
                    for (unsigned int j = PW; j < ny - PW; j++) {
                        y = pmin[1] + j * hy;
                        for (unsigned int i = PW; i < nx - PW; i++) {
                            x = pmin[0] + i * hx;
                            pp = i + nx * (j + ny * k);

#ifdef NLSM_NONLINEAR
                            r = sqrt(x * x + y * y + z * z);
                            if (r > 1.0e-6) {
                                chi_rhs[pp] = phi[pp];
                            } else {
                                chi_rhs[pp] = 0.0;
                            }
#else
                            chi_rhs[pp] = phi[pp];
#endif
                        }
                    }
                }
            })
            .name("chi_rhs_task");

    // set up the dependencies
    chi_xx_deriv.precede(phi_rhs_task);
    chi_yy_deriv.precede(phi_rhs_task);
    chi_zz_deriv.precede(phi_rhs_task);

    // NOTE: this could technically be optimized to avoid other
    // computations... NOTE: (2) additionally, these can be computed
    // at any time
    auto chi_x_deriv = taskflow
                           .emplace([&]() {
                               if (!bflag)
                                   deriv_x(grad_0_chi, chi, hx, sz, bflag);
                           })
                           .name("chi_x_deriv");
    auto chi_y_deriv = taskflow
                           .emplace([&]() {
                               if (!bflag)
                                   deriv_y(grad_1_chi, chi, hy, sz, bflag);
                           })
                           .name("chi_y_deriv");
    auto chi_z_deriv = taskflow
                           .emplace([&]() {
                               if (!bflag)
                                   deriv_z(grad_2_chi, chi, hz, sz, bflag);
                           })
                           .name("chi_z_deriv");

    auto phi_x_deriv = taskflow
                           .emplace([&]() {
                               if (!bflag)
                                   deriv_x(grad_0_phi, phi, hx, sz, bflag);
                           })
                           .name("phi_x_deriv");
    auto phi_y_deriv = taskflow
                           .emplace([&]() {
                               if (!bflag)
                                   deriv_y(grad_1_phi, phi, hy, sz, bflag);
                           })
                           .name("phi_y_deriv");
    auto phi_z_deriv = taskflow
                           .emplace([&]() {
                               if (!bflag)
                                   deriv_z(grad_2_phi, phi, hz, sz, bflag);
                           })
                           .name("phi_z_deriv");

    auto phi_rhs_boundary =
        taskflow
            .emplace([&]() {
                if (!bflag)
                    nlsm_bcs(phi_rhs, phi, grad_0_phi, grad_1_phi, grad_2_phi,
                             pmin, pmax, 1.0, 0.0, sz, bflag);
            })
            .name("phi_rhs_boundary");
    auto chi_rhs_boundary =
        taskflow
            .emplace([&]() {
                if (!bflag)
                    nlsm_bcs(chi_rhs, chi, grad_0_chi, grad_1_chi, grad_2_chi,
                             pmin, pmax, 1.0, 0.0, sz, bflag);
            })
            .name("chi_rhs_boundary");

    // gradient to boundary condition dependencies
    chi_x_deriv.precede(chi_rhs_boundary);
    chi_y_deriv.precede(chi_rhs_boundary);
    chi_z_deriv.precede(chi_rhs_boundary);
    phi_x_deriv.precede(phi_rhs_boundary);
    phi_y_deriv.precede(phi_rhs_boundary);
    phi_z_deriv.precede(phi_rhs_boundary);

    // now the rhs_boundary conditions must come after what they
    // computed
    phi_rhs_task.precede(phi_rhs_boundary);
    chi_rhs_task.precede(chi_rhs_boundary);

    // now we can do the KO_DISS_SIGMA
    // this must come after the boundary task finishes
    auto chi_ko_x_deriv =
        taskflow.emplace([&]() { ko_deriv_x(grad_0_chi, chi, hx, sz, bflag); })
            .name("chi_ko_x_deriv");
    auto chi_ko_y_deriv =
        taskflow.emplace([&]() { ko_deriv_y(grad_1_chi, chi, hy, sz, bflag); })
            .name("chi_ko_y_deriv");
    auto chi_ko_z_deriv =
        taskflow.emplace([&]() { ko_deriv_z(grad_2_chi, chi, hz, sz, bflag); })
            .name("chi_ko_z_deriv");
    chi_rhs_boundary.precede(chi_ko_x_deriv);
    chi_rhs_boundary.precede(chi_ko_y_deriv);
    chi_rhs_boundary.precede(chi_ko_z_deriv);

    auto phi_ko_x_deriv =
        taskflow.emplace([&]() { ko_deriv_x(grad_0_phi, phi, hx, sz, bflag); })
            .name("phi_ko_x_deriv");
    auto phi_ko_y_deriv =
        taskflow.emplace([&]() { ko_deriv_y(grad_1_phi, phi, hy, sz, bflag); })
            .name("phi_ko_y_deriv");
    auto phi_ko_z_deriv =
        taskflow.emplace([&]() { ko_deriv_z(grad_2_phi, phi, hz, sz, bflag); })
            .name("phi_ko_z_deriv");
    phi_rhs_boundary.precede(phi_ko_x_deriv);
    phi_rhs_boundary.precede(phi_ko_y_deriv);
    phi_rhs_boundary.precede(phi_ko_z_deriv);

    // then the KO calculation itself
    auto chi_ko_diss =
        taskflow
            .emplace([&]() {
                DendroRegister unsigned int pp;
                for (unsigned int k = PW; k < nz - PW; k++) {
                    for (unsigned int j = PW; j < ny - PW; j++) {
                        for (unsigned int i = PW; i < nx - PW; i++) {
                            pp = i + nx * (j + ny * k);

                            chi_rhs[pp] +=
                                sigma * (grad_0_chi[pp] + grad_1_chi[pp] +
                                         grad_2_chi[pp]);
                        }
                    }
                }
            })
            .name("chi_ko_diss");
    auto phi_ko_diss =
        taskflow
            .emplace([&]() {
                DendroRegister unsigned int pp;
                for (unsigned int k = PW; k < nz - PW; k++) {
                    for (unsigned int j = PW; j < ny - PW; j++) {
                        for (unsigned int i = PW; i < nx - PW; i++) {
                            pp = i + nx * (j + ny * k);

                            phi_rhs[pp] +=
                                sigma * (grad_0_phi[pp] + grad_1_phi[pp] +
                                         grad_2_phi[pp]);
                        }
                    }
                }
            })
            .name("phi_ko_diss");

    chi_ko_x_deriv.precede(chi_ko_diss);
    chi_ko_y_deriv.precede(chi_ko_diss);
    chi_ko_z_deriv.precede(chi_ko_diss);
    phi_ko_x_deriv.precede(phi_ko_diss);
    phi_ko_y_deriv.precede(phi_ko_diss);
    phi_ko_z_deriv.precede(phi_ko_diss);

    // make the executor based on the threads to use
    tf::Executor executor(tfdendro::threadsUse);
    nlsm::timer::t_rhs_a.stop();

    // taskflow.dump(std::cout);

    nlsm::timer::t_rhs.start();
    executor.run(taskflow).wait();
    nlsm::timer::t_rhs.stop();

    delete[] grad2_0_0_chi;
    delete[] grad2_1_1_chi;
    delete[] grad2_2_2_chi;

    delete[] grad_0_chi;
    delete[] grad_1_chi;
    delete[] grad_2_chi;

    delete[] grad_0_phi;
    delete[] grad_1_phi;
    delete[] grad_2_phi;
}
