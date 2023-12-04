#include "../../lib/taskflow/taskflow.hpp"
// #include "dendroTF_helpers.h"
#include <cstdlib>
#include <thread>

#include "hadrhs.h"
#include "rhs.h"

using namespace std;
using namespace bssn;

namespace tfdendro {
struct dendrotf_rhs_data {
    double *alpha;
    double *chi;
    double *K;
    double *gt0;
    double *gt1;
    double *gt2;
    double *gt3;
    double *gt4;
    double *gt5;
    double *beta0;
    double *beta1;
    double *beta2;
    double *At0;
    double *At1;
    double *At2;
    double *At3;
    double *At4;
    double *At5;
    double *Gt0;
    double *Gt1;
    double *Gt2;
    double *B0;
    double *B1;
    double *B2;

    double *a_rhs;
    double *chi_rhs;
    double *K_rhs;
    double *gt_rhs00;
    double *gt_rhs01;
    double *gt_rhs02;
    double *gt_rhs11;
    double *gt_rhs12;
    double *gt_rhs22;
    double *b_rhs0;
    double *b_rhs1;
    double *b_rhs2;
    double *At_rhs00;
    double *At_rhs01;
    double *At_rhs02;
    double *At_rhs11;
    double *At_rhs12;
    double *At_rhs22;
    double *Gt_rhs0;
    double *Gt_rhs1;
    double *Gt_rhs2;
    double *B_rhs0;
    double *B_rhs1;
    double *B_rhs2;

    unsigned int nx;
    unsigned int ny;
    unsigned int nz;

    double hx;
    double hy;
    double hz;

    unsigned int lambda[4];
    double lambda_f[2];
    double eta_power[2];

    int idx[3];

    unsigned int PW;

    unsigned int n;

    unsigned int BLK_SZ;

    double *deriv_base;

    double sigma;

    double **unzipVarsRHS;
    double **uZipVars;
    unsigned int offset;
    double *pmin;
    double *pmax;
    unsigned int *sz;
    unsigned int bflag;
};

char *executorThreads = std::getenv("DENDRO_TF_EXECUTOR_THREADS");
int threadsUse = executorThreads == NULL ? std::thread::hardware_concurrency()
                                         : atoi(executorThreads);

tf::Executor executor(threadsUse);
tf::Taskflow taskflow;

tf::Taskflow flow_deriv;

tf::Taskflow flow_kodiss;

dendrotf_rhs_data data;
}  // namespace tfdendro

void buildDerivativeGraph() {
#include "RHSINCLUDE_deriv_data_dependency.cpp"
}

void buildGraphVarsOnly() {
    // all of the tasks that require our actual data need to placeholders
    auto a_comp = tfdendro::taskflow.placeholder().name("a RHS");
    auto b_comp = tfdendro::taskflow.placeholder().name("b RHS");
    auto gt_comp = tfdendro::taskflow.placeholder().name("gt RHS");
    auto chi_comp = tfdendro::taskflow.placeholder().name("chi RHS");
    auto At_comp = tfdendro::taskflow.placeholder().name("At RHS");
    auto K_comp = tfdendro::taskflow.placeholder().name("K RHS");
    auto Gt_comp = tfdendro::taskflow.placeholder().name("Gt RHS");
    auto B_comp = tfdendro::taskflow.placeholder().name("B RHS");

    a_comp.data(&tfdendro::data).work([a_comp]() {
        auto d = *static_cast<tfdendro::dendrotf_rhs_data *>(a_comp.data());
        rhs_a_only(d.unzipVarsRHS, (const double **)d.uZipVars, d.offset,
                   d.pmin, d.pmax, d.sz, d.bflag);
    });
    b_comp.data(&tfdendro::data).work([b_comp]() {
        auto d = *static_cast<tfdendro::dendrotf_rhs_data *>(b_comp.data());
        rhs_b_only(d.unzipVarsRHS, (const double **)d.uZipVars, d.offset,
                   d.pmin, d.pmax, d.sz, d.bflag);
    });
    gt_comp.data(&tfdendro::data).work([gt_comp]() {
        auto d = *static_cast<tfdendro::dendrotf_rhs_data *>(gt_comp.data());
        rhs_gt_only(d.unzipVarsRHS, (const double **)d.uZipVars, d.offset,
                    d.pmin, d.pmax, d.sz, d.bflag);
    });
    chi_comp.data(&tfdendro::data).work([chi_comp]() {
        auto d = *static_cast<tfdendro::dendrotf_rhs_data *>(chi_comp.data());
        rhs_chi_only(d.unzipVarsRHS, (const double **)d.uZipVars, d.offset,
                     d.pmin, d.pmax, d.sz, d.bflag);
    });
    At_comp.data(&tfdendro::data).work([At_comp]() {
        auto d = *static_cast<tfdendro::dendrotf_rhs_data *>(At_comp.data());
        rhs_At_only(d.unzipVarsRHS, (const double **)d.uZipVars, d.offset,
                    d.pmin, d.pmax, d.sz, d.bflag);
    });
    K_comp.data(&tfdendro::data).work([K_comp]() {
        auto d = *static_cast<tfdendro::dendrotf_rhs_data *>(K_comp.data());
        rhs_K_only(d.unzipVarsRHS, (const double **)d.uZipVars, d.offset,
                   d.pmin, d.pmax, d.sz, d.bflag);
    });
    Gt_comp.data(&tfdendro::data).work([Gt_comp]() {
        auto d = *static_cast<tfdendro::dendrotf_rhs_data *>(Gt_comp.data());
        rhs_Gt_only(d.unzipVarsRHS, (const double **)d.uZipVars, d.offset,
                    d.pmin, d.pmax, d.sz, d.bflag);
    });
    B_comp.data(&tfdendro::data).work([B_comp]() {
        auto d = *static_cast<tfdendro::dendrotf_rhs_data *>(B_comp.data());
        rhs_B_only(d.unzipVarsRHS, (const double **)d.uZipVars, d.offset,
                   d.pmin, d.pmax, d.sz, d.bflag);
    });

    // tfdendro::taskflow.dump(std::cout);
}

void buildKOGraph() {
#include "RHSINCLUDE_kodiss_data.cpp"
}

void generateFullRHSGraph(int mpiRank) {
    if (!mpiRank)
        std::cout << "num workers: " << tfdendro::executor.num_workers()
                  << std::endl;

    // start by building derivative graph
    buildDerivativeGraph();

    // add this as a task to the main taskflow graph
    tf::Task derivTask = tfdendro::taskflow.composed_of(tfdendro::flow_deriv)
                             .name("Derivative Computation");

    // then establish the variables and data
    // all of the tasks that require our actual data need to placeholders
    auto a_comp = tfdendro::taskflow.placeholder().name("a RHS");
    auto b_comp = tfdendro::taskflow.placeholder().name("b RHS");
    auto gt_comp = tfdendro::taskflow.placeholder().name("gt RHS");
    auto chi_comp = tfdendro::taskflow.placeholder().name("chi RHS");
    auto At_comp = tfdendro::taskflow.placeholder().name("At RHS");
    auto K_comp = tfdendro::taskflow.placeholder().name("K RHS");
    auto Gt_comp = tfdendro::taskflow.placeholder().name("Gt RHS");
    auto B_comp = tfdendro::taskflow.placeholder().name("B RHS");

    a_comp.data(&tfdendro::data).work([a_comp]() {
        auto d = *static_cast<tfdendro::dendrotf_rhs_data *>(a_comp.data());
        rhs_a_only(d.unzipVarsRHS, (const double **)d.uZipVars, d.offset,
                   d.pmin, d.pmax, d.sz, d.bflag);
    });
    b_comp.data(&tfdendro::data).work([b_comp]() {
        auto d = *static_cast<tfdendro::dendrotf_rhs_data *>(b_comp.data());
        rhs_b_only(d.unzipVarsRHS, (const double **)d.uZipVars, d.offset,
                   d.pmin, d.pmax, d.sz, d.bflag);
    });
    gt_comp.data(&tfdendro::data).work([gt_comp]() {
        auto d = *static_cast<tfdendro::dendrotf_rhs_data *>(gt_comp.data());
        rhs_gt_only(d.unzipVarsRHS, (const double **)d.uZipVars, d.offset,
                    d.pmin, d.pmax, d.sz, d.bflag);
    });
    chi_comp.data(&tfdendro::data).work([chi_comp]() {
        auto d = *static_cast<tfdendro::dendrotf_rhs_data *>(chi_comp.data());
        rhs_chi_only(d.unzipVarsRHS, (const double **)d.uZipVars, d.offset,
                     d.pmin, d.pmax, d.sz, d.bflag);
    });
    At_comp.data(&tfdendro::data).work([At_comp]() {
        auto d = *static_cast<tfdendro::dendrotf_rhs_data *>(At_comp.data());
        rhs_At_only(d.unzipVarsRHS, (const double **)d.uZipVars, d.offset,
                    d.pmin, d.pmax, d.sz, d.bflag);
    });
    K_comp.data(&tfdendro::data).work([K_comp]() {
        auto d = *static_cast<tfdendro::dendrotf_rhs_data *>(K_comp.data());
        rhs_K_only(d.unzipVarsRHS, (const double **)d.uZipVars, d.offset,
                   d.pmin, d.pmax, d.sz, d.bflag);
    });
    Gt_comp.data(&tfdendro::data).work([Gt_comp]() {
        auto d = *static_cast<tfdendro::dendrotf_rhs_data *>(Gt_comp.data());
        rhs_Gt_only(d.unzipVarsRHS, (const double **)d.uZipVars, d.offset,
                    d.pmin, d.pmax, d.sz, d.bflag);
    });
    B_comp.data(&tfdendro::data).work([B_comp]() {
        auto d = *static_cast<tfdendro::dendrotf_rhs_data *>(B_comp.data());
        rhs_B_only(d.unzipVarsRHS, (const double **)d.uZipVars, d.offset,
                   d.pmin, d.pmax, d.sz, d.bflag);
    });

    // then make sure all derivatives are computed first
    derivTask.precede(a_comp);
    derivTask.precede(b_comp);
    derivTask.precede(gt_comp);
    derivTask.precede(chi_comp);
    derivTask.precede(At_comp);
    derivTask.precede(K_comp);
    derivTask.precede(Gt_comp);
    derivTask.precede(B_comp);

    // then build the ko diss graph
    buildKOGraph();

    // add this as a task to the main taskflow graph
    tf::Task koDissTask = tfdendro::taskflow.composed_of(tfdendro::flow_kodiss)
                              .name("KO Diss. Computation");

    // dependency so KO portion is done after
    a_comp.precede(koDissTask);
    b_comp.precede(koDissTask);
    gt_comp.precede(koDissTask);
    chi_comp.precede(koDissTask);
    At_comp.precede(koDissTask);
    K_comp.precede(koDissTask);
    Gt_comp.precede(koDissTask);
    B_comp.precede(koDissTask);

    // then we're good to go!

    // if (!mpiRank) {
    //     tfdendro::taskflow.dump(std::cout);
    // }
}

void bssnRHS_TF(double **uzipVarsRHS, const double **uZipVars,
                const ot::Block *blkList, unsigned int numBlocks) {
    // TODO: try task flowing out on this level
    tf::Taskflow taskflow;

    unsigned int offset;
    double ptmin[3], ptmax[3];
    unsigned int sz[3];
    unsigned int bflag;
    double dx, dy, dz;
    const Point pt_min(bssn::BSSN_COMPD_MIN[0], bssn::BSSN_COMPD_MIN[1],
                       bssn::BSSN_COMPD_MIN[2]);
    const Point pt_max(bssn::BSSN_COMPD_MAX[0], bssn::BSSN_COMPD_MAX[1],
                       bssn::BSSN_COMPD_MAX[2]);
    const unsigned int PW = bssn::BSSN_PADDING_WIDTH;

    for (unsigned int blk = 0; blk < numBlocks; blk++) {
        offset = blkList[blk].getOffset();
        sz[0] = blkList[blk].getAllocationSzX();
        sz[1] = blkList[blk].getAllocationSzY();
        sz[2] = blkList[blk].getAllocationSzZ();

        bflag = blkList[blk].getBlkNodeFlag();

        dx = blkList[blk].computeDx(pt_min, pt_max);
        dy = blkList[blk].computeDy(pt_min, pt_max);
        dz = blkList[blk].computeDz(pt_min, pt_max);

        ptmin[0] = GRIDX_TO_X(blkList[blk].getBlockNode().minX()) - PW * dx;
        ptmin[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().minY()) - PW * dy;
        ptmin[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ()) - PW * dz;

        ptmax[0] = GRIDX_TO_X(blkList[blk].getBlockNode().maxX()) + PW * dx;
        ptmax[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().maxY()) + PW * dy;
        ptmax[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ()) + PW * dz;

        bssnrhs_sep_taskflow_varsOnly(uzipVarsRHS, (const double **)uZipVars,
                                      offset, ptmin, ptmax, sz, bflag);
    }
}

void bssnrhs_sep_taskflow_FULL(double **unzipVarsRHS, const double **uZipVars,
                               const unsigned int &offset, const double *pmin,
                               const double *pmax, const unsigned int *sz,
                               const unsigned int &bflag,
                               const unsigned int threadID,
                               const unsigned int maxBlockSize) {
    tf::Executor executor;
    tf::Taskflow taskflow;

    // grab all of the pointers...

    const double *alpha = &uZipVars[VAR::U_ALPHA][offset];
    const double *chi = &uZipVars[VAR::U_CHI][offset];
    const double *K = &uZipVars[VAR::U_K][offset];
    const double *gt0 = &uZipVars[VAR::U_SYMGT0][offset];
    const double *gt1 = &uZipVars[VAR::U_SYMGT1][offset];
    const double *gt2 = &uZipVars[VAR::U_SYMGT2][offset];
    const double *gt3 = &uZipVars[VAR::U_SYMGT3][offset];
    const double *gt4 = &uZipVars[VAR::U_SYMGT4][offset];
    const double *gt5 = &uZipVars[VAR::U_SYMGT5][offset];
    const double *beta0 = &uZipVars[VAR::U_BETA0][offset];
    const double *beta1 = &uZipVars[VAR::U_BETA1][offset];
    const double *beta2 = &uZipVars[VAR::U_BETA2][offset];
    const double *At0 = &uZipVars[VAR::U_SYMAT0][offset];
    const double *At1 = &uZipVars[VAR::U_SYMAT1][offset];
    const double *At2 = &uZipVars[VAR::U_SYMAT2][offset];
    const double *At3 = &uZipVars[VAR::U_SYMAT3][offset];
    const double *At4 = &uZipVars[VAR::U_SYMAT4][offset];
    const double *At5 = &uZipVars[VAR::U_SYMAT5][offset];
    const double *Gt0 = &uZipVars[VAR::U_GT0][offset];
    const double *Gt1 = &uZipVars[VAR::U_GT1][offset];
    const double *Gt2 = &uZipVars[VAR::U_GT2][offset];
    const double *B0 = &uZipVars[VAR::U_B0][offset];
    const double *B1 = &uZipVars[VAR::U_B1][offset];
    const double *B2 = &uZipVars[VAR::U_B2][offset];

    double *const a_rhs = &unzipVarsRHS[VAR::U_ALPHA][offset];
    double *const chi_rhs = &unzipVarsRHS[VAR::U_CHI][offset];
    double *const K_rhs = &unzipVarsRHS[VAR::U_K][offset];
    double *const gt_rhs00 = &unzipVarsRHS[VAR::U_SYMGT0][offset];
    double *const gt_rhs01 = &unzipVarsRHS[VAR::U_SYMGT1][offset];
    double *const gt_rhs02 = &unzipVarsRHS[VAR::U_SYMGT2][offset];
    double *const gt_rhs11 = &unzipVarsRHS[VAR::U_SYMGT3][offset];
    double *const gt_rhs12 = &unzipVarsRHS[VAR::U_SYMGT4][offset];
    double *const gt_rhs22 = &unzipVarsRHS[VAR::U_SYMGT5][offset];
    double *const b_rhs0 = &unzipVarsRHS[VAR::U_BETA0][offset];
    double *const b_rhs1 = &unzipVarsRHS[VAR::U_BETA1][offset];
    double *const b_rhs2 = &unzipVarsRHS[VAR::U_BETA2][offset];
    double *const At_rhs00 = &unzipVarsRHS[VAR::U_SYMAT0][offset];
    double *const At_rhs01 = &unzipVarsRHS[VAR::U_SYMAT1][offset];
    double *const At_rhs02 = &unzipVarsRHS[VAR::U_SYMAT2][offset];
    double *const At_rhs11 = &unzipVarsRHS[VAR::U_SYMAT3][offset];
    double *const At_rhs12 = &unzipVarsRHS[VAR::U_SYMAT4][offset];
    double *const At_rhs22 = &unzipVarsRHS[VAR::U_SYMAT5][offset];
    double *const Gt_rhs0 = &unzipVarsRHS[VAR::U_GT0][offset];
    double *const Gt_rhs1 = &unzipVarsRHS[VAR::U_GT1][offset];
    double *const Gt_rhs2 = &unzipVarsRHS[VAR::U_GT2][offset];
    double *const B_rhs0 = &unzipVarsRHS[VAR::U_B0][offset];
    double *const B_rhs1 = &unzipVarsRHS[VAR::U_B1][offset];
    double *const B_rhs2 = &unzipVarsRHS[VAR::U_B2][offset];

    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];

    double hx = (pmax[0] - pmin[0]) / (nx - 1);
    double hy = (pmax[1] - pmin[1]) / (ny - 1);
    double hz = (pmax[2] - pmin[2]) / (nz - 1);

    const unsigned int lambda[4] = {BSSN_LAMBDA[0], BSSN_LAMBDA[1],
                                    BSSN_LAMBDA[2], BSSN_LAMBDA[3]};
    const double lambda_f[2] = {BSSN_LAMBDA_F[0], BSSN_LAMBDA_F[1]};
    const double eta_power[2] = {BSSN_ETA_POWER[0], BSSN_ETA_POWER[1]};

    int idx[3];
    const unsigned int PW = bssn::BSSN_PADDING_WIDTH;
    unsigned int n = sz[0] * sz[1] * sz[2];
    mem::memory_pool<double> *__mem_pool = &BSSN_MEM_POOL;

    bssn::timer::t_deriv.start();

    const unsigned int BLK_SZ = n;
    double *const deriv_base = bssn::BSSN_DERIV_WORKSPACE;

    const double sigma = KO_DISS_SIGMA;

    // clang-format off
// derivative pointers, which we'll keep
#include "bssnrhs_evar_derivs.h"
    // clang-format on

    tf::Taskflow flow_deriv;
    tf::Taskflow flow_kodiss;

    // clang-format off
#include "RHSINCLUDE_taskflow_derivs.cpp"
    // clang-format on

    tf::Task derivTask =
        taskflow.composed_of(flow_deriv).name("derivComputation");

    // then add in all of the RHS computations as tasks
    auto [a_comp, b_comp, gt_comp, chi_comp, At_comp, K_comp, Gt_comp, B_comp] =
        taskflow.emplace(
            [&]() {
                rhs_a_only(unzipVarsRHS, uZipVars, offset, pmin, pmax, sz,
                           bflag);
                // NOTE: temporary location
                if (bflag != 0) {
                    bssn_bcs(a_rhs, alpha, grad_0_alpha, grad_1_alpha,
                             grad_2_alpha, pmin, pmax, 1.0, 1.0, sz, bflag);
                }
            },
            [&]() {
                rhs_b_only(unzipVarsRHS, uZipVars, offset, pmin, pmax, sz,
                           bflag);
                // NOTE: temporary location
                if (bflag != 0) {
                    bssn_bcs(b_rhs0, beta0, grad_0_beta0, grad_1_beta0,
                             grad_2_beta0, pmin, pmax, 1.0, 0.0, sz, bflag);
                    bssn_bcs(b_rhs1, beta1, grad_0_beta1, grad_1_beta1,
                             grad_2_beta1, pmin, pmax, 1.0, 0.0, sz, bflag);
                    bssn_bcs(b_rhs2, beta2, grad_0_beta2, grad_1_beta2,
                             grad_2_beta2, pmin, pmax, 1.0, 0.0, sz, bflag);
                }
            },
            [&]() {
                rhs_gt_only(unzipVarsRHS, uZipVars, offset, pmin, pmax, sz,
                            bflag);
                // NOTE: temporary location
                if (bflag != 0) {
                    bssn_bcs(gt_rhs00, gt0, grad_0_gt0, grad_1_gt0, grad_2_gt0,
                             pmin, pmax, 1.0, 1.0, sz, bflag);
                    bssn_bcs(gt_rhs01, gt1, grad_0_gt1, grad_1_gt1, grad_2_gt1,
                             pmin, pmax, 1.0, 0.0, sz, bflag);
                    bssn_bcs(gt_rhs02, gt2, grad_0_gt2, grad_1_gt2, grad_2_gt2,
                             pmin, pmax, 1.0, 0.0, sz, bflag);
                    bssn_bcs(gt_rhs11, gt3, grad_0_gt3, grad_1_gt3, grad_2_gt3,
                             pmin, pmax, 1.0, 1.0, sz, bflag);
                    bssn_bcs(gt_rhs12, gt4, grad_0_gt4, grad_1_gt4, grad_2_gt4,
                             pmin, pmax, 1.0, 0.0, sz, bflag);
                    bssn_bcs(gt_rhs22, gt5, grad_0_gt5, grad_1_gt5, grad_2_gt5,
                             pmin, pmax, 1.0, 1.0, sz, bflag);
                }
            },
            [&]() {
                rhs_chi_only(unzipVarsRHS, uZipVars, offset, pmin, pmax, sz,
                             bflag);
                // NOTE: temporary location
                if (bflag != 0) {
                    bssn_bcs(chi_rhs, chi, grad_0_chi, grad_1_chi, grad_2_chi,
                             pmin, pmax, 1.0, 1.0, sz, bflag);
                }
            },
            [&]() {
                rhs_At_only(unzipVarsRHS, uZipVars, offset, pmin, pmax, sz,
                            bflag);
                // NOTE: temporary location
                if (bflag != 0) {
                    bssn_bcs(At_rhs00, At0, grad_0_At0, grad_1_At0, grad_2_At0,
                             pmin, pmax, 2.0, 0.0, sz, bflag);
                    bssn_bcs(At_rhs01, At1, grad_0_At1, grad_1_At1, grad_2_At1,
                             pmin, pmax, 2.0, 0.0, sz, bflag);
                    bssn_bcs(At_rhs02, At2, grad_0_At2, grad_1_At2, grad_2_At2,
                             pmin, pmax, 2.0, 0.0, sz, bflag);
                    bssn_bcs(At_rhs11, At3, grad_0_At3, grad_1_At3, grad_2_At3,
                             pmin, pmax, 2.0, 0.0, sz, bflag);
                    bssn_bcs(At_rhs12, At4, grad_0_At4, grad_1_At4, grad_2_At4,
                             pmin, pmax, 2.0, 0.0, sz, bflag);
                    bssn_bcs(At_rhs22, At5, grad_0_At5, grad_1_At5, grad_2_At5,
                             pmin, pmax, 2.0, 0.0, sz, bflag);
                }
            },
            [&]() {
                rhs_K_only(unzipVarsRHS, uZipVars, offset, pmin, pmax, sz,
                           bflag);
                // NOTE: temporary location
                if (bflag != 0) {
                    bssn_bcs(K_rhs, K, grad_0_K, grad_1_K, grad_2_K, pmin, pmax,
                             1.0, 0.0, sz, bflag);
                }
            },
            [&]() {
                rhs_Gt_only(unzipVarsRHS, uZipVars, offset, pmin, pmax, sz,
                            bflag);
                // NOTE: temporary location
                if (bflag != 0) {
                    bssn_bcs(Gt_rhs0, Gt0, grad_0_Gt0, grad_1_Gt0, grad_2_Gt0,
                             pmin, pmax, 2.0, 0.0, sz, bflag);
                    bssn_bcs(Gt_rhs1, Gt1, grad_0_Gt1, grad_1_Gt1, grad_2_Gt1,
                             pmin, pmax, 2.0, 0.0, sz, bflag);
                    bssn_bcs(Gt_rhs2, Gt2, grad_0_Gt2, grad_1_Gt2, grad_2_Gt2,
                             pmin, pmax, 2.0, 0.0, sz, bflag);
                }
            },
            [&]() {
                rhs_B_only(unzipVarsRHS, uZipVars, offset, pmin, pmax, sz,
                           bflag);
                // NOTE: temporary location
                if (bflag != 0) {
                    bssn_bcs(B_rhs0, B0, grad_0_B0, grad_1_B0, grad_2_B0, pmin,
                             pmax, 1.0, 0.0, sz, bflag);
                    bssn_bcs(B_rhs1, B1, grad_0_B1, grad_1_B1, grad_2_B1, pmin,
                             pmax, 1.0, 0.0, sz, bflag);
                    bssn_bcs(B_rhs2, B2, grad_0_B2, grad_1_B2, grad_2_B2, pmin,
                             pmax, 1.0, 0.0, sz, bflag);
                }
            });

    // make sure the derivative task comes before the RHS tasks
    derivTask.precede(a_comp);
    derivTask.precede(b_comp);
    derivTask.precede(gt_comp);
    derivTask.precede(chi_comp);
    derivTask.precede(At_comp);
    derivTask.precede(K_comp);
    derivTask.precede(Gt_comp);
    derivTask.precede(B_comp);

    // KO derivative portion
    // clang-format off
#include "RHSINCLUDE_taskflow_kodiss.cpp"
    // clang-format on

    tf::Task koDissTask =
        taskflow.composed_of(flow_kodiss).name("KODissComputation");

    a_comp.precede(koDissTask);
    b_comp.precede(koDissTask);
    gt_comp.precede(koDissTask);
    chi_comp.precede(koDissTask);
    At_comp.precede(koDissTask);
    K_comp.precede(koDissTask);
    Gt_comp.precede(koDissTask);
    B_comp.precede(koDissTask);

    bssn::timer::t_rhs.start();
    executor.run(taskflow).wait();
    bssn::timer::t_rhs.stop();

#if 0
        for (unsigned int m = 0; m < 24; m++) {
            std::cout<<"  || dtu("<<m<<")|| = "<<normLInfty(unzipVarsRHS[m] + offset, n)<<std::endl;
        }
#endif
}

void bssnrhs_sep_taskflow_varsOnly(
    double **unzipVarsRHS, const double **uZipVars, const unsigned int &offset,
    const double *pmin, const double *pmax, const unsigned int *sz,
    const unsigned int &bflag, const unsigned int threadID,
    const unsigned int maxBlockSize) {
    const double *alpha = &uZipVars[VAR::U_ALPHA][offset];
    const double *chi = &uZipVars[VAR::U_CHI][offset];
    const double *K = &uZipVars[VAR::U_K][offset];
    const double *gt0 = &uZipVars[VAR::U_SYMGT0][offset];
    const double *gt1 = &uZipVars[VAR::U_SYMGT1][offset];
    const double *gt2 = &uZipVars[VAR::U_SYMGT2][offset];
    const double *gt3 = &uZipVars[VAR::U_SYMGT3][offset];
    const double *gt4 = &uZipVars[VAR::U_SYMGT4][offset];
    const double *gt5 = &uZipVars[VAR::U_SYMGT5][offset];
    const double *beta0 = &uZipVars[VAR::U_BETA0][offset];
    const double *beta1 = &uZipVars[VAR::U_BETA1][offset];
    const double *beta2 = &uZipVars[VAR::U_BETA2][offset];
    const double *At0 = &uZipVars[VAR::U_SYMAT0][offset];
    const double *At1 = &uZipVars[VAR::U_SYMAT1][offset];
    const double *At2 = &uZipVars[VAR::U_SYMAT2][offset];
    const double *At3 = &uZipVars[VAR::U_SYMAT3][offset];
    const double *At4 = &uZipVars[VAR::U_SYMAT4][offset];
    const double *At5 = &uZipVars[VAR::U_SYMAT5][offset];
    const double *Gt0 = &uZipVars[VAR::U_GT0][offset];
    const double *Gt1 = &uZipVars[VAR::U_GT1][offset];
    const double *Gt2 = &uZipVars[VAR::U_GT2][offset];
    const double *B0 = &uZipVars[VAR::U_B0][offset];
    const double *B1 = &uZipVars[VAR::U_B1][offset];
    const double *B2 = &uZipVars[VAR::U_B2][offset];

    double *const a_rhs = &unzipVarsRHS[VAR::U_ALPHA][offset];
    double *const chi_rhs = &unzipVarsRHS[VAR::U_CHI][offset];
    double *const K_rhs = &unzipVarsRHS[VAR::U_K][offset];
    double *const gt_rhs00 = &unzipVarsRHS[VAR::U_SYMGT0][offset];
    double *const gt_rhs01 = &unzipVarsRHS[VAR::U_SYMGT1][offset];
    double *const gt_rhs02 = &unzipVarsRHS[VAR::U_SYMGT2][offset];
    double *const gt_rhs11 = &unzipVarsRHS[VAR::U_SYMGT3][offset];
    double *const gt_rhs12 = &unzipVarsRHS[VAR::U_SYMGT4][offset];
    double *const gt_rhs22 = &unzipVarsRHS[VAR::U_SYMGT5][offset];
    double *const b_rhs0 = &unzipVarsRHS[VAR::U_BETA0][offset];
    double *const b_rhs1 = &unzipVarsRHS[VAR::U_BETA1][offset];
    double *const b_rhs2 = &unzipVarsRHS[VAR::U_BETA2][offset];
    double *const At_rhs00 = &unzipVarsRHS[VAR::U_SYMAT0][offset];
    double *const At_rhs01 = &unzipVarsRHS[VAR::U_SYMAT1][offset];
    double *const At_rhs02 = &unzipVarsRHS[VAR::U_SYMAT2][offset];
    double *const At_rhs11 = &unzipVarsRHS[VAR::U_SYMAT3][offset];
    double *const At_rhs12 = &unzipVarsRHS[VAR::U_SYMAT4][offset];
    double *const At_rhs22 = &unzipVarsRHS[VAR::U_SYMAT5][offset];
    double *const Gt_rhs0 = &unzipVarsRHS[VAR::U_GT0][offset];
    double *const Gt_rhs1 = &unzipVarsRHS[VAR::U_GT1][offset];
    double *const Gt_rhs2 = &unzipVarsRHS[VAR::U_GT2][offset];
    double *const B_rhs0 = &unzipVarsRHS[VAR::U_B0][offset];
    double *const B_rhs1 = &unzipVarsRHS[VAR::U_B1][offset];
    double *const B_rhs2 = &unzipVarsRHS[VAR::U_B2][offset];

    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];

    double hx = (pmax[0] - pmin[0]) / (nx - 1);
    double hy = (pmax[1] - pmin[1]) / (ny - 1);
    double hz = (pmax[2] - pmin[2]) / (nz - 1);

    const unsigned int lambda[4] = {BSSN_LAMBDA[0], BSSN_LAMBDA[1],
                                    BSSN_LAMBDA[2], BSSN_LAMBDA[3]};
    const double lambda_f[2] = {BSSN_LAMBDA_F[0], BSSN_LAMBDA_F[1]};
    const double eta_power[2] = {BSSN_ETA_POWER[0], BSSN_ETA_POWER[1]};

    int idx[3];
    const unsigned int PW = bssn::BSSN_PADDING_WIDTH;
    unsigned int n = sz[0] * sz[1] * sz[2];
    mem::memory_pool<double> *__mem_pool = &BSSN_MEM_POOL;

    const unsigned int BLK_SZ = n;
    double *const deriv_base = bssn::BSSN_DERIV_WORKSPACE;

    const double sigma = KO_DISS_SIGMA;
    DendroRegister unsigned int pp;

#if 0
    bssn::timer::t_deriv.start();

// clang-format off
#include "bssnrhs_evar_derivs.h"
#include "bssnrhs_derivs.h"
#include "bssnrhs_derivs_adv.h"
    // clang-format on

    bssn::timer::t_deriv.stop();
#endif

#if 0
    tf::Taskflow taskflow;

    auto [a_comp, b_comp, gt_comp, chi_comp, At_comp, K_comp, Gt_comp, B_comp] =
        taskflow.emplace(
            [&]() {
                rhs_a_only(unzipVarsRHS, uZipVars, offset, pmin, pmax, sz,
                           bflag);
            },
            [&]() {
                rhs_b_only(unzipVarsRHS, uZipVars, offset, pmin, pmax, sz,
                           bflag);
            },
            [&]() {
                rhs_gt_only(unzipVarsRHS, uZipVars, offset, pmin, pmax, sz,
                            bflag);
            },
            [&]() {
                rhs_chi_only(unzipVarsRHS, uZipVars, offset, pmin, pmax, sz,
                             bflag);
            },
            [&]() {
                rhs_At_only(unzipVarsRHS, uZipVars, offset, pmin, pmax, sz,
                            bflag);
            },
            [&]() {
                rhs_K_only(unzipVarsRHS, uZipVars, offset, pmin, pmax, sz,
                           bflag);
            },
            [&]() {
                rhs_Gt_only(unzipVarsRHS, uZipVars, offset, pmin, pmax, sz,
                            bflag);
            },
            [&]() {
                rhs_B_only(unzipVarsRHS, uZipVars, offset, pmin, pmax, sz,
                           bflag);
            });
#endif

    // update the data that the graph will need
    tfdendro::data.unzipVarsRHS = unzipVarsRHS;
    tfdendro::data.uZipVars = (double **)uZipVars;
    tfdendro::data.offset = offset;
    tfdendro::data.pmin = (double *)pmin;
    tfdendro::data.pmax = (double *)pmax;
    tfdendro::data.sz = (unsigned int *)sz;
    tfdendro::data.bflag = bflag;

    tfdendro::data.BLK_SZ = BLK_SZ;
    tfdendro::data.deriv_base = deriv_base;
    tfdendro::data.hx = hx;
    tfdendro::data.hy = hy;
    tfdendro::data.hz = hz;

    tfdendro::data.sigma = sigma;
    tfdendro::data.PW = PW;
    tfdendro::data.nx = nx;
    tfdendro::data.ny = ny;
    tfdendro::data.nz = nz;

    tfdendro::data.alpha = (double *)alpha;
    tfdendro::data.chi = (double *)chi;
    tfdendro::data.K = (double *)K;
    tfdendro::data.gt0 = (double *)gt0;
    tfdendro::data.gt1 = (double *)gt1;
    tfdendro::data.gt2 = (double *)gt2;
    tfdendro::data.gt3 = (double *)gt3;
    tfdendro::data.gt4 = (double *)gt4;
    tfdendro::data.gt5 = (double *)gt5;
    tfdendro::data.beta0 = (double *)beta0;
    tfdendro::data.beta1 = (double *)beta1;
    tfdendro::data.beta2 = (double *)beta2;
    tfdendro::data.At0 = (double *)At0;
    tfdendro::data.At1 = (double *)At1;
    tfdendro::data.At2 = (double *)At2;
    tfdendro::data.At3 = (double *)At3;
    tfdendro::data.At4 = (double *)At4;
    tfdendro::data.At5 = (double *)At5;
    tfdendro::data.Gt0 = (double *)Gt0;
    tfdendro::data.Gt1 = (double *)Gt1;
    tfdendro::data.Gt2 = (double *)Gt2;
    tfdendro::data.B0 = (double *)B0;
    tfdendro::data.B1 = (double *)B1;
    tfdendro::data.B2 = (double *)B2;

    tfdendro::data.a_rhs = (double *)a_rhs;
    tfdendro::data.chi_rhs = (double *)chi_rhs;
    tfdendro::data.K_rhs = (double *)K_rhs;
    tfdendro::data.gt_rhs00 = (double *)gt_rhs00;
    tfdendro::data.gt_rhs01 = (double *)gt_rhs01;
    tfdendro::data.gt_rhs02 = (double *)gt_rhs02;
    tfdendro::data.gt_rhs11 = (double *)gt_rhs11;
    tfdendro::data.gt_rhs12 = (double *)gt_rhs12;
    tfdendro::data.gt_rhs22 = (double *)gt_rhs22;
    tfdendro::data.b_rhs0 = (double *)b_rhs0;
    tfdendro::data.b_rhs1 = (double *)b_rhs1;
    tfdendro::data.b_rhs2 = (double *)b_rhs2;
    tfdendro::data.At_rhs00 = (double *)At_rhs00;
    tfdendro::data.At_rhs01 = (double *)At_rhs01;
    tfdendro::data.At_rhs02 = (double *)At_rhs02;
    tfdendro::data.At_rhs11 = (double *)At_rhs11;
    tfdendro::data.At_rhs12 = (double *)At_rhs12;
    tfdendro::data.At_rhs22 = (double *)At_rhs22;
    tfdendro::data.Gt_rhs0 = (double *)Gt_rhs0;
    tfdendro::data.Gt_rhs1 = (double *)Gt_rhs1;
    tfdendro::data.Gt_rhs2 = (double *)Gt_rhs2;
    tfdendro::data.B_rhs0 = (double *)B_rhs0;
    tfdendro::data.B_rhs1 = (double *)B_rhs1;
    tfdendro::data.B_rhs2 = (double *)B_rhs2;

    bssn::timer::t_rhs.start();
    tfdendro::executor.run(tfdendro::taskflow).wait();
    bssn::timer::t_rhs.stop();

#if 0
    if (bflag != 0) {
        bssn::timer::t_bdyc.start();

        bssn_bcs(a_rhs, alpha, grad_0_alpha, grad_1_alpha, grad_2_alpha, pmin,
                 pmax, 1.0, 1.0, sz, bflag);
        bssn_bcs(chi_rhs, chi, grad_0_chi, grad_1_chi, grad_2_chi, pmin, pmax,
                 1.0, 1.0, sz, bflag);
        bssn_bcs(K_rhs, K, grad_0_K, grad_1_K, grad_2_K, pmin, pmax, 1.0, 0.0,
                 sz, bflag);

        bssn_bcs(b_rhs0, beta0, grad_0_beta0, grad_1_beta0, grad_2_beta0, pmin,
                 pmax, 1.0, 0.0, sz, bflag);
        bssn_bcs(b_rhs1, beta1, grad_0_beta1, grad_1_beta1, grad_2_beta1, pmin,
                 pmax, 1.0, 0.0, sz, bflag);
        bssn_bcs(b_rhs2, beta2, grad_0_beta2, grad_1_beta2, grad_2_beta2, pmin,
                 pmax, 1.0, 0.0, sz, bflag);

        bssn_bcs(Gt_rhs0, Gt0, grad_0_Gt0, grad_1_Gt0, grad_2_Gt0, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(Gt_rhs1, Gt1, grad_0_Gt1, grad_1_Gt1, grad_2_Gt1, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(Gt_rhs2, Gt2, grad_0_Gt2, grad_1_Gt2, grad_2_Gt2, pmin, pmax,
                 2.0, 0.0, sz, bflag);

        bssn_bcs(B_rhs0, B0, grad_0_B0, grad_1_B0, grad_2_B0, pmin, pmax, 1.0,
                 0.0, sz, bflag);
        bssn_bcs(B_rhs1, B1, grad_0_B1, grad_1_B1, grad_2_B1, pmin, pmax, 1.0,
                 0.0, sz, bflag);
        bssn_bcs(B_rhs2, B2, grad_0_B2, grad_1_B2, grad_2_B2, pmin, pmax, 1.0,
                 0.0, sz, bflag);

        bssn_bcs(At_rhs00, At0, grad_0_At0, grad_1_At0, grad_2_At0, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs01, At1, grad_0_At1, grad_1_At1, grad_2_At1, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs02, At2, grad_0_At2, grad_1_At2, grad_2_At2, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs11, At3, grad_0_At3, grad_1_At3, grad_2_At3, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs12, At4, grad_0_At4, grad_1_At4, grad_2_At4, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs22, At5, grad_0_At5, grad_1_At5, grad_2_At5, pmin, pmax,
                 2.0, 0.0, sz, bflag);

        bssn_bcs(gt_rhs00, gt0, grad_0_gt0, grad_1_gt0, grad_2_gt0, pmin, pmax,
                 1.0, 1.0, sz, bflag);
        bssn_bcs(gt_rhs01, gt1, grad_0_gt1, grad_1_gt1, grad_2_gt1, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(gt_rhs02, gt2, grad_0_gt2, grad_1_gt2, grad_2_gt2, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(gt_rhs11, gt3, grad_0_gt3, grad_1_gt3, grad_2_gt3, pmin, pmax,
                 1.0, 1.0, sz, bflag);
        bssn_bcs(gt_rhs12, gt4, grad_0_gt4, grad_1_gt4, grad_2_gt4, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(gt_rhs22, gt5, grad_0_gt5, grad_1_gt5, grad_2_gt5, pmin, pmax,
                 1.0, 1.0, sz, bflag);

        bssn::timer::t_bdyc.stop();
    }
#endif

#if 0
#include "bssnrhs_evar_derivs.h"

    bssn::timer::t_deriv.start();
#include "bssnrhs_ko_derivs.h"
    bssn::timer::t_deriv.stop();

    bssn::timer::t_rhs.start();

    for (unsigned int k = PW; k < nz - PW; k++) {
        for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = PW; i < nx - PW; i++) {
                pp = i + nx * (j + ny * k);

                a_rhs[pp] += sigma * (grad_0_alpha[pp] + grad_1_alpha[pp] +
                                      grad_2_alpha[pp]);
                b_rhs0[pp] += sigma * (grad_0_beta0[pp] + grad_1_beta0[pp] +
                                       grad_2_beta0[pp]);
                b_rhs1[pp] += sigma * (grad_0_beta1[pp] + grad_1_beta1[pp] +
                                       grad_2_beta1[pp]);
                b_rhs2[pp] += sigma * (grad_0_beta2[pp] + grad_1_beta2[pp] +
                                       grad_2_beta2[pp]);

                gt_rhs00[pp] +=
                    sigma * (grad_0_gt0[pp] + grad_1_gt0[pp] + grad_2_gt0[pp]);
                gt_rhs01[pp] +=
                    sigma * (grad_0_gt1[pp] + grad_1_gt1[pp] + grad_2_gt1[pp]);
                gt_rhs02[pp] +=
                    sigma * (grad_0_gt2[pp] + grad_1_gt2[pp] + grad_2_gt2[pp]);
                gt_rhs11[pp] +=
                    sigma * (grad_0_gt3[pp] + grad_1_gt3[pp] + grad_2_gt3[pp]);
                gt_rhs12[pp] +=
                    sigma * (grad_0_gt4[pp] + grad_1_gt4[pp] + grad_2_gt4[pp]);
                gt_rhs22[pp] +=
                    sigma * (grad_0_gt5[pp] + grad_1_gt5[pp] + grad_2_gt5[pp]);

                chi_rhs[pp] +=
                    sigma * (grad_0_chi[pp] + grad_1_chi[pp] + grad_2_chi[pp]);

                At_rhs00[pp] +=
                    sigma * (grad_0_At0[pp] + grad_1_At0[pp] + grad_2_At0[pp]);
                At_rhs01[pp] +=
                    sigma * (grad_0_At1[pp] + grad_1_At1[pp] + grad_2_At1[pp]);
                At_rhs02[pp] +=
                    sigma * (grad_0_At2[pp] + grad_1_At2[pp] + grad_2_At2[pp]);
                At_rhs11[pp] +=
                    sigma * (grad_0_At3[pp] + grad_1_At3[pp] + grad_2_At3[pp]);
                At_rhs12[pp] +=
                    sigma * (grad_0_At4[pp] + grad_1_At4[pp] + grad_2_At4[pp]);
                At_rhs22[pp] +=
                    sigma * (grad_0_At5[pp] + grad_1_At5[pp] + grad_2_At5[pp]);

                K_rhs[pp] +=
                    sigma * (grad_0_K[pp] + grad_1_K[pp] + grad_2_K[pp]);

                Gt_rhs0[pp] +=
                    sigma * (grad_0_Gt0[pp] + grad_1_Gt0[pp] + grad_2_Gt0[pp]);
                Gt_rhs1[pp] +=
                    sigma * (grad_0_Gt1[pp] + grad_1_Gt1[pp] + grad_2_Gt1[pp]);
                Gt_rhs2[pp] +=
                    sigma * (grad_0_Gt2[pp] + grad_1_Gt2[pp] + grad_2_Gt2[pp]);

                B_rhs0[pp] +=
                    sigma * (grad_0_B0[pp] + grad_1_B0[pp] + grad_2_B0[pp]);
                B_rhs1[pp] +=
                    sigma * (grad_0_B1[pp] + grad_1_B1[pp] + grad_2_B1[pp]);
                B_rhs2[pp] +=
                    sigma * (grad_0_B2[pp] + grad_1_B2[pp] + grad_2_B2[pp]);
            }
        }
    }

    bssn::timer::t_rhs.stop();

    bssn::timer::t_deriv.start();
    bssn::timer::t_deriv.stop();
#endif

#if 0
        for (unsigned int m = 0; m < 24; m++) {
            std::cout<<"  || dtu("<<m<<")|| = "<<normLInfty(unzipVarsRHS[m] + offset, n)<<std::endl;
        }
#endif
}
