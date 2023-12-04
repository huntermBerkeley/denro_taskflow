auto gt_rhs00_kodiss_g =
    tfdendro::flow_kodiss.placeholder().name("gt_rhs00_kodiss_g");
auto gt_rhs01_kodiss_g =
    tfdendro::flow_kodiss.placeholder().name("gt_rhs01_kodiss_g");
auto gt_rhs02_kodiss_g =
    tfdendro::flow_kodiss.placeholder().name("gt_rhs02_kodiss_g");
auto gt_rhs11_kodiss_g =
    tfdendro::flow_kodiss.placeholder().name("gt_rhs11_kodiss_g");
auto gt_rhs12_kodiss_g =
    tfdendro::flow_kodiss.placeholder().name("gt_rhs12_kodiss_g");
auto gt_rhs22_kodiss_g =
    tfdendro::flow_kodiss.placeholder().name("gt_rhs22_kodiss_g");
auto At_rhs00_kodiss_g =
    tfdendro::flow_kodiss.placeholder().name("At_rhs00_kodiss_g");
auto At_rhs01_kodiss_g =
    tfdendro::flow_kodiss.placeholder().name("At_rhs01_kodiss_g");
auto At_rhs02_kodiss_g =
    tfdendro::flow_kodiss.placeholder().name("At_rhs02_kodiss_g");
auto At_rhs11_kodiss_g =
    tfdendro::flow_kodiss.placeholder().name("At_rhs11_kodiss_g");
auto At_rhs12_kodiss_g =
    tfdendro::flow_kodiss.placeholder().name("At_rhs12_kodiss_g");
auto At_rhs22_kodiss_g =
    tfdendro::flow_kodiss.placeholder().name("At_rhs22_kodiss_g");
auto a_rhs_kodiss_g =
    tfdendro::flow_kodiss.placeholder().name("a_rhs_kodiss_g");
auto b_rhs0_kodiss_g =
    tfdendro::flow_kodiss.placeholder().name("b_rhs0_kodiss_g");
auto b_rhs1_kodiss_g =
    tfdendro::flow_kodiss.placeholder().name("b_rhs1_kodiss_g");
auto b_rhs2_kodiss_g =
    tfdendro::flow_kodiss.placeholder().name("b_rhs2_kodiss_g");
auto chi_rhs_kodiss_g =
    tfdendro::flow_kodiss.placeholder().name("chi_rhs_kodiss_g");
auto Gt_rhs0_kodiss_g =
    tfdendro::flow_kodiss.placeholder().name("Gt_rhs0_kodiss_g");
auto Gt_rhs1_kodiss_g =
    tfdendro::flow_kodiss.placeholder().name("Gt_rhs1_kodiss_g");
auto Gt_rhs2_kodiss_g =
    tfdendro::flow_kodiss.placeholder().name("Gt_rhs2_kodiss_g");
auto K_rhs_kodiss_g =
    tfdendro::flow_kodiss.placeholder().name("K_rhs_kodiss_g");
auto B_rhs0_kodiss_g =
    tfdendro::flow_kodiss.placeholder().name("B_rhs0_kodiss_g");
auto B_rhs1_kodiss_g =
    tfdendro::flow_kodiss.placeholder().name("B_rhs1_kodiss_g");
auto B_rhs2_kodiss_g =
    tfdendro::flow_kodiss.placeholder().name("B_rhs2_kodiss_g");

gt_rhs00_kodiss_g.data(&tfdendro::data).work([gt_rhs00_kodiss_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(gt_rhs00_kodiss_g.data());

    double* grad_0_gt0 = d.deriv_base + 36 * d.BLK_SZ;
    double* grad_1_gt0 = d.deriv_base + 37 * d.BLK_SZ;
    double* grad_2_gt0 = d.deriv_base + 38 * d.BLK_SZ;

    ko_deriv_x(grad_0_gt0, d.gt0, d.hx, d.sz, d.bflag);
    ko_deriv_y(grad_1_gt0, d.gt0, d.hy, d.sz, d.bflag);
    ko_deriv_z(grad_2_gt0, d.gt0, d.hz, d.sz, d.bflag);

    for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
        for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                DendroRegister unsigned int pp = i + d.nx * (j + d.ny * k);
                d.gt_rhs00[pp] += d.sigma * (grad_0_gt0[pp] + grad_1_gt0[pp] +
                                             grad_2_gt0[pp]);
            }
        }
    }
});
gt_rhs01_kodiss_g.data(&tfdendro::data).work([gt_rhs01_kodiss_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(gt_rhs01_kodiss_g.data());

    double* grad_0_gt1 = d.deriv_base + 39 * d.BLK_SZ;
    double* grad_1_gt1 = d.deriv_base + 40 * d.BLK_SZ;
    double* grad_2_gt1 = d.deriv_base + 41 * d.BLK_SZ;

    ko_deriv_x(grad_0_gt1, d.gt1, d.hx, d.sz, d.bflag);
    ko_deriv_y(grad_1_gt1, d.gt1, d.hy, d.sz, d.bflag);
    ko_deriv_z(grad_2_gt1, d.gt1, d.hz, d.sz, d.bflag);

    for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
        for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                DendroRegister unsigned int pp = i + d.nx * (j + d.ny * k);
                d.gt_rhs01[pp] += d.sigma * (grad_0_gt1[pp] + grad_1_gt1[pp] +
                                             grad_2_gt1[pp]);
            }
        }
    }
});
gt_rhs02_kodiss_g.data(&tfdendro::data).work([gt_rhs02_kodiss_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(gt_rhs02_kodiss_g.data());

    double* grad_0_gt2 = d.deriv_base + 42 * d.BLK_SZ;
    double* grad_1_gt2 = d.deriv_base + 43 * d.BLK_SZ;
    double* grad_2_gt2 = d.deriv_base + 44 * d.BLK_SZ;

    ko_deriv_x(grad_0_gt2, d.gt2, d.hx, d.sz, d.bflag);
    ko_deriv_y(grad_1_gt2, d.gt2, d.hy, d.sz, d.bflag);
    ko_deriv_z(grad_2_gt2, d.gt2, d.hz, d.sz, d.bflag);

    for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
        for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                DendroRegister unsigned int pp = i + d.nx * (j + d.ny * k);
                d.gt_rhs02[pp] += d.sigma * (grad_0_gt2[pp] + grad_1_gt2[pp] +
                                             grad_2_gt2[pp]);
            }
        }
    }
});
gt_rhs11_kodiss_g.data(&tfdendro::data).work([gt_rhs11_kodiss_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(gt_rhs11_kodiss_g.data());

    double* grad_0_gt3 = d.deriv_base + 45 * d.BLK_SZ;
    double* grad_1_gt3 = d.deriv_base + 46 * d.BLK_SZ;
    double* grad_2_gt3 = d.deriv_base + 47 * d.BLK_SZ;

    ko_deriv_x(grad_0_gt3, d.gt3, d.hx, d.sz, d.bflag);
    ko_deriv_y(grad_1_gt3, d.gt3, d.hy, d.sz, d.bflag);
    ko_deriv_z(grad_2_gt3, d.gt3, d.hz, d.sz, d.bflag);

    for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
        for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                DendroRegister unsigned int pp = i + d.nx * (j + d.ny * k);
                d.gt_rhs11[pp] += d.sigma * (grad_0_gt3[pp] + grad_1_gt3[pp] +
                                             grad_2_gt3[pp]);
            }
        }
    }
});
gt_rhs12_kodiss_g.data(&tfdendro::data).work([gt_rhs12_kodiss_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(gt_rhs12_kodiss_g.data());

    double* grad_0_gt4 = d.deriv_base + 48 * d.BLK_SZ;
    double* grad_1_gt4 = d.deriv_base + 49 * d.BLK_SZ;
    double* grad_2_gt4 = d.deriv_base + 50 * d.BLK_SZ;

    ko_deriv_x(grad_0_gt4, d.gt4, d.hx, d.sz, d.bflag);
    ko_deriv_y(grad_1_gt4, d.gt4, d.hy, d.sz, d.bflag);
    ko_deriv_z(grad_2_gt4, d.gt4, d.hz, d.sz, d.bflag);

    for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
        for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                DendroRegister unsigned int pp = i + d.nx * (j + d.ny * k);
                d.gt_rhs12[pp] += d.sigma * (grad_0_gt4[pp] + grad_1_gt4[pp] +
                                             grad_2_gt4[pp]);
            }
        }
    }
});
gt_rhs22_kodiss_g.data(&tfdendro::data).work([gt_rhs22_kodiss_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(gt_rhs22_kodiss_g.data());

    double* grad_0_gt5 = d.deriv_base + 51 * d.BLK_SZ;
    double* grad_1_gt5 = d.deriv_base + 52 * d.BLK_SZ;
    double* grad_2_gt5 = d.deriv_base + 53 * d.BLK_SZ;

    ko_deriv_x(grad_0_gt5, d.gt5, d.hx, d.sz, d.bflag);
    ko_deriv_y(grad_1_gt5, d.gt5, d.hy, d.sz, d.bflag);
    ko_deriv_z(grad_2_gt5, d.gt5, d.hz, d.sz, d.bflag);

    for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
        for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                DendroRegister unsigned int pp = i + d.nx * (j + d.ny * k);
                d.gt_rhs22[pp] += d.sigma * (grad_0_gt5[pp] + grad_1_gt5[pp] +
                                             grad_2_gt5[pp]);
            }
        }
    }
});
At_rhs00_kodiss_g.data(&tfdendro::data).work([At_rhs00_kodiss_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(At_rhs00_kodiss_g.data());

    double* grad_0_At0 = d.deriv_base + 54 * d.BLK_SZ;
    double* grad_1_At0 = d.deriv_base + 55 * d.BLK_SZ;
    double* grad_2_At0 = d.deriv_base + 56 * d.BLK_SZ;

    ko_deriv_x(grad_0_At0, d.At0, d.hx, d.sz, d.bflag);
    ko_deriv_y(grad_1_At0, d.At0, d.hy, d.sz, d.bflag);
    ko_deriv_z(grad_2_At0, d.At0, d.hz, d.sz, d.bflag);

    for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
        for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                DendroRegister unsigned int pp = i + d.nx * (j + d.ny * k);
                d.At_rhs00[pp] += d.sigma * (grad_0_At0[pp] + grad_1_At0[pp] +
                                             grad_2_At0[pp]);
            }
        }
    }
});
At_rhs01_kodiss_g.data(&tfdendro::data).work([At_rhs01_kodiss_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(At_rhs01_kodiss_g.data());

    double* grad_0_At1 = d.deriv_base + 57 * d.BLK_SZ;
    double* grad_1_At1 = d.deriv_base + 58 * d.BLK_SZ;
    double* grad_2_At1 = d.deriv_base + 59 * d.BLK_SZ;

    ko_deriv_x(grad_0_At1, d.At1, d.hx, d.sz, d.bflag);
    ko_deriv_y(grad_1_At1, d.At1, d.hy, d.sz, d.bflag);
    ko_deriv_z(grad_2_At1, d.At1, d.hz, d.sz, d.bflag);

    for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
        for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                DendroRegister unsigned int pp = i + d.nx * (j + d.ny * k);
                d.At_rhs01[pp] += d.sigma * (grad_0_At1[pp] + grad_1_At1[pp] +
                                             grad_2_At1[pp]);
            }
        }
    }
});
At_rhs02_kodiss_g.data(&tfdendro::data).work([At_rhs02_kodiss_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(At_rhs02_kodiss_g.data());

    double* grad_0_At2 = d.deriv_base + 60 * d.BLK_SZ;
    double* grad_1_At2 = d.deriv_base + 61 * d.BLK_SZ;
    double* grad_2_At2 = d.deriv_base + 62 * d.BLK_SZ;

    ko_deriv_x(grad_0_At2, d.At2, d.hx, d.sz, d.bflag);
    ko_deriv_y(grad_1_At2, d.At2, d.hy, d.sz, d.bflag);
    ko_deriv_z(grad_2_At2, d.At2, d.hz, d.sz, d.bflag);

    for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
        for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                DendroRegister unsigned int pp = i + d.nx * (j + d.ny * k);
                d.At_rhs02[pp] += d.sigma * (grad_0_At2[pp] + grad_1_At2[pp] +
                                             grad_2_At2[pp]);
            }
        }
    }
});
At_rhs11_kodiss_g.data(&tfdendro::data).work([At_rhs11_kodiss_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(At_rhs11_kodiss_g.data());

    double* grad_0_At3 = d.deriv_base + 63 * d.BLK_SZ;
    double* grad_1_At3 = d.deriv_base + 64 * d.BLK_SZ;
    double* grad_2_At3 = d.deriv_base + 65 * d.BLK_SZ;

    ko_deriv_x(grad_0_At3, d.At3, d.hx, d.sz, d.bflag);
    ko_deriv_y(grad_1_At3, d.At3, d.hy, d.sz, d.bflag);
    ko_deriv_z(grad_2_At3, d.At3, d.hz, d.sz, d.bflag);

    for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
        for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                DendroRegister unsigned int pp = i + d.nx * (j + d.ny * k);
                d.At_rhs11[pp] += d.sigma * (grad_0_At3[pp] + grad_1_At3[pp] +
                                             grad_2_At3[pp]);
            }
        }
    }
});
At_rhs12_kodiss_g.data(&tfdendro::data).work([At_rhs12_kodiss_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(At_rhs12_kodiss_g.data());

    double* grad_0_At4 = d.deriv_base + 66 * d.BLK_SZ;
    double* grad_1_At4 = d.deriv_base + 67 * d.BLK_SZ;
    double* grad_2_At4 = d.deriv_base + 68 * d.BLK_SZ;

    ko_deriv_x(grad_0_At4, d.At4, d.hx, d.sz, d.bflag);
    ko_deriv_y(grad_1_At4, d.At4, d.hy, d.sz, d.bflag);
    ko_deriv_z(grad_2_At4, d.At4, d.hz, d.sz, d.bflag);

    for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
        for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                DendroRegister unsigned int pp = i + d.nx * (j + d.ny * k);
                d.At_rhs12[pp] += d.sigma * (grad_0_At4[pp] + grad_1_At4[pp] +
                                             grad_2_At4[pp]);
            }
        }
    }
});
At_rhs22_kodiss_g.data(&tfdendro::data).work([At_rhs22_kodiss_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(At_rhs22_kodiss_g.data());

    double* grad_0_At5 = d.deriv_base + 69 * d.BLK_SZ;
    double* grad_1_At5 = d.deriv_base + 70 * d.BLK_SZ;
    double* grad_2_At5 = d.deriv_base + 71 * d.BLK_SZ;

    ko_deriv_x(grad_0_At5, d.At5, d.hx, d.sz, d.bflag);
    ko_deriv_y(grad_1_At5, d.At5, d.hy, d.sz, d.bflag);
    ko_deriv_z(grad_2_At5, d.At5, d.hz, d.sz, d.bflag);

    for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
        for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                DendroRegister unsigned int pp = i + d.nx * (j + d.ny * k);
                d.At_rhs22[pp] += d.sigma * (grad_0_At5[pp] + grad_1_At5[pp] +
                                             grad_2_At5[pp]);
            }
        }
    }
});
a_rhs_kodiss_g.data(&tfdendro::data).work([a_rhs_kodiss_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(a_rhs_kodiss_g.data());

    double* grad_0_alpha = d.deriv_base + 0 * d.BLK_SZ;
    double* grad_1_alpha = d.deriv_base + 1 * d.BLK_SZ;
    double* grad_2_alpha = d.deriv_base + 2 * d.BLK_SZ;

    ko_deriv_x(grad_0_alpha, d.alpha, d.hx, d.sz, d.bflag);
    ko_deriv_y(grad_1_alpha, d.alpha, d.hy, d.sz, d.bflag);
    ko_deriv_z(grad_2_alpha, d.alpha, d.hz, d.sz, d.bflag);

    for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
        for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                DendroRegister unsigned int pp = i + d.nx * (j + d.ny * k);
                d.a_rhs[pp] += d.sigma * (grad_0_alpha[pp] + grad_1_alpha[pp] +
                                          grad_2_alpha[pp]);
            }
        }
    }
});
b_rhs0_kodiss_g.data(&tfdendro::data).work([b_rhs0_kodiss_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(b_rhs0_kodiss_g.data());

    double* grad_0_beta0 = d.deriv_base + 3 * d.BLK_SZ;
    double* grad_1_beta0 = d.deriv_base + 4 * d.BLK_SZ;
    double* grad_2_beta0 = d.deriv_base + 5 * d.BLK_SZ;

    ko_deriv_x(grad_0_beta0, d.beta0, d.hx, d.sz, d.bflag);
    ko_deriv_y(grad_1_beta0, d.beta0, d.hy, d.sz, d.bflag);
    ko_deriv_z(grad_2_beta0, d.beta0, d.hz, d.sz, d.bflag);

    for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
        for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                DendroRegister unsigned int pp = i + d.nx * (j + d.ny * k);
                d.b_rhs0[pp] += d.sigma * (grad_0_beta0[pp] + grad_1_beta0[pp] +
                                           grad_2_beta0[pp]);
            }
        }
    }
});
b_rhs1_kodiss_g.data(&tfdendro::data).work([b_rhs1_kodiss_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(b_rhs1_kodiss_g.data());

    double* grad_0_beta1 = d.deriv_base + 6 * d.BLK_SZ;
    double* grad_1_beta1 = d.deriv_base + 7 * d.BLK_SZ;
    double* grad_2_beta1 = d.deriv_base + 8 * d.BLK_SZ;

    ko_deriv_x(grad_0_beta1, d.beta1, d.hx, d.sz, d.bflag);
    ko_deriv_y(grad_1_beta1, d.beta1, d.hy, d.sz, d.bflag);
    ko_deriv_z(grad_2_beta1, d.beta1, d.hz, d.sz, d.bflag);

    for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
        for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                DendroRegister unsigned int pp = i + d.nx * (j + d.ny * k);
                d.b_rhs1[pp] += d.sigma * (grad_0_beta1[pp] + grad_1_beta1[pp] +
                                           grad_2_beta1[pp]);
            }
        }
    }
});
b_rhs2_kodiss_g.data(&tfdendro::data).work([b_rhs2_kodiss_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(b_rhs2_kodiss_g.data());

    double* grad_0_beta2 = d.deriv_base + 9 * d.BLK_SZ;
    double* grad_1_beta2 = d.deriv_base + 10 * d.BLK_SZ;
    double* grad_2_beta2 = d.deriv_base + 11 * d.BLK_SZ;

    ko_deriv_x(grad_0_beta2, d.beta2, d.hx, d.sz, d.bflag);
    ko_deriv_y(grad_1_beta2, d.beta2, d.hy, d.sz, d.bflag);
    ko_deriv_z(grad_2_beta2, d.beta2, d.hz, d.sz, d.bflag);

    for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
        for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                DendroRegister unsigned int pp = i + d.nx * (j + d.ny * k);
                d.b_rhs2[pp] += d.sigma * (grad_0_beta2[pp] + grad_1_beta2[pp] +
                                           grad_2_beta2[pp]);
            }
        }
    }
});
chi_rhs_kodiss_g.data(&tfdendro::data).work([chi_rhs_kodiss_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(chi_rhs_kodiss_g.data());

    double* grad_0_chi = d.deriv_base + 21 * d.BLK_SZ;
    double* grad_1_chi = d.deriv_base + 22 * d.BLK_SZ;
    double* grad_2_chi = d.deriv_base + 23 * d.BLK_SZ;

    ko_deriv_x(grad_0_chi, d.chi, d.hx, d.sz, d.bflag);
    ko_deriv_y(grad_1_chi, d.chi, d.hy, d.sz, d.bflag);
    ko_deriv_z(grad_2_chi, d.chi, d.hz, d.sz, d.bflag);

    for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
        for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                DendroRegister unsigned int pp = i + d.nx * (j + d.ny * k);
                d.chi_rhs[pp] += d.sigma * (grad_0_chi[pp] + grad_1_chi[pp] +
                                            grad_2_chi[pp]);
            }
        }
    }
});
Gt_rhs0_kodiss_g.data(&tfdendro::data).work([Gt_rhs0_kodiss_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(Gt_rhs0_kodiss_g.data());

    double* grad_0_Gt0 = d.deriv_base + 24 * d.BLK_SZ;
    double* grad_1_Gt0 = d.deriv_base + 25 * d.BLK_SZ;
    double* grad_2_Gt0 = d.deriv_base + 26 * d.BLK_SZ;

    ko_deriv_x(grad_0_Gt0, d.Gt0, d.hx, d.sz, d.bflag);
    ko_deriv_y(grad_1_Gt0, d.Gt0, d.hy, d.sz, d.bflag);
    ko_deriv_z(grad_2_Gt0, d.Gt0, d.hz, d.sz, d.bflag);

    for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
        for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                DendroRegister unsigned int pp = i + d.nx * (j + d.ny * k);
                d.Gt_rhs0[pp] += d.sigma * (grad_0_Gt0[pp] + grad_1_Gt0[pp] +
                                            grad_2_Gt0[pp]);
            }
        }
    }
});
Gt_rhs1_kodiss_g.data(&tfdendro::data).work([Gt_rhs1_kodiss_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(Gt_rhs1_kodiss_g.data());

    double* grad_0_Gt1 = d.deriv_base + 27 * d.BLK_SZ;
    double* grad_1_Gt1 = d.deriv_base + 28 * d.BLK_SZ;
    double* grad_2_Gt1 = d.deriv_base + 29 * d.BLK_SZ;

    ko_deriv_x(grad_0_Gt1, d.Gt1, d.hx, d.sz, d.bflag);
    ko_deriv_y(grad_1_Gt1, d.Gt1, d.hy, d.sz, d.bflag);
    ko_deriv_z(grad_2_Gt1, d.Gt1, d.hz, d.sz, d.bflag);

    for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
        for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                DendroRegister unsigned int pp = i + d.nx * (j + d.ny * k);
                d.Gt_rhs1[pp] += d.sigma * (grad_0_Gt1[pp] + grad_1_Gt1[pp] +
                                            grad_2_Gt1[pp]);
            }
        }
    }
});
Gt_rhs2_kodiss_g.data(&tfdendro::data).work([Gt_rhs2_kodiss_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(Gt_rhs2_kodiss_g.data());

    double* grad_0_Gt2 = d.deriv_base + 30 * d.BLK_SZ;
    double* grad_1_Gt2 = d.deriv_base + 31 * d.BLK_SZ;
    double* grad_2_Gt2 = d.deriv_base + 32 * d.BLK_SZ;

    ko_deriv_x(grad_0_Gt2, d.Gt2, d.hx, d.sz, d.bflag);
    ko_deriv_y(grad_1_Gt2, d.Gt2, d.hy, d.sz, d.bflag);
    ko_deriv_z(grad_2_Gt2, d.Gt2, d.hz, d.sz, d.bflag);

    for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
        for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                DendroRegister unsigned int pp = i + d.nx * (j + d.ny * k);
                d.Gt_rhs2[pp] += d.sigma * (grad_0_Gt2[pp] + grad_1_Gt2[pp] +
                                            grad_2_Gt2[pp]);
            }
        }
    }
});
K_rhs_kodiss_g.data(&tfdendro::data).work([K_rhs_kodiss_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(K_rhs_kodiss_g.data());

    double* grad_0_K = d.deriv_base + 33 * d.BLK_SZ;
    double* grad_1_K = d.deriv_base + 34 * d.BLK_SZ;
    double* grad_2_K = d.deriv_base + 35 * d.BLK_SZ;

    ko_deriv_x(grad_0_K, d.K, d.hx, d.sz, d.bflag);
    ko_deriv_y(grad_1_K, d.K, d.hy, d.sz, d.bflag);
    ko_deriv_z(grad_2_K, d.K, d.hz, d.sz, d.bflag);

    for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
        for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                DendroRegister unsigned int pp = i + d.nx * (j + d.ny * k);
                d.K_rhs[pp] +=
                    d.sigma * (grad_0_K[pp] + grad_1_K[pp] + grad_2_K[pp]);
            }
        }
    }
});
B_rhs0_kodiss_g.data(&tfdendro::data).work([B_rhs0_kodiss_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(B_rhs0_kodiss_g.data());

    double* grad_0_B0 = d.deriv_base + 12 * d.BLK_SZ;
    double* grad_1_B0 = d.deriv_base + 13 * d.BLK_SZ;
    double* grad_2_B0 = d.deriv_base + 14 * d.BLK_SZ;

    ko_deriv_x(grad_0_B0, d.B0, d.hx, d.sz, d.bflag);
    ko_deriv_y(grad_1_B0, d.B0, d.hy, d.sz, d.bflag);
    ko_deriv_z(grad_2_B0, d.B0, d.hz, d.sz, d.bflag);

    for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
        for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                DendroRegister unsigned int pp = i + d.nx * (j + d.ny * k);
                d.B_rhs0[pp] +=
                    d.sigma * (grad_0_B0[pp] + grad_1_B0[pp] + grad_2_B0[pp]);
            }
        }
    }
});
B_rhs1_kodiss_g.data(&tfdendro::data).work([B_rhs1_kodiss_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(B_rhs1_kodiss_g.data());

    double* grad_0_B1 = d.deriv_base + 15 * d.BLK_SZ;
    double* grad_1_B1 = d.deriv_base + 16 * d.BLK_SZ;
    double* grad_2_B1 = d.deriv_base + 17 * d.BLK_SZ;

    ko_deriv_x(grad_0_B1, d.B1, d.hx, d.sz, d.bflag);
    ko_deriv_y(grad_1_B1, d.B1, d.hy, d.sz, d.bflag);
    ko_deriv_z(grad_2_B1, d.B1, d.hz, d.sz, d.bflag);

    for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
        for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                DendroRegister unsigned int pp = i + d.nx * (j + d.ny * k);
                d.B_rhs1[pp] +=
                    d.sigma * (grad_0_B1[pp] + grad_1_B1[pp] + grad_2_B1[pp]);
            }
        }
    }
});
B_rhs2_kodiss_g.data(&tfdendro::data).work([B_rhs2_kodiss_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(B_rhs2_kodiss_g.data());

    double* grad_0_B2 = d.deriv_base + 18 * d.BLK_SZ;
    double* grad_1_B2 = d.deriv_base + 19 * d.BLK_SZ;
    double* grad_2_B2 = d.deriv_base + 20 * d.BLK_SZ;

    ko_deriv_x(grad_0_B2, d.B2, d.hx, d.sz, d.bflag);
    ko_deriv_y(grad_1_B2, d.B2, d.hy, d.sz, d.bflag);
    ko_deriv_z(grad_2_B2, d.B2, d.hz, d.sz, d.bflag);

    for (unsigned int k = d.PW; k < d.nz - d.PW; k++) {
        for (unsigned int j = d.PW; j < d.ny - d.PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = d.PW; i < d.nx - d.PW; i++) {
                DendroRegister unsigned int pp = i + d.nx * (j + d.ny * k);
                d.B_rhs2[pp] +=
                    d.sigma * (grad_0_B2[pp] + grad_1_B2[pp] + grad_2_B2[pp]);
            }
        }
    }
});
