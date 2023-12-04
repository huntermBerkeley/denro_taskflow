auto [gt_rhs00_kodiss_g, gt_rhs01_kodiss_g, gt_rhs02_kodiss_g,
      gt_rhs11_kodiss_g, gt_rhs12_kodiss_g, gt_rhs22_kodiss_g,
      At_rhs00_kodiss_g, At_rhs01_kodiss_g, At_rhs02_kodiss_g,
      At_rhs11_kodiss_g, At_rhs12_kodiss_g, At_rhs22_kodiss_g, a_rhs_kodiss_g,
      b_rhs0_kodiss_g, b_rhs1_kodiss_g, b_rhs2_kodiss_g, chi_rhs_kodiss_g,
      Gt_rhs0_kodiss_g, Gt_rhs1_kodiss_g, Gt_rhs2_kodiss_g, K_rhs_kodiss_g,
      B_rhs0_kodiss_g, B_rhs1_kodiss_g, B_rhs2_kodiss_g] =
    flow_kodiss.emplace(
        [&]() {
            ko_deriv_x(grad_0_gt0, gt0, hx, sz, bflag);
            ko_deriv_y(grad_1_gt0, gt0, hy, sz, bflag);
            ko_deriv_z(grad_2_gt0, gt0, hz, sz, bflag);
            for (unsigned int k = PW; k < nz - PW; k++) {
                for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
                    for (unsigned int i = PW; i < nx - PW; i++) {
                        DendroRegister unsigned int pp = i + nx * (j + ny * k);
                        gt_rhs00[pp] +=
                            sigma *
                            (grad_0_gt0[pp] + grad_1_gt0[pp] + grad_2_gt0[pp]);
                    }
                }
            }
        },
        [&]() {
            ko_deriv_x(grad_0_gt1, gt1, hx, sz, bflag);
            ko_deriv_y(grad_1_gt1, gt1, hy, sz, bflag);
            ko_deriv_z(grad_2_gt1, gt1, hz, sz, bflag);
            for (unsigned int k = PW; k < nz - PW; k++) {
                for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
                    for (unsigned int i = PW; i < nx - PW; i++) {
                        DendroRegister unsigned int pp = i + nx * (j + ny * k);
                        gt_rhs01[pp] +=
                            sigma *
                            (grad_0_gt1[pp] + grad_1_gt1[pp] + grad_2_gt1[pp]);
                    }
                }
            }
        },
        [&]() {
            ko_deriv_x(grad_0_gt2, gt2, hx, sz, bflag);
            ko_deriv_y(grad_1_gt2, gt2, hy, sz, bflag);
            ko_deriv_z(grad_2_gt2, gt2, hz, sz, bflag);
            for (unsigned int k = PW; k < nz - PW; k++) {
                for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
                    for (unsigned int i = PW; i < nx - PW; i++) {
                        DendroRegister unsigned int pp = i + nx * (j + ny * k);
                        gt_rhs02[pp] +=
                            sigma *
                            (grad_0_gt2[pp] + grad_1_gt2[pp] + grad_2_gt2[pp]);
                    }
                }
            }
        },
        [&]() {
            ko_deriv_x(grad_0_gt3, gt3, hx, sz, bflag);
            ko_deriv_y(grad_1_gt3, gt3, hy, sz, bflag);
            ko_deriv_z(grad_2_gt3, gt3, hz, sz, bflag);
            for (unsigned int k = PW; k < nz - PW; k++) {
                for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
                    for (unsigned int i = PW; i < nx - PW; i++) {
                        DendroRegister unsigned int pp = i + nx * (j + ny * k);
                        gt_rhs11[pp] +=
                            sigma *
                            (grad_0_gt3[pp] + grad_1_gt3[pp] + grad_2_gt3[pp]);
                    }
                }
            }
        },
        [&]() {
            ko_deriv_x(grad_0_gt4, gt4, hx, sz, bflag);
            ko_deriv_y(grad_1_gt4, gt4, hy, sz, bflag);
            ko_deriv_z(grad_2_gt4, gt4, hz, sz, bflag);
            for (unsigned int k = PW; k < nz - PW; k++) {
                for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
                    for (unsigned int i = PW; i < nx - PW; i++) {
                        DendroRegister unsigned int pp = i + nx * (j + ny * k);
                        gt_rhs12[pp] +=
                            sigma *
                            (grad_0_gt4[pp] + grad_1_gt4[pp] + grad_2_gt4[pp]);
                    }
                }
            }
        },
        [&]() {
            ko_deriv_x(grad_0_gt5, gt5, hx, sz, bflag);
            ko_deriv_y(grad_1_gt5, gt5, hy, sz, bflag);
            ko_deriv_z(grad_2_gt5, gt5, hz, sz, bflag);
            for (unsigned int k = PW; k < nz - PW; k++) {
                for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
                    for (unsigned int i = PW; i < nx - PW; i++) {
                        DendroRegister unsigned int pp = i + nx * (j + ny * k);
                        gt_rhs22[pp] +=
                            sigma *
                            (grad_0_gt5[pp] + grad_1_gt5[pp] + grad_2_gt5[pp]);
                    }
                }
            }
        },
        [&]() {
            ko_deriv_x(grad_0_At0, At0, hx, sz, bflag);
            ko_deriv_y(grad_1_At0, At0, hy, sz, bflag);
            ko_deriv_z(grad_2_At0, At0, hz, sz, bflag);
            for (unsigned int k = PW; k < nz - PW; k++) {
                for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
                    for (unsigned int i = PW; i < nx - PW; i++) {
                        DendroRegister unsigned int pp = i + nx * (j + ny * k);
                        At_rhs00[pp] +=
                            sigma *
                            (grad_0_At0[pp] + grad_1_At0[pp] + grad_2_At0[pp]);
                    }
                }
            }
        },
        [&]() {
            ko_deriv_x(grad_0_At1, At1, hx, sz, bflag);
            ko_deriv_y(grad_1_At1, At1, hy, sz, bflag);
            ko_deriv_z(grad_2_At1, At1, hz, sz, bflag);
            for (unsigned int k = PW; k < nz - PW; k++) {
                for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
                    for (unsigned int i = PW; i < nx - PW; i++) {
                        DendroRegister unsigned int pp = i + nx * (j + ny * k);
                        At_rhs01[pp] +=
                            sigma *
                            (grad_0_At1[pp] + grad_1_At1[pp] + grad_2_At1[pp]);
                    }
                }
            }
        },
        [&]() {
            ko_deriv_x(grad_0_At2, At2, hx, sz, bflag);
            ko_deriv_y(grad_1_At2, At2, hy, sz, bflag);
            ko_deriv_z(grad_2_At2, At2, hz, sz, bflag);
            for (unsigned int k = PW; k < nz - PW; k++) {
                for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
                    for (unsigned int i = PW; i < nx - PW; i++) {
                        DendroRegister unsigned int pp = i + nx * (j + ny * k);
                        At_rhs02[pp] +=
                            sigma *
                            (grad_0_At2[pp] + grad_1_At2[pp] + grad_2_At2[pp]);
                    }
                }
            }
        },
        [&]() {
            ko_deriv_x(grad_0_At3, At3, hx, sz, bflag);
            ko_deriv_y(grad_1_At3, At3, hy, sz, bflag);
            ko_deriv_z(grad_2_At3, At3, hz, sz, bflag);
            for (unsigned int k = PW; k < nz - PW; k++) {
                for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
                    for (unsigned int i = PW; i < nx - PW; i++) {
                        DendroRegister unsigned int pp = i + nx * (j + ny * k);
                        At_rhs11[pp] +=
                            sigma *
                            (grad_0_At3[pp] + grad_1_At3[pp] + grad_2_At3[pp]);
                    }
                }
            }
        },
        [&]() {
            ko_deriv_x(grad_0_At4, At4, hx, sz, bflag);
            ko_deriv_y(grad_1_At4, At4, hy, sz, bflag);
            ko_deriv_z(grad_2_At4, At4, hz, sz, bflag);
            for (unsigned int k = PW; k < nz - PW; k++) {
                for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
                    for (unsigned int i = PW; i < nx - PW; i++) {
                        DendroRegister unsigned int pp = i + nx * (j + ny * k);
                        At_rhs12[pp] +=
                            sigma *
                            (grad_0_At4[pp] + grad_1_At4[pp] + grad_2_At4[pp]);
                    }
                }
            }
        },
        [&]() {
            ko_deriv_x(grad_0_At5, At5, hx, sz, bflag);
            ko_deriv_y(grad_1_At5, At5, hy, sz, bflag);
            ko_deriv_z(grad_2_At5, At5, hz, sz, bflag);
            for (unsigned int k = PW; k < nz - PW; k++) {
                for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
                    for (unsigned int i = PW; i < nx - PW; i++) {
                        DendroRegister unsigned int pp = i + nx * (j + ny * k);
                        At_rhs22[pp] +=
                            sigma *
                            (grad_0_At5[pp] + grad_1_At5[pp] + grad_2_At5[pp]);
                    }
                }
            }
        },
        [&]() {
            ko_deriv_x(grad_0_alpha, alpha, hx, sz, bflag);
            ko_deriv_y(grad_1_alpha, alpha, hy, sz, bflag);
            ko_deriv_z(grad_2_alpha, alpha, hz, sz, bflag);
            for (unsigned int k = PW; k < nz - PW; k++) {
                for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
                    for (unsigned int i = PW; i < nx - PW; i++) {
                        DendroRegister unsigned int pp = i + nx * (j + ny * k);
                        a_rhs[pp] +=
                            sigma * (grad_0_alpha[pp] + grad_1_alpha[pp] +
                                     grad_2_alpha[pp]);
                    }
                }
            }
        },
        [&]() {
            ko_deriv_x(grad_0_beta0, beta0, hx, sz, bflag);
            ko_deriv_y(grad_1_beta0, beta0, hy, sz, bflag);
            ko_deriv_z(grad_2_beta0, beta0, hz, sz, bflag);
            for (unsigned int k = PW; k < nz - PW; k++) {
                for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
                    for (unsigned int i = PW; i < nx - PW; i++) {
                        DendroRegister unsigned int pp = i + nx * (j + ny * k);
                        b_rhs0[pp] +=
                            sigma * (grad_0_beta0[pp] + grad_1_beta0[pp] +
                                     grad_2_beta0[pp]);
                    }
                }
            }
        },
        [&]() {
            ko_deriv_x(grad_0_beta1, beta1, hx, sz, bflag);
            ko_deriv_y(grad_1_beta1, beta1, hy, sz, bflag);
            ko_deriv_z(grad_2_beta1, beta1, hz, sz, bflag);
            for (unsigned int k = PW; k < nz - PW; k++) {
                for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
                    for (unsigned int i = PW; i < nx - PW; i++) {
                        DendroRegister unsigned int pp = i + nx * (j + ny * k);
                        b_rhs1[pp] +=
                            sigma * (grad_0_beta1[pp] + grad_1_beta1[pp] +
                                     grad_2_beta1[pp]);
                    }
                }
            }
        },
        [&]() {
            ko_deriv_x(grad_0_beta2, beta2, hx, sz, bflag);
            ko_deriv_y(grad_1_beta2, beta2, hy, sz, bflag);
            ko_deriv_z(grad_2_beta2, beta2, hz, sz, bflag);
            for (unsigned int k = PW; k < nz - PW; k++) {
                for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
                    for (unsigned int i = PW; i < nx - PW; i++) {
                        DendroRegister unsigned int pp = i + nx * (j + ny * k);
                        b_rhs2[pp] +=
                            sigma * (grad_0_beta2[pp] + grad_1_beta2[pp] +
                                     grad_2_beta2[pp]);
                    }
                }
            }
        },
        [&]() {
            ko_deriv_x(grad_0_chi, chi, hx, sz, bflag);
            ko_deriv_y(grad_1_chi, chi, hy, sz, bflag);
            ko_deriv_z(grad_2_chi, chi, hz, sz, bflag);
            for (unsigned int k = PW; k < nz - PW; k++) {
                for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
                    for (unsigned int i = PW; i < nx - PW; i++) {
                        DendroRegister unsigned int pp = i + nx * (j + ny * k);
                        chi_rhs[pp] +=
                            sigma *
                            (grad_0_chi[pp] + grad_1_chi[pp] + grad_2_chi[pp]);
                    }
                }
            }
        },
        [&]() {
            ko_deriv_x(grad_0_Gt0, Gt0, hx, sz, bflag);
            ko_deriv_y(grad_1_Gt0, Gt0, hy, sz, bflag);
            ko_deriv_z(grad_2_Gt0, Gt0, hz, sz, bflag);
            for (unsigned int k = PW; k < nz - PW; k++) {
                for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
                    for (unsigned int i = PW; i < nx - PW; i++) {
                        DendroRegister unsigned int pp = i + nx * (j + ny * k);
                        Gt_rhs0[pp] +=
                            sigma *
                            (grad_0_Gt0[pp] + grad_1_Gt0[pp] + grad_2_Gt0[pp]);
                    }
                }
            }
        },
        [&]() {
            ko_deriv_x(grad_0_Gt1, Gt1, hx, sz, bflag);
            ko_deriv_y(grad_1_Gt1, Gt1, hy, sz, bflag);
            ko_deriv_z(grad_2_Gt1, Gt1, hz, sz, bflag);
            for (unsigned int k = PW; k < nz - PW; k++) {
                for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
                    for (unsigned int i = PW; i < nx - PW; i++) {
                        DendroRegister unsigned int pp = i + nx * (j + ny * k);
                        Gt_rhs1[pp] +=
                            sigma *
                            (grad_0_Gt1[pp] + grad_1_Gt1[pp] + grad_2_Gt1[pp]);
                    }
                }
            }
        },
        [&]() {
            ko_deriv_x(grad_0_Gt2, Gt2, hx, sz, bflag);
            ko_deriv_y(grad_1_Gt2, Gt2, hy, sz, bflag);
            ko_deriv_z(grad_2_Gt2, Gt2, hz, sz, bflag);
            for (unsigned int k = PW; k < nz - PW; k++) {
                for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
                    for (unsigned int i = PW; i < nx - PW; i++) {
                        DendroRegister unsigned int pp = i + nx * (j + ny * k);
                        Gt_rhs2[pp] +=
                            sigma *
                            (grad_0_Gt2[pp] + grad_1_Gt2[pp] + grad_2_Gt2[pp]);
                    }
                }
            }
        },
        [&]() {
            ko_deriv_x(grad_0_K, K, hx, sz, bflag);
            ko_deriv_y(grad_1_K, K, hy, sz, bflag);
            ko_deriv_z(grad_2_K, K, hz, sz, bflag);
            for (unsigned int k = PW; k < nz - PW; k++) {
                for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
                    for (unsigned int i = PW; i < nx - PW; i++) {
                        DendroRegister unsigned int pp = i + nx * (j + ny * k);
                        K_rhs[pp] += sigma * (grad_0_K[pp] + grad_1_K[pp] +
                                              grad_2_K[pp]);
                    }
                }
            }
        },
        [&]() {
            ko_deriv_x(grad_0_B0, B0, hx, sz, bflag);
            ko_deriv_y(grad_1_B0, B0, hy, sz, bflag);
            ko_deriv_z(grad_2_B0, B0, hz, sz, bflag);
            for (unsigned int k = PW; k < nz - PW; k++) {
                for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
                    for (unsigned int i = PW; i < nx - PW; i++) {
                        DendroRegister unsigned int pp = i + nx * (j + ny * k);
                        B_rhs0[pp] += sigma * (grad_0_B0[pp] + grad_1_B0[pp] +
                                               grad_2_B0[pp]);
                    }
                }
            }
        },
        [&]() {
            ko_deriv_x(grad_0_B1, B1, hx, sz, bflag);
            ko_deriv_y(grad_1_B1, B1, hy, sz, bflag);
            ko_deriv_z(grad_2_B1, B1, hz, sz, bflag);
            for (unsigned int k = PW; k < nz - PW; k++) {
                for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
                    for (unsigned int i = PW; i < nx - PW; i++) {
                        DendroRegister unsigned int pp = i + nx * (j + ny * k);
                        B_rhs1[pp] += sigma * (grad_0_B1[pp] + grad_1_B1[pp] +
                                               grad_2_B1[pp]);
                    }
                }
            }
        },
        [&]() {
            ko_deriv_x(grad_0_B2, B2, hx, sz, bflag);
            ko_deriv_y(grad_1_B2, B2, hy, sz, bflag);
            ko_deriv_z(grad_2_B2, B2, hz, sz, bflag);
            for (unsigned int k = PW; k < nz - PW; k++) {
                for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
                    for (unsigned int i = PW; i < nx - PW; i++) {
                        DendroRegister unsigned int pp = i + nx * (j + ny * k);
                        B_rhs2[pp] += sigma * (grad_0_B2[pp] + grad_1_B2[pp] +
                                               grad_2_B2[pp]);
                    }
                }
            }
        });
