auto grad_0_alpha_g = tfdendro::flow_deriv.placeholder().name("grad_0_alpha_g");
auto grad_1_alpha_g = tfdendro::flow_deriv.placeholder().name("grad_1_alpha_g");
auto grad_2_alpha_g = tfdendro::flow_deriv.placeholder().name("grad_2_alpha_g");
auto grad_0_beta0_g = tfdendro::flow_deriv.placeholder().name("grad_0_beta0_g");
auto grad_1_beta0_g = tfdendro::flow_deriv.placeholder().name("grad_1_beta0_g");
auto grad_2_beta0_g = tfdendro::flow_deriv.placeholder().name("grad_2_beta0_g");
auto grad_0_beta1_g = tfdendro::flow_deriv.placeholder().name("grad_0_beta1_g");
auto grad_1_beta1_g = tfdendro::flow_deriv.placeholder().name("grad_1_beta1_g");
auto grad_2_beta1_g = tfdendro::flow_deriv.placeholder().name("grad_2_beta1_g");
auto grad_0_beta2_g = tfdendro::flow_deriv.placeholder().name("grad_0_beta2_g");
auto grad_1_beta2_g = tfdendro::flow_deriv.placeholder().name("grad_1_beta2_g");
auto grad_2_beta2_g = tfdendro::flow_deriv.placeholder().name("grad_2_beta2_g");
auto grad_0_B0_g = tfdendro::flow_deriv.placeholder().name("grad_0_B0_g");
auto grad_1_B0_g = tfdendro::flow_deriv.placeholder().name("grad_1_B0_g");
auto grad_2_B0_g = tfdendro::flow_deriv.placeholder().name("grad_2_B0_g");
auto grad_0_B1_g = tfdendro::flow_deriv.placeholder().name("grad_0_B1_g");
auto grad_1_B1_g = tfdendro::flow_deriv.placeholder().name("grad_1_B1_g");
auto grad_2_B1_g = tfdendro::flow_deriv.placeholder().name("grad_2_B1_g");
auto grad_0_B2_g = tfdendro::flow_deriv.placeholder().name("grad_0_B2_g");
auto grad_1_B2_g = tfdendro::flow_deriv.placeholder().name("grad_1_B2_g");
auto grad_2_B2_g = tfdendro::flow_deriv.placeholder().name("grad_2_B2_g");
auto grad_0_chi_g = tfdendro::flow_deriv.placeholder().name("grad_0_chi_g");
auto grad_1_chi_g = tfdendro::flow_deriv.placeholder().name("grad_1_chi_g");
auto grad_2_chi_g = tfdendro::flow_deriv.placeholder().name("grad_2_chi_g");
auto grad_0_Gt0_g = tfdendro::flow_deriv.placeholder().name("grad_0_Gt0_g");
auto grad_1_Gt0_g = tfdendro::flow_deriv.placeholder().name("grad_1_Gt0_g");
auto grad_2_Gt0_g = tfdendro::flow_deriv.placeholder().name("grad_2_Gt0_g");
auto grad_0_Gt1_g = tfdendro::flow_deriv.placeholder().name("grad_0_Gt1_g");
auto grad_1_Gt1_g = tfdendro::flow_deriv.placeholder().name("grad_1_Gt1_g");
auto grad_2_Gt1_g = tfdendro::flow_deriv.placeholder().name("grad_2_Gt1_g");
auto grad_0_Gt2_g = tfdendro::flow_deriv.placeholder().name("grad_0_Gt2_g");
auto grad_1_Gt2_g = tfdendro::flow_deriv.placeholder().name("grad_1_Gt2_g");
auto grad_2_Gt2_g = tfdendro::flow_deriv.placeholder().name("grad_2_Gt2_g");
auto grad_0_K_g = tfdendro::flow_deriv.placeholder().name("grad_0_K_g");
auto grad_1_K_g = tfdendro::flow_deriv.placeholder().name("grad_1_K_g");
auto grad_2_K_g = tfdendro::flow_deriv.placeholder().name("grad_2_K_g");
auto grad_0_gt0_g = tfdendro::flow_deriv.placeholder().name("grad_0_gt0_g");
auto grad_1_gt0_g = tfdendro::flow_deriv.placeholder().name("grad_1_gt0_g");
auto grad_2_gt0_g = tfdendro::flow_deriv.placeholder().name("grad_2_gt0_g");
auto grad_0_gt1_g = tfdendro::flow_deriv.placeholder().name("grad_0_gt1_g");
auto grad_1_gt1_g = tfdendro::flow_deriv.placeholder().name("grad_1_gt1_g");
auto grad_2_gt1_g = tfdendro::flow_deriv.placeholder().name("grad_2_gt1_g");
auto grad_0_gt2_g = tfdendro::flow_deriv.placeholder().name("grad_0_gt2_g");
auto grad_1_gt2_g = tfdendro::flow_deriv.placeholder().name("grad_1_gt2_g");
auto grad_2_gt2_g = tfdendro::flow_deriv.placeholder().name("grad_2_gt2_g");
auto grad_0_gt3_g = tfdendro::flow_deriv.placeholder().name("grad_0_gt3_g");
auto grad_1_gt3_g = tfdendro::flow_deriv.placeholder().name("grad_1_gt3_g");
auto grad_2_gt3_g = tfdendro::flow_deriv.placeholder().name("grad_2_gt3_g");
auto grad_0_gt4_g = tfdendro::flow_deriv.placeholder().name("grad_0_gt4_g");
auto grad_1_gt4_g = tfdendro::flow_deriv.placeholder().name("grad_1_gt4_g");
auto grad_2_gt4_g = tfdendro::flow_deriv.placeholder().name("grad_2_gt4_g");
auto grad_0_gt5_g = tfdendro::flow_deriv.placeholder().name("grad_0_gt5_g");
auto grad_1_gt5_g = tfdendro::flow_deriv.placeholder().name("grad_1_gt5_g");
auto grad_2_gt5_g = tfdendro::flow_deriv.placeholder().name("grad_2_gt5_g");
auto grad_0_At0_g = tfdendro::flow_deriv.placeholder().name("grad_0_At0_g");
auto grad_1_At0_g = tfdendro::flow_deriv.placeholder().name("grad_1_At0_g");
auto grad_2_At0_g = tfdendro::flow_deriv.placeholder().name("grad_2_At0_g");
auto grad_0_At1_g = tfdendro::flow_deriv.placeholder().name("grad_0_At1_g");
auto grad_1_At1_g = tfdendro::flow_deriv.placeholder().name("grad_1_At1_g");
auto grad_2_At1_g = tfdendro::flow_deriv.placeholder().name("grad_2_At1_g");
auto grad_0_At2_g = tfdendro::flow_deriv.placeholder().name("grad_0_At2_g");
auto grad_1_At2_g = tfdendro::flow_deriv.placeholder().name("grad_1_At2_g");
auto grad_2_At2_g = tfdendro::flow_deriv.placeholder().name("grad_2_At2_g");
auto grad_0_At3_g = tfdendro::flow_deriv.placeholder().name("grad_0_At3_g");
auto grad_1_At3_g = tfdendro::flow_deriv.placeholder().name("grad_1_At3_g");
auto grad_2_At3_g = tfdendro::flow_deriv.placeholder().name("grad_2_At3_g");
auto grad_0_At4_g = tfdendro::flow_deriv.placeholder().name("grad_0_At4_g");
auto grad_1_At4_g = tfdendro::flow_deriv.placeholder().name("grad_1_At4_g");
auto grad_2_At4_g = tfdendro::flow_deriv.placeholder().name("grad_2_At4_g");
auto grad_0_At5_g = tfdendro::flow_deriv.placeholder().name("grad_0_At5_g");
auto grad_1_At5_g = tfdendro::flow_deriv.placeholder().name("grad_1_At5_g");
auto grad_2_At5_g = tfdendro::flow_deriv.placeholder().name("grad_2_At5_g");
auto grad2_0_0_gt0_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_0_gt0_g");
auto grad2_0_1_gt0_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_1_gt0_g");
auto grad2_0_2_gt0_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_2_gt0_g");
auto grad2_1_1_gt0_g =
    tfdendro::flow_deriv.placeholder().name("grad2_1_1_gt0_g");
auto grad2_1_2_gt0_g =
    tfdendro::flow_deriv.placeholder().name("grad2_1_2_gt0_g");
auto grad2_2_2_gt0_g =
    tfdendro::flow_deriv.placeholder().name("grad2_2_2_gt0_g");
auto grad2_0_0_gt1_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_0_gt1_g");
auto grad2_0_1_gt1_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_1_gt1_g");
auto grad2_0_2_gt1_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_2_gt1_g");
auto grad2_1_1_gt1_g =
    tfdendro::flow_deriv.placeholder().name("grad2_1_1_gt1_g");
auto grad2_1_2_gt1_g =
    tfdendro::flow_deriv.placeholder().name("grad2_1_2_gt1_g");
auto grad2_2_2_gt1_g =
    tfdendro::flow_deriv.placeholder().name("grad2_2_2_gt1_g");
auto grad2_0_0_gt2_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_0_gt2_g");
auto grad2_0_1_gt2_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_1_gt2_g");
auto grad2_0_2_gt2_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_2_gt2_g");
auto grad2_1_1_gt2_g =
    tfdendro::flow_deriv.placeholder().name("grad2_1_1_gt2_g");
auto grad2_1_2_gt2_g =
    tfdendro::flow_deriv.placeholder().name("grad2_1_2_gt2_g");
auto grad2_2_2_gt2_g =
    tfdendro::flow_deriv.placeholder().name("grad2_2_2_gt2_g");
auto grad2_0_0_gt3_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_0_gt3_g");
auto grad2_0_1_gt3_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_1_gt3_g");
auto grad2_0_2_gt3_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_2_gt3_g");
auto grad2_1_1_gt3_g =
    tfdendro::flow_deriv.placeholder().name("grad2_1_1_gt3_g");
auto grad2_1_2_gt3_g =
    tfdendro::flow_deriv.placeholder().name("grad2_1_2_gt3_g");
auto grad2_2_2_gt3_g =
    tfdendro::flow_deriv.placeholder().name("grad2_2_2_gt3_g");
auto grad2_0_0_gt4_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_0_gt4_g");
auto grad2_0_1_gt4_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_1_gt4_g");
auto grad2_0_2_gt4_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_2_gt4_g");
auto grad2_1_1_gt4_g =
    tfdendro::flow_deriv.placeholder().name("grad2_1_1_gt4_g");
auto grad2_1_2_gt4_g =
    tfdendro::flow_deriv.placeholder().name("grad2_1_2_gt4_g");
auto grad2_2_2_gt4_g =
    tfdendro::flow_deriv.placeholder().name("grad2_2_2_gt4_g");
auto grad2_0_0_gt5_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_0_gt5_g");
auto grad2_0_1_gt5_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_1_gt5_g");
auto grad2_0_2_gt5_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_2_gt5_g");
auto grad2_1_1_gt5_g =
    tfdendro::flow_deriv.placeholder().name("grad2_1_1_gt5_g");
auto grad2_1_2_gt5_g =
    tfdendro::flow_deriv.placeholder().name("grad2_1_2_gt5_g");
auto grad2_2_2_gt5_g =
    tfdendro::flow_deriv.placeholder().name("grad2_2_2_gt5_g");
auto grad2_0_0_chi_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_0_chi_g");
auto grad2_0_1_chi_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_1_chi_g");
auto grad2_0_2_chi_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_2_chi_g");
auto grad2_1_1_chi_g =
    tfdendro::flow_deriv.placeholder().name("grad2_1_1_chi_g");
auto grad2_1_2_chi_g =
    tfdendro::flow_deriv.placeholder().name("grad2_1_2_chi_g");
auto grad2_2_2_chi_g =
    tfdendro::flow_deriv.placeholder().name("grad2_2_2_chi_g");
auto grad2_0_0_alpha_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_0_alpha_g");
auto grad2_0_1_alpha_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_1_alpha_g");
auto grad2_0_2_alpha_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_2_alpha_g");
auto grad2_1_1_alpha_g =
    tfdendro::flow_deriv.placeholder().name("grad2_1_1_alpha_g");
auto grad2_1_2_alpha_g =
    tfdendro::flow_deriv.placeholder().name("grad2_1_2_alpha_g");
auto grad2_2_2_alpha_g =
    tfdendro::flow_deriv.placeholder().name("grad2_2_2_alpha_g");
auto grad2_0_0_beta0_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_0_beta0_g");
auto grad2_0_1_beta0_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_1_beta0_g");
auto grad2_0_2_beta0_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_2_beta0_g");
auto grad2_1_1_beta0_g =
    tfdendro::flow_deriv.placeholder().name("grad2_1_1_beta0_g");
auto grad2_1_2_beta0_g =
    tfdendro::flow_deriv.placeholder().name("grad2_1_2_beta0_g");
auto grad2_2_2_beta0_g =
    tfdendro::flow_deriv.placeholder().name("grad2_2_2_beta0_g");
auto grad2_0_0_beta1_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_0_beta1_g");
auto grad2_0_1_beta1_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_1_beta1_g");
auto grad2_0_2_beta1_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_2_beta1_g");
auto grad2_1_1_beta1_g =
    tfdendro::flow_deriv.placeholder().name("grad2_1_1_beta1_g");
auto grad2_1_2_beta1_g =
    tfdendro::flow_deriv.placeholder().name("grad2_1_2_beta1_g");
auto grad2_2_2_beta1_g =
    tfdendro::flow_deriv.placeholder().name("grad2_2_2_beta1_g");
auto grad2_0_0_beta2_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_0_beta2_g");
auto grad2_0_1_beta2_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_1_beta2_g");
auto grad2_0_2_beta2_g =
    tfdendro::flow_deriv.placeholder().name("grad2_0_2_beta2_g");
auto grad2_1_1_beta2_g =
    tfdendro::flow_deriv.placeholder().name("grad2_1_1_beta2_g");
auto grad2_1_2_beta2_g =
    tfdendro::flow_deriv.placeholder().name("grad2_1_2_beta2_g");
auto grad2_2_2_beta2_g =
    tfdendro::flow_deriv.placeholder().name("grad2_2_2_beta2_g");

grad_0_alpha_g.data(&tfdendro::data).work([grad_0_alpha_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_0_alpha_g.data());
    double* grad_0_alpha = d.deriv_base + 0 * d.BLK_SZ;
    deriv_x(grad_0_alpha, d.alpha, d.hx, d.sz, d.bflag);
});
grad_1_alpha_g.data(&tfdendro::data).work([grad_1_alpha_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_1_alpha_g.data());
    double* grad_1_alpha = d.deriv_base + 1 * d.BLK_SZ;
    deriv_y(grad_1_alpha, d.alpha, d.hy, d.sz, d.bflag);
});
grad_2_alpha_g.data(&tfdendro::data).work([grad_2_alpha_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_2_alpha_g.data());
    double* grad_2_alpha = d.deriv_base + 2 * d.BLK_SZ;
    deriv_z(grad_2_alpha, d.alpha, d.hz, d.sz, d.bflag);
});
grad_0_beta0_g.data(&tfdendro::data).work([grad_0_beta0_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_0_beta0_g.data());
    double* grad_0_beta0 = d.deriv_base + 3 * d.BLK_SZ;
    deriv_x(grad_0_beta0, d.beta0, d.hx, d.sz, d.bflag);
});
grad_1_beta0_g.data(&tfdendro::data).work([grad_1_beta0_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_1_beta0_g.data());
    double* grad_1_beta0 = d.deriv_base + 4 * d.BLK_SZ;
    deriv_y(grad_1_beta0, d.beta0, d.hy, d.sz, d.bflag);
});
grad_2_beta0_g.data(&tfdendro::data).work([grad_2_beta0_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_2_beta0_g.data());
    double* grad_2_beta0 = d.deriv_base + 5 * d.BLK_SZ;
    deriv_z(grad_2_beta0, d.beta0, d.hz, d.sz, d.bflag);
});
grad_0_beta1_g.data(&tfdendro::data).work([grad_0_beta1_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_0_beta1_g.data());
    double* grad_0_beta1 = d.deriv_base + 6 * d.BLK_SZ;
    deriv_x(grad_0_beta1, d.beta1, d.hx, d.sz, d.bflag);
});
grad_1_beta1_g.data(&tfdendro::data).work([grad_1_beta1_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_1_beta1_g.data());
    double* grad_1_beta1 = d.deriv_base + 7 * d.BLK_SZ;
    deriv_y(grad_1_beta1, d.beta1, d.hy, d.sz, d.bflag);
});
grad_2_beta1_g.data(&tfdendro::data).work([grad_2_beta1_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_2_beta1_g.data());
    double* grad_2_beta1 = d.deriv_base + 8 * d.BLK_SZ;
    deriv_z(grad_2_beta1, d.beta1, d.hz, d.sz, d.bflag);
});
grad_0_beta2_g.data(&tfdendro::data).work([grad_0_beta2_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_0_beta2_g.data());
    double* grad_0_beta2 = d.deriv_base + 9 * d.BLK_SZ;
    deriv_x(grad_0_beta2, d.beta2, d.hx, d.sz, d.bflag);
});
grad_1_beta2_g.data(&tfdendro::data).work([grad_1_beta2_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_1_beta2_g.data());
    double* grad_1_beta2 = d.deriv_base + 10 * d.BLK_SZ;
    deriv_y(grad_1_beta2, d.beta2, d.hy, d.sz, d.bflag);
});
grad_2_beta2_g.data(&tfdendro::data).work([grad_2_beta2_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_2_beta2_g.data());
    double* grad_2_beta2 = d.deriv_base + 11 * d.BLK_SZ;
    deriv_z(grad_2_beta2, d.beta2, d.hz, d.sz, d.bflag);
});
grad_0_B0_g.data(&tfdendro::data).work([grad_0_B0_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_0_B0_g.data());
    double* grad_0_B0 = d.deriv_base + 12 * d.BLK_SZ;
    deriv_x(grad_0_B0, d.B0, d.hx, d.sz, d.bflag);
});
grad_1_B0_g.data(&tfdendro::data).work([grad_1_B0_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_1_B0_g.data());
    double* grad_1_B0 = d.deriv_base + 13 * d.BLK_SZ;
    deriv_y(grad_1_B0, d.B0, d.hy, d.sz, d.bflag);
});
grad_2_B0_g.data(&tfdendro::data).work([grad_2_B0_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_2_B0_g.data());
    double* grad_2_B0 = d.deriv_base + 14 * d.BLK_SZ;
    deriv_z(grad_2_B0, d.B0, d.hz, d.sz, d.bflag);
});
grad_0_B1_g.data(&tfdendro::data).work([grad_0_B1_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_0_B1_g.data());
    double* grad_0_B1 = d.deriv_base + 15 * d.BLK_SZ;
    deriv_x(grad_0_B1, d.B1, d.hx, d.sz, d.bflag);
});
grad_1_B1_g.data(&tfdendro::data).work([grad_1_B1_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_1_B1_g.data());
    double* grad_1_B1 = d.deriv_base + 16 * d.BLK_SZ;
    deriv_y(grad_1_B1, d.B1, d.hy, d.sz, d.bflag);
});
grad_2_B1_g.data(&tfdendro::data).work([grad_2_B1_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_2_B1_g.data());
    double* grad_2_B1 = d.deriv_base + 17 * d.BLK_SZ;
    deriv_z(grad_2_B1, d.B1, d.hz, d.sz, d.bflag);
});
grad_0_B2_g.data(&tfdendro::data).work([grad_0_B2_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_0_B2_g.data());
    double* grad_0_B2 = d.deriv_base + 18 * d.BLK_SZ;
    deriv_x(grad_0_B2, d.B2, d.hx, d.sz, d.bflag);
});
grad_1_B2_g.data(&tfdendro::data).work([grad_1_B2_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_1_B2_g.data());
    double* grad_1_B2 = d.deriv_base + 19 * d.BLK_SZ;
    deriv_y(grad_1_B2, d.B2, d.hy, d.sz, d.bflag);
});
grad_2_B2_g.data(&tfdendro::data).work([grad_2_B2_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_2_B2_g.data());
    double* grad_2_B2 = d.deriv_base + 20 * d.BLK_SZ;
    deriv_z(grad_2_B2, d.B2, d.hz, d.sz, d.bflag);
});
grad_0_chi_g.data(&tfdendro::data).work([grad_0_chi_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_0_chi_g.data());
    double* grad_0_chi = d.deriv_base + 21 * d.BLK_SZ;
    deriv_x(grad_0_chi, d.chi, d.hx, d.sz, d.bflag);
});
grad_1_chi_g.data(&tfdendro::data).work([grad_1_chi_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_1_chi_g.data());
    double* grad_1_chi = d.deriv_base + 22 * d.BLK_SZ;
    deriv_y(grad_1_chi, d.chi, d.hy, d.sz, d.bflag);
});
grad_2_chi_g.data(&tfdendro::data).work([grad_2_chi_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_2_chi_g.data());
    double* grad_2_chi = d.deriv_base + 23 * d.BLK_SZ;
    deriv_z(grad_2_chi, d.chi, d.hz, d.sz, d.bflag);
});
grad_0_Gt0_g.data(&tfdendro::data).work([grad_0_Gt0_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_0_Gt0_g.data());
    double* grad_0_Gt0 = d.deriv_base + 24 * d.BLK_SZ;
    deriv_x(grad_0_Gt0, d.Gt0, d.hx, d.sz, d.bflag);
});
grad_1_Gt0_g.data(&tfdendro::data).work([grad_1_Gt0_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_1_Gt0_g.data());
    double* grad_1_Gt0 = d.deriv_base + 25 * d.BLK_SZ;
    deriv_y(grad_1_Gt0, d.Gt0, d.hy, d.sz, d.bflag);
});
grad_2_Gt0_g.data(&tfdendro::data).work([grad_2_Gt0_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_2_Gt0_g.data());
    double* grad_2_Gt0 = d.deriv_base + 26 * d.BLK_SZ;
    deriv_z(grad_2_Gt0, d.Gt0, d.hz, d.sz, d.bflag);
});
grad_0_Gt1_g.data(&tfdendro::data).work([grad_0_Gt1_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_0_Gt1_g.data());
    double* grad_0_Gt1 = d.deriv_base + 27 * d.BLK_SZ;
    deriv_x(grad_0_Gt1, d.Gt1, d.hx, d.sz, d.bflag);
});
grad_1_Gt1_g.data(&tfdendro::data).work([grad_1_Gt1_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_1_Gt1_g.data());
    double* grad_1_Gt1 = d.deriv_base + 28 * d.BLK_SZ;
    deriv_y(grad_1_Gt1, d.Gt1, d.hy, d.sz, d.bflag);
});
grad_2_Gt1_g.data(&tfdendro::data).work([grad_2_Gt1_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_2_Gt1_g.data());
    double* grad_2_Gt1 = d.deriv_base + 29 * d.BLK_SZ;
    deriv_z(grad_2_Gt1, d.Gt1, d.hz, d.sz, d.bflag);
});
grad_0_Gt2_g.data(&tfdendro::data).work([grad_0_Gt2_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_0_Gt2_g.data());
    double* grad_0_Gt2 = d.deriv_base + 30 * d.BLK_SZ;
    deriv_x(grad_0_Gt2, d.Gt2, d.hx, d.sz, d.bflag);
});
grad_1_Gt2_g.data(&tfdendro::data).work([grad_1_Gt2_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_1_Gt2_g.data());
    double* grad_1_Gt2 = d.deriv_base + 31 * d.BLK_SZ;
    deriv_y(grad_1_Gt2, d.Gt2, d.hy, d.sz, d.bflag);
});
grad_2_Gt2_g.data(&tfdendro::data).work([grad_2_Gt2_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_2_Gt2_g.data());
    double* grad_2_Gt2 = d.deriv_base + 32 * d.BLK_SZ;
    deriv_z(grad_2_Gt2, d.Gt2, d.hz, d.sz, d.bflag);
});
grad_0_K_g.data(&tfdendro::data).work([grad_0_K_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_0_K_g.data());
    double* grad_0_K = d.deriv_base + 33 * d.BLK_SZ;
    deriv_x(grad_0_K, d.K, d.hx, d.sz, d.bflag);
});
grad_1_K_g.data(&tfdendro::data).work([grad_1_K_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_1_K_g.data());
    double* grad_1_K = d.deriv_base + 34 * d.BLK_SZ;
    deriv_y(grad_1_K, d.K, d.hy, d.sz, d.bflag);
});
grad_2_K_g.data(&tfdendro::data).work([grad_2_K_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_2_K_g.data());
    double* grad_2_K = d.deriv_base + 35 * d.BLK_SZ;
    deriv_z(grad_2_K, d.K, d.hz, d.sz, d.bflag);
});
grad_0_gt0_g.data(&tfdendro::data).work([grad_0_gt0_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_0_gt0_g.data());
    double* grad_0_gt0 = d.deriv_base + 36 * d.BLK_SZ;
    deriv_x(grad_0_gt0, d.gt0, d.hx, d.sz, d.bflag);
});
grad_1_gt0_g.data(&tfdendro::data).work([grad_1_gt0_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_1_gt0_g.data());
    double* grad_1_gt0 = d.deriv_base + 37 * d.BLK_SZ;
    deriv_y(grad_1_gt0, d.gt0, d.hy, d.sz, d.bflag);
});
grad_2_gt0_g.data(&tfdendro::data).work([grad_2_gt0_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_2_gt0_g.data());
    double* grad_2_gt0 = d.deriv_base + 38 * d.BLK_SZ;
    deriv_z(grad_2_gt0, d.gt0, d.hz, d.sz, d.bflag);
});
grad_0_gt1_g.data(&tfdendro::data).work([grad_0_gt1_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_0_gt1_g.data());
    double* grad_0_gt1 = d.deriv_base + 39 * d.BLK_SZ;
    deriv_x(grad_0_gt1, d.gt1, d.hx, d.sz, d.bflag);
});
grad_1_gt1_g.data(&tfdendro::data).work([grad_1_gt1_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_1_gt1_g.data());
    double* grad_1_gt1 = d.deriv_base + 40 * d.BLK_SZ;
    deriv_y(grad_1_gt1, d.gt1, d.hy, d.sz, d.bflag);
});
grad_2_gt1_g.data(&tfdendro::data).work([grad_2_gt1_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_2_gt1_g.data());
    double* grad_2_gt1 = d.deriv_base + 41 * d.BLK_SZ;
    deriv_z(grad_2_gt1, d.gt1, d.hz, d.sz, d.bflag);
});
grad_0_gt2_g.data(&tfdendro::data).work([grad_0_gt2_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_0_gt2_g.data());
    double* grad_0_gt2 = d.deriv_base + 42 * d.BLK_SZ;
    deriv_x(grad_0_gt2, d.gt2, d.hx, d.sz, d.bflag);
});
grad_1_gt2_g.data(&tfdendro::data).work([grad_1_gt2_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_1_gt2_g.data());
    double* grad_1_gt2 = d.deriv_base + 43 * d.BLK_SZ;
    deriv_y(grad_1_gt2, d.gt2, d.hy, d.sz, d.bflag);
});
grad_2_gt2_g.data(&tfdendro::data).work([grad_2_gt2_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_2_gt2_g.data());
    double* grad_2_gt2 = d.deriv_base + 44 * d.BLK_SZ;
    deriv_z(grad_2_gt2, d.gt2, d.hz, d.sz, d.bflag);
});
grad_0_gt3_g.data(&tfdendro::data).work([grad_0_gt3_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_0_gt3_g.data());
    double* grad_0_gt3 = d.deriv_base + 45 * d.BLK_SZ;
    deriv_x(grad_0_gt3, d.gt3, d.hx, d.sz, d.bflag);
});
grad_1_gt3_g.data(&tfdendro::data).work([grad_1_gt3_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_1_gt3_g.data());
    double* grad_1_gt3 = d.deriv_base + 46 * d.BLK_SZ;
    deriv_y(grad_1_gt3, d.gt3, d.hy, d.sz, d.bflag);
});
grad_2_gt3_g.data(&tfdendro::data).work([grad_2_gt3_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_2_gt3_g.data());
    double* grad_2_gt3 = d.deriv_base + 47 * d.BLK_SZ;
    deriv_z(grad_2_gt3, d.gt3, d.hz, d.sz, d.bflag);
});
grad_0_gt4_g.data(&tfdendro::data).work([grad_0_gt4_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_0_gt4_g.data());
    double* grad_0_gt4 = d.deriv_base + 48 * d.BLK_SZ;
    deriv_x(grad_0_gt4, d.gt4, d.hx, d.sz, d.bflag);
});
grad_1_gt4_g.data(&tfdendro::data).work([grad_1_gt4_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_1_gt4_g.data());
    double* grad_1_gt4 = d.deriv_base + 49 * d.BLK_SZ;
    deriv_y(grad_1_gt4, d.gt4, d.hy, d.sz, d.bflag);
});
grad_2_gt4_g.data(&tfdendro::data).work([grad_2_gt4_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_2_gt4_g.data());
    double* grad_2_gt4 = d.deriv_base + 50 * d.BLK_SZ;
    deriv_z(grad_2_gt4, d.gt4, d.hz, d.sz, d.bflag);
});
grad_0_gt5_g.data(&tfdendro::data).work([grad_0_gt5_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_0_gt5_g.data());
    double* grad_0_gt5 = d.deriv_base + 51 * d.BLK_SZ;
    deriv_x(grad_0_gt5, d.gt5, d.hx, d.sz, d.bflag);
});
grad_1_gt5_g.data(&tfdendro::data).work([grad_1_gt5_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_1_gt5_g.data());
    double* grad_1_gt5 = d.deriv_base + 52 * d.BLK_SZ;
    deriv_y(grad_1_gt5, d.gt5, d.hy, d.sz, d.bflag);
});
grad_2_gt5_g.data(&tfdendro::data).work([grad_2_gt5_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_2_gt5_g.data());
    double* grad_2_gt5 = d.deriv_base + 53 * d.BLK_SZ;
    deriv_z(grad_2_gt5, d.gt5, d.hz, d.sz, d.bflag);
});
grad_0_At0_g.data(&tfdendro::data).work([grad_0_At0_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_0_At0_g.data());
    double* grad_0_At0 = d.deriv_base + 54 * d.BLK_SZ;
    deriv_x(grad_0_At0, d.At0, d.hx, d.sz, d.bflag);
});
grad_1_At0_g.data(&tfdendro::data).work([grad_1_At0_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_1_At0_g.data());
    double* grad_1_At0 = d.deriv_base + 55 * d.BLK_SZ;
    deriv_y(grad_1_At0, d.At0, d.hy, d.sz, d.bflag);
});
grad_2_At0_g.data(&tfdendro::data).work([grad_2_At0_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_2_At0_g.data());
    double* grad_2_At0 = d.deriv_base + 56 * d.BLK_SZ;
    deriv_z(grad_2_At0, d.At0, d.hz, d.sz, d.bflag);
});
grad_0_At1_g.data(&tfdendro::data).work([grad_0_At1_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_0_At1_g.data());
    double* grad_0_At1 = d.deriv_base + 57 * d.BLK_SZ;
    deriv_x(grad_0_At1, d.At1, d.hx, d.sz, d.bflag);
});
grad_1_At1_g.data(&tfdendro::data).work([grad_1_At1_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_1_At1_g.data());
    double* grad_1_At1 = d.deriv_base + 58 * d.BLK_SZ;
    deriv_y(grad_1_At1, d.At1, d.hy, d.sz, d.bflag);
});
grad_2_At1_g.data(&tfdendro::data).work([grad_2_At1_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_2_At1_g.data());
    double* grad_2_At1 = d.deriv_base + 59 * d.BLK_SZ;
    deriv_z(grad_2_At1, d.At1, d.hz, d.sz, d.bflag);
});
grad_0_At2_g.data(&tfdendro::data).work([grad_0_At2_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_0_At2_g.data());
    double* grad_0_At2 = d.deriv_base + 60 * d.BLK_SZ;
    deriv_x(grad_0_At2, d.At2, d.hx, d.sz, d.bflag);
});
grad_1_At2_g.data(&tfdendro::data).work([grad_1_At2_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_1_At2_g.data());
    double* grad_1_At2 = d.deriv_base + 61 * d.BLK_SZ;
    deriv_y(grad_1_At2, d.At2, d.hy, d.sz, d.bflag);
});
grad_2_At2_g.data(&tfdendro::data).work([grad_2_At2_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_2_At2_g.data());
    double* grad_2_At2 = d.deriv_base + 62 * d.BLK_SZ;
    deriv_z(grad_2_At2, d.At2, d.hz, d.sz, d.bflag);
});
grad_0_At3_g.data(&tfdendro::data).work([grad_0_At3_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_0_At3_g.data());
    double* grad_0_At3 = d.deriv_base + 63 * d.BLK_SZ;
    deriv_x(grad_0_At3, d.At3, d.hx, d.sz, d.bflag);
});
grad_1_At3_g.data(&tfdendro::data).work([grad_1_At3_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_1_At3_g.data());
    double* grad_1_At3 = d.deriv_base + 64 * d.BLK_SZ;
    deriv_y(grad_1_At3, d.At3, d.hy, d.sz, d.bflag);
});
grad_2_At3_g.data(&tfdendro::data).work([grad_2_At3_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_2_At3_g.data());
    double* grad_2_At3 = d.deriv_base + 65 * d.BLK_SZ;
    deriv_z(grad_2_At3, d.At3, d.hz, d.sz, d.bflag);
});
grad_0_At4_g.data(&tfdendro::data).work([grad_0_At4_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_0_At4_g.data());
    double* grad_0_At4 = d.deriv_base + 66 * d.BLK_SZ;
    deriv_x(grad_0_At4, d.At4, d.hx, d.sz, d.bflag);
});
grad_1_At4_g.data(&tfdendro::data).work([grad_1_At4_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_1_At4_g.data());
    double* grad_1_At4 = d.deriv_base + 67 * d.BLK_SZ;
    deriv_y(grad_1_At4, d.At4, d.hy, d.sz, d.bflag);
});
grad_2_At4_g.data(&tfdendro::data).work([grad_2_At4_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_2_At4_g.data());
    double* grad_2_At4 = d.deriv_base + 68 * d.BLK_SZ;
    deriv_z(grad_2_At4, d.At4, d.hz, d.sz, d.bflag);
});
grad_0_At5_g.data(&tfdendro::data).work([grad_0_At5_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_0_At5_g.data());
    double* grad_0_At5 = d.deriv_base + 69 * d.BLK_SZ;
    deriv_x(grad_0_At5, d.At5, d.hx, d.sz, d.bflag);
});
grad_1_At5_g.data(&tfdendro::data).work([grad_1_At5_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_1_At5_g.data());
    double* grad_1_At5 = d.deriv_base + 70 * d.BLK_SZ;
    deriv_y(grad_1_At5, d.At5, d.hy, d.sz, d.bflag);
});
grad_2_At5_g.data(&tfdendro::data).work([grad_2_At5_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad_2_At5_g.data());
    double* grad_2_At5 = d.deriv_base + 71 * d.BLK_SZ;
    deriv_z(grad_2_At5, d.At5, d.hz, d.sz, d.bflag);
});
grad2_0_0_gt0_g.data(&tfdendro::data).work([grad2_0_0_gt0_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_0_gt0_g.data());
    double* grad2_0_0_gt0 = d.deriv_base + 72 * d.BLK_SZ;
    deriv_xx(grad2_0_0_gt0, d.gt0, d.hx, d.sz, d.bflag);
});
grad2_0_1_gt0_g.data(&tfdendro::data).work([grad2_0_1_gt0_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_1_gt0_g.data());
    double* grad2_0_1_gt0 = d.deriv_base + 73 * d.BLK_SZ;
    double* grad_0_gt0 = d.deriv_base + 36 * d.BLK_SZ;
    deriv_y(grad2_0_1_gt0, grad_0_gt0, d.hy, d.sz, d.bflag);
});
grad2_0_2_gt0_g.data(&tfdendro::data).work([grad2_0_2_gt0_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_2_gt0_g.data());
    double* grad2_0_2_gt0 = d.deriv_base + 74 * d.BLK_SZ;
    double* grad_0_gt0 = d.deriv_base + 36 * d.BLK_SZ;
    deriv_z(grad2_0_2_gt0, grad_0_gt0, d.hz, d.sz, d.bflag);
});
grad2_1_1_gt0_g.data(&tfdendro::data).work([grad2_1_1_gt0_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_1_1_gt0_g.data());
    double* grad2_1_1_gt0 = d.deriv_base + 75 * d.BLK_SZ;
    deriv_yy(grad2_1_1_gt0, d.gt0, d.hy, d.sz, d.bflag);
});
grad2_1_2_gt0_g.data(&tfdendro::data).work([grad2_1_2_gt0_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_1_2_gt0_g.data());
    double* grad2_1_2_gt0 = d.deriv_base + 76 * d.BLK_SZ;
    double* grad_1_gt0 = d.deriv_base + 37 * d.BLK_SZ;
    deriv_z(grad2_1_2_gt0, grad_1_gt0, d.hz, d.sz, d.bflag);
});
grad2_2_2_gt0_g.data(&tfdendro::data).work([grad2_2_2_gt0_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_2_2_gt0_g.data());
    double* grad2_2_2_gt0 = d.deriv_base + 77 * d.BLK_SZ;
    deriv_zz(grad2_2_2_gt0, d.gt0, d.hz, d.sz, d.bflag);
});
grad2_0_0_gt1_g.data(&tfdendro::data).work([grad2_0_0_gt1_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_0_gt1_g.data());
    double* grad2_0_0_gt1 = d.deriv_base + 78 * d.BLK_SZ;
    deriv_xx(grad2_0_0_gt1, d.gt1, d.hx, d.sz, d.bflag);
});
grad2_0_1_gt1_g.data(&tfdendro::data).work([grad2_0_1_gt1_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_1_gt1_g.data());
    double* grad2_0_1_gt1 = d.deriv_base + 79 * d.BLK_SZ;
    double* grad_0_gt1 = d.deriv_base + 39 * d.BLK_SZ;
    deriv_y(grad2_0_1_gt1, grad_0_gt1, d.hy, d.sz, d.bflag);
});
grad2_0_2_gt1_g.data(&tfdendro::data).work([grad2_0_2_gt1_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_2_gt1_g.data());
    double* grad2_0_2_gt1 = d.deriv_base + 80 * d.BLK_SZ;
    double* grad_0_gt1 = d.deriv_base + 39 * d.BLK_SZ;
    deriv_z(grad2_0_2_gt1, grad_0_gt1, d.hz, d.sz, d.bflag);
});
grad2_1_1_gt1_g.data(&tfdendro::data).work([grad2_1_1_gt1_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_1_1_gt1_g.data());
    double* grad2_1_1_gt1 = d.deriv_base + 81 * d.BLK_SZ;
    deriv_yy(grad2_1_1_gt1, d.gt1, d.hy, d.sz, d.bflag);
});
grad2_1_2_gt1_g.data(&tfdendro::data).work([grad2_1_2_gt1_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_1_2_gt1_g.data());
    double* grad2_1_2_gt1 = d.deriv_base + 82 * d.BLK_SZ;
    double* grad_1_gt1 = d.deriv_base + 40 * d.BLK_SZ;
    deriv_z(grad2_1_2_gt1, grad_1_gt1, d.hz, d.sz, d.bflag);
});
grad2_2_2_gt1_g.data(&tfdendro::data).work([grad2_2_2_gt1_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_2_2_gt1_g.data());
    double* grad2_2_2_gt1 = d.deriv_base + 83 * d.BLK_SZ;
    deriv_zz(grad2_2_2_gt1, d.gt1, d.hz, d.sz, d.bflag);
});
grad2_0_0_gt2_g.data(&tfdendro::data).work([grad2_0_0_gt2_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_0_gt2_g.data());
    double* grad2_0_0_gt2 = d.deriv_base + 84 * d.BLK_SZ;
    deriv_xx(grad2_0_0_gt2, d.gt2, d.hx, d.sz, d.bflag);
});
grad2_0_1_gt2_g.data(&tfdendro::data).work([grad2_0_1_gt2_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_1_gt2_g.data());
    double* grad2_0_1_gt2 = d.deriv_base + 85 * d.BLK_SZ;
    double* grad_0_gt2 = d.deriv_base + 42 * d.BLK_SZ;
    deriv_y(grad2_0_1_gt2, grad_0_gt2, d.hy, d.sz, d.bflag);
});
grad2_0_2_gt2_g.data(&tfdendro::data).work([grad2_0_2_gt2_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_2_gt2_g.data());
    double* grad2_0_2_gt2 = d.deriv_base + 86 * d.BLK_SZ;
    double* grad_0_gt2 = d.deriv_base + 42 * d.BLK_SZ;
    deriv_z(grad2_0_2_gt2, grad_0_gt2, d.hz, d.sz, d.bflag);
});
grad2_1_1_gt2_g.data(&tfdendro::data).work([grad2_1_1_gt2_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_1_1_gt2_g.data());
    double* grad2_1_1_gt2 = d.deriv_base + 87 * d.BLK_SZ;
    deriv_yy(grad2_1_1_gt2, d.gt2, d.hy, d.sz, d.bflag);
});
grad2_1_2_gt2_g.data(&tfdendro::data).work([grad2_1_2_gt2_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_1_2_gt2_g.data());
    double* grad2_1_2_gt2 = d.deriv_base + 88 * d.BLK_SZ;
    double* grad_1_gt2 = d.deriv_base + 43 * d.BLK_SZ;
    deriv_z(grad2_1_2_gt2, grad_1_gt2, d.hz, d.sz, d.bflag);
});
grad2_2_2_gt2_g.data(&tfdendro::data).work([grad2_2_2_gt2_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_2_2_gt2_g.data());
    double* grad2_2_2_gt2 = d.deriv_base + 89 * d.BLK_SZ;
    deriv_zz(grad2_2_2_gt2, d.gt2, d.hz, d.sz, d.bflag);
});
grad2_0_0_gt3_g.data(&tfdendro::data).work([grad2_0_0_gt3_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_0_gt3_g.data());
    double* grad2_0_0_gt3 = d.deriv_base + 90 * d.BLK_SZ;
    deriv_xx(grad2_0_0_gt3, d.gt3, d.hx, d.sz, d.bflag);
});
grad2_0_1_gt3_g.data(&tfdendro::data).work([grad2_0_1_gt3_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_1_gt3_g.data());
    double* grad2_0_1_gt3 = d.deriv_base + 91 * d.BLK_SZ;
    double* grad_0_gt3 = d.deriv_base + 45 * d.BLK_SZ;
    deriv_y(grad2_0_1_gt3, grad_0_gt3, d.hy, d.sz, d.bflag);
});
grad2_0_2_gt3_g.data(&tfdendro::data).work([grad2_0_2_gt3_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_2_gt3_g.data());
    double* grad2_0_2_gt3 = d.deriv_base + 92 * d.BLK_SZ;
    double* grad_0_gt3 = d.deriv_base + 45 * d.BLK_SZ;
    deriv_z(grad2_0_2_gt3, grad_0_gt3, d.hz, d.sz, d.bflag);
});
grad2_1_1_gt3_g.data(&tfdendro::data).work([grad2_1_1_gt3_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_1_1_gt3_g.data());
    double* grad2_1_1_gt3 = d.deriv_base + 93 * d.BLK_SZ;
    deriv_yy(grad2_1_1_gt3, d.gt3, d.hy, d.sz, d.bflag);
});
grad2_1_2_gt3_g.data(&tfdendro::data).work([grad2_1_2_gt3_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_1_2_gt3_g.data());
    double* grad2_1_2_gt3 = d.deriv_base + 94 * d.BLK_SZ;
    double* grad_1_gt3 = d.deriv_base + 46 * d.BLK_SZ;
    deriv_z(grad2_1_2_gt3, grad_1_gt3, d.hz, d.sz, d.bflag);
});
grad2_2_2_gt3_g.data(&tfdendro::data).work([grad2_2_2_gt3_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_2_2_gt3_g.data());
    double* grad2_2_2_gt3 = d.deriv_base + 95 * d.BLK_SZ;
    deriv_zz(grad2_2_2_gt3, d.gt3, d.hz, d.sz, d.bflag);
});
grad2_0_0_gt4_g.data(&tfdendro::data).work([grad2_0_0_gt4_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_0_gt4_g.data());
    double* grad2_0_0_gt4 = d.deriv_base + 96 * d.BLK_SZ;
    deriv_xx(grad2_0_0_gt4, d.gt4, d.hx, d.sz, d.bflag);
});
grad2_0_1_gt4_g.data(&tfdendro::data).work([grad2_0_1_gt4_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_1_gt4_g.data());
    double* grad2_0_1_gt4 = d.deriv_base + 97 * d.BLK_SZ;
    double* grad_0_gt4 = d.deriv_base + 48 * d.BLK_SZ;
    deriv_y(grad2_0_1_gt4, grad_0_gt4, d.hy, d.sz, d.bflag);
});
grad2_0_2_gt4_g.data(&tfdendro::data).work([grad2_0_2_gt4_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_2_gt4_g.data());
    double* grad2_0_2_gt4 = d.deriv_base + 98 * d.BLK_SZ;
    double* grad_0_gt4 = d.deriv_base + 48 * d.BLK_SZ;
    deriv_z(grad2_0_2_gt4, grad_0_gt4, d.hz, d.sz, d.bflag);
});
grad2_1_1_gt4_g.data(&tfdendro::data).work([grad2_1_1_gt4_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_1_1_gt4_g.data());
    double* grad2_1_1_gt4 = d.deriv_base + 99 * d.BLK_SZ;
    deriv_yy(grad2_1_1_gt4, d.gt4, d.hy, d.sz, d.bflag);
});
grad2_1_2_gt4_g.data(&tfdendro::data).work([grad2_1_2_gt4_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_1_2_gt4_g.data());
    double* grad2_1_2_gt4 = d.deriv_base + 100 * d.BLK_SZ;
    double* grad_1_gt4 = d.deriv_base + 49 * d.BLK_SZ;
    deriv_z(grad2_1_2_gt4, grad_1_gt4, d.hz, d.sz, d.bflag);
});
grad2_2_2_gt4_g.data(&tfdendro::data).work([grad2_2_2_gt4_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_2_2_gt4_g.data());
    double* grad2_2_2_gt4 = d.deriv_base + 101 * d.BLK_SZ;
    deriv_zz(grad2_2_2_gt4, d.gt4, d.hz, d.sz, d.bflag);
});
grad2_0_0_gt5_g.data(&tfdendro::data).work([grad2_0_0_gt5_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_0_gt5_g.data());
    double* grad2_0_0_gt5 = d.deriv_base + 102 * d.BLK_SZ;
    deriv_xx(grad2_0_0_gt5, d.gt5, d.hx, d.sz, d.bflag);
});
grad2_0_1_gt5_g.data(&tfdendro::data).work([grad2_0_1_gt5_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_1_gt5_g.data());
    double* grad2_0_1_gt5 = d.deriv_base + 103 * d.BLK_SZ;
    double* grad_0_gt5 = d.deriv_base + 51 * d.BLK_SZ;
    deriv_y(grad2_0_1_gt5, grad_0_gt5, d.hy, d.sz, d.bflag);
});
grad2_0_2_gt5_g.data(&tfdendro::data).work([grad2_0_2_gt5_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_2_gt5_g.data());
    double* grad2_0_2_gt5 = d.deriv_base + 104 * d.BLK_SZ;
    double* grad_0_gt5 = d.deriv_base + 51 * d.BLK_SZ;
    deriv_z(grad2_0_2_gt5, grad_0_gt5, d.hz, d.sz, d.bflag);
});
grad2_1_1_gt5_g.data(&tfdendro::data).work([grad2_1_1_gt5_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_1_1_gt5_g.data());
    double* grad2_1_1_gt5 = d.deriv_base + 105 * d.BLK_SZ;
    deriv_yy(grad2_1_1_gt5, d.gt5, d.hy, d.sz, d.bflag);
});
grad2_1_2_gt5_g.data(&tfdendro::data).work([grad2_1_2_gt5_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_1_2_gt5_g.data());
    double* grad2_1_2_gt5 = d.deriv_base + 106 * d.BLK_SZ;
    double* grad_1_gt5 = d.deriv_base + 52 * d.BLK_SZ;
    deriv_z(grad2_1_2_gt5, grad_1_gt5, d.hz, d.sz, d.bflag);
});
grad2_2_2_gt5_g.data(&tfdendro::data).work([grad2_2_2_gt5_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_2_2_gt5_g.data());
    double* grad2_2_2_gt5 = d.deriv_base + 107 * d.BLK_SZ;
    deriv_zz(grad2_2_2_gt5, d.gt5, d.hz, d.sz, d.bflag);
});
grad2_0_0_chi_g.data(&tfdendro::data).work([grad2_0_0_chi_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_0_chi_g.data());
    double* grad2_0_0_chi = d.deriv_base + 108 * d.BLK_SZ;
    deriv_xx(grad2_0_0_chi, d.chi, d.hx, d.sz, d.bflag);
});
grad2_0_1_chi_g.data(&tfdendro::data).work([grad2_0_1_chi_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_1_chi_g.data());
    double* grad2_0_1_chi = d.deriv_base + 109 * d.BLK_SZ;
    double* grad_0_chi = d.deriv_base + 21 * d.BLK_SZ;
    deriv_y(grad2_0_1_chi, grad_0_chi, d.hy, d.sz, d.bflag);
});
grad2_0_2_chi_g.data(&tfdendro::data).work([grad2_0_2_chi_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_2_chi_g.data());
    double* grad2_0_2_chi = d.deriv_base + 110 * d.BLK_SZ;
    double* grad_0_chi = d.deriv_base + 21 * d.BLK_SZ;
    deriv_z(grad2_0_2_chi, grad_0_chi, d.hz, d.sz, d.bflag);
});
grad2_1_1_chi_g.data(&tfdendro::data).work([grad2_1_1_chi_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_1_1_chi_g.data());
    double* grad2_1_1_chi = d.deriv_base + 111 * d.BLK_SZ;
    deriv_yy(grad2_1_1_chi, d.chi, d.hy, d.sz, d.bflag);
});
grad2_1_2_chi_g.data(&tfdendro::data).work([grad2_1_2_chi_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_1_2_chi_g.data());
    double* grad2_1_2_chi = d.deriv_base + 112 * d.BLK_SZ;
    double* grad_1_chi = d.deriv_base + 22 * d.BLK_SZ;
    deriv_z(grad2_1_2_chi, grad_1_chi, d.hz, d.sz, d.bflag);
});
grad2_2_2_chi_g.data(&tfdendro::data).work([grad2_2_2_chi_g]() {
    auto d = *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_2_2_chi_g.data());
    double* grad2_2_2_chi = d.deriv_base + 113 * d.BLK_SZ;
    deriv_zz(grad2_2_2_chi, d.chi, d.hz, d.sz, d.bflag);
});
grad2_0_0_alpha_g.data(&tfdendro::data).work([grad2_0_0_alpha_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_0_alpha_g.data());
    double* grad2_0_0_alpha = d.deriv_base + 114 * d.BLK_SZ;
    deriv_xx(grad2_0_0_alpha, d.alpha, d.hx, d.sz, d.bflag);
});
grad2_0_1_alpha_g.data(&tfdendro::data).work([grad2_0_1_alpha_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_1_alpha_g.data());
    double* grad2_0_1_alpha = d.deriv_base + 115 * d.BLK_SZ;
    double* grad_0_alpha = d.deriv_base + 0 * d.BLK_SZ;
    deriv_y(grad2_0_1_alpha, grad_0_alpha, d.hy, d.sz, d.bflag);
});
grad2_0_2_alpha_g.data(&tfdendro::data).work([grad2_0_2_alpha_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_2_alpha_g.data());
    double* grad2_0_2_alpha = d.deriv_base + 116 * d.BLK_SZ;
    double* grad_0_alpha = d.deriv_base + 0 * d.BLK_SZ;
    deriv_z(grad2_0_2_alpha, grad_0_alpha, d.hz, d.sz, d.bflag);
});
grad2_1_1_alpha_g.data(&tfdendro::data).work([grad2_1_1_alpha_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_1_1_alpha_g.data());
    double* grad2_1_1_alpha = d.deriv_base + 117 * d.BLK_SZ;
    deriv_yy(grad2_1_1_alpha, d.alpha, d.hy, d.sz, d.bflag);
});
grad2_1_2_alpha_g.data(&tfdendro::data).work([grad2_1_2_alpha_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_1_2_alpha_g.data());
    double* grad2_1_2_alpha = d.deriv_base + 118 * d.BLK_SZ;
    double* grad_1_alpha = d.deriv_base + 1 * d.BLK_SZ;
    deriv_z(grad2_1_2_alpha, grad_1_alpha, d.hz, d.sz, d.bflag);
});
grad2_2_2_alpha_g.data(&tfdendro::data).work([grad2_2_2_alpha_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_2_2_alpha_g.data());
    double* grad2_2_2_alpha = d.deriv_base + 119 * d.BLK_SZ;
    deriv_zz(grad2_2_2_alpha, d.alpha, d.hz, d.sz, d.bflag);
});
grad2_0_0_beta0_g.data(&tfdendro::data).work([grad2_0_0_beta0_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_0_beta0_g.data());
    double* grad2_0_0_beta0 = d.deriv_base + 120 * d.BLK_SZ;
    deriv_xx(grad2_0_0_beta0, d.beta0, d.hx, d.sz, d.bflag);
});
grad2_0_1_beta0_g.data(&tfdendro::data).work([grad2_0_1_beta0_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_1_beta0_g.data());
    double* grad2_0_1_beta0 = d.deriv_base + 121 * d.BLK_SZ;
    double* grad_0_beta0 = d.deriv_base + 3 * d.BLK_SZ;
    deriv_y(grad2_0_1_beta0, grad_0_beta0, d.hy, d.sz, d.bflag);
});
grad2_0_2_beta0_g.data(&tfdendro::data).work([grad2_0_2_beta0_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_2_beta0_g.data());
    double* grad2_0_2_beta0 = d.deriv_base + 122 * d.BLK_SZ;
    double* grad_0_beta0 = d.deriv_base + 3 * d.BLK_SZ;
    deriv_z(grad2_0_2_beta0, grad_0_beta0, d.hz, d.sz, d.bflag);
});
grad2_1_1_beta0_g.data(&tfdendro::data).work([grad2_1_1_beta0_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_1_1_beta0_g.data());
    double* grad2_1_1_beta0 = d.deriv_base + 123 * d.BLK_SZ;
    deriv_yy(grad2_1_1_beta0, d.beta0, d.hy, d.sz, d.bflag);
});
grad2_1_2_beta0_g.data(&tfdendro::data).work([grad2_1_2_beta0_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_1_2_beta0_g.data());
    double* grad2_1_2_beta0 = d.deriv_base + 124 * d.BLK_SZ;
    double* grad_1_beta0 = d.deriv_base + 4 * d.BLK_SZ;
    deriv_z(grad2_1_2_beta0, grad_1_beta0, d.hz, d.sz, d.bflag);
});
grad2_2_2_beta0_g.data(&tfdendro::data).work([grad2_2_2_beta0_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_2_2_beta0_g.data());
    double* grad2_2_2_beta0 = d.deriv_base + 125 * d.BLK_SZ;
    deriv_zz(grad2_2_2_beta0, d.beta0, d.hz, d.sz, d.bflag);
});
grad2_0_0_beta1_g.data(&tfdendro::data).work([grad2_0_0_beta1_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_0_beta1_g.data());
    double* grad2_0_0_beta1 = d.deriv_base + 126 * d.BLK_SZ;
    deriv_xx(grad2_0_0_beta1, d.beta1, d.hx, d.sz, d.bflag);
});
grad2_0_1_beta1_g.data(&tfdendro::data).work([grad2_0_1_beta1_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_1_beta1_g.data());
    double* grad2_0_1_beta1 = d.deriv_base + 127 * d.BLK_SZ;
    double* grad_0_beta1 = d.deriv_base + 6 * d.BLK_SZ;
    deriv_y(grad2_0_1_beta1, grad_0_beta1, d.hy, d.sz, d.bflag);
});
grad2_0_2_beta1_g.data(&tfdendro::data).work([grad2_0_2_beta1_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_2_beta1_g.data());
    double* grad2_0_2_beta1 = d.deriv_base + 128 * d.BLK_SZ;
    double* grad_0_beta1 = d.deriv_base + 6 * d.BLK_SZ;
    deriv_z(grad2_0_2_beta1, grad_0_beta1, d.hz, d.sz, d.bflag);
});
grad2_1_1_beta1_g.data(&tfdendro::data).work([grad2_1_1_beta1_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_1_1_beta1_g.data());
    double* grad2_1_1_beta1 = d.deriv_base + 129 * d.BLK_SZ;
    deriv_yy(grad2_1_1_beta1, d.beta1, d.hy, d.sz, d.bflag);
});
grad2_1_2_beta1_g.data(&tfdendro::data).work([grad2_1_2_beta1_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_1_2_beta1_g.data());
    double* grad2_1_2_beta1 = d.deriv_base + 130 * d.BLK_SZ;
    double* grad_1_beta1 = d.deriv_base + 7 * d.BLK_SZ;
    deriv_z(grad2_1_2_beta1, grad_1_beta1, d.hz, d.sz, d.bflag);
});
grad2_2_2_beta1_g.data(&tfdendro::data).work([grad2_2_2_beta1_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_2_2_beta1_g.data());
    double* grad2_2_2_beta1 = d.deriv_base + 131 * d.BLK_SZ;
    deriv_zz(grad2_2_2_beta1, d.beta1, d.hz, d.sz, d.bflag);
});
grad2_0_0_beta2_g.data(&tfdendro::data).work([grad2_0_0_beta2_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_0_beta2_g.data());
    double* grad2_0_0_beta2 = d.deriv_base + 132 * d.BLK_SZ;
    deriv_xx(grad2_0_0_beta2, d.beta2, d.hx, d.sz, d.bflag);
});
grad2_0_1_beta2_g.data(&tfdendro::data).work([grad2_0_1_beta2_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_1_beta2_g.data());
    double* grad2_0_1_beta2 = d.deriv_base + 133 * d.BLK_SZ;
    double* grad_0_beta2 = d.deriv_base + 9 * d.BLK_SZ;
    deriv_y(grad2_0_1_beta2, grad_0_beta2, d.hy, d.sz, d.bflag);
});
grad2_0_2_beta2_g.data(&tfdendro::data).work([grad2_0_2_beta2_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_0_2_beta2_g.data());
    double* grad2_0_2_beta2 = d.deriv_base + 134 * d.BLK_SZ;
    double* grad_0_beta2 = d.deriv_base + 9 * d.BLK_SZ;
    deriv_z(grad2_0_2_beta2, grad_0_beta2, d.hz, d.sz, d.bflag);
});
grad2_1_1_beta2_g.data(&tfdendro::data).work([grad2_1_1_beta2_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_1_1_beta2_g.data());
    double* grad2_1_1_beta2 = d.deriv_base + 135 * d.BLK_SZ;
    deriv_yy(grad2_1_1_beta2, d.beta2, d.hy, d.sz, d.bflag);
});
grad2_1_2_beta2_g.data(&tfdendro::data).work([grad2_1_2_beta2_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_1_2_beta2_g.data());
    double* grad2_1_2_beta2 = d.deriv_base + 136 * d.BLK_SZ;
    double* grad_1_beta2 = d.deriv_base + 10 * d.BLK_SZ;
    deriv_z(grad2_1_2_beta2, grad_1_beta2, d.hz, d.sz, d.bflag);
});
grad2_2_2_beta2_g.data(&tfdendro::data).work([grad2_2_2_beta2_g]() {
    auto d =
        *static_cast<tfdendro::dendrotf_rhs_data*>(grad2_2_2_beta2_g.data());
    double* grad2_2_2_beta2 = d.deriv_base + 137 * d.BLK_SZ;
    deriv_zz(grad2_2_2_beta2, d.beta2, d.hz, d.sz, d.bflag);
});

grad_0_gt0_g.precede(grad2_0_1_gt0_g);
grad_0_gt0_g.precede(grad2_0_2_gt0_g);
grad_1_gt0_g.precede(grad2_1_2_gt0_g);
grad_0_gt1_g.precede(grad2_0_1_gt1_g);
grad_0_gt1_g.precede(grad2_0_2_gt1_g);
grad_1_gt1_g.precede(grad2_1_2_gt1_g);
grad_0_gt2_g.precede(grad2_0_1_gt2_g);
grad_0_gt2_g.precede(grad2_0_2_gt2_g);
grad_1_gt2_g.precede(grad2_1_2_gt2_g);
grad_0_gt3_g.precede(grad2_0_1_gt3_g);
grad_0_gt3_g.precede(grad2_0_2_gt3_g);
grad_1_gt3_g.precede(grad2_1_2_gt3_g);
grad_0_gt4_g.precede(grad2_0_1_gt4_g);
grad_0_gt4_g.precede(grad2_0_2_gt4_g);
grad_1_gt4_g.precede(grad2_1_2_gt4_g);
grad_0_gt5_g.precede(grad2_0_1_gt5_g);
grad_0_gt5_g.precede(grad2_0_2_gt5_g);
grad_1_gt5_g.precede(grad2_1_2_gt5_g);
grad_0_chi_g.precede(grad2_0_1_chi_g);
grad_0_chi_g.precede(grad2_0_2_chi_g);
grad_1_chi_g.precede(grad2_1_2_chi_g);
grad_0_alpha_g.precede(grad2_0_1_alpha_g);
grad_0_alpha_g.precede(grad2_0_2_alpha_g);
grad_1_alpha_g.precede(grad2_1_2_alpha_g);
grad_0_beta0_g.precede(grad2_0_1_beta0_g);
grad_0_beta0_g.precede(grad2_0_2_beta0_g);
grad_1_beta0_g.precede(grad2_1_2_beta0_g);
grad_0_beta1_g.precede(grad2_0_1_beta1_g);
grad_0_beta1_g.precede(grad2_0_2_beta1_g);
grad_1_beta1_g.precede(grad2_1_2_beta1_g);
grad_0_beta2_g.precede(grad2_0_1_beta2_g);
grad_0_beta2_g.precede(grad2_0_2_beta2_g);
grad_1_beta2_g.precede(grad2_1_2_beta2_g);
