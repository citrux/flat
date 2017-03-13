Data one_particle_simulation(const Params & params,
        Data & data) {
    Particle p(seed);
    const float field_dimensionless_factor = e * v_f * params.dt / eV;
    const float E_xc = params.E_xc * field_dimensionless_factor;
    const float E_yc = params.E_yc * field_dimensionless_factor;
    const float E_x = params.E_x  * field_dimensionless_factor;
    const float E_y = params.E_y  * field_dimensionless_factor;
    const float H = params.H * v_f / c * field_dimensionless_factor;

    const float omega = params.omega * params.dt;
    const float phi = params.phi;
    const float photon_energy = hbar * params.omega;

    const int all_time = params.all_time / params.dt;
    const float deps = params.deps;
    const float T = params.T;

    const float wla_max = sqr(k * T * Dak) * eV / (2 * hbar *
        rho * sqr(v_s * v_f) * hbar * hbar);
    const float wlo_max = k * T * sqr(Dopt) * eV / (4 * hbar *
        optical_phonon_energy *
        rho * sqr(v_f));
    while () {}
}