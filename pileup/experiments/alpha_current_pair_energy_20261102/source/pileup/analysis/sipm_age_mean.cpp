// Deterministic acceleration of alpha_pair_deconvolution.AgeRecovery.forward.
// Same discrete mean-field model; no fast-math and no new physical parameters.
#include <algorithm>
#include <cmath>
#include <vector>

extern "C" int gagg_age_forward(const double* input, int samples, double cells,
    const double* recovery, const double* pde, const double* trigger,
    int ages, double xt, double* output) {
    if (samples < 1 || ages < 2 || cells <= 0) return 1;
    std::vector<double> state(ages, 0.0), hazard(ages), fire(ages), survived(ages);
    state.back() = 1.0;
    for (int i = 0; i < samples; ++i) {
        if (!std::isfinite(input[i]) || input[i] < 0) return 2;
        for (int j = 0; j < ages; ++j) hazard[j] = input[i] / cells * pde[j];
        double cross = 0.0, charge = 0.0;
        for (int iteration = 0; iteration < 6; ++iteration) {
            charge = 0.0;
            for (int j = 0; j < ages; ++j) {
                fire[j] = -std::expm1(-hazard[j] - cross * trigger[j]);
                charge += state[j] * (fire[j] * recovery[j]);
            }
            cross = xt * charge;
        }
        output[i] = cells * charge;
        double fired = 0.0;
        for (int j = 0; j < ages; ++j) {
            const double released = state[j] * fire[j];
            fired += released;
            survived[j] = state[j] - released;
        }
        for (int j = 1; j < ages; ++j) state[j] = survived[j - 1];
        state.back() += survived.back();
        state.front() = fired;
    }
    return 0;
}
