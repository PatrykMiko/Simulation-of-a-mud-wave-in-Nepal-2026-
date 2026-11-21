#pragma once

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <utility>
#include <vector>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace sciddica::kernels {

// --- constants mirroring sciddica_ss3/constants.py -------------------------

inline constexpr double GRAVITY = 9.81;
inline constexpr double EPS = 1.0e-9;
inline constexpr double MIN_THICKNESS = 1.0e-6;
inline constexpr double SQRT3 = 1.7320508075688772935274463415058723669;

// Axial coordinates (q, r) are stored as array[r, q].
// The six neighbours are separated by 60 degrees.
inline constexpr int NEIGHBOUR_DQ[6] = {1, 0, -1, -1, 0, 1};
inline constexpr int NEIGHBOUR_DR[6] = {0, 1, 1, 0, -1, -1};

inline constexpr double DIRECTION_X[6] = {
    1.0, 0.5, -0.5, -1.0, -0.5, 0.5};
inline constexpr double DIRECTION_Y[6] = {
    0.0, SQRT3 / 2.0, SQRT3 / 2.0, 0.0, -SQRT3 / 2.0, -SQRT3 / 2.0};

inline std::size_t idx2(std::size_t r, std::size_t q, std::size_t W) {
    return r * W + q;
}
inline std::size_t idx3(std::size_t d, std::size_t r, std::size_t q,
                        std::size_t H, std::size_t W) {
    return d * H * W + r * W + q;
}

// ---------------------------------------------------------------------------
// geometry: shifted / neighbour_stack
// ---------------------------------------------------------------------------

template <typename T>
void shifted(const T* src, T* dst, std::size_t H, std::size_t W, int dq,
             int dr, T fill) {
    std::fill(dst, dst + H * W, fill);
    const long long r0 = std::max<long long>(0, -dr);
    const long long r1 = std::min<long long>((long long)H, (long long)H - dr);
    const long long q0 = std::max<long long>(0, -dq);
    const long long q1 = std::min<long long>((long long)W, (long long)W - dq);
    if (r0 < r1 && q0 < q1) {
        for (long long r = r0; r < r1; ++r) {
            const long long sr = r + dr;
            for (long long q = q0; q < q1; ++q) {
                dst[idx2((std::size_t)r, (std::size_t)q, W)] =
                    src[idx2((std::size_t)sr, (std::size_t)(q + dq), W)];
            }
        }
    }
}

template <typename T>
void neighbour_stack(const T* src, T* out, std::size_t H, std::size_t W,
                     T fill) {
    for (int d = 0; d < 6; ++d) {
        shifted(src, out + idx3(d, 0, 0, H, W), H, W, NEIGHBOUR_DQ[d],
                NEIGHBOUR_DR[d], fill);
    }
}

// ---------------------------------------------------------------------------
// flow: minimization_algorithm
// ---------------------------------------------------------------------------

inline constexpr std::size_t OMP_GATE = 64 * 64;

template <typename T>
void minimization_algorithm(const T* mass, const T* centre,
                            const T* neighbours, const unsigned char* admissible,
                            T* outflow, std::size_t H, std::size_t W,
                            T eps) {
    const std::size_t n = H * W;
    const T inf = std::numeric_limits<T>::infinity();
    for (std::size_t i = 0; i < n; ++i) {
        const T m = mass[i];
        const T c = centre[i];

        T lv[7];
        bool fin[7];
        lv[0] = c;
        fin[0] = true;
        T lv_orig[7];
        lv_orig[0] = c;
        bool act[6];
        bool any_active = false;
        for (int d = 0; d < 6; ++d) {
            act[d] = admissible[idx3(d, 0, 0, H, W) + i] && (m > eps);
            any_active = any_active || act[d];
            const T nb = neighbours[idx3(d, 0, 0, H, W) + i];
            lv[d + 1] = act[d] ? nb : inf;
            fin[d + 1] = act[d];
            lv_orig[d + 1] = lv[d + 1];
        }

        T out_d[6] = {T(0), T(0), T(0), T(0), T(0), T(0)};
        if (!any_active || !(m > T(0))) {
            for (int d = 0; d < 6; ++d)
                outflow[idx3(d, 0, 0, H, W) + i] = T(0);
            continue;
        }

        // Insertion sort of (lv, fin) by lv ascending.
        for (int a = 1; a < 7; ++a) {
            T lv_key = lv[a];
            bool fin_key = fin[a];
            int b = a - 1;
            while (b >= 0 && lv[b] > lv_key) {
                lv[b + 1] = lv[b];
                fin[b + 1] = fin[b];
                --b;
            }
            lv[b + 1] = lv_key;
            fin[b + 1] = fin_key;
        }

        T cumsum[7] = {T(0)};
        T running = T(0);
        for (int k = 0; k < 7; ++k) {
            if (fin[k])
                running = running + lv[k];
            cumsum[k] = running;
        }

        T averages[7];
        bool filled[7];
        int k_filled = 0;
        for (int k = 0; k < 7; ++k) {
            const T count = T(k + 1);
            averages[k] = (m + cumsum[k]) / count;
            filled[k] = fin[k] && (averages[k] > lv[k]);
            k_filled += filled[k] ? 1 : 0;
        }

        const int gather = std::max(k_filled - 1, 0);
        const T fill_level = averages[gather];

        T total = T(0);
        for (int d = 0; d < 6; ++d) {
            T v = fill_level - lv_orig[d + 1];
            if (v < T(0))
                v = T(0);
            if (!act[d])
                v = T(0);
            out_d[d] = v;
            total = total + v;
        }

        T scale = T(1);
        if (total > eps)
            scale = std::min(T(1), m / total);

        const std::size_t base = idx3(0, 0, 0, H, W) + i;
        for (int d = 0; d < 6; ++d)
            outflow[base + d * H * W] = std::max(out_d[d] * scale, T(0));
    }
}

// ---------------------------------------------------------------------------
// flow: momentum_elevation_correction
// ---------------------------------------------------------------------------

template <typename T>
void momentum_elevation_correction(const T* thickness, const T* momentum_x,
                                   const T* momentum_y, const T* kinetic_head,
                                   const T* barycentre_x,
                                   const T* barycentre_y, T* out,
                                   std::size_t H, std::size_t W,
                                   double spacing, double pt, double g,
                                   double alpha, double eps) {
    const std::size_t n = H * W;
    const double apothem = spacing / 2.0;
    const double r_max = 0.6 * apothem;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(n >= OMP_GATE)
#endif
    for (std::size_t i = 0; i < n; ++i) {
        const double th = (double)thickness[i];
        const double px = (double)momentum_x[i];
        const double py = (double)momentum_y[i];
        double kh = (double)kinetic_head[i];
        if (kh < 0.0)
            kh = 0.0;
        const double bx0 = (double)barycentre_x[i];
        const double by0 = (double)barycentre_y[i];

        double w[6] = {0, 0, 0, 0, 0, 0};

        const double speed = std::min(std::sqrt(2.0 * g * kh), spacing / pt);
        const double magnitude = th * speed;

        const bool moving =
            (magnitude > eps) && (th > eps) &&
            (std::hypot(px, py) > eps);

        const double theta = std::atan2(py, px);

        double r = std::hypot(bx0, by0);
        double scale = 1.0;
        if (r > eps)
            scale = std::min(1.0, r_max / r);
        const double bx = bx0 * scale;
        const double by = by0 * scale;

        double face_angle[6];
        double rel[6];
        for (int d = 0; d < 6; ++d) {
            const double fx = apothem * DIRECTION_X[d] - bx;
            const double fy = apothem * DIRECTION_Y[d] - by;
            face_angle[d] = std::atan2(fy, fx);
            rel[d] = std::fmod(face_angle[d] - face_angle[0], 2.0 * M_PI);
            if (rel[d] < 0.0)
                rel[d] += 2.0 * M_PI;
        }
        double t_rel = std::fmod(theta - face_angle[0], 2.0 * M_PI);
        if (t_rel < 0.0)
            t_rel += 2.0 * M_PI;

        int lower = 0;
        for (int d = 0; d < 6; ++d)
            if (t_rel >= rel[d])
                lower = d;
        if (lower > 5)
            lower = 5;
        const int upper = (lower + 1) % 6;

        const double rel_lower = rel[lower];
        double rel_upper = rel[upper];
        if (lower == 5)
            rel_upper = 2.0 * M_PI;

        const double denom = std::max(rel_upper - rel_lower, eps);
        double fraction = (t_rel - rel_lower) / denom;
        if (fraction < 0.0)
            fraction = 0.0;
        if (fraction > 1.0)
            fraction = 1.0;

        double shares[6] = {0, 0, 0, 0, 0, 0};
        if (moving) {
            for (int d = 0; d < 6; ++d) {
                double s = 0.0;
                if (d == lower)
                    s = s + (1.0 - fraction);
                if (d == upper)
                    s = s + fraction;
                shares[d] = magnitude * s;
            }
        }

        for (int d = 0; d < 6; ++d) {
            const int rolled = (d + 3) % 6;
            double val = alpha * (shares[rolled] - shares[d]);
            double cap = 2.0 * th;
            if (val < -cap)
                val = -cap;
            if (val > cap)
                val = cap;
            w[d] = val;
        }

        const std::size_t base = idx3(0, 0, 0, H, W) + i;
        for (int d = 0; d < 6; ++d)
            out[base + d * H * W] = (T)w[d];
    }
}

// ---------------------------------------------------------------------------
// shift: circle segment scalar (double precision, mirrors numba kernel)
// ---------------------------------------------------------------------------

inline void circle_segment_scalar(double z, double r, double eps,
                                  double& frac, double& ec, double& ic) {
    double zc = z < -r ? -r : (z > r ? r : z);
    double rad = r * r - zc * zc;
    if (rad < 0.0)
        rad = 0.0;
    double root = std::sqrt(rad);
    double cap = r * r * std::acos(zc / r) - zc * root;
    double total = M_PI * r * r;
    double internal = total - cap;
    if (internal < 0.0)
        internal = 0.0;
    double f = cap / total;
    if (f < 0.0)
        f = 0.0;
    if (f > 1.0)
        f = 1.0;
    double rad15 = rad * std::sqrt(rad);
    double fm = (2.0 / 3.0) * rad15;
    double e = (cap > eps) ? (fm / cap) : 0.0;
    double ii = (internal > eps) ? (fm / internal) : 0.0;
    ii = -ii;
    if (f >= 1.0 - 1e-12)
        e = 0.0;
    if (f <= 1e-12)
        ii = 0.0;
    frac = f;
    ec = e;
    ic = ii;
}

// ---------------------------------------------------------------------------
// shift kernel (inert defaults identical to the numba version)
// ---------------------------------------------------------------------------

template <typename T>
void shift_kernel(const T* kh_in, const T* src, const T* nb, const T* bx,
                  const T* by, const T* fc, const unsigned char* sub,
                  const unsigned char* act,
                  double spacing, double dt, double g, double g_sub, double wr,
                  const double ux[6], const double uy[6], double eps,
                  T* fractions, T* travelled, T* heads, T* ibx, T* iby,
                  T* ebx, T* eby, std::size_t n) {
    const double apothem = spacing / 2.0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(n >= OMP_GATE)
#endif
    for (std::size_t i = 0; i < n; ++i) {
        double kh = (double)kh_in[i];
        if (kh < 0.0)
            kh = 0.0;

        // inert defaults
        for (int d = 0; d < 6; ++d) {
            fractions[idx3(d, 0, 0, 1, n) + i] = T(0);
            travelled[idx3(d, 0, 0, 1, n) + i] = T(0);
            heads[idx3(d, 0, 0, 1, n) + i] = (T)kh;
            ibx[idx3(d, 0, 0, 1, n) + i] = bx[i];
            iby[idx3(d, 0, 0, 1, n) + i] = by[i];
            ebx[idx3(d, 0, 0, 1, n) + i] = T(0);
            eby[idx3(d, 0, 0, 1, n) + i] = T(0);
        }

        const bool active = act[i];
        if (!active)
            continue;

        const double speed0 = std::sqrt(2.0 * g * kh);
        const bool is_sub = sub[i];

        for (int d = 0; d < 6; ++d) {
            const double pb = bx[i] * ux[d] + by[i] * uy[d];
            const double dx = spacing * ux[d] - bx[i];
            const double dy = spacing * uy[d] - by[i];
            double hdist = std::sqrt(dx * dx + dy * dy);
            if (hdist < eps)
                hdist = eps;

            const double drop = (double)src[i] - (double)nb[idx3(d, 0, 0, 1, n) + i];
            const double theta = std::atan2(drop, hdist);
            const bool can_move = (drop + kh) >= (double)fc[i] * hdist;

            const double gacc = is_sub ? g_sub : g;
            const double acc =
                gacc * (std::sin(theta) - (double)fc[i] * std::cos(theta));

            double sp = 0.0, sh = 0.0, new_speed = 0.0;
            if (is_sub && wr > eps) {
                const double decay = std::exp(-wr * dt);
                sp = speed0 * decay + (acc / wr) * (1.0 - decay);
                sh = (speed0 * (1.0 - decay) / wr) +
                     (acc / wr) * (dt - (1.0 - decay) / wr);
                const double den = wr * speed0 - acc;
                double ratio = (den > eps) ? (-acc / den) : 1.0;
                if (ratio < eps)
                    ratio = eps;
                if (ratio > 1.0)
                    ratio = 1.0;
                const double ts = -std::log(ratio) / wr;
                if (acc < 0.0 && speed0 >= 0.0 && ts < dt) {
                    const double sd = std::exp(-wr * ts);
                    sh = (speed0 * (1.0 - sd) / wr) +
                         (acc / wr) * (ts - (1.0 - sd) / wr);
                    sp = 0.0;
                }
                new_speed = sp < 0.0 ? 0.0 : sp;
            } else {
                sp = speed0 + acc * dt;
                if (sp < 0.0)
                    sp = 0.0;
                sh = speed0 * dt + 0.5 * acc * dt * dt;
                if (acc < 0.0 && (speed0 + acc * dt) < 0.0) {
                    double ts = (acc < -eps) ? (-speed0 / acc) : 0.0;
                    if (ts < 0.0)
                        ts = 0.0;
                    if (ts > dt)
                        ts = dt;
                    sh = speed0 * ts + 0.5 * acc * ts * ts;
                    sp = 0.0;
                }
                new_speed = sp;
            }

            double raw = sh;
            if (raw < 0.0)
                raw = 0.0;
            if (!can_move)
                raw = 0.0;

            double remaining = spacing - pb;
            if (remaining < 0.0)
                remaining = 0.0;
            double shift = raw < remaining ? raw : remaining;
            if (!can_move)
                shift = 0.0;

            double zc = apothem - pb - shift;
            double frac, cc, ic;
            circle_segment_scalar(zc, apothem, eps, frac, cc, ic);
            if (!can_move) {
                frac = 0.0;
                cc = 0.0;
                ic = 0.0;
            }

            const double perp_x = bx[i] - pb * ux[d];
            const double perp_y = by[i] - pb * uy[d];
            const double inner_along = pb + shift + ic;
            const double outer_along = pb + shift + cc - spacing;

            fractions[idx3(d, 0, 0, 1, n) + i] = (T)frac;
            travelled[idx3(d, 0, 0, 1, n) + i] = (T)shift;
            if (can_move) {
                heads[idx3(d, 0, 0, 1, n) + i] =
                    (T)(new_speed * new_speed / (2.0 * g));
                ibx[idx3(d, 0, 0, 1, n) + i] =
                    (T)(perp_x + inner_along * ux[d]);
                iby[idx3(d, 0, 0, 1, n) + i] =
                    (T)(perp_y + inner_along * uy[d]);
            }
            ebx[idx3(d, 0, 0, 1, n) + i] = (T)(perp_x + outer_along * ux[d]);
            eby[idx3(d, 0, 0, 1, n) + i] = (T)(perp_y + outer_along * uy[d]);
        }
    }
}

// ---------------------------------------------------------------------------
// _medium_arrays helper (cellwise air/water selection)
// ---------------------------------------------------------------------------

template <typename T>
void medium_arrays(const unsigned char* submerged, std::size_t n,
                   double adh_air, double pf_air, double pmt_air,
                   double per_air, double ptd_air, double ped_air,
                   double shift_friction_air,
                   double adh_w, double pf_w, double pmt_w, double per_w,
                   double ptd_w, double ped_w, double shift_friction_w,
                   T* adh, T* pf, T* pmt, T* per, T* ptd, T* ped,
                   T* friction) {
    for (std::size_t i = 0; i < n; ++i) {
        const bool s = submerged[i];
        adh[i] = (T)(s ? adh_w : adh_air);
        pf[i] = (T)(s ? pf_w : pf_air);
        pmt[i] = (T)(s ? pmt_w : pmt_air);
        per[i] = (T)(s ? per_w : per_air);
        ptd[i] = (T)(s ? ptd_w : ptd_air);
        ped[i] = (T)(s ? ped_w : ped_air);
        friction[i] = (T)(s ? shift_friction_w : shift_friction_air);
    }
}

// ---------------------------------------------------------------------------
// _record_impact
// ---------------------------------------------------------------------------

template <typename T>
void record_impact(unsigned char* mask, const T* thickness, std::size_t n,
                   T threshold) {
    for (std::size_t i = 0; i < n; ++i) {
        if (thickness[i] > threshold)
            mask[i] = 1;
    }
}

template <typename T>
struct StepWorkspace {
    std::size_t H = 0, W = 0, n = 0;
    std::vector<unsigned char> sub, sub_nb, active_region, active_mask, act_all, allowed1, allowed2;
    std::vector<T> adh, pf, pmt, per, ptd, ped, friction;
    std::vector<unsigned char> plane_u8, planes_u8;
    std::vector<T> above_threshold, erosion, head_after_erosion_loss, dilution, new_head_mob, mass_ratio, velocity_scale, momentum_scale;
    std::vector<T> a_nb, th_nb, bed_surface, neighbour_bed_surface;
    std::vector<T> movable_mass, central_unmovable, requested1, th_after1;
    std::vector<T> stationary_buf, sum_x_buf, sum_y_buf, bx1, by1, head1, px1, py1;
    std::vector<T> phase1_inflow, plane_buf, surface_after1, neighbour_surface_after1;
    std::vector<T> residual_mass_after1, movable_mass2, central_unmovable2, requested2;
    std::vector<T> retained, retained_head, px_retained, py_retained;
    std::vector<T> incoming_mass, incoming_head_mass, incoming_px, incoming_py, incoming_bx_mass, incoming_by_mass;
    std::vector<T> mass_buf, head_buf, r1_kh1_buf, r2_kh2_buf, mass_k_arr;
    std::vector<T> new_thickness, new_head_vals, new_px, new_py;
    std::vector<T> bx2, by2, retained_barycentre_mass, new_bx, new_by;
    std::vector<T> correction, phase1_neighbour_heights, energy_proxy_res;
    struct Pbuf {
        std::vector<T> requested, external, received, kinetic_head, shifted_kinetic_head, fraction, travelled, internal;
        std::vector<T> internal_barycentre_x, internal_barycentre_y, external_barycentre_x, external_barycentre_y;
    };
    Pbuf ph1, ph2;
    void resize_if_needed(std::size_t h, std::size_t w) {
        if (H==h && W==w && n==h*w) return;
        H=h; W=w; n=h*w;
        std::size_t n6 = 6*n;
        sub.assign(n,0); sub_nb.assign(n6,0); active_region.assign(n,0); active_mask.assign(n,0); act_all.assign(n,1); allowed1.assign(n6,0); allowed2.assign(n6,0);
        adh.assign(n,T(0)); pf.assign(n,T(0)); pmt.assign(n,T(0)); per.assign(n,T(0)); ptd.assign(n,T(0)); ped.assign(n,T(0)); friction.assign(n,T(0));
        plane_u8.assign(n,0); planes_u8.assign(n6,0);
        above_threshold.assign(n,T(0)); erosion.assign(n,T(0)); head_after_erosion_loss.assign(n,T(0)); dilution.assign(n,T(0)); new_head_mob.assign(n,T(0)); mass_ratio.assign(n,T(0)); velocity_scale.assign(n,T(0)); momentum_scale.assign(n,T(0));
        a_nb.assign(n6,T(0)); th_nb.assign(n6,T(0)); bed_surface.assign(n,T(0)); neighbour_bed_surface.assign(n6,T(0));
        movable_mass.assign(n,T(0)); central_unmovable.assign(n,T(0)); requested1.assign(n6,T(0)); th_after1.assign(n,T(0));
        stationary_buf.assign(n,T(0)); sum_x_buf.assign(n,T(0)); sum_y_buf.assign(n,T(0)); bx1.assign(n,T(0)); by1.assign(n,T(0)); head1.assign(n,T(0)); px1.assign(n,T(0)); py1.assign(n,T(0));
        phase1_inflow.assign(n,T(0)); plane_buf.assign(n,T(0)); surface_after1.assign(n,T(0)); neighbour_surface_after1.assign(n6,T(0));
        residual_mass_after1.assign(n,T(0)); movable_mass2.assign(n,T(0)); central_unmovable2.assign(n,T(0)); requested2.assign(n6,T(0));
        retained.assign(n,T(0)); retained_head.assign(n,T(0)); px_retained.assign(n,T(0)); py_retained.assign(n,T(0));
        incoming_mass.assign(n,T(0)); incoming_head_mass.assign(n,T(0)); incoming_px.assign(n,T(0)); incoming_py.assign(n,T(0)); incoming_bx_mass.assign(n,T(0)); incoming_by_mass.assign(n,T(0));
        mass_buf.assign(n,T(0)); head_buf.assign(n,T(0)); r1_kh1_buf.assign(n,T(0)); r2_kh2_buf.assign(n,T(0)); mass_k_arr.assign(n,T(0));
        new_thickness.assign(n,T(0)); new_head_vals.assign(n,T(0)); new_px.assign(n,T(0)); new_py.assign(n,T(0));
        bx2.assign(n,T(0)); by2.assign(n,T(0)); retained_barycentre_mass.assign(n,T(0)); new_bx.assign(n,T(0)); new_by.assign(n,T(0));
        correction.assign(n6,T(0)); phase1_neighbour_heights.assign(n6,T(0)); energy_proxy_res.assign(n,T(0));
        auto resize_p = [&](Pbuf& p){ p.requested.assign(n6,T(0)); p.external.assign(n6,T(0)); p.received.assign(n6,T(0)); p.kinetic_head.assign(n6,T(0)); p.shifted_kinetic_head.assign(n6,T(0)); p.fraction.assign(n6,T(0)); p.travelled.assign(n6,T(0)); p.internal.assign(n6,T(0)); p.internal_barycentre_x.assign(n6,T(0)); p.internal_barycentre_y.assign(n6,T(0)); p.external_barycentre_x.assign(n6,T(0)); p.external_barycentre_y.assign(n6,T(0)); };
        resize_p(ph1); resize_p(ph2);
    }
};
template <typename T> inline StepWorkspace<T>& get_ws() { thread_local StepWorkspace<T> ws; return ws; }


// ---------------------------------------------------------------------------
// step_full: one complete synchronous cellular-automaton step.
// Mirrors SciddicaSS3._step_full exactly (windowing stays on the Python side).
// Returns the mass-lost increment [m]; Python updates time/steps/mass_lost.
// ---------------------------------------------------------------------------

template <typename T>
T step_full(T* altitude, T* thickness, T* soil_depth, T* kinetic_head,
            T* momentum_x, T* momentum_y, T* barycentre_x, T* barycentre_y,
            const unsigned char* inside,
            // parameters
            double pc, double pt, double rho, double alpha_scale, double pml,
            double pwr,
            double adh_air, double pf_air, double pmt_air, double per_air,
            double ptd_air, double ped_air, double shift_friction_air,
            double adh_w, double pf_w, double pmt_w, double per_w,
            double ptd_w, double ped_w, double shift_friction_w,
            double rho_water, double water_level, int has_water_level,
            std::size_t H, std::size_t W) {
    auto& ws = get_ws<T>();
    ws.resize_if_needed(H,W);
    const std::size_t n = H * W;
    const T eps = (T)EPS;
    const T min_th = (T)MIN_THICKNESS;
    const double g = GRAVITY;
    const double spacing = SQRT3 * pc;
    const double cell_area = 1.5 * SQRT3 * pc * pc;
    const double g_sub = g * std::max(rho - rho_water, 0.0) / rho;
    const double alpha = alpha_scale * 2.0 * pt / spacing;

    // work buffers
    auto& sub = ws.sub; auto& sub_nb = ws.sub_nb; auto& active_region = ws.active_region; auto& active_mask = ws.active_mask;
    auto& adh = ws.adh; auto& pf = ws.pf; auto& pmt = ws.pmt; auto& per = ws.per; auto& ptd = ws.ptd; auto& ped = ws.ped; auto& friction = ws.friction;
    auto submerged_mask = [&](const T* alt) {
        for (std::size_t i = 0; i < n; ++i)
            sub[i] = has_water_level && (double)alt[i] < water_level ? 1 : 0;
    };

    auto fill_medium = [&]() {
        medium_arrays(sub.data(), n,
                      adh_air, pf_air, pmt_air, per_air, ptd_air, ped_air,
                      shift_friction_air,
                      adh_w, pf_w, pmt_w, per_w, ptd_w, ped_w,
                      shift_friction_w,
                      adh.data(), pf.data(), pmt.data(), per.data(),
                      ptd.data(), ped.data(), friction.data());
    };

    // neighbour_stack of a uchar mask into ws.planes_u8 (6,n)
    auto neighbour_bool_planes = [&](const unsigned char* mask, unsigned char* planes_out) {
        for (int d = 0; d < 6; ++d) {
            std::fill(ws.plane_u8.begin(), ws.plane_u8.end(), (unsigned char)0);
            shifted<unsigned char>(mask, ws.plane_u8.data(), H, W,
                                       NEIGHBOUR_DQ[d], NEIGHBOUR_DR[d],
                                       (unsigned char)0);
            const std::size_t base = (std::size_t)d * n;
            for (std::size_t i = 0; i < n; ++i)
                planes_out[base + i] = ws.plane_u8[i];
        }
    };

    auto any_active_region = [&](const std::vector<unsigned char>& act_mask) {
        bool any = false;
        for (std::size_t i = 0; i < n; ++i) {
            if (act_mask[i]) { any = true; break; }
        }
        return any;
    };

    // ---- mobilize soil (first pass) ----
    submerged_mask(altitude);
    fill_medium();

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(n >= OMP_GATE)
#endif
    for (std::size_t i = 0; i < n; ++i) {
        T v = kinetic_head[i] - pmt[i];
        T above_threshold = v > T(0) ? v : T(0);
        T e = per[i] * above_threshold;
        if (e > soil_depth[i]) e = soil_depth[i];
        if (!(thickness[i] > eps)) e = T(0);
        if (e < T(0)) e = T(0);
        T v2 = kinetic_head[i] - ped[i] * above_threshold;
        T head_after_erosion_loss = v2 < T(0) ? T(0) : v2;
        T denom = thickness[i] + e;
        T dilution = (denom > eps) ? (thickness[i] / denom) : T(1);
        T new_head_v = (e > eps) ? (head_after_erosion_loss * dilution) : kinetic_head[i];
        T mass_ratio = (thickness[i] > eps) ? ((thickness[i] + e) / thickness[i]) : T(1);
        T ratio = (kinetic_head[i] > eps) ? (new_head_v / kinetic_head[i]) : T(1);
        T velocity_scale = (T)std::sqrt((double)ratio);
        T momentum_scale = mass_ratio * velocity_scale;

        altitude[i] = altitude[i] - e;
        T sd = soil_depth[i] - e;
        soil_depth[i] = sd < T(0) ? T(0) : sd;
        thickness[i] = thickness[i] + e;
        kinetic_head[i] = new_head_v;
        momentum_x[i] = momentum_x[i] * momentum_scale;
        momentum_y[i] = momentum_y[i] * momentum_scale;
    }

    // submerged recomputed after altitude change
    submerged_mask(altitude);
    fill_medium();

    // ---- active region ----
    std::fill(ws.act_all.begin(), ws.act_all.end(), 1u); auto& act_all = ws.act_all;
    const T act_thresh = min_th * T(10.0);
    for (std::size_t i = 0; i < n; ++i)
        active_mask[i] = (thickness[i] > act_thresh)
                       || (kinetic_head[i] > eps)
                       || (std::hypot(momentum_x[i], momentum_y[i]) > eps);

    const bool any_active = any_active_region(active_mask);
    if (any_active) {
        neighbour_bool_planes(active_mask.data(), ws.planes_u8.data());
        for (std::size_t i = 0; i < n; ++i) {
            bool any = active_mask[i];
            for (int d = 0; d < 6 && !any; ++d)
                any = ws.planes_u8[idx3(d, 0, 0, 1, n) + i] > T(0.5);
            active_region[i] = any;
        }
    } else {
        std::fill(active_region.begin(), active_region.end(), false);
    }

    // sub_nb[d] = submerged neighbour in direction d
    neighbour_bool_planes(sub.data(), ws.planes_u8.data());
    for (std::size_t i = 0; i < n; ++i)
        for (int d = 0; d < 6; ++d)
            sub_nb[idx3(d, 0, 0, 1, n) + i] = ws.planes_u8[idx3(d, 0, 0, 1, n) + i] > T(0.5) ? 1 : 0;

    // a_nb = neighbour_stack(a, inf); th_nb = neighbour_stack(th, 0)
    auto& a_nb = ws.a_nb; auto& th_nb = ws.th_nb;
    neighbour_stack<T>(altitude, a_nb.data(), H, W,
                       (T)std::numeric_limits<T>::infinity());
    neighbour_stack<T>(thickness, th_nb.data(), H, W, T(0));

    // bed_surface = a + th ; neighbour_bed_surface = a_nb + th_nb
    auto& bed_surface = ws.bed_surface; auto& neighbour_bed_surface = ws.neighbour_bed_surface;
    for (std::size_t i = 0; i < n; ++i)
        bed_surface[i] = altitude[i] + thickness[i];
    for (int d = 0; d < 6; ++d)
        for (std::size_t i = 0; i < n; ++i)
            neighbour_bed_surface[idx3(d, 0, 0, 1, n) + i] =
                a_nb[idx3(d, 0, 0, 1, n) + i] +
                th_nb[idx3(d, 0, 0, 1, n) + i];

    // =======================================================================
    // Phase helper: _make_phase (writes into ws.ph1/ws.ph2 to avoid per-step allocations)
    // =======================================================================
    using Pbuf = typename StepWorkspace<T>::Pbuf;
    auto make_phase = [&](const std::vector<T>& requested,
                          const T* kinetic_head_src,
                          const T* source_surface,
                          const T* neighbour_surface,
                          const unsigned char* neighbour_submerged,
                          const T* barycentre_x_src,
                          const T* barycentre_y_src,
                          Pbuf& ph) -> void {
        std::copy(requested.begin(), requested.end(), ph.requested.begin());
        const unsigned char* act_mask = act_all.data();

        // shift_fraction
        shift_kernel(kinetic_head_src, source_surface, neighbour_surface,
                        barycentre_x_src, barycentre_y_src, friction.data(),
                        sub.data(), act_mask, spacing, pt, g, g_sub, pwr,
                        DIRECTION_X, DIRECTION_Y, (double)eps,
                        ph.fraction.data(), ph.travelled.data(),
                        ph.shifted_kinetic_head.data(),
                        ph.internal_barycentre_x.data(),
                        ph.internal_barycentre_y.data(),
                        ph.external_barycentre_x.data(),
                        ph.external_barycentre_y.data(), n);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(n >= OMP_GATE)
#endif
        for (std::size_t i = 0; i < n; ++i)
            for (int d = 0; d < 6; ++d) {
                const std::size_t base = idx3(d, 0, 0, 1, n) + i;
                ph.external[base] = requested[base] * ph.fraction[base];
                bool enters_water = (!sub[i]) && neighbour_submerged[base];
                T loss = enters_water ? (T)pml : T(0);
                ph.received[base] = ph.external[base] * (T(1) - loss);
                T eh = ph.shifted_kinetic_head[base] * (T(1) - ptd[i]);
                if (!(ph.external[base] > eps)) eh = T(0);
                ph.kinetic_head[base] = eh;
                T v = requested[base] - ph.external[base];
                ph.internal[base] = v < T(0) ? T(0) : v;
            }
    };

    // =======================================================================
    // Phase 1
    // =======================================================================
    // correction = momentum_elevation_correction(...)
    auto& correction = ws.correction;
    momentum_elevation_correction(thickness, momentum_x, momentum_y,
                                     kinetic_head, barycentre_x, barycentre_y,
                                     correction.data(), H, W, spacing, pt, g,
                                     alpha, (double)eps);

    auto& phase1_neighbour_heights = ws.phase1_neighbour_heights;
    for (int d = 0; d < 6; ++d)
        for (std::size_t i = 0; i < n; ++i)
            phase1_neighbour_heights[idx3(d, 0, 0, 1, n) + i] =
                neighbour_bed_surface[idx3(d, 0, 0, 1, n) + i] +
                correction[idx3(d, 0, 0, 1, n) + i];

    // allowed1
    auto& allowed1 = ws.allowed1;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(n >= OMP_GATE)
#endif
    for (int di = 0; di < 6 * (int)n; ++di) {
        int d = di / (int)n; std::size_t i = (std::size_t)(di % (int)n);
        if (any_active && !active_region[i]) { allowed1[idx3(d, 0, 0, 1, n) + i] = 0; continue; }
        bool v = (inside[idx3(d, 0, 0, 1, n) + i] != 0) &&
                 ((altitude[i] + thickness[i] -
                   phase1_neighbour_heights[idx3(d, 0, 0, 1, n) + i]) >
                  pf[i]);
        allowed1[idx3(d, 0, 0, 1, n) + i] = v ? 1 : 0;
    }

    // movable_mass = max(th - adh, 0); centre = a + adh + kh
    auto& movable_mass = ws.movable_mass; auto& central_unmovable = ws.central_unmovable;
    for (std::size_t i = 0; i < n; ++i) {
        T mm = thickness[i] - adh[i];
        movable_mass[i] = mm < T(0) ? T(0) : mm;
        central_unmovable[i] = altitude[i] + adh[i] + kinetic_head[i];
    }

    auto& requested1 = ws.requested1;
    minimization_algorithm(movable_mass.data(), central_unmovable.data(),
                              phase1_neighbour_heights.data(),
                              allowed1.data(), requested1.data(), H, W, eps);

    auto& phase1 = ws.ph1;
    make_phase(requested1, kinetic_head, bed_surface.data(),
               neighbour_bed_surface.data(), sub_nb.data(),
               barycentre_x, barycentre_y, phase1);

    // th_after1 = max(th - phase1.external.sum(0), 0)
    auto& th_after1 = ws.th_after1;
    for (std::size_t i = 0; i < n; ++i) {
        T sum_ext = T(0);
        for (int d = 0; d < 6; ++d)
            sum_ext = sum_ext + phase1.external[idx3(d, 0, 0, 1, n) + i];
        T v = thickness[i] - sum_ext;
        th_after1[i] = v < T(0) ? T(0) : v;
    }

    // _internal_barycentre (writes into ws.stationary_buf/sum_x_buf/sum_y_buf)
    auto internal_barycentre = [&](const T* mass_before, const T* bx_src,
                                   const T* by_src,
                                   const std::vector<T>& requested,
                                   const std::vector<T>& external,
                                   const std::vector<T>& ibx,
                                   const std::vector<T>& iby,
                                   T* bx_out, T* by_out) {
        auto& stationary = ws.stationary_buf; auto& sum_x = ws.sum_x_buf; auto& sum_y = ws.sum_y_buf;
        for (std::size_t i = 0; i < n; ++i) {
            T sum_req = T(0);
            for (int d = 0; d < 6; ++d)
                sum_req = sum_req + requested[idx3(d, 0, 0, 1, n) + i];
            T v = mass_before[i] - sum_req;
            stationary[i] = v < T(0) ? T(0) : v;
            sum_x[i] = stationary[i] * bx_src[i];
            sum_y[i] = stationary[i] * by_src[i];
        }
        for (std::size_t i = 0; i < n; ++i)
            for (int d = 0; d < 6; ++d) {
                const std::size_t base = idx3(d, 0, 0, 1, n) + i;
                T iv = requested[base] - external[base];
                if (iv < T(0))
                    iv = T(0);
                sum_x[i] = sum_x[i] + iv * ibx[base];
                sum_y[i] = sum_y[i] + iv * iby[base];
            }
        for (std::size_t i = 0; i < n; ++i) {
            T sum_ext = T(0);
            for (int d = 0; d < 6; ++d)
                sum_ext = sum_ext + external[idx3(d, 0, 0, 1, n) + i];
            T remaining = mass_before[i] - sum_ext;
            if (remaining < T(0))
                remaining = T(0);
            bx_out[i] = (remaining > eps) ? (sum_x[i] / remaining) : bx_src[i];
            by_out[i] = (remaining > eps) ? (sum_y[i] / remaining) : by_src[i];
        }
    };

    auto& bx1 = ws.bx1; auto& by1 = ws.by1;
    internal_barycentre(thickness, barycentre_x, barycentre_y,
                        phase1.requested, phase1.external,
                        phase1.internal_barycentre_x,
                        phase1.internal_barycentre_y,
                        bx1.data(), by1.data());

    // head1 = _head_after_internal_flow(th, th_after1, kh, phase1, ptd, False)
    auto head_after_internal_flow = [&](const T* mass_before,
                                        const T* mass_after,
                                        const T* kinetic_head_src,
                                        const Pbuf& ph,
                                        const T* ptd_arr,
                                        bool apply_turbulence,
                                        T* result) {
        auto& stationary = ws.stationary_buf; auto& energy_proxy = ws.energy_proxy_res;
        for (std::size_t i = 0; i < n; ++i) {
            T sum_req = T(0);
            for (int d = 0; d < 6; ++d)
                sum_req = sum_req + ph.requested[idx3(d, 0, 0, 1, n) + i];
            T v = mass_before[i] - sum_req;
            stationary[i] = v < T(0) ? T(0) : v;
            T khv = kinetic_head_src[i] < T(0) ? T(0) : kinetic_head_src[i];
            energy_proxy[i] = stationary[i] * khv;
        }
        for (std::size_t i = 0; i < n; ++i)
            for (int d = 0; d < 6; ++d) {
                const std::size_t base = idx3(d, 0, 0, 1, n) + i;
                energy_proxy[i] =
                    energy_proxy[i] + ph.internal[base] * ph.shifted_kinetic_head[base];
            }
        for (std::size_t i = 0; i < n; ++i) {
            T res;
            if (mass_after[i] > eps)
                res = energy_proxy[i] / mass_after[i];
            else {
                T khv = kinetic_head_src[i] < T(0) ? T(0) : kinetic_head_src[i];
                res = khv;
            }
            if (apply_turbulence)
                res = res * (T(1) - ptd_arr[i]);
            result[i] = res < T(0) ? T(0) : res;
        }
    };

    auto& head1 = ws.head1;
    head_after_internal_flow(thickness, th_after1.data(), kinetic_head, phase1, ptd.data(), false, head1.data());

    // px1, py1 = _momentum_after_internal_flow(th, mx, my, phase1, True)
    auto momentum_after_internal_flow = [&](const T* mass_before, const T* mx,
                                            const T* my, const Pbuf& ph,
                                            bool keep_stationary,
                                            T* px_out, T* py_out) {
        auto& stationary = ws.stationary_buf;
        for (std::size_t i = 0; i < n; ++i) {
            T sum_req = T(0);
            for (int d = 0; d < 6; ++d)
                sum_req = sum_req + ph.requested[idx3(d, 0, 0, 1, n) + i];
            T v = mass_before[i] - sum_req;
            stationary[i] = v < T(0) ? T(0) : v;
        }
        for (std::size_t i = 0; i < n; ++i) {
            T sf = (mass_before[i] > eps) ? (stationary[i] / mass_before[i])
                                          : T(0);
            px_out[i] = keep_stationary ? (mx[i] * sf) : T(0);
            py_out[i] = keep_stationary ? (my[i] * sf) : T(0);
        }
        for (std::size_t i = 0; i < n; ++i)
            for (int d = 0; d < 6; ++d) {
                const std::size_t base = idx3(d, 0, 0, 1, n) + i;
                double skh = (double)ph.shifted_kinetic_head[base];
                if (skh < 0.0)
                    skh = 0.0;
                double speed = std::sqrt(2.0 * g * skh);
                double im =
                    (double)rho * cell_area * (double)ph.internal[base] * speed;
                if (!(ph.travelled[base] > eps))
                    im = 0.0;
                px_out[i] = px_out[i] + (T)(im * DIRECTION_X[d]);
                py_out[i] = py_out[i] + (T)(im * DIRECTION_Y[d]);
            }
    };

    auto& px1 = ws.px1; auto& py1 = ws.py1;
    momentum_after_internal_flow(thickness, momentum_x, momentum_y, phase1, true, px1.data(), py1.data());

    // =======================================================================
    // Phase 2
    // =======================================================================
    // phase1_inflow = sum_d shifted(phase1.received[d], -dq, -dr, 0)
    auto& phase1_inflow = ws.phase1_inflow; auto& plane_buf = ws.plane_buf; std::fill(phase1_inflow.begin(), phase1_inflow.end(), T(0)); std::fill(plane_buf.begin(), plane_buf.end(), T(0));
    for (int d = 0; d < 6; ++d) {
        const T* src = phase1.received.data() + idx3(d, 0, 0, 1, n);
        shifted(src, plane_buf.data(), H, W, -NEIGHBOUR_DQ[d],
                   -NEIGHBOUR_DR[d], T(0));
        for (std::size_t i = 0; i < n; ++i)
            phase1_inflow[i] = phase1_inflow[i] + plane_buf[i];
    }

    // surface_after1 = a + th_after1 + phase1_inflow
    auto& surface_after1 = ws.surface_after1; auto& neighbour_surface_after1 = ws.neighbour_surface_after1;
    for (std::size_t i = 0; i < n; ++i)
        surface_after1[i] = altitude[i] + th_after1[i] + phase1_inflow[i];
    neighbour_stack<T>(surface_after1.data(), neighbour_surface_after1.data(),
                       H, W, (T)std::numeric_limits<T>::infinity());

    // allowed2
    auto& allowed2 = ws.allowed2;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(n >= OMP_GATE)
#endif
    for (int di = 0; di < 6 * (int)n; ++di) {
        int d = di / (int)n; std::size_t i = (std::size_t)(di % (int)n);
        if (any_active && !active_region[i]) { allowed2[idx3(d, 0, 0, 1, n) + i] = 0; continue; }
        bool v = (inside[idx3(d, 0, 0, 1, n) + i] != 0) &&
                 ((surface_after1[i] -
                   neighbour_surface_after1[idx3(d, 0, 0, 1, n) + i]) >
                  pf[i]);
        allowed2[idx3(d, 0, 0, 1, n) + i] = v ? 1 : 0;
    }

    // residual_mass_after1 = max(th - phase1.requested.sum(0), 0)
    auto& residual_mass_after1 = ws.residual_mass_after1;
    for (std::size_t i = 0; i < n; ++i) {
        T sum_req = T(0);
        for (int d = 0; d < 6; ++d)
            sum_req = sum_req + phase1.requested[idx3(d, 0, 0, 1, n) + i];
        T v = thickness[i] - sum_req;
        residual_mass_after1[i] = v < T(0) ? T(0) : v;
    }

    auto& movable_mass2 = ws.movable_mass2; auto& central_unmovable2 = ws.central_unmovable2;
    for (std::size_t i = 0; i < n; ++i) {
        T mm = residual_mass_after1[i] - adh[i];
        movable_mass2[i] = mm < T(0) ? T(0) : mm;
        central_unmovable2[i] = altitude[i] + adh[i] + head1[i];
    }
    auto& requested2 = ws.requested2;
    minimization_algorithm(movable_mass2.data(), central_unmovable2.data(),
                              neighbour_surface_after1.data(), allowed2.data(),
                              requested2.data(), H, W, eps);

    auto& phase2 = ws.ph2;
    make_phase(requested2, head1.data(), surface_after1.data(),
               neighbour_surface_after1.data(), sub_nb.data(),
               bx1.data(), by1.data(), phase2);

    // retained = max(th_after1 - phase2.external.sum(0), 0)
    auto& retained = ws.retained;
    for (std::size_t i = 0; i < n; ++i) {
        T sum_ext = T(0);
        for (int d = 0; d < 6; ++d)
            sum_ext = sum_ext + phase2.external[idx3(d, 0, 0, 1, n) + i];
        T v = th_after1[i] - sum_ext;
        retained[i] = v < T(0) ? T(0) : v;
    }

    // retained_head = _head_after_internal_flow(th_after1, retained, head1,
    //                                           phase2, ptd, True)
    auto& retained_head = ws.retained_head;
    head_after_internal_flow(th_after1.data(), retained.data(), head1.data(), phase2, ptd.data(), true, retained_head.data());

    // px_retained, py_retained = _momentum_after_internal_flow(th_after1,
    // px1, py1, phase2, False)
    auto& px_retained = ws.px_retained; auto& py_retained = ws.py_retained;
    momentum_after_internal_flow(th_after1.data(), px1.data(), py1.data(), phase2, true, px_retained.data(), py_retained.data());

    // turbulence scale
    for (std::size_t i = 0; i < n; ++i) {
        double v = 1.0 - (double)ptd[i];
        if (v < 0.0)
            v = 0.0;
        px_retained[i] = px_retained[i] * (T)std::sqrt(v);
        py_retained[i] = py_retained[i] * (T)std::sqrt(v);
    }

    // ---- _inflows(phase1, phase2) ----
    auto& incoming_mass = ws.incoming_mass; auto& incoming_head_mass = ws.incoming_head_mass; auto& incoming_px = ws.incoming_px; auto& incoming_py = ws.incoming_py; auto& incoming_bx_mass = ws.incoming_bx_mass; auto& incoming_by_mass = ws.incoming_by_mass; std::fill(incoming_mass.begin(),incoming_mass.end(),0); std::fill(incoming_head_mass.begin(),incoming_head_mass.end(),0); std::fill(incoming_px.begin(),incoming_px.end(),0); std::fill(incoming_py.begin(),incoming_py.end(),0); std::fill(incoming_bx_mass.begin(),incoming_bx_mass.end(),0); std::fill(incoming_by_mass.begin(),incoming_by_mass.end(),0);
    auto& mass_buf = ws.mass_buf; auto& head_buf = ws.head_buf; auto& r1_kh1_buf = ws.r1_kh1_buf; auto& r2_kh2_buf = ws.r2_kh2_buf; auto& mass_k_arr = ws.mass_k_arr;

    // Recompute inflows cleanly (direction-by-direction, accumulating).
    std::fill(incoming_mass.begin(), incoming_mass.end(), T(0));
    std::fill(incoming_head_mass.begin(), incoming_head_mass.end(), T(0));
    std::fill(incoming_px.begin(), incoming_px.end(), T(0));
    std::fill(incoming_py.begin(), incoming_py.end(), T(0));
    std::fill(incoming_bx_mass.begin(), incoming_bx_mass.end(), T(0));
    std::fill(incoming_by_mass.begin(), incoming_by_mass.end(), T(0));
    for (int d = 0; d < 6; ++d) {
        const int ndq = -NEIGHBOUR_DQ[d];
        const int ndr = -NEIGHBOUR_DR[d];
        const T* r1 = phase1.received.data() + idx3(d, 0, 0, 1, n);
        const T* r2 = phase2.received.data() + idx3(d, 0, 0, 1, n);
        const T* kh1 = phase1.kinetic_head.data() + idx3(d, 0, 0, 1, n);
        const T* kh2 = phase2.kinetic_head.data() + idx3(d, 0, 0, 1, n);
        const T* ebx1 = phase1.external_barycentre_x.data() +
                        idx3(d, 0, 0, 1, n);
        const T* eby1 = phase1.external_barycentre_y.data() +
                        idx3(d, 0, 0, 1, n);
        const T* ebx2 = phase2.external_barycentre_x.data() +
                        idx3(d, 0, 0, 1, n);
        const T* eby2 = phase2.external_barycentre_y.data() +
                        idx3(d, 0, 0, 1, n);

        // mass1 = shifted(r1), mass2 = shifted(r2)
        shifted(r1, mass_buf.data(), H, W, ndq, ndr, T(0));
        shifted(r2, head_buf.data(), H, W, ndq, ndr, T(0));

        // Save mass_k for all cells before overwriting buffers
        for (std::size_t i = 0; i < n; ++i) {
            mass_k_arr[i] = mass_buf[i] + head_buf[i];
            incoming_mass[i] = incoming_mass[i] + mass_k_arr[i];
        }

        // r1*kh1 and r2*kh2 products (separate buffers)
        for (std::size_t j = 0; j < n; ++j) {
            r1_kh1_buf[j] = r1[j] * kh1[j];
            r2_kh2_buf[j] = r2[j] * kh2[j];
        }

        // Shift r*kh products
        shifted(r1_kh1_buf.data(), plane_buf.data(), H, W, ndq, ndr, T(0));
        shifted(r2_kh2_buf.data(), mass_buf.data(), H, W, ndq, ndr, T(0));
        for (std::size_t i = 0; i < n; ++i) {
            T head_mass_k = plane_buf[i] + mass_buf[i];
            incoming_head_mass[i] = incoming_head_mass[i] + head_mass_k;

            T mean_head =
                (mass_k_arr[i] > eps) ? (head_mass_k / mass_k_arr[i]) : T(0);
            double speed_k = std::sqrt(2.0 * g * (double)mean_head);
            double momentum_k =
                (double)rho * cell_area * (double)mass_k_arr[i] * speed_k;
            incoming_px[i] =
                incoming_px[i] + (T)(momentum_k * DIRECTION_X[d]);
            incoming_py[i] =
                incoming_py[i] + (T)(momentum_k * DIRECTION_Y[d]);
        }

        // source_bx_mass = r1*ebx1[d] + r2*ebx2[d] ; shifted
        for (std::size_t j = 0; j < n; ++j) {
            mass_buf[j] = r1[j] * ebx1[j] + r2[j] * ebx2[j];
            head_buf[j] = r1[j] * eby1[j] + r2[j] * eby2[j];
        }
        shifted(mass_buf.data(), plane_buf.data(), H, W, ndq, ndr, T(0));
        for (std::size_t i = 0; i < n; ++i)
            incoming_bx_mass[i] = incoming_bx_mass[i] + plane_buf[i];
        shifted(head_buf.data(), plane_buf.data(), H, W, ndq, ndr, T(0));
        for (std::size_t i = 0; i < n; ++i)
            incoming_by_mass[i] = incoming_by_mass[i] + plane_buf[i];
    }

    // new_thickness = retained + incoming_mass
    // new_head = divide(retained*retained_head + incoming_head_mass,
    //                   new_thickness, where>eps) ; max(.,0)
    // new_px = px_retained + incoming_px ; new_py = ...
    auto& new_thickness = ws.new_thickness; auto& new_head_vals = ws.new_head_vals; auto& new_px = ws.new_px; auto& new_py = ws.new_py;
    for (std::size_t i = 0; i < n; ++i)
        new_thickness[i] = retained[i] + incoming_mass[i];
    for (std::size_t i = 0; i < n; ++i) {
        T res;
        if (new_thickness[i] > eps)
            res = (retained[i] * retained_head[i] + incoming_head_mass[i]) /
                  new_thickness[i];
        else
            res = T(0);
        new_head_vals[i] = res < T(0) ? T(0) : res;
    }
    for (std::size_t i = 0; i < n; ++i) {
        new_px[i] = px_retained[i] + incoming_px[i];
        new_py[i] = py_retained[i] + incoming_py[i];
    }

    // bx2, by2, retained_barycentre_mass = _internal_barycentre(th_after1,
    // bx1, by1, phase2.requested, phase2.external, ...)
    auto& bx2 = ws.bx2; auto& by2 = ws.by2;
    internal_barycentre(th_after1.data(), bx1.data(), by1.data(),
                        phase2.requested, phase2.external,
                        phase2.internal_barycentre_x,
                        phase2.internal_barycentre_y,
                        bx2.data(), by2.data());

    // retained_barycentre_mass = remaining from internal_barycentre (mass_before -
    // external.sum). Recompute:
    auto& retained_barycentre_mass = ws.retained_barycentre_mass;
    for (std::size_t i = 0; i < n; ++i) {
        T sum_ext = T(0);
        for (int d = 0; d < 6; ++d)
            sum_ext = sum_ext + phase2.external[idx3(d, 0, 0, 1, n) + i];
        T v = th_after1[i] - sum_ext;
        retained_barycentre_mass[i] = v < T(0) ? T(0) : v;
    }

    // new_bx = divide(retained_barycentre_mass*bx2 + incoming_bx_mass,
    //                 new_thickness, where>eps)
    auto& new_bx = ws.new_bx; auto& new_by = ws.new_by;
    for (std::size_t i = 0; i < n; ++i) {
        T vx = retained_barycentre_mass[i] * bx2[i] + incoming_bx_mass[i];
        T vy = retained_barycentre_mass[i] * by2[i] + incoming_by_mass[i];
        new_bx[i] =
            (new_thickness[i] > eps) ? (vx / new_thickness[i]) : T(0);
        new_by[i] =
            (new_thickness[i] > eps) ? (vy / new_thickness[i]) : T(0);
    }

    // lost_depth = sum(phase1.external - phase1.received) +
    //              sum(phase2.external - phase2.received)
    T lost_depth = T(0);
    for (std::size_t i = 0; i < n; ++i)
        for (int d = 0; d < 6; ++d) {
            const std::size_t base = idx3(d, 0, 0, 1, n) + i;
            lost_depth = lost_depth +
                         (phase1.external[base] - phase1.received[base]) +
                         (phase2.external[base] - phase2.received[base]);
        }

    // dead = new_thickness < MIN_THICKNESS (reuse ws.plane_u8)
    auto& dead = ws.plane_u8;
    for (std::size_t i = 0; i < n; ++i)
        dead[i] = new_thickness[i] < min_th ? 1 : 0;

    for (std::size_t i = 0; i < n; ++i) {
        if (dead[i]) {
            thickness[i] = T(0);
            kinetic_head[i] = T(0);
            momentum_x[i] = T(0);
            momentum_y[i] = T(0);
            barycentre_x[i] = T(0);
            barycentre_y[i] = T(0);
        } else {
            T nv = new_thickness[i];
            if (nv < T(0))
                nv = T(0);
            thickness[i] = nv;
            kinetic_head[i] = new_head_vals[i];
            momentum_x[i] = new_px[i];
            momentum_y[i] = new_py[i];
            barycentre_x[i] = new_bx[i];
            barycentre_y[i] = new_by[i];
        }
    }

    return lost_depth * (T)cell_area;
}

} // namespace sciddica::kernels
