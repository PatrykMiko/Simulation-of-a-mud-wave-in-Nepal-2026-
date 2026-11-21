#include "kernels.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <cstring>

namespace py = pybind11;

namespace {

template <typename T>
T* ptr(py::array_t<T> a) {
    return static_cast<T*>(a.request().ptr);
}

template <typename T>
void check_2d(py::array_t<T> a, std::size_t H, std::size_t W) {
    auto r = a.request();
    if (r.ndim != 2 || r.shape[0] != (py::ssize_t)H ||
        r.shape[1] != (py::ssize_t)W)
        throw std::runtime_error("array must have shape (H, W)");
}

template <typename T>
py::array_t<T> make6(std::size_t H, std::size_t W) {
    return py::array_t<T>({(py::ssize_t)6, (py::ssize_t)H, (py::ssize_t)W});
}

template <typename T>
py::array_t<T> make2(std::size_t H, std::size_t W) {
    return py::array_t<T>({(py::ssize_t)H, (py::ssize_t)W});
}

py::array_t<unsigned char> make_uchar2(std::size_t H, std::size_t W) {
    return py::array_t<unsigned char>({(py::ssize_t)H, (py::ssize_t)W});
}

} // namespace

// ---------------------------------------------------------------------------
// step_full: one full synchronous step. Overloaded for float32 / float64.
// Python-API validation/array access stays under the GIL; the compute kernel
// runs with the GIL released (no Python API/object access, only raw buffers).
// ---------------------------------------------------------------------------

template <typename T>
T step_full_impl(py::array_t<T> altitude, py::array_t<T> thickness,
                 py::array_t<T> soil_depth, py::array_t<T> kinetic_head,
                 py::array_t<T> momentum_x, py::array_t<T> momentum_y,
                 py::array_t<T> barycentre_x, py::array_t<T> barycentre_y,
                 py::array inside, double pc, double pt, double rho,
                 double alpha_scale, double pml, double pwr,
                 double adh_air, double pf_air, double pmt_air, double per_air,
                 double ptd_air, double ped_air, double shift_friction_air,
                 double adh_w, double pf_w, double pmt_w, double per_w,
                 double ptd_w, double ped_w, double shift_friction_w,
                 double rho_water, double water_level, int has_water_level) {
    auto ra = altitude.request();
    std::size_t H = (std::size_t)ra.shape[0];
    std::size_t W = (std::size_t)ra.shape[1];
    check_2d(thickness, H, W);
    check_2d(soil_depth, H, W);
    check_2d(kinetic_head, H, W);
    check_2d(momentum_x, H, W);
    check_2d(momentum_y, H, W);
    check_2d(barycentre_x, H, W);
    check_2d(barycentre_y, H, W);

    auto ri = inside.request();
    if (ri.ndim != 3 || ri.shape[0] != 6 || ri.shape[1] != (py::ssize_t)H ||
        ri.shape[2] != (py::ssize_t)W)
        throw std::runtime_error("inside must have shape (6, H, W)");

    // Resolve raw pointers / values while the GIL is held (uses pybind11 accessors).
    T* p_alt = ptr(altitude); T* p_th = ptr(thickness); T* p_soil = ptr(soil_depth);
    T* p_kh = ptr(kinetic_head); T* p_mx = ptr(momentum_x); T* p_my = ptr(momentum_y);
    T* p_bx = ptr(barycentre_x); T* p_by = ptr(barycentre_y);
    auto* p_in = static_cast<unsigned char*>(ri.ptr);

    T ret;
    {
        py::gil_scoped_release release;
        ret = sciddica::kernels::step_full<T>(
            p_alt, p_th, p_soil, p_kh, p_mx, p_my, p_bx, p_by, p_in,
            pc, pt, rho, alpha_scale, pml, pwr,
            adh_air, pf_air, pmt_air, per_air, ptd_air, ped_air,
            shift_friction_air, adh_w, pf_w, pmt_w, per_w, ptd_w, ped_w,
            shift_friction_w, rho_water, water_level, has_water_level, H, W);
    }
    return ret;
}

#define BIND_STEP_FULL(T)                                                   \
    m.def("step_full", &step_full_impl<T>,                                  \
          py::arg("altitude"), py::arg("thickness"),                        \
          py::arg("soil_depth"), py::arg("kinetic_head"),                   \
          py::arg("momentum_x"), py::arg("momentum_y"),                     \
          py::arg("barycentre_x"), py::arg("barycentre_y"),                 \
          py::arg("inside"), py::arg("pc"), py::arg("pt"),                  \
          py::arg("rho"), py::arg("alpha_scale"), py::arg("pml"),           \
          py::arg("pwr"), py::arg("adh_air"), py::arg("pf_air"),            \
          py::arg("pmt_air"), py::arg("per_air"), py::arg("ptd_air"),       \
          py::arg("ped_air"), py::arg("shift_friction_air"),                \
          py::arg("adh_w"), py::arg("pf_w"), py::arg("pmt_w"),              \
          py::arg("per_w"), py::arg("ptd_w"), py::arg("ped_w"),             \
          py::arg("shift_friction_w"), py::arg("rho_water"),                \
          py::arg("water_level"), py::arg("has_water_level"),               \
          "Advance the cellular automaton one synchronous step (full grid).")

// ---------------------------------------------------------------------------
// record_impact
// ---------------------------------------------------------------------------

template <typename T>
void record_impact_impl(py::array_t<T> thickness, py::array mask,
                        double impact_threshold) {
    auto r = thickness.request();
    std::size_t H = (std::size_t)r.shape[0];
    std::size_t W = (std::size_t)r.shape[1];
    std::size_t n = H * W;
    auto rm = mask.request();
    if (rm.ndim != 2 || rm.shape[0] != (py::ssize_t)H ||
        rm.shape[1] != (py::ssize_t)W)
        throw std::runtime_error("mask must have shape (H, W)");
    const T* p_th = ptr(thickness);
    auto* p_mask = static_cast<unsigned char*>(rm.ptr);
    {
        py::gil_scoped_release release;
        sciddica::kernels::record_impact(p_mask, p_th, n, (T)impact_threshold);
    }
}

#define BIND_RECORD_IMPACT(T)                                               \
    m.def("record_impact", &record_impact_impl<T>,                          \
          py::arg("thickness"), py::arg("mask"), py::arg("impact_threshold"),\
          "Extend the maximum impact mask with currently impacted cells.")

// ---------------------------------------------------------------------------
// shift_fraction
// ---------------------------------------------------------------------------

template <typename T>
py::tuple shift_fraction_impl(
    py::array_t<T> kinetic_head, py::array_t<T> source_surface,
    py::array_t<T> neighbour_surface, py::array_t<T> barycentre_x,
    py::array_t<T> barycentre_y, py::array_t<T> shift_friction,
    py::array submerged, py::array active, double spacing, double pt,
    double g, double g_sub, double wr, double eps) {
    auto rk = kinetic_head.request();
    std::size_t H = (std::size_t)rk.shape[0];
    std::size_t W = (std::size_t)rk.shape[1];
    std::size_t n = H * W;

    auto rs = submerged.request();
    if (rs.ndim != 2 || rs.shape[0] != (py::ssize_t)H ||
        rs.shape[1] != (py::ssize_t)W)
        throw std::runtime_error("submerged must have shape (H, W)");

    std::vector<unsigned char> act(n, 1);
    if (!active.is_none()) {
        auto ra = active.request();
        if (ra.ndim != 2 || ra.shape[0] != (py::ssize_t)H ||
            ra.shape[1] != (py::ssize_t)W)
            throw std::runtime_error("active must have shape (H, W)");
        const unsigned char* ap = static_cast<const unsigned char*>(ra.ptr);
        for (std::size_t i = 0; i < n; ++i)
            act[i] = ap[i] ? (unsigned char)1 : (unsigned char)0;
    }

    py::array_t<T> fractions = make6<T>(H, W);
    py::array_t<T> travelled = make6<T>(H, W);
    py::array_t<T> heads = make6<T>(H, W);
    py::array_t<T> ibx = make6<T>(H, W);
    py::array_t<T> iby = make6<T>(H, W);
    py::array_t<T> ebx = make6<T>(H, W);
    py::array_t<T> eby = make6<T>(H, W);

    sciddica::kernels::shift_kernel(ptr(kinetic_head), ptr(source_surface),
                          ptr(neighbour_surface), ptr(barycentre_x),
                          ptr(barycentre_y), ptr(shift_friction),
                          static_cast<const unsigned char*>(rs.ptr),
                          act.data(), spacing, pt, g, g_sub, wr,
                          sciddica::kernels::DIRECTION_X, sciddica::kernels::DIRECTION_Y, eps,
                          ptr(fractions), ptr(travelled), ptr(heads),
                          ptr(ibx), ptr(iby), ptr(ebx), ptr(eby), n);

    return py::make_tuple(fractions, travelled, heads, ibx, iby, ebx, eby);
}

#define BIND_SHIFT_FRACTION(T)                                              \
    m.def("shift_fraction", &shift_fraction_impl<T>,                        \
          py::arg("kinetic_head"), py::arg("source_surface"),               \
          py::arg("neighbour_surface"), py::arg("barycentre_x"),            \
          py::arg("barycentre_y"), py::arg("shift_friction"),               \
          py::arg("submerged"), py::arg("active"), py::arg("spacing"),       \
          py::arg("pt"), py::arg("g"), py::arg("g_sub"), py::arg("wr"),      \
          py::arg("eps"),                                                   \
          "Shift kernel: external fraction, travelled, heads, barycentres.")

// ---------------------------------------------------------------------------
// minimization_algorithm
// ---------------------------------------------------------------------------

template <typename T>
py::array_t<T> minimization_algorithm_impl(
    py::array_t<T> movable_mass, py::array_t<T> central_unmovable_height,
    py::array_t<T> neighbour_heights, py::array allowed, double eps_d) {
    auto rm = movable_mass.request();
    std::size_t H = (std::size_t)rm.shape[0];
    std::size_t W = (std::size_t)rm.shape[1];
    std::size_t n = H * W;
    T eps = (T)eps_d;

    auto rc = central_unmovable_height.request();
    auto rn = neighbour_heights.request();
    auto ra = allowed.request();
    if (rc.ndim != 2 || rc.shape[0] != (py::ssize_t)H ||
        rc.shape[1] != (py::ssize_t)W)
        throw std::runtime_error("central_unmovable_height must be (H, W)");
    if (rn.ndim != 3 || rn.shape[0] != 6 || rn.shape[1] != (py::ssize_t)H ||
        rn.shape[2] != (py::ssize_t)W)
        throw std::runtime_error("neighbour_heights must be (6, H, W)");
    if (ra.ndim != 3 || ra.shape[0] != 6 || ra.shape[1] != (py::ssize_t)H ||
        ra.shape[2] != (py::ssize_t)W)
        throw std::runtime_error("allowed must be (6, H, W)");

    py::array_t<T> outflow = make6<T>(H, W);
    sciddica::kernels::minimization_algorithm(ptr(movable_mass),
                                    ptr(central_unmovable_height),
                                    ptr(neighbour_heights),
                                    static_cast<const unsigned char*>(ra.ptr),
                                    ptr(outflow), H, W, eps);
    return outflow;
}

#define BIND_MINIMIZATION(T)                                                \
    m.def("minimization_algorithm", &minimization_algorithm_impl<T>,        \
          py::arg("movable_mass"), py::arg("central_unmovable_height"),     \
          py::arg("neighbour_heights"), py::arg("allowed"),                 \
          py::arg("eps"),                                                   \
          "Minimization Algorithm of the Differences.")

// ---------------------------------------------------------------------------
// momentum_elevation_correction
// ---------------------------------------------------------------------------

template <typename T>
py::array_t<T> momentum_elevation_correction_impl(
    py::array_t<T> thickness, py::array_t<T> momentum_x,
    py::array_t<T> momentum_y, py::array_t<T> kinetic_head,
    py::array_t<T> barycentre_x, py::array_t<T> barycentre_y,
    double spacing, double pt, double g, double alpha_scale, double eps_d) {
    auto r = thickness.request();
    std::size_t H = (std::size_t)r.shape[0];
    std::size_t W = (std::size_t)r.shape[1];
    double alpha = alpha_scale * 2.0 * pt / spacing;
    T eps = (T)eps_d;
    py::array_t<T> out = make6<T>(H, W);
    sciddica::kernels::momentum_elevation_correction(ptr(thickness), ptr(momentum_x),
                                           ptr(momentum_y), ptr(kinetic_head),
                                           ptr(barycentre_x), ptr(barycentre_y),
                                           ptr(out), H, W, spacing, pt, g,
                                           alpha, (double)eps);
    return out;
}

#define BIND_MOMENTUM(T)                                                    \
    m.def("momentum_elevation_correction", &momentum_elevation_correction_impl<T>,\
          py::arg("thickness"), py::arg("momentum_x"), py::arg("momentum_y"),\
          py::arg("kinetic_head"), py::arg("barycentre_x"),                 \
          py::arg("barycentre_y"), py::arg("spacing"), py::arg("pt"),       \
          py::arg("g"), py::arg("alpha_scale"), py::arg("eps"),             \
          "SS3 fictitious-slope elevation correction.")

// ---------------------------------------------------------------------------
// geometry: shifted / neighbour_stack
// ---------------------------------------------------------------------------

template <typename T>
py::array_t<T> shifted_impl(py::array_t<T> array, int dq, int dr,
                            T fill_value) {
    auto r = array.request();
    std::size_t H = (std::size_t)r.shape[0];
    std::size_t W = (std::size_t)r.shape[1];
    py::array_t<T> out = make2<T>(H, W);
    sciddica::kernels::shifted(ptr(array), ptr(out), H, W, dq, dr, fill_value);
    return out;
}

template <typename T>
py::array_t<T> neighbour_stack_impl(py::array_t<T> array, T fill_value) {
    auto r = array.request();
    std::size_t H = (std::size_t)r.shape[0];
    std::size_t W = (std::size_t)r.shape[1];
    py::array_t<T> out = make6<T>(H, W);
    sciddica::kernels::neighbour_stack(ptr(array), ptr(out), H, W, fill_value);
    return out;
}

#define BIND_SHIFTED(T)                                                     \
    m.def("shifted", &shifted_impl<T>, py::arg("array"), py::arg("dq"),      \
          py::arg("dr"), py::arg("fill_value"),                             \
          "Shift a 2-D array by (dq, dr) without wrapping (geometry).")

#define BIND_NEIGHBOUR_STACK(T)                                             \
    m.def("neighbour_stack", &neighbour_stack_impl<T>,                      \
          py::arg("array"), py::arg("fill_value"),                          \
          "Stack the six neighbours into a (6, H, W) array.")

// ---------------------------------------------------------------------------
// mobilize_soil
// ---------------------------------------------------------------------------

template <typename T>
py::tuple mobilize_soil_impl(py::array_t<T> altitude, py::array_t<T> thickness,
                             py::array_t<T> soil_depth,
                             py::array_t<T> kinetic_head,
                             py::array_t<T> momentum_x,
                             py::array_t<T> momentum_y,
                             double water_level, int has_water_level,
                             double eps_d,
                             double adh_air, double pf_air, double pmt_air,
                             double per_air, double ptd_air, double ped_air,
                             double shift_friction_air,
                             double adh_w, double pf_w, double pmt_w,
                             double per_w, double ptd_w, double ped_w,
                             double shift_friction_w) {
    auto r = altitude.request();
    std::size_t H = (std::size_t)r.shape[0];
    std::size_t W = (std::size_t)r.shape[1];
    std::size_t n = H * W;
    T eps = (T)eps_d;

    std::vector<unsigned char> sub(n);
    for (std::size_t i = 0; i < n; ++i)
        sub[i] = has_water_level && (double)ptr(altitude)[i] < water_level
                     ? 1
                     : 0;

    std::vector<T> adh(n), pf(n), pmt(n), per(n), ptd(n), ped(n),
        friction(n);
    sciddica::kernels::medium_arrays(sub.data(), n,
                           adh_air, pf_air, pmt_air, per_air, ptd_air, ped_air,
                           shift_friction_air,
                           adh_w, pf_w, pmt_w, per_w, ptd_w, ped_w,
                           shift_friction_w,
                           adh.data(), pf.data(), pmt.data(), per.data(),
                           ptd.data(), ped.data(), friction.data());

    std::vector<T> above_threshold(n), erosion(n), head_after_erosion(n),
        dilution(n), new_head(n), mass_ratio(n), velocity_scale(n),
        momentum_scale(n);
    for (std::size_t i = 0; i < n; ++i) {
        T v = ptr(kinetic_head)[i] - pmt[i];
        above_threshold[i] = v > T(0) ? v : T(0);
    }
    for (std::size_t i = 0; i < n; ++i) {
        T e = per[i] * above_threshold[i];
        if (e > ptr(soil_depth)[i])
            e = ptr(soil_depth)[i];
        if (!(ptr(thickness)[i] > eps))
            e = T(0);
        if (e < T(0))
            e = T(0);
        erosion[i] = e;
    }
    for (std::size_t i = 0; i < n; ++i) {
        T v = ptr(kinetic_head)[i] - ped[i] * above_threshold[i];
        head_after_erosion[i] = v < T(0) ? T(0) : v;
    }
    for (std::size_t i = 0; i < n; ++i) {
        T denom = ptr(thickness)[i] + erosion[i];
        dilution[i] = (denom > eps) ? (ptr(thickness)[i] / denom) : T(1);
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (erosion[i] > eps)
            new_head[i] = head_after_erosion[i] * dilution[i];
        else
            new_head[i] = ptr(kinetic_head)[i];
    }
    for (std::size_t i = 0; i < n; ++i) {
        mass_ratio[i] = (ptr(thickness)[i] > eps)
                            ? ((ptr(thickness)[i] + erosion[i]) /
                               ptr(thickness)[i])
                            : T(1);
        T ratio = (ptr(kinetic_head)[i] > eps)
                      ? (new_head[i] / ptr(kinetic_head)[i])
                      : T(1);
        velocity_scale[i] = (T)std::sqrt((double)ratio);
        momentum_scale[i] = mass_ratio[i] * velocity_scale[i];
    }
    for (std::size_t i = 0; i < n; ++i) {
        ptr(altitude)[i] = ptr(altitude)[i] - erosion[i];
        T sd = ptr(soil_depth)[i] - erosion[i];
        ptr(soil_depth)[i] = sd < T(0) ? T(0) : sd;
        ptr(thickness)[i] = ptr(thickness)[i] + erosion[i];
        ptr(kinetic_head)[i] = new_head[i];
        ptr(momentum_x)[i] = ptr(momentum_x)[i] * momentum_scale[i];
        ptr(momentum_y)[i] = ptr(momentum_y)[i] * momentum_scale[i];
    }

    for (std::size_t i = 0; i < n; ++i)
        sub[i] = has_water_level && (double)ptr(altitude)[i] < water_level
                     ? 1
                     : 0;
    sciddica::kernels::medium_arrays(sub.data(), n,
                           adh_air, pf_air, pmt_air, per_air, ptd_air, ped_air,
                           shift_friction_air,
                           adh_w, pf_w, pmt_w, per_w, ptd_w, ped_w,
                           shift_friction_w,
                           adh.data(), pf.data(), pmt.data(), per.data(),
                           ptd.data(), ped.data(), friction.data());

    py::array_t<unsigned char> out_sub = make_uchar2(H, W);
    std::memcpy(ptr(out_sub), sub.data(), n);
    py::array_t<T> out_adh = make2<T>(H, W);
    py::array_t<T> out_pf = make2<T>(H, W);
    py::array_t<T> out_ptd = make2<T>(H, W);
    py::array_t<T> out_friction = make2<T>(H, W);
    std::memcpy(ptr(out_adh), adh.data(), n * sizeof(T));
    std::memcpy(ptr(out_pf), pf.data(), n * sizeof(T));
    std::memcpy(ptr(out_ptd), ptd.data(), n * sizeof(T));
    std::memcpy(ptr(out_friction), friction.data(), n * sizeof(T));

    return py::make_tuple(out_sub, out_adh, out_pf, out_ptd, out_friction);
}

#define BIND_MOBILIZE_SOIL(T)                                               \
    m.def("mobilize_soil", &mobilize_soil_impl<T>,                          \
          py::arg("altitude"), py::arg("thickness"), py::arg("soil_depth"),  \
          py::arg("kinetic_head"), py::arg("momentum_x"),                   \
          py::arg("momentum_y"), py::arg("water_level"),                    \
          py::arg("has_water_level"), py::arg("eps"),                       \
          py::arg("adh_air"), py::arg("pf_air"), py::arg("pmt_air"),         \
          py::arg("per_air"), py::arg("ptd_air"), py::arg("ped_air"),        \
          py::arg("shift_friction_air"),                                    \
          py::arg("adh_w"), py::arg("pf_w"), py::arg("pmt_w"),               \
          py::arg("per_w"), py::arg("ptd_w"), py::arg("ped_w"),              \
          py::arg("shift_friction_w"),                                      \
          "Soil mobilization (mutates state; returns submerged, adh, pf,"   \
          " ptd, friction).")

PYBIND11_MODULE(sciddica_ss3_cpp, m) {
    m.doc() = "SCIDDICA-SS3 C++ kernels";
    m.attr("__version__") = "0.2.0";
    m.attr("ABI_VERSION") = 2;

    BIND_STEP_FULL(float);
    BIND_STEP_FULL(double);
    BIND_RECORD_IMPACT(float);
    BIND_RECORD_IMPACT(double);
    BIND_SHIFT_FRACTION(float);
    BIND_SHIFT_FRACTION(double);
    BIND_MINIMIZATION(float);
    BIND_MINIMIZATION(double);
    BIND_MOMENTUM(float);
    BIND_MOMENTUM(double);
    BIND_SHIFTED(float);
    BIND_SHIFTED(double);
    BIND_NEIGHBOUR_STACK(float);
    BIND_NEIGHBOUR_STACK(double);
    BIND_MOBILIZE_SOIL(float);
    BIND_MOBILIZE_SOIL(double);
}
