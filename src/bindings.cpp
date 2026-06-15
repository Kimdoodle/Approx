#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../include/parse.h"
#include "../include/error_bound.h"
#include "../include/fakeciphertext.h"
#include "../include/eval_plain.h"

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace py = pybind11;

class DecompCache {
public:
    nlohmann::json cache_root;

    DecompCache(const std::string& cache_file) {
        cache_root = load_decomp_cache(cache_file);
        if (cache_root.is_null()) {
            throw std::runtime_error("Failed to load decomp cache: " + cache_file);
        }
    }

    EvalStep restore_EvalStep(const std::vector<double>& coeff) const {
        const std::string key = get_poly_type_key(coeff);
        if (!cache_root.contains(key)) {
            throw std::runtime_error("Cache key not found: " + key);
        }

        Poly poly(coeff);
        std::shared_ptr<Decomp> dcmp =
            reconstruct_decomp_from_cache(cache_root.at(key), poly, 0);

        if (!dcmp) {
            throw std::runtime_error("Failed to reconstruct Decomp from cache.");
        }

        return EvalStep(dcmp);
    }
};

EvalStep restore_EvalStep(const std::vector<double>& coeff, const std::string& cache_file) {
    DecompCache cache(cache_file);
    return cache.restore_EvalStep(coeff);
}

fct::Ciphertext evaluate(fct::Ciphertext x, const EvalStep& es, ErrBound& eb) {
    EvalStep es_copy = es;
    return fct::evaluate(x, es_copy, eb);
}

std::pair<double, double> estimate_range(ErrBound& eb, std::pair<double, double> interval, const EvalStep& es, double interval_size, bool lipschitz) {
    EvalStep es_copy = es;
    return fct::estimate_range(eb, interval, es_copy, interval_size, lipschitz);
}

std::pair<double, double> estimate_composite_range(ErrBound& eb, std::pair<double, double> interval, const std::vector<EvalStep>& eval_steps, double interval_size, bool lipschitz) {
    std::vector<EvalStep> es_copy = eval_steps;
    return fct::estimate_composite_range(eb, interval, es_copy, interval_size, lipschitz);
}

void bind_ciphertext(py::module_& m) {
    py::class_<fct::Ciphertext>(m, "Ciphertext")
        .def(py::init<>())
        .def(py::init<double, int, double, double>(),
            py::arg("x"), py::arg("size"),
            py::arg("err") = 0.0, py::arg("scale") = 1.0)

        .def(py::init<double, double, double, double, double, bool>(),
            py::arg("start"), py::arg("end"),
            py::arg("interval_size") = 1.0 / 16384.0,
            py::arg("err") = 0.0, py::arg("scale") = 1.0,
            py::arg("encode") = true)

        .def(py::init([](std::vector<double> x, double err, double scale, bool encode) {
                return fct::Ciphertext(x, err, scale, encode);
            }),
            py::arg("x"), py::arg("err") = 0.0,
            py::arg("scale") = 1.0, py::arg("encode") = true)

        .def(py::init([](std::vector<std::pair<double, double>> intervals,
                        double err,
                        double scale,
                        bool encode) {
                return fct::Ciphertext(intervals, err, scale, encode);
            }),
            py::arg("intervals"), py::arg("err") = 0.0,
            py::arg("scale") = 1.0, py::arg("encode") = true)

        .def(py::init([](std::vector<double> x,
                        std::vector<double> ct_high,
                        std::vector<double> ct_low,
                        double err,
                        double scale) {
                return fct::Ciphertext(x, ct_high, ct_low, err, scale);
            }),
            py::arg("x"), py::arg("ct_high"), py::arg("ct_low"),
            py::arg("err") = 0.0, py::arg("scale") = 1.0)

        .def_readwrite("pt", &fct::Ciphertext::pt)
        .def_readwrite("err", &fct::Ciphertext::err)
        .def_readwrite("ct_high", &fct::Ciphertext::ct_high)
        .def_readwrite("ct_low", &fct::Ciphertext::ct_low)
        .def_readwrite("scale", &fct::Ciphertext::scale)
        .def("set_error", &fct::Ciphertext::set_error, py::arg("err"));
}

void bind_ciphertext_ops(py::module_& m) {
    m.def("add", &fct::add,
        py::arg("ct1"), py::arg("ct2"));

    m.def("add_pt", &fct::add_pt,
        py::arg("ct"), py::arg("scalar"));

    m.def("mult", &fct::mult,
        py::arg("ct1"), py::arg("ct2"), py::arg("eb"));

    m.def("square", &fct::square,
        py::arg("ct"), py::arg("eb"));

    m.def("mult_pt", &fct::mult_pt,
        py::arg("scalar"), py::arg("ct"), py::arg("eb"),
        py::arg("rescale") = true);

    m.def("dec", &fct::dec, py::arg("ct"));
    m.def("dec_interval", &fct::dec_interval, py::arg("ct"));
    m.def("dec_intervals", &fct::dec_intervals, py::arg("ct"));
}

PYBIND11_MODULE(EBModule, m) {
    m.doc() = "pybind11 bindings for error_bound, parse, and fake ciphertext";

    py::class_<DecompCache>(m, "DecompCache")
        .def(py::init<const std::string&>(), py::arg("cache_file"))
        .def("restore_EvalStep", &DecompCache::restore_EvalStep, py::arg("coeff"));

    m.def("restore_EvalStep", &restore_EvalStep,
          py::arg("coeff"), py::arg("cache_file"));

    m.def("evaluate", &evaluate,
          py::arg("x"), py::arg("es"), py::arg("eb"));

    m.def("estimate_range", &estimate_range,
          py::arg("eb"), py::arg("interval"), py::arg("es"), py::arg("interval_size"), py::arg("lipschitz") = false);

    m.def("estimate_lipschitz_range", &estimate_composite_range,
          py::arg("eb"), py::arg("interval"), py::arg("eval_steps"), py::arg("interval_size"), py::arg("lipschitz"));
          
    m.def("bit_diff", &bit_diff, py::arg("a"), py::arg("b"));
    m.def("clearNum", &clearNum, py::arg("a"));

    py::class_<ErrBound>(m, "ErrBound")
        .def(py::init<double, int, int, int>(),
             py::arg("sigma"), py::arg("N"), py::arg("h"), py::arg("s"))
        .def_readwrite("Bc", &ErrBound::Bc)
        .def_readwrite("Bs", &ErrBound::Bs)
        .def_readwrite("scale", &ErrBound::scale)
        .def("B_clean", &ErrBound::B_clean,
             py::arg("sigma"), py::arg("N"), py::arg("h"))
        .def("B_scale", &ErrBound::B_scale,
             py::arg("N"), py::arg("h"));

    bind_ciphertext(m);
    bind_ciphertext_ops(m);

    py::class_<step>(m, "step")
        .def(py::init<>())
        .def_readwrite("op", &step::op)
        .def_readwrite("key1", &step::key1)
        .def_readwrite("key2", &step::key2)
        .def_readwrite("save_key", &step::save_key);

    py::class_<EvalStep>(m, "EvalStep")
        .def(py::init<>())
        .def_readwrite("term_count", &EvalStep::term_count)
        .def_readwrite("coeff_count", &EvalStep::coeff_count)
        .def_readwrite("powers", &EvalStep::powers)
        .def_readwrite("coeffs", &EvalStep::coeffs)
        .def_readwrite("terms", &EvalStep::terms)
        .def_readwrite("eval_step", &EvalStep::eval_step)
        .def("get_value", &EvalStep::get_value, py::arg("key"))
        .def("add_value", &EvalStep::add_value,
             py::arg("key_header"), py::arg("key_number"), py::arg("data"))
        .def("add_step", &EvalStep::add_step,
             py::arg("key1"), py::arg("key2"), py::arg("op"))
        .def("print_step", &EvalStep::print_step);
}
