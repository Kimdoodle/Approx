/* #include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../include/error_bound.h"
#include "../include/parse.h"

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace py = pybind11;

PYBIND11_MODULE(EBModule, m) {
    m.doc() = "pybind11 bindings for error_bound";

    m.def("bit_diff", &bit_diff,
          py::arg("a"), py::arg("b"));

    m.def("clearNum", &clearNum,
          py::arg("a"));

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
        // cal_bound 는 제외

    py::class_<fct::Ciphertext>(m, "Ciphertext")
        .def(py::init<>())
        .def(py::init<double, int, double, double>(),
             py::arg("x"), py::arg("size"),
             py::arg("err") = 0.0, py::arg("scale") = 1.0)
        .def(py::init<std::vector<double>&, double, double>(),
             py::arg("x"), py::arg("err") = 0.0, py::arg("scale") = 1.0)
        .def(py::init<std::vector<double>&, std::vector<double>&, std::vector<double>&, double, double>(),
             py::arg("x"), py::arg("ct_high"), py::arg("ct_low"),
             py::arg("err") = 0.0, py::arg("scale") = 1.0)
        .def_readwrite("pt", &fct::Ciphertext::pt)
        .def_readwrite("err", &fct::Ciphertext::err)
        .def_readwrite("ct_high", &fct::Ciphertext::ct_high)
        .def_readwrite("ct_low", &fct::Ciphertext::ct_low)
        .def_readwrite("scale", &fct::Ciphertext::scale)
        .def("set_error", &fct::Ciphertext::set_error,
             py::arg("err"));

    m.def("ct_add", &fct::add,
          py::arg("ct1"), py::arg("ct2"));

    m.def("ct_add_pt", &fct::add_pt,
          py::arg("ct"), py::arg("scalar"));

    m.def("ct_mult", &fct::mult,
          py::arg("ct1"), py::arg("ct2"), py::arg("eb"));

    m.def("ct_square", &fct::square,
          py::arg("ct"), py::arg("eb"));

    m.def("ct_mult_pt", &fct::mult_pt,
          py::arg("scalar"), py::arg("ct"), py::arg("eb"),
          py::arg("rescale") = true);

    m.def("ct_dec", &fct::dec,
          py::arg("ct"));
} */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../include/parse.h"
#include "../include/error_bound.h"

#include <memory>
#include <stdexcept>
#include <string>
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

fct::Ciphertext evaluate_fct(fct::Ciphertext x, const EvalStep& es, ErrBound& eb) {
    EvalStep es_copy = es;
    return fct::evaluate(x, es_copy, eb);
}

fct::Ciphertext evaluate_fct_cl(fct::Ciphertext x, ErrBound& eb) {
    return fct::evaluate_cl(x, eb);
}

PYBIND11_MODULE(EBModule, m) {
    m.doc() = "pybind11 bindings for error_bound + parse";

    py::class_<DecompCache>(m, "DecompCache")
        .def(py::init<const std::string&>(), py::arg("cache_file"))
        .def("restore_EvalStep", &DecompCache::restore_EvalStep, py::arg("coeff"));

    m.def("restore_EvalStep", &restore_EvalStep,
          py::arg("coeff"), py::arg("cache_file"));

    m.def("evaluate_fct", &evaluate_fct,
          py::arg("x"), py::arg("es"), py::arg("eb"));

    m.def("evaluate_fct_cl", &evaluate_fct_cl,
          py::arg("x"), py::arg("eb"));
          
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

    py::class_<fct::Ciphertext>(m, "Ciphertext")
        .def(py::init<>())
        .def(py::init<double, int, double, double>(),
             py::arg("x"), py::arg("size"),
             py::arg("err") = 0.0, py::arg("scale") = 1.0)
        .def(py::init([](std::vector<double> x, double err, double scale, bool encode) {
                return fct::Ciphertext(x, err, scale, encode);
             }),
             py::arg("x"), py::arg("err") = 0.0,
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

    m.def("ct_add", &fct::add, py::arg("ct1"), py::arg("ct2"));
    m.def("ct_add_pt", &fct::add_pt, py::arg("ct"), py::arg("scalar"));
    m.def("ct_mult", &fct::mult, py::arg("ct1"), py::arg("ct2"), py::arg("eb"));
    m.def("ct_square", &fct::square, py::arg("ct"), py::arg("eb"));
    m.def("ct_mult_pt", &fct::mult_pt,
          py::arg("scalar"), py::arg("ct"), py::arg("eb"),
          py::arg("rescale") = true);
    m.def("ct_dec", &fct::dec, py::arg("ct"));

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