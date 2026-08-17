#include <iostream>
#include <array>
#include <cmath>
#include <cinttypes>
#include <vector>
#include "model.hpp"

struct error_record {
    size_t model_index;
    double error;
};

int main() {

    //base input variable
    double _height_cm = 175.0;
    double _weight_kg = 70.0;
    body_param body {_height_cm, _weight_kg};
    model::set_body_param(body);

    std::vector<model> models;
    models.resize(490);

    std::vector<error_record> _errors;
    _errors.reserve(32);

    //model prepare logic...
    size_t idx = 0;
    for (const auto pr3 : pr3_values) {
        for (const auto pr4 : pr4_values) {
            for (const auto pr7 : pr7_values) {
                models[idx++] = model{param{pr3, pr4, pr7}};
            }
        }
    }

    bool _first_flag = true;
    double _prev_temp_skin = 0;
    size_t _step_count = 0;

    while (true) {

        //input variable
        double _temp_air = 0;
        double _temp_skin = 0;
        double _extern_hum = 0;
        double _basal_kcal_per_min = 0;
        double _exercise_kcal_per_min = 0;



        const auto kcal_min_to_W_m2 =
               [area = body.area](double kcal_per_min) -> double {
                   constexpr double kcal_to_j = 4184.0;
                   return kcal_per_min * kcal_to_j / 60.0 / area;
        };

        //Kcal/min 기준으로 계산.. f: Kcal/min -> (W/m^2)
        double _m_basal = kcal_min_to_W_m2(_basal_kcal_per_min); // W/m²
        double _m_ex    = kcal_min_to_W_m2(_exercise_kcal_per_min); // W/m²

        double _m_total =
            _m_basal + _m_ex;

        constexpr double _w_eff = 0.4;

        double _w =
            _m_ex * _w_eff;

        auto normalize = [](double percentage) -> double {
            return percentage / 100.0;
        };

        const double _normalized_extern_hum = std::clamp(normalize(_extern_hum), 0.0,1.0);

        //호흡으로 빠져나가는 열
        double _q_res = 0.0023 * _m_total * (44 - _normalized_extern_hum * sat_vapor_pressure_mmHg(_temp_air));

        if (_first_flag) {
            _first_flag = false;
            _prev_temp_skin = _temp_skin;
            continue;
        }

        const double measured_skin_delta = _temp_skin  - _prev_temp_skin;

        //run models....
        for (auto& _model : models) {
            _model.step(_temp_air, _normalized_extern_hum, _q_res, _m_total, _w);
            _model.accumulate_error(measured_skin_delta);
        }
        _step_count++;
        size_t _candidate_count = 0;
        size_t try_count = 0;
        constexpr double _default_distance = 0.001;
        constexpr double _distance_increment = 0.0001;
        double _select_distance = _default_distance;
        _errors.clear();
        do {
            try_count++;
            size_t index = 0;
            for (auto& _model : models) {
                double _error = _model.mean_error(_step_count);
                if (_error < _select_distance) {
                    _candidate_count++;
                    _errors.push_back({.model_index = index, .error = _error});
                }
                index++;
            }
            _select_distance += _distance_increment;

        } while (_candidate_count == 0);

        double _estimated_core_temperature = 0;
        for (const auto& _error : _errors) {
            _estimated_core_temperature += models[_error.model_index].core_temperature();
        }

        _estimated_core_temperature /= static_cast<double>(_candidate_count);

        //report result
        std::cout
        << "estimated core temperature: "
        << _estimated_core_temperature << " C\n"
        << "selected models: "
        << _candidate_count << " / "
        << 490 << '\n'
        << "effective threshold: "
        << (_select_distance - _distance_increment) << '\n';

        _prev_temp_skin = _temp_skin;

        std::cin.get();
    }
    return 0;
}
