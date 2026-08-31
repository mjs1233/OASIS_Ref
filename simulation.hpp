#ifndef OASISREF_SIMULATION_HPP
#define OASISREF_SIMULATION_HPP
#include <iostream>

#include "model.hpp"
#include <vector>

class simulation {
private:
    struct error_record {
        size_t model_index;
        double error;
    };

public:
    simulation(float height,float weight) :
    body(height,weight) {

        models.reserve(model_count);
        for (double p3 : pr3_values) {
            for (double p4 : pr4_values) {
                for (double p7 : pr7_values) {
                    models.emplace_back(param {p3, p4, p7});
                }
            }
        }

        for (auto& m : models) {
            m.set_body_param(body);
        }

        errors.reserve(32);
    }

    double step(double t_air,double t_skin, double hum_ext, double basal_kmin, double exercise_kmin) {

        double _temp_air = t_air;
        double _temp_skin = t_skin;
        double _extern_hum = hum_ext;
        double _basal_kcal_per_min = basal_kmin;
        double _exercise_kcal_per_min = exercise_kmin;

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

        if (first_flag) {
            first_flag = false;
            prev_temp_skin = _temp_skin;
            return 0;
        }

        const double measured_skin_delta = _temp_skin  - prev_temp_skin;

        //run models....
        for (auto& _model : models) {
            _model.step(_temp_air, _normalized_extern_hum, _q_res, _m_total, _w);
            _model.accumulate_error(measured_skin_delta);
        }
        step_count++;
        size_t _candidate_count = 0;
        size_t try_count = 0;
        constexpr double _default_distance = 0.001;
        constexpr double _distance_increment = 0.0001;
        double _select_distance = _default_distance;
        errors.clear();
        do {
            try_count++;
            size_t index = 0;
            for (auto& _model : models) {
                double _error = _model.mean_error(step_count);

                if (std::isfinite(_error) &&_error <= _select_distance) {
                    _candidate_count++;
                    errors.push_back({.model_index = index, .error = _error});
                }
                index++;
            }
            _select_distance += _distance_increment;

        } while (_candidate_count == 0);

        double _estimated_core_temperature = 0;
        for (const auto& _error : errors) {
            _estimated_core_temperature += models[_error.model_index].core_temperature();
        }

        _estimated_core_temperature /= static_cast<double>(_candidate_count);

        //report result
        std::cout
        << "estimated core temperature: "
        << _estimated_core_temperature << " C\n"
        << "selected models: "
        << _candidate_count << " / "
        << model_count << '\n'
        << "effective threshold: "
        << (_select_distance - _distance_increment) << '\n';

        prev_temp_skin = _temp_skin;
        return _estimated_core_temperature;
    }

private:

    static constexpr std::size_t model_count = pr3_values.size() * pr4_values.size() * pr7_values.size();
    body_param body {};
    std::vector<model> models;

    std::vector<error_record> errors;
    bool first_flag = true;
    double prev_temp_skin = 0;
    size_t step_count = 0;

};
#endif //OASISREF_SIMULATION_HPP