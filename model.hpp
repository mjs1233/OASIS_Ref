#ifndef OASISREF_MODEL_HPP
#define OASISREF_MODEL_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>

struct body_param {

    body_param() = default;

    /*
     * @param : height (cm)
     * @param : weight (kg)
     */
    body_param(double height, double weight) {
        //Kurazumi 1994
        area = 100.315 * std::pow(weight, 0.383) * std::pow(height, 0.693) * 1e-4;


        mass_skin = weight * 0.05;
        mass_core = weight * 0.95;
    }

    double area {};
    double mass_skin {};
    double mass_core {};
};


static constexpr double pr1 = 36.6;
static constexpr double pr2 = 34.1;
static constexpr double pr5 = 75.0;
static constexpr double pr6 = 0.0;

constexpr std::array pr3_values {
    100.0, 80.0, 60.0, 40.0, 20.0, 10.0, 5.0
};

constexpr std::array pr4_values {
    12.6, 10.08, 7.56, 5.04, 2.52, 1.26, 0.63, 0.315, 0.1575, 0.07875
};

constexpr std::array pr7_values {
    250.0, 200.0, 150.0, 100.0, 50.0, 25.0, 12.5
};

struct param {
public:
    param() = default;

    param(double pr3,double pr4, double pr7) :
    value {pr1,pr2,pr3,pr4,pr5,pr6,pr7} {}

    std::array<double,7> value{};

    double operator[](std::size_t idx) const {
        if (idx < 1 || idx > value.size()) {
            throw std::out_of_range("param index must be 1..7");
        }
        return value[idx - 1];
    }
};

inline auto sat_vapor_pressure_mmHg = [](double temp)-> double {
    return std::exp(18.6686 - (4030.183 / (temp + 235.0)));
};

class model {
public:
    model() = default;

    model(param pr) :
    m_pr(pr),
    m_temp_core{pr[1]},
    m_temp_skin{pr[2]} {}

    //@param : hum_air (normalized 0-1)
    double step(double temp_air,double hum_air, double q_res, double m_total, double w,double dt = 60.0) {
        hum_air = std::clamp(hum_air, 0.0, 1.0);
        const double _core_signal =
            std::max(0.0, m_temp_core - m_pr[1]);

        const double _skin_signal =
            std::max(0.0, m_temp_skin - m_pr[2]);

        //전도 열전달
        constexpr double _k_min = 5.28;
        double _q_cond = _k_min * (m_temp_core - m_temp_skin);
        double _v_blo = (m_pr[4] + m_pr[5] * _core_signal) * (1.0 / 60.0);

        //혈류 열전달
        constexpr double _c_blo = 1.163;
        double _q_blo = _c_blo * _v_blo * (m_temp_core - m_temp_skin);

        //대류 + 복사
        constexpr double _h_conv = 4.3;
        constexpr double _h_rad = 5.23;
        constexpr double _f_cl = 0.53;
        double _q_conv_q_rad = (_h_conv + _h_rad) * (m_temp_skin - temp_air) * _f_cl;

        //땀 분비량
        double _m_rsw =
            m_pr[7] * _core_signal + m_pr[3] * _core_signal * _skin_signal * 1.0/1000.0 * 1.0/60.0;

        double _p_skin = sat_vapor_pressure_mmHg(m_temp_skin);
        constexpr double _f_pcl = 0.73;

        double _e_max = 2.2 * _h_conv * (_p_skin - hum_air * sat_vapor_pressure_mmHg(temp_air)) * _f_pcl;
        _e_max = std::max(_e_max, 0.0);

        double _q_rsw_raw = 0.7 * _m_rsw * std::pow(2.0,_skin_signal / 3.0);
        _q_rsw_raw = std::max(_q_rsw_raw, 0.0);

        double _q_diff;
        double _q_rsw;

        if (_q_rsw_raw > _e_max) {
            _q_rsw  = _e_max;
            _q_diff = 0.0;
        }
        else {
            _q_rsw  = _q_rsw_raw;
            _q_diff = 0.06 * _e_max;
        }

        auto kcal_to_J = [](double kcal) -> double {
            return kcal * 4184.0;
        };

        constexpr double _c_skin = kcal_to_J(0.97);
        constexpr double _c_core = kcal_to_J(0.97);

        //skin_temp, core_temp gradient
        double _d_temp_skin = (_q_cond + _q_blo - _q_diff - _q_rsw - _q_conv_q_rad) * body.area * dt / (body.mass_skin * _c_skin);
        double _d_temp_core = (m_total - w - q_res- _q_cond- _q_blo) * body.area * dt / (body.mass_core * _c_core);

        //update
        m_temp_skin_delta = _d_temp_skin;
        m_temp_core += _d_temp_core;
        m_temp_skin += _d_temp_skin;

        return m_temp_skin;
    }

    void accumulate_error(double measured_temp_skin_delta) {
        m_error_sum += std::abs(measured_temp_skin_delta - m_temp_skin_delta);
    }

    [[nodiscard]]
    double mean_error(size_t step_count) const {
        if (step_count > 0) {
            return m_error_sum / step_count;
        }
        return std::numeric_limits<double>::infinity();
    }



    [[nodiscard]]
    double core_temperature() const {
        return m_temp_core;
    }

    [[nodiscard]]
    double skin_temperature() const {
        return m_temp_skin;
    }

    void set_body_param(body_param param) {
        body = param;
    }

private:
    param m_pr {};
    double m_temp_core {};
    double m_temp_skin {};

    double m_temp_skin_delta{};
    double m_error_sum {};

    body_param body {};

};

#endif //OASISREF_MODEL_HPP
