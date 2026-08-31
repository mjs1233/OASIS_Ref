//
// Created by user on 26. 8. 25..
//

#ifndef OASISREF_TESTHISTORY_HPP
#define OASISREF_TESTHISTORY_HPP
#include <filesystem>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <stdexcept>
#include <cmath>
#include <cctype>
#include <cstdint>

class test_history {
public:

    struct record {
        double exercise_kcal_per_min {};
        double temp_air {};
        double hum {};
        double temp_skin {};
        double temp_core {};
    };

    explicit test_history(std::filesystem::path const& path) {
        const std::string stem = path.stem().string();
        id        = parse_prefixed_number(stem, 'P');
        condition = parse_prefixed_number(stem, 'C');

        std::ifstream file(path);
        if (!file.is_open()) {
            throw std::runtime_error("TestHistory: failed to open file: " + path.string());
        }

        std::string header_line;
        if (!std::getline(file, header_line)) {
            throw std::runtime_error("TestHistory: empty file: " + path.string());
        }

        const auto col_index = parse_header(header_line);

        auto require_col = [&](const char* name) -> std::size_t {
            auto it = col_index.find(name);
            if (it == col_index.end()) {
                throw std::runtime_error(
                    std::string("TestHistory: missing column '") + name + "' in " + path.string());
            }
            return it->second;
        };

        const std::size_t idx_t_min      = require_col("t_min");
        const std::size_t idx_age        = require_col("age");
        const std::size_t idx_weight     = require_col("weight_kg");
        const std::size_t idx_height     = require_col("height_cm");
        const std::size_t idx_basal      = require_col("basal_kcal_per_min");
        const std::size_t idx_exercise   = require_col("exercise_kcal_per_min");
        const std::size_t idx_air_temp   = require_col("air_temp_c");
        const std::size_t idx_hum        = require_col("ext_hum_pct");
        const std::size_t idx_tcore      = require_col("tcore_c");
        const std::size_t idx_tskin      = require_col("tskin_head_c");

        struct raw_row {
            uint32_t min_index;
            record rec;
        };
        std::vector<raw_row> rows;

        bool first_row = true;
        std::string line;
        while (std::getline(file, line)) {
            if (line.empty()) {
                continue;
            }
            const auto fields = split_csv_line(line);

            if (first_row) {
                age    = static_cast<uint32_t>(std::llround(to_double(fields, idx_age)));
                weight = to_double(fields, idx_weight);
                height = to_double(fields, idx_height);
                basal_kcal_per_min = to_double(fields, idx_basal);
                first_row = false;
            }

            const double t_min = to_double(fields, idx_t_min);
            const auto min_index = static_cast<uint32_t>(std::llround(t_min));

            record rec {};
            rec.exercise_kcal_per_min = to_double(fields, idx_exercise);
            rec.temp_air              = to_double(fields, idx_air_temp);
            rec.hum                   = to_double(fields, idx_hum);
            rec.temp_skin             = to_double(fields, idx_tskin);
            rec.temp_core             = to_double(fields, idx_tcore);

            rows.push_back({min_index, rec});
        }

        if (rows.empty()) {
            return;
        }

        uint32_t max_index = 0;
        for (const auto& r : rows) {
            max_index = std::max(max_index, r.min_index);
        }

        records.assign(static_cast<std::size_t>(max_index) + 1, record{});
        for (const auto& r : rows) {
            records[r.min_index] = r.rec;
        }
    }


    record get(uint32_t min) {
        if (min >= records.size()) {
            return {};
        }
        return records[min];
    }

    uint32_t id {};
    uint32_t condition {};
    uint32_t age {};
    float weight{};
    float height{};
    double basal_kcal_per_min {};

    std::vector<record> records;

private:
    static std::unordered_map<std::string, std::size_t> parse_header(std::string const& line) {
        std::unordered_map<std::string, std::size_t> map;
        const auto fields = split_csv_line(line);
        for (std::size_t i = 0; i < fields.size(); ++i) {
            map[fields[i]] = i;
        }
        return map;
    }

    static std::vector<std::string> split_csv_line(std::string const& line) {
        std::vector<std::string> fields;
        std::stringstream ss(line);
        std::string field;
        while (std::getline(ss, field, ',')) {
            fields.push_back(field);
        }
        return fields;
    }

    static double to_double(std::vector<std::string> const& fields, std::size_t idx) {
        if (idx >= fields.size() || fields[idx].empty()) {
            return 0.0;
        }
        try {
            return std::stod(fields[idx]);
        } catch (...) {
            return 0.0;
        }
    }

    // 파일명 stem(예: "P17_C2")에서 주어진 prefix 문자 바로 뒤에 오는
    // 연속된 숫자를 찾아 파싱한다. 예: prefix='P' -> 17, prefix='C' -> 2
    static uint32_t parse_prefixed_number(std::string const& stem, char prefix) {
        for (std::size_t i = 0; i < stem.size(); ++i) {
            if (stem[i] == prefix && i + 1 < stem.size() &&
                std::isdigit(static_cast<unsigned char>(stem[i + 1]))) {
                std::size_t j = i + 1;
                while (j < stem.size() && std::isdigit(static_cast<unsigned char>(stem[j]))) {
                    ++j;
                }
                return static_cast<uint32_t>(std::stoul(stem.substr(i + 1, j - i - 1)));
            }
        }
        throw std::runtime_error(
            "TestHistory: cannot parse '" + std::string(1, prefix) + "<number>' from filename: " + stem);
    }
};

#endif //OASISREF_TESTHISTORY_HPP