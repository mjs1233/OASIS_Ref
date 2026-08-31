#include "simulation.hpp"
#include "test_history.hpp"

int main() {
    test_history history{"../res/dataset0/P1_C1.csv"};
    simulation sim{history.height,history.weight};

    uint32_t min = 0;
    std::cout << history.height << std::endl;
    for (min = 0; min < history.records.size(); min++) {
        const test_history::record rec = history.get(min);
        double est = sim.step(rec.temp_air,rec.temp_skin,rec.hum,history.basal_kcal_per_min,rec.exercise_kcal_per_min);
        std::cout << "estimated " <<est << " // measured " << rec.temp_core << "\n";

    }
}