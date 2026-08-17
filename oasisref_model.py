"""Two-node core-temperature estimator based on the OASISRef C++ model.

Inputs passed to ``CoreTemperatureEstimator.update`` are assumed to be sampled
at a fixed interval (60 seconds by default). Metabolic inputs are in kcal/min.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import product
import math


PR1 = 36.6
PR2 = 34.1
PR5 = 75.0
PR6 = 0.0

PR3_VALUES = (100.0, 80.0, 60.0, 40.0, 20.0, 10.0, 5.0)
PR4_VALUES = (
    12.6,
    10.08,
    7.56,
    5.04,
    2.52,
    1.26,
    0.63,
    0.315,
    0.1575,
    0.07875,
)
PR7_VALUES = (250.0, 200.0, 150.0, 100.0, 50.0, 25.0, 12.5)

DEFAULT_DISTANCE = 0.001
DISTANCE_INCREMENT = 0.0001
KCAL_TO_J = 4184.0


@dataclass(frozen=True, slots=True)
class BodyParam:
    height_cm: float
    weight_kg: float
    area: float = field(init=False)
    mass_skin: float = field(init=False)
    mass_core: float = field(init=False)

    def __post_init__(self) -> None:
        if not math.isfinite(self.height_cm) or self.height_cm <= 0.0:
            raise ValueError("height_cm must be a positive finite value")
        if not math.isfinite(self.weight_kg) or self.weight_kg <= 0.0:
            raise ValueError("weight_kg must be a positive finite value")

        # Kurazumi, 1994. Height is in cm and weight is in kg.
        area = (
            100.315
            * self.weight_kg**0.383
            * self.height_cm**0.693
            * 1.0e-4
        )
        object.__setattr__(self, "area", area)
        object.__setattr__(self, "mass_skin", self.weight_kg * 0.05)
        object.__setattr__(self, "mass_core", self.weight_kg * 0.95)


@dataclass(frozen=True, slots=True)
class Param:
    pr3: float
    pr4: float
    pr7: float
    pr1: float = PR1
    pr2: float = PR2
    pr5: float = PR5
    pr6: float = PR6


def sat_vapor_pressure_mmhg(temp_c: float) -> float:
    return math.exp(18.6686 - 4030.183 / (temp_c + 235.0))


@dataclass(slots=True)
class ThermalModel:
    body: BodyParam
    param: Param
    temp_core: float = field(init=False)
    temp_skin: float = field(init=False)
    temp_skin_delta: float = 0.0
    error_sum: float = 0.0
    error_count: int = 0

    def __post_init__(self) -> None:
        self.temp_core = self.param.pr1
        self.temp_skin = self.param.pr2

    def step(
        self,
        temp_air: float,
        hum_air: float,
        q_res: float,
        m_total: float,
        work: float,
        dt_seconds: float = 60.0,
    ) -> float:
        if not math.isfinite(dt_seconds) or dt_seconds <= 0.0:
            raise ValueError("dt_seconds must be a positive finite value")

        hum_air = min(max(hum_air, 0.0), 1.0)

        core_signal = max(0.0, self.temp_core - self.param.pr1)
        skin_signal = max(0.0, self.temp_skin - self.param.pr2)

        k_min = 5.28
        q_cond = k_min * (self.temp_core - self.temp_skin)

        # pr6 is removed in the referenced method, so its denominator is 1.
        v_blo = (self.param.pr4 + self.param.pr5 * core_signal) / 60.0

        c_blo = 1.163
        q_blo = c_blo * v_blo * (self.temp_core - self.temp_skin)

        h_conv = 4.3
        h_rad = 5.23
        f_cl = 0.53
        q_conv_rad = (
            (h_conv + h_rad) * (self.temp_skin - temp_air) * f_cl
        )

        m_rsw = (
            self.param.pr7 * core_signal
            + self.param.pr3
            * core_signal
            * skin_signal
            / 1000.0
            / 60.0
        )

        p_skin = sat_vapor_pressure_mmhg(self.temp_skin)
        p_air = hum_air * sat_vapor_pressure_mmhg(temp_air)
        f_pcl = 0.73
        e_max = max(0.0, 2.2 * h_conv * (p_skin - p_air) * f_pcl)

        q_rsw_raw = max(
            0.0,
            0.7 * m_rsw * 2.0 ** (skin_signal / 3.0),
        )

        if q_rsw_raw > e_max:
            q_rsw = e_max
            q_diff = 0.0
        else:
            q_rsw = q_rsw_raw
            q_diff = 0.06 * e_max

        c_skin = 0.97 * KCAL_TO_J
        c_core = 0.97 * KCAL_TO_J

        d_temp_skin = (
            (q_cond + q_blo - q_diff - q_rsw - q_conv_rad)
            * self.body.area
            * dt_seconds
            / (self.body.mass_skin * c_skin)
        )
        d_temp_core = (
            (m_total - work - q_res - q_cond - q_blo)
            * self.body.area
            * dt_seconds
            / (self.body.mass_core * c_core)
        )

        # Both deltas were calculated from the same previous state.
        self.temp_skin_delta = d_temp_skin
        self.temp_core += d_temp_core
        self.temp_skin += d_temp_skin
        return self.temp_skin

    def accumulate_error(self, measured_temp_skin_delta: float) -> None:
        self.error_sum += abs(
            measured_temp_skin_delta - self.temp_skin_delta
        )
        self.error_count += 1

    @property
    def mean_error(self) -> float:
        if self.error_count == 0:
            return math.inf
        return self.error_sum / self.error_count


@dataclass(frozen=True, slots=True)
class EstimationResult:
    estimated_core_temperature: float
    selected_model_count: int
    total_model_count: int
    effective_threshold: float
    selected_model_indices: tuple[int, ...]


class CoreTemperatureEstimator:
    """Maintains all 490 model states and filters them at each sample."""

    def __init__(
        self,
        height_cm: float,
        weight_kg: float,
        *,
        dt_seconds: float = 60.0,
        default_distance: float = DEFAULT_DISTANCE,
        distance_increment: float = DISTANCE_INCREMENT,
    ) -> None:
        if not math.isfinite(dt_seconds) or dt_seconds <= 0.0:
            raise ValueError("dt_seconds must be a positive finite value")
        if not math.isfinite(default_distance) or default_distance < 0.0:
            raise ValueError("default_distance must be finite and nonnegative")
        if not math.isfinite(distance_increment) or distance_increment <= 0.0:
            raise ValueError("distance_increment must be positive and finite")

        self.body = BodyParam(height_cm, weight_kg)
        self.dt_seconds = dt_seconds
        self.default_distance = default_distance
        self.distance_increment = distance_increment
        self.previous_measured_skin: float | None = None

        self.models = [
            ThermalModel(self.body, Param(pr3, pr4, pr7))
            for pr3, pr4, pr7 in product(
                PR3_VALUES,
                PR4_VALUES,
                PR7_VALUES,
            )
        ]
        if len(self.models) != 490:
            raise RuntimeError("unexpected parameter combination count")

    def _kcal_min_to_w_m2(self, kcal_per_min: float) -> float:
        return kcal_per_min * KCAL_TO_J / 60.0 / self.body.area

    def update(
        self,
        *,
        temp_air: float,
        measured_skin_temp: float,
        humidity_percent: float,
        basal_kcal_per_min: float,
        exercise_kcal_per_min: float,
    ) -> EstimationResult | None:
        """Consume one sample; the first sample establishes the skin baseline."""

        values = (
            temp_air,
            measured_skin_temp,
            humidity_percent,
            basal_kcal_per_min,
            exercise_kcal_per_min,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("all sample inputs must be finite")

        humidity = min(max(humidity_percent / 100.0, 0.0), 1.0)

        if self.previous_measured_skin is None:
            self.previous_measured_skin = measured_skin_temp
            return None

        measured_skin_delta = (
            measured_skin_temp - self.previous_measured_skin
        )

        m_basal = self._kcal_min_to_w_m2(basal_kcal_per_min)
        m_ex = self._kcal_min_to_w_m2(exercise_kcal_per_min)
        m_total = m_basal + m_ex
        work = m_ex * 0.4
        q_res = (
            0.0023
            * m_total
            * (44.0 - humidity * sat_vapor_pressure_mmhg(temp_air))
        )

        for model in self.models:
            model.step(
                temp_air,
                humidity,
                q_res,
                m_total,
                work,
                self.dt_seconds,
            )
            model.accumulate_error(measured_skin_delta)

        self.previous_measured_skin = measured_skin_temp

        finite_errors = [
            model.mean_error
            for model in self.models
            if math.isfinite(model.mean_error)
        ]
        if not finite_errors:
            raise RuntimeError("all model errors are non-finite")

        threshold = self.default_distance
        selected: list[int] = []
        while not selected:
            selected = [
                index
                for index, model in enumerate(self.models)
                if math.isfinite(model.mean_error)
                and model.mean_error <= threshold
            ]
            if not selected:
                threshold += self.distance_increment

        estimated_core = sum(
            self.models[index].temp_core for index in selected
        ) / len(selected)

        return EstimationResult(
            estimated_core_temperature=estimated_core,
            selected_model_count=len(selected),
            total_model_count=len(self.models),
            effective_threshold=threshold,
            selected_model_indices=tuple(selected),
        )


if __name__ == "__main__":
    estimator = CoreTemperatureEstimator(
        height_cm=175.0,
        weight_kg=70.0,
    )

    print(
        "Enter: air_C skin_C humidity_percent "
        "basal_kcal_per_min exercise_kcal_per_min"
    )
    while True:
        try:
            line = input("> ").strip()
        except EOFError:
            break
        if not line:
            continue

        try:
            air, skin, humidity, basal, exercise = map(float, line.split())
            result = estimator.update(
                temp_air=air,
                measured_skin_temp=skin,
                humidity_percent=humidity,
                basal_kcal_per_min=basal,
                exercise_kcal_per_min=exercise,
            )
        except ValueError as error:
            print(f"invalid input: {error}")
            continue

        if result is None:
            print("initial skin-temperature sample stored")
            continue

        print(
            f"estimated core temperature: "
            f"{result.estimated_core_temperature:.6f} C\n"
            f"selected models: {result.selected_model_count} / "
            f"{result.total_model_count}\n"
            f"effective threshold: {result.effective_threshold:.6f}"
        )
