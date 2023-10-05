import numpy as np
import pandas as pd
from rocketpy import Environment, SolidMotor, Rocket, Flight
from sympy import *

MOTOR_LENGTH = 1.239  # m

# The number of grains in the rocket motor
NUM_GRAINS = 6

OUTER_GRAIN_DIAMETER = 0.088  # m
OUTER_GRAIN_RADIUS = OUTER_GRAIN_DIAMETER / 2  # m

INNER_DIAMETER_LIST = [0.058, 0.058, 0.058, 0.058, 0.058, 0.058]  # m
INNER_DIAMETER_AVG = sum(INNER_DIAMETER_LIST) / len(INNER_DIAMETER_LIST)  # m
INNER_GRAIN_RADIUS = INNER_DIAMETER_AVG / 2  # m

GRAIN_HEIGHT = 0.20483  # m

PROPELLANT_MASS = 9.021  # kg
Empty_Mass = 5.4
Motor_length = 1239
Theory_thinckness = 0.005  # unit:m
out_radius = 0.049
Dry_inertia = [0.33, 0.33, 0.012, 0, 0, 0]
Distance_each_grain = 0.001

NOZZLE_DIAMETER = 0.0275  # m
NOZZLE_RADIUS = NOZZLE_DIAMETER / 2  # m
grains_center_of_mass_position = -0.6
center_of_dry_mass_position = -0.7
# The radius of the throat of the motor nozzle in metres (inner part of the nozzle)
THROAT_DIAMETER = 0.0146  # m
THROAT_RADIUS = THROAT_DIAMETER / 2  # m

POWER_OFF_DRAG = "power-off-drag-ork-N5800.csv"
POWER_ON_DRAG = "power-on-drag-ork-N5800.csv"

files = [POWER_ON_DRAG, POWER_OFF_DRAG]

# Iterate through each of the drag files, and ensure that all duplicate Mach number lines are removed prior to being passed to the Rocket object
for file in files:
    df = pd.read_csv(file)
    df.columns = ['MACH', 'CD']
    df = df.drop_duplicates(subset='MACH')
    df.to_csv(file, index=False, header=False)


def drogueTrigger(p, h, y):
    # p = pressure
    # h= apogee
    # y = [x, y, z, vx, vy, vz, e0, e1, e2, e3, w1, w2, w3]
    # activate drogue when vz < 0 m/s.
    return True if y[5] < 0 else False


# Define the parachute trigger that deploys the main parachute at a certain height
def mainTrigger(p, h, y):
    # p = pressure
    # h = apogee
    # y = [x, y, z, vx, vy, vz, e0, e1, e2, e3, w1, w2, w3]
    # activate main when vz < 0 m/s and z < 800 m.
    return True if y[5] < 0 and y[2] < 500 else False


# Define a objective funtion

# It considers only launch stage(the journey after reaching apogee).
# The meric paraemter inlcuding distance CG and CP, horizontal distance and apogee.
# the varing parameter is ballast weight(0 to 3000g) and poistion from nose to body cube.
# The inital weight for apogee, stabilityMargin and hozrion distance will be 50, 30 and 20 ( sum of weight will 100)

# initialization
weight_apg = 70
weight_hzor = 10
weight_CPACG = 20

target_appogee = 9144
target_hozri = 0
target_stability = 1

# Define a CG calculator function
DRAG_COEFF_DROGUE = 1.55
DRAG_COEFF_MAIN = 2.20
DROGUE_RADIUS = 0.61
MAIN_RADIUS = 1.22
Nose_cone_w = 0.642  # mass: unit=kg
Nose_cone_l = 38.5 / 100  # distance from nose to current objectve center of mass : unit: m

Body_tub1_w = 1.379
Body_tub1_l = 65.8 / 100

Body_tub2_w = 0.839
Body_tub2_l = 85.8 / 100

Body_tub3_w = 0.191
Body_tub3_l = 90.8 / 100

Body_tub4_w = 0.394
Body_tub4_l = 100.8 / 100

Body_tub5_w = 1.73
Body_tub5_l = 168.3 / 100

Body_tub6_w = 1.49
Body_tub6_l = 221.8 / 100

Fin_w = 0.805
Fin_l = 392.5 / 100

BoatTail_w = 0.1
BoatTail_l = 407.5 / 100

Drog_chute_w = 0.039
Drog_chute_l = 115 / 100

Avonic_w = 2
Avonic_l = 159 / 100

Main_chute_w = 0.12
Main_chute_l = 220 / 100

Motor_w = 5.4
Motor_l = 340 / 100

Ring1_w = 0.28
Ring1_l = 408.5 / 100

Ring2_w = 0.114
Ring2_l = 375.5 / 100

Ring3_w = 0.114
Ring3_l = 351.5 / 100

Ring4_w = 0.114
Ring4_l = 325.5 / 100


def Cen_mass_cal(W, L):
    Z = symbols("Z")  # set center of mass distance as Z

    equation = Nose_cone_w * (Z - Nose_cone_l) + Body_tub1_w * (Z - Body_tub1_l) + Body_tub2_w * (
            Z - Body_tub2_l) + Body_tub3_w * (Z - Body_tub3_l) + Body_tub4_w * (Z - Body_tub4_l) + Body_tub5_w * (
                       Z - Body_tub5_l) + Body_tub6_w * (Z - Body_tub6_l) + Fin_w * (Z - Fin_l) + BoatTail_w * (
                       Z - BoatTail_l) + Drog_chute_w * (Z - Drog_chute_l) + Avonic_w * (
                       Z - Avonic_l) + Main_chute_w * (Z - Main_chute_l) + Motor_w * (Z - Motor_l) + Ring1_w * (
                       Z - Ring1_l) + Ring2_w * (Z - Ring2_l) + Ring3_w * (Z - Ring3_l) + Ring4_w * (
                       Z - Ring4_l) + W * (Z - L)

    sol = solve(equation, Z)
    return float(sol[0])


def Rokcetweigh_len_benchmark(Ballastweigh, position):
    Env = Environment(latitude=32.990254, longitude=-106.974998, elevation=1401, date=(2023, 8, 20, 12))

    N5800 = SolidMotor(
        thrust_source="Cesaroni_20146N5800-P.eng",
        burn_time=3.5,
        grain_number=NUM_GRAINS,
        dry_mass=Empty_Mass,
        dry_inertia=Dry_inertia,
        grain_separation=Distance_each_grain,
        grain_density=(PROPELLANT_MASS / NUM_GRAINS) / (np.pi * (OUTER_GRAIN_RADIUS ** 2) * GRAIN_HEIGHT - np.pi * (
                INNER_GRAIN_RADIUS ** 2) * GRAIN_HEIGHT),
        grain_outer_radius=OUTER_GRAIN_RADIUS,
        grain_initial_inner_radius=INNER_GRAIN_RADIUS,
        grain_initial_height=GRAIN_HEIGHT,
        throat_radius=THROAT_RADIUS,
        nozzle_radius=NOZZLE_RADIUS,
        interpolation_method="linear",
        grains_center_of_mass_position=grains_center_of_mass_position,
        center_of_dry_mass_position=center_of_dry_mass_position,
    )

    MOTOR_MASS = 14.826  # kg

    MOTOR_CASING_MASS = MOTOR_MASS - PROPELLANT_MASS
    LOADED_ROCKET = 20.2 + Ballastweigh  # kg
    UNLOADED_ROCKET = LOADED_ROCKET - MOTOR_MASS + MOTOR_CASING_MASS  # kg
    DIAMETER = 15.4 / 100  # m - note that this is the diameter of the widest point of the rocket body
    RADIUS = DIAMETER / 2
    UNLOADED_CG = Cen_mass_cal(Ballastweigh, position)
    # UNLOADED_CG = 220 / 100
    TOTAL_LENGTH = 398 / 100
    unloaded_i = 0.049 * 0.5 + UNLOADED_ROCKET * (
            1 / 12) * TOTAL_LENGTH * TOTAL_LENGTH * 1.0767  # there is offset paraemter for reducing error from curent function
    UNLOADED_INERTIA_I = unloaded_i  # yaw axis moment after motor burnout

    unloaded_z = 0.049
    UNLOADED_INERTIA_Z = unloaded_z  # roll axis moment after motor burnout

    Inertia = [UNLOADED_INERTIA_Z, UNLOADED_INERTIA_Z, UNLOADED_INERTIA_I, 0, 0, 0]

    MOTOR_CG_DIST_FROM_CONE = TOTAL_LENGTH - MOTOR_LENGTH / 2
    CG_TO_MOTOR_CG = UNLOADED_CG - MOTOR_CG_DIST_FROM_CONE
    CONE_LENGTH = 0.75  # m
    CONE_KIND = 'ogive'
    CONE_DIST_TO_CG = UNLOADED_CG - CONE_LENGTH  # from cone base to CG in m
    DRAG_COEFF_DROGUE = 1.55
    DRAG_COEFF_MAIN = 2.20
    DROGUE_RADIUS = 0.61
    MAIN_RADIUS = 1.22

    Deimos = Rocket(
        inertia=Inertia,
        radius=RADIUS,
        mass=UNLOADED_ROCKET,
        center_of_mass_without_motor=(UNLOADED_CG - CONE_LENGTH),
        power_off_drag=POWER_OFF_DRAG,
        power_on_drag=POWER_ON_DRAG,
        coordinate_system_orientation='nose_to_tail',
    )
    FIN_NUM = 4
    SPAN_HEIGHT = 0.141  # m
    FIN_ROOT_CHORD = 0.344  # m
    FIN_TIP_CHORD = 0.012  # m
    FIN_OFFSET = 0.004  # m
    TOP_OF_FIN = TOTAL_LENGTH - (FIN_ROOT_CHORD + FIN_OFFSET)
    DIST_TO_CG = UNLOADED_CG - TOP_OF_FIN  # m

    Motor_te = Deimos.add_motor(N5800, MOTOR_CG_DIST_FROM_CONE - CONE_LENGTH)
    NoseCone = Deimos.add_nose(length=CONE_LENGTH, kind=CONE_KIND, position=-CONE_LENGTH)
    FinSet = Deimos.add_trapezoidal_fins(
        FIN_NUM, span=SPAN_HEIGHT, root_chord=FIN_ROOT_CHORD, tip_chord=FIN_TIP_CHORD, position=TOP_OF_FIN - CONE_LENGTH
    )
    Tail = Deimos.add_tail(
        top_radius=0.0635, bottom_radius=0.0435, length=0.060, position=TOTAL_LENGTH - CONE_LENGTH - 0.44
    )
    Main = Deimos.add_parachute(
        "Main",
        cd_s=DRAG_COEFF_MAIN * (MAIN_RADIUS) ** 2 * np.pi,  # Times reference area
        trigger=mainTrigger,
        sampling_rate=105,
        lag=1.5,
        noise=(0, 8.3, 0.5),
    )
    Drogue = Deimos.add_parachute(
        "Drogue",
        cd_s=DRAG_COEFF_DROGUE * (DROGUE_RADIUS) ** 2 * np.pi,  # Times reference area
        trigger=drogueTrigger,
        sampling_rate=105,
        lag=1.5,
        noise=(0, 8.3, 0.5),
    )
    TestFlight = Flight(rocket=Deimos, environment=Env, inclination=88, heading=0, rail_length=5.2,
                        terminate_on_apogee=True)
    return TestFlight


position = 1.78
# Number of simualtion
N = 60
# initial matrix
apogee_record = np.zeros(N)
Horizo_record = np.zeros(N)
Stabiliy_record = np.zeros(N)
Ballsat_weight = np.zeros(N)
for i in range(0, N):
    Ballastweigh = 0 + i * 0.04
    Ballsat_weight[i] = Ballastweigh
    results = Rokcetweigh_len_benchmark(Ballastweigh, position)
    apogee_record[i] = results.apogee - target_appogee
    Horizo_record[i] = np.sqrt(results.apogee_x ** 2 + results.apogee_y ** 2)
    Stabiliy_record[i] = results.static_margin.get_source().sum()

apogee_record_pr = apogee_record / (np.max(apogee_record) - np.min(apogee_record))
Horizo_record_pr = Horizo_record / (np.max(Horizo_record) - np.min(Horizo_record))
Stabiliy_record_pr = Stabiliy_record / (np.max(Stabiliy_record) - np.min(Stabiliy_record))
Output_pd = pd.DataFrame({"weight": Ballsat_weight, "apogee": apogee_record_pr, "horizon": Horizo_record_pr,
                          "Stability": Stabiliy_record_pr})
Output_pd.to_excel("Ballst_te1_result1.xlsx")
