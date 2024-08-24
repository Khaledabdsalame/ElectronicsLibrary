
Electrical Engineering Calculations
This repository contains a collection of C++ functions for performing various electrical engineering calculations. The code includes functions for calculating voltage, current, resistance, power, and more, suitable for both DC and AC circuits.

Overview
The code provides functionality for:

Voltage calculations
Current calculations
Resistance calculations
Power calculations
AC and DC circuit analysis
Parallel and series resistors
Temperature effects on resistance
Features
Voltage Calculations: Compute voltage differences, momentary voltage, and RMS voltage.
Current Calculations: Determine electric current, momentary current, and RMS current.
Resistance Calculations: Calculate resistance for series and parallel resistors, and effects of temperature on resistance.
Power Calculations: Compute power based on current, resistance, and circuit type.
AC Circuit Analysis: Handle calculations for single-phase and three-phase AC circuits.
Installation
Clone the Repository:

bash
Copier le code
git clone https://github.com/yourusername/electrical-engineering-calculations.git
Navigate to the Project Directory:

bash
Copier le code
cd electrical-engineering-calculations
Compile the Code:

Ensure you have a C++ compiler installed. You can compile the code using g++:

bash
Copier le code
g++ -o electrical_calculations main.cpp
Replace main.cpp with the name of your source file if different.

Usage
Voltage Calculations
float E_Voltage_v1(float φ2, float φ1)
Calculates the voltage difference between two points.

Parameters:
φ2: Electric potential at point 2 (V)
φ1: Electric potential at point 1 (V)
Returns: Voltage difference 
𝑉
=
𝜑
2
−
𝜑
1
V=φ2−φ1 (V)
float E_Voltage_v2(float E, float Q)
Calculates the voltage from energy and charge.

Parameters:
E: Energy in joules (J)
Q: Electric charge in coulombs (C)
Returns: Voltage 
𝑉
=
𝐸
/
𝑄
V=E/Q (V)
Current Calculations
float Electric_Current(float q, float t)
Calculates the electric current from charge and time.

Parameters:
q: Charge (C)
t: Time (s)
Returns: Current 
𝐼
I (A)
float Mumentary_currant(float Ipeak, float w, float temp, float θ)
Calculates the momentary current in an AC circuit.

Parameters:
Ipeak: Peak current (A)
w: Angular frequency (rad/s)
temp: Time (s)
θ: Phase angle (rad)
Returns: Momentary current 
𝑖
(
𝑡
)
i(t) (A)
Resistance Calculations
float Resistance(float ρ, float l, float a)
Calculates the resistance of a conductor.

Parameters:
ρ: Resistivity (Ω·m)
l: Length of the conductor (m)
a: Cross-sectional area (m²)
Returns: Resistance 
𝑅
R (Ω)
float Temperature_eff(float r1, float α, float T1, float T2)
Calculates the resistance of a resistor at a different temperature.

Parameters:
r1: Resistance at temperature 
𝑇
1
T1 (Ω)
α: Temperature coefficient (1/°C)
T1: Initial temperature (°C)
T2: Final temperature (°C)
Returns: Resistance 
𝑅
2
R2 at temperature 
𝑇
2
T2 (Ω)
Power Calculations
float Power(string i, string r, float angle, int circuit, int type)
Calculates power based on current, resistance, and circuit type.

Parameters:
i: Current value with unit (A, C)
r: Resistance value with unit (Ω, J)
angle: Phase angle (rad)
circuit: Identifier for circuit type (DC_circuit, AC_circuit_single_phase, AC_cicuit_3_phase)
type: Identifier for power type (line_to_line, line_to_zero)
Returns: Power (W)
Contributing
Feel free to fork the repository and submit pull requests with improvements or additional features. Please ensure that all code adheres to the existing style and includes appropriate documentation.
Constant used to convert peak voltage/current to RMS (Root Mean Square) value.

R_t: 0xE9864CD9
Identifier for calculating resistance in parallel circuits.

I_1: 0x4478917C
Identifier for calculating current in parallel circuits.

Pi: 3.14
Value of π (pi).

electron_charge: -1.602e-19
Charge of an electron in coulombs (C).

proton_charge: 1.602e-19
Charge of a proton in coulombs (C).

coulomb_constant: 8.988e9
Coulomb's constant in N·m²/C².

line_to_zero: 0x0004
Identifier for line-to-zero voltage in three-phase AC circuits.

Functions
Voltage
float E_Voltage_v1(float φ2, float φ1)
Calculates the voltage difference between two points.

Parameters:
φ2: Electric potential at point 2 (V)
φ1: Electric potential at point 1 (V)
Returns: Voltage difference 
𝑉
=
𝜑
2
−
𝜑
1
V=φ2−φ1 (V)
float E_Voltage_v2(float E, float Q)
Calculates the voltage from energy and charge.

Parameters:
E: Energy in joules (J)
Q: Electric charge in coulombs (C)
Returns: Voltage 
𝑉
=
𝐸
/
𝑄
V=E/Q (V)
Voltage Divider
float Voltage_Div(float vt, float Ri, vector<float> R)
Calculates the voltage drop across a resistor in a series circuit.

Parameters:
vt: Total voltage across the series (V)
Ri: Resistance of the specific resistor (Ω)
R: Vector of resistances in the series (Ω)
Returns: Voltage drop across resistor 
𝑅
𝑖
Ri (V)
Voltage Calculation with Ohm's Law
float Voltage_ohm(float i, float r, float Iz, float Z, int circuit)
Calculates the voltage drop using Ohm's law.

Parameters:
i: Current (A) for DC circuit
r: Resistance (Ω) for DC circuit
Iz: Current (A) for AC circuit
Z: Impedance (Ω) for AC circuit
circuit: Identifier for the type of circuit (DC_circuit or AC_circuit)
Returns: Voltage drop (V)
Momentary Voltage
float momentry_Voltage(float MaximalVoltage, float w, float θ, float t)
Calculates the momentary voltage in an AC circuit.

Parameters:
MaximalVoltage: Peak voltage (V)
w: Angular frequency (rad/s)
θ: Phase angle (rad)
t: Time (s)
Returns: Momentary voltage 
𝑣
(
𝑡
)
v(t) (V)
RMS Voltage
float RMS_Voltage(float Vmax)
Calculates the RMS voltage from peak voltage.

Parameters:
Vmax: Peak voltage (V)
Returns: RMS voltage 
𝑉
𝑟
𝑚
𝑠
V 
rms
​
  (V)
Peak-to-Peak Voltage
float VoltageP_P(float Vmax)
Calculates the peak-to-peak voltage from peak voltage.

Parameters:
Vmax: Peak voltage (V)
Returns: Peak-to-peak voltage 
𝑉
𝑝
−
𝑝
V 
p−p
​
  (V)
Electric Current
float Electric_Current(float q, float t)
Calculates the electric current from charge and time.

Parameters:
q: Charge (C)
t: Time (s)
Returns: Current 
𝐼
I (A)
float Electric_current_Ohm(float v, float R)
Calculates the current using Ohm's law.

Parameters:
v: Voltage (V)
R: Resistance (Ω)
Returns: Current 
𝐼
I (A)
Current in Parallel Circuits
float Current_parallel(vector<float> c_parallel)
Calculates the total current in parallel circuits.

Parameters:
c_parallel: Vector of currents through each parallel load (A)
Returns: Total current 
𝐼
𝑡
𝑜
𝑡
𝑎
𝑙
I 
total
​
  (A)
Current Division
float current_div(float R1, float R2, float R3, float It, float Rt, int type)
Calculates the equivalent resistance or current division in parallel circuits.

Parameters:
R1, R2, R3: Resistances of parallel resistors (Ω)
It: Total current (A)
Rt: Equivalent resistance (Ω)
type: Identifier for resistance (R_t) or current (I_1)
Returns: Equivalent resistance or current (Ω or A)
Alternating Current (AC)
float Alterent_Current(float Vz, float Z)
Calculates the current in an AC circuit using voltage and impedance.

Parameters:
Vz: Voltage drop (V)
Z: Impedance (Ω)
Returns: Current 
𝐼
𝑧
I 
z
​
  (A)
Angular Frequency
float Angular_Frequancy(float frequancy)
Calculates the angular frequency from the frequency.

Parameters:
frequancy: Frequency (Hz)
Returns: Angular frequency 
𝜔
ω (rad/s)
Momentary Current
float Mumentary_currant(float Ipeak, float w, float temp, float θ)
Calculates the momentary current in an AC circuit.

Parameters:
Ipeak: Peak current (A)
w: Angular frequency (rad/s)
temp: Time (s)
θ: Phase angle (rad)
Returns: Momentary current 
𝑖
(
𝑡
)
i(t) (A)
RMS Current
float RMS_current(float Ipeak)
Calculates the RMS current from peak current.

Parameters:
Ipeak: Peak current (A)
Returns: RMS current 
𝐼
𝑟
𝑚
𝑠
I 
rms
​
  (A)
Peak-to-Peak Current
float Ip_p(float Ipeak)
Calculates the peak-to-peak current from peak current.

Parameters:
Ipeak: Peak current (A)
Returns: Peak-to-peak current 
𝐼
𝑝
−
𝑝
I 
p−p
​
  (A)
Resistance
float Resistance(float ρ, float l, float a)
Calculates the resistance of a conductor.

Parameters:
ρ: Resistivity (Ω·m)
l: Length of the conductor (m)
a: Cross-sectional area (m²)
Returns: Resistance 
𝑅
R (Ω)
float Resistance_Ohm(float v, float i)
Calculates the resistance using Ohm's law.

Parameters:
v: Voltage (V)
i: Current (A)
Returns: Resistance 
𝑅
R (Ω)
Temperature Effects of Resistance
float Temperature_eff(float r1, float α, float T1, float T2)
Calculates the resistance of a resistor at a different temperature.

Parameters:
r1: Resistance at temperature 
𝑇
1
T1 (Ω)
α: Temperature coefficient (1/°C)
T1: Initial temperature (°C)
T2: Final temperature (°C)
Returns: Resistance 
𝑅
2
R2 at temperature 
𝑇
2
T2 (Ω)
Resistance in Series
float Resistance_series(vector<float> r_series)
Calculates the total resistance of resistors in series.

Parameters:
r_series: Vector of resistances in series (Ω)
Returns: Total resistance 
𝑅
𝑡
𝑜
𝑡
𝑎
𝑙
R 
total
​
  (Ω)
Resistance in Parallel
float Resistance_parallel(vector<float> r_series)
Calculates the total resistance of resistors in parallel.

Parameters:
r_series: Vector of resistances in parallel (Ω)
Returns: Total resistance 
𝑅
𝑡
𝑜
𝑡
𝑎
𝑙
R 
total
​
  (Ω)
Power
float Power(string i, string r, float angle, int circuit, int type)
Calculates power based on current, resistance, and circuit type.

Parameters:
i: Current value with unit (A, C)
r: Resistance value with unit (Ω, J)
angle: Phase angle (rad)
circuit: Identifier for circuit type (DC_circuit, AC_circuit_single_phase, AC_cicuit_3_phase)
type: Identifier for power type (line_to_line, line_to_zero)
Returns: Power (W)
