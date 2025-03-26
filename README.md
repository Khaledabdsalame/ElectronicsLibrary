# Electronics.h - C++ Library for Electrical Engineering Calculations

## Overview

`Electronics.h` is a C++ header file designed to assist with various electrical engineering calculations. This library provides easy-to-use functions and constants that simplify the computation of key electrical properties such as voltage, current, charge, and more. It is ideal for students, engineers, or hobbyists working with electrical circuits and related systems.

## Features

- **Calculate Electrical Properties**: Compute voltage, current, resistance, power, and other essential electrical quantities.
- **Constants & Macros**: Includes commonly used constants like the speed of light, electron charge, and more.
- **Unit Conversion**: Facilitates conversions between different electrical units (e.g., kilowatt-hours to joules).
- **Capacitor and Inductor Laws**: Contains formulas for various calculations involving capacitors and inductors.
  
## Functions

The library offers a range of functions to help you perform essential calculations:

### Voltage Calculation
- **Function**: `E_Voltage_v1(float φ2, float φ1)`
  - **Description**:Calculates the voltage as the difference between two electric potentials.
  
### Current Calculation
- **Function**: `Electric_Current(float q, float t)`
  - **Description**: Description: Calculates the voltage using energy and electric charge.

### Power Calculation
- **Function**: `power(string i, string r, float angle, int circuit, int type)`
  - **Description**: Calculates the electrical power in different types of circuits, including AC and DC, based on the provided current, resistance, and phase angle..


## Installation

To use `Electronics.h` in your project, simply download the header file and include it in your C++ code:

```cpp
#include "Electronics.h"
```
Ensure that the file is located in a directory accessible to your project, or set up your build system accordingly.
# Example Usage
Here is an example demonstrating how to use the library for a basic calculation:
```cpp
#include <iostream>

int main() {
    vector <float> voltages = {3.0, 4.0, 5.0};
    float total_voltage = Voltage_Series(voltages);
    std::cout << "The voltage is: " << total_voltage << " V" << std::endl;
    
    return 0;
}
```
Documentation
For detailed documentation on how each function works and other utilities, visit the **[official documentation](https://amber-carlynn-12.tiiny.site/)**.
# License 
GNU License Copyright (c) Khaled Abdsalame And Logbibi Mabrouk From [CSA Team](https://github.com/CSA-club)


