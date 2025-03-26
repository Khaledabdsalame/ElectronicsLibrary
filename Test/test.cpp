

#include <iostream>
#include "Electronics.h"
int main() {
    vector <float>voltages = { 3.0, 4.0, 5.0 };
    float total_voltage = Voltage_Series(voltages);
    std::cout << "The voltage is: " << total_voltage << " V" << std::endl;

    return 0;
}