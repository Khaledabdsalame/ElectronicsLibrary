#include "Electronics.h"

float E_Voltage_v1(float φ2, float φ1) {
    return φ2 - φ1;
}

float E_Voltage_v2(float E, float Q) {
    return  E / Q;
}

float Voltage_Series(vector<float> v_series) {
    float sum = 0;
    for (const float& i : v_series) {
        sum = sum + i;
    }
    return sum;

}


float Voltage_Div(float vt, float Ri, vector<float>R) {
    float sum_r = 0;

    for (const float& i : R)
    {
        sum_r = sum_r + i;
    }

    float vi = (vt) * (Ri / sum_r);

    return vi;

}


float Voltage_ohm(float i, float r, float Iz, float Z, int circuit) {
    if (circuit == DC_circuit)
    {
        float vr = i * r;
        return vr;
    }

    else if (circuit == AC_circuit)
    {
        float vz = Iz * Z;
        return vz;
    }
    else
    {
        return 0;
    }

}

float momentry_Voltage(float MaximalVoltage, float w, float θ, float t) {
    return  MaximalVoltage * sin((w * t) + θ);

}

float RMS_Voltage(float Vmax) {
    return Vmax * RMS_Const;
}

float VoltageP_P(float Vmax) {
    return  2 * Vmax;
}

float Electric_Current(float q, float t) {
    return  q / t;

};

float Electric_current_Ohm(float v, float R) {
    return  v / R;
}


float Current_parallel(vector <float> c_parallel) {
    float sum = 0;
    for (const float& i : c_parallel) {
        sum = sum + i;
    }
    return sum;
}

float current_div(float R1, float R2, float R3, float It, float Rt, int type) {
    if (type == R_t)
    {
        float rt = 1 / (1 / R2 + 1 / R3);
        return rt;
    }
    else if (type == I_1)
    {
        float I1 = It * Rt / (R1 + Rt);
        return I1;
    }
    else
    {
        return 0;
    }

}

float Alterent_Current(float Vz, float Z) {
    return  Vz / Z;

};

float  Angular_Frequancy(float frequancy) {
    return  2 * frequancy * Pi;

}

float Mumentary_currant(float Ipeak, float w, float temp, float θ) {

    return  Ipeak * sin((w * temp) + θ);
}

float RMS_current(float Ipeak) {

    return Ipeak * RMS_Const;

}

float Ip_p(float Ipeak) {
    return Ipeak * 2;
}

float Resistance(float ρ, float l, float a) {
    return  ρ * (l / a);
}

float Resistance_Ohm(float v, float i) {
    return v / i;
}

float Temperature_eff(float r1, float α, float T1, float T2) {
    return r1 * (1 + α * (T2 - T1));
}

float Resistance_series(vector<float>r_series) {
    float sum = 0;
    for (const float& i : r_series) {
        sum = sum + i;
    }
    return sum;

}

float Resistance_parallel(vector<float>r_series) {
    float sum = 0;
    for (const float& i : r_series) {
        sum = sum + (1 / i);
    }
    return sum;
}

float power(string i, string r, float angle, int circuit, int type) {
    if (circuit == AC_circuit_single_phase)
    {
        string types[5] = { "jou", "sec", "vol", "amp", "ohm" };
        string lti = i.substr(i.length() - 3);
        string ltr = r.substr(r.length() - 3);
        string numberParti = i.substr(0, i.length() - 3);
        float numberi = stof(numberParti);
        string numberPartr = r.substr(0, r.length() - 3);
        float numberr = stof(numberPartr);
        bool iis = false;
        bool ris = false;
        for (int j = 0; j < 5; j++)
        {
            if (lti == types[j])
            {
                iis = true;
            }
        }
        for (int j = 0; j < 5; j++)
        {
            if (ltr == types[j])
            {
                ris = true;
            }
        }
        if (!iis || !ris)
        {
            return -1;
        }


        if (lti == types[0] && ltr == types[1])
        {
            return numberi / numberr;
        }
        else if (lti == types[1] && ltr == types[0])
        {
            return numberr / numberi;
        }
        else if (lti == types[2] && ltr == types[3])
        {
            return numberi * numberr * cos(angle);
        }
        else if (lti == types[3] && ltr == types[2])
        {
            return numberr * numberi * cos(angle);
        }
        else if (lti == types[3] && ltr == types[4])
        {
            return pow(numberi, 2) * numberr * cos(angle);
        }
        else if (lti == types[4] && ltr == types[3])
        {
            return pow(numberr, 2) * numberi * cos(angle);
        }
        else if (lti == types[2] && ltr == types[4])
        {
            return pow(numberi, 2) / numberr * cos(angle);
        }
        else if (lti == types[4] && ltr == types[2])
        {
            return pow(numberr, 2) / numberi * cos(angle);
        }

    }
    else if (circuit == DC_circuit)
    {
        string types[5] = { "jou", "sec", "vol", "amp", "ohm" };
        string lti = i.substr(i.length() - 3);
        string ltr = r.substr(r.length() - 3);
        string numberParti = i.substr(0, i.length() - 3);
        float numberi = stof(numberParti);
        string numberPartr = r.substr(0, r.length() - 3);
        float numberr = stof(numberPartr);
        bool iis = false;
        bool ris = false;
        for (int j = 0; j < 5; j++)
        {
            if (lti == types[j])
            {
                iis = true;
            }
        }
        for (int j = 0; j < 5; j++)
        {
            if (ltr == types[j])
            {
                ris = true;
            }
        }
        if (!iis || !ris)
        {
            return -1;
        }


        if (lti == types[0] && ltr == types[1])
        {
            return numberi / numberr;
        }
        else if (lti == types[1] && ltr == types[0])
        {
            return numberr / numberi;
        }
        else if (lti == types[2] && ltr == types[3])
        {
            return numberi * numberr;
        }
        else if (lti == types[3] && ltr == types[2])
        {
            return numberr * numberi;
        }
        else if (lti == types[3] && ltr == types[4])
        {
            return pow(numberi, 2) * numberr;
        }
        else if (lti == types[4] && ltr == types[3])
        {
            return pow(numberr, 2) * numberi;
        }
        else if (lti == types[2] && ltr == types[4])
        {
            return pow(numberi, 2) / numberr;
        }
        else if (lti == types[4] && ltr == types[2])
        {
            return pow(numberr, 2) / numberi;
        }

    }
    else if (circuit == AC_cicuit_3_phase && type == line_to_line)
    {
        string types[5] = { "jou", "sec", "vol", "amp", "ohm" };
        string lti = i.substr(i.length() - 3);
        string ltr = r.substr(r.length() - 3);
        string numberParti = i.substr(0, i.length() - 3);
        float numberi = stof(numberParti);
        string numberPartr = r.substr(0, r.length() - 3);
        float numberr = stof(numberPartr);
        bool iis = false;
        bool ris = false;
        for (int j = 0; j < 5; j++)
        {
            if (lti == types[j])
            {
                iis = true;
            }
        }
        for (int j = 0; j < 5; j++)
        {
            if (ltr == types[j])
            {
                ris = true;
            }
        }
        if (!iis || !ris)
        {
            return -1;
        }


        if (lti == types[0] && ltr == types[1])
        {
            return numberi / numberr;
        }
        else if (lti == types[1] && ltr == types[0])
        {
            return numberr / numberi;
        }
        else if (lti == types[2] && ltr == types[3])
        {
            return sqrt(3) * numberi * numberr * cos(angle);
        }
        else if (lti == types[3] && ltr == types[2])
        {
            return sqrt(3) * numberr * numberi * cos(angle);
        }
        else if (lti == types[3] && ltr == types[4])
        {
            return 3 * pow(numberi, 2) * numberr * cos(angle);
        }
        else if (lti == types[4] && ltr == types[3])
        {
            return 3 * pow(numberr, 2) * numberi * cos(angle);
        }
        else if (lti == types[2] && ltr == types[4])
        {
            return 3 * (pow(numberi, 2) / numberr) * cos(angle);
        }
        else if (lti == types[4] && ltr == types[2])
        {
            return 3 * (pow(numberr, 2) / numberi) * cos(angle);
        }


    }
    else if (circuit == AC_cicuit_3_phase && type == line_to_zero)
    {
        string types[5] = { "jou", "sec", "vol", "amp", "ohm" };
        string lti = i.substr(i.length() - 3);
        string ltr = r.substr(r.length() - 3);
        string numberParti = i.substr(0, i.length() - 3);
        float numberi = stof(numberParti);
        string numberPartr = r.substr(0, r.length() - 3);
        float numberr = stof(numberPartr);
        bool iis = false;
        bool ris = false;
        for (int j = 0; j < 5; j++)
        {
            if (lti == types[j])
            {
                iis = true;
            }
        }
        for (int j = 0; j < 5; j++)
        {
            if (ltr == types[j])
            {
                ris = true;
            }
        }
        if (!iis || !ris)
        {
            return -1;
        }


        if (lti == types[0] && ltr == types[1])
        {
            return numberi / numberr;
        }
        else if (lti == types[1] && ltr == types[0])
        {
            return numberr / numberi;
        }
        else if (lti == types[2] && ltr == types[3])
        {
            return sqrt(3) * numberi * numberr * cos(angle);
        }
        else if (lti == types[3] && ltr == types[2])
        {
            return sqrt(3) * numberr * numberi * cos(angle);
        }
        else if (lti == types[3] && ltr == types[4])
        {
            return 3 * pow(numberi, 2) * numberr * cos(angle);
        }
        else if (lti == types[4] && ltr == types[3])
        {
            return 3 * pow(numberr, 2) * numberi * cos(angle);
        }
        else if (lti == types[2] && ltr == types[4])
        {
            return 3 * (pow(numberi, 2) / numberr) * cos(angle);
        }
        else if (lti == types[4] && ltr == types[2])
        {
            return 3 * (pow(numberr, 2) / numberi) * cos(angle);
        }
    }

    else
    {
        return -1;

    }

}

float Power_Real(float Vrms, float Irms, float φ) {
    return Vrms * Irms * cos(φ);
}

float Reactive_power(float Vrms, float Irms, float φ) {
    return Vrms * Irms * sin(φ);

}

float Apperent_power(float vrms, float Irms) {
    return vrms * Irms;


}


float relation_real_reactive_apperant(string i, string r) {
    string types[5] = { "wat", "var", "vam" };
    string lti = i.substr(i.length() - 3);
    string ltr = r.substr(r.length() - 3);
    string numberParti = i.substr(0, i.length() - 3);
    float numberi = stof(numberParti);
    string numberPartr = r.substr(0, r.length() - 3);
    float numberr = stof(numberPartr);
    bool iis = false;
    bool ris = false;
    for (int j = 0; j < 5; j++)
    {
        if (lti == types[j])
        {
            iis = true;
        }
    }
    for (int j = 0; j < 5; j++)
    {
        if (ltr == types[j])
        {
            ris = true;
        }
    }
    if (!iis || !ris)
    {
        return -1;
    }


    if (lti == types[0] && ltr == types[1])
    {
        return sqrt(numberi + numberr);
    }
    else if (lti == types[1] && ltr == types[0])
    {
        return sqrt(numberr + numberi);
    }
    else if (lti == types[1] && ltr == types[2])
    {
        return sqrt(numberi - numberr);
    }
    else if (lti == types[2] && ltr == types[1])
    {
        return sqrt(numberr - numberi);
    }
    else if (lti == types[0] && ltr == types[2])
    {
        return  sqrt(numberi - numberr);
    }
    else if (lti == types[2] && ltr == types[0])
    {
        return sqrt(numberr - numberi);
    }

}

float electric_charge(string i, string r) {
    string types[3] = { "sec", "amp", "cou" };

    string lti = i.substr(i.length() - 3);
    string ltr = r.substr(r.length() - 3);

    string numberParti = i.substr(0, i.length() - 3);
    float numberi = stof(numberParti);

    string numberPartr = r.substr(0, r.length() - 3);
    float numberr = stof(numberPartr);

    bool iis = false;
    bool ris = false;

    for (int j = 0; j < 3; j++)
    {
        if (lti == types[j])
        {
            iis = true;
        }
    }
    for (int j = 0; j < 3; j++)
    {
        if (ltr == types[j])
        {
            ris = true;
        }
    }
    if (!iis || !ris)
    {
        return -1;
    }


    if (lti == types[0] && ltr == types[1])
    {
        return numberi / numberr;
    }
    else if (lti == types[0] && ltr == types[2])
    {
        return numberi / numberr;
    }
    else if (lti == types[2] && ltr == types[0])
    {
        return numberr / numberi;
    }
    else if (lti == types[1] && ltr == types[2])
    {
        return numberi * numberr;
    }
    else if (lti == types[2] && ltr == types[1])
    {
        return numberi * numberr;
    }
    else if (lti == types[1] && ltr == types[0])
    {
        return numberr / numberi;
    }
    else {
        return -1;
    }
}

float coulomb_law(float q1, float q2, float r) {
    return (q1 * q2 * coulomb_constant) / pow(r, 2);
}

float efficiency(float in, float out) {
    return 100 * (out / in);
}


float Power_factor(float power_real, float apperant_power) {
    return power_real * abs(apperant_power);

}

namespace units_conv {

    // Voltage conversions
    double VoltsToMillivolts(double volts) { return volts * 1000.0; }
    double VoltsToKilovolts(double volts) { return volts / 1000.0; }
    double VoltsToMicrovolts(double volts) { return volts * 1e6; }
    double MillivoltsToVolts(double millivolts) { return millivolts / 1000.0; }
    double KilovoltsToVolts(double kilovolts) { return kilovolts * 1000.0; }
    double MicrovoltsToVolts(double microvolts) { return microvolts / 1e6; }

    // Current conversions
    double AmperesToMilliamperes(double amperes) { return amperes * 1000.0; }
    double AmperesToKiloamperes(double amperes) { return amperes / 1000.0; }
    double AmperesToMicroamperes(double amperes) { return amperes * 1e6; }
    double MilliamperesToAmperes(double milliamperes) { return milliamperes / 1000.0; }
    double KiloamperesToAmperes(double kiloamperes) { return kiloamperes * 1000.0; }
    double MicroamperesToAmperes(double microamperes) { return microamperes / 1e6; }

    // Resistance conversions
    double OhmsToKiloohms(double ohms) { return ohms / 1000.0; }
    double OhmsToMegaohms(double ohms) { return ohms / 1e6; }
    double OhmsToMilliohms(double ohms) { return ohms * 1000.0; }
    double KiloohmsToOhms(double kiloohms) { return kiloohms * 1000.0; }
    double MegaohmsToOhms(double megaohms) { return megaohms * 1e6; }
    double MilliohmsToOhms(double milliohms) { return milliohms / 1000.0; }

    // Power conversions
    double WattsToMilliwatts(double watts) { return watts * 1000.0; }
    double WattsToKilowatts(double watts) { return watts / 1000.0; }
    double WattsToMicrowatts(double watts) { return watts * 1e6; }
    double MilliwattsToWatts(double milliwatts) { return milliwatts / 1000.0; }
    double KilowattsToWatts(double kilowatts) { return kilowatts * 1000.0; }
    double MicrowattsToWatts(double microwatts) { return microwatts / 1e6; }

    // Capacitance conversions
    double FaradsToMicrofarads(double farads) { return farads * 1e6; }
    double FaradsToNanofarads(double farads) { return farads * 1e9; }
    double FaradsToPicofarads(double farads) { return farads * 1e12; }
    double MicrofaradsToFarads(double microfarads) { return microfarads / 1e6; }
    double NanofaradsToFarads(double nanofarads) { return nanofarads / 1e9; }
    double PicofaradsToFarads(double picofarads) { return picofarads / 1e12; }

    // Inductance conversions
    double HenriesToMillihenries(double henries) { return henries * 1000.0; }
    double HenriesToMicrohenries(double henries) { return henries * 1e6; }
    double MillihenriesToHenries(double millihenries) { return millihenries / 1000.0; }
    double MicrohenriesToHenries(double microhenries) { return microhenries / 1e6; }

    // Conductance conversions
    double SiemensToMillisiemens(double siemens) { return siemens * 1000.0; }
    double MillisiemensToSiemens(double millisiemens) { return millisiemens / 1000.0; }
    // Energy conversions
    double KilowattHoursToWattHours(double kWh) { return kWh * 1000.0; }
    double WattHoursToKilowattHours(double Wh) { return Wh / 1000.0; }
    double KilowattHoursToJoules(double kWh) { return kWh * 3.6e6; }
    double JoulesToKilowattHours(double joules) { return joules / 3.6e6; }
    double WattHoursToJoules(double Wh) { return Wh * 3600.0; }
    double JoulesToWattHours(double joules) { return joules / 3600.0; }
    double JoulesToCalories(double joules) { return joules / 4.184; }
    double CaloriesToJoules(double calories) { return calories * 4.184; }
    double JoulesToKilocalories(double joules) { return joules / 4184.0; }
    double KilocaloriesToJoules(double kilocalories) { return kilocalories * 4184.0; }
    double JoulesToElectronvolts(double joules) { return joules / 1.60218e-19; }
    double ElectronvoltsToJoules(double electronvolts) { return electronvolts * 1.60218e-19; }
}
namespace elctronics_components {
    namespace Capacitor {
        // Function to calculate capacitance (C = Q / V)
        double CalculateCapacitance(double charge, double voltage) {
            return charge / voltage;
        }

        // Function to calculate charge (Q = C * V)
        double CalculateCharge(double capacitance, double voltage) {
            return capacitance * voltage;
        }

        // Function to calculate voltage (V = Q / C)
        double CalculateVoltage(double charge, double capacitance) {
            return charge / capacitance;
        }

        // Function to calculate energy stored in the capacitor (E = 0.5 * C * V^2)
        double CalculateEnergy(double capacitance, double voltage) {
            return 0.5 * capacitance * pow(voltage, 2);
        }

        // Function to calculate equivalent capacitance in series (1/C_eq = 1/C1 + 1/C2 + ...)
        double CalculateSeriesCapacitance(const std::vector<double>& capacitances) {
            double reciprocalSum = 0.0;
            for (double capacitance : capacitances) {
                reciprocalSum += 1.0 / capacitance;
            }
            return 1.0 / reciprocalSum;
        }

        // Function to calculate equivalent capacitance in parallel (C_eq = C1 + C2 + ...)
        double CalculateParallelCapacitance(const std::vector<double>& capacitances) {
            double sum = 0.0;
            for (double capacitance : capacitances) {
                sum += capacitance;
            }
            return sum;
        }
    }
    namespace indactor {
        // Function to calculate inductance (L = V * t / I)
        double CalculateInductance(double voltage, double time, double current) {
            return (voltage * time) / current;
        }

        // Function to calculate voltage across the inductor (V = L * dI/dt)
        double CalculateInductorVoltage(double inductance, double currentChangeRate) {
            return inductance * currentChangeRate;
        }

        // Function to calculate current through the inductor (I = V * t / L)
        double CalculateInductorCurrent(double voltage, double time, double inductance) {
            return (voltage * time) / inductance;
        }

        // Function to calculate energy stored in the inductor (E = 0.5 * L * I^2)
        double CalculateInductorEnergy(double inductance, double current) {
            return 0.5 * inductance * pow(current, 2);
        }

        // Function to calculate equivalent inductance in series (L_eq = L1 + L2 + ...)
        double CalculateSeriesInductance(const std::vector<double>& inductances) {
            double sum = 0.0;
            for (double inductance : inductances) {
                sum += inductance;
            }
            return sum;
        }

        // Function to calculate equivalent inductance in parallel (1/L_eq = 1/L1 + 1/L2 + ...)
        double CalculateParallelInductance(const std::vector<double>& inductances) {
            double reciprocalSum = 0.0;
            for (double inductance : inductances) {
                reciprocalSum += 1.0 / inductance;
            }
            return 1.0 / reciprocalSum;
        }

    }
    namespace Risistor {
        // Function to calculate resistance of a 4-band resistor
        double Calculate4BandResistor(int digit1, int digit2, double multiplier) {
            return (10 * digit1 + digit2) * multiplier;
        }

        // Function to calculate resistance of a 5-band or 6-band resistor
        double Calculate5Or6BandResistor(int digit1, int digit2, int digit3, double multiplier) {
            return (100 * digit1 + 10 * digit2 + digit3) * multiplier;
        }

    }
}


