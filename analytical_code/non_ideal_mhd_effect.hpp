#ifndef NON_IDEAL_MHD_EFFECT_HPP_
#define NON_IDEAL_MHD_EFFECT_HPP_

#include "constants.hpp"
#include "array.hpp"

#include <iostream>
#include <cmath>

/* constants */

constexpr double POLARIZABILITY_H  = 0.667;
constexpr double POLARIZABILITY_H2 = 0.804;
constexpr double POLARIZABILITY_He = 0.207;

constexpr double ME   = cst::ELECTRON_MASS;      // electron mass
constexpr double MUI  = 29.0;                    // ion mean molecular weight
constexpr double MION = MUI * cst::PROTON_MASS;  // ion mass
constexpr double SI   = 1.0;                     // sticking probablity of ions
constexpr double SE   = 0.6;                     // sticking probablity of electorns

constexpr double M_H2 = 2.0 * cst::PROTON_MASS;  // H2 mass
constexpr double M_He = 4.0 * cst::PROTON_MASS;  // He mass

/* prototype */
struct GasChargeState;
class  DustState;
class  CalculationOfIonizationDegree;
struct Conductivity;
struct Resistivity;
class  NonIdealMHDeffects;



/*
   Gas Charge State
*/
struct GasChargeState
{
    double ne_; // electron number density
    double ni_; // ion number density
};

/*
   Dust State
*/
class DustState
{
private:

    bool is_allocate_arrays_;

public:

    DustState() { is_allocate_arrays_ = false; }
    ~DustState() { DeleteArrays(); }

    bool IsAllocateArrays() { return is_allocate_arrays_; }
    void AllocateArrays();
    void AllocateArrays(int nbin);
    void DeleteArrays();

    double GetTotalNeutralDustNumberDensity();
    double GetTotalPositivelySingleChargedDustNumberDensity();
    double GetTotalNegativelySingleChargedDustNumberDensity();


    int nbin_;                     // number of bins
    array::DoubleArray1D ad_;      // dust size
    array::DoubleArray1D sd_;      // dust cross section
    array::DoubleArray1D md_;      // dust mass
    array::DoubleArray1D nd_;      // dust number density
    array::DoubleArray1D Zd_;      // dust charge number
    array::DoubleArray1D tau_;     // the normalized temperature
    array::DoubleArray1D Omega_I_; // TO2022 eq.(15)
    array::DoubleArray1D J_I_0_;   // J(tau(I), nu = 0)
    double nd_total_;              // total dust number density
    double psi_;                   // the parameter which determines the dust charge distribution TO2022 eq.(23) - (25)
    double eps_;                   // the ratio of ions to elecrons density
};

/*
    CalculationOfIonizationDegree
*/
class CalculationOfIonizationDegree
/*
    Tsukamoto and Okuzumi 2022 (TO2022)
*/
{
public:

    CalculationOfIonizationDegree();
    ~CalculationOfIonizationDegree();

    void SetParameters(double number_density, double temperature, double ionization_rate, DustState *pdust);
    void SetGasNumberDensity(double number_density);
    void SetTemperature(double temperature);
    void SetIonizationRate(double ionization_rate);
    void SetDustState(DustState *pdust);

    void CalculateIonizationDegree(GasChargeState *pgcs);
    
private:

    void CalculateDustStateParameters();
    double CalculateGasPhaseRecombinationRateCoefficient(double T) {
        return 2.4e-7 * std::pow(T_/300.0, -0.69);
    };
    double CalculatePsi(double eps);

    double ng_;         // gas number density [cm^-3]
    double T_;          // temperature [K]
    double zeta_;       // ionization rate [s^-1]
    DustState *pdust_;  // pointer of dust state

    double ui_;         // thermal velocity of ion
    double ue_;         // thermal velocity of electron
    double beta_;       // recombination rate between ion and electron
};

/* Conductivity */
struct Conductivity
{
    double sigmaO_;    // Ohmic conductivity
    double sigmaH_;    // Hall conductivity
    double sigmaP_;    // Pedersen conductivity
    double sigmaO_i_;  // Ohmic conductivity of ion
    double sigmaH_i_;  // Hall conductivity of ion
    double sigmaP_i_;  // Pedersen conductivity pf ion
    double sigmaO_e_;  // Ohmic conductivity of electron
    double sigmaH_e_;  // Hall conductivity of electron
    double sigmaP_e_;  // Pedersen conducvtivity of electron
    double sigmaO_d_;  // Ohmic conductivity of dust
    double sigmaH_d_;  // Hall conductivity of dust
    double sigmaP_d_;  // Pedersen conductivity of dust
};

/* Resistivity */
struct Resistivity
{
    double etaO_;  // Ohmic resistivity
    double etaH_;  // Hall resistivity
    double etaA_;  // Ambipoar diffusion resistivity
};

/* NonIdealMHDeffects  */
class  NonIdealMHDeffects
{
public:

    NonIdealMHDeffects();
    ~NonIdealMHDeffects();

    Conductivity cond_;
    Resistivity  res_;

    // void CalculateConductivity(double nH, double T, double B, GasChargeState *pgcs, DustState *pdust);
    void CalculateResistivity(double nH, double T, double B, GasChargeState *pgcs, DustState *pdust);

private:

    double xH2_;            // abundance of H2 
    double xHe_;            // abundance of He

    double sigmavH2_coef_;
    double sigmavHe_coef_;
    double sigmavd_coef_;
};


#endif /* NON_IDEAL_MHD_EFFECT_HPP_ */