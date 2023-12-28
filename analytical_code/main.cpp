#include "array.hpp"
#include "constants.hpp"
#include "non_ideal_mhd_effect.hpp"

#include <iostream>
#include <fstream>
#include <string>

// g++ -std=c++17 chemi_and_res.cpp main.cpp

void SetDustSingleSize(DustState &dust, double ad, double fdg, double rhog, double rhodi);

void SetDustMRNSizeDistribution(DustState &dust, int nbin, double admin, double admax, double q, double fdg, double rhog, double rhodi);

void Test1();

void Test2();

void Test3();

int main()
{
    Test1();
    // Test2();
    // Test3();

    return 0;
}

void Test1()
{
    int nbin;
    double ng, rhog, T, zeta, B, ad, admin, admax, q, fdg, rhodi;

    ng    = 1.0e4;
    rhog  = 2.34 * cst::PROTON_MASS * ng;
    T     = 10.0;
    zeta  = 1.0e-17;
    ad    = 1.0e-5;
    admin = 5.0e-7;
    admax = 2.5e-5;
    q     = -3.5;
    nbin  = 30;
    fdg   = 1.0e-2;
    rhodi = 1.4;
    B     = 2.0e-7 * std::sqrt(ng);

    DustState dust;
    // SetDustSingleSize(dust, ad, fdg, rhog, rhodi);
    SetDustMRNSizeDistribution(dust, nbin, admin, admax, q, fdg, rhog, rhodi);

    std::cout << std::scientific;
    // for (int i = 0; i < dust.nbin_; ++i) {
    //     std::cout << std::setw(5) << i << " "
    //               << dust.ad_[i] << " "
    //               << dust.sd_[i] << " "
    //               << dust.md_[i] << " "
    //               << dust.nd_[i]
    //               << std::endl;
    // }

    // return;

    GasChargeState gcs;

    CalculationOfIonizationDegree cid;
    cid.SetParameters(ng, T, zeta, &dust);
    cid.CalculateIonizationDegree(&gcs);


    double nd0, ndp1, ndm1;
    nd0  = dust.GetTotalNeutralDustNumberDensity();
    ndp1 = dust.GetTotalPositivelySingleChargedDustNumberDensity();
    ndm1 = dust.GetTotalNegativelySingleChargedDustNumberDensity();

    NonIdealMHDeffects nih;

    nih.CalculateResistivity(ng, T, B, &gcs, &dust);

    std::cout << std::scientific;
    std::cout << "ni         = " << gcs.ni_ << std::endl;
    std::cout << "ne         = " << gcs.ne_ << std::endl;
    std::cout << "nd(Z = 0)  = " << nd0  << std::endl;
    std::cout << "nd(Z = 1)  = " << ndp1 << std::endl;
    std::cout << "nd(Z = -1) = " << ndm1 << std::endl;

    std::cout << std::endl;
    std::cout << "xi         = " << (gcs.ni_/ng) << std::endl;
    std::cout << "xe         = " << (gcs.ne_/ng) << std::endl;
    std::cout << "xd(Z = 0)  = " << (nd0/ng)  << std::endl;
    std::cout << "xd(Z = 1)  = " << (ndp1/ng) << std::endl;
    std::cout << "xd(Z = -1) = " << (ndm1/ng) << std::endl;

    std::cout << std::endl;
    std::cout << nih.res_.etaO_ << std::endl;
    std::cout << nih.res_.etaH_ << std::endl;
    std::cout << nih.res_.etaA_ << std::endl;

    return;
}


void Test2()
{
    int N = 150;
    int nbin = 30;
    double ng, rhog, T, zeta, B, ad, admin, admax, q, fdg, rhodi;
    double ng1, ng2, dlnng, ngc, gammaT;
    double nd0, ndp1, ndm1;
    std::string filename ="test2_MRN.txt";

    ng1    = 1.0e4;
    ng2    = 1.0e15;
    dlnng  = (std::log10(ng2) - std::log10(ng1)) / (double(N-1));
    gammaT = 7.0 / 5.0;
    ngc    = 2.0e10;

    zeta  = 1.0e-17;
    ad    = 1.0e-5;
    admin = 5.0e-7;
    // admin = 1.0e-5;
    admax = 2.5e-5;
    q     = -3.5;
    fdg   = 1.0e-2;
    rhodi = 2.0;

    DustState dust;
    dust.AllocateArrays(nbin);
    GasChargeState gcs;
    CalculationOfIonizationDegree cid;
    NonIdealMHDeffects nih;

    std::ofstream file(filename, std::ios::trunc | std::ios::out);
    if (!file.is_open()) {
        std::cout << "#Error : file is not open " << std::endl;
        return;
    }

    file << std::scientific;
    for (int i = 0; i < N; ++i){
        ng   = ng1 * std::pow(10.0, dlnng*double(i));
        rhog = 2.34 * cst::PROTON_MASS * ng;
        T    = 10.0 * (1.0 + gammaT * std::pow(ng / ngc, gammaT - 1.0));
        B    = 2.0e-7 * std::sqrt(ng);
        // SetDustSingleSize(dust, ad, fdg, rhog, rhodi);
        SetDustMRNSizeDistribution(dust, nbin, admin, admax, q, fdg, rhog, rhodi);
        cid.SetParameters(ng, T, zeta, &dust);
        cid.CalculateIonizationDegree(&gcs);
        nd0  = dust.GetTotalNeutralDustNumberDensity();
        ndp1 = dust.GetTotalPositivelySingleChargedDustNumberDensity();
        ndm1 = dust.GetTotalNegativelySingleChargedDustNumberDensity();
        nih.CalculateResistivity(ng, T, B, &gcs, &dust);
        // std::cout << std::scientific << ng << " " << T << " " << (gcs.ne_/ng) << " " << (gcs.ni_/ng) << std::endl;
        file << ng << " " << T << " " << (gcs.ne_/ng) << " " << (gcs.ni_/ng) << " "
             << (nd0/ng) << " " << (ndp1/ng) << " " << (ndm1/ng) << " "
             << nih.cond_.sigmaO_e_ << " " << nih.cond_.sigmaH_e_ << " " << nih.cond_.sigmaP_e_ << " "
             << nih.cond_.sigmaO_i_ << " " << nih.cond_.sigmaH_i_ << " " << nih.cond_.sigmaP_i_ << " "
             << nih.cond_.sigmaO_d_ << " " << nih.cond_.sigmaH_d_ << " " << nih.cond_.sigmaP_d_ << " "
             << nih.cond_.sigmaO_   << " " << nih.cond_.sigmaH_   << " " << nih.cond_.sigmaP_   << " "
             << nih.res_.etaO_ << " " << nih.res_.etaH_ << " " << nih.res_.etaA_ << " "
             << std::endl;
    }

    file.close();
}


void Test3()
{
    int Nz = 200;
    double zHmax = 5.0;
    double dzH = zHmax / double(Nz - 1.0);
    double zH;

    double Ms    = 0.5 * cst::SOLAR_MASS; 
    double rAU   = 50.0;
    double rcm   = 1.0 * cst::ASTRONOMICAL_UNIT;
    double Sigma = 1.7e3 * std::pow(rAU, -1.5);
    double T     = 2.8e2 * std::pow(rAU, -0.5);
    // double T     = 1.5e2 * std::pow(rAU, -3.0/7.0);
    double cs    = std::sqrt(cst::BOLTZMANN_CONSTANT * T / cst::GAS_MOLECULAR_MASS);
    double Omega = std::sqrt(cst::GRAVITATIONAL_CONSTANT*Ms / (rcm * rcm * rcm));
    double Hg    = cs / Omega;
    double rhog0 = Sigma / (std::sqrt(2.0*M_PI) * Hg);
    double Bz    = 1.0e-1;
    double ng, rhog, zeta, nd0, ndp1, ndm1, vA, elsO, elsH, elsA;

    double ad = 1.0e-5;
    double admin = 1.0e-5;
    // double admax = 2.5e-5;
    double admax = 1.0e-2;
    double q     = -3.5;
    double fdg   = 1.0e-1;
    double rhodi = 2.0;
    int    nbin  = 150;

    std::string filename ="test3.txt";
    std::ofstream file(filename, std::ios::trunc | std::ios::out);
    if (!file.is_open()) {
        std::cout << "#Error : file is not open " << std::endl;
        return;
    }

    DustState dust;
    dust.AllocateArrays(nbin);
    GasChargeState gcs;
    CalculationOfIonizationDegree cid;
    NonIdealMHDeffects nih;

    file << std::scientific;
    for (int i = 0; i < Nz; ++i) {
        zH   = dzH * i;
        rhog = rhog0 * std::exp(- 0.5 * zH * zH);
        ng   = rhog / cst::GAS_MOLECULAR_MASS;
        zeta = 1.0e-17;
        SetDustMRNSizeDistribution(dust, nbin, admin, admax, q, fdg, rhog, rhodi);
        cid.SetParameters(ng, T, zeta, &dust);
        cid.CalculateIonizationDegree(&gcs);
        nd0  = dust.GetTotalNeutralDustNumberDensity();
        ndp1 = dust.GetTotalPositivelySingleChargedDustNumberDensity();
        ndm1 = dust.GetTotalNegativelySingleChargedDustNumberDensity();
        nih.CalculateResistivity(ng, T, Bz, &gcs, &dust);
        vA = Bz / std::sqrt(4.0*M_PI*rhog);
        elsO = vA*vA / (Omega * nih.res_.etaO_);
        elsH = vA*vA / (Omega * nih.res_.etaH_);
        elsA = vA*vA / (Omega * nih.res_.etaA_);
        file << zH << " " << rhog << " " << (gcs.ne_/ng) << " " << (gcs.ni_/ng) << " "
             << (nd0/ng) << " " << (ndp1/ng) << " " << (ndm1/ng) << " "
             << nih.cond_.sigmaO_e_ << " " << nih.cond_.sigmaH_e_ << " " << nih.cond_.sigmaP_e_ << " "
             << nih.cond_.sigmaO_i_ << " " << nih.cond_.sigmaH_i_ << " " << nih.cond_.sigmaP_i_ << " "
             << nih.cond_.sigmaO_d_ << " " << nih.cond_.sigmaH_d_ << " " << nih.cond_.sigmaP_d_ << " "
             << nih.cond_.sigmaO_   << " " << nih.cond_.sigmaH_   << " " << nih.cond_.sigmaP_   << " "
             << nih.res_.etaO_ << " " << nih.res_.etaH_ << " " << nih.res_.etaA_ << " "
             << elsO <<  " " << elsH << " " << elsA
             << std::endl;
    }

    file.close();
}


void SetDustSingleSize(DustState &dust, double ad, double fdg, double rhog, double rhodi)
{
    dust.nbin_ = 1;
    if (dust.IsAllocateArrays() == false) {
        dust.AllocateArrays();
    }
    dust.ad_[0] = ad;
    dust.sd_[0] = M_PI * dust.ad_[0] * dust.ad_[0];
    dust.md_[0] = (4.0*M_PI/3.0) * dust.ad_[0] * dust.ad_[0] * dust.ad_[0] * rhodi;
    dust.nd_[0] = fdg * rhog / dust.md_[0];
}

void SetDustMRNSizeDistribution(DustState &dust, int nbin, double admin, double admax, double q, double fdg, double rhog, double rhodi)
{
    dust.nbin_ = nbin;
    if (dust.IsAllocateArrays() == false) {
        dust.AllocateArrays(nbin);
    }

    double C1 = 3.0 * (4.0 + q) * fdg * rhog / (4.0 * M_PI * rhodi * (std::pow(admax, 4.0 + q) - std::pow(admin, 4.0 + q)));
    double C2;
    double dlna = (std::log10(admax) - std::log10(admin)) / double(nbin - 1);
    double admin_i, admax_i;

    for (int i = 0; i < nbin; ++i) {

        admin_i = admin * std::pow(10.0, dlna*double(i));
        admax_i = admin * std::pow(10.0, dlna*double(i+1));

        C2 = (1.0 + q) / (std::pow(admax_i, 1.0 + q) - std::pow(admin_i, 1.0 + q));

        dust.ad_[i] = C2 * (std::pow(admax_i, 2.0 + q) - std::pow(admin_i, 2.0 + q)) / (2.0 + q);
        dust.sd_[i] = M_PI * C2 * (std::pow(admax_i, 3.0 + q) - std::pow(admin_i, 3.0 + q)) / (3.0 + q);
        dust.md_[i] = (4.0 * M_PI / 3.0) * C2 * (std::pow(admax_i, 4.0 + q) - std::pow(admin_i, 4.0 + q)) / (4.0 + q);
        dust.nd_[i] = C1 * (std::pow(admax_i, 1.0 + q) - std::pow(admin_i, 1.0 + q)) / (1.0 + q);

    }
}