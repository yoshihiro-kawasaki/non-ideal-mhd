#include "non_ideal_mhd_effect.hpp"


double IonizationRate(double number_density, double temperature, double zeta0, bool attenuation)
{
    double zeta, sigmaH2, sigma0, rho, zetaRA;

    if (attenuation) {
        zetaRA = 7.3e-19;
        // zetaRA = 1.1e-22;
        sigma0 = 96.0;
        rho = 1.4 * cst::PROTON_MASS * number_density;
        sigmaH2 = std::sqrt(cst::BOLTZMANN_CONSTANT*temperature*rho / (M_PI*cst::GRAVITATIONAL_CONSTANT*cst::GAS_MOLECULAR_MASS));
        zeta = zeta0 * std::exp(-sigmaH2/sigma0) + zetaRA;
    } else {
        zeta = zeta0;
    }
    return zeta;
}

double MagneticField(double number_density)
{
    double B;
    B = 1.43e-7 * std::sqrt(number_density);
    // B = 30.0e-6 * std::sqrt(number_density/1.0e4);
    return B;
}

// ########################################################################### //
// class DustState
// ########################################################################### //

void DustState::AllocateArrays()
{
    if (is_allocate_arrays_ == true) {
        std::cout << "DustState is already allocated" << std::endl;
        return;
    }

    if (nbin_ == 0) {
        std::cout << "DustState : nbin = 0" << std::endl;
        return;
    }

    ad_.Resize(nbin_);
    sd_.Resize(nbin_);
    md_.Resize(nbin_);
    nd_.Resize(nbin_);
    Zd_.Resize(nbin_);
    tau_.Resize(nbin_);
    Omega_I_.Resize(nbin_);
    J_I_0_.Resize(nbin_);

    is_allocate_arrays_ = true;
}


void DustState::AllocateArrays(int nbin)
{
    if (is_allocate_arrays_ == true) {
        std::cout << "DustState is already allocated" << std::endl;
        return;
    }

    nbin_ = nbin;

    ad_.Resize(nbin_);
    sd_.Resize(nbin_);
    md_.Resize(nbin_);
    nd_.Resize(nbin_);
    Zd_.Resize(nbin_);
    tau_.Resize(nbin_);
    Omega_I_.Resize(nbin_);
    J_I_0_.Resize(nbin_);

    is_allocate_arrays_ = true;
}


void DustState::DeleteArrays()
{
    if (is_allocate_arrays_ == false) {
        return;
    }

    is_allocate_arrays_ = false;
}

double DustState::GetTotalNeutralDustNumberDensity()
{
    double eps = eps_, eps_inv = 1.0 / eps_;
    double nd0 = 0;
    double Xi, nd0_low, nd0_high, dZ2, coef_dZ2, Zd_high;

    coef_dZ2 = (1.0 - psi_) / (2.0 - psi_);

    for (int i = 0; i < nbin_; ++i) {
        Xi       = Omega_I_[i] * (eps + eps_inv) + 1.0;
        nd0_low  = nd_[i] / Xi;
        Zd_high  = psi_ * tau_[i];
        dZ2      = coef_dZ2 * tau_[i];
        nd0_high = 0.0; //(nd_[i] / std::sqrt(2.0*M_PI*dZ2)) * std::exp(-0.5*Zd_[i]*Zd_[i] / dZ2);
        // nd0_high = (nd_[i] / std::sqrt(2.0*M_PI*dZ2)) * std::exp(-0.5*Zd_high*Zd_high / dZ2);
        nd0     += nd0_low + nd0_high;
    }

    return nd0;
}

double DustState::GetTotalPositivelySingleChargedDustNumberDensity()
{
    double eps = eps_, eps_inv = 1.0 / eps_;
    double ndp1 = 0;
    double Xi, ndp1_low, ndp1_high, dZ2, coef_dZ2, Zd_high;

    coef_dZ2 = (1.0 - psi_) / (2.0 - psi_);

    for (int i = 0; i < nbin_; ++i) {
        Xi        = Omega_I_[i] * (eps + eps_inv) + 1.0;
        ndp1_low  = eps * Omega_I_[i] * nd_[i] / Xi;
        Zd_high   = psi_ * tau_[i];
        dZ2       = coef_dZ2 * tau_[i];
        ndp1_high = 0.0; // (nd_[i] / std::sqrt(2.0*M_PI*dZ2)) * std::exp(- 0.5*(1.0 -Zd_high)*(1.0 - Zd_high) / dZ2);
        ndp1     += ndp1_low + ndp1_high;
    }

    return ndp1;
}

double DustState::GetTotalNegativelySingleChargedDustNumberDensity()
{
    double eps = eps_, eps_inv = 1.0 / eps_;
    double ndm1 = 0;
    double Xi, ndm1_low, ndm1_high, dZ2, coef_dZ2, Zd_high;

    coef_dZ2 = (1.0 - psi_) / (2.0 - psi_);

    for (int i = 0; i < nbin_; ++i) {
        Xi        = Omega_I_[i] * (eps + eps_inv) + 1.0;
        ndm1_low  = eps_inv * Omega_I_[i] * nd_[i] / Xi;
        Zd_high   = psi_ * tau_[i]; 
        dZ2       = coef_dZ2 * tau_[i];
        ndm1_high = 0.0; //(nd_[i] / std::sqrt(2.0*M_PI*dZ2)) * std::exp(- 0.5*(-1.0 -Zd_high)*(-1.0 - Zd_high) / dZ2);
        ndm1     += ndm1_low + ndm1_high;
    }

    return ndm1;
}

// ########################################################################### //
// class CalculationOfIonizationDegree
// ########################################################################### //

CalculationOfIonizationDegree::CalculationOfIonizationDegree()
{

}

CalculationOfIonizationDegree::~CalculationOfIonizationDegree()
{

}

void CalculationOfIonizationDegree::SetParameters(double gas_number_density, double temperature, double ionization_rate, DustState *pdust)
{
    ng_    = gas_number_density;
    T_     = temperature;
    zeta_  = ionization_rate;
    pdust_ = pdust;

    // SetIonizationRate(zeta_);

    // thermal velocity
    ui_ = std::sqrt(8.0 * cst::BOLTZMANN_CONSTANT * T_ / (M_PI * MION));
    ue_ = std::sqrt(8.0 * cst::BOLTZMANN_CONSTANT * T_ / (M_PI * ME));
    // gas-phase recombination rate
    beta_ = CalculateGasPhaseRecombinationRateCoefficient(T_);

    CalculateDustStateParameters();
    
    return;
}

void CalculationOfIonizationDegree::SetGasNumberDensity(double gas_number_density)
{
    ng_ = gas_number_density;
}

void CalculationOfIonizationDegree::SetTemperature(double temperature)
{
    T_ = temperature;
    // thermal velocity
    ui_ = std::sqrt(8.0 * cst::BOLTZMANN_CONSTANT * T_ / (M_PI * MION));
    ue_ = std::sqrt(8.0 * cst::BOLTZMANN_CONSTANT * T_ / (M_PI * ME));
    // gas-phase recombination rate
    beta_ = CalculateGasPhaseRecombinationRateCoefficient(T_);
    // dust
    CalculateDustStateParameters();
}

void CalculationOfIonizationDegree::SetIonizationRate(double ionization_rate)
{
    // ionization rate (see Okuzumi 2009 3.1.2)
    double xH2 = 0.5, xHe = 9.75e-2;
    double zetaH2, zetaHe, zetaH;
    zetaH2 = 0.97 * ionization_rate * xH2;
    zetaHe = 0.84 * ionization_rate * xHe;
    zetaH  = 0.03 * ionization_rate * xH2;
    zeta_  = zetaH2 + zetaHe + zetaH;
    return;
}


void CalculationOfIonizationDegree::SetDustState(DustState *pdust)
{
    pdust_ = pdust;
    CalculateDustStateParameters();
    return;
}


void CalculationOfIonizationDegree::CalculateDustStateParameters()
{
    if (pdust_->IsAllocateArrays() == false) {
        std::cerr << __func__ <<  " : DustState is not allocated" << std::endl;
        return;
    }

    int nbin = pdust_->nbin_;
    double J_I_m1;
    pdust_->nd_total_ = 0.0;
    for (int i = 0; i < nbin; ++i) {
        pdust_->Zd_[i]      = 0.0;
        pdust_->tau_[i]     = pdust_->ad_[i] * cst::BOLTZMANN_CONSTANT * T_ / (cst::CHARGE_UNIT*cst::CHARGE_UNIT);
        pdust_->J_I_0_[i]   = 1.0 + std::sqrt(M_PI / (2.0 * pdust_->tau_[i]));
        J_I_m1              = (1.0 + 1.0 / pdust_->tau_[i])  * (1.0 + std::sqrt(2.0 / (pdust_->tau_[i] + 2.0)));
        pdust_->Omega_I_[i] = pdust_->J_I_0_[i] / J_I_m1;
        pdust_->nd_total_  += pdust_->nd_[i];
    }

    return;
}


double CalculationOfIonizationDegree::CalculatePsi(double eps)
{
    int iter, iter_max;
    double psi, psi_old, err, f, df, tol;

    psi = 0.0;
    iter = 0;
    iter_max = 50;
    tol = 1.0e-7;

    while (1)
    {
        iter++;
        f  = std::exp(psi) - eps * (1.0 - psi);
        df = std::exp(psi) + eps;
        psi_old = psi;
        psi = psi_old - f / df;
        err = std::abs(psi - psi_old);
        if (err < tol) break;
        if (iter == iter_max) {
            std::cerr << __func__  << " warning : iter_max : eps = " << eps << std::endl;
            break;
        }
    }

    return psi;
}

void CalculationOfIonizationDegree::CalculateIonizationDegree(GasChargeState *pgcs)
{
    const int max_iter = 100;
    const double prec  = 1.0e-7;

    int nbin = pdust_->nbin_;
    int iter = 0;
    double eps, eps_old, err_eps, psi, dpsi_deps, eps_inv, Xi, dXi;
    double Z_low, Z_high, Z_mean, dZ_low, dZ_high, dZ, dZ_mean;
    double Ji_low, Ji_high, Je_low, Je_high, Ji, Je, sdi_mean, sde_mean;
    double dJi, dJi_low, dJi_high, dJe, dJe_low, dJe_high;
    double dsdi_mean, dsde_mean;
    double ni, ne, dni, dne, Ai, Bi, Ae, Be;
    double f, df;
    double inv_ndtotal;

    inv_ndtotal = 1.0 / pdust_->nd_total_;

    Ai = ue_ * SE * pdust_->nd_total_ / (2.0 * beta_);
    Bi = 4.0 * beta_ * zeta_ * ng_ / (SI * ui_ * SE * ue_ * pdust_->nd_total_ * pdust_->nd_total_);

    Ae = ui_ * SI * pdust_->nd_total_ / (2.0 * beta_);
    Be = Bi; //4.0 * beta_ * zeta_ * ng_ / (SI * ui_ * SE * ue_ * pdust_->nd_total_ * pdust_->nd_total_);

    // initial guess eps (asume ni = ne)
    eps     = SI * ui_ / (SE * ue_);  // TO2022 eq.(8) 
    eps_old = 0.0;
    
    while (1) 
    {
        iter++;

        // calculate psi, dpsi
        psi       = CalculatePsi(eps);
        dpsi_deps = (1.0 - psi) / (eps + std::exp(psi));

        // calculate mean Z, Ji, Je and dZ_eps, dJi_deps, dJe_deps 
        Z_mean    = 0.0;
        dZ_mean   = 0.0;
        sdi_mean  = 0.0;
        sde_mean  = 0.0;
        dsdi_mean = 0.0;
        dsde_mean = 0.0;
        eps_inv   = 1.0 / eps;

        // ion-dust/electorn-dust collsion cross section for high tau case
        Ji_high  = 1.0 - psi;      // TO2022 eq.(27)
        Je_high  = std::exp(psi);  // TO2022 eq.(28)
        dJi_high = - dpsi_deps;
        dJe_high = std::exp(psi) * dpsi_deps;

        for (int i = 0; i < nbin; ++i) {

            Xi     = pdust_->Omega_I_[i] * (eps + eps_inv) + 1.0;          // TO2022 eq.(16)
            dXi    = pdust_->Omega_I_[i] * (1.0 - eps_inv * eps_inv);

            // dust charge number
            Z_low           = pdust_->Omega_I_[i] * (eps - eps_inv) / Xi;  // TO2022 eq.(17)
            Z_high          = psi * pdust_->tau_[i];                       // TO2022 eq.(24)
            pdust_->Zd_[i]  = Z_low + Z_high;                              // TO2022 eq.(29)
            Z_mean         += pdust_->Zd_[i] * pdust_->nd_[i];

            dZ_low          = - pdust_->Omega_I_[i] * dXi * (eps - eps_inv) / (Xi*Xi) + pdust_->Omega_I_[i] * (1.0 + eps_inv * eps_inv) / Xi;
            dZ_high         = dpsi_deps * pdust_->tau_[i];
            dZ              = dZ_low + dZ_high;
            dZ_mean        += dZ * pdust_->nd_[i];

            // ion-dust/electorn-dust collsion cross section
            Ji_low    = pdust_->J_I_0_[i] * (1.0 + eps_inv) / Xi; // TO2022 eq.(18)
            Ji        = Ji_low + Ji_high;                         // TO2022 eq.(30)
            Je_low    = pdust_->J_I_0_[i] * (1.0 + eps) / Xi;     // TO2022 eq.(19)
            Je        = Je_low + Je_high;                         // TO2022 eq.(31)

            sdi_mean += pdust_->sd_[i] * pdust_->nd_[i] * Ji;
            sde_mean += pdust_->sd_[i] * pdust_->nd_[i] * Je;

            dJi_low   = - pdust_->J_I_0_[i] * dXi * (1.0 + eps_inv) / (Xi*Xi) - pdust_->J_I_0_[i] * eps_inv * eps_inv / Xi;
            dJi       = dJi_low + dJi_high;
            dJe_low   = - pdust_->J_I_0_[i] * dXi * (1.0 + eps) / (Xi*Xi) + pdust_->J_I_0_[i] / Xi;
            dJe       = dJe_low + dJe_high;

            dsdi_mean += pdust_->sd_[i] * pdust_->nd_[i] * dJi;
            dsde_mean += pdust_->sd_[i] * pdust_->nd_[i] * dJe;
        }

        Z_mean    *= inv_ndtotal;
        dZ_mean   *= inv_ndtotal;
        sdi_mean  *= inv_ndtotal;
        sde_mean  *= inv_ndtotal;
        dsdi_mean *= inv_ndtotal;
        dsde_mean *= inv_ndtotal;

        // calculate ni, ne, dni_deps, dne_deps
        ni = Ai * sde_mean * (std::sqrt(1.0 + Bi / (sde_mean * sdi_mean)) - 1.0);  // TO2022 eq.(32)
        ne = Ae * sdi_mean * (std::sqrt(1.0 + Be / (sde_mean * sdi_mean)) - 1.0);  // TO2022 eq.(33)

        dni = ni * dsde_mean / sde_mean 
              - 0.5 * Ai * Bi * (dsdi_mean + (sdi_mean/sde_mean)*dsde_mean) 
                / (dsdi_mean*dsdi_mean * std::sqrt(1.0 + Bi / (sde_mean*sdi_mean)));

        dne = ne * dsdi_mean / sdi_mean 
              - 0.5 * Ae * Be * (dsde_mean + (sde_mean/sdi_mean)*dsdi_mean) 
                / (dsde_mean*dsde_mean * std::sqrt(1.0 + Be / (sde_mean*sdi_mean)));

        // calculate feps, df_deps
        f  = ni  - ne  + Z_mean  * pdust_->nd_total_;    // TO2022 eq.(34)
        df = dni - dne + dZ_mean * pdust_->nd_total_;

        // update eps
        eps_old = eps;
        eps = eps_old - f / df;
        if (std::isnan(eps)) {
            std::cerr << __func__  << " warning : eps is nan : iter = " << iter << std::endl;
            break;
        }

        // convergence check
        err_eps = std::abs(eps - eps_old);
        if (err_eps < prec) break;
        if (iter == max_iter) {
            std::cerr << __func__  << " warning : max_iter = " << max_iter << " : eps = " << eps << std::endl;
            break;
        }

    } // end while

    pgcs->ne_    = ne;
    pgcs->ni_    = ni;
    pdust_->psi_ = psi;
    pdust_->eps_ = eps;

    return;
}



// ############################################################################################################ //
// class NonIdealMHDeffects
// ############################################################################################################ // 


NonIdealMHDeffects::NonIdealMHDeffects()
{
    xH2_           = 0.5;
    xHe_           = 0.0975;
    sigmavH2_coef_ = 2.81e-9 * std::sqrt(POLARIZABILITY_H2);                     // Pinto and Galli 2008, eq.(A.5)
    sigmavHe_coef_ = 2.81e-9 * std::sqrt(POLARIZABILITY_He);                     // Pinto and Galli 2008, eq.(A.5)
    sigmavd_coef_  = 1.3 * std::sqrt(128.0*cst::BOLTZMANN_CONSTANT/(9.0*M_PI));  // Pinto and Galli 2008, eq.(25)
}

NonIdealMHDeffects::~NonIdealMHDeffects()
{

}


void NonIdealMHDeffects::CalculateResistivity(double ng, double T, double B, GasChargeState *pgcs, DustState *pdust)
{
    int nbin;
    double nH, rho_H2, rho_He, sqrtT, logT, sigmavd_coef, cB, cc4pi;
    double sigma_eH2, sigma_eHe, sigma_iH2, sigma_iHe, sigma_dH2, sigma_dHe;
    double tauH2inv, tauHeinv, tau, omega_cyc, reduced_mass, beta, Z, dZ2, ndZ, Xi, eps, eps_inv;
    double sigmaPerp;

    nH           = (cst::MEAN_MOLECULAR_WEIGHT / 1.4) * ng; // convert nH (hydrogen nuclrei number density) to ng (number density of gas particles)
    rho_H2       = M_H2 * xH2_ * nH;
    // rho_H2       = cst::GAS_MOLECULAR_MASS * ng;
    rho_He       = M_He * xHe_ * nH;
    sqrtT        = std::sqrt(T);
    logT         = std::log(T);
    sigmavd_coef = sigmavd_coef_ * sqrtT;
    cB           = cst::SPEED_OF_LIGHT / B;

    /* 
        electron, sigmav from Pinto & Galli (2008) table1
    */
    Z = -1.0;
    // electron - H2
    sigma_eH2 = 1.0e-9 * sqrtT * (0.535 + 0.203*logT - 0.136*logT*logT + 0.050*std::pow(logT, 3.0));
    // sigma_eH2 = 1.0e-15 * std::sqrt(8.0 * cst::BOLTZMANN_CONSTANT * T / (M_PI * cst::ELECTRON_MASS));
    tauH2inv  = sigma_eH2 * rho_H2 / (cst::ELECTRON_MASS + M_H2);
    // electron - He
    sigma_eHe = 1.0e-9 * 0.428 * sqrtT;
    tauHeinv  = sigma_eHe * rho_He / (cst::ELECTRON_MASS + M_He);
    // tau, cyclotron freq
    tau       = 1.0 / (tauH2inv + tauHeinv);
    omega_cyc = Z * cst::CHARGE_UNIT * B /(cst::ELECTRON_MASS * cst::SPEED_OF_LIGHT);
    // hall parameter
    beta      = tau * omega_cyc;
    // conductivity 
    cond_.sigmaO_e_ =   cB * Z * cst::CHARGE_UNIT * pgcs->ne_ * beta;
    cond_.sigmaH_e_ = - cB * Z * cst::CHARGE_UNIT * pgcs->ne_ * beta*beta / (1.0 + beta*beta);
    cond_.sigmaP_e_ =   cB * Z * cst::CHARGE_UNIT * pgcs->ne_ * beta / (1.0 + beta*beta);

    /* 
        ion (Mg+), sigmav from Pinto & Galli (2008) Appendix (A.5)
    */
    Z = 1.0;
    // collision with H2
    reduced_mass = MION * M_H2 / (MION + M_H2);
    sigma_iH2    = sigmavH2_coef_ * std::pow(reduced_mass/cst::PROTON_MASS, -0.5);
    tauH2inv     = sigma_iH2 * rho_H2 / (MION + M_H2);
    // collision with He
    reduced_mass = MION * M_He / (MION + M_He);
    sigma_iHe    = sigmavHe_coef_ * std::pow(reduced_mass/cst::PROTON_MASS, -0.5);
    tauHeinv     = sigma_iHe * rho_He / (MION + M_He);
    // tau (collision(drag) time), cyclotron freq
    tau          = 1.0 / (tauH2inv + tauHeinv);
    omega_cyc    = Z * cst::CHARGE_UNIT * B /(MION * cst::SPEED_OF_LIGHT);
    // hall parameter
    beta         = tau * omega_cyc;
    // conductivity
    cond_.sigmaO_i_ =   cB * Z * cst::CHARGE_UNIT * pgcs->ni_ * beta;
    cond_.sigmaH_i_ = - cB * Z * cst::CHARGE_UNIT * pgcs->ni_ * beta*beta / (1.0 + beta*beta);
    cond_.sigmaP_i_ =   cB * Z * cst::CHARGE_UNIT * pgcs->ni_ * beta / (1.0 + beta*beta);

    /*
        dust 
        The dust conductivities are separately calculated from low-temperarure and high-temperature dust size distribution and
        then summed up.
    */
    nbin    = pdust->nbin_;
    eps     = pdust->eps_;
    eps_inv = 1.0 / eps;
    cond_.sigmaO_d_ = 0.0;
    cond_.sigmaH_d_ = 0.0;
    cond_.sigmaP_d_ = 0.0;
    for (int i = 0; i < nbin; ++i) {

        // colision with H2
        sigma_dH2 = pdust->sd_[i] * sigmavd_coef / std::sqrt(M_H2);
        tauH2inv  = sigma_dH2 * rho_H2 / (pdust->md_[i] + M_H2);
        // collision with He
        sigma_dHe = pdust->sd_[i] * sigmavd_coef / std::sqrt(M_He);
        tauHeinv  = sigma_dHe * rho_He / (pdust->md_[i] + M_He);
        // tau
        tau       = 1.0 / (tauH2inv + tauHeinv);

        omega_cyc = cst::CHARGE_UNIT * B / (pdust->md_[i] * cst::SPEED_OF_LIGHT);
        Xi        = pdust->Omega_I_[i] * (eps + eps_inv) + 1.0;

        /* low-temperature  (Z = -1, +1) */ 
        Z         = -1.0;
        beta      = tau * Z * omega_cyc;
        ndZ       = eps_inv * pdust->Omega_I_[i] * pdust->nd_[i] / Xi;
        cond_.sigmaO_d_ +=   cB * Z * cst::CHARGE_UNIT * ndZ * beta;
        cond_.sigmaH_d_ += - cB * Z * cst::CHARGE_UNIT * ndZ * beta*beta / (1.0 + beta*beta);
        cond_.sigmaP_d_ +=   cB * Z * cst::CHARGE_UNIT * ndZ * beta / (1.0 + beta*beta);

        Z         = 1.0;
        beta      = tau * Z * omega_cyc;
        ndZ       = eps * pdust->Omega_I_[i] * pdust->nd_[i] / Xi;
        cond_.sigmaO_d_ +=   cB * Z * cst::CHARGE_UNIT * ndZ * beta;
        cond_.sigmaH_d_ += - cB * Z * cst::CHARGE_UNIT * ndZ * beta*beta / (1.0 + beta*beta);
        cond_.sigmaP_d_ +=   cB * Z * cst::CHARGE_UNIT * ndZ * beta / (1.0 + beta*beta);

        /* high-temperature (Z = <Z>_high) */ 
        Z         = pdust->psi_ * pdust->tau_[i];
        beta      = tau * Z * omega_cyc;
        dZ2       = (1.0 - pdust->psi_) * pdust->tau_[i] / (2.0 - pdust->psi_);
        ndZ       = pdust->nd_[i] / std::sqrt(2.0 * M_PI * dZ2);
        cond_.sigmaO_d_ +=   cB * Z * cst::CHARGE_UNIT * ndZ * beta;
        cond_.sigmaH_d_ += - cB * Z * cst::CHARGE_UNIT * ndZ * beta*beta / (1.0 + beta*beta);
        cond_.sigmaP_d_ +=   cB * Z * cst::CHARGE_UNIT * ndZ * beta / (1.0 + beta*beta);
    }

    

    /* total conductivity */
    cond_.sigmaO_ = cond_.sigmaO_e_ + cond_.sigmaO_i_ + cond_.sigmaO_d_;
    cond_.sigmaH_ = cond_.sigmaH_e_ + cond_.sigmaH_i_ + cond_.sigmaH_d_;
    cond_.sigmaP_ = cond_.sigmaP_e_ + cond_.sigmaP_i_ + cond_.sigmaP_d_;

    sigmaPerp = std::sqrt(cond_.sigmaH_*cond_.sigmaH_ + cond_.sigmaP_*cond_.sigmaP_);

    /* resistivity */
    cc4pi = cst::SPEED_OF_LIGHT * cst::SPEED_OF_LIGHT / (4.0*M_PI);
    res_.etaO_ = cc4pi / cond_.sigmaO_;
    res_.etaH_ = cc4pi * cond_.sigmaH_ / (sigmaPerp * sigmaPerp);
    res_.etaA_ = cc4pi * cond_.sigmaP_ / (sigmaPerp * sigmaPerp) - res_.etaO_;
}

