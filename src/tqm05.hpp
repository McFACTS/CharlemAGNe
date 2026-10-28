/*
 *  File: disk1d.hpp
 *
 *  BSD 3-Clause License
 *
 *  Copyright (c) 2026, CharlemAGNe developers
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the copyright holder nor the names of its
 *     contributors may be used to endorse or promote products derived from
 *     this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
#pragma once

//////// Include ////////
//// Local ////
#include "constants.hpp"
#include "opacity.hpp"
#include "disk1d.hpp"
//// Third Party ////
//// Standard Library ////
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <utility>
#include <string>


//////// Setup /////////
using namespace constants;

//////// Objects ////////
class TQM05 : public Disk1D {
public:
  /* Constructor for TQM05
   * @param smbh_mass Schwarzschild radius
   * @param eddington_ratio Fraction of Eddington Luminosity
   * @param alpha Shakura-Sunyaev viscosity parameter
   * @param X_mass_fraction Schwarzschild radius
   * @param radiative_efficiency The radiative efficiency
  */
  TQM05(
      std::string opac_name,
      std::string opac_path,
      int opac_id,
      double smbh_mass = 1.e8,
      double velocity_dispersion = 0.0,
      double radiative_efficiency = 1.e-3,
      double transport_efficiency = 0.2,
      double supernova_pressure = 1.0,
      double X_mass_fraction = 0.7,
      double radius_inner = 6.0,
      double radius_outer = 4.178e+7,
      double mdot_outer = 2.0169e+25
  ): smbh_mass_(smbh_mass * SOLAR_MASS_KG),
        velocity_dispersion_(velocity_dispersion),
        radiative_efficiency_(radiative_efficiency),
        transport_efficiency_(transport_efficiency),
        supernova_pressure_(supernova_pressure),
        X_mass_fraction_(X_mass_fraction),
        radius_outer_(radius_outer),
        radius_inner_(radius_inner),
        mdot_outer_(mdot_outer),
        opac_(opac_name, opac_path, opac_id) {}

  const PAGNCombined& opac() const { return opac_; }

  /* const is a promise that it does not modify anything in the object */
  /* accessors will have it. setters will not. */
  double GetSMBHMass() const { return smbh_mass_; }
  double GetSMBHMassSolar() const { return smbh_mass_ / SOLAR_MASS_KG; }
  double GetVelocityDispersion() const {
    return (velocity_dispersion_ == 0.0) ? 
        2.e5 * std::pow(smbh_mass_ / (1.3e8 * SOLAR_MASS_KG), 1.0/4.24):
        velocity_dispersion_;}
  double GetRadiativeEfficiency() const { return radiative_efficiency_; }
  double GetTransportEfficiency() const { return transport_efficiency_; }
  double GetSupernovaPressure() const { return supernova_pressure_; }
  double GetX() const { return X_mass_fraction_; }
  double GetMDotOuter() const { return mdot_outer_; }

  /* Rg */
  double GetRg() const {
    return GetSMBHMass() * G_SI / (C_SI * C_SI);
  }
  double GetRadiusInner() const { return radius_inner_ * GetRg(); }
  double GetRadiusOuter() const { return radius_outer_ * GetRg(); }

  /* angular_velocity (omega; from Kepler) */
  double GetAngularVelocity(double r) const {
    return std::sqrt(
        (G_SI * GetSMBHMass() / (r * r * r)) +
        (2.0 * GetVelocityDispersion() * GetVelocityDispersion() / (r * r))
    );
  }

  /* Effective temperature to the fourth power */
  double GetEffectiveTemperature4(double r, double angular_velocity) const {
    return 3.0 * GetMDotOuter() * (1.0 - std::sqrt(GetRadiusInner() / r)) * angular_velocity * angular_velocity / (8.0 * M_PI * STEFAN_BOLTZMANN);
  }
  /* Effective temperature */
  double GetEffectiveTemperature(double r, double angular_velocity) const {
    return std::sqrt(
        std::sqrt(
            GetEffectiveTemperature4(r, angular_velocity)
        )
    );
  }

  /// Solve for T and rho in a 2D rootfind
  ///
  /// For a 2D Newton's method, X_{i + 1} = X_{i} - J^{-1} f(X_{i})
  /// where here f(X_{i}) has rows which are the pressure term and the radiation term
  /// which are both being minimized.
  Parcel solve_T_rho(
        double r,
        double log10T,
        double log10rho,
        const int iter_max = 30,
        const double tol = 1.e-5
    ) {
    // Construct Parcel
    Parcel out = Parcel(r);

    // Fill things outside of the hot loop
    out.omega = GetAngularVelocity(r);
    out.Teff = GetEffectiveTemperature(r, out.omega);
    out.Mdotprime = GetMDotOuter() * (1. - std::sqrt(GetRadiusInner() / r));

    double fn_pressure, fn_radiation;
    double jac_inv[4] = {0.};
    for (int iter=0; iter < iter_max; iter++) {
      // Call the function here so we don't do it twice
      calc_inv_jacobian<&TQM05::T_rho_eqn>(
        *this,
        out,
        jac_inv,
        log10T, log10rho,
        fn_pressure, fn_radiation
      );
      double rhs_T = jac_inv[0] * fn_pressure + jac_inv[1] * fn_radiation;
      double rhs_rho  = jac_inv[2] * fn_pressure + jac_inv[3] * fn_radiation;

      // stop. We're all done here.
      if (std::abs(fn_pressure) < tol && std::abs(fn_radiation) < tol) {
        out.converged = true;
        out.steps = iter;
        break;
      }
      // continue iterating
      log10T -= rhs_T;
      log10rho -= rhs_rho;
    }
    // Update steps
    if (!out.converged) {out.steps = iter_max;}

    // Fill Q / h
    //out.Q = out.omega * out.omega / (2.0 * M_PI * G_SI * out.rho);
    return out;
    
  }
private:
  /* Privates of things in constructor */
  double smbh_mass_;
  double velocity_dispersion_;
  double radiative_efficiency_;
  double transport_efficiency_;
  double supernova_pressure_;
  double X_mass_fraction_;
  double radius_inner_;
  double radius_outer_;
  double mdot_outer_;
  PAGNCombined opac_;

  // Equation for T / rho rootfind
  void T_rho_eqn(
          Parcel &parcel,
          double log10T,
          double log10rho,
          double &fn_pressure,
          double &fn_radiation
  ) {

        parcel.T = std::exp(log10T * M_LN10);
        parcel.rho = std::exp(log10rho * M_LN10);
        double Teff4 = parcel.Teff * parcel.Teff * parcel.Teff * parcel.Teff;
        double T4 = parcel.T * parcel.T * parcel.T * parcel.T;

        parcel.kappa = opac_.eval(parcel.rho, parcel.T);
        //printf("kappa: %g\n", kappa);

        parcel.h = std::sqrt(GetMDotOuter() / (
            4.0 * M_PI * parcel.r_cyl * parcel.omega * 
            parcel.rho * GetTransportEfficiency()
          )
        );

        parcel.cs = parcel.h * parcel.omega;

        parcel.tauV = parcel.kappa * parcel.rho * parcel.h;
        //printf("tauV: %g\n", tauV);

        parcel.opacfac = 0.75 * (parcel.tauV + (4./3.) + ((2./3.) / parcel.tauV));
        //printf("opacfac: %g\n", opacfac);
        
        // Handle rad cost
        double T4ratio = Teff4 / T4;
        double Rr = T4ratio * parcel.opacfac;
        fn_radiation = Rr - 1;

        parcel.PPgas = (KBOLTZMANN / AMU) * parcel.T / 
            (parcel.cs * parcel.cs);
        //printf("PPgas: %g\n", PPgas);
        parcel.PPrad = parcel.tauV * STEFAN_BOLTZMANN * Teff4 / (
            C_SI * parcel.rho * parcel.cs * parcel.cs
        );
        //printf("PPrad: %g\n", PPrad);

        // Rewritten as dimensionless, dividing by cs2
        fn_pressure = parcel.PPgas + parcel.PPrad - 1;
        //printf("fn_pressure: %g; fn_rad: %g\n", fn_pressure, fn_radiation);
        //throw std::runtime_error("Crash!\n");
  }

  // Equation for T / eta rootfind
  void T_eta_eqn(
          Parcel &parcel,
          double log10T,
          double log10eta,
          double &fn_pressure,
          double &fn_radiation
  ) {

        parcel.T = std::exp(log10T * M_LN10);
        double eta = std::exp(log10eta * M_LN10);
        double T4 = parcel.T * parcel.T * parcel.T * parcel.T;
        // TODO move out of solve
        parcel.rho = parcel.omega * parcel.omega / (M_SQRT2 * M_PI * G_SI);
        parcel.h = std::sqrt(
            GetMDotOuter() / (
                4.0 * M_PI * parcel.r_cyl * parcel.omega * parcel.rho * 
                    GetTransportEfficiency()
            )
        );

        parcel.kappa = opac_.eval(parcel.rho, parcel.T);
        //printf("kappa: %g\n", kappa);

        // TODO move out of loop
        parcel.cs = parcel.h * parcel.omega;

        parcel.tauV = parcel.kappa * parcel.rho * parcel.h;
        //printf("tauV: %g\n", tauV);

        parcel.opacfac = 0.75 * (parcel.tauV + (4./3.) + ((2./3.) / parcel.tauV));
        //printf("opacfac: %g\n", opacfac);
        
        // Handle rad cost
        double Teff4 = T4 / parcel.opacfac;
        parcel.Teff = std::sqrt(std::sqrt(Teff4));

        // Follow pAGN conventions
        double factorT = parcel.T * KBOLTZMANN / AMU;
        fn_pressure = 1.0 + 
            (2.0 * parcel.h * eta * parcel.omega * GetRadiativeEfficiency() * C_SI)
            *
            (((parcel.tauV / 2.0) + GetSupernovaPressure()) / factorT) -
            (parcel.cs * parcel.cs / factorT);
        double factor_Teff = STEFAN_BOLTZMANN * Teff4;
        fn_radiation = 1.0 - parcel.rho * parcel.h * parcel.omega * eta * GetRadiativeEfficiency()
            * C_SI * C_SI / factor_Teff -
            3.0 * parcel.omega * parcel.omega * parcel.Mdotprime / (8.  * M_PI * factor_Teff);
        //printf("fn_radiation: %g; fn_rad: %g\n", fn_pressure, fn_radiation);
        //printf("fn_pressure: %g; fn_rad: %g\n", fn_pressure, fn_radiation);
        //throw std::runtime_error("Crash!\n");

        // TODO set eta?
  }
};
