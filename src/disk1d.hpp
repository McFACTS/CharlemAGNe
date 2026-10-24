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
//// Parcel ////
/// A parcel of gas
struct Parcel {
    //// Members 
    // The distance from the axis of the disk to the point
    double r_cyl;
    // The angular velocity
    double omega = 0.0;
    // Effective temperature
    double Teff4 = 0.0;
    double Teff = 0.0;
    // Mdotprime
    double Mdotprime = 0.0;
    // Midplane temperature
    double T4 = 0.0;
    double T = 0.0;
    double log10T = 0.0;
    // Sound speed
    double cs = 0.0;
    double log10cs = 0.0;
    // density
    double rho = 0.0;
    double log10rho = 0.0;
    // opacity
    double kappa = 0.0;
    double tauV = 0.0;
    // pressure
    double PPgas = 0.0;
    double PPrad = 0.0;
    // opacfac
    double opacfac = 0.0;
    // Toomre instability
    double Q = 0.0;
    // Scale height
    double h = 0.0;

    //// Constructor
    Parcel(double r_cyl) : r_cyl(r_cyl) {};
};

//// Disks ////
class Disk1D {
public:
  /* Define this in daughter classes */
  Disk1D() {}
};

class SG03 : public Disk1D {
public:
  /* Constructor for SG03
   * Fun fact, Cecilia Payne-Gaposchkin came up with the X/Y/Z
   * mass abundance labels, and she
   * determined that stars are mostly hydrogen and helium.
   * @param smbh_mass Schwarzschild radius
   * @param eddington_ratio Fraction of Eddington Luminosity
   * @param alpha Shakura-Sunyaev viscosity parameter
   * @param X_mass_fraction Schwarzschild radius
   * @param radiative_efficiency The radiative efficiency
  */
  SG03(
      std::string opac_name,
      std::string opac_path,
      int opac_id,
      double smbh_mass = 1.0e8,
      double eddington_ratio = 0.5,
      double alpha = 0.01,
      double X_mass_fraction = 0.7,
      double radiative_efficiency = 0.1)
      : smbh_mass_(smbh_mass * SOLAR_MASS_KG),
        eddington_ratio_(eddington_ratio),
        alpha_(alpha),
        X_mass_fraction_(X_mass_fraction),
        radiative_efficiency_(radiative_efficiency),
        opac_(opac_name, opac_path, opac_id) {}

  const PAGNCombined& opac() const { return opac_; }

  /* const is a promise that it does not modify anything in the object */
  /* accessors will have it. setters will not. */
  double GetSMBHMass() const { return smbh_mass_; }
  double GetSMBHMassSolar() const { return smbh_mass_ / SOLAR_MASS_KG; }
  double GetEddingtonRatio() const { return eddington_ratio_; }
  double GetAlpha() const { return alpha_; }
  double GetX() const { return X_mass_fraction_; }
  double GetRadiativeEfficiency() const { return radiative_efficiency_; }

  /* Mass accretion rate */
  double GetMDot() const {
    /* Electron scattering opacity */
    /* cm^2 g^-1 * (1 m / 100 cm)^2 * (1 kg / 1000 g)^-1 */
    double electron_scattering_opacity = 0.2 * (1.0 + GetX()) * 0.1;
    /* Eddington Luminosity */
    double eddington_luminosity =
        4.0 * M_PI * G_SI * C_SI * GetSMBHMass() /
        electron_scattering_opacity;
    /* Mdot at Eddington */
    double mdot_eddington =
        eddington_luminosity /
        (GetRadiativeEfficiency() * C_SI * C_SI);
    /* Mass accretion rate */
    return mdot_eddington * GetEddingtonRatio();
  }

  /* Rg */
  double GetRg() const {
    return GetSMBHMass() * G_SI / (C_SI * C_SI);
  }

  /* Rmin */
  double GetRmin() const {
    /*  Schwarzschild radius / (4 * rad_eff) 
     *  Schwarzschild radius = 2 Rg */
    return GetRg() / (2.0 * GetRadiativeEfficiency());
  }

  /* These are the things
   *
  */

  /* angular_velocity (omega; from Kepler) */
  double GetAngularVelocity(double r) const {
    return std::sqrt(G_SI * GetSMBHMass() / (r * r * r));
  }

  /* Effective temperature to the fourth power */
  double GetEffectiveTemperature4(double r, double angular_velocity) const {
    return 3.0 * GetMDot() * (1.0 - std::sqrt(GetRmin() / r)) * angular_velocity * angular_velocity / (8.0 * M_PI * STEFAN_BOLTZMANN);
  }
  /* Effective temperature */
  double GetEffectiveTemperature(double r, double angular_velocity) const {
    return std::sqrt(
        std::sqrt(
            GetEffectiveTemperature4(r, angular_velocity)
        )
    );
  }

  /// Solve for cs and T in a 2D rootfind
  /// We are trying to solve these two equations:
  /// $c_s^2 = \frac{p_{\rm tot}}{\rho}$
  /// $T^4 = (\frac{3}{8}\tau_V + \frac{1}{2} + \frac{1}{4\tau_V}) T_{\rm eff}^4
  /// by finding the sound speed and temperature which satisfy them.
  ///
  /// For a 2D Newton's method, X_{i + 1} = X_{i} - J^{-1} f(X_{i})
  /// where here f(X_{i}) has rows which are the pressure term and the radiation term
  /// which are both being minimized.
  Parcel solve_cs_T(
        double r,
        double log10cs = 0.0,
        double log10T = 0.0
    ) {
    // Construct Parcel
    Parcel out = Parcel(r);

    // Fill things outside of the hot loop
    out.omega = GetAngularVelocity(r);
    out.Teff4 = GetEffectiveTemperature4(r, out.omega);
    out.Teff = std::sqrt(std::sqrt(out.Teff4));
    out.Mdotprime = GetMDot() * (1. - std::sqrt(GetRmin() / r));
    // Check arguments for rootfind guess
    if (log10cs == 0.0) {
        log10cs = std::log10(out.omega * r * 1.e-2);
    }
    if (log10T == 0.0) {
        log10T = std::log10(out.Teff);
    }

    // TODO args or properties
    constexpr int iter_max = 30;
    constexpr double tol = 1e-5;
    double fn_pressure, fn_radiation;
    cs_T_eqn(out, log10cs, log10T, fn_pressure, fn_radiation);
    
    double jac_inv[4] = {0.};
    for (int iter=0; iter < iter_max; iter++) {
      // Call the function here so we don't do it twice
      cs_T_eqn(out, log10cs, log10T, fn_pressure, fn_radiation);
      calc_inv_jacobian(out, jac_inv, log10cs, log10T, fn_pressure, fn_radiation);
      double rhs_cs = jac_inv[0] * fn_pressure + jac_inv[1] * fn_radiation;
      double rhs_T  = jac_inv[2] * fn_pressure + jac_inv[3] * fn_radiation;

      // stop. We're all done here.
      if (std::abs(fn_pressure) < tol && std::abs(fn_radiation) < tol)
        break;
      // continue iterating
      log10cs -= rhs_cs;
      log10T -= rhs_T;
    }

    // TODO fill Q / h
    return out;
    
  }


private:
  /* Privates of things in constructor */
  double smbh_mass_;
  double eddington_ratio_;
  double alpha_;
  double X_mass_fraction_;
  double radiative_efficiency_;
  PAGNCombined opac_;

  void calc_inv_jacobian(Parcel &parcel, double jac_inv[4], double log10cs, double log10T, double p, double rad) {

    constexpr double eps = 1e-5;
    double p_mc, rad_mc, p_mT, rad_mT;
    cs_T_eqn(parcel, log10cs - eps, log10T, p_mc, rad_mc);
    cs_T_eqn(parcel, log10cs, log10T - eps, p_mT, rad_mT);

    double jac[4] = {(p - p_mc) / eps, (p - p_mT) / eps,
                     (rad - rad_mc) / eps, (rad - rad_mT) / eps};
    // Analytic inverse
    jac_inv[0] = jac[3];
    jac_inv[1] = -jac[1];
    jac_inv[2] = -jac[2];
    jac_inv[3] = jac[0];
    double inv_det = 1. / (jac[0]*jac[3] - jac[1]*jac[2]);
    for (int i=0; i < 4; i++)
      jac_inv[i] *= inv_det;
  }

  // Equation for cs/T rootfind
  void cs_T_eqn(Parcel &parcel, double log10cs, double log10T, double &fn_pressure, double &fn_radiation) {

        parcel.log10cs = log10cs;
        parcel.log10T = log10T;
        parcel.cs = std::exp(log10cs * M_LN10);
        double cs2 = parcel.cs * parcel.cs;
        parcel.T = std::exp(log10T * M_LN10);

        // Sound speed has a different expression depending on b=0,1
        parcel.rho = (1) ?
            parcel.Mdotprime * parcel.omega * parcel.omega / (6 * M_PI * GetAlpha() * cs2 * parcel.cs) :
            parcel.Mdotprime * parcel.omega * parcel.omega * AMU / (6 * M_PI * GetAlpha() * parcel.cs * parcel.T * KBOLTZMANN);
        parcel.log10rho = std::log10(parcel.rho);
        //printf("rho: %g\n", rho);

        parcel.kappa = opac_.eval(parcel.rho, parcel.T);
        //printf("kappa: %g\n", kappa);

        parcel.tauV = parcel.kappa * parcel.rho * parcel.cs / parcel.omega;
        //printf("tauV: %g\n", tauV);

        parcel.PPgas = (KBOLTZMANN / AMU) * parcel.T / cs2;
        //printf("PPgas: %g\n", PPgas);
        parcel.PPrad = parcel.tauV * STEFAN_BOLTZMANN / (2 * C_SI) * parcel.Teff4 / (parcel.rho * cs2);
        //printf("PPrad: %g\n", PPrad);

        // Rewritten as dimensionless, dividing by cs2
        fn_pressure = parcel.PPgas + parcel.PPrad - 1;
        // Same, dividing by Teff4
        parcel.opacfac = (0.375 * parcel.tauV + 0.5 + 0.25 / parcel.tauV);
        //printf("opacfac: %g\n", opacfac);
        parcel.T4 = parcel.T * parcel.T * parcel.T * parcel.T;
        double T4ratio = parcel.Teff4 / parcel.T4;
        double Rr = T4ratio * parcel.opacfac;
        fn_radiation = Rr - 1;
        //printf("fn_pressure: %g; fn_rad: %g\n", fn_pressure, fn_radiation);
        //throw std::runtime_error("Crash!\n");
  }

};
