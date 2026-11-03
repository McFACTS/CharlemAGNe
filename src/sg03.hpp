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
      double radiative_efficiency = 0.1,
      bool viscosity_switch = false)
      : smbh_mass_(smbh_mass * SOLAR_MASS_KG),
        eddington_ratio_(eddington_ratio),
        alpha_(alpha),
        X_mass_fraction_(X_mass_fraction),
        radiative_efficiency_(radiative_efficiency),
        viscosity_switch_(viscosity_switch),
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
        double log10T = 0.0,
        const SolverOpts& opts = {}
    ) {
    // Construct Parcel
    Parcel out = Parcel(r);

    // Fill things outside of the hot loop
    out.omega = GetAngularVelocity(r);
    out.Teff = GetEffectiveTemperature(r, out.omega);
    out.Mdotprime = GetMDot() * (1. - std::sqrt(GetRmin() / r));
    // Check arguments for rootfind guess
    if (log10cs == 0.0) {
        log10cs = std::log10(out.omega * r * 1.e-2);
    }
    if (log10T == 0.0) {
        log10T = std::log10(out.Teff);
    }

    // Current guess
    double cur[2] = {log10cs, log10T};
    int nfev = 0; // number of evaluations
    // Run the solver
    solve_lmdif<&SG03::cs_T_eqn>(*this, out, cur, opts, nfev);
    // evaluate on the point where it lands
    double fn_pressure, fn_radiation;
    cs_T_eqn(out, cur[0], cur[1], fn_pressure, fn_radiation);
    // Check for errors
    out.converged = std::isfinite(fn_pressure) && std::isfinite(fn_radiation)
      && std::abs(fn_pressure) < opts.conv_tol
      && std::abs(fn_radiation) < opts.conv_tol;
    out.steps = nfev;

    // Fill Q / h
    out.Q = out.omega * out.omega / (2.0 * M_PI * G_SI * out.rho);
    out.h = out.cs / out.omega;
    return out;
    
  }
  /// Solve for Teff and T in a 2D rootfind
  /// We are trying to solve these two equations:
  /// $c_s^2 = \frac{p_{\rm tot}}{T}$
  /// $T^4 = (\frac{3}{8}\tau_V + \frac{1}{2} + \frac{1}{4\tau_V}) T_{\rm eff}^4
  /// by finding the effective temperature and midplane temperature which satisfy them
  ///
  /// For a 2D Newton's method, X_{i + 1} = X_{i} - J^{-1} f(X_{i})
  /// where here f(X_{i}) has rows which are the pressure term and the radiation term
  /// which are both being minimized.
  Parcel solve_Teff_T(
        double r,
        double log10Teff,
        double log10T,
        const SolverOpts& opts = {}
    ) {
    // Construct Parcel
    Parcel out = Parcel(r);

    // Fill things outside of the hot loop
    out.omega = GetAngularVelocity(r);
    out.Teff = GetEffectiveTemperature(r, out.omega);
    out.Mdotprime = GetMDot() * (1. - std::sqrt(GetRmin() / r));
    // The guess is not optional for this one.

    // We know Q and rho by fiat
    out.Q = 1.0;
    out.rho = out.omega * out.omega / (2.0 * M_PI * G_SI);

    // We may or may not know cs by fiat
    if (!viscosity_switch_) {
        out.cs = std::pow(
          out.Mdotprime * out.omega * out.omega / (
            6.0 * M_PI * GetAlpha() * out.rho),
          1.0/3.0
        );
    }

    // current guess
    double cur[2] = {log10Teff, log10T};
    int nfev = 0;
    // Solve
    solve_lmdif<&SG03::Teff_T_eqn>(*this, out, cur, opts, nfev);
    // evaluate on the point where it lands
    double fn_pressure, fn_radiation;
    Teff_T_eqn(out, cur[0], cur[1], fn_pressure, fn_radiation);
    // Check for errors
    out.converged = std::isfinite(fn_pressure) && std::isfinite(fn_radiation)
      && std::abs(fn_pressure) < opts.conv_tol
      && std::abs(fn_radiation) < opts.conv_tol;
    out.steps = nfev;

    // Fill Q / h
    if (out.converged) {
        out.h = out.cs / out.omega;
    }
    return out;
  }

  Profile WalkDisk(
    std::vector<double> r_cyl,
    const double xtol = 1.e-4,
    const int iter_max = 200,
    const int iter_max_first = 1000
  ) {
    // Initialize out
    Profile out = Profile(r_cyl);
    // Initialize things we might reuse
    bool star_formation = false;
    bool switched_solve = false;
    double omega = GetAngularVelocity(r_cyl[0]);
    double Teff = GetEffectiveTemperature(r_cyl[0], omega);
    Parcel prev;
    Parcel cur;
    const SolverOpts opts_first{iter_max_first, xtol};
    const SolverOpts opts_warm{iter_max, xtol};
    // Enter loop
    for (size_t iter=0; iter < r_cyl.size(); iter++) {
      if (iter == 0) {
        cur = solve_cs_T(
            r_cyl[iter],
            std::log10(omega * r_cyl[iter] * 1.e-2),
            std::log10(Teff),
            opts_first
        );
      } else {
        if (!star_formation) {
          cur = solve_cs_T(
              r_cyl[iter],
              std::log10(prev.cs),
              std::log10(prev.T),
              opts_warm
          );
          // 
          if (!(std::isfinite(cur.Q) && cur.Q > 1.0)) {
            star_formation = true;
            switched_solve = true;
          }
        }
        // Not a bug that this is separate check.
        if (star_formation) {
          cur = solve_Teff_T(
            r_cyl[iter],
            std::log10(switched_solve ? cur.Teff : prev.Teff),
            std::log10(switched_solve ? cur.T : prev.T),
            switched_solve ? opts_first : opts_warm
          );
        }
      }
      // Catch non-converged
      if (!cur.converged) {
        // Re-evaluate the residuals to throw only on hard failures
        double fn_p, fn_r;
        if (star_formation) {
          Teff_T_eqn(cur, std::log10(cur.Teff), std::log10(cur.T), fn_p, fn_r);
        } else {
          cs_T_eqn(cur, std::log10(cur.cs), std::log10(cur.T), fn_p, fn_r);
        }
        // Check the maximum residual
        double residual = std::max(std::fabs(fn_p), std::fabs(fn_r));
        fprintf(stderr, "non-converged i=%zu r=%.6gRg sf=%d residual=%.3e\n",
            iter, r_cyl[iter] / GetRg(), (int)star_formation, residual);
        if (!std::isfinite(residual) || residual > 1.0) {
          throw std::runtime_error("Failed to converge disk!\n");
        }
      }
      // We have constructed cur; let's update our profiles
      out.converged[iter]       = cur.converged;
      out.steps[iter]           = cur.steps;
      out.star_formation[iter]  = star_formation;
      out.omega[iter]           = cur.omega;
      out.Teff[iter]            = cur.Teff;
      out.Mdotprime[iter]       = cur.Mdotprime;
      out.Mdot[iter]            = cur.Mdot;
      out.T[iter]               = cur.T;
      out.cs[iter]              = cur.cs;
      out.rho[iter]             = cur.rho;
      out.kappa[iter]           = cur.kappa;
      out.tauV[iter]            = cur.tauV;
      out.PPgas[iter]           = cur.PPgas;
      out.PPrad[iter]           = cur.PPrad;
      out.opacfac[iter]         = cur.opacfac;
      out.Q[iter]               = cur.Q;
      out.h[iter]               = cur.h;
      out.eta[iter]             = cur.eta;

      // Profile specific things
      out.dr_cyl[iter] = (iter == 0) ? 0.0 : cur.r_cyl - out.r_cyl[iter-1];
      out.dA[iter] = out.dr_cyl[iter] * M_PI * 4.0 * cur.r_cyl;
      out.dV[iter] = out.dA[iter] * cur.h;
      out.mass_disk[iter] = out.dV[iter] * cur.rho;
      out.Mdot_sf[iter] = 0.0;

      // Update loop
      prev = cur;
    }
    return out;
  }


private:
  /* Privates of things in constructor */
  double smbh_mass_;
  double eddington_ratio_;
  double alpha_;
  bool viscosity_switch_;
  double X_mass_fraction_;
  double radiative_efficiency_;
  PAGNCombined opac_;

protected:
  void cs_T_eqn(
          Parcel &parcel,
          double log10cs,
          double log10T,
          double &fn_pressure,
          double &fn_radiation
  ) {

        parcel.cs = std::exp(log10cs * M_LN10);
        double cs2 = parcel.cs * parcel.cs;
        parcel.T = std::exp(log10T * M_LN10);
        double Teff4 = parcel.Teff * parcel.Teff * parcel.Teff * parcel.Teff;

        // Density / sound speed have different expression depending on b=0,1
        parcel.rho = (viscosity_switch_) ?
            parcel.Mdotprime * parcel.omega * parcel.omega * AMU / (6 * M_PI * GetAlpha() * parcel.cs * parcel.T * KBOLTZMANN):
            parcel.Mdotprime * parcel.omega * parcel.omega / (6 * M_PI * GetAlpha() * cs2 * parcel.cs);
        //printf("rho: %g\n", rho);

        parcel.kappa = opac_.eval(parcel.rho, parcel.T);
        //printf("kappa: %g\n", kappa);

        parcel.tauV = parcel.kappa * parcel.rho * parcel.cs / parcel.omega;
        //printf("tauV: %g\n", tauV);

        parcel.PPgas = (KBOLTZMANN / AMU) * parcel.T / cs2;
        //printf("PPgas: %g\n", PPgas);
        parcel.PPrad = parcel.tauV * STEFAN_BOLTZMANN / (2 * C_SI) * Teff4 / (parcel.rho * cs2);
        //printf("PPrad: %g\n", PPrad);

        // Rewritten as dimensionless, dividing by cs2
        fn_pressure = parcel.PPgas + parcel.PPrad - 1;
        // Same, dividing by Teff4
        parcel.opacfac = (0.375 * parcel.tauV + 0.5 + 0.25 / parcel.tauV);
        //printf("opacfac: %g\n", opacfac);
        double T4 = parcel.T * parcel.T * parcel.T * parcel.T;
        double T4ratio = Teff4 / T4;
        double Rr = T4ratio * parcel.opacfac;
        fn_radiation = Rr - 1;
        //printf("fn_pressure: %g; fn_rad: %g\n", fn_pressure, fn_radiation);
        //throw std::runtime_error("Crash!\n");
  }

  // Equation for Teff/T rootfind
  void Teff_T_eqn(
          Parcel &parcel,
          double log10Teff,
          double log10T,
          double &fn_pressure,
          double &fn_radiation
  ) {

        parcel.Teff = std::exp(log10Teff * M_LN10);
        parcel.T = std::exp(log10T * M_LN10);
        double Teff4 = parcel.Teff * parcel.Teff * parcel.Teff * parcel.Teff;
        parcel.cs = (viscosity_switch_) ?
            parcel.Mdotprime * parcel.omega * parcel.omega * AMU / (
              6.0 * M_PI * GetAlpha() * parcel.rho * KBOLTZMANN * parcel.T)
              : parcel.cs;
        double cs2 = parcel.cs * parcel.cs;

        parcel.kappa = opac_.eval(parcel.rho, parcel.T);
        //printf("kappa: %g\n", kappa);

        parcel.tauV = parcel.kappa * parcel.rho * parcel.cs / parcel.omega;
        //printf("tauV: %g\n", tauV);

        parcel.PPgas = (KBOLTZMANN / AMU) * parcel.T / cs2;
        //printf("PPgas: %g\n", PPgas);
        parcel.PPrad = parcel.tauV * STEFAN_BOLTZMANN / (2 * C_SI) * Teff4 / (parcel.rho * cs2);
        //printf("PPrad: %g\n", PPrad);

        // Rewritten as dimensionless, dividing by cs2
        fn_pressure = parcel.PPgas + parcel.PPrad - 1;
        // Same, dividing by Teff4
        parcel.opacfac = (0.375 * parcel.tauV + 0.5 + 0.25 / parcel.tauV);
        //printf("opacfac: %g\n", opacfac);
        double T4 = parcel.T * parcel.T * parcel.T * parcel.T;
        double T4ratio = Teff4 / T4;
        double Rr = T4ratio * parcel.opacfac;
        fn_radiation = Rr - 1;
        //printf("fn_pressure: %g; fn_rad: %g\n", fn_pressure, fn_radiation);
        //throw std::runtime_error("Crash!\n");
  }
};
