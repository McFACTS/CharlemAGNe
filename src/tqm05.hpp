/*
 *  File: tqm05.hpp
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
#include "utils.hpp"
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
   * @param opac_name The name of the opacity table in the sp5 file
   * @param opac_path The path to the sp5 file
   * @param opac_id The ID field of the opacity table in the sp5 file
   * @param smbh_mass Schwarzschild radius
   * @param velocity_dispersion the 'sigma' describing stellar velocity disp
   * @param radiative_efficiency  radiative efficiency of star formation
   * @param transport_efficiency The fraction of the local sound speed which conveys angular momentum transport efficiency.
   * @param supernova_pressure Fraction of pressure from supernovae in rad pressure
   * @param X_mass_fraction Hydrogen abundance
   * @param radius_inner The inner radius (R_g) of the disk
   * @param radius_outer The largest radius for which accretion rate may be constant
   * @param mdot_outer The accretion rate at radius_outer
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
        opac_(opac_name, opac_path, opac_id) {
    // Cache the velocity dispersion
    velocity_dispersion_ = GetVelocityDispersion();
    // Generate quantities for the outer disk radius
    angular_velocity_outer_ = GetAngularVelocity(GetRadiusOuter());
    star_formation_outer_ = std::sqrt((2. * M_SQRT2) * GetMDotOuter() *
      G_SI / (GetTransportEfficiency() * (
        angular_velocity_outer_ * angular_velocity_outer_ * angular_velocity_outer_ *
        GetRadiusOuter() * GetRadiusOuter() * GetRadiusOuter())));
    temperature_outer_ = std::sqrt(
      star_formation_outer_ * GetVelocityDispersion() * GetVelocityDispersion() /
      GetRadiusOuter()
    ) * std::sqrt(std::sqrt(
      (3. * C_SI) / (
        8. * M_SQRT2 * M_PI * G_SI * STEFAN_BOLTZMANN
      )
    ));
    rho_outer_ = angular_velocity_outer_ * angular_velocity_outer_ / (
      M_SQRT2 * M_PI * G_SI);
    aspect_ratio_outer_ = GetRadiusOuter() * star_formation_outer_ / (
      2.0 * M_SQRT2);
    const double MDotMax = 8. * M_PI * rho_outer_ * aspect_ratio_outer_ * 
      GetVelocityDispersion() * GetVelocityDispersion()  * GetRadiusOuter() /
      (C_SI * GetRadiativeEfficiency());
    
    // Check mdot against mdot max
    if (GetMDotOuter() > MDotMax) {
        printf(
            "Cannot set MDotOuter (%g) greater than MDotMax (%g)!\n",
            GetMDotOuter(), MDotMax
        );
        throw std::runtime_error("See stdout message!");
    }
  }

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

  double GetAngularVelocityOuter() const { return angular_velocity_outer_; }
  double GetStarFormationOuter() const { return star_formation_outer_; }
  double GetTemperatureOuter() const { return temperature_outer_; }
  double GetRhoOuter() const { return rho_outer_; }
  double GetAspectRatioOuter() const { return aspect_ratio_outer_; }
  double GetKappaOuter() const {
    return opac_.eval(GetRhoOuter(),GetTemperatureOuter());
  }
  double GetTauVOuter() const {
    return GetKappaOuter() * GetVelocityDispersion() * GetVelocityDispersion()
      * GetStarFormationOuter() / (2.0 * M_PI * G_SI * GetRadiusOuter());
  }
  double GetTeff4Outer() const {
    return (4./3.) * GetTemperatureOuter() * GetTemperatureOuter()
      * GetTemperatureOuter() * GetTemperatureOuter() / (
        GetTauVOuter() + (2. / (3. * GetTauVOuter())) + (4. / 3.)
      );
  }
  double GetSFEfficiencyOuter() const {
    return 2.0 * STEFAN_BOLTZMANN * GetTeff4Outer() / (
      C_SI * C_SI * GetRadiativeEfficiency()
      * 2. * GetRhoOuter() * GetAspectRatioOuter() * GetAngularVelocityOuter()
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
        double mdot_local,
        const SolverOpts& opts = {}
    ) {
    // Construct Parcel
    Parcel out = Parcel(r);

    // Fill things outside of the hot loop
    out.omega = GetAngularVelocity(r);
    out.Mdot = mdot_local;
    out.Mdotprime = out.Mdot * (1. - std::sqrt(GetRadiusInner() / r));
    out.Teff = std::sqrt(
      std::sqrt(
        3.0 * out.Mdotprime * out.omega * out.omega / 
          (8.0 * M_PI * STEFAN_BOLTZMANN)
      )
    );

    // current guess
    double cur[2] = {log10T, log10rho};
    int nfev = 0;
    // Run solver
    solve_lmdif<&TQM05::T_rho_eqn>(*this, out, cur, opts, nfev);

    // Re-evaluate the final point
    double fn_pressure, fn_radiation;
    T_rho_eqn(out, cur[0], cur[1], fn_pressure, fn_radiation);
    // Assign outputs
    if (std::isfinite(fn_pressure) && std::isfinite(fn_radiation)) {
      out.max_residual = std::max(std::fabs(fn_pressure), std::fabs(fn_radiation));
      if (out.max_residual < opts.conv_tol) {
        out.converged = true;
      }
    }
    out.steps = nfev;

    // Fill Q / h
    out.Q = out.omega * out.omega / (M_SQRT2 * M_PI * G_SI * out.rho);
    out.star_formation = !(std::isfinite(out.Q)) || (out.Q <= 1.0);
    return out;
    
  }
  /// Solve for T and eta in a 2D rootfind
  ///
  /// For a 2D Newton's method, X_{i + 1} = X_{i} - J^{-1} f(X_{i})
  /// where here f(X_{i}) has rows which are the pressure term and the radiation term
  /// which are both being minimized.
  Parcel solve_T_eta(
        double r,
        double log10T,
        double log10eta,
        double mdot_local,
        const SolverOpts& opts = {}
    ) {
    // Construct Parcel
    Parcel out = Parcel(r);

    // Fill things outside of the hot loop
    out.star_formation = true;
    out.omega = GetAngularVelocity(r);
    out.Mdot = mdot_local;
    out.Mdotprime = out.Mdot * (1. - std::sqrt(GetRadiusInner() / r));
    out.rho = out.omega * out.omega / (M_SQRT2 * M_PI * G_SI);
    out.h = std::sqrt(
      out.Mdot / (
        4.0 * M_PI * r * out.omega * out.rho * 
          GetTransportEfficiency()
      )
    );
    // Sound speed
    out.cs = out.h * out.omega;

    // current guess
    double cur[2] = {log10T, log10eta};
    int nfev = 0;
    // Run solver
    solve_lmdif<&TQM05::T_eta_eqn>(*this, out, cur, opts, nfev);

    // Re-evaluate final point
    double fn_pressure, fn_radiation;
    T_eta_eqn(out, cur[0], cur[1], fn_pressure, fn_radiation);
    // Assign outputs
    if (std::isfinite(fn_pressure) && std::isfinite(fn_radiation)) {
      out.max_residual = std::max(std::fabs(fn_pressure), std::fabs(fn_radiation));
      if (out.max_residual < opts.conv_tol) {
        out.converged = true;
      }
    }
    out.steps = nfev;

    // Fill Q / h
    out.Q = out.omega * out.omega / (M_SQRT2 * M_PI * G_SI * out.rho);
    return out;
    
  }
  /// Analytic prescription outside Rout
  Parcel solve_outer(double r) {
    // Construct Parcel
    Parcel out = Parcel(r);
    // Angular velocity
    out.omega = GetAngularVelocity(r);
    // We define Q this way here
    out.Q = 1.0;
    // Density
    out.rho = out.omega * out.omega / (M_SQRT2 * M_PI * G_SI);
    // Mass loss
    out.Mdot = GetMDotOuter();
    // star formation parameter
    double fg = std::sqrt((2.0 * M_SQRT2 * out.Mdot * G_SI) /
      (GetTransportEfficiency() * out.omega * out.omega * out.omega * r * r * r));
    // scale height / aspect ratio
    out.h = r * fg  / (2.0 * M_SQRT2);
    // sound speed
    out.cs = GetVelocityDispersion() * fg / 2.0;
    // Temperature
    out.T = std::sqrt(fg * GetVelocityDispersion() * GetVelocityDispersion() / r) *
	  std::sqrt(std::sqrt(3. * C_SI / (8. * M_SQRT2 * M_PI * G_SI * STEFAN_BOLTZMANN)));
    // Opacity
    out.kappa = opac_.eval(out.rho, out.T);
    // Optical depth
    out.tauV = out.kappa * GetVelocityDispersion() * GetVelocityDispersion() * fg /
      (2. * M_PI * G_SI * r);
    // Effective temperature
    out.Teff = std::sqrt(
      std::sqrt(
        (4.0 / 3.0) * (out.T * out.T * out.T * out.T)  / (out.tauV + (2. / (3. * out.tauV)) + (4./3.))
      )
    );
    // temporary product
    double Sigmastar = 2. * STEFAN_BOLTZMANN * (out.Teff * out.Teff * out.Teff * out.Teff) / (C_SI * C_SI * GetRadiativeEfficiency());
    // Star formation efficiency
    out.eta = Sigmastar / (2 * out.rho * out.h * out.omega);
    // Don't crash please
    out.converged = true;
    return out;
  }


  // Construct a 1D TQM disk and solve along it to return profiles
  Profile SolveDisk(
    double rmin = 6.0, // Rg
    double rmax = 0.0, // Rg
    int npts = 10000,
    const SolverOpts& opts = {}
  ) {
    // Convert rmin to SI
    rmin = rmin * GetRg();
    // Check for default rmax
    if (rmax == 0.0) {
      rmax = GetRadiusOuter();
    } else {
      rmax = rmax * GetRg();
    }
    // Check rmax is >= outer radius
    if (rmax < GetRadiusOuter()) {
      throw std::runtime_error("rmax must be >= outer disk radius!\n");
    }
    // Check number of points
    if (npts < 2) {
      throw std::runtime_error("Must have at least two shells!\n");
    }
    // Get r_edges
    // Initialize radii
    const std::vector<double> r_edge = geomspace(
      rmin,
      rmax,
      npts + 1
    );
    // Initialize centers
    std::vector<double> r_centers(npts);
    for (std::size_t i = 0; i + 1 < r_edge.size(); i++) {
       r_centers[i] = 0.5 * (r_edge[i] + r_edge[i + 1]);
    }
    // Initialize output profile
    Profile out = Profile(r_centers);

    // State tracking
    bool star_formation = false;
    bool first_inside = true;
    Parcel prev;
    Parcel cur;

    // Walk from outer to inner radius
    for (int iter = npts - 1; iter >= 0; iter--) {
      const double r = r_centers[iter];
      const double mdot_local = first_inside ? GetMDotOuter() :
        prev.Mdot - out.Mdot_sf[iter + 1];
      // Check for outer disk
      if (r > GetRadiusOuter()) {
        cur = solve_outer(r);
      } else {
        // Always try no-SF solver first
        cur = solve_T_rho(
          r,
          std::log10(first_inside ? GetTemperatureOuter() : prev.T),
          std::log10(first_inside ? GetRhoOuter() : prev.rho),
          mdot_local,
          opts
        );
        // Check for Star formation
        if (!(std::isfinite(cur.Q)) || (cur.Q <= 1.0)) {
          star_formation = true;
          // Try SF solver
          cur = solve_T_eta(
            r,
            std::log10(first_inside ? GetTemperatureOuter() : prev.T),
            std::log10(first_inside ? GetSFEfficiencyOuter() : (prev.eta > 1e-8) ? prev.eta : 1e-8),
            mdot_local,
            opts
          );
        } else {
          // No star formation. No crash.
          star_formation = false;
        }
        // Catch non-converged
        if (!cur.converged) {
          // Check the maximum residual
          fprintf(stderr, "non-converged i=%d r=%.6gRg sf=%d residual=%.3e\n",
            iter, cur.r_cyl / GetRg(), (int)star_formation, cur.max_residual);
          if (!std::isfinite(cur.max_residual) || cur.max_residual > 10.0) {
            throw std::runtime_error("Failed to converge disk!\n");
          }
        }
      }
      // We have constructed cur; let's update our profiles
      out.converged[iter]       = cur.converged;
      out.max_residual[iter]    = cur.max_residual;
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
      out.dr_cyl[iter] = r_edge[iter+1] - r_edge[iter];
      out.dA[iter] = out.dr_cyl[iter] * M_PI * 4.0 * cur.r_cyl;
      out.dV[iter] = out.dA[iter] * cur.h;
      out.mass_disk[iter] = out.dV[iter] * cur.rho;
      out.Mdot_sf[iter] = (r_centers[iter] > GetRadiusOuter()) ? 0.0 : out.mass_disk[iter] * cur.omega * cur.eta;

      // Update loop
      first_inside = false;
      prev = cur;
    }
  return out;
  }

  // Write disk to file
  void WriteDisk(
    const char* fname,
    double rmin = 6.0,
    double rmax = 0.0,
    int npts = 10000,
    const SolverOpts& opts = {}
  ) {
    // Solve the disk
    Profile tqm05 = SolveDisk(rmin, rmax, npts, opts);
    // open a file
    std::FILE* tqm05_out = std::fopen(fname, "w");
    if (!tqm05_out) {
      throw::std::runtime_error(std::string("WriteDisk cannot open '")
        + fname + "': " + std::strerror(errno));
    }
    fprintf(tqm05_out, "#iter, converged, max_residual, steps, star_formation, r_cyl[R_g], omega, Teff, Mdotprime, Mdot, T, cs, rho, kappa, tauV, PPgas, PPrad, opacfac, Q, h, eta, dr, dA, dV, mass_disk, Mdot_sf\n");
    for (size_t i = 0; i < tqm05.r_cyl.size(); i++) {
      fprintf(
        tqm05_out,
        "%lu %d %g %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
        i,
        (int)tqm05.converged[i],
        tqm05.max_residual[i],
        tqm05.steps[i],
        (int)tqm05.star_formation[i],
        tqm05.r_cyl[i] / GetRg(),
        tqm05.omega[i],
        tqm05.Teff[i],
        tqm05.Mdotprime[i],
        tqm05.Mdot[i],
        tqm05.T[i],
        tqm05.cs[i],
        tqm05.rho[i],
        tqm05.kappa[i],
        tqm05.tauV[i],
        tqm05.PPgas[i],
        tqm05.PPrad[i],
        tqm05.opacfac[i],
        tqm05.Q[i],
        tqm05.h[i],
        tqm05.eta[i],
        tqm05.dr_cyl[i],
        tqm05.dA[i],
        tqm05.dV[i],
        tqm05.mass_disk[i],
        tqm05.Mdot_sf[i]
      );
    }
  std::fclose(tqm05_out);
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
  double angular_velocity_outer_;
  double star_formation_outer_;
  double temperature_outer_;
  double rho_outer_;
  double aspect_ratio_outer_;
  PAGNCombined opac_;

protected:
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

        parcel.h = std::sqrt(parcel.Mdot / (
            4.0 * M_PI * parcel.r_cyl * parcel.omega * 
            parcel.rho * GetTransportEfficiency()
          )
        );

        parcel.cs = parcel.h * parcel.omega;

        parcel.tauV = parcel.kappa * parcel.rho * parcel.h;

        parcel.opacfac = 0.75 * (parcel.tauV + (4./3.) + ((2./3.) / parcel.tauV));
        
        // Handle rad cost
        double T4ratio = Teff4 / T4;
        double Rr = T4ratio * parcel.opacfac;
        fn_radiation = Rr - 1;

        parcel.PPgas = (KBOLTZMANN / AMU) * parcel.T / 
            (parcel.cs * parcel.cs);
        parcel.PPrad = parcel.tauV * STEFAN_BOLTZMANN * Teff4 / (
            C_SI * parcel.rho * parcel.cs * parcel.cs
        );

        // Rewritten as dimensionless, dividing by cs2
        fn_pressure = parcel.PPgas + parcel.PPrad - 1;
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
        parcel.eta = std::exp(log10eta * M_LN10);
        double T4 = parcel.T * parcel.T * parcel.T * parcel.T;

        parcel.kappa = opac_.eval(parcel.rho, parcel.T);

        parcel.tauV = parcel.kappa * parcel.rho * parcel.h;

        parcel.opacfac = 0.75 * (parcel.tauV + (4./3.) + ((2./3.) / parcel.tauV));
        
        // Handle rad cost
        double Teff4 = T4 / parcel.opacfac;
        parcel.Teff = std::sqrt(std::sqrt(Teff4));

        // Follow pAGN conventions
        double factorT = parcel.T * KBOLTZMANN / AMU;
        fn_pressure = 1.0 + 
            (2.0 * parcel.h * parcel.eta * parcel.omega * GetRadiativeEfficiency() * C_SI)
            *
            (((parcel.tauV / 2.0) + GetSupernovaPressure()) / factorT) -
            (parcel.cs * parcel.cs / factorT);
        double factor_Teff = STEFAN_BOLTZMANN * Teff4;
        fn_radiation = 1.0 - parcel.rho * parcel.h * parcel.omega
            * parcel.eta * GetRadiativeEfficiency()
            * C_SI * C_SI / factor_Teff -
            3.0 * parcel.omega * parcel.omega * parcel.Mdotprime / (8.  * M_PI * factor_Teff);
  }
};
