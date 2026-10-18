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
//// Third Party ////
//// Standard Library ////
#include <cmath>
#include <cstdio>


//////// Setup /////////
using namespace constants;

//////// Objects /////////
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
      double smbh_mass = 1.0e8,
      double eddington_ratio = 0.5,
      double alpha = 0.01,
      double X_mass_fraction = 0.7,
      double radiative_efficiency = 0.1)
      : smbh_mass_(smbh_mass * SOLAR_MASS_KG),
        eddington_ratio_(eddington_ratio),
        alpha_(alpha),
        X_mass_fraction_(X_mass_fraction),
        radiative_efficiency_(radiative_efficiency) {}

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
    return 3.0 * GetSMBHMass() * (1 - std::sqrt(GetRmin() / r)) * angular_velocity * angular_velocity / (8.0 * M_PI);
  }
  /* Effective temperature */
  double GetEffectiveTemperature(double r, double angular_velocity) const {
    return std::sqrt(
        std::sqrt(
            GetEffectiveTemperature4(r, angular_velocity)
        )
    );
  }


private:
  /* Privates of things in constructor */
  double smbh_mass_;
  double eddington_ratio_;
  double alpha_;
  double X_mass_fraction_;
  double radiative_efficiency_;
};
