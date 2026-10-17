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

class Disk1D {
public:
  /* Define this in daughter classes */
  Disk1D() {}
};

class SG03 : public Disk1D {
public:
  /* Constructor for SG03 */
  SG03(
      /* SMBH mass in solar masses */
      double smbh_mass = 1.0e8,
      /* Fraction of Eddington Luminosity */
      double eddington_ratio = 0.5,
      /* Shakura-Sunyaev parameter */
      double alpha = 0.01,
      /* Hydrogen abundance of gas in the disk */
      /* Fun fact, Cecilia Payne-Gaposchkin came up with this, and she
       * determined that stars are mostly hydrogen and helium. */
      double X_mass_fraction = 0.7,
      /* Radiative efficiency */
      double radiative_efficiency = 0.1)
      : smbh_mass_(smbh_mass), eddington_ratio_(eddington_ratio), alpha_(alpha),
        X_mass_fraction_(X_mass_fraction),
        radiative_efficiency_(radiative_efficiency) {}

  /* const is a promise that it does not modify anything in the object */
  /* accessors will have it. setters will not. */
  double GetSMBHMass() const { return smbh_mass_; }
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
        4.0 * M_PI * constants::G_SI * constants::C_SI * GetSMBHMass() *
        constants::SOLAR_MASS_KG / electron_scattering_opacity;
    /* Mdot at Eddington */
    double mdot_eddington =
        eddington_luminosity /
        (GetRadiativeEfficiency() * constants::C_SI * constants::C_SI);
    /* Mass accretion rate */
    return mdot_eddington * GetEddingtonRatio();
  }

private:
  /* Privates of things in constructor */
  double smbh_mass_;
  double eddington_ratio_;
  double alpha_;
  double X_mass_fraction_;
  double radiative_efficiency_;
};
