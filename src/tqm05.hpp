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
        supernovae_pressure_(supernovae_pressure),
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
        2.e5 * std::pow(smbh_mass_ / (1.3e8 * SOLAR_MASS_KG), 1.0/4.24):;"
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

};
