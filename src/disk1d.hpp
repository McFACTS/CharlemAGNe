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
    // If we successfully solved the disk at this radius
    bool converged = false;
    // Number of steps iterated
    int steps = 0;
    // Star formation flag
    bool star_formation = false;
    // The distance from the axis of the disk to the point
    double r_cyl;
    // The angular velocity
    double omega = 0.0;
    // Effective temperature
    double Teff = 0.0;
    // Mdotprime
    double Mdotprime = 0.0;
    // Midplane temperature
    double T = 0.0;
    // Sound speed
    double cs = 0.0;
    // density
    double rho = 0.0;
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

#define HOT_LOOP_ARGS Parcel& state, double logX0, double logX1, double &cost0, double &cost1

//// Disks ////
class Disk1D {
public:
  /* Define this in daughter classes */
  Disk1D() {}

protected:
  template<auto HotLoop, class D>
  void calc_inv_jacobian(
    D& self,
    Parcel &parcel,
    double jac_inv[4],
    double log10X0,
    double log10X1,
    double &cost0, 
    double &cost1
  ) {
    constexpr double eps = 1e-5;
    double cost0_m0, cost1_m0, cost0_m1, cost1_m1;
    (self.*HotLoop)(parcel, log10X0 - eps, log10X1, cost0_m0, cost1_m0);
    (self.*HotLoop)(parcel, log10X0, log10X1 - eps, cost0_m1, cost1_m1);
    (self.*HotLoop)(parcel, log10X0, log10X1, cost0, cost1);

    double jac[4] = {(cost0 - cost0_m0) / eps, (cost0 - cost0_m1) / eps,
                     (cost1 - cost1_m0) / eps, (cost1 - cost1_m1) / eps};
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
};

