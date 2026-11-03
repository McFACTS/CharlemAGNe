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
#include <cminpack.h> // for lmdif/lmder root finders
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
    double r_cyl = 0.0;
    // The angular velocity
    double omega = 0.0;
    // Effective temperature
    double Teff = 0.0;
    // Mdotprime
    double Mdotprime = 0.0;
    // Mdot (for TQM05)
    double Mdot = 0.0;
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
    // Star formation efficiency
    double eta = 0.0;

    //// Constructor
    Parcel() = default;
    explicit Parcel(double r_cyl) : r_cyl(r_cyl) {};
};

/// A profile of gas
struct Profile {
    //// Members 
    // If we successfully solved the disk at this radius
    std::vector<bool> converged;
    // Number of steps iterated
    std::vector<int> steps;
    // Star formation flag
    std::vector<bool> star_formation;
    // The distance from the axis of the disk to the point
    std::vector<double> r_cyl;
    // The angular velocity
    std::vector<double> omega;
    // Effective temperature
    std::vector<double> Teff;
    // Mdotprime
    std::vector<double> Mdotprime;
    // Mdot (for TQM05)
    std::vector<double> Mdot;
    // Midplane temperature
    std::vector<double> T;
    // Sound speed
    std::vector<double> cs;
    // density
    std::vector<double> rho;
    // opacity
    std::vector<double> kappa;
    std::vector<double> tauV;
    // pressure
    std::vector<double> PPgas;
    std::vector<double> PPrad;
    // opacfac
    std::vector<double> opacfac;
    // Toomre instability
    std::vector<double> Q;
    // Scale height
    std::vector<double> h;
    // Star formation efficiency
    std::vector<double> eta;

    //// Profile specific things
    std::vector<double> dr_cyl;
    std::vector<double> dA;
    std::vector<double> dV;
    std::vector<double> mass_disk;
    std::vector<double> Mdot_sf;

    //// Constructor
    Profile(std::vector<double> r_cyl) : r_cyl(r_cyl) {
      const size_t n = r_cyl.size();
      // Non doubles
      converged.assign(n, false);
      steps.assign(n, 0);
      star_formation.assign(n, false);
      // Doubles
      omega.resize(n);
      Teff.resize(n);
      Mdotprime.resize(n);
      Mdot.resize(n);
      T.resize(n);
      cs.resize(n);
      rho.resize(n);
      kappa.resize(n);
      tauV.resize(n);
      PPgas.resize(n);
      PPrad.resize(n);
      opacfac.resize(n);
      Q.resize(n);
      h.resize(n);
      eta.resize(n);
      dr_cyl.resize(n);
      dA.resize(n);
      dV.resize(n);
      mass_disk.resize(n);
      Mdot_sf.resize(n);
    };
};

//// Solver options ////
struct SolverOpts {
  int maxfev = 200;         // Function evaluations
  double ftol = 1.e-5;        // relative f(x) reduction
  double xtol = 1.e-5;        // relative step size (dx)
  double gtol = 0.0;        // Orthogonality test
  double epsfcn = 1.e-5;    // Step size for estimating Jacobian
  double factor = 0.1;      // initial "trust-region" boundary.
  double conv_tol = 1.e-1; // Check if we converged
};

//// wrapper ////
// This is for passing the parcel through to the hot loop function
template<auto HotLoop, class D>
struct LmdifCtx {
    D* self; // D as in disk
    Parcel* parcel;
};

template<auto HotLoop, class D>
int lmdif_fcn(
    void* p,
    int /*m*/,
    int /*n*/,
    const double* x,
    double* fvec,
    int /*iflag*/) {
  // Define a LmdifCtx using the templated types
  auto* ctx = static_cast<LmdifCtx<HotLoop, D>*>(p);
  // This is a call of the HotLoop function
  (ctx->self->*HotLoop)(*ctx->parcel, x[0], x[1], fvec[0], fvec[1]);
  return 0;
}

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

  // Run cminpack lmdif on the 2D hot loop
  template<auto HotLoop, class D>
  int solve_lmdif(
    D& self,
    Parcel& parcel,
    double x[2],
    const SolverOpts& opts,
    int& nfev) {
    // These are all 2D problems
    constexpr int m = 2, n = 2;
    // Variables that cminpack needs
    double fvec[m];
    double diag[n], fjac[m * n], qtf[n];
    double wa1[n], wa2[n], wa3[n], wa4[n];
    int ipvt[n];
    diag[0] = diag[1] = 1.0; // mode=1 user-supplied scaling
    // Create the initial state by running LmdifCtx onece for the guess
    LmdifCtx<HotLoop, D> ctx{&self, &parcel};
    // Do the solve
    int info = lmdif(
      lmdif_fcn<HotLoop, D>,
      &ctx,
      m, n, x, fvec,
      opts.ftol, opts.xtol, opts.gtol, opts.maxfev, opts.epsfcn,
      diag, 1 /*mode*/, opts.factor, 0 /*nprint*/, &nfev,
      fjac, m /*ldfjac*/, ipvt, qtf, wa1, wa2, wa3, wa4
    );
    return info;
    }
};

