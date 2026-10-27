/*
 *  File: main.cpp
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

//
//////// Include ////////
//// Local ////
#include "constants.hpp"
#include "disk1d.hpp"
#include "opacity.hpp"
//// Third Party ////
//// Standard Library ////
#include <cstdio>
#include <iostream>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <chrono>

constexpr char CHARLEMAGNE_OPACITY_NAME[]   = "CHARLEMAGNE_OPACITY_NAME";
//std::string name = "pagn_combined",
constexpr char CHARLEMAGNE_OPACITY_PATH[]   = "CHARLEMAGNE_OPACITY_PATH";
//std::string path = "./pagn_combined.sp5",
constexpr char CHARLEMAGNE_OPACITY_ID[]     = "CHARLEMAGNE_OPACITY_ID";
//int opacid = 1001

int main(int argc, char *argv[]) {
  // Get arguments
  if (!std::getenv(CHARLEMAGNE_OPACITY_NAME)) {
      throw std::runtime_error("Must define CHARLEMAGNE_OPACITY_NAME");
  }
  if (!std::getenv(CHARLEMAGNE_OPACITY_PATH)) {
      throw std::runtime_error("Must define CHARLEMAGNE_OPACITY_PATH");
  }
  if (!std::getenv(CHARLEMAGNE_OPACITY_ID)) {
      throw std::runtime_error("Must define CHARLEMAGNE_OPACITY_ID");
  }
  std::string opac_name = std::getenv(CHARLEMAGNE_OPACITY_NAME);
  std::string opac_path = std::getenv(CHARLEMAGNE_OPACITY_PATH);
  int opac_id = std::stoi(std::getenv(CHARLEMAGNE_OPACITY_ID));

  // setup disk
  auto disk = SG03(
    opac_name,
    opac_path,
    opac_id,
    1.0e8,
    0.5,
    0.01,
    0.7,
    0.1,
    false
  );

  //tests
  printf("disk smbh_mass: %g\n", disk.GetSMBHMassSolar());
  printf("disk eddington_ratio: %g\n", disk.GetEddingtonRatio());
  printf("disk alpha: %g\n", disk.GetAlpha());
  printf("disk X: %g\n", disk.GetX());
  printf("disk Mdot: %g\n", disk.GetMDot());
  double r_800 = 233981685308525.47;
  double omega = disk.GetAngularVelocity(r_800);
  double Teff4 = disk.GetEffectiveTemperature4(r_800, omega);
  double Teff = disk.GetEffectiveTemperature(r_800, omega);
  printf("omega (%g r_g) = %g\n", r_800, disk.GetAngularVelocity(r_800));
  printf("Teff4 (%g r_g) = %g\n", r_800, Teff4);
  
  /* setup xguess */
  double xguess_a = std::log10(omega * r_800 * 1.e-2);
  double xguess_b = std::log10(Teff);// + 1.6;
  printf("solve_cs_T:");
  printf("x_guess (%g r_g) = %g, %g\n", r_800, xguess_a, xguess_b);
  /* do the solve */
  auto tic = std::chrono::high_resolution_clock::now();
  Parcel parcel = disk.solve_cs_T(
    r_800,
    xguess_a,
    xguess_b
  );
  auto toc = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(toc-tic);
  printf("x_solve (%g r_g) = %g, %g\n", r_800, std::log10(parcel.cs), std::log10(parcel.T));
  printf("time: %g ms\n", duration.count() / 1000000.0);
  printf("converged: %d\n", parcel.converged);
  printf("r_cyl: %g\n", parcel.r_cyl);
  printf("cs: %g\n", parcel.cs);
  printf("T: %g\n", parcel.T);
  printf("Teff: %g\n", parcel.Teff);
  printf("rho: %g\n", parcel.rho);
  printf("kappa: %g\n", parcel.kappa);
  printf("tauV: %g\n", parcel.tauV);
  printf("Q: %g\n", parcel.Q);
  printf("h: %g\n", parcel.h);

  /* setup xguess */
  xguess_a = std::log10(Teff);// + 1.6;
  xguess_b = std::log10(Teff) + 1.6;
  printf("solve_Teff_T:");
  printf("x_guess (%g r_g) = %g, %g\n", r_800, xguess_a, xguess_b);
  /* do the solve */
  tic = std::chrono::high_resolution_clock::now();
  parcel = disk.solve_Teff_T(
    r_800,
    xguess_a,
    xguess_b
  );
  toc = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::nanoseconds>(toc-tic);
  printf("x_solve (%g r_g) = %g, %g\n", r_800, std::log10(parcel.Teff), std::log10(parcel.T));
  printf("time: %g ms\n", duration.count() / 1000000.0);
  printf("converged: %d\n", parcel.converged);
  printf("r_cyl: %g\n", parcel.r_cyl);
  printf("cs: %g\n", parcel.cs);
  printf("T: %g\n", parcel.T);
  printf("Teff: %g\n", parcel.Teff);
  printf("rho: %g\n", parcel.rho);
  printf("kappa: %g\n", parcel.kappa);
  printf("tauV: %g\n", parcel.tauV);
  printf("Q: %g\n", parcel.Q);
  printf("h: %g\n", parcel.h);
  // R: 233981685308525.47
  // cs: 55589.63200337918
  // T: 70383.14284818874
  // Teff: 3603.107000477657
  // rho: 2.470550e-06
  // kappa: 9.099875e-02
  // tauV: 3.882695e+05
  // Q: 1.000000e+00
  // h: 1.727047e+12

  return 0;
}
