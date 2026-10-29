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
#include "sg03.hpp"
#include "tqm05.hpp"
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
  auto disk = TQM05(
    opac_name,
    opac_path,
    opac_id
  );

  //tests
  /* accessors will have it. setters will not. */
  printf("TQM05 disk modl\n");
  printf("disk smbh_mass: %g\n", disk.GetSMBHMassSolar());
  printf("velocity_dispersion (sigma): %g\n", disk.GetVelocityDispersion());
  printf("radiative efficiency (epsilon): %g\n", disk.GetRadiativeEfficiency());
  printf("transport efficiency (m): %g\n", disk.GetTransportEfficiency());
  printf("supernova pressure (xi): %g\n", disk.GetSupernovaPressure());
  printf("disk X: %g\n", disk.GetX());
  printf("MDotOuter: %g\n", disk.GetMDotOuter());
  printf("Inner radius: %g\n", disk.GetRadiusInner());
  printf("Outer radius: %g\n", disk.GetRadiusOuter());
  printf("Angular velocity outer (omega): %g\n", disk.GetAngularVelocityOuter());
  printf("Star Formation (fg0): %g\n", disk.GetStarFormationOuter());
  printf("Temperature outer: %g\n", disk.GetTemperatureOuter());
  printf("Density outer: %g\n", disk.GetRhoOuter());
  printf("Aspect ratio outer: %g\n", disk.GetAspectRatioOuter());
  // printf("disk eddington_ratio: %g\n", disk.GetEddingtonRatio());
  // printf("disk alpha: %g\n", disk.GetAlpha());
  // printf("disk Mdot: %g\n", disk.GetMDot());
  // double r_800 = 233981685308525.47;
  // double omega = disk.GetAngularVelocity(r_800);
  // double Teff4 = disk.GetEffectiveTemperature4(r_800, omega);
  // double Teff = disk.GetEffectiveTemperature(r_800, omega);
  // printf("omega (%g r_g) = %g\n", r_800, disk.GetAngularVelocity(r_800));
  // printf("Teff4 (%g r_g) = %g\n", r_800, Teff4);
  
  // Find point
  Parcel first_no_sf = disk.solve_T_rho(
        6.167141351355161e+18,
        1.93088819, // log10_T
        -15.04406584, // log10_rho
        disk.GetMDotOuter(),
        1,
        1.e-5
  );
  printf("first_no_sf converged: %d\n", first_no_sf.converged);
  printf("first_no_sf steps: %d\n", first_no_sf.steps);
  printf("r_cyl: %g\n", first_no_sf.r_cyl);
  printf("cs: %g\n", first_no_sf.cs);
  printf("T: %g\n", first_no_sf.T);
  printf("Teff: %g\n", first_no_sf.Teff);
  printf("rho: %g\n", first_no_sf.rho);
  printf("kappa: %g\n", first_no_sf.kappa);
  printf("tauV: %g\n", first_no_sf.tauV);
  printf("Q: %g\n", first_no_sf.Q);
  printf("h: %g\n", first_no_sf.h);
  // /* setup xguess */
  // double xguess_a = std::log10(omega * r_800 * 1.e-2);
  // double xguess_b = std::log10(Teff);// + 1.6;
  // printf("solve_cs_T:");
  // printf("x_guess (%g r_g) = %g, %g\n", r_800, xguess_a, xguess_b);
  // /* do the solve */
  // auto tic = std::chrono::high_resolution_clock::now();
  // Parcel parcel = disk.solve_cs_T(
  //   r_800,
  //   xguess_a,
  //   xguess_b
  // );
  // auto toc = std::chrono::high_resolution_clock::now();
  // auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(toc-tic);
  // printf("x_solve (%g r_g) = %g, %g\n", r_800, std::log10(parcel.cs), std::log10(parcel.T));
  // printf("time: %g ms\n", duration.count() / 1000000.0);
  // printf("converged: %d\n", parcel.converged);
  // printf("r_cyl: %g\n", parcel.r_cyl);
  // printf("cs: %g\n", parcel.cs);
  // printf("T: %g\n", parcel.T);
  // printf("Teff: %g\n", parcel.Teff);
  // printf("rho: %g\n", parcel.rho);
  // printf("kappa: %g\n", parcel.kappa);
  // printf("tauV: %g\n", parcel.tauV);
  // printf("Q: %g\n", parcel.Q);
  // printf("h: %g\n", parcel.h);

  // /* setup xguess */
  // xguess_a = std::log10(Teff);// + 1.6;
  // xguess_b = std::log10(Teff) + 1.6;
  // printf("solve_Teff_T:");
  // printf("x_guess (%g r_g) = %g, %g\n", r_800, xguess_a, xguess_b);
  // /* do the solve */
  // tic = std::chrono::high_resolution_clock::now();
  // parcel = disk.solve_Teff_T(
  //   r_800,
  //   xguess_a,
  //   xguess_b
  // );
  // toc = std::chrono::high_resolution_clock::now();
  // duration = std::chrono::duration_cast<std::chrono::nanoseconds>(toc-tic);
  // printf("x_solve (%g r_g) = %g, %g\n", r_800, std::log10(parcel.Teff), std::log10(parcel.T));
  // printf("time: %g ms\n", duration.count() / 1000000.0);
  // printf("converged: %d\n", parcel.converged);
  // printf("r_cyl: %g\n", parcel.r_cyl);
  // printf("cs: %g\n", parcel.cs);
  // printf("T: %g\n", parcel.T);
  // printf("Teff: %g\n", parcel.Teff);
  // printf("rho: %g\n", parcel.rho);
  // printf("kappa: %g\n", parcel.kappa);
  // printf("tauV: %g\n", parcel.tauV);
  // printf("Q: %g\n", parcel.Q);
  // printf("h: %g\n", parcel.h);

  return 0;
}
