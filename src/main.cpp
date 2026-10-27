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
  auto disk = SG03(opac_name, opac_path, opac_id);

  //tests
  printf("disk smbh_mass: %g\n", disk.GetSMBHMassSolar());
  printf("disk eddington_ratio: %g\n", disk.GetEddingtonRatio());
  printf("disk alpha: %g\n", disk.GetAlpha());
  printf("disk X: %g\n", disk.GetX());
  printf("disk Mdot: %g\n", disk.GetMDot());
  double r_12 = 12.0 * disk.GetRg();
  double omega = disk.GetAngularVelocity(r_12);
  double Teff4 = disk.GetEffectiveTemperature4(r_12, omega);
  double Teff = disk.GetEffectiveTemperature(r_12, omega);
  printf("omega (%g r_g) = %g\n", r_12, disk.GetAngularVelocity(r_12));
  printf("Teff4 (%g r_g) = %g\n", r_12, Teff4);
  
  /* setup xguess */
  double xguess_a = std::log10(omega * r_12 * 1.e-2);
  double xguess_b = std::log10(Teff);
  printf("x_guess (%g r_g) = %g, %g\n", r_12, xguess_a, xguess_b);
  /* do the solve */
  auto tic = std::chrono::high_resolution_clock::now();
  Parcel parcel = disk.solve_cs_T(r_12);//, xguess_a, xguess_b);
  auto toc = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(toc-tic);
  printf("x_solve (%g r_g) = %g, %g\n", r_12, std::log10(parcel.cs), std::log10(parcel.T));
  printf("time: %g ms\n", duration.count() / 1000000.0);
  printf("r_cyl: %g\n", parcel.r_cyl);
  printf("rho: %g\n", parcel.rho);
  printf("kappa: %g\n", parcel.kappa);
  printf("tauV: %g\n", parcel.tauV);
  printf("Q: %g\n", parcel.Q);
  printf("h: %g\n", parcel.h);


  return 0;
}
