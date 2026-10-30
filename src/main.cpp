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
#include "utils.hpp"
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
  std::FILE* sg03_csv = std::fopen("./charlemagne_sg03.csv", "w");

  // setup disk
  auto disk = SG03(
    opac_name,
    opac_path,
    opac_id
  );
  // r_space
  const double rmin = 6.0 * disk.GetRg();
  const double rmax = 1.e7 * disk.GetRg();
  const double npts = 1.e4;
  std::vector<double> r_space = geomspace(rmin, rmax, npts);
  auto tic = std::chrono::high_resolution_clock::now();
  Profile sg03_profile = disk.WalkDisk(r_space);
  auto toc = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(toc-tic);
  printf("Solved disk in %g s\n", duration.count() / 1.e9);
  fprintf(sg03_csv, "#iter, converged, steps, star_formation, omega, Teff, Mdotprime, Mdot, T, cs, rho, kappa, tauV, PPgas, PPrad, opacfac, Q, h, eta, dr, dA, dV, mass_disk, Mdot_sf\n");
  for (size_t i = 0; i < npts; i++) {
    fprintf(
        sg03_csv,
        "%lu %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
        i,
        sg03_profile.converged[i],
        sg03_profile.steps[i],
        sg03_profile.star_formation[i],
        sg03_profile.omega[i],
        sg03_profile.Teff[i],
        sg03_profile.Mdotprime[i],
        sg03_profile.Mdot[i],
        sg03_profile.T[i],
        sg03_profile.cs[i],
        sg03_profile.rho[i],
        sg03_profile.kappa[i],
        sg03_profile.tauV[i],
        sg03_profile.PPgas[i],
        sg03_profile.PPrad[i],
        sg03_profile.opacfac[i],
        sg03_profile.Q[i],
        sg03_profile.h[i],
        sg03_profile.eta[i],
        sg03_profile.dr_cyl[i],
        sg03_profile.dA[i],
        sg03_profile.dV[i],
        sg03_profile.mass_disk[i],
        sg03_profile.Mdot_sf[i]
    );
  }
  std::fclose(sg03_csv);
  
  return 0;
}
