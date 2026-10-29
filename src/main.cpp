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
  printf("Kappa outer: %g\n", disk.GetKappaOuter());
  printf("TauV outer: %g\n", disk.GetTauVOuter());
  printf("Teff4 outer: %g\n", disk.GetTeff4Outer());
  printf("SFEfficiency outer: %g\n", disk.GetSFEfficiencyOuter());
  
  // Find point
  Parcel first_sf = disk.solve_T_eta(
        6.167141351355161e+18,
        std::log10(disk.GetTemperatureOuter()),
        std::log10(disk.GetSFEfficiencyOuter()),
        disk.GetMDotOuter(),
        30,
        1.e-5
  );
  printf("first_sf converged: %d\n", first_sf.converged);
  printf("first_sf steps: %d\n", first_sf.steps);
  printf("r_cyl: %g\n", first_sf.r_cyl);
  printf("cs: %g\n", first_sf.cs);
  printf("T: %g\n", first_sf.T);
  printf("Teff: %g\n", first_sf.Teff);
  printf("rho: %g\n", first_sf.rho);
  printf("kappa: %g\n", first_sf.kappa);
  printf("tauV: %g\n", first_sf.tauV);
  printf("Q: %g\n", first_sf.Q);
  printf("h: %g\n", first_sf.h);

  // Find second point
  Parcel second_sf = disk.solve_T_eta(
    6.157431703601016e+18,
    std::log10(first_sf.T),
    std::log10(first_sf.eta),
    2.0163459832179886e+25,
    30,
    1.e-5
  );
  printf("second_sf converged: %d\n", second_sf.converged);
  printf("second_sf steps: %d\n", second_sf.steps);
  printf("r_cyl: %g\n", second_sf.r_cyl);
  printf("cs: %g\n", second_sf.cs);
  printf("T: %g\n", second_sf.T);
  printf("Teff: %g\n", second_sf.Teff);
  printf("rho: %g\n", second_sf.rho);
  printf("kappa: %g\n", second_sf.kappa);
  printf("tauV: %g\n", second_sf.tauV);
  printf("Q: %g\n", second_sf.Q);
  printf("h: %g\n", second_sf.h);

  // Try a solve with no star formation
  Parcel no_sf = disk.solve_T_rho(
    2285812867742.8877,
    std::log10(240555.264243547),
    std::log10(6.558288583390488e-10),
    1.406323456972877e+23,
    30,
    1.e-5
  );
  printf("no_sf converged: %d\n", no_sf.converged);
  printf("no_sf steps: %d\n", no_sf.steps);
  printf("r_cyl: %g\n", no_sf.r_cyl);
  printf("cs: %g\n", no_sf.cs);
  printf("T: %g\n", no_sf.T);
  printf("Teff: %g\n", no_sf.Teff);
  printf("rho: %g\n", no_sf.rho);
  printf("kappa: %g\n", no_sf.kappa);
  printf("tauV: %g\n", no_sf.tauV);
  printf("Q: %g\n", no_sf.Q);
  printf("h: %g\n", no_sf.h);

  return 0;
}
