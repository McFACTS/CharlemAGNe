/*
 *  File: test_tqm05.cpp
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

//////// Include ////////
//// Local ////
#include "disk1d.hpp"
#include "tqm05.hpp"

//// Third Party ////
#include <catch2/catch_test_macros.hpp>

//// Standard Library ////

//////// Setup ////////
constexpr char CHARLEMAGNE_OPACITY_NAME[]   = "CHARLEMAGNE_OPACITY_NAME";
//std::string name = "pagn_combined",
constexpr char CHARLEMAGNE_OPACITY_PATH[]   = "CHARLEMAGNE_OPACITY_PATH";
//std::string path = "./pagn_combined.sp5",
constexpr char CHARLEMAGNE_OPACITY_ID[]     = "CHARLEMAGNE_OPACITY_ID";
//int opacid = 1001

//////// Tests ////////

//// TQM05 Defaults ////
TEST_CASE("TQM05 Disk object defaults", "[tqm05]") {
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
  auto disk = TQM05(opac_name, opac_path, opac_id);
  // pAGN
  // sigma: 187999.4434187186
  // epsilon: 0.001
  // m: 0.2
  // xi: 1.0
  // Mdot_out: 2.0163459832179886e+25
  // Rin: 886001814620.0636
  // Rout: 6.172e+18
  // Omega: 4.376226558685041e-14
  // fg0: 0.9839371036124078
  // T0: 120.75098750027597
  // rho0: 6.4584603958072986e-18
  // h0: 2.1470801741235356e+18

  // main.cpp
  // velocity_dispersion (sigma): 187999
  // radiative efficiency (epsilon): 0.001
  // transport efficiency (m): 0.2
  // supernova pressure (xi): 1
  // disk X: 0.7
  // MDotOuter: 2.0169e+25
  // Inner radius: 8.85975e+11
  // Outer radius: 6.16934e+18
  // Angular velocity outer (omega): 4.37464e-14
  // Star Formation (f0): 0.984081
  // Temperature outer: 120.738
  // Density outer: 6.45379e-18
  // Aspect ratio outer: 2.14647e+18

  REQUIRE((std::abs(disk.GetSMBHMassSolar() - 1.e8) / 1.e8) < 1.e-7);
  REQUIRE((std::abs(disk.GetVelocityDispersion() - 187999.4434187186) / 187999.4434187186) < 1.e-5);
  REQUIRE((std::abs(disk.GetRadiativeEfficiency() - 0.001) / 0.001) < 1.e-6);
  REQUIRE((std::abs(disk.GetTransportEfficiency() - 0.2) / 0.2) < 1.e-6);
  REQUIRE((std::abs(disk.GetSupernovaPressure() - 1.0) / 1.0) < 1.e-6);
  REQUIRE((std::abs(disk.GetX() - 0.7) / 0.7) < 1.e-8);
  REQUIRE((std::abs(disk.GetMDotOuter() - 2.0163459832179886e+25) / 2.0163459832179886e+25) < 1.e-3);
  REQUIRE((std::abs(disk.GetRadiusInner() - 886001814620.0636) / 886001814620.0636) < 1.e-3);
  REQUIRE((std::abs(disk.GetRadiusOuter() - 6.172e+18) / 6.172e+18) < 1.e-3);
  REQUIRE((std::abs(disk.GetAngularVelocityOuter() - 4.376226558685041e-14) / 4.376226558685041e-14) < 1.e-3);
  REQUIRE((std::abs(disk.GetStarFormationOuter() - 0.9839371036124078) / 0.9839371036124078) < 1.e-3);
  REQUIRE((std::abs(disk.GetTemperatureOuter() - 120.75098750027597) / 120.75098750027597) < 1.e-3);
  REQUIRE((std::abs(disk.GetRhoOuter() - 6.4584603958072986e-18) / 6.4584603958072986e-18) < 1.e-3);
  REQUIRE((std::abs(disk.GetAspectRatioOuter() - 2.1470801741235356e+18) / 2.1470801741235356e+18) < 1.e-3);
  REQUIRE((std::abs(disk.GetKappaOuter() - 0.17684028668423085) / 0.17684028668423085) < 1.e-1);
  REQUIRE((std::abs(disk.GetTauVOuter() - 2.37789229981872) / 2.37789229981872) < 1.e-1);
  REQUIRE((std::abs(disk.GetTeff4Outer() - 71015968.52973816) / 71015968.52973816) < 1.e-1);
  REQUIRE((std::abs(disk.GetSFEfficiencyOuter() - 0.07389095846905776) / 0.07389095846905776) < 1.e-1);
}

TEST_CASE("TQM05 Disk object first steps", "[tqm05]") {
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
  auto disk = TQM05(opac_name, opac_path, opac_id);
  // Find point
  Parcel first_sf = disk.solve_T_eta(
        6.167141351355161e+18,
        std::log10(disk.GetTemperatureOuter()),
        std::log10(disk.GetSFEfficiencyOuter()),
        disk.GetMDotOuter()
  );
  REQUIRE(first_sf.converged);
  // pAGN
  // cs: 93887.12587907449
  // T: 119.79945768699307
  // Teff: 90.83386987173469
  // rho: 6.4584603958072986e-18
  // kappa: 0.17510017298606326
  // tauV: 2.426173321390316
  // eta: 0.07074697050876151
  // h: 2.1453899751315773e+18
  //
  // main.cpp
  // r_cyl: 6.16714e+18
  // cs: 93900
  // T: 119.786
  // Teff: 89.7629
  // rho: 6.45845e-18
  // kappa: 0.190705
  // tauV: 2.64276
  // Q: 1
  // h: 2.14569e+18
  REQUIRE((std::abs(first_sf.cs - 93887.12587907449) / 93887.12587907449) < 1.e-3);
  REQUIRE((std::abs(first_sf.T - 119.79945768699307) / 119.79945768699307) < 1.e-3);
  REQUIRE((std::abs(first_sf.Teff - 90.83386987173469) / 90.83386987173469) < 1.e-1);
  REQUIRE((std::abs(first_sf.rho - 6.4584603958072986e-18) / 6.4584603958072986e-18) < 1.e-3);
  REQUIRE((std::abs(first_sf.kappa - 0.17510017298606326) / 0.17510017298606326) < 1.e-1);
  REQUIRE((std::abs(first_sf.tauV - 2.426173321390316) / 2.426173321390316) < 1.e-1);
  REQUIRE((std::abs(first_sf.h - 2.1453899751315773e+18) / 2.1453899751315773e+18) < 1.e-3);
  REQUIRE((std::abs(first_sf.eta - 0.07074697050876151) / 0.07074697050876151) < 1.e-1);
  // Well, it agrees in each quantity to at least 10%, 
  // even using an entirely different opacity interpolation scheme.
  // Find second point
  Parcel second_sf = disk.solve_T_eta(
    6.157431703601016e+18,
    std::log10(first_sf.T),
    std::log10(first_sf.eta),
    2.0163459832179886e+25
  );
  /* Second  test
   * pAGN
   * R: 6.157431703601016e+18
   * Mdot: 2.0163459832179886e+25
   * cs: 93886.03239932859
   * T: 119.89355934263655
   * Teff: 90.87602945979515
   * rho: 6.479146973372726e-18
   * kappa: 0.17524372784116451
   * tauV: 2.432019697365853
   * eta: 0.07065281725968725
   * h: 2.1419373966624404e+18
   *
   * main.cpp
   * second_sf converged: 1
   * second_sf steps: 2
   * r_cyl: 6.15743e+18
   * cs: 93886.1
   * T: 119.872
   * Teff: 89.8078
   * rho: 6.47914e-18
   * kappa: 0.190731
   * tauV: 2.64696
   * Q: 1
   * h: 2.14194e+18
  */
  REQUIRE(second_sf.converged);
  REQUIRE((std::abs(second_sf.cs - 93886.03239932859) / 93886.03239932859) < 1.e-3);
  REQUIRE((std::abs(second_sf.T - 119.89355934263655) / 119.89355934263655) < 1.e-3);
  REQUIRE((std::abs(second_sf.Teff - 90.87602945979515) / 90.87602945979515) < 1.e-1);
  REQUIRE((std::abs(second_sf.rho - 6.479146973372726e-18) / 6.479146973372726e-18) < 1.e-3);
  REQUIRE((std::abs(second_sf.kappa - 0.17524372784116451) / 0.17524372784116451) < 1.e-1);
  REQUIRE((std::abs(second_sf.tauV - 2.432019697365853) /2.4320196973658536) < 1.e-1);
  REQUIRE((std::abs(second_sf.h - 2.1419373966624404e+18) / 2.1419373966624404e+18) < 1.e-4);
  REQUIRE((std::abs(first_sf.eta - 0.07065281725968725) / 0.07065281725968725) < 1.e-1);
}

TEST_CASE("TQM05 Disk object no sf step", "[tqm05]") {
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
  auto disk = TQM05(opac_name, opac_path, opac_id);
  // Find point
  // Try a solve with no star formation
  Parcel no_sf = disk.solve_T_rho(
    2285812867742.8877,
    std::log10(240555.264243547),
    std::log10(6.558288583390488e-10),
    1.406323456972877e+23
  );
  REQUIRE(no_sf.converged);
  /*
  pAGN:
  R: 2285812867742.8877
  Mdot: 1.406323456972877e+23
  cs: 35177765.141558014
  T: 240791.1076947123
  Teff: 105559.22962529745
  rho: 6.594353356182405e-10
  kappa: 0.049933921917359105
  tauV: 34.74818464029484
  eta: 0.0
  h: 1055271545158.8575

  main.cpp:
  r_cyl: 2.28581e+12
  cs: 1.51148e+07
  T: 240790
  Teff: 105559
  rho: 3.57187e-09
  kappa: 0.0214548
  tauV: 34.7476
  Q: 1.04913e+09
  h: 4.53425e+11
  */
  REQUIRE((std::abs(no_sf.T - 240791.1076947123) / 240791.1076947123) < 1.e-3);
  REQUIRE((std::abs(no_sf.Teff - 105559.22962529745) / 105559.22962529745) < 1.e-3);
  REQUIRE((std::abs(no_sf.tauV - 34.74818464029484) / 34.74818464029484) < 1.e-3);
  // This is way off the edge of the opacity table, 
  //  so this is how pAGN and singularity-opac extrapolate differently.
  //REQUIRE((std::abs(no_sf.cs - 35177765.141558014) / 35177765.141558014) < 1.e-3);
  //REQUIRE((std::abs(no_sf.rho - 6.594353356182405e-10) / 6.594353356182405e-10) < 1.e-3);
  //REQUIRE((std::abs(no_sf.kappa - 0.049933921917359105) / 0.049933921917359105) < 1.e-1);
  //REQUIRE((std::abs(no_sf.h - 1055271545158.8575) / 1055271545158.8575) < 1.e-3);
}

TEST_CASE("TQM05 Default solve doesn't crash", "[tqm05]") {
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
  auto tqm05_disk = TQM05(opac_name ,opac_path, opac_id);
  auto tqm05_profile = tqm05_disk.SolveDisk();
}

TEST_CASE("TQM05 Standard options converges", "[tqm05]") {
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
  auto tqm05_disk = TQM05(opac_name ,opac_path, opac_id);
  auto opts = SolverOpts(
    200 /*maxfev*/,
    1.e-5 /*xtol*/,
    1.e-8 /*ftol*/,
    0.0 /*gtol*/,
    1.e-7 /*epsfcn*/,
    0.1 /*factor*/,
    2.0 /*conv_tol*/
  );
  auto tqm05_profile = tqm05_disk.SolveDisk(6.0, 0.0, 10000, opts);
  for (size_t i = 0; i < tqm05_profile.r_cyl.size(); i++) {
    REQUIRE(tqm05_profile.max_residual[i] < 10.0);
  }
}

TEST_CASE("TQM05 Outer disk", "[tqm05]") {
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
  auto tqm05_disk = TQM05(opac_name ,opac_path, opac_id);
  auto opts = SolverOpts(
    200 /*maxfev*/,
    1.e-5 /*xtol*/,
    1.e-8 /*ftol*/,
    0.0 /*gtol*/,
    1.e-7 /*epsfcn*/,
    0.1 /*factor*/,
    10.0 /*conv_tol*/
  );
  auto tqm05_profile = tqm05_disk.SolveDisk(6.0, 1.e8, 10000, opts);
  for (size_t i = 0; i < tqm05_profile.r_cyl.size(); i++) {
    REQUIRE(tqm05_profile.max_residual[i] < 10.0);
  }
}
