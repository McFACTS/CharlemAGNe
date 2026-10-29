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

TEST_CASE("TQM05 Disk object outer disk radius", "[tqm05]") {
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
  Parcel first_no_sf = disk.solve_T_rho(
        6.167141351355161e+18,
        1.93088819, // log10_T
        -15.04406584, // log10_rho
        disk.GetMDotOuter(),
        1,
        1.e-5
  );
  REQUIRE(first_no_sf.converged);
}
