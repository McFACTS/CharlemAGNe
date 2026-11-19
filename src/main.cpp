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
#include "sg03.hpp"
#include "tqm05.hpp"
#include "utils.hpp"
//// Third Party ////
//// Standard Library ////
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

constexpr char CHARLEMAGNE_OPACITY_NAME[] = "CHARLEMAGNE_OPACITY_NAME";
// std::string name = "pagn_combined",
constexpr char CHARLEMAGNE_OPACITY_PATH[] = "CHARLEMAGNE_OPACITY_PATH";
// std::string path = "./pagn_combined.sp5",
constexpr char CHARLEMAGNE_OPACITY_ID[] = "CHARLEMAGNE_OPACITY_ID";
// int opacid = 1001

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
  auto sg03_disk = SG03(opac_name, opac_path, opac_id);
  auto tic = std::chrono::high_resolution_clock::now();
  sg03_disk.WriteDisk("./charlemagne_sg03.csv");
  auto toc = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::nanoseconds>(toc - tic);
  printf("Solved SG03 disk in %g s\n", duration.count() / 1.e9);

  // Run TQM disk
  auto tqm05_disk = TQM05(opac_name, opac_path, opac_id);
  tic = std::chrono::high_resolution_clock::now();
  tqm05_disk.WriteDisk("./charlemagne_tqm05.csv");
  toc = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::nanoseconds>(toc - tic);
  printf("Solved TQM05 disk in %g s\n", duration.count() / 1.e9);
  return 0;
}
