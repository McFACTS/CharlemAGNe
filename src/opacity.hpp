/*
 *  File: opacity.hpp
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
//// Third Party ////
#include <singularity-opac/photons/mean_opacity_photons.hpp>
//// Standard Library ////
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <utility>


//////// Setup /////////
using namespace constants;
constexpr char CHARLEMAGNE_OPACITY_NAME[]   = "CHARLEMAGNE_OPACITY_NAME";
//std::string name = "pagn_combined",
constexpr char CHARLEMAGNE_OPACITY_PATH[]   = "CHARLEMAGNE_OPACITY_PATH";
//std::string path = "./pagn_combined.sp5",
constexpr char CHARLEMAGNE_OPACITY_ID[]     = "CHARLEMAGNE_OPACITY_ID";
//int opacid = 1001

//////// Objects /////////
class OpacityHandler {
public:
  /* Define this in daughter classes */
  OpacityHandler() {}
};

class PAGNCombined : public OpacityHandler {
public:
  /* Constructor for PAGNCombined OpacityHandler */
  PAGNCombined(
    const std::string name,
    const std::string path,
    int opacid
  ) :
    name_(std::move(name)),
    path_(std::move(path)),
    opacid_(std::move(opacid)),
    mean_opac_cgs_(path_, name_)
  {
    // Make sure the filepath checks out
    if (!std::filesystem::is_regular_file(path_)) {
        throw std::runtime_error("File not found: " + path_);
    }
  }

  /* Constructor from env */
  PAGNCombined() {
    const char *opacid_char = std::getenv(CHARLEMAGNE_OPACITY_ID);
    if (!std::getenv(CHARLEMAGNE_OPACITY_NAME)) {
        throw std::runtime_error("Must define CHARLEMAGNE_OPACITY_NAME");
    }
    if (!std::getenv(CHARLEMAGNE_OPACITY_PATH)) {
        throw std::runtime_error("Must define CHARLEMAGNE_OPACITY_PATH");
    }
    if (!std::getenv(CHARLEMAGNE_OPACITY_ID)) {
        throw std::runtime_error("Must define CHARLEMAGNE_OPACITY_ID");
    }
    name_ = std::getenv(CHARLEMAGNE_OPACITY_NAME);
    path_ = std::getenv(CHARLEMAGNE_OPACITY_PATH);
    opacid_ = std::stoi(std::getenv(CHARLEMAGNE_OPACITY_ID));
    //MeanOpacityBase mean_opac_cgs_(path, name, false);
    mean_opac_cgs_ = singularity::photons::MeanOpacityBase(path_, name_);
    //printf("name_: %s; path_: %s; opacid: %d\n", name_, path_, opacid_);
    // Make sure at least one char
    if (name_.empty()) {
        throw std::runtime_error("Cannot load from env: " + std::string(CHARLEMAGNE_OPACITY_NAME));
    }
    if (path_.empty()) {
        throw std::runtime_error("Cannot load from env: " + std::string(CHARLEMAGNE_OPACITY_PATH));
    }
    // Make sure the filepath checks out
    if (!std::filesystem::is_regular_file(path_)) {
        throw std::runtime_error("File not found: " + path_);
    }
  }

  // Return the name of the opacity table
  std::string GetName() const { return name_; }

  // Return the value of the opacity
  double eval(double rho_si, double temp) {
    const double logRopac = std::log10(rho_si) - 3.0 * std::log10(temp) + 15.0;
    const double fake_rho = (rho_si / (temp * temp * temp)) * 1.e15;
    double value = 0.1 * mean_opac_cgs_.RosselandMeanAbsorptionCoefficient(fake_rho, temp) / fake_rho;
    return value;
  }

private:
    std::string name_;
    std::string path_;
    int opacid_;
    singularity::photons::MeanOpacityBase mean_opac_cgs_;
};
