/*
 *  File: utils.hpp
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
//// Standard Library ////
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <string>

//////// Setup /////////
using namespace constants;

// We are trying to replicate np.geomspace
inline std::vector<double> geomspace(
    double start,
    double stop,
    size_t num,
    bool endpoint = true) {
  // Initialize output
  std::vector<double> out;
  if (num == 0) return out;
  out.reserve(num);
  // Some of the loop logic
  const double log_start = std::log(start);
  const double log_stop  = std::log(stop);
  const size_t div       = endpoint ? num - 1 : num;
  const double step = (div > 0) ? (log_stop - log_start) / div : 0.0;
  // Loop through
  for (size_t i = 0; i < num; i++) {
    out.push_back(std::exp(log_start + i * step));
  }
  // Pin the endpoints 
  out.front() = start;
  if (endpoint && num > 1) out.back() = stop;
  return out;
};
