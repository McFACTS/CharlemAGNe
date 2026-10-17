/*
 *  File: test_sg03.cpp
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

//// Third Party ////
#include <catch2/catch_test_macros.hpp>

//// Standard Library ////

//////// Tests ////////

//// SG03 Defaults ////
TEST_CASE("SG03 Disk object defaults", "[sg03]") {
  auto A = SG03();
  REQUIRE(A.GetSMBHMass() == 1e8);
  REQUIRE(A.GetEddingtonRatio() == 0.5);
  REQUIRE(A.GetAlpha() == 0.01);
  REQUIRE(A.GetX() == 0.7);
  REQUIRE(A.GetRadiativeEfficiency() == 0.1);
}

//// SG03 Assign ////
TEST_CASE("SG03 Disk object value assignments", "[sg03]") {
  auto B = SG03(1.e7,  /* SMBH mass */
                1.0,   /* Eddington ratio */
                0.02,  /* Alpha */
                0.715, /* X mass abundance */
                0.2    /* radiative efficiency */
  );
  REQUIRE(B.GetSMBHMass() == 1e7);
  REQUIRE(B.GetEddingtonRatio() == 1.0);
  REQUIRE(B.GetAlpha() == 0.02);
  REQUIRE(B.GetX() == 0.715);
  REQUIRE(B.GetRadiativeEfficiency() == 0.2);
}
