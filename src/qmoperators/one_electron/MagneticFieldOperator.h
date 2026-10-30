/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2023 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

#pragma once

#include "tensor/RankZeroOperator.h"

#include "PositionOperator.h"

/** @class MagneticFieldOperator
 *
 * @brief Real component of the magnetic field operator
 *
 * Diamagnetic interaction between the external magnetic field
 * and the electrons.
 *
 */

namespace mrchem {

class MagneticFieldOperator final : public RankZeroOperator {
public:
    MagneticFieldOperator(const Eigen::Vector3d &f, const mrcpp::Coord<3> &o)
            : MagneticFieldOperator(std::array<double, 3>{f[0], f[1], f[2]}, PositionOperator(o)) {}

    MagneticFieldOperator(const Eigen::Vector3d &f, PositionOperator r)
            : MagneticFieldOperator(std::array<double, 3>{f[0], f[1], f[2]}, r) {}

    MagneticFieldOperator(const std::array<double, 3> &f, const mrcpp::Coord<3> &o)
            : MagneticFieldOperator(f, PositionOperator(o)) {}

    MagneticFieldOperator(const std::array<double, 3> &f, PositionOperator r) {
        /*auto B2 = f[0] * f[0] + f[1] * f[1] + f[2] * f[2];
        auto r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
        auto Br = f[0] * r[0] + f[1] * r[1] + f[2] * r[2];
        RankZeroOperator MF_e = 0.125 * (B2 * r2 - Br * Br);*/

        /*auto Bx2 = f[0] * f[0];
        auto By2 = f[1] * f[1];
        auto Bz2 = f[2] * f[2];
        auto rx2 = r[0] * r[0];
        auto ry2 = r[1] * r[1];
        auto rz2 = r[2] * r[2];

        auto squares = Bx2 * ry2 + Bx2 * rz2 + By2 * rx2 + By2 * rz2 + Bz2 * rx2 + Bz2 * ry2;
        auto overlaps = -2.0 * (f[0] * f[1] * r[0] * r[1] + f[0] * f[2] * r[0] * r[2] + f[1] * f[2] * r[1] * r[2]);
        RankZeroOperator MF_e = 0.125 * (squares + overlaps);*/

        /*MagneticVectorPotentialOperator A(f, r);

        RankZeroOperator MF_e = 0.5 * (A[0] * A[0] + A[1] * A[1] + A[2] * A[2]);*/

        auto Byz = (f[1] * r[2] - f[2] * r[1]);
        auto Bxz = (f[2] * r[0] - f[0] * r[2]);
        auto Bxy = (f[0] * r[1] - f[1] * r[0]);

        RankZeroOperator MF_e = 0.125 * (Bxy * Bxy + Bxz * Bxz + Byz * Byz);

        // Invoke operator= to assign *this operator
        RankZeroOperator &HMF = (*this);
        HMF = MF_e;
        HMF.name() = "B . r";
    }
};

} // namespace mrchem
