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
#include "H_B_dip.h"

/** @class MagneticFieldOperatorImag
 *
 * @brief Imaginary component of the external magnetic field operator
 *
 * An external magnetic field interacts with the molecular dipole
 * moment. The operator is simply implemented as a scalar product of
 * the dipole moment operator with the magnetic field vector.
 *
 */

namespace mrchem {

class MagneticFieldOperatorImag final : public RankZeroOperator {
public:
    MagneticFieldOperatorImag(const Eigen::Vector3d &f, std::shared_ptr<mrcpp::DerivativeOperator<3>> D, const mrcpp::Coord<3> &o)
            : MagneticFieldOperatorImag(std::array<double, 3>{f[0], f[1], f[2]}, D, o) {}

    MagneticFieldOperatorImag(const std::array<double, 3> &f, std::shared_ptr<mrcpp::DerivativeOperator<3>> D, const mrcpp::Coord<3> &o) {
        auto nabla = NablaOperator(D, true);
        auto r = PositionOperator(o);

        auto Lyz = (r[1] * nabla[2] - r[2] * nabla[1]);
        auto Lzx = (r[2] * nabla[0] - r[0] * nabla[2]);
        auto Lxy = (r[0] * nabla[1] - r[1] * nabla[0]);

        // -i * B * L_0
        RankZeroOperator MF_dip = -0.5 * (f[0] * Lyz + f[1] * Lzx + f[2] * Lxy);

        /*auto Lyz = (f[1] * r[2] - f[2] * r[1]);
        auto Lzx = (f[2] * r[0] - f[0] * r[2]);
        auto Lxy = (f[0] * r[1] - f[1] * r[0]);

        // -i * B * L_0
        RankZeroOperator MF_dip = -0.5 * (nabla[0] * Lyz + nabla[1] * Lzx + nabla[2] * Lxy);*/

        /*auto r = PositionOperator(o);
        auto p = MomentumOperator(D);
        MagneticVectorPotentialOperator A(f, r);

        // -i * B * L_0
        RankZeroOperator MF_dip = (A[0] * p[0] + A[1] * p[1] + A[2] * p[2]);*/
        MF_dip.setImag(true);

        // Invoke operator= to assign *this operator
        RankZeroOperator &HMF = (*this);
        HMF = MF_dip;
        HMF.name() = "B . mu_B";
    }
};

} // namespace mrchem
