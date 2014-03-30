/*
 * code_archive_cb.cpp
 *
 *  Created on: Mar 30, 2014
 *      Author: max
 */

//This contains code fragments that are no longer in use in
//the main codebase, but should rather not get lost.

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//           Related to cb methods apart assaad_berg    //
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////


//////////////////
// maindet.cpp
//////////////////
            ("checkerboardMethod", po::value<std::string>(&modelpar.checkerboardMethod)->default_value("assaad_berg"), "method to use for the checkerboard decomposition: santos, assaad or assaad_berg")

//////////////////
// parameters.h
//////////////////
    std::string checkerboardMethod;		//SDW if checkerboard: "santos" or "assaad" or "assaad_berg"
checkerboardMethod(),
   & checkerboardMethod


//////////////////
// detsdw.h
//////////////////


enum CheckerboardMethod {
	CB_NONE,				//regular, dense matrix products
	CB_SANTOS,				//checkerboard, four break-ups in e^{K_*} as described in: R. R. dos Santos, Braz. J. Phys 33, 36 (2003).
	CB_ASSAAD,				//checkerboard, two break-ups in e^{K_*} as described in: F. F. Assaad, in Quantum Simulations Complex Many-Body Syst. From Theory to Algorithms, edited by J. Grotendorst, D. Marx, and A. Muramatsu (FZ-J�lich, J�lich, Germany, 2002).
	CB_ASSAAD_BERG,			//checkerboard, two break-ups, making sure all multiplications are symmetric, as described by Erez Berg
};
std::string cbmToString(CheckerboardMethod cbm);


    const std::string checkerboardMethod;

    template <class Matrix>
    MatCpx cbLMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_SANTOS>,
    							  const Matrix& A, Band band, int sign, bool invertedCbOrder = false);
    template <class Matrix>
    MatCpx cbRMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_SANTOS>,
    							  const Matrix& A, Band band, int sign, bool invertedCbOrder = false);
    template <class Matrix>
    MatCpx cbLMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_ASSAAD>,
    							  const Matrix& A, Band band, int sign, bool invertedCbOrder = false);
    template <class Matrix>
    MatCpx cbRMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_ASSAAD>,
    							  const Matrix& A, Band band, int sign, bool invertedCbOrder = false);


    template<class Matrix>
    void cb_santos_applyBondFactorsLeft(Matrix& result, NeighDir neigh, uint32_t subgroup, num ch, num sh);
    template<class Matrix>
    void cb_santos_applyBondFactorsRight(Matrix& result, NeighDir neigh, uint32_t subgroup, num ch, num sh);



////////////////////
// detsdw.cpp
////////////////////

std::string cbmToString(CheckerboardMethod cbm) {
    switch (cbm) {
    case CB_NONE: return "NONE";
    case CB_SANTOS: return "santos";
    case CB_ASSAAD: return "assaad";
    case CB_ASSAAD_BERG: return "assaad_berg";
    default: return "INVALID_CHECKERBOARD_METHOD";
    }
}

    CheckerboardMethod cbm = CB_NONE;
    if (pars.checkerboard) {
        if (pars.checkerboardMethod == "santos") {
            cbm = CB_SANTOS;
        } else if (pars.checkerboardMethod == "assaad") {
            cbm = CB_ASSAAD;
        } else if (pars.checkerboardMethod == "assaad_berg") {
            cbm = CB_ASSAAD_BERG;
        } else {
            throw ParameterWrong("checkerboardMethod", pars.checkerboardMethod);
        }
    }


        if (cbm == CB_SANTOS) {
        return std::unique_ptr<DetModel>(new DetSDW<false,CB_SANTOS>(rng, pars));
    } else
    if (cbm == CB_ASSAAD) {
        return std::unique_ptr<DetModel>(new DetSDW<false,CB_ASSAAD>(rng, pars));
    } else

        checkerboardMethod(pars.checkerboardMethod),



    if (CB) {
        meta["checkerboardMethod"] = checkerboardMethod;
    }




//neigh == XNEIGH:
//   subgroup == 0:  bonds (2*i_x, i_y)--(2*i_x + 1, i_y)
//   subgroup == 1:  bonds (2*i_x + 1, i_y)--(2*i_x + 2, i_y)
//neigh == YNEIGH:
//   subgroup == 0:  bonds (i_x, 2*i_y)--(i_x, 2*i_y + 1)
//   subgroup == 1:  bonds (i_x, 2*i_y + 1)--(i_x, 2*i_y + 2)
template<bool TD, CheckerboardMethod CB>
template<class Matrix>
void DetSDW<TD,CB>::cb_santos_applyBondFactorsLeft(Matrix& result, const NeighDir neigh, const uint32_t subgroup, const num ch, const num sh) {
    assert(subgroup == 0 or subgroup == 1);
    assert(neigh == XPLUS or neigh == YPLUS);
    arma::Row<cpx> new_row_i(N);
    for (uint32_t i1 = subgroup; i1 < L; i1 += 2) {
        for (uint32_t i2 = 0; i2 < L; ++i2) {
            uint32_t i;
            switch (neigh) {
            case XPLUS:
                i = this->coordsToSite(i1, i2);
                break;
            case YPLUS:
                i = this->coordsToSite(i2, i1);
                break;
            default: //should not be reached
                break;
            }
            uint32_t j = spaceNeigh(neigh, i);
            //change rows i and j of result
            num b_sh = sh;
            if ((bc == APBC_X or bc == APBC_XY) and neigh == XPLUS and i1 == L-1) {
                //crossed antiperiodic boundary
                b_sh *= -1;
            }
            else if ((bc == APBC_Y or bc == APBC_XY) and neigh == YPLUS and i1 == L-1) {
                //crossed antiperiodic boundary
                b_sh *= -1;
            }
            new_row_i     = ch * result.row(i) + b_sh * result.row(j);
            result.row(j) = b_sh * result.row(i) + ch * result.row(j);
            result.row(i) = new_row_i;
        }
    }
}



// with sign = +/- 1, band = XBAND|YBAND: set R := E^(sign * dtau * K_band) * A
// using the method described in R. R. dos Santos, Braz. J. Phys 33, 36 (2003).
template<bool TD, CheckerboardMethod CB>
template<class Matrix> inline
MatCpx DetSDW<TD,CB>::cbLMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_SANTOS>,
                                             const Matrix& A, Band band, int sign, bool invertedCbOrder) {
    MatCpx result = A;      //can't avoid this copy

    if (not invertedCbOrder) {
        cb_santos_applyBondFactorsLeft(result, XPLUS, 0, coshHopHor[band], sign * sinhHopHor[band]);
        cb_santos_applyBondFactorsLeft(result, YPLUS, 0, coshHopVer[band], sign * sinhHopVer[band]);
        cb_santos_applyBondFactorsLeft(result, XPLUS, 1, coshHopHor[band], sign * sinhHopHor[band]);
        cb_santos_applyBondFactorsLeft(result, YPLUS, 1, coshHopVer[band], sign * sinhHopVer[band]);
    } else {
        cb_santos_applyBondFactorsLeft(result, YPLUS, 1, coshHopVer[band], sign * sinhHopVer[band]);
        cb_santos_applyBondFactorsLeft(result, XPLUS, 1, coshHopHor[band], sign * sinhHopHor[band]);
        cb_santos_applyBondFactorsLeft(result, YPLUS, 0, coshHopVer[band], sign * sinhHopVer[band]);
        cb_santos_applyBondFactorsLeft(result, XPLUS, 0, coshHopHor[band], sign * sinhHopHor[band]);
    }
    return result;
}



// with sign = +/- 1, band = XBAND|YBAND: set R := E^(sign * dtau * K_band) * A
// using the method described in F. F. Assaad, in Quantum Simulations Complex Many-Body Syst. From Theory to Algorithms, edited by J. Grotendorst, D. Marx, and A. Muramatsu (FZ-Jülich, Jülich, Germany, 2002).
template<bool TD, CheckerboardMethod CB>
template<class Matrix> inline
MatCpx DetSDW<TD,CB>::cbLMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_ASSAAD>,
                                             const Matrix& A, Band band, int sign, bool invertedCbOrder) {
    MatCpx result = A;      //can't avoid this copy

    if (not invertedCbOrder) {
        cb_assaad_applyBondFactorsLeft(result, 0, coshHopHor[band], sign * sinhHopHor[band], coshHopVer[band], sign * sinhHopVer[band]);
        cb_assaad_applyBondFactorsLeft(result, 1, coshHopHor[band], sign * sinhHopHor[band], coshHopVer[band], sign * sinhHopVer[band]);
    } else {
        cb_assaad_applyBondFactorsLeft(result, 1, coshHopHor[band], sign * sinhHopHor[band], coshHopVer[band], sign * sinhHopVer[band]);
        cb_assaad_applyBondFactorsLeft(result, 0, coshHopHor[band], sign * sinhHopHor[band], coshHopVer[band], sign * sinhHopVer[band]);
    }
    return result;
}



//neigh == XNEIGH:
//   subgroup == 0:  bonds (2*i_x, i_y)--(2*i_x + 1, i_y)
//   subgroup == 1:  bonds (2*i_x + 1, i_y)--(2*i_x + 2, i_y)
//neigh == YNEIGH:
//   subgroup == 0:  bonds (i_x, 2*i_y)--(i_x, 2*i_y + 1)
//   subgroup == 1:  bonds (i_x, 2*i_y + 1)--(i_x, 2*i_y + 2)
template<bool TD, CheckerboardMethod CB>
template<class Matrix>
void DetSDW<TD,CB>::cb_santos_applyBondFactorsRight(Matrix& result, const NeighDir neigh, const uint32_t subgroup, const num ch, const num sh) {
    assert(subgroup == 0 or subgroup == 1);
    assert(neigh == XPLUS or neigh == YPLUS);
    arma::Col<cpx> new_col_i(N);
    for (uint32_t i1 = subgroup; i1 < L; i1 += 2) {
        for (uint32_t i2 = 0; i2 < L; ++i2) {
            uint32_t i;
            switch (neigh) {
            case XPLUS:
                i = this->coordsToSite(i1, i2);
                break;
            case YPLUS:
                i = this->coordsToSite(i2, i1);
                break;
            default: //should not be reached
                break;
            }
            uint32_t j = spaceNeigh(neigh, i);
            //change columns i and j of result
            num b_sh = sh;
            if ((bc == APBC_X or bc == APBC_XY) and neigh == XPLUS and i1 == L-1) {
                //crossed antiperiodic boundary
                b_sh *= -1;
            }
            else if ((bc == APBC_Y or bc == APBC_XY) and neigh == YPLUS and i1 == L-1) {
                //crossed antiperiodic boundary
                b_sh *= -1;
            }
            new_col_i     = ch * result.col(i) + b_sh * result.col(j);
            result.col(j) = b_sh * result.col(i) + ch * result.col(j);
            result.col(i) = new_col_i;
        }
    }
}



// with sign = +/- 1, band = XBAND|YBAND: return A * E^(sign * dtau * K_band)
// using the method described in R. R. dos Santos, Braz. J. Phys 33, 36 (2003).
template<bool TD, CheckerboardMethod CB>
template <class Matrix> inline
MatCpx DetSDW<TD,CB>::cbRMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_SANTOS>,
                                             const Matrix& A, Band band, int sign, bool invertedCbOrder) {
    MatCpx result = A;      //can't avoid this copy

    //order reversed wrt cbLMultHoppingExp
    if (not invertedCbOrder) {
        cb_santos_applyBondFactorsRight(result, YPLUS, 1, coshHopVer[band], sign * sinhHopVer[band]);
        cb_santos_applyBondFactorsRight(result, XPLUS, 1, coshHopHor[band], sign * sinhHopHor[band]);
        cb_santos_applyBondFactorsRight(result, YPLUS, 0, coshHopVer[band], sign * sinhHopVer[band]);
        cb_santos_applyBondFactorsRight(result, XPLUS, 0, coshHopHor[band], sign * sinhHopHor[band]);
    } else {
        cb_santos_applyBondFactorsRight(result, XPLUS, 0, coshHopHor[band], sign * sinhHopHor[band]);
        cb_santos_applyBondFactorsRight(result, YPLUS, 0, coshHopVer[band], sign * sinhHopVer[band]);
        cb_santos_applyBondFactorsRight(result, XPLUS, 1, coshHopHor[band], sign * sinhHopHor[band]);
        cb_santos_applyBondFactorsRight(result, YPLUS, 1, coshHopVer[band], sign * sinhHopVer[band]);
    }

    return result;
}



template<bool TD, CheckerboardMethod CB>
template<class Matrix> inline
MatCpx DetSDW<TD,CB>::cbRMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_ASSAAD>,
                                             const Matrix& A, Band band, int sign, bool invertedCbOrder) {
    MatCpx result = A;      //can't avoid this copy

    //order reversed wrt cbLMultHoppingExp
    if (not invertedCbOrder) {
        cb_assaad_applyBondFactorsRight(result, 1, coshHopHor[band], sign * sinhHopHor[band], coshHopVer[band], sign * sinhHopVer[band]);
        cb_assaad_applyBondFactorsRight(result, 0, coshHopHor[band], sign * sinhHopHor[band], coshHopVer[band], sign * sinhHopVer[band]);
    } else {
        cb_assaad_applyBondFactorsRight(result, 1, coshHopHor[band], sign * sinhHopHor[band], coshHopVer[band], sign * sinhHopVer[band]);
        cb_assaad_applyBondFactorsRight(result, 0, coshHopHor[band], sign * sinhHopHor[band], coshHopVer[band], sign * sinhHopVer[band]);
    }
    return result;
}



//shiftGreenSymmetric:
    // unclear to me: why do I need to explicitly qualify 'this->' in the lambdas?
    else if (CB == CB_SANTOS) {
        return shiftGreenSymmetric_impl(
                        //rightMultiply
            // output and input are NxN blocks of a complex matrix
            // this effectively multiplies e^{+ dtau K^band_b / 2} e^{+ dtau K^band_a / 2}
            // to the right of input and stores the result in output
            [this](SubMatCpx output, SubMatCpx input, Band band) -> void {
                output = input;            //copy
                this->cb_santos_applyBondFactorsRight(output, YPLUS, 1, coshHopVerHalf[band], +sinhHopVerHalf[band]);
                this->cb_santos_applyBondFactorsRight(output, XPLUS, 1, coshHopHorHalf[band], +sinhHopHorHalf[band]);
                this->cb_santos_applyBondFactorsRight(output, YPLUS, 0, coshHopVerHalf[band], +sinhHopVerHalf[band]);
                this->cb_santos_applyBondFactorsRight(output, XPLUS, 0, coshHopHorHalf[band], +sinhHopHorHalf[band]);
            },
            //leftMultiply
            // output and input are NxN blocks of a complex matrix
            // this effectively multiplies e^{- dtau K^band_a / 2} e^{- dtau K^band_b / 2}
            // to the left of input and stores the result in output
            [this](SubMatCpx output, SubMatCpx input, Band band) -> void {
                output = input;            //copy
                this->cb_santos_applyBondFactorsLeft(output, XPLUS, 0, coshHopHorHalf[band], -sinhHopHorHalf[band]);
                this->cb_santos_applyBondFactorsLeft(output, YPLUS, 0, coshHopVerHalf[band], -sinhHopVerHalf[band]);
                this->cb_santos_applyBondFactorsLeft(output, XPLUS, 1, coshHopHorHalf[band], -sinhHopHorHalf[band]);
                this->cb_santos_applyBondFactorsLeft(output, YPLUS, 1, coshHopVerHalf[band], -sinhHopVerHalf[band]);
            }
        );
    }
    else if (CB == CB_ASSAAD) {
        return shiftGreenSymmetric_impl(
                        //rightMultiply
            // output and input are NxN blocks of a complex matrix
            // this effectively multiplies [Input] * e^{+ dtau K^band_b / 2} e^{+ dtau K^band_a / 2}
            // to the right of input and stores the result in output
            [this](SubMatCpx output, SubMatCpx input, Band band) -> void {
                output = input;      //copy
                this->cb_assaad_applyBondFactorsRight(output, 1, coshHopHorHalf[band], +sinhHopHorHalf[band],
                                                                 coshHopVerHalf[band], +sinhHopVerHalf[band]);
                this->cb_assaad_applyBondFactorsRight(output, 0, coshHopHorHalf[band], +sinhHopHorHalf[band],
                                                                       coshHopVerHalf[band], +sinhHopVerHalf[band]);
            },
            //leftMultiply
            // output and input are NxN blocks of a complex matrix
            // this effectively multiplies e^{- dtau K^band_b / 2} e^{- dtau K^band_a / 2} * [Input]
            // to the left of input and stores the result in output
            [this](SubMatCpx output, SubMatCpx input, Band band) -> void {
                output = input;      //copy
                this->cb_assaad_applyBondFactorsLeft(output, 0, coshHopHorHalf[band], -sinhHopHorHalf[band],
                                                                       coshHopVerHalf[band], -sinhHopVerHalf[band]);
                this->cb_assaad_applyBondFactorsLeft(output, 1, coshHopHorHalf[band], -sinhHopHorHalf[band],
                                                                coshHopVerHalf[band], -sinhHopVerHalf[band]);
            }
        );
    }




//explicit template instantiations:
//template class DetSDW<true,CB_NONE>;
template class DetSDW<false,CB_NONE>;
//template class DetSDW<true,CB_SANTOS>;
template class DetSDW<false,CB_SANTOS>;
//template class DetSDW<true,CB_ASSAAD>;
template class DetSDW<false,CB_ASSAAD>;
//template class DetSDW<true,CB_ASSAAD_BERG>;
template class DetSDW<false,CB_ASSAAD_BERG>;





////////////////////
// detqmc.h
////////////////////
        // if (DetSDW<true, CB_ASSAAD>* p = dynamic_cast<DetSDW<true, CB_ASSAAD>*>(replica.get())) {
        //     p->loadContents(SerializeContentsKey(), ar);
        // } else

        if (DetSDW<false, CB_ASSAAD>* p = dynamic_cast<DetSDW<false, CB_ASSAAD>*>(replica.get())) {
            p->loadContents(SerializeContentsKey(), ar);
        } else


        // if (DetSDW<true, CB_ASSAAD>* p = dynamic_cast<DetSDW<true, CB_ASSAAD>*>(replica.get())) {
        //     p->saveContents(SerializeContentsKey(), ar);
        // } else

        if (DetSDW<false, CB_ASSAAD>* p = dynamic_cast<DetSDW<false, CB_ASSAAD>*>(replica.get())) {
            p->saveContents(SerializeContentsKey(), ar);
        } else
