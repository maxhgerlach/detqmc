/*
 * code_archive.cpp
 *
 *  Created on: Mar 25, 2014
 *      Author: max
 */

//This contains code fragments that are no longer in use in
//the main codebase, but should rather not get lost.

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//           Related to the "rescale move"              //
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

/////////////////////
//  from detsdw.h  //
/////////////////////
    const bool rescale;
    const uint32_t rescaleInterval;
    const num rescaleGrowthFactor;
    const num rescaleShrinkFactor;
    uint32_t acceptedRescales;
    uint32_t attemptedRescales;

    //Try a timeslice-global move, where all the phi-fields of a timeslice are multiplied
    //by a common factor.
    void attemptTimesliceRescaleMove(uint32_t timeslice, num factor);
    num deltaSPhiTimesliceRescale(uint32_t timeslice, num factor);

    	ar & acceptedRescales & attemptedRescales;

///////////////////////
//  from detsdw.cpp  //
///////////////////////
, "rescale"

        rescale(pars.rescale), rescaleInterval(pars.rescaleInterval),
        rescaleGrowthFactor(pars.rescaleGrowthFactor), rescaleShrinkFactor(pars.rescaleShrinkFactor),
        acceptedRescales(0), attemptedRescales(0),

    META_INSERT(rescale);
    if (rescale) {
        META_INSERT(rescaleInterval);
        META_INSERT(rescaleGrowthFactor);
        META_INSERT(rescaleShrinkFactor);
    }

    if (rescale) {
        if (performedSweeps % rescaleInterval == 0) {
            num rnd = rng.rand01();
            if (rnd <= 0.5) {
                attemptTimesliceRescaleMove(timeslice, rescaleGrowthFactor);
            } else {
                //attemptGlobalRescaleMove(timeslice, rescaleShrinkFactor);
                attemptTimesliceRescaleMove(timeslice, 1.0 / rescaleGrowthFactor);
            }
        }
    }

template<bool TD, CheckerboardMethod CB>
inline void DetSDW<TD,CB>::attemptTimesliceRescaleMove(uint32_t timeslice, num factor) {
    timing.start("sdw-attemptTimesliceRescaleMove");

    //see hand-written notes and Ipython notebook sdw-rescale-move to understand these formulas

    //original fields
    //TODO: unnecessary copies
    const VecNum a  {phi2.col(timeslice)};
    const VecCpx b  {phi0.col(timeslice), -phi1.col(timeslice)};
    const VecCpx bc {phi0.col(timeslice), +phi1.col(timeslice)};
    const VecNum x  {phiSinh.col(timeslice)};
    const VecNum c  {phiCosh.col(timeslice)};

    // //DEBUG
    // VecNum oldphi0 = phi0.col(timeslice);
    // VecNum oldphi1 = phi1.col(timeslice);
    // VecNum oldphi2 = phi2.col(timeslice);
    // debugSaveMatrix(oldphi0, "old_phi0");
    // debugSaveMatrix(oldphi1, "old_phi1");
    // debugSaveMatrix(oldphi2, "old_phi2");

    //rescaled fields
    const VecNum rphi0 {factor * phi0.col(timeslice)};
    const VecNum rphi1 {factor * phi1.col(timeslice)};
    const VecNum rphi2 {factor * phi2.col(timeslice)};
    const VecNum ra  {rphi2};
    const VecCpx rb  {rphi0, -rphi1};
    const VecCpx rbc {rphi0, +rphi1};
    using arma::pow; using arma::sqrt; using arma::sinh; using arma::cosh;
    const VecNum rnorm { sqrt(pow(rphi0,2) + pow(rphi1,2) + pow(rphi2,2)) };
    const VecNum rx    { sinh(dtau * rnorm) / rnorm };
    const VecNum rc    { cosh(dtau * rnorm) };

    //DEBUG
    // debugSaveMatrix(rphi0, "new_phi0");
    // debugSaveMatrix(rphi1, "new_phi1");
    // debugSaveMatrix(rphi2, "new_phi2");

    // 1) Calculate Delta = exp(-dtau V(a',b',c'))*exp(+dtau V(a,b,c)) - 1
    // Delta is setup by 4x4 blocks of size NxN, each being diagonal.
    // 4 blocks are zero, apart from that there are 5 different blocks:
    const VecNum delta_a  { rc % a % x - ra % rx % c };
    const VecNum delta_ma { -delta_a };
    const VecNum delta_c  { rc % c - ra % rx % a % x - rx % arma::real(rb % bc) % x - arma::ones<VecNum>(N) };    //Note: rb % bc will result in a purely real result
    // the block diagonals that are complex are stored with real and imaginary parts
    // separated:
    const VecNum delta_b_r  { rc % arma::real(b) % x - arma::real(rb) % rx % c };
    const VecNum delta_b_i  { rc % arma::imag(b) % x - arma::imag(rb) % rx % c };
    const VecNum delta_bc_r { delta_b_r };
    const VecNum delta_bc_i { -delta_b_i };

    //DEBUG
//    debugSaveMatrix(delta_a, "delta_a");
//    debugSaveMatrix(delta_ma, "delta_ma");
//    debugSaveMatrix(delta_c, "delta_c");
//    debugSaveMatrix(delta_b_r , "delta_b_r");
//    debugSaveMatrix(delta_b_i , "delta_b_i");
//    debugSaveMatrix(delta_bc_r, "delta_bc_r");
//    debugSaveMatrix(delta_bc_i, "delta_bc_i");

    // real part of matrix represented by 4x4 array of pointers to our vectors
    using std::array; using std::cref;
    array< array<const VecNum*, 4>, 4> delta_r;
    delta_r[0][0] = &delta_c;
    delta_r[0][1] = 0;
    delta_r[0][2] = &delta_a;
    delta_r[0][3] = &delta_b_r;
    delta_r[1][0] = 0;
    delta_r[1][1] = &delta_c;
    delta_r[1][2] = &delta_bc_r;
    delta_r[1][3] = &delta_ma;
    delta_r[2][0] = &delta_a;
    delta_r[2][1] = &delta_b_r;
    delta_r[2][2] = &delta_c;
    delta_r[2][3] = 0;
    delta_r[3][0] = &delta_bc_r;
    delta_r[3][1] = &delta_ma;
    delta_r[3][2] = 0;
    delta_r[3][3] = &delta_c;

    // imaginary part of matrix: only the antidiagonal blocks
    array< array<const VecNum*, 4>, 4> delta_i {{}};        //init with nulls
    delta_i[0][3] = &delta_b_i;
    delta_i[1][2] = &delta_bc_i;
    delta_i[2][1] = &delta_b_i;
    delta_i[3][0] = &delta_bc_i;

    // 2) Compute the matrix M = I + Delta * (I - G(timeslice))
    MatCpx oneMinusG { arma::eye(4*N,4*N) - g };
    //DEBUG
    // for (uint32_t r = 0; r < 4; ++r) {
    //     for (uint32_t c = 0; c < 4; ++c) {
    //         const VecNum* ptr = delta_r[r][c];
    //         std::string basename = "delta_r"+numToString(r)+"_c"+numToString(c);
    //         if (ptr) {
    //             debugSaveMatrix(*ptr, basename);
    //         } else {
    //             debugSaveMatrix(VecNum(arma::zeros<VecNum>(N)), basename);
    //         }
    //         const VecNum* ptr2 = delta_i[r][c];
    //         std::string basename2 = "delta_i"+numToString(r)+"_c"+numToString(c);
    //         if (ptr2) {
    //             debugSaveMatrix(*ptr2, basename2);
    //         } else {
    //             debugSaveMatrix(VecNum(arma::zeros<VecNum>(N)), basename2);
    //         }
    //     }
    // }
    MatCpx M { arma::eye(4*N,4*N), arma::zeros(4*N,4*N) };
#define block(matrix, row, col) matrix.submat((row) * N, (col) * N, ((row) + 1) * N - 1, ((col) + 1) * N - 1)
    //real parts
    for (uint32_t row = 0; row < 4; ++row) {
        for (uint32_t col = 0; col < 4; ++col) {
            //skip the zero blocks of delta:
            uint32_t skip_i;
            switch (row) {
            case 0: skip_i = 1; break;
            case 1: skip_i = 0; break;
            case 2: skip_i = 3; break;
            case 3: skip_i = 2; break;
            }
            uint32_t start_i = 0;
            for (uint32_t i = start_i; i < 4; ++i) {
                if (i == skip_i) continue;
                block(M,row,col) += arma::diagmat(*(delta_r[row][i])) *
                        block(oneMinusG, i, col);
            }
        }
    }
    //imaginary parts
    for (uint32_t row = 0; row < 4; ++row) {
        VecCpx temp { arma::zeros<VecNum>(N), *(delta_i[row][3 - row]) };   //antidiagonal: col = 3 - row
        for (uint32_t i = 0; i < 4; ++i) {
            block(M,row,i) += arma::diagmat(temp) * block(oneMinusG, 3 - row, i);
        }
    }
#undef block

    // //DEBUG
    // debugSaveMatrix(MatNum(arma::real(g)), "gslice_old_real");
    // debugSaveMatrix(MatNum(arma::imag(g)), "gslice_old_imag");
    // debugSaveMatrixCpx(M, "M");
    // MatCpx g_new = trans(solve(trans(M), trans(g)));
    // debugSaveMatrix(MatNum(arma::real(g_new)), "gslice_new_real");
    // debugSaveMatrix(MatNum(arma::imag(g_new)), "gslice_new_imag");
    // //END_DEBUG

    // 3) Compute probability of accepting the global rescale move
    num probFermion = arma::det(M).real();
    num probBoson = std::exp(-deltaSPhiTimesliceRescale(timeslice, factor));
    num prob = probFermion * probBoson;

    // //DEBUG check probBoson
    // num sphi_old = phiAction();

    // //DEBUG info
    // std::cout << "Rescale factor " << factor << " -> probFermion = " << probFermion
    //           << " \tprobBoson = " << probBoson << '\n';

//    // DEBUG-CHECK
//    VecNum debug_phi0_before = phi0.col(timeslice);
//    VecNum debug_phi1_before = phi1.col(timeslice);
//    VecNum debug_phi2_before = phi2.col(timeslice);
//    num sphi_old = phiAction();
//    phi0.col(timeslice) = rphi0;
//    phi1.col(timeslice) = rphi1;
//    phi2.col(timeslice) = rphi2;
//    num sphi_new = phiAction();
//    phi0.col(timeslice) = debug_phi0_before;
//    phi1.col(timeslice) = debug_phi1_before;
//    phi2.col(timeslice) = debug_phi2_before;
//    num prob_check = std::exp(-(sphi_new - sphi_old));
//    assert((prob_check - probBoson) / probBoson < 1E-10);
//
//    MatCpx debug_Delta(4*N,4*N);
//    debug_Delta = computePotentialExponential(-1, rphi0, rphi1, rphi2)
//            * computePotentialExponential(+1, debug_phi0_before, debug_phi1_before, debug_phi2_before)
//            - arma::eye<MatCpx>(4*N, 4*N);
//    MatCpx debug_M(4*N,4*N);
//    debug_M = arma::eye<MatCpx>(4*N, 4*N) + debug_Delta *
//            (arma::eye<MatCpx>(4*N, 4*N) - g);
//    num debug_det = arma::det(debug_M).real();
//    assert(debug_det > 0);
//    assert(probFermion > 0);
//    assert((debug_det - probFermion) / probFermion < 1E-10);
//    // END-DEBUG-CHECK

//    // DEBUG - count shrink/Grow
//    static int countShrink = 0;
//    static int countGrow = 0;
//    // ENd-DEBUG - count shrink/Grow

    if (prob > 1.0 or rng.rand01() < prob) {
//        // DEBUG - count shrink/Grow
//        if (factor < 1.0) {
//            ++countShrink;
//        } else if (factor > 1.0) {
//            ++countGrow;
//        }
//        // END-DEBUG - count shrink/Grow

        // //DEBUG info
        // std::cout << "Accepted!" << std::endl;

        //count accepted update
        ++acceptedRescales;

        //update phi-fields and dependent quantities
        phi0.col(timeslice) = rphi0;
        phi1.col(timeslice) = rphi1;
        phi2.col(timeslice) = rphi2;
        phiCosh.col(timeslice) = rc;
        phiSinh.col(timeslice) = rx;

        // //DEBUG check probBoson
        // num sphi_new = phiAction();
        // num delta_sphi = sphi_new - sphi_old;
        // num probCheck = std::exp(-delta_sphi);
        // std::cout << "Check probBoson = " << probCheck << std::endl;

        using arma::trans; using arma::solve;

        //update Green function
        g = trans(solve(trans(M), trans(g)));
        //TODO: the three transpositions here bug me

        // //DEBUG
        // std::cout << MatNum(arma::abs(g - g_new)).max() << std::endl;
        // //END DEBUG
    } else {
        // //DEBUG info
        // std::cout << "Rejected!" << std::endl;

        // //DEBUG check probBoson
        // phi0.col(timeslice) = rphi0;
        // phi1.col(timeslice) = rphi1;
        // phi2.col(timeslice) = rphi2;

        // //DEBUG check probBoson
        // num sphi_new = phiAction();
        // num delta_sphi = sphi_new - sphi_old;
        // num probCheck = std::exp(-delta_sphi);

        // phi0.col(timeslice) = oldphi0;
        // phi1.col(timeslice) = oldphi1;
        // phi2.col(timeslice) = oldphi2;

        // std::cout << "Check probBoson = " << probCheck << std::endl;
    }

    // //DEBUG
    // if (rng.rand01() < 0.15) {
    //     exit(1);
    // }
    // //END DEBUG

    ++attemptedRescales;

//    // DEBUG - count shrink/Grow
//    std::cout << "Shrink: " << countShrink << " Grow: " << countGrow << "\n";
//    // END-DEBUG - count shrink/Grow

    timing.stop("sdw-attemptTimesliceRescaleMove");
}


template<bool TD, CheckerboardMethod CB>
num DetSDW<TD,CB>::deltaSPhiTimesliceRescale(uint32_t timeslice, num factor) {
    using std::pow;
    num delta1 = 0;
    for (uint32_t site_i = 0; site_i < N; ++site_i) {
        uint32_t neighSites[] = {spaceNeigh(XPLUS, site_i), spaceNeigh(YPLUS, site_i)};   //for Icpc
        for (uint32_t site_j : neighSites) {
            delta1 += pow(phi0(site_i, timeslice) - phi0(site_j, timeslice), 2)
                    + pow(phi1(site_i, timeslice) - phi1(site_j, timeslice), 2)
                    + pow(phi2(site_i, timeslice) - phi2(site_j, timeslice), 2);
        }
    }

    num delta2 = 0;
    for (uint32_t site_i = 0; site_i < N; ++site_i) {
        delta2 += pow(phi0(site_i, timeslice), 2)
                + pow(phi1(site_i, timeslice), 2)
                + pow(phi2(site_i, timeslice), 2);
    }

    num delta3 = 0;
    for (uint32_t site_i = 0; site_i < N; ++site_i) {
        delta3 += pow(pow(phi0(site_i, timeslice), 2)
                    + pow(phi1(site_i, timeslice), 2)
                    + pow(phi2(site_i, timeslice), 2), 2);
    }

    num delta4 = 0;
    uint32_t timeslicePlus = timeNeigh(ChainDir::PLUS, timeslice);
    uint32_t timesliceMinus = timeNeigh(ChainDir::MINUS, timeslice);
    for (uint32_t site_i = 0; site_i < N; ++site_i) {
        delta4 += (pow(factor,2) - 1.0) * (
                      pow(phi0(site_i, timeslice), 2)
                    + pow(phi1(site_i, timeslice), 2)
                    + pow(phi2(site_i, timeslice), 2)
                );
        delta4 -= (factor - 1.0) * (
                      phi0(site_i, timeslice) * (phi0(site_i, timesliceMinus) + phi0(site_i, timeslicePlus))
                    + phi1(site_i, timeslice) * (phi1(site_i, timesliceMinus) + phi1(site_i, timeslicePlus))
                    + phi2(site_i, timeslice) * (phi2(site_i, timesliceMinus) + phi2(site_i, timeslicePlus))
                );
    }

    num delta = (dtau/2.0) * (pow(factor,2) - 1.0) * delta1
              + (dtau*r/2.0) * (pow(factor,2) - 1.0) * delta2
              + (dtau*u/4.0) * (pow(factor,4) - 1.0) * delta3
              + (1.0/(c*dtau)) * delta4;

    return delta;
}

    if (rescale) {
        num ratio = num(acceptedRescales) / num(attemptedRescales);
        std::cout << "Timeslice rescale move acceptance ratio = " << ratio
                  << std::endl;
    }



/////////////////////////
//  from parameters.h  //
/////////////////////////
    bool rescale;	//perform timeslice rescale move?
    uint32_t rescaleInterval;		// attempt global rescale move every # sweeps
    num rescaleGrowthFactor;		// factor by which to size up the fields
    num rescaleShrinkFactor;		// factor by which to size down the fields

            rescale(), rescaleInterval(), rescaleGrowthFactor(), rescaleShrinkFactor(),

           & rescale & rescaleInterval & rescaleGrowthFactor & rescaleShrinkFactor

////////////////////////
//  from maindet.cpp  //
////////////////////////
            ("rescale", po::value<bool>(&modelpar.rescale)->default_value(false), "SDW: perform timeslice rescale move?")
            ("rescaleInterval", po::value<uint32_t>(&modelpar.rescaleInterval)->default_value(100), "attempt timeslice rescale move every # sweeps")
            ("rescaleGrowthFactor", po::value<num>(&modelpar.rescaleGrowthFactor)->default_value(1.05), "factor by which to attempt to grow the fields")
            ("rescaleShrinkFactor", po::value<num>(&modelpar.rescaleShrinkFactor)->default_value(0.95), "factor by which to attempt to shrink the fields")
