//Code by Max Henner Gerlach, 2010--2012,
//used for the diploma thesis "Directional Ordering in the Classical Compass Model in Two and Three Dimensions"
//contact: maxgerlach@gmail.com

/*
 * multireweighthistosptjk.cpp
 *
 *  Created on: Jul 25, 2011
 *      Author: gerlach
 */

// generalized for SDW DQMC (2015-02-06 - )

#include "mrpt-jk.h"
#include "tools.h"
#include "statistics.h"
#include "numerics.h"

using namespace std;

MultireweightHistosPTJK::MultireweightHistosPTJK(unsigned jkTotalBlocks, std::ostream& outStream)
    : MultireweightHistosPT(outStream), blockCount(jkTotalBlocks) {
    out << "Jackknifing: " << jkTotalBlocks << " blocks" << std::endl;
}


MultireweightHistosPTJK::~MultireweightHistosPTJK() {
}

void MultireweightHistosPTJK::setUpHistogramsHelperJK() {
    lZ_bl.resize(boost::extents[blockCount][numReplicas]);
    initArray(lZ_bl, LogVal(1.0));

    H_bkm.resize(boost::extents[blockCount][numReplicas][binCount]);
    initArray(H_bkm, 0);
    g_bkm.resize(boost::extents[blockCount][numReplicas][binCount]);
    initArray(g_bkm, 0);
    H_blm.resize(boost::extents[blockCount][numReplicas][binCount]);
    initArray(H_blm, 0);
    N_bkl.resize(boost::extents[blockCount][numReplicas][numReplicas]);
    initArray(N_bkl, 0);

    H_bm.resize(boost::extents[blockCount][binCount]);
    initArray(H_bm, 0);

    lOmega_bm.resize(boost::extents[blockCount][binCount]);
    initArray(lOmega_bm, LogVal(1.0));
}

void MultireweightHistosPTJK::setUpHistogramsJK(int binCount_) {
    MultireweightHistosPT::setUpHistograms(binCount_);      //whole data set
    setUpHistogramsHelperJK();
}

void MultireweightHistosPTJK::setUpHistogramsIsingJK() {
    MultireweightHistosPT::setUpHistogramsIsing();      //whole data set
    setUpHistogramsHelperJK();
}

void MultireweightHistosPTJK::createHistogramsHelperJK() {
    out << "Creating energy histograms etc., jackknife ready. minEnergyNormalized=" << minEnergyNormalized
            << " maxEnergyNormalized=" << maxEnergyNormalized
            << " binSize=" << binSize
            << " binCount=" << binCount << endl;

    for (unsigned k = 0; k < numReplicas; ++k) {
        unsigned N_k = energyTimeSeries[k]->size();
        for (unsigned n = 0; n < N_k; ++n) {
            //get bin number, correct due to rounding in cast:
            int m = static_cast<int>(((*energyTimeSeries[k])[n] - minEnergyNormalized) / binSize);
            //get beta index
            int l = (*(betaIndexTimeSeries[k]))[n];

            //update quantities over the whole data set:
            (*m_kn[k])[n] = m;
            N_kl[k][l] += 1;
            ++H_km[k][m];
            ++H_m[m];               //!!!
            ++H_lm[l][m];           //!!!

            unsigned curBlock = (n * blockCount) / N_k;     //integer division rounds down
            for (signed b = 0; b < (signed)curBlock; ++b) {
                N_bkl[b][k][l] += 1;
                ++H_bkm[b][k][m];
                ++H_bm[b][m];       //!!!
                ++H_blm[b][l][m];   //!!!
            }
            for (signed b = curBlock + 1; b < (signed)blockCount; ++b) {
                N_bkl[b][k][l] += 1;
                ++H_bkm[b][k][m];
                ++H_bm[b][m];       //!!!
                ++H_blm[b][l][m];   //!!!
            }
        }
        out << "." << flush;
    }
    destroyAll(betaIndexTimeSeries);                    //not needed any more
    out << " done" << endl;
}

void MultireweightHistosPTJK::createHistogramsHelperDiscreteJK() {
    out << "Creating energy histograms etc., jackknife ready. minEnergyNormalized=" << minEnergyNormalized
            << " maxEnergyNormalized=" << maxEnergyNormalized
            << " binSize=" << binSize
            << " binCount=" << binCount
            << " for originally discrete bins " << endl;

    for (unsigned k = 0; k < numReplicas; ++k) {
        unsigned N_k = energyTimeSeries[k]->size();
        for (unsigned n = 0; n < N_k; ++n) {
            //get bin number, correct due to rounding in cast:
//            int m = static_cast<int>(((*energyTimeSeries[k])[n] - minEnergyNormalized) / binSize);
//            int m = int(round(((*energyTimeSeries[k])[n] - minEnergyNormalized) / binSize));
            int m = int(round(((*energyTimeSeries[k])[n] * systemN - minEnergy) / deltaU));
            //get beta index
            int l = (*(betaIndexTimeSeries[k]))[n];

            //update quantities over the whole data set:
            (*m_kn[k])[n] = m;
            N_kl[k][l] += 1;
            ++H_km[k][m];
            ++H_m[m];               //!!!
            ++H_lm[l][m];           //!!!

            unsigned curBlock = (n * blockCount) / N_k;     //integer division rounds down
            for (signed b = 0; b < (signed)curBlock; ++b) {
                N_bkl[b][k][l] += 1;
                ++H_bkm[b][k][m];
                ++H_bm[b][m];       //!!!
                ++H_blm[b][l][m];   //!!!
            }
            for (signed b = curBlock + 1; b < (signed)blockCount; ++b) {
                N_bkl[b][k][l] += 1;
                ++H_bkm[b][k][m];
                ++H_bm[b][m];       //!!!
                ++H_blm[b][l][m];   //!!!
            }
        }
        out << "." << flush;
    }
    destroyAll(betaIndexTimeSeries);                    //not needed any more
    out << " done" << endl;
}

void MultireweightHistosPTJK::createHistograms(int binCount_) {
    setUpHistogramsJK(binCount_);
    createHistogramsHelperJK();
}

void MultireweightHistosPTJK::createHistogramsIsing() {
    setUpHistogramsIsingJK();
    createHistogramsHelperDiscreteJK();
}

MultireweightHistosPT::ResultsMap* MultireweightHistosPTJK::
        directNoReweighting() {
    out << "Computing estimates from time series without any reweighting, jackknifed... "
        << flush;
    ResultsMap* results = new ResultsMap;

    //map (jk-block, beta-index) -> value
    Double2Array meanEnergy_bl(boost::extents[blockCount][numReplicas]);
    initArray(meanEnergy_bl, 0);
    Double2Array meanEnergySquared_bl(boost::extents[blockCount][numReplicas]);
    initArray(meanEnergySquared_bl, 0);
    Double2Array meanObs_bl(boost::extents[blockCount][numReplicas]);
    initArray(meanObs_bl, 0);
    Double2Array meanObsSquared_bl(boost::extents[blockCount][numReplicas]);
    initArray(meanObsSquared_bl, 0);
    Double2Array meanObsToTheFourth_bl(boost::extents[blockCount][numReplicas]);
    initArray(meanObsToTheFourth_bl, 0);

    //calculate for each replica separately:
    #pragma omp parallel for
    for (signed k = 0; k < (signed)numReplicas; ++k) {
        //arrays for this replica
        Double2Array k_meanEnergy_bl(boost::extents[blockCount][numReplicas]);
        initArray(k_meanEnergy_bl, 0);
        Double2Array k_meanEnergySquared_bl(boost::extents[blockCount][numReplicas]);
        initArray(k_meanEnergySquared_bl, 0);
        Double2Array k_meanObs_bl(boost::extents[blockCount][numReplicas]);
        initArray(k_meanObs_bl, 0);
        Double2Array k_meanObsSquared_bl(boost::extents[blockCount][numReplicas]);
        initArray(k_meanObsSquared_bl, 0);
        Double2Array k_meanObsToTheFourth_bl(boost::extents[blockCount][numReplicas]);
        initArray(k_meanObsToTheFourth_bl, 0);
        unsigned N_k = energyTimeSeries[k]->size();
        unsigned jkBlockSize = N_k / blockCount;
        //sum up:
        for (unsigned n = 0; n < N_k; ++n) {
            int bi = (*betaIndexTimeSeries[k])[n];
            double e = (*energyTimeSeries[k])[n];
            double o = (*observableTimeSeries[k])[n];
            const unsigned curBlock = n / jkBlockSize;
            for (unsigned b = 0; b < blockCount; ++b) {
                //TODO: optimize away this if: (two for loops -- should be
                //done automatically anyway?)
                if (b == curBlock) {
                    //leave out for this block
                    continue;
                } else {
                    k_meanEnergy_bl[b][bi] += e;
                    k_meanEnergySquared_bl[b][bi] += e*e;
                    k_meanObs_bl[b][bi] += o;
                    k_meanObsSquared_bl[b][bi] += o*o;
                    k_meanObsToTheFourth_bl[b][bi] += o*o*o*o;
                }
            }
        }
        double samplesPerBlock = double(N_k - jkBlockSize)
                / double(numReplicas);      //divisor: number of temperatures
        //divide to form averages:
        for (unsigned b = 0; b < blockCount; ++b) {
            for (unsigned bi = 0; bi < numReplicas; ++bi) {
                k_meanEnergy_bl[b][bi] /= samplesPerBlock;
                k_meanEnergySquared_bl[b][bi] /= samplesPerBlock;
                k_meanObs_bl[b][bi] /= samplesPerBlock;
                k_meanObsSquared_bl[b][bi] /= samplesPerBlock;
                k_meanObsToTheFourth_bl[b][bi] /= samplesPerBlock;
            }
        }
        //average over replicas -- sum up everything from all replicas
        #pragma omp critical
        {
            for (unsigned b = 0; b < blockCount; ++b) {
                for (unsigned bi = 0; bi < numReplicas; ++bi) {
                    meanEnergy_bl[b][bi] += k_meanEnergy_bl[b][bi];
                    meanEnergySquared_bl[b][bi] += k_meanEnergySquared_bl[b][bi];
                    meanObs_bl[b][bi] += k_meanObs_bl[b][bi];
                    meanObsSquared_bl[b][bi] += k_meanObsSquared_bl[b][bi];
                    meanObsToTheFourth_bl[b][bi] +=
                            k_meanObsToTheFourth_bl[b][bi];
                }
            }
        }
    }
    //average over replicas -- divide, compute quantities from higher moments
    Double2Array heatCapacity_bl(boost::extents[blockCount][numReplicas]);
    Double2Array suscObs_bl(boost::extents[blockCount][numReplicas]);
    Double2Array binderObs_bl(boost::extents[blockCount][numReplicas]);
    for (unsigned b = 0; b < blockCount; ++b) {
        for (unsigned bi = 0; bi < numReplicas; ++bi) {
            double beta = betas[bi];
            meanEnergy_bl[b][bi] /= double(numReplicas);
            meanEnergySquared_bl[b][bi] /= double(numReplicas);
            meanObs_bl[b][bi] /= double(numReplicas);
            meanObsSquared_bl[b][bi] /= double(numReplicas);
            meanObsToTheFourth_bl[b][bi] /= double(numReplicas);
            heatCapacity_bl[b][bi] = systemN * beta*beta *
                    (meanEnergySquared_bl[b][bi] - pow(meanEnergy_bl[b][bi], 2));
            suscObs_bl[b][bi] = systemN *
                    (meanObsSquared_bl[b][bi] - pow(meanObs_bl[b][bi], 2));
            binderObs_bl[b][bi] = 1.0 - (meanObsToTheFourth_bl[b][bi] /
                    (3 * pow(meanObsSquared_bl[b][bi], 2)));
        }
    }

    //combine Jackknife blocks, for each beta...
    for (unsigned bi = 0; bi < numReplicas; ++bi) {
        double beta = betas[bi];
        //sum over jk blocks --> averages computed this way have reduced bias.
        double jkSumEnergy = 0;
        double jkSumHeatCapacity = 0;
        double jkSumObs = 0;
        double jkSumSusc = 0;
        double jkSumBinder = 0;
        for (unsigned b = 0; b < blockCount; ++b) {
            jkSumEnergy += meanEnergy_bl[b][bi];
            jkSumHeatCapacity += heatCapacity_bl[b][bi];
            jkSumObs += meanObs_bl[b][bi];
            jkSumSusc += suscObs_bl[b][bi];
            jkSumBinder += binderObs_bl[b][bi];
        }
        double resEnergy = jkSumEnergy / blockCount;
        double resHeatCapacity = jkSumHeatCapacity / blockCount;
        double resObs = jkSumObs / blockCount;
        double resSusc = jkSumSusc / blockCount;
        double resBinder = jkSumBinder / blockCount;

        //error estimation
        double sqDevEnergy = 0;
        double sqDevHeatCapacity = 0;
        double sqDevObs = 0;
        double sqDevSusc = 0;
        double sqDevBinder = 0;
        for (unsigned b = 0; b < blockCount; ++b) {
            sqDevEnergy += pow(resEnergy - meanEnergy_bl[b][bi], 2);
            sqDevHeatCapacity +=
                    pow(resHeatCapacity- heatCapacity_bl[b][bi], 2);
            sqDevObs += pow(resObs - meanObs_bl[b][bi], 2);
            sqDevSusc += pow(resSusc - suscObs_bl[b][bi], 2);
            sqDevBinder += pow(resBinder - binderObs_bl[b][bi], 2);
        }
        //set values and errors:
        (*results)[beta] = ReweightingResult(resEnergy, sqrt(sqDevEnergy),
                resHeatCapacity, sqrt(sqDevHeatCapacity),
                resObs, sqrt(sqDevObs),
                resSusc, sqrt(sqDevSusc),
                resBinder, sqrt(sqDevBinder));
    }

    out << "done" << endl;
    return results;
}


void MultireweightHistosPTJK::findDensityOfStatesNonIteratively() {
    MultireweightHistosPT::findDensityOfStatesNonIteratively();

    out << "Non-iterative estimate of the density of states, jackknife blocks... " << endl;

    const double deltaU = binSize * systemN;

    for (int b = 0; b < (signed)blockCount; ++b) {
        LogVal2Array x_lm(boost::extents[numReplicas][binCount - 1]);       // <-> difference of microcanonical entropies at bins m, m+1
        Double2Array w_lm(boost::extents[numReplicas][binCount - 1]);       //weights at temperature l, energy bin m
        vector<double> w_m(binCount - 1, 0);                                //sum_l { w_lm }
        for (unsigned l = 0; l < numReplicas; ++l) {
            double temperatureExponenet = +betas[l] * deltaU;
            for (unsigned m = 0; m < binCount - 1; ++m) {
                double curHist = H_blm[b][l][m];
                double nextHist = H_blm[b][l][m + 1];
                if (curHist > 0 and nextHist > 0) {
                    x_lm[l][m] = (LogVal(nextHist * deltaU) / LogVal(curHist * deltaU)) * toLogValExp(temperatureExponenet);
                    w_lm[l][m] = curHist * nextHist / (curHist + nextHist);
                    w_m[m] += w_lm[l][m];
                } else {
                    x_lm[l][m] = LogVal();          //very small
                    w_lm[l][m] = 0;
                }
            }
        }

        lOmega_bm[b][0] = LogVal(1.0);
        for (unsigned m = 0; m < binCount - 1; ++m) {
            double ln_x_m = 0;
            for (unsigned l = 0; l < numReplicas; ++l) {
                if (w_m[m] > 0) ln_x_m += (w_lm[l][m] / w_m[m]) * x_lm[l][m].lnx;
            }

            lOmega_bm[b][m + 1].lnx = lOmega_bm[b][m].lnx + ln_x_m;
        }

        out << "." << flush;
    }

    out << "Done" << endl;
}

void MultireweightHistosPTJK::updateEffectiveCountsJK() {
    //jackknife-blocked:

    Heff_bm.resize(boost::extents[blockCount][binCount]);
    initArray(Heff_bm, 0);
    Neff_blm.resize(boost::extents[blockCount][numReplicas][binCount]);
    initArray(Neff_blm, 0);
    lHeff_bm.resize(boost::extents[blockCount][binCount]);
    lNeff_blm.resize(boost::extents[blockCount][numReplicas][binCount]);
    lPrecalc_blm.resize(boost::extents[blockCount][numReplicas][binCount]);

    #pragma omp parallel for
    for (signed b = 0; b < (signed)blockCount; ++b) {
        for (unsigned k = 0; k < numReplicas; ++k) {
            for (unsigned m = 0; m < binCount; ++m) {
                Heff_bm[b][m] += double(H_bkm[b][k][m]) / g_km[k][m];
                for (unsigned l = 0; l < numReplicas; ++l) {
                    Neff_blm[b][l][m] += double(N_bkl[b][k][l]) / g_km[k][m];
                }
            }
        }

        for (unsigned m = 0; m < binCount; ++m) {
            lHeff_bm[b][m] = (Heff_bm[b][m] != 0 ? LogVal(Heff_bm[b][m]) : LogVal(LogVal::LogZero));
            for (unsigned l = 0; l < numReplicas; ++l) {
                lNeff_blm[b][l][m] = (Neff_blm[b][l][m] != 0 ? LogVal(Neff_blm[b][l][m]) : LogVal(LogVal::LogZero));
                lPrecalc_blm[b][l][m] = lNeff_blm[b][l][m] * lBinSize * toLogValExp(-betas[l] * U_m[m]);    //for updateDensityOfStates
            }
        }
    }
}

inline void MultireweightHistosPTJK::updateDensityOfStatesJK(unsigned b) {
    //precalculated the part that does not depend on the estimates of the partition functions
    //lPrecalc_blm[b][l][m] = lNeff_blm[b][l][m] * lBinSize * toLogValExp(-betas[l] * U_m[m]);
    for (int m = 0; m < (signed)binCount; ++m) {
        lOmega_bm[b][m] = lHeff_bm[b][m];
        LogVal denominator = lPrecalc_blm[b][0][m] / lZ_bl[b][0];
        for (unsigned l = 1; l < numReplicas; ++l) {
            denominator += lPrecalc_blm[b][l][m] / lZ_bl[b][l];
        }
        lOmega_bm[b][m] /= denominator;
    }
}



void MultireweightHistosPTJK::findPartitionFunctionsAndDensityOfStates(double tolerance, int maxIterations) {
    //iteration for whole data set:
    MultireweightHistosPT::findPartitionFunctionsAndDensityOfStates(tolerance, maxIterations);

    //jackknife data sets:
    out << "Updating effective counts (JK)... " << flush;
    updateEffectiveCountsJK();
    out << "Done." << endl;

    out << "Starting iteration to estimate density of states (jackknife), tolerance=" << tolerance << " maxIterations=" << maxIterations << endl;

    #pragma omp parallel for
    for (signed b = 0; b < (signed)blockCount; ++b) {
        //starting value: result for the whole data set
        for (unsigned l = 0; l < numReplicas; ++l) {
            lZ_bl[b][l] = lZ_l[l];
        }
        updateDensityOfStatesJK(b);

        boost::multi_array<LogVal, 1> lZ_l_lastIteration = lZ_bl[b];

        int iterations = 0;

        double deltaSquared = 0;    //this is sum_l { ((Z_l - Z_l_old) / (Z_l)) ** 2 }, stop iteration once deltaSquared < tolerance ** 2
        do {
            ++iterations;

            deltaSquared = 0;

            //update estimates of partition functions
            lZ_bl[b][0] = lOmega_bm[b][0] * lBinSize * toLogValExp(-betas[0] * U_m[0]);
            for (unsigned m = 1; m < binCount; ++m) {
                lZ_bl[b][0] += lOmega_bm[b][m] * lBinSize * toLogValExp(-betas[0] * U_m[m]);
            }
            for (int l = 1; l < (signed)numReplicas; ++l) {
                lZ_l_lastIteration[l] = lZ_bl[b][l];            //store old value
                lZ_bl[b][l] = lOmega_bm[b][0] * lBinSize * toLogValExp(-betas[l] * U_m[0]);
                for (unsigned m = 1; m < binCount; ++m) {
                    lZ_bl[b][l] += lOmega_bm[b][m] * lBinSize * toLogValExp(-betas[l] * U_m[m]);
                }
                lZ_bl[b][l] /= lZ_bl[b][0];             //normalize

                deltaSquared += pow(expm1(lZ_bl[b][l].lnx - lZ_l_lastIteration[l].lnx), 2);     //gauge change from last iteration
            }
            lZ_bl[b][0].lnx = 0;            //sets to 1 (normalize)

            if (iterations % 10 == 0) {
                out << "block " << b << ", Iteration " << iterations << ", deltaSquared=" << deltaSquared << endl;
            }

            updateDensityOfStatesJK(b);
        } while (iterations < maxIterations and deltaSquared >= tolerance * tolerance);
    }
    out << "Done." << endl;
}


void MultireweightHistosPTJK::saveLogDensityOfStates(const std::string& filename) {
    //determine error
    vector<double> omegaError_m(binCount);
    for (unsigned m = 0; m < binCount; ++m) {
        double avg_lOmega_bm_lnx = 0;
        for (unsigned b = 0; b < blockCount; ++b) {
            avg_lOmega_bm_lnx += lOmega_bm[b][m].lnx;
        }
        avg_lOmega_bm_lnx /= blockCount;
        double sqDev = 0;
        for (unsigned b = 0; b < blockCount; ++b) {
//            sqDev += pow(lOmega_m[m].lnx - lOmega_bm[b][m].lnx, 2);
            sqDev += pow(avg_lOmega_bm_lnx - lOmega_bm[b][m].lnx, 2);
        }
        omegaError_m[m] = sqrt((double(blockCount - 1) / blockCount) * sqDev);
    }

    ofstream output(filename.c_str());
    output.precision(15);
    output.setf(std::ios::scientific, std::ios::floatfield);
    output  << "## logarithm of density of states, mrpt estimation\n"
            << "## energy (normalized)\t ln(d.o.s.)\t error (jk)"
            << endl;
    for (unsigned m = 0; m < binCount; ++m) {
        output  << U_m[m] / systemN
                << "\t" << lOmega_m[m] / lOmega_m[binCount / 2]
                << "\t" << omegaError_m[m]
                << '\n';
    }
}

void MultireweightHistosPTJK::saveLogDensityOfStatesIsing(const std::string& filename) {
    //determine error
    vector<double> omegaError_m(binCount);
    for (unsigned m = 0; m < binCount; ++m) {
        double avg_lOmega_bm_lnx = 0;
        for (unsigned b = 0; b < blockCount; ++b) {
            avg_lOmega_bm_lnx += lOmega_bm[b][m].lnx;
        }
        avg_lOmega_bm_lnx /= blockCount;
        double sqDev = 0;
        for (unsigned b = 0; b < blockCount; ++b) {
//            sqDev += pow(lOmega_m[m].lnx - lOmega_bm[b][m].lnx, 2);
            sqDev += pow(avg_lOmega_bm_lnx - lOmega_bm[b][m].lnx, 2);
        }
        omegaError_m[m] = sqrt((double(blockCount - 1) / blockCount) * sqDev);
    }

    ofstream output(filename.c_str());
    output.precision(15);
    output.setf(std::ios::scientific, std::ios::floatfield);
    output  << "## logarithm of density of states, mrpt estimation\n"
            << "## ln(d.o.s.) normalized to ln(2) for the lowest energy entry\n"
            << "## energy (normalized)\t ln(d.o.s.)\t error (jk)"
            << endl;
    LogVal norm = lOmega_m[0] / LogVal(2);
    //Normierung ist additiv im Logarithmus -- Fehler des
    //logarithmierten Werts bleibt also gleich

    //Normierung ist nur dann korrekt und stimmig mit Beale,
    //wenn U[0] == -2 * systemN
    //d.h., wenn Grundzustand wirklich erreicht!

    for (unsigned m = 0; m < binCount; ++m) {
        output  << U_m[m] / systemN
                << "\t" << lOmega_m[m] / norm
                << "\t" << omegaError_m[m]
                << '\n';
    }
}



ReweightingResult MultireweightHistosPTJK::reweightJackknifeInternal(
        double targetBeta, unsigned jkBlock,
        const DoubleSeriesCollection& w_kn) {
    //everything normalized by system volume:
    double meanEnergy = 0;
    double meanEnergySquared = 0;
    double meanObservable = 0;
    double meanObservableSquared = 0;
    double meanObservableToTheFourth = 0;

    for (int k = 0; k < (signed)numReplicas; ++k) {
        unsigned N_k = energyTimeSeries[k]->size();
        unsigned jkBlockSize = N_k / blockCount;
        for (unsigned n = 0; n < N_k; ++n) {
            const unsigned curBlock = n / jkBlockSize;
            if (curBlock == jkBlock) {
                n += jkBlockSize;
                continue;
            } else {
                double e = (*energyTimeSeries[k])[n];
                double o = (*observableTimeSeries[k])[n];
                meanEnergy += e * (*w_kn[k])[n];
                meanEnergySquared += e*e * (*w_kn[k])[n];
                meanObservable += o * (*w_kn[k])[n];
                meanObservableSquared += o*o * (*w_kn[k])[n];
                meanObservableToTheFourth += o*o*o*o * (*w_kn[k])[n];
            }
        }
    }

    double heatCapacity = double(systemN) * targetBeta * targetBeta
        * (meanEnergySquared - pow(meanEnergy, 2));
    double suscObservable = double(systemN) * (meanObservableSquared - pow(meanObservable, 2));
    double binderObservable = 1.0 - (meanObservableToTheFourth / (3 * pow(meanObservableSquared, 2)));

    return ReweightingResult(meanEnergy, heatCapacity, meanObservable, suscObservable, binderObservable);
}


void MultireweightHistosPTJK::reweight1stMomentInternalJK(const DoubleSeriesCollection& timeSeries, const DoubleSeriesCollection& w_kn, unsigned jkBlock, double& firstMoment) {
    double first = 0;
    for (int k = 0; k < (signed)numReplicas; ++k) {
        unsigned N_k = timeSeries[k]->size();
        unsigned jkBlockSize = N_k / blockCount;
        for (unsigned n = 0; n < N_k; ++n) {
            const unsigned curBlock = n / jkBlockSize;
            if (curBlock == jkBlock) {
                n += jkBlockSize;
                continue;
            } else {
                double v = (*timeSeries[k])[n];
                first += v * (*w_kn[k])[n];
            }
        }
    }
    firstMoment = first;
}

void MultireweightHistosPTJK::reweight1stMoment2ndMomentInternalJK(
        const DoubleSeriesCollection& timeSeries,
        const DoubleSeriesCollection& w_kn, unsigned jkBlock,
        double& firstMoment, double& secondMoment) {
    double first = 0;
    double second = 0;
    for (int k = 0; k < (signed)numReplicas; ++k) {
        unsigned N_k = timeSeries[k]->size();
        unsigned jkBlockSize = N_k / blockCount;
        for (unsigned n = 0; n < N_k; ++n) {
            const unsigned curBlock = n / jkBlockSize;
            if (curBlock == jkBlock) {
                n += jkBlockSize;
                continue;
            } else {
                double v = (*timeSeries[k])[n];
                first += v * (*w_kn[k])[n];
                second += v*v * (*w_kn[k])[n];
            }
        }
    }
    firstMoment = first;
    secondMoment = second;
}

void MultireweightHistosPTJK::reweight2ndMoment4thMomentInternalJK(
        const DoubleSeriesCollection& timeSeries,
        const DoubleSeriesCollection& w_kn, unsigned jkBlock,
        double& secondMoment, double& fourthMoment) {
    double second = 0;
    double fourth = 0;
    for (int k = 0; k < (signed)numReplicas; ++k) {
        unsigned N_k = timeSeries[k]->size();
        unsigned jkBlockSize = N_k / blockCount;
        for (unsigned n = 0; n < N_k; ++n) {
            const unsigned curBlock = n / jkBlockSize;
            if (curBlock == jkBlock) {
                n += jkBlockSize;
                continue;
            } else {
                double v = (*timeSeries[k])[n];
                second += v*v * (*w_kn[k])[n];
                fourth += v*v*v*v * (*w_kn[k])[n];
            }
        }
    }
    fourthMoment = fourth;
    secondMoment = second;
}


double MultireweightHistosPTJK::reweightEnergyJK(double targetBeta, unsigned jkBlock) {
    DoubleSeriesCollection w_kn = computeWeightsJK(targetBeta, jkBlock);

    double result = 0;
    reweight1stMomentInternalJK(energyTimeSeries, w_kn, jkBlock, result);

    destroyAll(w_kn);

    return result;
}

double MultireweightHistosPTJK::reweightSpecificHeatJK(double targetBeta, unsigned jkBlock) {
    DoubleSeriesCollection w_kn = computeWeightsJK(targetBeta, jkBlock);

    double firstMoment = 0;
    double secondMoment = 0;
    reweight1stMoment2ndMomentInternalJK(energyTimeSeries, w_kn, jkBlock, firstMoment, secondMoment);
    double result = systemN * pow(targetBeta, 2) *
            (secondMoment - pow(firstMoment, 2));;

    destroyAll(w_kn);

    return result;
}

double MultireweightHistosPTJK::reweightObservableJK(double targetBeta, unsigned jkBlock) {
    DoubleSeriesCollection w_kn = computeWeightsJK(targetBeta, jkBlock);

    double result = 0;
    reweight1stMomentInternalJK(observableTimeSeries, w_kn, jkBlock, result);

    destroyAll(w_kn);

    return result;
}

double MultireweightHistosPTJK::reweightObservableSusceptibilityJK(
        double targetBeta, unsigned jkBlock) {
    DoubleSeriesCollection w_kn = computeWeightsJK(targetBeta, jkBlock);

    double firstMoment = 0;
    double secondMoment = 0;
    reweight1stMoment2ndMomentInternalJK(observableTimeSeries, w_kn, jkBlock, firstMoment, secondMoment);
    double result = systemN * (secondMoment - pow(firstMoment, 2));;

    destroyAll(w_kn);

    return result;
}

double MultireweightHistosPTJK::reweightObservableBinderJK(double targetBeta,
        unsigned jkBlock) {
    DoubleSeriesCollection w_kn = computeWeightsJK(targetBeta, jkBlock);

    double secondMoment = 0;
    double fourthMoment = 0;
    reweight2ndMoment4thMomentInternalJK(observableTimeSeries, w_kn, jkBlock,
            secondMoment, fourthMoment);
    double result = 1.0 - (fourthMoment / (3 * pow(secondMoment, 2)));

    destroyAll(w_kn);

    return result;
}

MultireweightHistosPTJK::DoubleSeriesCollection MultireweightHistosPTJK::computeWeightsJK(double targetBeta, unsigned jkBlock) {
//  out << "Computing weights w_kn at beta=" << targetBeta << ", jackknife block:" << jkBlock << endl;

    vector<LogVal> arguments(binCount);
    arguments[0] = lOmega_bm[jkBlock][0] * toLogValExp(-targetBeta * U_m[0]);
    LogVal normalization = arguments[0];
    for (unsigned m = 1; m < binCount; ++m) {
        arguments[m] = lOmega_bm[jkBlock][m] * toLogValExp(-targetBeta * U_m[m]);
        normalization += arguments[m];
    }
    //normalize arguments, calculate weight corresponding to bin
    vector<double> weightFromBin(binCount);
    for (unsigned m = 0; m < binCount; ++m) {
        if (H_bm[jkBlock][m] == 0) {
            weightFromBin[m] = 0;
        } else {
            arguments[m] /= normalization;
            weightFromBin[m] = toDouble(arguments[m]) / double(H_bm[jkBlock][m]);
        }
    }

    DoubleSeriesCollection w_kn(numReplicas, 0);

    //calculate weight for each sample
    for (int k = 0; k < (signed)numReplicas; ++k) {
        unsigned N_k = m_kn[k]->size();
        w_kn[k] = new vector<double>(N_k);
        unsigned jkBlockSize = N_k / blockCount;
        for (unsigned n = 0; n < N_k; ++n) {
            unsigned curBlock = (n * blockCount) / N_k;
            if (curBlock == jkBlock) {
                n += jkBlockSize;
                continue;
            } else {
                unsigned m = (*m_kn[k])[n];
                (*w_kn[k])[n] = weightFromBin[m];
            }
        }
    }

    return w_kn;
}

ReweightingResult MultireweightHistosPTJK::reweight(double targetBeta) {
    out << "Jackknife estimation of energy and " << observable << " moments at beta=" << targetBeta << endl;

    //jackknifed dataset:
    vector<double> jkEnergy_b(blockCount);
    vector<double> jkHeatCapacity_b(blockCount);
    vector<double> jkObs_b(blockCount);
    vector<double> jkSusc_b(blockCount);
    vector<double> jkBinder_b(blockCount);
    #pragma omp parallel for
    for (signed b = 0; b < (signed)blockCount; ++b) {
        DoubleSeriesCollection w_kn = computeWeightsJK(targetBeta, b);

        ReweightingResult results =
                reweightJackknifeInternal(targetBeta, b, w_kn);
        jkEnergy_b[b] = results.energyAvg;
        jkHeatCapacity_b[b] = results.heatCapacity;
        jkObs_b[b] = results.obsAvg;
        jkSusc_b[b] = results.obsSusc;
        jkBinder_b[b] = results.obsBinder;

        destroyAll(w_kn);
        out << '#' << flush;
    }

    //Jackknife for error and average (--> bias correction)
    ReweightingResult totalResult;
    jackknife(totalResult.energyAvg, totalResult.energyError, jkEnergy_b);
    jackknife(totalResult.heatCapacity, totalResult.heatCapacityError,
            jkHeatCapacity_b);
    jackknife(totalResult.obsAvg, totalResult.obsError, jkObs_b);
    jackknife(totalResult.obsSusc, totalResult.obsSuscError, jkSusc_b);
    jackknife(totalResult.obsBinder, totalResult.obsBinderError, jkBinder_b);

    out << " Done." << endl;

    return totalResult;
}


HistogramDouble* MultireweightHistosPTJK::reweightEnergyHistogram(double targetBeta) {
    out << "Jackknife estimation of energy histogram at beta=" << targetBeta << ", " << binCount << " bins... " << endl;

    //helper array:
    vector<LogVal> arguments(binCount);
    //histogram from whole data set
    vector<double> histo_m(binCount);
    arguments[0] = lOmega_m[0] * toLogValExp(-targetBeta * U_m[0]);
    LogVal normalization = arguments[0];
    for (unsigned m = 1; m < binCount; ++m) {
        arguments[m] = lOmega_m[m] * toLogValExp(-targetBeta * U_m[m]);
        normalization += arguments[m];
    }
    for (unsigned m = 0; m < binCount; ++m) {
        arguments[m] /= normalization;
        histo_m[m] = toDouble(arguments[m]);
    }

    //histograms created from jacknifed sub-sets
    Double2Array histo_bm(boost::extents[blockCount][binCount]);
    initArray(histo_bm, 0);
    #pragma omp parallel for
    for (signed b = 0; b < (signed)blockCount; ++ b) {
        vector<LogVal> args(binCount);
        args[0] = lOmega_bm[b][0] * toLogValExp(-targetBeta * U_m[0]);
        LogVal normalization = args[0];
        for (unsigned m = 1; m < binCount; ++m) {
            args[m] = lOmega_bm[b][m] * toLogValExp(-targetBeta * U_m[m]);
            normalization += args[m];
        }
        for (unsigned m = 0; m < binCount; ++m) {
            args[m] /= normalization;
            histo_bm[b][m] = toDouble(args[m]);
        }
    }

    //error estimation
    vector<double> sqDev(binCount, 0);
    for (unsigned b = 0; b < blockCount; ++b) {
        for (unsigned m = 0; m < binCount; ++m) {
            sqDev[m] += pow(histo_bm[b][m] - histo_m[m], 2);
        }
    }
    vector<double> errors(binCount);
    double factor = double(blockCount - 1) / double(blockCount);
    for (unsigned m = 0; m < binCount; ++m) {
        errors[m] = sqrt(factor * sqDev[m]);
    }

    //prepare instance of nice histogram class to return
    HistogramDouble* result = new HistogramDouble;
    result->assignVector(histo_m,
            minEnergyNormalized, maxEnergyNormalized, targetBeta, systemN);
    result->assignErrorBarVector(errors);
    result->headerLines +=
            "## MRPT reweighted histogram of normalized energy\n";
    result->updateMeta();

    return result;
}


HistogramDouble* MultireweightHistosPTJK::reweightObservableHistogram(double targetBeta, unsigned obsBinCount) {
    out << "Jackknife estimation of " << observable << " histogram at beta=" << targetBeta << ", " << obsBinCount << " bins... " << endl;

    //histogram created from the whole data-set
    vector<double> obsHisto(obsBinCount, 0.0);
    DoubleSeriesCollection total_w_kn = computeWeights(targetBeta);
    MultireweightHistosPT::reweightObservableHistogramInternal(targetBeta, obsBinCount, total_w_kn, obsHisto);
    destroyAll(total_w_kn);

    HistogramDouble* result = new HistogramDouble;
    result->assignVector(obsHisto, minObservableNormalized, maxObservableNormalized, targetBeta, systemN);
    result->headerLines += "## MRPT reweighted histogram of normalized " + observable + "\n";
    result->updateMeta();

    const double SMALL = 1e-10;     //to fit maxObservableNormalized into the highest bin
    double obsBinSize = (maxObservableNormalized - minObservableNormalized + SMALL) / obsBinCount;

    Double2Array obsHisto_bm(boost::extents[blockCount][obsBinCount]);
    initArray(obsHisto_bm, 0);

    #pragma omp parallel for
    for (signed b = 0; b < (signed)blockCount; ++ b) {
        DoubleSeriesCollection w_kn = computeWeightsJK(targetBeta, b);
        for (unsigned k = 0; k < numReplicas; ++k) {
            unsigned N_k = observableTimeSeries[k]->size();
            unsigned jkBlockSize = N_k / blockCount;
            for (unsigned n = 0; n < N_k; ++n) {
                const int curBlock = n / jkBlockSize;
                if (curBlock == b) {
                    n += jkBlockSize;
                    continue;
                } else {
                    unsigned curBin = static_cast<int>(((*observableTimeSeries[k])[n] - minObservableNormalized) / obsBinSize);
                    obsHisto_bm[b][curBin] += (*w_kn[k])[n];
                }
            }
        }
        destroyAll(w_kn);
        cout << "#" << flush;
    }

    vector<double> sqDev(obsBinCount, 0);
    for (unsigned b = 0; b < blockCount; ++b) {
        for (unsigned m = 0; m < obsBinCount; ++m) {
            sqDev[m] += pow(obsHisto_bm[b][m] - obsHisto[m], 2);
        }
    }
    vector<double> errors(obsBinCount);
    double factor = double(blockCount - 1) / double(blockCount);
    for (unsigned m = 0; m < obsBinCount; ++m) {
        errors[m] = sqrt(factor * sqDev[m]);
    }

    result->assignErrorBarVector(errors);

    cout << "Done." << endl;

    return result;
}

HistogramDouble* MultireweightHistosPTJK::reweightEnergyHistogramJK(
        double targetBeta, unsigned jkBlock) {
    vector<double> histo_m(binCount, 0.0);

    vector<LogVal> args(binCount);
    args[0] = lOmega_bm[jkBlock][0] * toLogValExp(-targetBeta * U_m[0]);
    LogVal normalization = args[0];
    for (unsigned m = 1; m < binCount; ++m) {
        args[m] = lOmega_bm[jkBlock][m] * toLogValExp(-targetBeta * U_m[m]);
        normalization += args[m];
    }
    for (unsigned m = 0; m < binCount; ++m) {
        args[m] /= normalization;
        histo_m[m] = toDouble(args[m]);
    }

    //prepare instance of nice histogram class to return
    HistogramDouble* result = new HistogramDouble;
    result->assignVector(histo_m,
            minEnergyNormalized, maxEnergyNormalized, targetBeta, systemN);
    result->headerLines +=
            "## MRPT reweighted histogram of normalized energy\n";
    result->updateMeta();

    return result;
}

HistogramDouble* MultireweightHistosPTJK::reweightObservableHistogramJK(
        double targetBeta, unsigned obsBinCount, unsigned jkBlock) {
    const double SMALL = 1e-10;     //to fit maxObservableNormalized into the highest bin
    double obsBinSize = (maxObservableNormalized - minObservableNormalized + SMALL) / obsBinCount;
    vector<double> obsHisto_m(obsBinCount, 0.0);

    DoubleSeriesCollection w_kn = computeWeightsJK(targetBeta, jkBlock);
    for (unsigned k = 0; k < numReplicas; ++k) {
        unsigned N_k = observableTimeSeries[k]->size();
        unsigned jkBlockSize = N_k / blockCount;
        for (unsigned n = 0; n < N_k; ++n) {
            const int curBlock = n / jkBlockSize;
            if (curBlock == jkBlock) {
                n += jkBlockSize;
                continue;
            } else {
                unsigned curBin = static_cast<int>(
                        ((*observableTimeSeries[k])[n] -
                                minObservableNormalized) / obsBinSize);
                obsHisto_m[curBin] += (*w_kn[k])[n];
            }
        }
    }
    destroyAll(w_kn);

    HistogramDouble* result = new HistogramDouble;
    result->assignVector(obsHisto_m, minObservableNormalized, maxObservableNormalized,
            targetBeta, systemN);
    result->updateMeta();

    return result;
}



//TODO:unnecessary weight-recalculcation in reweightObservableHistogram
ReweightingResult MultireweightHistosPTJK::reweightWithHistograms(double targetBeta, unsigned  obsBinCount) {
    ReweightingResult results = reweight(targetBeta);
    results.energyHistogram = reweightEnergyHistogram(targetBeta);
    results.obsHistogram = reweightObservableHistogram(targetBeta, obsBinCount);

    return results;
}


class SuscMinCallableJK {
    MultireweightHistosPTJK* master; map<double, double>& pointsEvaluated; unsigned jkBlock;
public:
    SuscMinCallableJK(MultireweightHistosPTJK* master, map<double, double>& pointsEvaluated, unsigned jkBlock) :
        master(master), pointsEvaluated(pointsEvaluated), jkBlock(jkBlock) { }
    double operator()(double beta) {
        double susc = master->reweightObservableSusceptibilityJK(beta, jkBlock);
        pointsEvaluated[beta] = susc;
//      master->out << " b:" << jkBlock << ", beta: " << beta << " => " << susc << endl;
        return -susc;
    }
};

class BinderMinCallableJK {
    MultireweightHistosPTJK* master; map<double, double>& pointsEvaluated; unsigned jkBlock;
public:
    BinderMinCallableJK(MultireweightHistosPTJK* master,
            map<double, double>& pointsEvaluated, unsigned jkBlock) :
        master(master), pointsEvaluated(pointsEvaluated), jkBlock(jkBlock) { }
    double operator()(double beta) {
        double binder = master->reweightObservableBinderJK(beta, jkBlock);
        pointsEvaluated[beta] = binder;
        return binder;
    }
};

class SpecificHeatDiscreteMinCallableJK {
    MultireweightHistosPTJK* master; map<double, double>& pointsEvaluated; unsigned jkBlock;
public:
    SpecificHeatDiscreteMinCallableJK (MultireweightHistosPTJK* master, map<double, double>& pointsEvaluated, unsigned jkBlock) :
        master(master), pointsEvaluated(pointsEvaluated), jkBlock(jkBlock) { }
    double operator()(double beta) {
        double specHeat = master->reweightSpecificHeatDiscreteJK(beta, jkBlock);
        pointsEvaluated[beta] = specHeat;
//      master->out << " b:" << jkBlock << ", beta: " << beta << " => " << specHeat << endl;
        return -specHeat;
    }
};

//take only data from one jackknife block
//these histograms are no longer needed later on --> destructor frees them
class EnergyHistogramPeakDiffMinCallableJK {
    MultireweightHistosPTJK* master;
    HistogramDouble* foundHistogram;
    double tolerance;
    unsigned jkBlock;
public:
    EnergyHistogramPeakDiffMinCallableJK(double tolerance,
            MultireweightHistosPTJK* master, unsigned jkBlock) :
        master(master), foundHistogram(0), jkBlock(jkBlock),
        tolerance(tolerance)
    {}
    double operator()(double beta) {
        destroy(foundHistogram);
        foundHistogram = master->reweightEnergyHistogramJK(beta, jkBlock);
        double peakDiff = histogramPeakDiff(foundHistogram, tolerance);
        return peakDiff;
    }
    HistogramDouble* getHistogram() {
        return foundHistogram;
    }
    ~EnergyHistogramPeakDiffMinCallableJK() {
        destroy(foundHistogram);
    }
};

class ObsHistogramPeakDiffMinCallableJK {
    MultireweightHistosPTJK* master;
    HistogramDouble* foundHistogram;
    double tolerance;
    unsigned numBins;
    unsigned jkBlock;
public:
    ObsHistogramPeakDiffMinCallableJK(double tolerance,
            unsigned numBins, MultireweightHistosPTJK* master,
            unsigned jkBlock) :
        master(master), foundHistogram(0),
        tolerance(tolerance), numBins(numBins), jkBlock(jkBlock)
    {}
    double operator()(double beta) {
        destroy(foundHistogram);
        foundHistogram = master->
                reweightObservableHistogramJK(beta, numBins, jkBlock);
        double peakDiff = histogramPeakDiff(foundHistogram, tolerance);
        return peakDiff;
    }
    HistogramDouble* getHistogram() {
        return foundHistogram;
    }
    ~ObsHistogramPeakDiffMinCallableJK() {
        destroy(foundHistogram);
    }
};

class EnergyHistogramWeightDiffMinCallableJK {
    MultireweightHistosPTJK* master;
    HistogramDouble* foundHistogram;
    double tolerance;
    double cutOff;
    unsigned jkBlock;
public:
    EnergyHistogramWeightDiffMinCallableJK(double tolerance,
            double cutOff, MultireweightHistosPTJK* master, unsigned jkBlock) :
        master(master), foundHistogram(0), tolerance(tolerance),
        cutOff(cutOff), jkBlock(jkBlock)
    {}
    double operator()(double beta) {
        destroy(foundHistogram);
        foundHistogram = master->reweightEnergyHistogramJK(beta, jkBlock);
        double weightDiff = histogramWeightDiff(foundHistogram,
                cutOff);
        return weightDiff;
    }
    HistogramDouble* getHistogram() {
        return foundHistogram;
    }
    ~EnergyHistogramWeightDiffMinCallableJK() {
        destroy(foundHistogram);
    }
};

class ObsHistogramWeightDiffMinCallableJK {
    MultireweightHistosPTJK* master;
    HistogramDouble* foundHistogram;
    double tolerance;
    double cutOff;
    unsigned numBins;
    unsigned jkBlock;
public:
    ObsHistogramWeightDiffMinCallableJK(double tolerance,
            double cutOff, unsigned numBins, MultireweightHistosPTJK* master,
            unsigned jkBlock) :
        master(master), foundHistogram(0), tolerance(tolerance), cutOff(cutOff),
        numBins(numBins), jkBlock(jkBlock)
    {}
    double operator()(double beta) {
        destroy(foundHistogram);
        foundHistogram = master->
                reweightObservableHistogramJK(beta, numBins, jkBlock);
        double weightDiff = histogramWeightDiff(foundHistogram,
                cutOff);
        return weightDiff;
    }
    HistogramDouble* getHistogram() {
        return foundHistogram;
    }
    ~ObsHistogramWeightDiffMinCallableJK() {
        destroy(foundHistogram);
    }
};


void MultireweightHistosPTJK::findMaxObservableSusceptibility(double& betaMax, double& betaMaxError, double& suscMax, double& suscMaxError,
        std::map<double, double>& pointsEvaluated, std::vector<std::map<double,double> >& pointsEvaluatedJK,
        double betaStart, double betaEnd) {
    //whole data-set:
    MultireweightHistosPT::findMaxObservableSusceptibility(betaMax, suscMax, pointsEvaluated, betaStart, betaEnd);

    out << "Jackknifed search for maximum of " << observable << " susceptibility" << endl;

    //jackknife blocking
    vector<double> jkSusc_b(blockCount);
    vector<double> jkBeta_b(blockCount);
    #pragma omp parallel for
    for (signed b = 0; b < (signed)blockCount; ++b) {
        SuscMinCallableJK f(this, pointsEvaluatedJK[b], b);
        brentMinimize(jkBeta_b[b], jkSusc_b[b], f, betaStart, betaEnd);
        jkSusc_b[b] *= -1;
        out << "final b:" << b << ", beta: " << jkBeta_b[b] << " => " << jkSusc_b[b] << endl;
    }

    double jkSuscBlockAverage;
    double jkBetaBlockAverage;
    jackknife(jkSuscBlockAverage, suscMaxError, jkSusc_b);
    jackknife(jkBetaBlockAverage, betaMaxError, jkBeta_b);
    out << "averaged beta: " << jkBetaBlockAverage << endl;
    out << "averaged susc: " << jkSuscBlockAverage << endl;

    out << "Final result: maximum of " << observable << " susceptibility: " << suscMax << ", error: " << suscMaxError << '\n'
        << "at beta=" << betaMax << ", error: " << betaMaxError << endl;
}

void MultireweightHistosPTJK::findMinBinder(double& betaMin,
        double& betaMinError, double& binderMin, double& binderMinError,
        std::map<double, double>& pointsEvaluated,
        std::vector<std::map<double, double> >& pointsEvaluatedJK,
        double betaStart, double betaEnd) {
    //whole data-set:
    MultireweightHistosPT::findMinBinder(betaMin, binderMin,
            pointsEvaluated, betaStart, betaEnd);

    out << "Jackknifed search for minimum of "
        << observable << " binder cumulant" << endl;

    //jackknife blocking
    vector<double> jkBinder_b(blockCount);
    vector<double> jkBeta_b(blockCount);
    #pragma omp parallel for
    for (signed b = 0; b < (signed)blockCount; ++b) {
        BinderMinCallableJK f(this, pointsEvaluatedJK[b], b);
        brentMinimize(jkBeta_b[b], jkBinder_b[b], f, betaStart, betaEnd);
        out << "final b:" << b << ", beta: " << jkBeta_b[b] << " => " << jkBinder_b[b] << endl;
    }

    double jkBinderBlockAverage;
    double jkBetaBlockAverage;
    jackknife(jkBinderBlockAverage, binderMinError, jkBinder_b);
    jackknife(jkBetaBlockAverage, betaMinError, jkBeta_b);
    out << "averaged beta:   " << jkBetaBlockAverage << endl;
    out << "averaged binder: " << jkBinderBlockAverage << endl;

    out << "Final result: minimum of " << observable << " binder cumulant: "
        << binderMin << ", error: " << binderMinError << '\n'
        << "at beta=" << betaMin << ", error: " << betaMinError << endl;
}


void MultireweightHistosPTJK::findMaxSpecificHeatDiscrete(
        double& betaMax, double& betaMaxError,
        double& specificHeatMax, double& specificHeatMaxError,
        std::map<double, double>& pointsEvaluated, std::vector<std::map<double,double> >& pointsEvaluatedJK,
        double betaStart, double betaEnd) {
    //whole data-set:
    MultireweightHistosPT::findMaxSpecificHeatDiscrete(betaMax, specificHeatMax, pointsEvaluated, betaStart, betaEnd);

    out << "Jackknifed search for maximum of specific heat (discrete)" << endl;

    //jackknife blocking
    vector<double> jkSH_b(blockCount);
    vector<double> jkBeta_b(blockCount);
    #pragma omp parallel for
    for (signed b = 0; b < (signed)blockCount; ++b) {
        SpecificHeatDiscreteMinCallableJK f(this, pointsEvaluatedJK[b], b);
        brentMinimize(jkBeta_b[b], jkSH_b[b], f, betaStart, betaEnd);
        jkSH_b[b] *= -1;
        out << "final b:" << b << ", beta: " << jkBeta_b[b] << " => " << jkSH_b[b] << endl;
    }

    double jkSHBlockAverage;
    double jkBetaBlockAverage;
    jackknife(jkSHBlockAverage, specificHeatMaxError, jkSH_b);
    jackknife(jkBetaBlockAverage, betaMaxError, jkBeta_b);
    out << "averaged beta: " << jkBetaBlockAverage << endl;
    out << "averaged s.h.: " << jkSHBlockAverage << endl;

    out << "Final result: maximum of specific heat: " << specificHeatMax << ", error: " << specificHeatMaxError << '\n'
        << "at beta=" << betaMax << ", error: " << betaMaxError << endl;
}


void MultireweightHistosPTJK::findEnergyEqualHeight(
        double& betaDouble, double& betaDoubleError,
        double& relDip, double& relDipError, HistogramDouble*& histo,
        double betaStart, double betaEnd, double tolerance) {
    //whole data-set
    MultireweightHistosPT::findEnergyEqualHeight(betaDouble, relDip,
            histo, betaStart, betaEnd, tolerance);

    out << "Jackknifed search for energy histogram with equal "
            "height double peak "
        << "on " << blockCount << " blocks" << endl;

    //jackknife blocking
    vector<double> jkRelDip_b(blockCount);
    vector<double> jkBeta_b(blockCount);
    #pragma omp parallel for
    for (signed b = 0; b < (signed)blockCount; ++b) {
        EnergyHistogramPeakDiffMinCallableJK f(tolerance, this, b);
        double peakDiff = -1.;
        brentMinimize(jkBeta_b[b], peakDiff, f, betaStart, betaEnd);
        jkRelDip_b[b] = histogramRelativeDip(f.getHistogram(), tolerance);
        out << "final b:" << b << ", beta: " << jkBeta_b[b]
                << " => peakDiff: " << peakDiff << ", relativeDip: "
                << jkRelDip_b[b] << endl;
    }

    double jkRelDipBlockAverage;
    double jkBetaBlockAverage;
    jackknife(jkRelDipBlockAverage, relDipError, jkRelDip_b);
    jackknife(jkBetaBlockAverage, betaDoubleError, jkBeta_b);
    out << "averaged beta:   " << jkBetaBlockAverage << endl;
    out << "averaged relDip: " << jkRelDipBlockAverage << endl;

    out << "Final result: relative dip: " << relDip << ", error: " << relDipError << '\n'
        << "at beta=" << betaDouble << ", error: " << betaDoubleError << endl;
}

void MultireweightHistosPTJK::findObservableEqualHeight(
        double& betaDouble, double& betaDoubleError,
        double& relDip, double& relDipError, HistogramDouble*& histo,
        double betaStart, double betaEnd,
        unsigned numBins, double tolerance) {
    //whole data-set
    MultireweightHistosPT::findObsEqualHeight(betaDouble, relDip,
            histo, betaStart, betaEnd, numBins, tolerance);

    out << "Jackknifed search for " << observable <<
           " histogram with equal height double peak "
        << "on " << blockCount << " blocks" << endl;

    //jackknife blocking
    vector<double> jkRelDip_b(blockCount);
    vector<double> jkBeta_b(blockCount);
    #pragma omp parallel for
    for (signed b = 0; b < (signed)blockCount; ++b) {
        ObsHistogramPeakDiffMinCallableJK f(tolerance, numBins, this, b);
        double peakDiff = -1.;
        brentMinimize(jkBeta_b[b], peakDiff, f, betaStart, betaEnd);
        jkRelDip_b[b] = histogramRelativeDip(f.getHistogram(), tolerance);
        out << "final b:" << b << ", beta: " << jkBeta_b[b]
                << " => peakDiff: " << peakDiff << ", relativeDip: "
                << jkRelDip_b[b] << endl;
    }

    double jkRelDipBlockAverage = average(&jkRelDip_b);
    double jkBetaDoubleBlockAverage = average(&jkBeta_b);
    jackknife(jkRelDipBlockAverage, relDipError, jkRelDip_b);
    jackknife(jkBetaDoubleBlockAverage, betaDoubleError, jkBeta_b);
    out << "averaged beta:   " << jkBetaDoubleBlockAverage << endl;
    out << "averaged relDip: " << jkRelDipBlockAverage << endl;

    out << "Final result: relative dip: " << relDip << ", error: " << relDipError << '\n'
        << "at beta=" << betaDouble << ", error: " << betaDoubleError << endl;
}

void MultireweightHistosPTJK::findEnergyEqualHeightWeight(
        double& betaDoubleEH, double& betaDoubleErrorEH,
        double& relDipEH, double& relDipErrorEH,
        HistogramDouble*& histoResultEH,
        double& betaDoubleEW, double& betaDoubleErrorEW,
        double& relDipEW, double& relDipErrorEW,
        HistogramDouble*& histoResultEW,
        double betaStart, double betaEnd, double tolerance) {
    //whole data-set
    MultireweightHistosPT::findEnergyEqualHeight(betaDoubleEH, relDipEH,
            histoResultEH, betaStart, betaEnd, tolerance);
    MultireweightHistosPT::findEnergyEqualWeight(betaDoubleEW, relDipEW,
            histoResultEW, histoResultEH, betaStart, betaEnd, tolerance);

    out << "Jackknifed search for energy histograms with equal "
           "height and weight double peak "
        << "on " << blockCount << " blocks" << endl;

    vector<double> jkRelDipEH_b(blockCount);
    vector<double> jkBetaEH_b(blockCount);
    vector<double> jkRelDipEW_b(blockCount);
    vector<double> jkBetaEW_b(blockCount);

    #pragma omp parallel for
    for (signed b = 0; b < (signed)blockCount; ++b) {
        EnergyHistogramPeakDiffMinCallableJK fEH(tolerance, this, b);
        double peakDiff = -1.;
        brentMinimize(jkBetaEH_b[b], peakDiff, fEH, betaStart, betaEnd);
        jkRelDipEH_b[b] = histogramRelativeDip(fEH.getHistogram(), tolerance);
        out << "final equal height, b:" << b << ", beta: " << jkBetaEH_b[b]
                << " => peakDiff: " << peakDiff << ", relativeDip: "
                << jkRelDipEH_b[b] << endl;
        double cut = histogramMinimumLocation(fEH.getHistogram(), tolerance);
        EnergyHistogramWeightDiffMinCallableJK fEW(tolerance, cut, this, b);
        double weightDiff = -1.;
        brentMinimize(jkBetaEW_b[b], weightDiff, fEW, betaStart, betaEnd);
        jkRelDipEW_b[b] = histogramRelativeDip(fEW.getHistogram(), tolerance);
        out << "final equal weight, b:" << b << ", beta: " << jkBetaEW_b[b]
                << " => weightDiff: " << weightDiff << ", relativeDip: "
                << jkRelDipEW_b[b] << endl;
    }
    double jkRelDipBlockAverageEH;
    double jkBetaBlockAverageEH;
    jackknife(jkRelDipBlockAverageEH, relDipErrorEH, jkRelDipEH_b);
    jackknife(jkBetaBlockAverageEH, betaDoubleErrorEH, jkBetaEH_b);
    out << "EH averaged beta:   " << jkBetaBlockAverageEH << endl;
    out << "EH averaged relDip: " << jkRelDipBlockAverageEH << endl;
    double jkRelDipBlockAverageEW;
    double jkBetaBlockAverageEW;
    jackknife(jkRelDipBlockAverageEW, relDipErrorEW, jkRelDipEW_b);
    jackknife(jkBetaBlockAverageEW, betaDoubleErrorEW, jkBetaEW_b);
    out << "EW averaged beta:   " << jkBetaBlockAverageEW << endl;
    out << "EW averaged relDip: " << jkRelDipBlockAverageEW << endl;

    out << "Equal Height Final result: relative dip: " << relDipEH
        << ", error: " << relDipErrorEH << '\n'
        << "at beta=" << betaDoubleEH << ", error: "
        << betaDoubleErrorEH << endl;
    out << "Equal Weight Final result: relative dip: " << relDipEW
        << ", error: " << relDipErrorEW << '\n'
        << "at beta=" << betaDoubleEW << ", error: "
        << betaDoubleErrorEW << endl;

}

void MultireweightHistosPTJK::findObservableEqualHeightWeight(
        double& betaDoubleEH, double& betaDoubleErrorEH,
        double& relDipEH, double& relDipErrorEH,
        HistogramDouble*& histoResultEH,
        double& betaDoubleEW, double& betaDoubleErrorEW,
        double& relDipEW, double& relDipErrorEW,
        HistogramDouble*& histoResultEW,
        double betaStart, double betaEnd, unsigned numBins,
        double tolerance) {
    //whole data-set
    MultireweightHistosPT::findObsEqualHeight(betaDoubleEH, relDipEH,
            histoResultEH, betaStart, betaEnd, numBins, tolerance);
    MultireweightHistosPT::findObsEqualWeight(betaDoubleEW, relDipEW,
            histoResultEW, histoResultEH, betaStart, betaEnd, numBins,
            tolerance);

    out << "Jackknifed search for " << observable << "histograms with equal "
           "height and weight double peak "
        << "on " << blockCount << " blocks" << endl;

    vector<double> jkRelDipEH_b(blockCount);
    vector<double> jkBetaEH_b(blockCount);
    vector<double> jkRelDipEW_b(blockCount);
    vector<double> jkBetaEW_b(blockCount);

    #pragma omp parallel for
    for (signed b = 0; b < (signed)blockCount; ++b) {
        ObsHistogramPeakDiffMinCallableJK fEH(tolerance, numBins, this, b);
        double peakDiff = -1.;
        brentMinimize(jkBetaEH_b[b], peakDiff, fEH, betaStart, betaEnd);
        jkRelDipEH_b[b] = histogramRelativeDip(fEH.getHistogram(), tolerance);
        out << "final equal height, b:" << b << ", beta: " << jkBetaEH_b[b]
                << " => peakDiff: " << peakDiff << ", relativeDip: "
                << jkRelDipEH_b[b] << endl;
        double cut = histogramMinimumLocation(fEH.getHistogram(), tolerance);
        ObsHistogramWeightDiffMinCallableJK fEW(tolerance, cut, numBins, this, b);
        double weightDiff = -1.;
        brentMinimize(jkBetaEW_b[b], weightDiff, fEH, betaStart, betaEnd);
        jkRelDipEW_b[b] = histogramRelativeDip(fEH.getHistogram(), tolerance);
        out << "final equal weight, b:" << b << ", beta: " << jkBetaEW_b[b]
                << " => weightDiff: " << weightDiff << ", relativeDip: "
                << jkRelDipEW_b[b] << endl;
    }
    double jkRelDipBlockAverageEH;
    double jkBetaBlockAverageEH;
    jackknife(jkRelDipBlockAverageEH, relDipErrorEH, jkRelDipEH_b);
    jackknife(jkBetaBlockAverageEH, betaDoubleErrorEH, jkBetaEH_b);
    out << "EH averaged beta:   " << jkBetaBlockAverageEH << endl;
    out << "EH averaged relDip: " << jkRelDipBlockAverageEH << endl;
    double jkRelDipBlockAverageEW;
    double jkBetaBlockAverageEW;
    jackknife(jkRelDipBlockAverageEW, relDipErrorEW, jkRelDipEW_b);
    jackknife(jkBetaBlockAverageEW, betaDoubleErrorEW, jkBetaEW_b);
    out << "EW averaged beta:   " << jkBetaBlockAverageEW << endl;
    out << "EW averaged relDip: " << jkRelDipBlockAverageEW << endl;

    out << "Equal Height Final result: relative dip: " << relDipEH
        << ", error: " << relDipErrorEH << '\n'
        << "at beta=" << betaDoubleEH << ", error: "
        << betaDoubleErrorEH << endl;
    out << "Equal Weight Final result: relative dip: " << relDipEW
        << ", error: " << relDipErrorEW << '\n'
        << "at beta=" << betaDoubleEW << ", error: "
        << betaDoubleErrorEW << endl;
}




void MultireweightHistosPTJK::energyRelDip(
        double& relDip, double& relDipError, HistogramDouble*& histo,
        double targetBeta, double tolerance) {
    //whole data-set
    MultireweightHistosPT::energyRelDip(relDip,
            histo, targetBeta, tolerance);

    //jackknife blocking
    vector<double> jkRelDip_b(blockCount);
    #pragma omp parallel for
    for (signed b = 0; b < (signed)blockCount; ++b) {
        HistogramDouble* histojk = reweightEnergyHistogramJK(targetBeta, b);
        jkRelDip_b[b] = histogramRelativeDip(histojk, tolerance);
        destroy(histojk);
    }

    double jkRelDipBlockAverage;
    jackknife(jkRelDipBlockAverage, relDipError, jkRelDip_b);
}

void MultireweightHistosPTJK::obsRelDip(
        double& relDip, double& relDipError, HistogramDouble*& histo,
        double targetBeta, unsigned numBins, double tolerance) {
    //whole data-set
    MultireweightHistosPT::obsRelDip(relDip,
            histo, targetBeta, numBins, tolerance);

    //jackknife blocking
    vector<double> jkRelDip_b(blockCount);
    #pragma omp parallel for
    for (signed b = 0; b < (signed)blockCount; ++b) {
        HistogramDouble* histojk = reweightObservableHistogramJK(
                targetBeta, numBins, b);
        jkRelDip_b[b] = histogramRelativeDip(histojk, tolerance);
        destroy(histojk);
    }

    double jkRelDipBlockAverage;
    jackknife(jkRelDipBlockAverage, relDipError, jkRelDip_b);
}


//TODO: error estimation really should be put into a more general
//function
void MultireweightHistosPTJK::saveH_km_errors(const std::string & filename) {
    ofstream outHkmErrors(filename.c_str());
    double factor = double(blockCount - 1) / double(blockCount);
    for (unsigned k = 0; k < numReplicas; ++k) {
        for (unsigned m = 0; m < binCount; ++m) {
            //error estimation
            double sqDev = 0;
            double avgH_bkm = 0;
            for (unsigned b = 0; b < blockCount; ++b) {
                avgH_bkm += H_bkm[b][k][m];
            }
            avgH_bkm /= blockCount;
            for (unsigned b = 0; b < blockCount; ++b) {
                //the following original calculation was wrong:
                //sqDev += pow(H_bkm[b][k][m] - H_km[k][m], 2);
                //the totals of H_bkm are always smaller than the
                //totals of H_km !
                sqDev += pow(H_bkm[b][k][m] - avgH_bkm, 2);
            }
            double error = sqrt(factor * sqDev);

            outHkmErrors << error << '\t';
        }
        outHkmErrors << '\n';
    }
}


ReweightingResult MultireweightHistosPTJK::reweightDiscrete(double beta) {
//    //estimates from whole data set:
//    ReweightingResult result = MultireweightHistosPT::reweightDiscrete(beta);

    //estimates for jackknifed sub sets:
    vector<double> energy_b(blockCount);
    vector<double> heatCapacity_b(blockCount);

    //estimation on jacknifed sub-sets
    #pragma omp parallel for
    for (signed b = 0; b < (signed)blockCount; ++b) {
        //helper array:
        vector<LogVal> args(binCount);
        args[0] = lOmega_bm[b][0] * toLogValExp(-beta * U_m[0]);
        LogVal normalization = args[0];
        for (unsigned m = 1; m < binCount; ++m) {
            args[m] = lOmega_bm[b][m] * toLogValExp(-beta * U_m[m]);
            normalization += args[m];
        }
        double estEnergyNorm = 0;
        double estEnergySqNorm = 0;
        for (unsigned m = 0; m < binCount; ++m) {
            args[m] /= normalization;
            double prob = toDouble(args[m]);
            double energyNorm = U_m[m] / systemN;
            double energySqNorm = energyNorm*energyNorm;
            estEnergyNorm += prob * energyNorm;
            estEnergySqNorm += prob * energySqNorm;
        }
        energy_b[b] = estEnergyNorm;
        heatCapacity_b[b] = beta*beta * systemN *
            (estEnergySqNorm - estEnergyNorm*estEnergyNorm);
    }

    //total result: average of jackknife blocks
    ReweightingResult result;
    jackknife(result.energyAvg, result.energyError, energy_b);
    jackknife(result.heatCapacity, result.heatCapacityError, heatCapacity_b);

    return result;
}


double MultireweightHistosPTJK::reweightSpecificHeatDiscreteJK(
        double targetBeta, unsigned jkBlock) {
    //helper array:
    vector<LogVal> args(binCount);
    args[0] = lOmega_bm[jkBlock][0] * toLogValExp(-targetBeta * U_m[0]);
    LogVal normalization = args[0];
    for (unsigned m = 1; m < binCount; ++m) {
        args[m] = lOmega_bm[jkBlock][m] * toLogValExp(-targetBeta * U_m[m]);
        normalization += args[m];
    }
    double estEnergyNorm = 0;
    double estEnergySqNorm = 0;
    for (unsigned m = 0; m < binCount; ++m) {
        args[m] /= normalization;
        double prob = toDouble(args[m]);
        double energyNorm = U_m[m] / systemN;
        double energySqNorm = energyNorm*energyNorm;
        estEnergyNorm += prob * energyNorm;
        estEnergySqNorm += prob * energySqNorm;
    }
    return targetBeta*targetBeta * systemN *
            (estEnergySqNorm - estEnergyNorm*estEnergyNorm);
}

