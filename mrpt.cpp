//Code by Max Henner Gerlach, 2010--2012,
//used for the diploma thesis "Directional Ordering in the Classical Compass Model in Two and Three Dimensions"
//contact: maxgerlach@gmail.com

/*
 * multireweighthisto-pt.cpp
 *
 *  Created on: Jun 28, 2011
 *      Author: gerlach
 */

// generalized for SDW DQMC (2015-02-06 - )


#include <algorithm>
#include <omp.h>
#include <dlib/misc_api.h>          //directory handling
#include "tools.h"
#include "mrpt.h"
#include "dataseriesloader.h"
#include "dataserieswriter.h"
#include "datamapwriter.h"
#include "statistics.h"
#include "numerics.h"

using namespace std;

template void printVector<LogVal>(const vector<LogVal>&);
template void printVector<double>(const vector<double>&);

MultireweightHistosPT::MultireweightHistosPT(ostream& outStream) :
        numReplicas(0), systemN(0),
        minEnergyNormalized(0), maxEnergyNormalized(0), binCount(0), binSize(0), lBinSize(1.0),
        basicConfig(false), out(outStream),
{
    out << "max threads: " << omp_get_max_threads() << endl;
}

MultireweightHistosPT::~MultireweightHistosPT() {
    destroyAll(energyTimeSeries);
    destroyAll(observableTimeSeries);
    destroyAll(betaIndexTimeSeries);
    destroyAll(m_kn);
}

void MultireweightHistosPT::addSimulationInfo(const std::string& filename) {
    DataSeriesLoader<double> info;
    info.readFromFile(filename);
    numReplicas = 0;
    info.getMeta("numTemperatures", numReplicas);
    systemN = 0;
    info.getMeta("N", systemN);
    systemL = 0;
    info.getMeta("L", systemL);
    infoNumSamples = 0;
    try {
        info.getMeta("sweepsBetweenMeasurements", infoSweepsBetweenMeasurements);
    } catch (KeyUndefined& exc) {
        infoSweepsBetweenMeasurements = 1;
    }
    try {
        info.getMeta("curSamples", infoNumSamples);
    } catch (KeyUndefined& exc) {
        try {
            info.getMeta("totalSamples", infoNumSamples);
        } catch (KeyUndefined& exc) {
            try {
                info.getMeta("samples", infoNumSamples);
            } catch (KeyUndefined& exc) {
                try {
                    info.getMeta("totalSweeps", infoNumSamples);
                } catch (KeyUndefined& exc) {
                    infoNumSamples = 0;
                }
            }
        }
    }

    info.deleteData();

    out << "Simulation info loaded from " << filename << " numReplicas=" << numReplicas << " systemN=" << systemN << endl;

    addedEnergyTimeSeries = 0;
    addedObservableTimeSeries = 0;
    energyTimeSeries.resize(numReplicas, 0);
    observableTimeSeries.resize(numReplicas, 0);
    betaIndexTimeSeries.resize(numReplicas, 0);
    lZ_l.resize(numReplicas, LogVal(1.0));

    betas.loadBetas(filename);
    out << "betas: ";
    for (unsigned b = 0; b < betas.getSize(); ++b) {
        out << betas[b] << " ";
    }
    out << endl;
    basicConfig = true;
}


void MultireweightHistosPT::addInputTimeSeries(const std::string& filename, unsigned subsample, unsigned discardEntries) {
    if (not basicConfig) {
        throw GeneralException("tried to add input time series before basic configuration from simulation info was done (need to parse info.dat first)");
    }
    DoubleSeriesLoader input;
    out << "adding time-series " << filename << " - " << flush;
    input.readFromFile(filename, subsample, discardEntries, infoNumSamples);
    unsigned replicaIndex;
    input.getMeta("r", replicaIndex);
    string obs;
    input.getMeta("observable", obs);
    out << obs << " - " << input.getData(0)->size() << endl;
    if (obs == "energy") {
//      if (energyTimeSeries.size() < replicaIndex + 1) {
//          //content preserving resizes:
//          unsigned newNumReplicas = replicaIndex + 1;
//          energyTimeSeries.resize(newNumReplicas, 0);
//          betaIndexTimeSeries.resize(newNumReplicas, 0);
//      }
        if (energyTimeSeries[replicaIndex] != 0) {
            throw GeneralException("two time series added for " + obs + " from replica-index " + numToString(replicaIndex));
        }
        energyTimeSeries[replicaIndex] = input.getData(input.getColumns() - 1);
        if (input.getColumns() != 2) {
            throw GeneralException(filename + ": expected a two column time series file, but got " + numToString(input.getColumns()) + " column(s)");
        }
        ++addedEnergyTimeSeries;
        //update betaIndexTimeSeries:
        std::vector<double>* dBetaIndexTimeSeries = input.getData(0);           //values loaded in as doubles...
        betaIndexTimeSeries[replicaIndex] = new vector<int>(energyTimeSeries[replicaIndex]->size());
        for (unsigned n = 0; n < dBetaIndexTimeSeries->size(); ++n) {
            int betaIndex = static_cast<int>((*dBetaIndexTimeSeries)[n]);       //TODO:useless cast!
            (*(betaIndexTimeSeries[replicaIndex]))[n] = betaIndex;
        }
        delete dBetaIndexTimeSeries;        //don't need this in memory anymore
    } else if (observable == "" or obs == observable) {
        observable = obs;
//      if (observableTimeSeries.size() < replicaIndex + 1) {
//          observableTimeSeries.resize(replicaIndex + 1, 0);
//      }
        if (observableTimeSeries[replicaIndex] != 0) {
            throw GeneralException("two time series added for " + obs + " from replica-index " + numToString(replicaIndex));
        }
        observableTimeSeries[replicaIndex] = input.getData(input.getColumns() - 1);
        ++addedObservableTimeSeries;
        if (input.getColumns() == 2) {
            delete input.getData(0);        //don't need this in memory
        }
    } else {
        throw GeneralException("in " + filename + ": Observable doesn't match previous\n" +
                obs + " vs. " + observable + "\n");
    }
    unsigned N;
    input.getMeta("N", N);
    if (systemN == 0) {
        systemN = N;
    } else if (systemN != N) {
        throw GeneralException("in " + filename + ": Non matching system sizes: " + numToString(N) + " vs. " + numToString(systemN));
    }
}

void MultireweightHistosPT::sortTimeSeriesByBeta() {
    out << "Sorting time series by beta..." << flush;

    if (addedObservableTimeSeries > 0) {
        //all time series need to have the same length
        unsigned M = energyTimeSeries[0]->size();
        for (unsigned r = 1; r < numReplicas; ++r) {
            if (energyTimeSeries[r]->size() != M or
                    observableTimeSeries[r]->size() != M) {
                throw GeneralException("Time series length mismatch: for beta-sorted "
                        "timeseries all need to have the same length");
            }
        }
        //for each sample:
        //first these vectors hold the original ordering by replica index
        //then they are sorted back into the time series, sorted by beta index
        vector<double> tempEnergies(numReplicas);
        vector<double> tempObs(numReplicas);
        vector<int> tempBetaIndices(numReplicas);
        for (unsigned m = 0; m < M; ++m) {
            for (unsigned r = 0; r < numReplicas; ++r) {
                tempEnergies[r] = (*energyTimeSeries[r])[m];
                tempObs[r] = (*observableTimeSeries[r])[m];
                tempBetaIndices[r] = (*betaIndexTimeSeries[r])[m];
            }
            for (unsigned r = 0; r < numReplicas; ++r) {
                int bi = tempBetaIndices[r];
                (*energyTimeSeries[bi])[m] = tempEnergies[r];
                (*observableTimeSeries[bi])[m] = tempObs[r];
                //the following is now trivial, but later code depends on it
                (*betaIndexTimeSeries[bi])[m] = bi;
            }
        }
    } else {
        //all time series need to have the same length
        unsigned M = energyTimeSeries[0]->size();
        for (unsigned r = 1; r < numReplicas; ++r) {
            if (energyTimeSeries[r]->size() != M) {
                throw GeneralException("Time series length mismatch: for beta-sorted "
                        "timeseries all need to have the same length");
            }
        }
        //for each sample:
        //first these vectors hold the original ordering by replica index
        //then they are sorted back into the time series, sorted by beta index
        vector<double> tempEnergies(numReplicas);
        vector<int> tempBetaIndices(numReplicas);
        for (unsigned m = 0; m < M; ++m) {
            for (unsigned r = 0; r < numReplicas; ++r) {
                tempEnergies[r] = (*energyTimeSeries[r])[m];
                tempBetaIndices[r] = (*betaIndexTimeSeries[r])[m];
            }
            for (unsigned r = 0; r < numReplicas; ++r) {
                int bi = tempBetaIndices[r];
                (*energyTimeSeries[bi])[m] = tempEnergies[r];
                //the following is now trivial, but later code depends on it
                (*betaIndexTimeSeries[bi])[m] = bi;
            }
        }
    }
    out << endl;
}

MultireweightHistosPT::ResultsMap* MultireweightHistosPT::
        directNoReweighting() {
    out << "Computing estimates from time series without any reweighting... "
        << flush;

    ResultsMap* results = new ResultsMap;

    //vectors betaIndex (of this->betas) -> value at that temperature
    //everything normalized by system volume
    vector<double> meanEnergy_l(numReplicas, 0);
    vector<double> meanEnergySquared_l(numReplicas, 0);
    vector<double> meanObs_l(numReplicas, 0);
    vector<double> meanObsSquared_l(numReplicas);
    vector<double> meanObsToTheFourth_l(numReplicas, 0);

    //calculate for each replica separately:
    #pragma omp parallel for
    for (signed k = 0; k < (signed)numReplicas; ++k) {
        //vectors for this replica (still map betaindex -> values):
        vector<double> k_meanEnergy_l(numReplicas, 0);
        vector<double> k_meanEnergySquared_l(numReplicas, 0);
        vector<double> k_meanObs_l(numReplicas, 0);
        vector<double> k_meanObsSquared_l(numReplicas);
        vector<double> k_meanObsToTheFourth_l(numReplicas, 0);
        unsigned N_k = energyTimeSeries[k]->size();
        //sum up:
        for (unsigned n = 0; n < N_k; ++n) {
            int bi = (*betaIndexTimeSeries[k])[n];
            double e = (*energyTimeSeries[k])[n];
            double o = (*observableTimeSeries[k])[n];
            k_meanEnergy_l[bi] += e;
            k_meanEnergySquared_l[bi] += e*e;
            k_meanObs_l[bi] += o;
            k_meanObsSquared_l[bi] += o*o;
            k_meanObsToTheFourth_l[bi] += o*o*o*o;
        }
        //divide to form averages, divisor: no of samples for this temperature
        //in the time series (numReplicas == #temperatures)
        double divisor = double(N_k) / double(numReplicas);
        for (unsigned bi = 0; bi < numReplicas; ++bi) {
            k_meanEnergy_l[bi] /= divisor;
            k_meanEnergySquared_l[bi] /= divisor;
            k_meanObs_l[bi] /= divisor;
            k_meanObsSquared_l[bi] /= divisor;
            k_meanObsToTheFourth_l[bi] /= divisor;
        }
        //average over replicas -- sum up everything from all replicas
        #pragma omp critical
        {
            for (unsigned bi = 0; bi < numReplicas; ++bi) {
                meanEnergy_l[bi] += k_meanEnergy_l[bi];
                meanEnergySquared_l[bi] += k_meanEnergySquared_l[bi];
                meanObs_l[bi] += k_meanObs_l[bi];
                meanObsSquared_l[bi] += k_meanObsSquared_l[bi];
                meanObsToTheFourth_l[bi] +=
                        k_meanObsToTheFourth_l[bi];
            }
        }
    }
    //average over replicas -- divide
    for (unsigned bi = 0; bi < numReplicas; ++bi) {
        meanEnergy_l[bi] /= double(numReplicas);
        meanEnergySquared_l[bi] /= double(numReplicas);
        meanObs_l[bi] /= double(numReplicas);
        meanObsSquared_l[bi] /= double(numReplicas);
        meanObsToTheFourth_l[bi] /= double(numReplicas);
    }
    //compute final results and put them into the map
    for (unsigned bi = 0; bi < numReplicas; ++bi) {
        double beta = betas[bi];
        double heatCapacity = systemN * beta*beta *
                (meanEnergySquared_l[bi] - pow(meanEnergy_l[bi], 2));
        double suscObs = systemN *
                (meanObsSquared_l[bi] - pow(meanObs_l[bi], 2));
        double binderObs = 1.0 - (meanObsToTheFourth_l[bi] /
                        (3 * pow(meanObsSquared_l[bi], 2)));
        //set results without specification of errors:
        (*results)[beta] = ReweightingResult(meanEnergy_l[bi],
                heatCapacity,
                meanObs_l[bi],
                suscObs,
                binderObs);
    }
    out << "done." << endl;
    return results;
}

void MultireweightHistosPT::setUpHistograms(int binCount_) {
    if (addedEnergyTimeSeries == 0) {
        throw GeneralException("No energy time series added!");
    }
    if (addedObservableTimeSeries == 0) {
        cout << "Warning: No observable time series added." << endl;
    }
    if (numReplicas != betas.getSize() or numReplicas != addedEnergyTimeSeries
            or (addedObservableTimeSeries > 0 and
                numReplicas != addedObservableTimeSeries)) {
//        throw GeneralException("Mismatched number of replicas (forgot to add time series? wrong betas?)");
        cout << "Warning: Number of added time series does not match number "
                "of replicas from simulation info file." << endl;
    }
    //replace null-pointer with zero-length vectors at the spots, where
    //no time series were added
    for (unsigned k = 0; k < numReplicas; ++k) {
        if (energyTimeSeries[k] == 0) {
            energyTimeSeries[k] = new vector<double>(0);
        }
        if (observableTimeSeries[k] == 0) {
            observableTimeSeries[k] = new vector<double>(0);
        }
        if (betaIndexTimeSeries[k] == 0) {
            betaIndexTimeSeries[k] = new vector<int>(0);
        }
    }

    binCount = binCount_;

    findMinMax(energyTimeSeries, minEnergyNormalized, maxEnergyNormalized);
    minEnergy = minEnergyNormalized * systemN;
    maxEnergy = maxEnergyNormalized * systemN;
    deltaU = (maxEnergyNormalized - minEnergyNormalized) / binCount;
    const double SMALL = 1e-10;     //to fit maxEnergy into the highest bin
    binSize = (maxEnergyNormalized - minEnergyNormalized + SMALL) / binCount;
    lBinSize = LogVal(binSize);

    if (addedObservableTimeSeries > 0) {
        findMinMax(observableTimeSeries, minObservableNormalized, maxObservableNormalized);
    }

    U_m.resize(binCount);
    for (unsigned m = 0; m < binCount; ++m) {
        //take normalized energy at bin center, then undo normalization
        U_m[m] = (minEnergyNormalized + binSize * (m + 0.5)) * systemN;
    }

    lZ_l = vector<LogVal>(numReplicas, LogVal(1.0));

    H_km.resize(boost::extents[numReplicas][binCount]);
    initArray(H_km, 0);
    g_km.resize(boost::extents[numReplicas][binCount]);
    initArray(g_km, 0);
    H_lm.resize(boost::extents[numReplicas][binCount]);
    initArray(H_lm, 0);
    N_kl.resize(boost::extents[numReplicas][numReplicas]);
    initArray(N_kl, 0);

    H_m = Histo(binCount, 0);

    m_kn = IntSeriesCollection(numReplicas, 0);
    for (unsigned k = 0; k < numReplicas; ++k) {
        m_kn[k] = new std::vector<int>(energyTimeSeries[k]->size());
    }

    lOmega_m.resize(binCount, LogVal(1.0));
}

void MultireweightHistosPT::setUpHistogramsIsing() {
    if (addedEnergyTimeSeries == 0) {
        throw GeneralException("No energy time series added!");
    }
    if (addedObservableTimeSeries == 0) {
        cout << "Warning: No observable time series added." << endl;
    }
    if (numReplicas != betas.getSize() or numReplicas != addedEnergyTimeSeries
            or (addedObservableTimeSeries > 0 and
                numReplicas != addedObservableTimeSeries)) {
//        throw GeneralException("Mismatched number of replicas (forgot to add time series? wrong betas?)");
        cout << "Warning: Number of added time series does not match number "
                "of replicas from simulation info file." << endl;
    }
    //replace null-pointer with zero-length vectors at the spots, where
    //no time series were added
    for (unsigned k = 0; k < numReplicas; ++k) {
        if (energyTimeSeries[k] == 0) {
            energyTimeSeries[k] = new vector<double>(0);
        }
        if (observableTimeSeries[k] == 0) {
            observableTimeSeries[k] = new vector<double>(0);
        }
        if (betaIndexTimeSeries[k] == 0) {
            betaIndexTimeSeries[k] = new vector<int>(0);
        }
    }

    findMinMax(energyTimeSeries, minEnergyNormalized, maxEnergyNormalized);

    //actual min and max energies:
    minEnergy = round(minEnergyNormalized * double(systemN));
    maxEnergy = round(maxEnergyNormalized * double(systemN));

    //special knowledge of 2D Ising model:
    //offset by half a bin width, else the binning routines won't work due to
    //rounding errors:
//    minEnergyNormalized -= double(2) / systemN;
    //the energies -2N+4 and 2N-4 aren't physically accessible, but we'll
    //include them in the histograms anyway
    binCount = int(maxEnergy - minEnergy) / 4 + 1;

    deltaU = 4.0;
    const double SMALL = 1e-10;     //to fit maxEnergy into the highest bin
    binSize = deltaU / double(systemN) + SMALL;
    lBinSize = LogVal(binSize);

    if (addedObservableTimeSeries) {
        findMinMax(observableTimeSeries, minObservableNormalized, maxObservableNormalized);
    }

    U_m.resize(binCount);
    for (int m = 0; m < int(binCount); ++m) {
          U_m[m] = double(minEnergy) + m * deltaU;
    }

    lZ_l = vector<LogVal>(numReplicas, LogVal(1.0));

    H_km.resize(boost::extents[numReplicas][binCount]);
    initArray(H_km, 0);
    g_km.resize(boost::extents[numReplicas][binCount]);
    initArray(g_km, 0);
    H_lm.resize(boost::extents[numReplicas][binCount]);
    initArray(H_lm, 0);
    N_kl.resize(boost::extents[numReplicas][numReplicas]);
    initArray(N_kl, 0);

    H_m = Histo(binCount, 0);

    m_kn = IntSeriesCollection(numReplicas, 0);
    for (unsigned k = 0; k < numReplicas; ++k) {
        m_kn[k] = new std::vector<int>(energyTimeSeries[k]->size());
    }

    lOmega_m.resize(binCount, LogVal(1.0));
}

void MultireweightHistosPT::createHistogramsHelper() {
    out << "Creating energy histograms etc. minEnergyNormalized=" << minEnergyNormalized
            << " maxEnergyNormalized=" << maxEnergyNormalized
            << " binSize=" << binSize
            << " binCount=" << binCount << endl;

//  #pragma omp parallel for
    for (int k = 0; k < (signed)numReplicas; ++k) {
        for (unsigned n = 0; n < energyTimeSeries[k]->size(); ++n) {
            //get bin number, correct due to rounding in cast:
            int m = static_cast<int>(((*energyTimeSeries[k])[n] - minEnergyNormalized) / binSize);
            //get beta index
            int l = (*(betaIndexTimeSeries[k]))[n];

            N_kl[k][l] += 1;

            (*m_kn[k])[n] = m;
            ++H_km[k][m];
//          #pragma omp atomic
            ++H_m[m];
//          #pragma omp atomic
            ++H_lm[l][m];
        }
        out << "." << flush;
    }

    destroyAll(betaIndexTimeSeries);                    //not needed any more

    out << " done" << endl;
}

void MultireweightHistosPT::createHistogramsHelperDiscrete() {
    out << "Creating energy histograms etc. minEnergyNormalized=" << minEnergyNormalized
            << " maxEnergyNormalized=" << maxEnergyNormalized
            << " binSize=" << binSize
            << " binCount=" << binCount
            << " for originally discrete bins " << endl;

//  #pragma omp parallel for
    for (int k = 0; k < (signed)numReplicas; ++k) {
        for (unsigned n = 0; n < energyTimeSeries[k]->size(); ++n) {
            //get bin number, correct due to rounding in cast:
//            int m = static_cast<int>(((*energyTimeSeries[k])[n] - minEnergyNormalized) / binSize);
            int m = int(round(((*energyTimeSeries[k])[n] * systemN - minEnergy) / deltaU));
            //get beta index
            int l = (*(betaIndexTimeSeries[k]))[n];

            N_kl[k][l] += 1;

            (*m_kn[k])[n] = m;
            ++H_km[k][m];
//          #pragma omp atomic
            ++H_m[m];
//          #pragma omp atomic
            ++H_lm[l][m];
        }
        out << "." << flush;
    }

    destroyAll(betaIndexTimeSeries);                    //not needed any more

    out << " done" << endl;
}

void MultireweightHistosPT::createHistograms(int binCount_) {
    setUpHistograms(binCount_);
    createHistogramsHelper();
}

void MultireweightHistosPT::createHistogramsIsing() {
    setUpHistogramsIsing();
    createHistogramsHelperDiscrete();
}


void MultireweightHistosPT::measureGlobalInefficiencies(bool saveAutocorr) {
    out << "Measuring energy time series inefficiencies g_k, "
           "ignoring differences in individual bins" << endl;

    g_km.resize(boost::extents[numReplicas][binCount]);

    if (saveAutocorr) {
        dlib::create_directory("autocorr");
        dlib::set_current_dir("autocorr");
    }

    #pragma omp parallel for
    for (int k = 0; k < (signed)numReplicas; ++k) {
        AutoCorrMap* autocorr = 0;
        if (saveAutocorr) {
            autocorr = new AutoCorrMap;
        }
        double tauint = tauint_adaptive(energyTimeSeries[k], autocorr);
        double g_k = 1 + 2 * tauint;
        for (unsigned m = 0; m < binCount; ++m) {
            g_km[k][m] = g_k;
        }
        if (saveAutocorr) {
            IntDoubleMapWriter writeAutocorr;
            writeAutocorr.setData(autocorr);
            writeAutocorr.addHeaderText("Autocorrelation function of the energy");
            writeAutocorr.addHeaderText("Replica-index (corresponds to beta-index if time series were sorted):");
            writeAutocorr.addMeta("k", k);
            writeAutocorr.addHeaderText("Estimated integrated autocorrelation time");
            writeAutocorr.addMeta("tau-int", tauint);
            writeAutocorr.addMeta("g", g_k);
            writeAutocorr.addHeaderText("gap t \t autocorr");
            writeAutocorr.writeToFile("autocorr-k" + numToString(k) + ".dat");
            destroy(autocorr);
        }
        out << "." << flush;
    }

    if (saveAutocorr) {
        dlib::set_current_dir("..");
    }

    out << " done" << endl;
}

void MultireweightHistosPT::writeOutEnergyTauInt(const std::string& filename) {
    out << "Estimating energy tau-ints... " << endl;
    map<double, double> tauint;
    #pragma omp parallel for
    for (int k = 0; k < (signed)numReplicas; ++k) {
        tauint.insert(make_pair(betas[k],
                tauint_adaptive(energyTimeSeries[k]) *
                infoSweepsBetweenMeasurements));
        cout << "." << flush;
    }
    DoubleMapWriter w;
    w.addHeaderText("Integrated autocorrelation time of the energy, used fast adaptive method");
    w.addHeaderText("Unit: Monte Carlo Sweeps (even if less measured samples were taken)");
    w.addMeta("sweepsBetweenMeasurements", infoSweepsBetweenMeasurements);
    w.addMeta("observable", "energy");
    w.addMeta("L", systemL);
    w.addMeta("N", systemN);
    w.addHeaderText("beta \t tau-int");
    w.setData(&tauint);
    w.writeToFile(filename);
    out << " Done." << endl;
}

void MultireweightHistosPT::writeOutObsTauInt(const std::string& filename) {
    out << "Estimating " << observable << " tau-ints... " << endl;
    map<double, double> tauint;
    #pragma omp parallel for
    for (int k = 0; k < (signed)numReplicas; ++k) {
        tauint.insert(make_pair(betas[k],
                tauint_adaptive(observableTimeSeries[k]) *
                infoSweepsBetweenMeasurements));
        cout << "." << flush;
    }
    DoubleMapWriter w;
    w.addHeaderText("Integrated autocorrelation time of " + observable +
            ", used fast adaptive method");
    w.addHeaderText("Unit: Monte Carlo Sweeps (even if less measured samples were taken)");
    w.addMeta("sweepsBetweenMeasurements", infoSweepsBetweenMeasurements);
    w.addMeta("observable", observable);
    w.addMeta("L", systemL);
    w.addMeta("N", systemN);
    w.addHeaderText("beta \t tau-int");
    w.setData(&tauint);
    w.writeToFile(filename);
    out << " Done." << endl;
}


void MultireweightHistosPT::setBinInefficienciesToUnity() {
    out << "Setting all bin inefficiencies g_km to 1";

    g_km.resize(boost::extents[numReplicas][binCount]);
    std::fill(g_km.data(), g_km.data() + g_km.num_elements(), 1.0);

    out << endl;
}


void MultireweightHistosPT::measureBinInefficiencies(bool saveAutocorr) {
    out << "Measuring bin inefficiencies g_km" << endl;

    g_km.resize(boost::extents[numReplicas][binCount]);
    std::fill(g_km.data(), g_km.data() + g_km.num_elements(), 1.0);

    if (saveAutocorr) {
        dlib::create_directory("autocorr");
        dlib::set_current_dir("autocorr");
        for (unsigned k = 0; k < numReplicas; ++k) {
            if (energyTimeSeries[k]->size() > 0) {
                string dir = "k" + numToString(k);
                dlib::create_directory(dir);
            }
        }
    }

    #pragma omp parallel for
    for (int k = 0; k < (signed)numReplicas; ++k) {
        vector<AutoCorrMap*> autocorr_m(binCount, 0);
        if (saveAutocorr) {
            for (unsigned m = 0; m < binCount; ++m) {
                autocorr_m[m] = new AutoCorrMap;
            }
        }

        //adaptive integration of autocorrelation function
        //high time resolution for small lag times, lower resolution for higher times in the
        //slowly decaying tail of the autocorrelation function
        vector<char> zeroCrossed(binCount, false);
        unsigned countZeroCrossed = 0;

        unsigned N = energyTimeSeries[k]->size();
        if (N == 0) {
            //skip empty time series
            continue;
        }

        vector<double> binMeanOccupation(binCount, 0);
        vector<double> binMeanSquared(binCount, 0);
        for (unsigned m = 0; m < binCount; ++m) {
            binMeanOccupation[m] = double(H_km[k][m]) / double(N);
            binMeanSquared[m] = pow(binMeanOccupation[m], 2);

            if (binMeanOccupation[m] == 0) {
                //bin never occupied in this replica --> g_km has no real sense
                //==> don't take this bin into account at all in the following, leave g_km at 1.0
                zeroCrossed[m] = true;
                ++countZeroCrossed;
            }
        }

        unsigned i = 1;
        unsigned t_i = 1;
        do {
            //compute bin autocorrelation functions
            vector<double> binAutoCorr(binCount, 0);
            for (unsigned n = 0; n < N - t_i; ++n) {
                int bin1 = (*m_kn[k])[n];
                int bin2 = (*m_kn[k])[n + t_i];
                if (bin1 == bin2) {
                    binAutoCorr[bin1] += 1.0;
                }
            }
            for (unsigned m = 0; m < binCount; ++m) {
                if (not zeroCrossed[m]) {
                    binAutoCorr[m] /= double(N - t_i);
                    binAutoCorr[m] -= binMeanSquared[m];
                    binAutoCorr[m] /= (binMeanOccupation[m] - binMeanSquared[m]);   //normalize by variance
                    if (binAutoCorr[m] <= 0) {
                        zeroCrossed[m] = true;
                        ++countZeroCrossed;
                    } else {
                        //weighted addition to estimate integrated autocorrelation time
                        double t_next = 1.0 + (i + 1.0) * (i) / 2.0;
                        g_km[k][m] += 2.0 * binAutoCorr[m] * (t_next - t_i);
                        if (saveAutocorr) {
                            autocorr_m[m]->insert(make_pair(t_i, binAutoCorr[m]));
                        }
                    }
                }
            }

            i += 1;
            //lag time for next step of the iteration
            t_i = 1.0 + i * (i - 1.0) / 2.0;
        } while (t_i < N - 1 and countZeroCrossed < binCount);

        if (saveAutocorr) {
            for (unsigned m = 0; m < binCount; ++m) {
                IntDoubleMapWriter wa;
                wa.setData(autocorr_m[m]);
                wa.addHeaderText("Autocorrelation function of the characteristic function of one energy bin");
                wa.addHeaderText("Replica-index (corresponds to beta-index if time series were sorted):");
                wa.addMeta("k", k);
                wa.addHeaderText("Bin-index");
                wa.addMeta("m", m);
                wa.addHeaderText("Corresponding energy");
                wa.addMeta("U_m", U_m[m]);
                wa.addHeaderText("Corresponding energy normalized by system N");
                wa.addMeta("e_m", double(U_m[m]) / double(systemN));
                wa.addHeaderText("Statistical inefficiency estimated from integrated autocorrelation time");
                wa.addMeta("g_km", g_km[k][m]);
                wa.addHeaderText("gap t \t autocorr");
                wa.writeToFile("k" + numToString(k) + "/autocorr-k" + numToString(k)
                        + "m-" + numToString(m) + ".dat");
                destroy(autocorr_m[m]);
            }
        }

        out << "." << flush;
    }

    if (saveAutocorr) {
        dlib::set_current_dir("..");
    }

    out << " done" << endl;
}

void MultireweightHistosPT::saveg_km(const std::string& filename) {
    ofstream outgkm(filename.c_str());
    for (unsigned k = 0; k < numReplicas; ++k) {
        for (unsigned m = 0; m < binCount; ++m) {
            outgkm << g_km[k][m] << '\t';
        }
        outgkm << '\n';
    }
}

void MultireweightHistosPT::saveH_km(const std::string& filename) {
    ofstream outHkm(filename.c_str());
    for (unsigned k = 0; k < numReplicas; ++k) {
        for (unsigned m = 0; m < binCount; ++m) {
            outHkm << H_km[k][m] << '\t';
        }
        outHkm << '\n';
    }
}


void MultireweightHistosPT::saveU_m(const std::string& filename) {
    ofstream outUm(filename.c_str());
    for (unsigned m = 0; m < binCount; ++m) {
        outUm << U_m[m] << '\n';
    }
}


void MultireweightHistosPT::saveNeff_lm(const std::string& filename) {
    ofstream output(filename.c_str());
    for (unsigned l = 0; l < numReplicas; ++l) {
        for (unsigned m = 0; m < binCount; ++m) {
            output << Neff_lm[l][m] << '\t';
        }
        output << '\n';
    }
}

void MultireweightHistosPT::saveNeff_l(const std::string& filename) {
    ofstream output(filename.c_str());
    for (unsigned l = 0; l < numReplicas; ++l) {
        output << Neff_l[l] << '\n';
    }
}

void MultireweightHistosPT::savelNeff_lm(const std::string& filename) {
    ofstream output(filename.c_str());
    for (unsigned l = 0; l < numReplicas; ++l) {
        for (unsigned m = 0; m < binCount; ++m) {
            output << lNeff_lm[l][m] << '\t';
        }
        output << '\n';
    }
}



inline void MultireweightHistosPT::updateEffectiveCounts() {
    Heff_m.resize(binCount, 0);
    Neff_lm.resize(boost::extents[numReplicas][binCount]);
    initArray(Neff_lm, 0);
    for (unsigned k = 0; k < numReplicas; ++k) {
        for (unsigned m = 0; m < binCount; ++m) {
            Heff_m[m] += double(H_km[k][m]) / g_km[k][m];
            for (unsigned l = 0; l < numReplicas; ++l) {
                Neff_lm[l][m] += double(N_kl[k][l]) / g_km[k][m];
            }
        }
    }
    lHeff_m.resize(binCount);
    lNeff_lm.resize(boost::extents[numReplicas][binCount]);
    lPrecalc_lm.resize(boost::extents[numReplicas][binCount]);
    Neff_l.resize(numReplicas, 0);
    for (unsigned m = 0; m < binCount; ++m) {
        lHeff_m[m] = (Heff_m[m] != 0 ? LogVal(Heff_m[m]) : LogVal(LogVal::LogZero));
        for (unsigned l = 0; l < numReplicas; ++l) {
            lNeff_lm[l][m] = (Neff_lm[l][m] != 0 ? LogVal(Neff_lm[l][m]) : LogVal(LogVal::LogZero));
            lPrecalc_lm[l][m] = lNeff_lm[l][m] * lBinSize * toLogValExp(-betas[l] * U_m[m]);    //for updateDensityOfStates
            Neff_l[l] += Neff_lm[l][m];
        }
    }
}

inline void MultireweightHistosPT::updateDensityOfStates() {
    //precalculated the part that does not depend on the estimates of the partition functions
    //lPrecalc_lm[l][m] = lNeff_lm[l][m] * lBinSize * toLogValExp(-betas[l] * U_m[m]);
    for (int m = 0; m < (signed)binCount; ++m) {
        lOmega_m[m] = lHeff_m[m];
        LogVal denominator = lPrecalc_lm[0][m] / lZ_l[0];
        for (unsigned l = 1; l < numReplicas; ++l) {
            denominator += lPrecalc_lm[l][m] / lZ_l[l];
        }
        lOmega_m[m] /= denominator;
    }
}

void MultireweightHistosPT::findDensityOfStatesNonIteratively() {
    out << "Non-iterative estimate of the density of states... " << flush;

    const double deltaU = binSize * systemN;
    LogVal2Array x_lm(boost::extents[numReplicas][binCount - 1]);       // <-> difference of microcanonical entropies at bins m, m+1
    Double2Array w_lm(boost::extents[numReplicas][binCount - 1]);       //weights at temperature l, energy bin m
    vector<double> w_m(binCount - 1, 0);                                //sum_l { w_lm }
    for (unsigned l = 0; l < numReplicas; ++l) {
        double temperatureExponent = +betas[l] * deltaU;
        for (unsigned m = 0; m < binCount - 1; ++m) {
//          double curHist = max(H_lm[l][m], 1);
//          double nextHist = max(H_lm[l][m + 1], 1);
            double curHist = H_lm[l][m];
            double nextHist = H_lm[l][m + 1];
            if (curHist > 0 and nextHist > 0) {
//              x_lm[l][m].lnx = std::log(nextHist) - std::log(curHist) + temperatureInfluence;
                x_lm[l][m] = (LogVal(nextHist * deltaU) / LogVal(curHist * deltaU)) * toLogValExp(temperatureExponent);
                w_lm[l][m] = curHist * nextHist / (curHist + nextHist);
                w_m[m] += w_lm[l][m];
            } else {
                x_lm[l][m] = LogVal();          //very small
                w_lm[l][m] = 0;
            }
        }
    }

    lOmega_m[0] = LogVal(1.0);
//  LogVal normalization = lOmega_m[0];
    for (unsigned m = 0; m < binCount - 1; ++m) {
//      LogVal x_m;         //very small
        double ln_x_m = 0;
        for (unsigned l = 0; l < numReplicas; ++l) {
//          cout << w_lm[l][m] << '\t' << deltaS_lm[l][m] << endl;
//          LogVal add;
//          add = pow(x_lm[l][m], w_lm[l][m]);
//          x_m += add;
            if (w_m[m] > 0) ln_x_m += (w_lm[l][m] / w_m[m]) * x_lm[l][m].lnx;
        }
//      if (w_m[m] > 0) x_m /= w_m[m];
//      else x_m = LogVal();
//      cout << m << "\t" << lOmega_m[m] << '\t' << deltaS_m << "\n";

//      lOmega_m[m + 1] = lOmega_m[m] * x_m;
        lOmega_m[m + 1].lnx = lOmega_m[m].lnx + ln_x_m;
//      normalization += lOmega_m[m + 1];
    }
//  for (unsigned m = 0; m < binCount; ++m) {
//      lOmega_m[m] /= normalization;
//  }

    out << "Done" << endl;
}


ReweightingResult MultireweightHistosPT::reweightDiscrete(double targetBeta) {
    vector<LogVal> arguments(binCount);
    arguments[0] = lOmega_m[0] * toLogValExp(-targetBeta * U_m[0]);
    LogVal normalization = arguments[0];
    for (unsigned m = 1; m < binCount; ++m) {
        arguments[m] = lOmega_m[m] * toLogValExp(-targetBeta * U_m[m]);
        normalization += arguments[m];
    }

    double estEnergyNorm = 0;
    double estEnergySqNorm = 0;

    for (unsigned m = 0; m < binCount; ++m) {
        double energyNorm = U_m[m] / systemN;
        double energySqNorm = energyNorm*energyNorm;
        arguments[m] /= normalization;
        double prob = toDouble(arguments[m]);
        estEnergyNorm += energyNorm * prob;
        estEnergySqNorm += energySqNorm * prob;
    }

    ReweightingResult result;
    result.energyAvg = estEnergyNorm;
    result.heatCapacity = targetBeta*targetBeta * systemN *
            (estEnergySqNorm - estEnergyNorm*estEnergyNorm);

    return result;
}

void MultireweightHistosPT::findPartitionFunctionsAndDensityOfStates(double tolerance, int maxIterations) {
    out << "Updating effective counts... " << flush;
    updateEffectiveCounts();
    out << "Done." << endl;

    out << "Starting iteration to estimate density of states, tolerance=" << tolerance << " maxIterations=" << maxIterations << endl;

//  updateDensityOfStates();

    vector<LogVal> lZ_l_lastIteration = lZ_l;

    int iterations = 0;

    double deltaSquared = 0;    //this is sum_l { ((Z_l - Z_l_old) / (Z_l)) ** 2 }, stop iteration once deltaSquared < tolerance ** 2
    do {
        ++iterations;

        deltaSquared = 0;

        //update estimates of partition functions
        lZ_l[0] = lOmega_m[0] * lBinSize * toLogValExp(-betas[0] * U_m[0]);
        for (unsigned m = 1; m < binCount; ++m) {
            lZ_l[0] += lOmega_m[m] * lBinSize * toLogValExp(-betas[0] * U_m[m]);
        }
        #pragma omp parallel for reduction(+: deltaSquared)
        for (int l = 1; l < (signed)numReplicas; ++l) {
            lZ_l_lastIteration[l] = lZ_l[l];            //store old value
            lZ_l[l] = lOmega_m[0] * lBinSize * toLogValExp(-betas[l] * U_m[0]);
            for (unsigned m = 1; m < binCount; ++m) {
                lZ_l[l] += lOmega_m[m] * lBinSize * toLogValExp(-betas[l] * U_m[m]);
            }
            lZ_l[l] /= lZ_l[0];             //normalize

            deltaSquared += pow(expm1(lZ_l[l].lnx - lZ_l_lastIteration[l].lnx), 2);     //gauge change from last iteration
        }
        lZ_l[0].lnx = 0;            //sets to 1 (normalize)

        if (iterations % 10 == 0) {
            out << " Iteration " << iterations << " deltaSquared=" << deltaSquared << endl;
        }

        updateDensityOfStates();
    } while (iterations < maxIterations and deltaSquared >= tolerance * tolerance);
    out << "Done." << endl;
}

MultireweightHistosPT::DoubleSeriesCollection MultireweightHistosPT::computeWeights(double targetBeta) {
//    out << "Computing weights w_kn at beta=" << targetBeta << endl;
    out << " " << targetBeta << flush;

    vector<LogVal> arguments(binCount);
    arguments[0] = lOmega_m[0] * toLogValExp(-targetBeta * U_m[0]);
    LogVal normalization = arguments[0];
    for (unsigned m = 1; m < binCount; ++m) {
        arguments[m] = lOmega_m[m] * toLogValExp(-targetBeta * U_m[m]);
        normalization += arguments[m];
    }
    //normalize arguments, calculate weight corresponding to bin
    vector<double> weightFromBin(binCount);
    for (unsigned m = 0; m < binCount; ++m) {
        if (H_m[m] == 0) {
            weightFromBin[m] = 0;
        } else {
            arguments[m] /= normalization;
            weightFromBin[m] = toDouble(arguments[m]) / double(H_m[m]);
        }
    }

    DoubleSeriesCollection w_kn(numReplicas, 0);

    //calculate weight for each sample
    #pragma omp parallel for
    for (int k = 0; k < (signed)numReplicas; ++k) {
        unsigned N_k = m_kn[k]->size();
        w_kn[k] = new vector<double>(N_k);
        for (unsigned n = 0; n < N_k; ++n) {
            unsigned m = (*m_kn[k])[n];
            (*w_kn[k])[n] = weightFromBin[m];
        }
    }
    out << "." << flush;

    return w_kn;
}

void MultireweightHistosPT::saveLogDensityOfStates(const string& filename) {
    ofstream output(filename.c_str());
    output.precision(15);
    output.setf(std::ios::scientific, std::ios::floatfield);
    output  << "## logarithm of density of states, mrpt estimation\n"
            << "## ln(d.o.s.) normalized to zero at center bin\n"
            << "## energy (normalized by volume)\t ln(d.o.s.)"
            << endl;
    for (unsigned m = 0; m < binCount; ++m) {
        output  << U_m[m] / systemN << "\t" << lOmega_m[m] / lOmega_m[binCount / 2] << '\n';
    }
}

void MultireweightHistosPT::saveLogDensityOfStatesIsing(const string& filename) {
    ofstream output(filename.c_str());
    output.precision(15);
    output.setf(std::ios::scientific, std::ios::floatfield);
    output  << "## logarithm of density of states, mrpt estimation\n"
            << "## ln(d.o.s.) normalized to ln(2) for the lowest energy entry\n"
            << "## energy (normalized by volume)\t ln(d.o.s.)"
            << endl;
    LogVal norm = lOmega_m[0] / LogVal(2);
    //Normierung ist nur dann korrekt und stimmig mit Beale,
    //wenn U[0] == -2 * systemN
    //d.h., wenn Grundzustand wirklich erreicht!
    for (unsigned m = 0; m < binCount; ++m) {
        output  << U_m[m] / systemN << "\t" << lOmega_m[m] / norm << '\n';
    }
}

void MultireweightHistosPT::savePartitionFunctions(const string& filename) {
    out << "Saving estimated partition functions to file " << filename << "... " << flush;
    DataSeriesWriter<vector<LogVal> > saveLogs;
    saveLogs.setData(&lZ_l);
    saveLogs.addHeaderText("estimated partition functions ordered by beta-index, logarithmized");
    saveLogs.writeToFile(filename, 16);
    out << "Done." << endl;
}

void MultireweightHistosPT::loadPartitionFunctions(const string &filename) {
    out << "Loading ln(Z) estimates from file " << filename << endl;
    try {
        DataSeriesLoader<LogVal> loadLogs;
        loadLogs.readFromFile(filename);
        lZ_l = *(loadLogs.getData(0));  //copy
        loadLogs.deleteData();

        updateEffectiveCounts();
        updateDensityOfStates();            //else lZ_l would be overridden in findPartitionFunctionsAndDensityOfStates
    } catch (ReadError& err) {
        cerr << "Error loading file with ln(Z)." << endl;
    }
}

void MultireweightHistosPT::reweight1stMomentInternalWithoutErrors(const DoubleSeriesCollection& timeSeries, const DoubleSeriesCollection& w_kn, double& firstMoment) {
    double first = 0;
    #pragma omp parallel for reduction( + : first)
    for (int k = 0; k < (signed)numReplicas; ++k) {
        unsigned N_k = timeSeries[k]->size();
        for (unsigned n = 0; n < N_k; ++n) {
            double v = (*timeSeries[k])[n];
            first += v * (*w_kn[k])[n];
        }
        out << ".";
    }
    out << '\n';
    firstMoment = first;
}

void MultireweightHistosPT::reweight1stMoment2ndMomentInternalWithoutErrors(
        const DoubleSeriesCollection& timeSeries,
        const DoubleSeriesCollection& w_kn,
        double& firstMoment, double& secondMoment) {
    double first = 0;
    double second = 0;
    #pragma omp parallel for reduction( + : first, second)
    for (int k = 0; k < (signed)numReplicas; ++k) {
        unsigned N_k = timeSeries[k]->size();
        for (unsigned n = 0; n < N_k; ++n) {
            double v = (*timeSeries[k])[n];
            first += v * (*w_kn[k])[n];
            second += v*v * (*w_kn[k])[n];
        }
    }
    out << '#' << flush;
    firstMoment = first;
    secondMoment = second;
}

void MultireweightHistosPT::reweight2ndMoment4thMomentInternalWithoutErrors(
        const DoubleSeriesCollection& timeSeries,
        const DoubleSeriesCollection& w_kn,
        double& secondMoment, double& fourthMoment) {
    double second = 0;
    double fourth = 0;
    #pragma omp parallel for reduction( + : second, fourth)
    for (int k = 0; k < (signed)numReplicas; ++k) {
        unsigned N_k = timeSeries[k]->size();
        for (unsigned n = 0; n < N_k; ++n) {
            double v = (*timeSeries[k])[n];
            second += v*v * (*w_kn[k])[n];
            fourth += v*v*v*v * (*w_kn[k])[n];
        }
    }
    out << '#' << flush;
    fourthMoment = fourth;
    secondMoment = second;
}

ReweightingResult MultireweightHistosPT::reweightWithoutErrorsInternal
        (double targetBeta, const DoubleSeriesCollection& w_kn) {
    out << "Reweighting without errors to beta=" << targetBeta << endl;

    //everything normalized by system volume
    double meanEnergy = 0;
    double meanEnergySquared = 0;
    double meanObservable = 0;
    double meanObservableSquared = 0;
    double meanObservableToTheFourth = 0;

    #pragma omp parallel for reduction( + : meanEnergy, meanEnergySquared, meanObservable, meanObservableSquared, meanObservableToTheFourth)
    for (int k = 0; k < (signed)numReplicas; ++k) {
        unsigned N_k = energyTimeSeries[k]->size();
        for (unsigned n = 0; n < N_k; ++n) {
            double e = (*energyTimeSeries[k])[n];
            double o = (*observableTimeSeries[k])[n];
            meanEnergy += e * (*w_kn[k])[n];
            meanEnergySquared += e*e * (*w_kn[k])[n];
            meanObservable += o * (*w_kn[k])[n];
            meanObservableSquared += o*o * (*w_kn[k])[n];
            meanObservableToTheFourth += o*o*o*o * (*w_kn[k])[n];
        }
        out << ".";
    }

    double heatCapacity = systemN * targetBeta * targetBeta * (meanEnergySquared - pow(meanEnergy, 2));
    double suscObservable = systemN * (meanObservableSquared - pow(meanObservable, 2));
    double binderObservable = 1.0 - (meanObservableToTheFourth / (3 * pow(meanObservableSquared, 2)));

    out << " Done." << endl;

    return ReweightingResult(meanEnergy, heatCapacity, meanObservable, suscObservable, binderObservable);
}

double MultireweightHistosPT::reweightEnergy(double targetBeta) {
    DoubleSeriesCollection w_kn = computeWeights(targetBeta);

    double result = 0;
    reweight1stMomentInternalWithoutErrors(energyTimeSeries, w_kn, result);

    destroyAll(w_kn);

    return result;
}

double MultireweightHistosPT::reweightSpecificHeat(double targetBeta) {
    DoubleSeriesCollection w_kn = computeWeights(targetBeta);

    double firstMoment = 0;
    double secondMoment = 0;
    reweight1stMoment2ndMomentInternalWithoutErrors(energyTimeSeries, w_kn, firstMoment, secondMoment);
    double result = systemN * pow(targetBeta, 2) *
            (secondMoment - pow(firstMoment, 2));;

    destroyAll(w_kn);

    return result;
}

double MultireweightHistosPT::reweightObservable(double targetBeta) {
    DoubleSeriesCollection w_kn = computeWeights(targetBeta);

    double result = 0;
    reweight1stMomentInternalWithoutErrors(observableTimeSeries, w_kn, result);

    destroyAll(w_kn);

    return result;
}

double MultireweightHistosPT::reweightObservableSusceptibility(double targetBeta) {
    DoubleSeriesCollection w_kn = computeWeights(targetBeta);

    double firstMoment = 0;
    double secondMoment = 0;
    reweight1stMoment2ndMomentInternalWithoutErrors(observableTimeSeries, w_kn, firstMoment, secondMoment);
    double result = systemN * (secondMoment - pow(firstMoment, 2));

    destroyAll(w_kn);

    return result;
}

double MultireweightHistosPT::reweightObservableBinder(double targetBeta) {
    DoubleSeriesCollection w_kn = computeWeights(targetBeta);

    double secondMoment = 0;
    double fourthMoment = 0;
    reweight2ndMoment4thMomentInternalWithoutErrors(
            observableTimeSeries, w_kn, secondMoment, fourthMoment);
    double result = 1.0 - (fourthMoment / (3 * pow(secondMoment, 2)));

    destroyAll(w_kn);

    return result;
}

ReweightingResult MultireweightHistosPT::reweight(double targetBeta) {
    DoubleSeriesCollection w_kn = computeWeights(targetBeta);

    ReweightingResult results = reweightWithoutErrorsInternal(targetBeta, w_kn);

    destroyAll(w_kn);

    return results;
}

ReweightingResult MultireweightHistosPT::reweightWithHistograms(double targetBeta, unsigned obsBinCount) {
    DoubleSeriesCollection w_kn = computeWeights(targetBeta);

    ReweightingResult results = reweightWithoutErrorsInternal(targetBeta, w_kn);
    results.energyHistogram = reweightEnergyHistogram(targetBeta);
    results.obsHistogram = reweightObservableHistogramUsingWeights(
            targetBeta, obsBinCount, w_kn);

    destroyAll(w_kn);

    return results;
}

HistogramDouble* MultireweightHistosPT::reweightEnergyHistogramWithoutErrors(
        double targetBeta) {
    out << "Reweighting energy histogram... ";
    HistogramDouble* result = new HistogramDouble;
    result->minBin = minEnergyNormalized;
    result->maxBin = maxEnergyNormalized;
    result->spacing = binSize;
    result->binCount = binCount;
    result->N = systemN;
    result->beta = targetBeta;
    result->total = 0;

    vector<LogVal> arguments(binCount);
    arguments[0] = lOmega_m[0] * toLogValExp(-targetBeta * U_m[0]);
    LogVal normalization = arguments[0];
    for (unsigned m = 1; m < binCount; ++m) {
        arguments[m] = lOmega_m[m] * toLogValExp(-targetBeta * U_m[m]);
        normalization += arguments[m];
    }
    for (unsigned m = 0; m < binCount; ++m) {
        arguments[m] /= normalization;
        double prob = toDouble(arguments[m]);
        double energyNormalized = U_m[m] / systemN;
        result->histo[energyNormalized] = prob;
        result->total += prob;
    }
    result->headerLines += "## MRPT reweighted histogram of normalized energy\n";
    result->updateMeta();
    out << "Done." << endl;
    return result;
}

HistogramDouble* MultireweightHistosPT::reweightObservableHistogramUsingWeights(double targetBeta, unsigned obsBinCount,
            const DoubleSeriesCollection& w_kn) {
    HistogramDouble* result = new HistogramDouble;

    vector<double> obsHisto(obsBinCount, 0.0);
    reweightObservableHistogramInternal(targetBeta, obsBinCount, w_kn, obsHisto);

    result->assignVector(obsHisto, minObservableNormalized, maxObservableNormalized, targetBeta, systemN);
    result->headerLines += "## MRPT reweighted histogram of normalized " + observable + "\n";
    result->updateMeta();

    return result;
}

void MultireweightHistosPT::computeAndSaveHistogramCrossCorr() {
    out << "Estimating and saving histogram cross correlation tables... " << flush;
    typedef boost::multi_array<double, 3> Double3Array;
    Double3Array rho_kmm;      //cross correlation coefficients
    rho_kmm.resize(boost::extents[numReplicas][binCount][binCount]);
    //in this case:
    // rho_kmn = (-p_km * p_kn) / sqrt(p_km(1-p_km)*p_kn(1-p_kn))
    #pragma omp parallel for
    for (int k = 0; k < (signed)numReplicas; ++k) {
        unsigned N_k = energyTimeSeries[k]->size();
        if (N_k == 0) {
            //skip non-added time series
            continue;
        }
        for (unsigned m1 = 0; m1 < binCount; ++m1) {
            for (unsigned m2 = 0; m2 < binCount; ++m2) {
                if (m1 == m2) {
                    rho_kmm[k][m1][m2] = 1.0;
                } else {
                    double p_km = double(H_km[k][m1]) / N_k;
                    double p_kn = double(H_km[k][m2]) / N_k;
                    rho_kmm[k][m1][m2] =
                            -p_km*p_kn / sqrt(p_km*(1-p_km)*p_kn*(1-p_kn));
                }
            }
        }
    }

    dlib::create_directory("crosscorr");
    dlib::set_current_dir("crosscorr");
    for (unsigned k = 0; k < numReplicas; ++k) {
        unsigned N_k = energyTimeSeries[k]->size();
        if (N_k == 0) {
            //skip non-added time series
            continue;
        }
        ofstream outCrossCorr(
                ("crosscorr-k" + numToString(k) + ".table").c_str());
        outCrossCorr <<
             "## rho_kmn = (-p_km * p_kn) / sqrt(p_km(1-p_km)*p_kn(1-p_kn))\n";
        outCrossCorr << "# k = " << k << "\n";
        for (unsigned m1 = 0; m1 < binCount; ++m1) {
            for (unsigned m2 = 0; m2 < binCount; ++m2) {
                outCrossCorr << rho_kmm[k][m1][m2] << '\t';
            }
            outCrossCorr << '\n';
        }
    }
    dlib::set_current_dir("..");

    out << "done" << endl;
}

void MultireweightHistosPT::computeAndSaveHistogramCrossCorrAlt() {
    out << "Estimating and saving histogram cross correlation tables (alternative version) " << endl;
    typedef boost::multi_array<double, 3> Double3Array;
    Double3Array rho_kmm;      //cross correlation coefficients
    rho_kmm.resize(boost::extents[numReplicas][binCount][binCount]);
    std::fill(rho_kmm.data(), rho_kmm.data() + rho_kmm.num_elements(), 0.0);
    //in this case:
    // rho_kmn = <(\psi_m - <psi_m>)(\psi_n - <psi_n>)> /
    //            sqrt(Var(psi_m)Var(psi_n))>
    #pragma omp parallel for
    for (int k = 0; k < (signed)numReplicas; ++k) {
        unsigned N_k = energyTimeSeries[k]->size();
        if (N_k == 0) {
            //skip non-added time series
            continue;
        }
        //normalized histogram:
        vector<double> p_m(binCount);
        for (unsigned m = 0; m < binCount; ++m) {
            p_m[m] = double(H_km[k][m]) / N_k;
        }
        for (unsigned n = 0; n < N_k; ++n) {
            unsigned m = (*m_kn[k])[n];        //occupied bin
            //contribution for this and the other bins:
            for (unsigned m1 = 0; m1 < m; ++m1) {
                for (unsigned m2 = 0; m2 < m; ++m2) {
                    rho_kmm[k][m1][m2] += (.0 - p_m[m1]) * (.0 - p_m[m2]);
                }
                rho_kmm[k][m1][m] += (.0 - p_m[m1]) * (1.0 - p_m[m]);
                for (unsigned m2 = m + 1; m2 < binCount; ++m2) {
                    rho_kmm[k][m1][m2] += (.0 - p_m[m1]) * (.0 - p_m[m2]);
                }
            }
            rho_kmm[k][m][m] += (1.0 - p_m[m]) * (1.0 - p_m[m]);
            for (unsigned m1 = m + 1; m1 < binCount; ++m1) {
                for (unsigned m2 = 0; m2 < m; ++m2) {
                    rho_kmm[k][m1][m2] += (.0 - p_m[m1]) * (.0 - p_m[m2]);
                }
                rho_kmm[k][m1][m] += (.0 - p_m[m1]) * (1.0 - p_m[m]);
                for (unsigned m2 = m + 1; m2 < binCount; ++m2) {
                    rho_kmm[k][m1][m2] += (.0 - p_m[m1]) * (.0 - p_m[m2]);
                }
            }
        }
        //normalize by number of samples and bin variances
        for (unsigned m1 = 0; m1 < binCount; ++m1) {
            for (unsigned m2 = 0; m2 < binCount; ++m2) {
                rho_kmm[k][m1][m2] /=
                        N_k * sqrt(p_m[m1]*(1-p_m[m1])*p_m[m2]*(1-p_m[m2]));
            }
        }
        out << "." << flush;
    }

    dlib::create_directory("crosscorr-alt");
    dlib::set_current_dir("crosscorr-alt");
    for (unsigned k = 0; k < numReplicas; ++k) {
        unsigned N_k = energyTimeSeries[k]->size();
        if (N_k == 0) {
            //skip non-added time series
            continue;
        }
        ofstream outCrossCorr(
                ("crosscorr-k" + numToString(k) + ".table").c_str());
        outCrossCorr <<
                "## rho_kmn = <(psi_m - <psi_m>)(psi_n - <psi_n>)> / "
                "sqrt(Var(psi_m)Var(psi_n))>\n";
        outCrossCorr << "# k = " << k << "\n";
        for (unsigned m1 = 0; m1 < binCount; ++m1) {
            for (unsigned m2 = 0; m2 < binCount; ++m2) {
                outCrossCorr << rho_kmm[k][m1][m2] << '\t';
            }
            outCrossCorr << '\n';
        }
    }
    dlib::set_current_dir("..");

    out << " done" << endl;
}

void MultireweightHistosPT::reweightObservableHistogramInternal(double targetBeta, unsigned obsBinCount,
        const DoubleSeriesCollection& w_kn, std::vector<double>& obsHisto) {
    const double SMALL = 1e-10;     //to fit maxObservableNormalized into the highest bin
    double obsBinSize = (maxObservableNormalized - minObservableNormalized + SMALL) / obsBinCount;

    for (unsigned k = 0; k < numReplicas; ++k) {
        unsigned N_k = observableTimeSeries[k]->size();
        for (unsigned n = 0; n < N_k; ++n) {
            unsigned curBin = static_cast<int>(((*observableTimeSeries[k])[n] - minObservableNormalized) / obsBinSize);
            obsHisto[curBin] += (*w_kn[k])[n];
        }
    }
}

HistogramDouble* MultireweightHistosPT::reweightObservableHistogramWithoutErrors(
        double targetBeta, unsigned obsBinCount) {
    out << "Reweighting " << observable << " to generate histogram at beta=" << targetBeta << ", " << obsBinCount << " bins... " << flush;
    DoubleSeriesCollection w_kn = computeWeights(targetBeta);
    HistogramDouble* result = reweightObservableHistogramUsingWeights(targetBeta, obsBinCount, w_kn);
    destroyAll(w_kn);
    out << "Done." << endl;
    return result;
}






class SuscMinCallable {
    MultireweightHistosPT* master; map<double, double>& pointsEvaluated;
public:
    SuscMinCallable(MultireweightHistosPT* master, map<double, double>& pointsEvaluated) :
        master(master), pointsEvaluated(pointsEvaluated) { }
    double operator()(double beta) {
        double susc = master->reweightObservableSusceptibility(beta);
        pointsEvaluated[beta] = susc;
        master->out << " beta: " << beta << " => " << susc << endl;
        return -susc;
    }
};

class BinderMinCallable {
    MultireweightHistosPT* master; map<double, double>& pointsEvaluated;
public:
    BinderMinCallable(MultireweightHistosPT* master, map<double, double>& pointsEvaluated) :
        master(master), pointsEvaluated(pointsEvaluated) { }
    double operator()(double beta) {
        double binder = master->reweightObservableBinder(beta);
        pointsEvaluated[beta] = binder;
        master->out << " beta: " << beta << " => " << binder << endl;
        return binder;
    }
};

class SpecificHeatDiscreteMinCallable {
    MultireweightHistosPT* master;
    map<double, double>& pointsEvaluated;
public:
    SpecificHeatDiscreteMinCallable (MultireweightHistosPT* master, map<double, double>& pointsEvaluated) :
        master(master), pointsEvaluated(pointsEvaluated) { }
    double operator()(double beta) {
        //double specHeat = master->reweightSpecificHeat(beta);
        double specHeat = static_cast<MultireweightHistosPT*>(master)->reweightDiscrete(beta).heatCapacity;
        pointsEvaluated[beta] = specHeat;
        master->out << " beta: " << beta << " => " << specHeat << endl;
        return -specHeat;
    }
};

void MultireweightHistosPT::findMaxObservableSusceptibility(
        double& betaMax, double& suscMax,
        map<double, double>& pointsEvaluated,
        double betaStart, double betaEnd) {
    SuscMinCallable f(this, pointsEvaluated);
    out << "Searching max of " << observable << " susceptibility...\n"
        << " intermediate points:" << endl;
    brentMinimize(betaMax, suscMax, f, betaStart, betaEnd);
    suscMax = -suscMax;
    out << "final result: beta: " << betaMax << " => " << suscMax << endl;
}

void MultireweightHistosPT::findMaxSpecificHeatDiscrete(
        double& betaMax, double& specificHeatMax,
        map<double, double>& pointsEvaluated,
        double betaStart, double betaEnd) {
    SpecificHeatDiscreteMinCallable f(this, pointsEvaluated);
    out << "Searching max of specific heat (discrete!)...\n"
        << " intermediate points:" << endl;
    brentMinimize(betaMax, specificHeatMax, f, betaStart, betaEnd);
    specificHeatMax = -specificHeatMax;
    out << "final result: beta: " << betaMax << " => " << specificHeatMax << endl;
}

void MultireweightHistosPT::findMinBinder(double& betaMin, double& binderMin,
        std::map<double, double>& pointsEvaluated,
        double betaStart, double betaEnd) {
    BinderMinCallable f(this, pointsEvaluated);
    out << "Searching min of " << observable << " binder cumulant...\n"
        << " intermediate points:" << endl;
    brentMinimize(betaMin, binderMin, f, betaStart, betaEnd);
    out << "final result: beta: " << betaMin << " => " << binderMin << endl;
}


class EnergyHistogramPeakDiffMinCallable {
    MultireweightHistosPT* master;
    HistogramDouble* foundHistogram;
    double tolerance;
public:
    EnergyHistogramPeakDiffMinCallable(double tolerance,
            MultireweightHistosPT* master) :
        master(master), foundHistogram(0), tolerance(tolerance)
    {}
    double operator()(double beta) {
        destroy(foundHistogram);
        foundHistogram = master->reweightEnergyHistogramWithoutErrors(
                beta);
        double peakDiff = histogramPeakDiff(foundHistogram, tolerance);
        master->out
            << " beta: " << beta << " => peakDiff:" << peakDiff << endl;
        return peakDiff;
    }
    HistogramDouble* getHistogram() {
        return foundHistogram;
    }
};

class ObsHistogramPeakDiffMinCallable {
    MultireweightHistosPT* master;
    HistogramDouble* foundHistogram;
    double tolerance;
    unsigned numBins;
public:
    ObsHistogramPeakDiffMinCallable(double tolerance, unsigned numBins,
            MultireweightHistosPT* master) :
        master(master), foundHistogram(0),
        tolerance(tolerance), numBins(numBins)
    {}
    double operator()(double beta) {
        destroy(foundHistogram);
        foundHistogram = master->
                reweightObservableHistogramWithoutErrors(
                        beta, numBins);
        double peakDiff = histogramPeakDiff(foundHistogram, tolerance);
        master->out
            << " beta: " << beta << " => peakDiff:" << peakDiff << endl;
        return peakDiff;
    }
    HistogramDouble* getHistogram() {
        return foundHistogram;
    }
};

class EnergyHistogramWeightDiffMinCallable {
    MultireweightHistosPT* master;
    HistogramDouble* foundHistogram;
    double tolerance;
    double cutOff;
public:
    EnergyHistogramWeightDiffMinCallable(double tolerance,
            double cutOff, MultireweightHistosPT* master) :
        master(master), foundHistogram(0), tolerance(tolerance), cutOff(cutOff)
    {}
    double operator()(double beta) {
        destroy(foundHistogram);
        foundHistogram = master->reweightEnergyHistogramWithoutErrors(
                beta);
        double weightDiff = histogramWeightDiff(foundHistogram,
                cutOff);
        master->out
            << " beta: " << beta << " => weightDiff:" << weightDiff << endl;
        return weightDiff;
    }
    HistogramDouble* getHistogram() {
        return foundHistogram;
    }
};

class ObsHistogramWeightDiffMinCallable {
    MultireweightHistosPT* master;
    HistogramDouble* foundHistogram;
    double tolerance;
    double cutOff;
    unsigned numBins;
public:
    ObsHistogramWeightDiffMinCallable(double tolerance,
            double cutOff, unsigned numBins,
            MultireweightHistosPT* master) :
        master(master), foundHistogram(0),
        tolerance(tolerance), cutOff(cutOff), numBins(numBins)
    {}
    double operator()(double beta) {
        destroy(foundHistogram);
        foundHistogram = master->
                reweightObservableHistogramWithoutErrors(
                        beta, numBins);
        double weightDiff = histogramWeightDiff(foundHistogram,
                cutOff);
        master->out << " beta: " << beta
                << " => weightDiff:" << weightDiff << endl;
        return weightDiff;
    }
    HistogramDouble* getHistogram() {
        return foundHistogram;
    }
};


void MultireweightHistosPT::findEnergyEqualHeight(
        double& betaDouble, double& relDip, HistogramDouble*& histoResult,
        double betaStart, double betaEnd, double tolerance) {
    out << "Searching energy histogram with equal height double peak\n"
        << " intermediate points:" << endl;
    EnergyHistogramPeakDiffMinCallable f(tolerance, this);
    double peakDiff = -1.;
    brentMinimize(betaDouble, peakDiff, f, betaStart, betaEnd);
    histoResult = f.getHistogram();
    relDip = histogramRelativeDip(histoResult, tolerance);
    out << "final result: beta: " << betaDouble
        << " => peakDiff:" << peakDiff
        << ", relativeDip: " << relDip << endl;
}

void MultireweightHistosPT::findObsEqualHeight(
        double& betaDouble, double& relDip, HistogramDouble*& histoResult,
        double betaStart, double betaEnd,
        unsigned numBins, double tolerance) {
    out << "Searching " << observable
        << " histogram with " << numBins
        << " bins with equal height double peak\n"
        << " intermediate points:" << endl;
    ObsHistogramPeakDiffMinCallable f(tolerance, numBins, this);
    double peakDiff = -1.;
    brentMinimize(betaDouble, peakDiff, f, betaStart, betaEnd);
    histoResult = f.getHistogram();
    relDip = histogramRelativeDip(histoResult, tolerance);
    out << "final result: beta: " << betaDouble
        << " => peakDiff:" << peakDiff
        << ", relativeDip: " << relDip << endl;
}

void MultireweightHistosPT::findEnergyEqualWeight(
        double& betaDouble, double& relDip, HistogramDouble*& histoResult,
        const HistogramDouble* equalHeightHisto,
        double betaStart, double betaEnd, double tolerance) {
    out << "Searching energy histogram with equal weight double peak\n"
        << " intermediate points:" << endl;
    double cutOff = histogramMinimumLocation(equalHeightHisto, tolerance);
    EnergyHistogramWeightDiffMinCallable f(tolerance, cutOff, this);
    double weightDiff = -1.;
    brentMinimize(betaDouble, weightDiff, f, betaStart, betaEnd);
    histoResult = f.getHistogram();
    relDip = histogramRelativeDip(histoResult, tolerance);
    out << "final result: beta: " << betaDouble
        << " => weightDiff:" << weightDiff
        << ", relativeDip: " << relDip << endl;
}

void MultireweightHistosPT::findObsEqualWeight(
        double& betaDouble, double& relDip, HistogramDouble*& histoResult,
        const HistogramDouble* equalHeightHisto,
        double betaStart, double betaEnd, unsigned numBins, double tolerance) {
    out << "Searching " << observable
        << " histogram with " << numBins
        << " bins with equal weight double peak\n"
        << " intermediate points:" << endl;
    double cutOff = histogramMinimumLocation(equalHeightHisto, tolerance);
    ObsHistogramWeightDiffMinCallable f(tolerance, cutOff, numBins, this);
    double weightDiff = -1.;
    brentMinimize(betaDouble, weightDiff, f, betaStart, betaEnd);
    histoResult = f.getHistogram();
    relDip = histogramRelativeDip(histoResult, tolerance);
    out << "final result: beta: " << betaDouble
        << " => weightDiff:" << weightDiff
        << ", relativeDip: " << relDip << endl;
}

void MultireweightHistosPT::energyRelDip(double& relDip, HistogramDouble*& histoResult,
            double targetBeta, double tolerance) {
    histoResult = reweightEnergyHistogramWithoutErrors(targetBeta);
    relDip = histogramRelativeDip(histoResult, tolerance);
}

void MultireweightHistosPT::obsRelDip(double& relDip, HistogramDouble*& histoResult,
            double targetBeta, unsigned numBins, double tolerance) {
    histoResult = reweightObservableHistogramWithoutErrors(targetBeta,
            numBins);
    relDip = histogramRelativeDip(histoResult, tolerance);
}
