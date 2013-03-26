/*
 * detqmc.h
 *
 * Handling of determinantal quantum Monte Carlo simulations
 *
 *  Created on: Dec 10, 2012
 *      Author: gerlach
 */

#ifndef DETQMC_H_
#define DETQMC_H_

#include <vector>
#include <functional>
#include <memory>
#include "metadata.h"
#include "parameters.h"
#include "detmodel.h"
#include "observablehandler.h"
#include "rngwrapper.h"



// Class handling the simulation
class DetQMC {
public:
	DetQMC(const ModelParams& parsmodel, const MCParams& parsmc);

	void run();			//carry out simulation determined by parsmc given in construction

	// do numSweeps thermalization sweeps
	void thermalize(unsigned numSweeps);
	// do numSweeps sweeps taking measurements
	void measure(unsigned numSweeps, unsigned measureInterval);
	// update results stored on disk
	void saveResults();

	virtual ~DetQMC();
protected:
	ModelParams parsmodel;
	MCParams parsmc;

	enum class GreenUpdateType {Simple, Stabilized};
	GreenUpdateType greenUpdateType;
	std::function<void()> sweepFunc;	//the replica member function that will be called
										//to perform a sweep (depending on greenUpdate).
										//adds a function pointer layer
	std::function<void()> sweepThermalizationFunc;		//during thermalization this may
														//be a different one

	MetadataMap modelMeta;
	MetadataMap mcMeta;
	RngWrapper rng;
	std::unique_ptr<DetModel> replica;
	typedef std::unique_ptr<ScalarObservableHandler> ObsPtr;
	typedef std::unique_ptr<VectorObservableHandler> VecObsPtr;
	std::vector<ObsPtr> obsHandlers;
	std::vector<VecObsPtr> vecObsHandlers;
	unsigned sweepsDone;						//Measurement sweeps done

	MetadataMap prepareMCMetadataMap() const;
};





#endif /* DETQMC_H_ */
