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

#include <boost/serialization/split_member.hpp>
#include <boost/serialization/map.hpp>				//for MetadataMap
#include "boost_serialize_uniqueptr.h"
#include "boost_serialize_vector_uniqueptr.h"

// Class handling the simulation
class DetQMC {
public:
	//constructor to init a new simulation:
	DetQMC(const ModelParams& parsmodel, const MCParams& parsmc);

	//constructor to resume a simulation from a dumped state file:
	//if newSweeps > 0 is specified, adjust the number of target sweeps
	DetQMC(const std::string& stateFileName, unsigned newSweeps = 0);


	void run();			//carry out simulation determined by parsmc given in construction

	// do numSweeps thermalization sweeps
	void thermalize(unsigned numSweeps);
	// do numSweeps sweeps taking measurements
	void measure(unsigned numSweeps, unsigned measureInterval);
	// update results stored on disk
	void saveResults();

	virtual ~DetQMC();
protected:
	//helper for constructors:
	void initFromParameters(const ModelParams& parsmodel, const MCParams& parsmc);

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
	std::vector<VecObsPtr> vecObsHandlers;		//need to be pointers: holds both KeyValueObservableHandlers and VectorObservableHandlers
	unsigned sweepsDone;						//Measurement sweeps done

	MetadataMap prepareMCMetadataMap() const;

protected:
	//TODO: use boost nvp (name-value-pair) serialization
	friend class boost::serialization::access;

	//serialization code that is equal for save/load
	template<class Archive>
	void serialize(Archive& ar, const unsigned int version) {
		//callbacks etc. need not be serialized as they are
		//determined via parsmodel & parsmc upon construction.
    	ar & parsmodel & parsmc;
//    	ar & modelMeta & mcMeta;
    	ar & rng;
//    	ar & greenUpdateType;
    	ar & replica;
    	ar & obsHandlers & vecObsHandlers;
    	ar & sweepsDone;
	}
};





#endif /* DETQMC_H_ */
