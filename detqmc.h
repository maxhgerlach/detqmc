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
#include "dethubbard.h"
#include "detsdw.h"
#include "observablehandler.h"
#include "rngwrapper.h"
#include "exceptions.h"

#include "boost/serialization/split_member.hpp"
//#include "boost/serialization/map.hpp"				//for MetadataMap
//#include "boost_serialize_uniqueptr.h"
//#include "boost_serialize_vector_uniqueptr.h"

class SerializeContentsKey;

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
	// dump simulation parameters and the current state to a Boost::S11n archive
	void saveState();

	virtual ~DetQMC();
protected:
	//helper for constructors -- set all parameters and initialize contained objects
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
	unsigned sweepsDoneThermalization;			//thermalization sweeps done

	MetadataMap prepareMCMetadataMap() const;

private:
	//Serialize only the content data that has changed after construction.
	//Only call for "deserialization" after DetQMC has already been constructed and initialized!
	template<class Archive>
	void serializeContents(Archive& ar) {
		//callbacks etc. need not be serialized as they are
		//determined via parsmodel & parsmc upon construction.
//    	ar & parsmodel & parsmc;
//    	ar & modelMeta & mcMeta;
    	ar & rng;					//serialize completely
//    	ar & greenUpdateType;

    	//The template member functions serializeContents(Archive&) cannot be virtual,
    	//so we have to resort to RTTI to serialize the right object.
    	if (DetHubbard* p = dynamic_cast<DetHubbard*>(replica.get())) {
    		p->serializeContents(SerializeContentsKey(), ar);
    	} else if (DetSDW* p = dynamic_cast<DetSDW*>(replica.get())) {
    		p->serializeContents(SerializeContentsKey(), ar);
    	} else {
    		throw SerializationError("Tried to serialize contents of unsupported replica");
    	}

    	for (auto p = obsHandlers.begin(); p != obsHandlers.end(); ++p) {
    		//ATM no further derived classes of ScalarObservableHandler have a method serializeContents
    		(*p)->serializeContents(SerializeContentsKey(), ar);
    	}
    	for (auto p = vecObsHandlers.begin(); p != vecObsHandlers.end(); ++p) {
    		//ATM no further derived classes of VectorObservableHandler have a method serializeContents
    		(*p)->serializeContents(SerializeContentsKey(), ar);
    	}
    	ar & sweepsDone & sweepsDoneThermalization;
	}
};


//BOOST_CLASS_EXPORT_GUID(DetQMC, "DetQMC")


//Only one member function of DetQMC is allowed to make instances of this class.
//In this way access to the member function serializeContents() of other classes
//is restricted.
// compare to http://stackoverflow.com/questions/6310720/declare-a-member-function-of-a-forward-declared-class-as-friend
class SerializeContentsKey {
  SerializeContentsKey() {} // default ctor private
  SerializeContentsKey(const SerializeContentsKey&) {} // copy ctor private

  // grant access to one method
  template<class Archive>
  friend void DetQMC::serializeContents(Archive& ar);
};




#endif /* DETQMC_H_ */
