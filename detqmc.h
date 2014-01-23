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

#include "boost/preprocessor/comma.hpp"
#include "boost/timer/timer.hpp"
#include "boost/serialization/split_member.hpp"


class SerializeContentsKey;

// Class handling the simulation
class DetQMC {
public:
    //constructor to init a new simulation:
    DetQMC(const ModelParams& parsmodel, const MCParams& parsmc);

    //constructor to resume a simulation from a dumped state file:
    //we allow to change some MC parameters at this point:
    //  sweeps & saveInterval
    //if values > than the old values are specified, change them
    DetQMC(const std::string& stateFileName, const MCParams& newParsmc);


    //carry out simulation determined by parsmc given in construction,
    //- handle thermalization & measurement stages as necessary
    //- save state and results periodically
    //- if granted walltime is almost over, save state & results
    //  and exit gracefully
    void run();

    // update results stored on disk
    void saveResults();
    // dump simulation parameters and the current state to a Boost::S11n archive,
    // also write out information about the current simulation state to info.dat
    void saveState();

    virtual ~DetQMC();
protected:
    //helper for constructors -- set all parameters and initialize contained objects
    void initFromParameters(const ModelParams& parsmodel, const MCParams& parsmc);

    ModelParams parsmodel;
    MCParams parsmc;

    enum GreenUpdateType {GreenUpdateTypeSimple, GreenUpdateTypeStabilized};
    GreenUpdateType greenUpdateType;
//  std::function<void()> sweepFunc;    //the replica member function that will be called
//                                      //to perform a sweep (depending on greenUpdate).
//                                      //adds a function pointer layer
//  std::function<void()> sweepThermalizationFunc;      //during thermalization this may
//                                                      //be a different one

    MetadataMap modelMeta;
    MetadataMap mcMeta;
    RngWrapper rng;
    std::unique_ptr<DetModel> replica;
    typedef std::unique_ptr<ScalarObservableHandler> ObsPtr;
    typedef std::unique_ptr<VectorObservableHandler> VecObsPtr;
    std::vector<ObsPtr> obsHandlers;
    std::vector<VecObsPtr> vecObsHandlers;      //need to be pointers: holds both KeyValueObservableHandlers and VectorObservableHandlers
    uint32_t sweepsDone;                        //Measurement sweeps done
    uint32_t sweepsDoneThermalization;          //thermalization sweeps done

    uint32_t swCounter;         //helper counter in run() -- e.g. sweeps between measurements -- should also be serialized

    boost::timer::cpu_timer elapsedTimer;           //during this simulation run
    uint32_t curWalltimeSecs() {
        return elapsedTimer.elapsed().wall / 1000 / 1000 / 1000; // ns->mus->ms->s
    }
    uint32_t totalWalltimeSecs;             //this is serialized and carries the elapsed walltime in seconds
                                                    //accumulated over all runs, updated on call of saveResults()
    uint32_t walltimeSecsLastSaveResults;       //timer seconds at previous saveResults() call --> used to update totalWalltimeSecs
    uint32_t grantedWalltimeSecs;               //walltime the simulation is allowed to run


    MetadataMap prepareMCMetadataMap() const;

private:
    //Serialize only the content data that has changed after construction.
    //Only call for deserialization after DetQMC has already been constructed and initialized!

    //separate functions loadContents, saveContents; both employ serializeContentsCommon
    template<class Archive>
    void loadContents(Archive& ar) {
        serializeContentsCommon(ar);

        //The template member functions saveContents(Archive&) cannot be virtual,
        //so we have to resort to RTTI to serialize the right object.
        //Unfortunately this is fugly.
        if (DetHubbard<true, true>* p = dynamic_cast<DetHubbard<true, true>*>(replica.get())) {
            p->loadContents(SerializeContentsKey(), ar);
        } else
        if (DetHubbard<true, false>* p = dynamic_cast<DetHubbard<true, false>*>(replica.get())) {
            p->loadContents(SerializeContentsKey(), ar);
        } else
        if (DetHubbard<false, true>* p = dynamic_cast<DetHubbard<false, true>*>(replica.get())) {
            p->loadContents(SerializeContentsKey(), ar);
        } else
        if (DetHubbard<false, false>* p = dynamic_cast<DetHubbard<false, false>*>(replica.get())) {
            p->loadContents(SerializeContentsKey(), ar);
        } else
        if (DetSDW<true, CB_NONE>* p = dynamic_cast<DetSDW<true, CB_NONE>*>(replica.get())) {
            p->loadContents(SerializeContentsKey(), ar);
        } else
        if (DetSDW<true, CB_ASSAAD>* p = dynamic_cast<DetSDW<true, CB_ASSAAD>*>(replica.get())) {
            p->loadContents(SerializeContentsKey(), ar);
        } else
        if (DetSDW<true, CB_ASSAAD_BERG>* p = dynamic_cast<DetSDW<true, CB_ASSAAD_BERG>*>(replica.get())) {
            p->loadContents(SerializeContentsKey(), ar);
        } else
        if (DetSDW<false, CB_NONE>* p = dynamic_cast<DetSDW<false, CB_NONE>*>(replica.get())) {
            p->loadContents(SerializeContentsKey(), ar);
        } else
        if (DetSDW<false, CB_ASSAAD>* p = dynamic_cast<DetSDW<false, CB_ASSAAD>*>(replica.get())) {
            p->loadContents(SerializeContentsKey(), ar);
        } else
        if (DetSDW<false, CB_ASSAAD_BERG>* p = dynamic_cast<DetSDW<false, CB_ASSAAD_BERG>*>(replica.get())) {
            p->loadContents(SerializeContentsKey(), ar);
        }
        else {
            throw SerializationError("Tried to load contents of unsupported replica");
        }
    }

    template<class Archive>
    void saveContents(Archive& ar) {
        serializeContentsCommon(ar);

        //The template member functions saveContents(Archive&) cannot be virtual,
        //so we have to resort to RTTI to serialize the right object.
        //Unfortunately this is fugly.
        if (DetHubbard<true, true>* p = dynamic_cast<DetHubbard<true, true>*>(replica.get())) {
            p->saveContents(SerializeContentsKey(), ar);
        } else
        if (DetHubbard<true, false>* p = dynamic_cast<DetHubbard<true, false>*>(replica.get())) {
            p->saveContents(SerializeContentsKey(), ar);
        } else
        if (DetHubbard<false, true>* p = dynamic_cast<DetHubbard<false, true>*>(replica.get())) {
            p->saveContents(SerializeContentsKey(), ar);
        } else
        if (DetHubbard<false, false>* p = dynamic_cast<DetHubbard<false, false>*>(replica.get())) {
            p->saveContents(SerializeContentsKey(), ar);
        } else
        if (DetSDW<true, CB_NONE>* p = dynamic_cast<DetSDW<true, CB_NONE>*>(replica.get())) {
            p->saveContents(SerializeContentsKey(), ar);
        } else
        if (DetSDW<true, CB_ASSAAD>* p = dynamic_cast<DetSDW<true, CB_ASSAAD>*>(replica.get())) {
            p->saveContents(SerializeContentsKey(), ar);
        } else
        if (DetSDW<true, CB_ASSAAD_BERG>* p = dynamic_cast<DetSDW<true, CB_ASSAAD_BERG>*>(replica.get())) {
            p->saveContents(SerializeContentsKey(), ar);
        } else
        if (DetSDW<false, CB_NONE>* p = dynamic_cast<DetSDW<false, CB_NONE>*>(replica.get())) {
            p->saveContents(SerializeContentsKey(), ar);
        } else
        if (DetSDW<false, CB_ASSAAD>* p = dynamic_cast<DetSDW<false, CB_ASSAAD>*>(replica.get())) {
            p->saveContents(SerializeContentsKey(), ar);
        } else
        if (DetSDW<false, CB_ASSAAD_BERG>* p = dynamic_cast<DetSDW<false, CB_ASSAAD_BERG>*>(replica.get())) {
            p->saveContents(SerializeContentsKey(), ar);
        }
        else {
            throw SerializationError("Tried to save contents of unsupported replica");
        }
    }


    template<class Archive>
    void serializeContentsCommon(Archive& ar) {
        ar & rng;                   //serialize completely

        for (auto p = obsHandlers.begin(); p != obsHandlers.end(); ++p) {
            //ATM no further derived classes of ScalarObservableHandler have a method serializeContents
            (*p)->serializeContents(SerializeContentsKey(), ar);
        }
        for (auto p = vecObsHandlers.begin(); p != vecObsHandlers.end(); ++p) {
            //ATM no further derived classes of VectorObservableHandler have a method serializeContents
            (*p)->serializeContents(SerializeContentsKey(), ar);
        }
        ar & sweepsDone & sweepsDoneThermalization;
        ar & swCounter;

        ar & totalWalltimeSecs;
    }
};


//Only one member function of DetQMC is allowed to make instances of this class.
//In this way access to the member function serializeContents() of other classes
//is restricted.
// compare to http://stackoverflow.com/questions/6310720/declare-a-member-function-of-a-forward-declared-class-as-friend
class SerializeContentsKey {
  SerializeContentsKey() {} // default ctor private
  SerializeContentsKey(const SerializeContentsKey&) {} // copy ctor private

  // grant access to few methods
  template<class Archive>
  friend void DetQMC::saveContents(Archive& ar);
  template<class Archive>
  friend void DetQMC::loadContents(Archive& ar);
  template<class Archive>
  friend void DetQMC::serializeContentsCommon(Archive& ar);
};




#endif /* DETQMC_H_ */
