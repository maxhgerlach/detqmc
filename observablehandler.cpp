/*
 * observablehandler.cpp
 *
 *  Created on: Dec 13, 2012
 *      Author: gerlach
 */

#include "observablehandler.h"

void outputResults(const std::vector<std::unique_ptr<ScalarObservableHandler>>& obsHandlers) {
	typedef std::map<std::string, num> StringNumMap;
	typedef std::shared_ptr<StringNumMap> StringNumMapPtr;

	StringNumMapPtr values(new StringNumMap);
	StringNumMapPtr errors(new StringNumMap);
	for (auto p = obsHandlers.cbegin(); p != obsHandlers.cend(); ++p) {
		num val, err;
		std::tie(val, err) = (*p)->evaluateJackknife();
		std::string obsname = (*p)->name;
		(*values)[obsname] = val;
		(*errors)[obsname] = err;
	}
	DataMapWriter<std::string, num> output;
	output.setData(values);
	output.setErrors(errors);
	output.addHeaderText("Monte Carlo results for observable expectation values");
	output.addMetadataMap((*obsHandlers.begin())->metaModel);
	output.addMetadataMap((*obsHandlers.begin())->metaMC);
	output.addMeta("key", "observable");
	output.addHeaderText("observable\t value \t error");
	output.writeToFile("results.values");
}

void outputResults(const std::vector<std::unique_ptr<VectorObservableHandler>>& obsHandlers) {
	typedef std::map<num, num> NumMap;
	typedef std::shared_ptr<NumMap> NumMapPtr;
	typedef DataMapWriter<num,num> NumMapWriter;

	for (auto p = obsHandlers.cbegin(); p != obsHandlers.cend(); ++p) {
		const std::unique_ptr<VectorObservableHandler>& obsptr = *p;
		unsigned N = obsptr->getVectorSize();
		arma::Col<num> values, errors;
		std::tie(values, errors) = obsptr->evaluateJackknife();
		NumMapPtr valmap(new NumMap);
		NumMapPtr errmap(new NumMap);
		for (unsigned site = 0; site < N; ++site) {
			valmap->insert(std::make_pair(num(site), values[site]));
			errmap->insert(std::make_pair(num(site), errors[site]));
		}
		NumMapWriter output;
		output.setData(valmap);
		output.setErrors(errmap);
		output.addHeaderText("Monte Carlo results for vector observable " + obsptr->name +
				" expectation values");
		output.addMetadataMap(obsptr->metaModel);
		output.addMetadataMap(obsptr->metaMC);
		output.addMeta("key", "site");
		output.addMeta("observable", obsptr->name);
		output.addHeaderText("site \t value \t error");
		output.writeToFile("results-" + obsptr->name + ".values");
	}
}





