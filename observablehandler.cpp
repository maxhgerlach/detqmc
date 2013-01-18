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
	typedef DataMapWriter<std::string, num> StringNumWriter;

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




