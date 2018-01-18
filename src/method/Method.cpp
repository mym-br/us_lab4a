#include "Method.h"

#include "Exception.h"
#include "Project.h"
#include "SectorialScanMethod.h"
#include "ShowImageMethod.h"
#include "SingleAcquisitionMethod.h"
#include "STAMethod.h"
#include "TestMethod.h"



namespace Lab {

Method::NameMap Method::nameMap_;



void
Method::fillNameMap()
{
	nameMap_["single_acquisition"                                       ] = SINGLE_ACQUISITION;
	nameMap_["sectorial_scan:single_precision:network"                  ] = SECTORIAL_SCAN_SP_NETWORK;
	nameMap_["sectorial_scan:single_precision:network:continuous"       ] = SECTORIAL_SCAN_SP_NETWORK_CONTINUOUS;
	nameMap_["sectorial_scan:single_precision:network:trigger"          ] = SECTORIAL_SCAN_SP_NETWORK_TRIGGER;
	nameMap_["sectorial_scan:single_precision:saved"                    ] = SECTORIAL_SCAN_SP_SAVED;
	nameMap_["sta:sectorial:simple:simulated"                           ] = STA_SECTORIAL_SIMPLE_SIMULATED;
	nameMap_["sta:sectorial:simple:saved"                               ] = STA_SECTORIAL_SIMPLE_SAVED;
	nameMap_["sta:sectorial:simulated"                                  ] = STA_SECTORIAL_SIMULATED;
	nameMap_["sta:sectorial:double_precision:network"                   ] = STA_SECTORIAL_DP_NETWORK;
	nameMap_["sta:sectorial:double_precision:saved"                     ] = STA_SECTORIAL_DP_SAVED;
	nameMap_["sta:sectorial:vectorial:double_precision:saved"           ] = STA_SECTORIAL_VECTORIAL_DP_SAVED;
	nameMap_["sta:sectorial:single_precision:saved"                     ] = STA_SECTORIAL_SP_SAVED;
	nameMap_["sta:sectorial:vectorial:single_precision:saved"           ] = STA_SECTORIAL_VECTORIAL_SP_SAVED;
	nameMap_["sta:save_signals"                                         ] = STA_SAVE_SIGNALS;
	nameMap_["show_image"                                               ] = SHOW_IMAGE;
	nameMap_["test"                                                     ] = TEST;
}

Method*
Method::get(Project& project)
{
	switch (project.method()) {
	case SINGLE_ACQUISITION:
		return new SingleAcquisitionMethod(project);
	case SECTORIAL_SCAN_SP_NETWORK:                           // falls through
	case SECTORIAL_SCAN_SP_NETWORK_CONTINUOUS:                // falls through
	case SECTORIAL_SCAN_SP_NETWORK_TRIGGER:                   // falls through
	case SECTORIAL_SCAN_SP_SAVED:
		return new SectorialScanMethod<float>(project);
	case STA_SECTORIAL_SIMPLE_SIMULATED:                      // falls through
	case STA_SECTORIAL_SIMPLE_SAVED:                          // falls through
	case STA_SECTORIAL_SIMULATED:                             // falls through
	case STA_SECTORIAL_DP_NETWORK:                            // falls through
	case STA_SECTORIAL_DP_SAVED:                              // falls through
	case STA_SECTORIAL_VECTORIAL_DP_SAVED:                    // falls through
	case STA_SAVE_SIGNALS:
		return new STAMethod<double>(project);
	case STA_SECTORIAL_SP_SAVED:
	case STA_SECTORIAL_VECTORIAL_SP_SAVED:
		return new STAMethod<float>(project);
	case SHOW_IMAGE:
		return new ShowImageMethod(project);
	case TEST:
		return new TestMethod(project);
	default:
		THROW_EXCEPTION(InvalidParameterException, "Invalid method: " << project.method() << '.');
	}
}

Method::Type
Method::findByName(const std::string& name)
{
	NameMap::const_iterator iter = nameMap_.find(name);
	if (iter == nameMap_.end()) {
		THROW_EXCEPTION(InvalidParameterException, "Could not find a method with name \"" << name << "\".");
	}
	return iter->second;
}

// This search is slow.
std::string
Method::findByType(Method::Type type)
{
	for (NameMap::const_iterator iter = nameMap_.begin(); iter != nameMap_.end(); ++iter) {
		if (iter->second == type) {
			return iter->first;
		}
	}
	THROW_EXCEPTION(InvalidParameterException, "Could not find a method with type " << static_cast<int>(type) << '.');
}

} // namespace Lab
