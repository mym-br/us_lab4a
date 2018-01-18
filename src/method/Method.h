#ifndef METHOD_H_
#define METHOD_H_

#include <string>

#include <boost/unordered_map.hpp>



namespace Lab {

class Project;

// The static functions in this class can only be accessed by one thread at a time.
class Method {
public:
	enum Type {
		INVALID,
		SINGLE_ACQUISITION,
		SECTORIAL_SCAN_SP_NETWORK,
		SECTORIAL_SCAN_SP_NETWORK_CONTINUOUS,
		SECTORIAL_SCAN_SP_NETWORK_TRIGGER,
		SECTORIAL_SCAN_SP_SAVED,
		STA_SECTORIAL_SIMPLE_SIMULATED,
		STA_SECTORIAL_SIMPLE_SAVED,
		STA_SECTORIAL_SIMULATED,
		STA_SECTORIAL_DP_NETWORK,
		STA_SECTORIAL_DP_SAVED,
		STA_SECTORIAL_VECTORIAL_DP_SAVED,
		STA_SECTORIAL_SP_SAVED,
		STA_SECTORIAL_VECTORIAL_SP_SAVED,
		STA_SAVE_SIGNALS,
		SHOW_IMAGE,
		TEST
	};

	Method() {}
	virtual ~Method() {}

	virtual void execute() = 0;

	// This function must be called in main().
	static void fillNameMap();

	static Type findByName(const std::string& name);
	static std::string findByType(Type type);
	static Method* get(Project& project);
private:
	typedef boost::unordered_map<std::string, Type> NameMap;

	static NameMap nameMap_;
};

} // namespace Lab

#endif /* METHOD_H_ */
