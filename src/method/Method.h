#ifndef METHOD_H_
#define METHOD_H_

#include <string>

#include <boost/unordered_map.hpp>



namespace Lab {

enum class MethodType {
	invalid,
	single_acquisition,
	sectorial_scan_sp_network,
	sectorial_scan_sp_network_continuous,
	sectorial_scan_sp_network_trigger,
	sectorial_scan_sp_saved,
	sta_sectorial_simple_simulated,
	sta_sectorial_simple_saved,
	sta_sectorial_simulated,
	sta_sectorial_dp_network,
	sta_sectorial_dp_saved,
	sta_sectorial_vectorial_dp_saved,
	sta_sectorial_sp_saved,
	sta_sectorial_vectorial_sp_saved,
	sta_save_signals,
	show_image,
	test
};

class Project;

class MethodNameMap {
public:
	MethodNameMap();
	~MethodNameMap();

	MethodType findByName(const std::string& name);
private:
	typedef boost::unordered_map<std::string, MethodType> Map;

	Map map_;
};

class Method {
public:
	Method() {}
	virtual ~Method() {}

	virtual void execute() = 0;

	static MethodType findByName(const std::string& name) { return nameMap_.findByName(name); }
	static Method* get(Project& project);
private:
	static MethodNameMap nameMap_;
};

} // namespace Lab

#endif /* METHOD_H_ */
