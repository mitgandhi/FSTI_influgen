#ifndef _FSTINETLIC_
#define _FSTINETLIC_

#include <string>

namespace fsti_netlic
{
	enum license_type
	{
		NODE,
		NETWORK,
		NODE_THEN_NETWORK
	};

	void * checkout_license(license_type type, std::string interface_name);
	void release_license(void * lic);
};

#endif