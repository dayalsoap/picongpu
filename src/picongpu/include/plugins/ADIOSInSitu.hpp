#pragma once

#include "plugins/ILightweightPlugin.hpp"

namespace picongpu
{
    namespace adiosinsitu
    {
	
	using namespace PMacc;

	class ADIOSInSitu : public ILightweightPlugin
	{
	public:
	    ADIOSInSitu()
		{
		    /* register our plugin during creation */
		    Environment<>::get().PluginConnector().registerPlugin(this);
		}

	    std::string pluginGetName() const
		{
		    return "ADIOSInSitu";
		}

	    void notify(uint32_t currentStep)
		{
		    /* notification callback for simulation step currentStep
		     * called every notifyPeriod steps */
		    fprintf(stderr, "Jai's plugin notified.\n");
		}

	    void pluginRegisterHelp(po::options_description& desc)
		{
		    /* register command line parameters for your plugin */
		    desc.add_options()
			("adiosinsitu.period", po::value<uint32_t > (&notifyPeriod)->default_value(0),
			 "Enable ADIOSInSitu [for each n-th step]");
		}

	    void setMappingDescription(MappingDesc *cellDescription)
		{
		}

	private:
	    uint32_t notifyPeriod;

	    void pluginLoad()
		{
		    /* called when plugin is loaded, command line flags are available here
		     * set notification period for our plugin at the PluginConnector */
		    Environment<>::get().PluginConnector().setNotificationPeriod(this, notifyPeriod);
		}

	    void pluginUnload()
		{
		    /* called when plugin is unloaded, cleanup here */
		}
	};
    } // adiosinsitu
} // picongpu
