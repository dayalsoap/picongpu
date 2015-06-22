#pragma once

#include "plugins/ILightweightPlugin.hpp"

#include <pthread.h>
#include <cassert>
#include <sstream>
#include <list>
#include <vector>

#include "types.h"
#include "simulation_types.hpp"
#include "plugins/adios/ADIOSInSitu.def"

#include "particles/frame_types.hpp"

#include <adios.h>
#include <adios_read.h>
#include <adios_error.h>

#include "fields/FieldB.hpp"
#include "fields/FieldE.hpp"
#include "fields/FieldJ.hpp"
#include "fields/FieldTmp.hpp"
#include "particles/operations/CountParticles.hpp"

#include "dataManagement/DataConnector.hpp"
#include "mappings/simulation/GridController.hpp"
#include "mappings/simulation/SubGrid.hpp"
#include "dimensions/GridLayout.hpp"
#include "pluginSystem/PluginConnector.hpp"
#include "simulationControl/MovingWindow.hpp"
#include "math/Vector.hpp"

#include "plugins/ILightweightPlugin.hpp"
#include <boost/mpl/vector.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/begin_end.hpp>
#include <boost/mpl/find.hpp>
#include <boost/filesystem.hpp>

#include <boost/type_traits.hpp>
#if !defined(_WIN32)
#include <unistd.h>
#endif

#define INITIAL_STEP 0
#define FINAL_STEP 10


namespace picongpu
{
    namespace adiosinsitu
    {
	
	using namespace PMacc;

	/**
	 * Helper function to define an adios variable.
	 **/
	template <unsigned DIM>
	int64_t defineAdiosVar(int64_t group_id,
			       const char * name,
			       const char * path,
			       enum ADIOS_DATATYPES type,
			       PMacc::math::UInt64<DIM> dimensions,
			       PMacc::math::UInt64<DIM> globalDimensions,
			       PMacc::math::UInt64<DIM> offset)
	{
	    int64_t var_id = 0;
	    if ((DIM == 1) && (globalDimensions.productOfComponents() == 1)) {
		/* scalars need empty size strings */
		var_id = adios_define_var(
		    group_id, name, path, type, 0, 0, 0);
	    } else {
		// first we have to define the variables that hold the array dimensions.		
		var_id = adios_define_var(
		    group_id, name, path, type,
		    dimensions.revert().toString(",", "").c_str(),
		    globalDimensions.revert().toString(",", "").c_str(),
		    offset.revert().toString(",", "").c_str());
	    }
	    
	    log<picLog::INPUT_OUTPUT > ("ADIOS: Defined varID=%1% for '%2%' at %3% for %4%/%5% elements") %
		var_id % std::string(name) % offset.toString() % dimensions.toString() % globalDimensions.toString();
	    return var_id;
	}

// struct ThreadParams
// {
//     uint32_t currentStep;                   /** current simulation step */
//     std::string fullFilename;

//     /** current dump is a checkpoint */
//     bool isCheckpoint;
//     ADIOS_FILE* fp;                          /* file pointer for checkpoint file */

//     MPI_Comm adiosComm;                     /* MPI communicator for adios lib */
//     bool adiosBufferInitialized;            /* set if ADIOS buffer has been allocated */
//     int64_t adiosFileHandle;                /* ADIOS file handle */
//     int64_t adiosGroupHandle;               /* ADIOS group handle */
//     uint64_t adiosGroupSize;                /* size of ADIOS group in bytes */
//     uint32_t adiosAggregators;              /* number of ADIOS aggregators for MPI_AGGREGATE */
//     uint32_t adiosOST;                      /* number of ADIOS OST for MPI_AGGREGATE */
//     std::string adiosBasePath;              /* base path for the current step */
//     std::string adiosCompression;           /* ADIOS data transform compression method */

//     PMacc::math::UInt64<simDim> fieldsSizeDims;
//     PMacc::math::UInt64<simDim> fieldsGlobalSizeDims;
//     PMacc::math::UInt64<simDim> fieldsOffsetDims;

//     std::list<int64_t> adiosFieldVarIds;        /* var IDs for fields in order of appearance */
//     std::list<int64_t> adiosParticleAttrVarIds; /* var IDs for particle attributes in order of appearance */
//     std::list<int64_t> adiosSpeciesIndexVarIds; /* var IDs for species index tables in order of appearance */

//     GridLayout<simDim> gridLayout;
//     MappingDesc *cellDescription;

//     float *fieldBfr;                                /* temp. buffer for fields */

//     Window window;                                  /* window describing the volume to be dumped */

//     DataSpace<simDim> localWindowToDomainOffset;    /** offset from local moving window to local domain */
// };
	
	class ADIOSInSitu : public ILightweightPlugin
	{
	public:
	    ADIOSInSitu()
		{
		    /* register our plugin during creation */
		    Environment<>::get().PluginConnector().registerPlugin(this);
		    fprintf(stderr, "Plugin instantiated.\n");
		}

	    std::string pluginGetName() const
		{
		    return "ADIOSInSitu";
		}

	    void notify(uint32_t currentStep)
		{
		    const PMacc::Selection<simDim>& localDomain =
			Environment<simDim>::get().SubGrid().getLocalDomain();
		    mThreadParams.currentStep = currentStep;
		    mThreadParams.cellDescription = this->cellDescription;
		    double start = MPI_Wtime();
		    __getTransactionEvent().waitForFinished();
		    double end = MPI_Wtime();
		    /* notification callback for simulation step currentStep
		     * called every notifyPeriod steps */
		    fprintf(stderr, "Jai's plugin notified %lf.\n", end-start);
		    if (currentStep == INITIAL_STEP) {
			// define field vars.
			// define species vars.
			// define particle vars. 
		    }
		    // pull data from GPU
		    // write it.
		    
		    if (currentStep == FINAL_STEP) {
			//finalize.
		    }
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
	    
	    ThreadParams2 mThreadParams;
	    MappingDesc *cellDescription;
	    //uint32_t notifyPeriod;
	    std::string filename;
	    std::string checkpointFilename;
	    std::string restartFilename;
	    std::string checkpointDirectory;
	    
	    /* select MPI method, #OSTs and #aggregators */
	    std::string mpiTransportParams;
	    
	    uint32_t restartChunkSize;
	    
	    DataSpace<simDim> mpi_pos;
	    DataSpace<simDim> mpi_size;
	};
    } // adiosinsitu
} // picongpu
