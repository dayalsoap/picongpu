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
	
	class ADIOSInSitu : public ILightweightPlugin
	{
	private:

	    template<typename UnitType>
	    static std::vector<float_64> createUnit(UnitType unit, uint32_t numComponents)
		{
		    std::vector<float_64> tmp(numComponents);
		    for (uint32_t i = 0; i < numComponents; ++i)
			tmp[i] = unit[i];
		    return tmp;
		}

	    /**
	     * Write calculated fields to adios file.
	     */
	    template< typename T >
	    struct GetFields
	    {
	    private:
		typedef typename T::ValueType ValueType;
		typedef typename GetComponentsType<ValueType>::type ComponentType;

	    public:

		HDINLINE void operator()(ThreadParams* params)
		    {
#ifndef __CUDA_ARCH__
			DataConnector &dc = Environment<simDim>::get().DataConnector();
			
			T* field = &(dc.getData<T > (T::getName()));
			params->gridLayout = field->getGridLayout();

			PICToAdios<ComponentType> adiosType;
			writeField(params,
				   sizeof(ComponentType),
				   adiosType.type,
				   GetNComponents<ValueType>::value,
				   T::getName(),
				   field->getHostDataBox().getPointer());

			dc.releaseData(T::getName());
#endif
		    }

	    }; // struct GetFields
	    
    /** Calculate FieldTmp with given solver and particle species
     * and write them to adios.
     *
     * FieldTmp is calculated on device and than dumped to adios.
     */
    template< typename Solver, typename Species >
    struct GetFields<FieldTmpOperation<Solver, Species> >
    {

        /*
         * This is only a wrapper function to allow disable nvcc warnings.
         * Warning: calling a __host__ function from __host__ __device__
         * function.
         * Use of PMACC_NO_NVCC_HDWARNING is not possible if we call a virtual
         * method inside of the method were we disable the warnings.
         * Therefore we create this method and call a new method were we can
         * call virtual functions.
         */
        PMACC_NO_NVCC_HDWARNING
        HDINLINE void operator()(ThreadParams* tparam)
        {
            this->operator_impl(tparam);
        }
    private:
        typedef typename FieldTmp::ValueType ValueType;
        typedef typename GetComponentsType<ValueType>::type ComponentType;

        /** Create a name for the adios identifier.
         */
        static std::string getName()
        {
            std::stringstream str;
            str << Solver().getName();
            str << "_";
            str << Species::FrameType::getName();
            return str.str();
        }

        HINLINE void operator_impl(ThreadParams* params)
        {
            DataConnector &dc = Environment<>::get().DataConnector();

            /*## update field ##*/

            /*load FieldTmp without copy data to host*/
            FieldTmp* fieldTmp = &(dc.getData<FieldTmp > (FieldTmp::getName(), true));
            /*load particle without copy particle data to host*/
            Species* speciesTmp = &(dc.getData<Species >(Species::FrameType::getName(), true));

            fieldTmp->getGridBuffer().getDeviceBuffer().setValue(ValueType::create(0.0));
            /*run algorithm*/
            fieldTmp->computeValue < CORE + BORDER, Solver > (*speciesTmp, params->currentStep);

            EventTask fieldTmpEvent = fieldTmp->asyncCommunication(__getTransactionEvent());
            __setTransactionEvent(fieldTmpEvent);
            /* copy data to host that we can write same to disk*/
            fieldTmp->getGridBuffer().deviceToHost();
            dc.releaseData(Species::FrameType::getName());
            /*## finish update field ##*/

            const uint32_t components = GetNComponents<ValueType>::value;
            PICToAdios<ComponentType> adiosType;

            params->gridLayout = fieldTmp->getGridLayout();
            /*write data to ADIOS file*/
            writeField(params,
                       sizeof(ComponentType),
                       adiosType.type,
                       components,
                       getName(),
                       fieldTmp->getHostDataBox().getPointer());

            dc.releaseData(FieldTmp::getName());
	    
        }

    };

	    static void defineFieldVar(ThreadParams* params,
				       uint32_t nComponents, ADIOS_DATATYPES adiosType, const std::string name,
				       std::vector<float_64> unit)
		{
		    const std::string name_lookup_tpl[] = {"x", "y", "z", "w"};

		    for (uint32_t c = 0; c < nComponents; c++)
		    {
			std::stringstream datasetName;
			datasetName << params->adiosBasePath << ADIOS_PATH_FIELDS << name;
			if (nComponents > 1)
			    datasetName << "/" << name_lookup_tpl[c];

			/* define adios var for field, e.g. field_FieldE_y */
			const char* path = NULL;
			int64_t adiosFieldVarId = defineAdiosVar<simDim>(
			    params->adiosGroupHandle,
			    datasetName.str().c_str(),
			    path,
			    adiosType,
			    params->fieldsSizeDims,
			    params->fieldsGlobalSizeDims,
			    params->fieldsOffsetDims,
			    true,
			    params->adiosCompression);
	    
			params->adiosFieldVarIds.push_back(adiosFieldVarId);

			/* already add the sim_unit attribute so `adios_group_size` calculates
			 * the reservation for the buffer correctly */
			AdiosDoubleType adiosDoubleType;

			ADIOS_CMD(adios_define_attribute(params->adiosGroupHandle,
							 "sim_unit", datasetName.str().c_str(), adiosDoubleType.type,
							 flt2str(unit.at(c)).c_str(), ""));
		    }
		}

    /**
     * Collect field sizes to set adios group size.
     */
	    template< typename T >
	    struct CollectFieldsSizes
	    {
	    public:
		typedef typename T::ValueType ValueType;
		typedef typename T::UnitValueType UnitType;
		typedef typename GetComponentsType<ValueType>::type ComponentType;
		
		static std::vector<float_64> getUnit()
		    {
			UnitType unit = T::getUnit();
			return createUnit(unit, T::numComponents);
		    }

		HDINLINE void operator()(ThreadParams* params)
		    {
#ifndef __CUDA_ARCH__
			const uint32_t components = T::numComponents;

			// adios buffer size for this dataset (all components)
			uint64_t localGroupSize =
			    params->window.localDimensions.size.productOfComponents() *
			    sizeof(ComponentType) *
			    components;

			params->adiosGroupSize += localGroupSize;

			PICToAdios<ComponentType> adiosType;
			defineFieldVar(params, components, adiosType.type, T::getName(), getUnit());
#endif
		    }
	    };

	        /**
     * Collect field sizes to set adios group size.
     * Specialization.
     */
	    template< typename Solver, typename Species >
	    struct CollectFieldsSizes<FieldTmpOperation<Solver, Species> >
	    {
	    public:

		PMACC_NO_NVCC_HDWARNING
		HDINLINE void operator()(ThreadParams* tparam)
		    {
			this->operator_impl(tparam);
		    }

	    private:
		typedef typename FieldTmp::ValueType ValueType;
		typedef typename FieldTmp::UnitValueType UnitType;
		typedef typename GetComponentsType<ValueType>::type ComponentType;

		/** Create a name for the adios identifier.
		 */
		static std::string getName()
		    {
			std::stringstream str;
			str << Solver().getName();
			str << "_";
			str << Species::FrameType::getName();
			return str.str();
		    }

		/** Get the unit for the result from the solver*/
		static std::vector<float_64> getUnit()
		    {
			UnitType unit = FieldTmp::getUnit<Solver>();
			const uint32_t components = GetNComponents<ValueType>::value;
			return createUnit(unit, components);
		    }

		HINLINE void operator_impl(ThreadParams* params)
		    {
			const uint32_t components = GetNComponents<ValueType>::value;

			// adios buffer size for this dataset (all components)
			uint64_t localGroupSize =
			    params->window.localDimensions.size.productOfComponents() *
			    sizeof(ComponentType) *
			    components;

			params->adiosGroupSize += localGroupSize;

			PICToAdios<ComponentType> adiosType;
			defineFieldVar(params, components, adiosType.type, getName(), getUnit());
		    }

	    };
	    
	public:
	    ADIOSInSitu() :
		filename("adiospicongpu.bp"),
		flexpathTransportParams("QUEUE_SIZE=10"),
		notifyPeriod(0),
		lastSpeciesSyncStep(PMacc::traits::limits::Max<uint32_t>::value)
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

    static void writeField(ThreadParams *params, const uint32_t sizePtrType,
                           ADIOS_DATATYPES adiosType,
                           const uint32_t nComponents, const std::string name,
                           void *ptr)
    {
        log<picLog::INPUT_OUTPUT > ("ADIOS: write field: %1% %2% %3%") %
            name % nComponents % ptr;

        /* data to describe source buffer */
        GridLayout<simDim> field_layout = params->gridLayout;
        DataSpace<simDim> field_full = field_layout.getDataSpace();
        DataSpace<simDim> field_no_guard = params->window.localDimensions.size;
        DataSpace<simDim> field_guard = field_layout.getGuard() + params->localWindowToDomainOffset;

        /* write the actual field data */
        for (uint32_t d = 0; d < nComponents; d++)
        {
            const size_t plane_full_size = field_full[1] * field_full[0] * nComponents;
            const size_t plane_no_guard_size = field_no_guard[1] * field_no_guard[0];

            /* copy strided data from source to temporary buffer
             *
             * \todo use d1Access as in `include/plugins/hdf5/writer/Field.hpp`
             */
            const int maxZ = simDim == DIM3 ? field_no_guard[2] : 1;
            const int guardZ = simDim == DIM3 ? field_guard[2] : 0;
            for (int z = 0; z < maxZ; ++z)
            {
                for (int y = 0; y < field_no_guard[1]; ++y)
                {
                    const size_t base_index_src =
                                (z + guardZ) * plane_full_size +
                                (y + field_guard[1]) * field_full[0] * nComponents;

                    const size_t base_index_dst =
                                z * plane_no_guard_size +
                                y * field_no_guard[0];

                    for (int x = 0; x < field_no_guard[0]; ++x)
                    {
                        size_t index_src = base_index_src + (x + field_guard[0]) * nComponents + d;
                        size_t index_dst = base_index_dst + x;

                        params->fieldBfr[index_dst] = ((float_32*)ptr)[index_src];
                    }
                }
            }

            /* Write the actual field data. The id is on the front of the list. */
            if (params->adiosFieldVarIds.empty())
                throw std::runtime_error("Cannot write field (var id list is empty)");

            int64_t adiosFieldVarId = *(params->adiosFieldVarIds.begin());
            params->adiosFieldVarIds.pop_front();
            ADIOS_CMD(adios_write_byid(params->adiosFileHandle, adiosFieldVarId, params->fieldBfr));
        }
    }
	    
	    ThreadParams mThreadParams;
	    MappingDesc *cellDescription;
	    //uint32_t notifyPeriod;
	    std::string filename;
	    std::string checkpointDirectory;
	    
	    /* Flexpath Queuesize */
	    std::string flexpathTransportParams;	   
	    uint32_t lastSpeciesSyncStep;
	    
	    DataSpace<simDim> mpi_pos;
	    DataSpace<simDim> mpi_size;
	};
    } // adiosinsitu
} // picongpu
