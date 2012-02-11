/**
 * @file Observation.h
 * @brief Extending the Likelihood Observation to work with SolarSystemTools
 * @author G. Johannesson
 *
 * $Header: $
 */

#ifndef SolarSystemTools_Observation_h
#define SolarSystemTools_Observation_h

#include "Likelihood/Observation.h"
#include "SolarSystemTools/ExposureCubeSun.h"

namespace SolarSystemTools {

/**
 * @class Observation
 * @brief A container class composed of all of the data-related
 * classes associated with a particular observation.
 *
 * The contained classes were all originally implemented as
 * Singletons.  This class encapsulates the observation information
 * comprising those classes, providing a single access point whilst
 * allowing multiple observations to be considered within the same
 * program instance.
 *
 * @author G. Johannesson
 *
 * $Header: $
 */

class Observation : public Likelihood::Observation {

public:

//    Observation() : m_respFuncs(0), m_scData(0), m_roiCuts(0), m_expCube(0),
//       m_expMap(0), m_eventCont(0) {}

   Observation(Likelihood::ResponseFunctions * respFuncs, 
               Likelihood::ScData * scData,
               Likelihood::RoiCuts * roiCuts,
               Likelihood::ExposureCube * expCube,
               ExposureCubeSun * expCubeSun,
               Likelihood::ExposureMap * expMap,
               Likelihood::EventContainer * eventCont) :
		 Likelihood::Observation(respFuncs,scData,roiCuts,expCube,expMap,eventCont),
      m_expCubeSun(expCubeSun) {}

	 Observation(Likelihood::Observation &obs, ExposureCubeSun *expCubeSun) :
		 Likelihood::Observation(obs), m_expCubeSun(expCubeSun) {}

   const ExposureCubeSun & expCubeSun() const {
      return *m_expCubeSun;
   }

   ExposureCubeSun & expCubeSun() {
      return *m_expCubeSun;
   }

private:

   /// @todo Assert ownership of these pointers and delete in destructor.
   ExposureCubeSun * m_expCubeSun;

};

} //namespace SolarSystemTools

#endif // SolarSystemTools_Observation_h
