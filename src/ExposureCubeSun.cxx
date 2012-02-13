/**
 * @file ExposureCubeSun.cxx
 * @brief Implementation for ExposureCubeSun wrapper class of SolarSystemTools::Exposure
 * @author G. Johannesson
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/SolarSystemTools/src/ExposureCubeSun.cxx,v 1.1.1.1 2012/02/11 02:26:40 gudlaugu Exp $
 */

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "SolarSystemTools/ExposureCubeSun.h"
#include "SolarSystemTools/Observation.h"

namespace SolarSystemTools {

ExposureCubeSun::ExposureCubeSun(const ExposureCubeSun & other) 
   : m_exposure(new ExposureSun(*(other.m_exposure))),
     m_weightedExposure(0), m_efficiencyFactor(0),
     m_haveFile(other.m_haveFile),
     m_fileName(other.m_fileName), 
     m_hasPhiDependence(other.m_hasPhiDependence) {
   if (other.m_weightedExposure) {
      m_weightedExposure = new ExposureSun(*(other.m_weightedExposure));
   }
   if (other.m_efficiencyFactor) {
      m_efficiencyFactor = other.m_efficiencyFactor->clone();
   }
}

void ExposureCubeSun::readExposureCubeSun(std::string filename) {
   facilities::Util::expandEnvVar(&filename);
   m_fileName = filename;
   m_exposure = new ExposureSun(filename);
   try {
      m_weightedExposure = new ExposureSun(filename,
                                                   "WEIGHTED_EXPOSURESUN");
   } catch(tip::TipException &) {
      m_weightedExposure = 0;
   }
   m_haveFile = true;
   m_hasPhiDependence = phiDependence(filename);
}

double ExposureCubeSun::livetime(const astro::SkyDir & dir,
                              double costheta, double costhetasun, double phi) const {
   return m_exposure->data()[dir](costheta, costhetasun, phi);
}

void ExposureCubeSun::
costhetaBinsSun(std::vector<double> &muSunbounds) const {
   bool sqrtbins(CosineBinner2D::thetaBinning() == "SQRT(1-COSTHETA)");
   double cosminSun(CosineBinner2D::cosmin2());
   size_t nbinsSun(CosineBinner2D::nbins2());
   muSunbounds.clear();
//   for (int i(nbins); i >= 0; i--) {
   for (size_t i(0); i < nbinsSun+1; i++) {
      double factor(static_cast<double>(i)/nbinsSun);
      if (sqrtbins) {
         factor *= factor;
      }
      muSunbounds.push_back(1. - factor*(1. - cosminSun));
   }
}

bool ExposureCubeSun::phiDependence(const std::string & filename) const {
   const tip::Table * table 
      = tip::IFileSvc::instance().readTable(filename, "EXPOSURESUN");
   const tip::Header & header(table->getHeader());
   long nphibins;
   try {
      header["PHIBINS"].get(nphibins);
   } catch (tip::TipException &) {
      nphibins = 0;
   }
   return nphibins > 0;
}

ExposureCubeSun::Aeff::Aeff(double energy, 
                         const Observation & observation) 
   : m_energy(energy), m_observation(observation) {
// Turn off phi-dependence if omitted from livetime cube.
   bool phi_dependence(m_observation.expCubeSun().hasPhiDependence());
   if (!phi_dependence) {
      std::map<unsigned int, irfInterface::Irfs *>::const_iterator respIt 
         = m_observation.respFuncs().begin();
      for ( ; respIt != m_observation.respFuncs().end(); ++respIt) {
         respIt->second->aeff()->setPhiDependence(false);
      }
   }
}

double ExposureCubeSun::Aeff::operator()(double cosTheta, double phi) const {
   double inclination = acos(cosTheta)*180./M_PI;
// Check if we need to truncate the inclination past 70 deg as in 
// MeanPsf::Aeff::operator()(...)
   // if (inclination > 70.) {
   //    return 0;
   // }
   std::map<unsigned int, irfInterface::Irfs *>::const_iterator respIt 
      = m_observation.respFuncs().begin();
	 double aeff_val(0);
   for ( ; respIt != m_observation.respFuncs().end(); ++respIt) {
         irfInterface::IAeff * aeff = respIt->second->aeff();
         aeff_val += aeff->value(m_energy, inclination, phi);
   }
   return aeff_val;
   return 0;
}

} // namespace SolarSystemTools
