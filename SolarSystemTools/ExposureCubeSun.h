/**
 * @file ExposureCubeSun.h
 * @brief Exposure time hypercube.
 * @author G. Johannesson <gudlaugu@glast2.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/SolarSystemTools/ExposureCubeSun.h,v 1.2 2012/02/13 22:24:36 gudlaugu Exp $
 */

#ifndef SolarSystemTools_ExposureCubeSun_h
#define SolarSystemTools_ExposureCubeSun_h

#include <stdexcept>
#include <vector>

#include "facilities/Util.h"

#include "astro/SkyDir.h"

#include "irfInterface/IEfficiencyFactor.h"

#include "SolarSystemTools/ExposureSun.h"

#include "SolarSystemTools/CosineBinner2D.h"

namespace SolarSystemTools {

   class Observation;

/**
 * @class ExposureCubeSun
 * @brief Exposure time as a function of sky position and inclination
 *        wrt the instrument z-axis
 *
 * @author G. Johannesson
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/SolarSystemTools/ExposureCubeSun.h,v 1.2 2012/02/13 22:24:36 gudlaugu Exp $
 */

class ExposureCubeSun {

public:

   ExposureCubeSun() : m_exposure(0), m_weightedExposure(0), 
                    m_efficiencyFactor(0),
                    m_haveFile(false), m_fileName(""),
                    m_hasPhiDependence(false) {}

   ExposureCubeSun(const ExposureCubeSun & other);

   ~ExposureCubeSun() {
      delete m_exposure;
      delete m_weightedExposure;
      delete m_efficiencyFactor;
   }

   void readExposureCubeSun(std::string filename);

   double livetime(const astro::SkyDir & dir, double costheta,
                   double costhetasun, double phi=1) const;

   void setEfficiencyFactor(const irfInterface::IEfficiencyFactor * eff) {
      if (eff) {
         m_efficiencyFactor = eff->clone();
      }
   }

#ifndef SWIG
   template<class T>
   double value(const astro::SkyDir & dir, double costhetasun, const T & aeff,
                bool weighted_lt=false) const {
      if (m_hasPhiDependence) {
         AeffWrapper<T> myAeff(aeff);
         if (weighted_lt && m_weightedExposure) {
            return m_weightedExposure->integral(dir, costhetasun, myAeff);
         }
         return m_exposure->integral(dir, costhetasun, myAeff);
      }
      if (weighted_lt && m_weightedExposure) {
         return m_weightedExposure->operator()(dir, costhetasun, aeff);
      }
      return (*m_exposure)(dir, costhetasun, aeff);
   }

   // Compute the exposure with trigger rate- and energy-dependent
   // efficiency corrections.
   template<class T>
   double value(const astro::SkyDir & dir, double costhetasun, const T & aeff, 
                double energy) const {
      double factor1(1), factor2(0);
      if (m_efficiencyFactor) {
         m_efficiencyFactor->getLivetimeFactors(energy, factor1, factor2);
      }
      double value1(value(dir, costhetasun, aeff));
      double value2(0);
      double exposure(factor1*value1);
      if (factor2 != 0) {
         value2 = value(dir, costhetasun, aeff, true);
         exposure += factor2*value2;
      }
      if (exposure < 0) {
         throw std::runtime_error("ExposureCubeSun::value: exposure < 0");
      }
      return exposure;
   }
#endif // SWIG

   bool haveFile() const {
      return m_haveFile;
   }

   const std::string & fileName() const {
      return m_fileName;
   }

   bool hasPhiDependence() const {
      return m_hasPhiDependence;
   }

	 healpix::Healpix getHealpix() const {
		 return m_exposure->data().healpix();
	 }

	 //Returns the boundaries
	 void costhetaBinsSun(std::vector<double> &muSunbounds) const;

	 size_t ncosthetaBinsSun() const {
		  return CosineBinner2D::nbins2();
	 }

   class Aeff {
   public:
      Aeff(double energy, const Observation & observation); 
      virtual ~Aeff() {}
      virtual double operator()(double cosTheta, double phi=0) const;
      virtual double integral(double cosTheta, double phi=0) const {
         return operator()(cosTheta, phi);
      }
   protected:
      double m_energy;
      const Observation & m_observation;
   };

private:

#ifndef SWIG
// healpix::CosineBinner::integral assumes that phi is in radians instead of
// degrees, contrary to established conventions.  This class wraps the 
// integral(costh, phi) method to do the conversion.
   template <class Aeff>
   class AeffWrapper {
   public:
      AeffWrapper(const Aeff & aeff) : m_aeff(aeff) {}
      double integral(double costh, double phi) const {
         phi *= 180./M_PI;
         return m_aeff.integral(costh, phi);
      }
   private:
      const Aeff & m_aeff;
   };
#endif

   ExposureSun * m_exposure;
   ExposureSun * m_weightedExposure;

   irfInterface::IEfficiencyFactor * m_efficiencyFactor;

   bool m_haveFile;

   std::string m_fileName;

   bool m_hasPhiDependence;

   bool phiDependence(const std::string & filename) const;

};

} // namespace SolarSystemTools

#endif // SolarSystemTools_ExposureCubeSun_h
