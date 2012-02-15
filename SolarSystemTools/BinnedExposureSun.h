/**
 * @file BinnedExposureSun.h
 * @brief All-Sky exposure map for use by SourceMap for DiffuseSource 
 * integrations
 * @author G. Johannesson
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/SolarSystemTools/BinnedExposureSun.h,v 1.2 2012/02/13 22:24:36 gudlaugu Exp $
 */

#ifndef SolarSystemTools_BinnedExposureSun_h
#define SolarSystemTools_BinnedExposureSun_h

#include <stdexcept>
#include <string>
#include <vector>

#include "SolarSystemTools/ExposureCubeSun.h"

namespace st_app {
   class AppParGroup;
}

namespace astro {
   class SkyProj;
}

namespace Likelihood {
   class CountsMap;
}

namespace SolarSystemTools {

   class Observation;

/**
 * @class BinnedExposureSun
 * @brief This class encapsulates the calculation of and access to 
 * the integral of the effective area over live time.
 *
 * @author G. Johannesson
 */

class BinnedExposureSun {

public:

   BinnedExposureSun();

   BinnedExposureSun(const Likelihood::CountsMap & cmap, 
                  const Observation & observation, 
                  bool useEbounds=true,
                  const st_app::AppParGroup * pars=0);

   BinnedExposureSun(const std::vector<double> & energies,
                  const Observation & observation,
                  const st_app::AppParGroup * pars=0);

   BinnedExposureSun(const std::string & filename);

   ~BinnedExposureSun();

   /// @return ExposureSun (effective area integrated over time) (cm^2-s)
   /// @param energy True photon energy (MeV)
   /// @param ra Right Ascension of desired sky location (degrees)
   /// @param dec Declination of desired sky location (degrees)
   /// @param costhetasun Cosine of the angle from sun
   double operator()(double energy, double ra, double dec, double costhetasun) const;

   void writeOutput(const std::string & filename) const;

   const std::vector<double> & energies() const {
      return m_energies;
   }

   const std::vector<double> & costhetasun() const {
      return m_costhetasun;
   }

   void setBoundaryFlag(bool enforce_boundaries) {
      m_enforce_boundaries = enforce_boundaries;
   }

   class Aeff : public ExposureCubeSun::Aeff {
   public:
      Aeff(double energy, const Observation & observation,
           double costhmin, double costhmax) 
         : ExposureCubeSun::Aeff(energy, observation),
           m_costhmin(costhmin), m_costhmax(costhmax) {}
      virtual double operator()(double cosTheta, double phi=0) const;
   private:
      double m_costhmin;
      double m_costhmax;
			mutable std::map< std::pair<double,double>, double> m_cache;
   };
   

protected:

// Disable copy constructor and copy assignment operator
   BinnedExposureSun(const BinnedExposureSun &) {
      throw std::runtime_error("BinnedExposureSun copy constructor is disabled");
   }

   BinnedExposureSun & operator=(const BinnedExposureSun &) {
      throw std::runtime_error("BinnedExposureSun copy assignment operator "
                               "is disabled");
      return *this;
   }

   void setMapGeometry(const Likelihood::CountsMap & cmap);

   void setMapGeometry(const st_app::AppParGroup & pars);

   void setMapGeometry();

private:

   const Observation * m_observation;

   std::vector<float> m_exposureMap;

   std::vector<double> m_energies;
   std::vector<double> m_costhetasun;

   // Coordinate system parameters to be fed to SkyProj.  These are
   // either set to all-sky values (CAR, CEL) or to match the geometry
   // of the input counts map.
   std::string m_proj_name;
   double m_crpix[2];
   double m_crval[2];
   double m_cdelt[2];
   double m_crota2;
   bool m_isGalactic;

   astro::SkyProj * m_proj;

   std::vector<long> m_naxes;

   double m_costhmin;
   double m_costhmax;

   bool m_enforce_boundaries;

   void setCosThetaBounds(const st_app::AppParGroup & pars);

   void computeMap();

};

} // namespace SolarSystemTools

#endif // SolarSystemTools_BinnedExposure_h
