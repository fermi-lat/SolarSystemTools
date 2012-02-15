/**
 * @file SolarTemplate.h
 * @brief Calculates an average template for the sun given a healpix exposure
 * map binnen in solar angles and a binned exposure map
 * @author G. Johannesson
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/SolarSystemTools/SolarTemplate.h,v 1.2 2012/02/13 22:24:36 gudlaugu Exp $
 */

#ifndef SolarSystemTools_SolarTemplate_h
#define SolarSystemTools_SolarTemplate_h

#include <stdexcept>
#include <string>
#include <vector>

#include "SolarSystemTools/HealpixExposureSun.h"

namespace st_app {
   class AppParGroup;
}

namespace astro {
   class SkyProj;
}

namespace Likelihood {
	class CountsMap;
	 class BinnedExposure;
}

namespace SolarSystemTools {

   class Observation;
	 class SolarProfile;

/**
 * @class SolarTemplate
 * @brief This class encapsulates the calculation of a solar template for use
 * in likelihood analysis
 *
 * @author G. Johannesson
 */

class SolarTemplate {

public:

   SolarTemplate(const Likelihood::CountsMap & cmap, 
                  const HealpixExposureSun & hpexpsun, 
                  const Likelihood::BinnedExposure & avgexp, 
									const SolarProfile & prof,
                  bool useEbounds=true,
                  const st_app::AppParGroup * pars=0);

   SolarTemplate(const std::vector<double> & energies,
                  const HealpixExposureSun & hpexpsun, 
                  const Likelihood::BinnedExposure & avgexp, 
									const SolarProfile & prof,
                  const st_app::AppParGroup * pars=0);

   ~SolarTemplate();

   void writeOutput(const std::string & filename) const;

   const std::vector<double> & energies() const {
      return m_energies;
   }

protected:

// Disable copy constructor and copy assignment operator
   SolarTemplate & operator=(const SolarTemplate &) {
      throw std::runtime_error("SolarTemplate copy assignment operator "
                               "is disabled");
      return *this;
   }

   void setMapGeometry(const Likelihood::CountsMap & cmap);

   void setMapGeometry(const st_app::AppParGroup & pars);

   void setMapGeometry();

private:

	 const HealpixExposureSun & m_expsun;
	 const Likelihood::BinnedExposure & m_avgexp;
	 const SolarProfile & m_profile;

	 std::vector<float> m_template;

   std::vector<double> m_energies;

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

   void computeMap();

	 class ProfFun : public HealpixExposureSun::Fun {
		 public:
			 ProfFun ( const std::vector<double> & costhsun, 
					       const std::vector<double> & energies, 
					       const std::vector<double> & cache, size_t ie );
			 virtual double operator() (double costhetasun) const;
		 private:
			 const std::vector<double> & m_costhsun, m_energies, m_cache;
			 const size_t m_ie;
	 };

};

} // namespace SolarSystemTools

#endif // SolarSystemTools_SolarTemplate_h
