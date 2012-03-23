/**
 * @file HealpixExposureSun.h
 * @brief All-Sky exposure map for use by SourceMap for DiffuseSource 
 * integrations
 * @author G. Johannesson
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/SolarSystemTools/HealpixExposureSun.h,v 1.2 2012/03/21 22:50:19 gudlaugu Exp $
 */

#ifndef SolarSystemTools_HealpixExposureSun_h
#define SolarSystemTools_HealpixExposureSun_h

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

namespace healpix {
	template<typename T>
	class HealpixArray;
}

namespace SolarSystemTools {

   class Observation;

/**
 * @class HealpixExposureSun
 * @brief This class encapsulates the calculation of and access to 
 * the integral of the effective area over live time.
 *
 * @author G. Johannesson
 */

class HealpixExposureSun {

public:

	 class Fun {
		 public:
			 virtual double operator() (double thetasun) const = 0;
			 virtual ~Fun() {}
	 };

   HealpixExposureSun();

   HealpixExposureSun(const std::vector<double> & energies,
                  const Observation & observation,
                  const st_app::AppParGroup * pars=0);

   HealpixExposureSun(const std::string & filename);

   /// @return ExposureSun (effective area integrated over time) (cm^2-s)
   /// @param energy True photon energy (MeV)
   /// @param ra Right Ascension of desired sky location (degrees)
   /// @param dec Declination of desired sky location (degrees)
   /// @param thetasun the angle from sun
   double operator()(double energy, double ra, double dec, double thetasun) const;

   /// @return Integral of function f integrated over solar angle
   /// @param energy True photon energy (MeV)
   /// @param ra Right Ascension of desired sky location (degrees)
   /// @param dec Declination of desired sky location (degrees)
   double integrate(double energy, double ra, double dec, const Fun &f) const;

   void writeOutput(const std::string & filename) const;

   const std::vector<double> & energies() const {
      return m_energies;
   }

   const std::vector<double> & thetasun() const {
      return m_thetasun;
   }

   void setBoundaryFlag(bool enforce_boundaries) {
      m_enforce_boundaries = enforce_boundaries;
   }

protected:

// Disable copy constructor and copy assignment operator
   HealpixExposureSun(const HealpixExposureSun &) {
      throw std::runtime_error("HealpixExposureSun copy constructor is disabled");
   }

   HealpixExposureSun & operator=(const HealpixExposureSun &) {
      throw std::runtime_error("HealpixExposureSun copy assignment operator "
                               "is disabled");
      return *this;
   }

   void setMapGeometry(const st_app::AppParGroup & pars);

   void setMapGeometry();

private:

   const Observation * m_observation;

	 // The array is sparse in thetasun only
	 class pixel : public std::vector<float> {
		 public:
			 float & value (size_t index);
			 float value (size_t index) const;
			 
			 std::vector<size_t> indices() const;
			 std::vector<size_t> sparseIndices() const;
			 
			 static void setStride(size_t stride);

		 private:
			 // The data is stored sparse in first index i
			 // index = i*stride + j
	     typedef std::vector< std::pair<size_t, size_t> > pixelIndex;
			 pixelIndex m_ifaketoreal, m_irealtofake;
			 
			 static size_t m_stride;

		   static bool less_than(const std::pair<size_t,size_t> &a, const std::pair<size_t,size_t> &b) { return a.first < b.first; }
	 };
	 healpix::HealpixArray<pixel> m_exposureMap;

   std::vector<double> m_energies;
   std::vector<double> m_thetasun;

   double m_costhmin;
   double m_costhmax;

   bool m_enforce_boundaries;

   void setCosThetaBounds(const st_app::AppParGroup & pars);

   void computeMap();

};

} // namespace SolarSystemTools

#endif // SolarSystemTools_HealpixExposureSun_h
