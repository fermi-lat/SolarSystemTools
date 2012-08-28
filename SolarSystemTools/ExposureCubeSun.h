/**
 * @file ExposureCubeSun.h
 * @brief Exposure time hypercube.
 * @author G. Johannesson <gudlaugu@glast2.stanford.edu>
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/SolarSystemTools/ExposureCubeSun.h,v 1.9 2012/08/24 13:48:46 gudlaugu Exp $
 */

#ifndef SolarSystemTools_ExposureCubeSun_h
#define SolarSystemTools_ExposureCubeSun_h

#include <stdexcept>
#include <vector>

#include "facilities/Util.h"

#include "astro/SkyDir.h"

#include "Likelihood/RoiCuts.h"

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
 * $Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/SolarSystemTools/ExposureCubeSun.h,v 1.9 2012/08/24 13:48:46 gudlaugu Exp $
 */

class ExposureCubeSun {

public:

   ExposureCubeSun() : m_exposure(0), m_weightedExposure(0), 
                    m_efficiencyFactor(0),
                    m_haveFile(false), m_fileName(""),
                    m_hasPhiDependence(false),
										m_distCosCut(-1),
										m_tstart(0), m_tstop(0),
		   m_timeDist(0), m_time(0),
       m_timeCuts(std::vector<std::pair<double, double> >(0)),
       m_gtis(std::vector<std::pair<double, double> >(0))
	{}

   ExposureCubeSun(double skybin, double costhetabin, double thetabinSun, double thetamax, double powerbinsun,
                const std::vector< std::pair<double, double> > & timeCuts,
                const std::vector< std::pair<double, double> > & gtis,
                astro::SolarSystem::Body=astro::SolarSystem::SUN,
								double zenmax=180.);

   ExposureCubeSun(const ExposureCubeSun & other);

   ~ExposureCubeSun() {
      delete m_exposure;
      delete m_weightedExposure;
      delete m_efficiencyFactor;
   }

   void readExposureCubeSun(std::string filename);

   double livetime(const astro::SkyDir & dir, double costheta,
                   double thetasun, double phi=1) const;

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

	 bool hasCosthetasun(const astro::SkyDir & dir, double costhetasun) const {
		 return m_exposure->hasCosthetasun(dir, costhetasun);
	 }

	 healpix::Healpix getHealpix() const {
		 return m_exposure->data().healpix();
	 }

	 //Returns the boundaries
	 void thetaBinsSun(std::vector<double> &thSunbounds) const;

	 size_t nthetaBinsSun() const {
		  return CosineBinner2D::nbins2();
	 }

	 ExposureCubeSun & operator += (const ExposureCubeSun &other);

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

   void load(const tip::Table * tuple, bool verbose=true);

   tip::Index_t numIntervals() const {
      return m_numIntervals;
   }

	 double distCosCut() const { return m_distCosCut; }
	 double avgDist() const;

	 double tstart() const { return m_tstart; }
	 double tstop() const { return m_tstop; }

   /// @brief Normally one would re-implement the
   /// map_tools::Exposure::write(...) member function from the base
   /// class, but it is not virtual, so we add this method instead to
   /// avoid possible confusion if these classes are used
   /// polymorphically.
   void writeFile(const std::string & outfile, double start, double stop, const Likelihood::RoiCuts &cuts) const;

   /// @param start MET start time of interval (seconds)
   /// @param stop MET stop time of interval (seconds)
   /// @param timeCuts Time range cuts
   /// @param gtis Good Time Intervals
   /// @param fraction Fraction of the interval to use in exposure
   ///        calculation.  This is a return value.
   static bool 
   acceptInterval(double start, double stop, 
                  const std::vector< std::pair<double, double> > & timeCuts,
                  const std::vector< std::pair<double, double> > & gtis,
                  double & fraction);

   static double overlap(const std::pair<double, double> & interval1,
                         const std::pair<double, double> & interval2);


	 static std::string bodyToString(astro::SolarSystem::Body body);
	 static astro::SolarSystem::Body stringToBody(const std::string &body);

protected:

   ExposureCubeSun & operator=(const ExposureCubeSun &) {
      return *this;
   }


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
	 void writeKeywords(const std::string &outfile, const std::string &extname, const Likelihood::RoiCuts &cuts) const;
	 void readKeywords(const std::string &outfile, const std::string &extname);
   double m_costhetabin;
   double m_thetabin;
   double m_thetamax;

   ExposureSun * m_exposure;
   ExposureSun * m_weightedExposure;

   irfInterface::IEfficiencyFactor * m_efficiencyFactor;

   bool m_haveFile;

   std::string m_fileName;

   bool m_hasPhiDependence;

   const std::vector< std::pair<double, double> > & m_timeCuts;
   const std::vector< std::pair<double, double> > & m_gtis;

	 /// Start of mission in Julian date
		static const double s_mjd_missionStart;
		astro::SolarSystem m_source_dir;
		astro::SolarSystem::Body m_body;

		/// cosine of maximum angle to apply the distance weighting
		double m_distCosCut;
		/// time * distance to the source (time in seconds, distance 
		/// in light seconds).  Used to calculate the time weighted average distance.
		double m_timeDist;
		/// Total livetime accumulated
		double m_time;

   /// Minimum time to be considered given GTIs (MET s)
   double m_tmin;

   /// Maximum time to be considered given GTIs (MET s)
   double m_tmax;

	 mutable double m_tstart, m_tstop;

   /// Number of FT2 intervals that have been loaded.
   tip::Index_t m_numIntervals;

   static bool overlaps(const std::pair<double, double> & interval1,
                        std::pair<double, double> & interval2);

   void writeBins(const std::string & outfile) const;

   void setCosbinsFieldFormat(const std::string & outfile,
                              const ExposureSun * self,
                              const std::string & extname) const;

   void fitsReportError(int status, const std::string & routine) const;

   void computeBins(std::vector<double> & mubounds, std::vector<double> &thSunbounds) const;

   bool phiDependence(const std::string & filename) const;

};

} // namespace SolarSystemTools

#endif // SolarSystemTools_ExposureCubeSun_h
