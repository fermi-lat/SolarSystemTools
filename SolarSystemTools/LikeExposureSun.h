/**
 * @file LikeExposureSun.h
 * @brief Exposure class for moving sources for use by the SolarSystemTools.
 * @author G. Johannesson
 *
 * $Header: $
 */

#ifndef SolarSystemTools_LikeExposureSun_h
#define SolarSystemTools_LikeExposureSun_h

#include <utility>

#include "tip/tip_types.h"

#include "astro/SkyDir.h"
#include "astro/SolarSystem.h"
#include "healpix/HealpixArray.h"
#include "SolarSystemTools/CosineBinner2D.h"
#include "SolarSystemTools/ExposureSun.h"
#include "map_tools/Exposure.h"


namespace tip {
   class Table;
}

namespace SolarSystemTools {

/**
 * @class LikeExposureSun
 *
 * @brief Class to aid in computing an exposure time hypercube that
 * includes the ROI time range cuts and GTIs.
 *
 * It is pixelated using Healpix binning
 *
 * @author J. Chiang
 *
 * $Header: $
 */

class LikeExposureSun : public ExposureSun {

public:

   LikeExposureSun(double skybin, double costhetabin, double costhetabinSun,
                const std::vector< std::pair<double, double> > & timeCuts,
                const std::vector< std::pair<double, double> > & gtis,
                double zenmax=180.);

   LikeExposureSun(const LikeExposureSun & other);

   void load(const tip::Table * tuple, bool verbose=true);

   tip::Index_t numIntervals() const {
      return m_numIntervals;
   }

   /// @brief Normally one would re-implement the
   /// map_tools::Exposure::write(...) member function from the base
   /// class, but it is not virtual, so we add this method instead to
   /// avoid possible confusion if these classes are used
   /// polymorphically.
   void writeFile(const std::string & outfile) const;

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

protected:

   LikeExposureSun & operator=(const LikeExposureSun &) {
      return *this;
   }

private:

   double m_costhetabin;
   double m_costhetabinSun;

	 ExposureSun * m_weightedExposure;

   const std::vector< std::pair<double, double> > & m_timeCuts;
   const std::vector< std::pair<double, double> > & m_gtis;

	 /// Start of mission in Julian date
		static const double s_mjd_missionStart;
		astro::SolarSystem m_solar_dir;
   /// Minimum time to be considered given GTIs (MET s)
   double m_tmin;

   /// Maximum time to be considered given GTIs (MET s)
   double m_tmax;

   /// Number of FT2 intervals that have been loaded.
   tip::Index_t m_numIntervals;

   static bool overlaps(const std::pair<double, double> & interval1,
                        std::pair<double, double> & interval2);

   void writeFilename(const std::string & outfile) const;

   void writeLivetimes(const std::string & outfile,
                       const ExposureSun * self=0,
                       const std::string & extname="EXPOSURESUN") const;

   void writeCosbins(const std::string & outfile) const;

   void setCosbinsFieldFormat(const std::string & outfile,
                              const ExposureSun * self,
                              const std::string & extname) const;

   void fitsReportError(int status, const std::string & routine) const;

   void computeCosbins(std::vector<double> & mubounds, std::vector<double> &muSunbounds) const;

};

} // namespace SolarSystemTools

#endif // SolarSystemTools_LikeExposureSun_h
