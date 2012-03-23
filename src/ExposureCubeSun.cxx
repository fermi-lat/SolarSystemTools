/**
 * @file ExposureCubeSun.cxx
 * @brief Implementation for ExposureCubeSun wrapper class of SolarSystemTools::Exposure
 * @author G. Johannesson
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/src/ExposureCubeSun.cxx,v 1.3 2012/03/21 22:50:20 gudlaugu Exp $
 */

#include <iomanip>

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"
#include "st_stream/StreamFormatter.h"
#include "facilities/commonUtilities.h"

#include "fitsio.h"

#include "SolarSystemTools/ExposureCubeSun.h"
#include "SolarSystemTools/Observation.h"

namespace {
   bool compareFirst(const std::pair<double, double> & a, 
                     const std::pair<double, double> & b) {
      return a.first < b.first;
   }
   bool compareSecond(const std::pair<double, double> & a, 
                      const std::pair<double, double> & b) {
      return a.second < b.second;
   }
}

namespace SolarSystemTools {

const double ExposureCubeSun::s_mjd_missionStart(astro::JulianDate::missionStart());

ExposureCubeSun::ExposureCubeSun(const ExposureCubeSun & other) 
   : m_exposure(new ExposureSun(*(other.m_exposure))),
     m_weightedExposure(0), m_efficiencyFactor(0),
     m_haveFile(other.m_haveFile),
     m_fileName(other.m_fileName), 
		 m_gtis(other.m_gtis),
		 m_timeCuts(other.m_timeCuts),
		 m_solar_dir(other.m_solar_dir),
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

ExposureCubeSun::
ExposureCubeSun(double skybin, double costhetabin, double thetabin, double thetamax, double powerbinsun,
             const std::vector< std::pair<double, double> > & timeCuts,
             const std::vector< std::pair<double, double> > & gtis,
             double zenmax)
   : m_efficiencyFactor(0),
	   m_haveFile(false), m_fileName(""),
	   m_hasPhiDependence(false),
     m_costhetabin(costhetabin), m_timeCuts(timeCuts), m_gtis(gtis),
     m_thetabin(thetabin), m_thetamax(thetamax), m_numIntervals(0), 
		 m_solar_dir(astro::SolarSystem::SUN),
     m_exposure(new ExposureSun(skybin, costhetabin, thetabin, thetamax, powerbinsun,
                                                std::cos(zenmax*M_PI/180.))),
     m_weightedExposure(new ExposureSun(skybin, costhetabin, thetabin, thetamax, powerbinsun,
                                                std::cos(zenmax*M_PI/180.)))
     {
   if (!gtis.empty()) {
      for (size_t i = 0; i < gtis.size(); i++) {
         if (i == 0 || gtis.at(i).first < m_tmin) {
            m_tmin = gtis.at(i).first;
         }
         if (i == 0 || gtis.at(i).second > m_tmax) {
            m_tmax = gtis.at(i).second;
         }
      }
   } else {
      throw std::runtime_error("ExposureCubeSun::ExposureCubeSun: GTIs are empty.\n"
                               "Cannot proceed with livetime calculation.");
   }
}

void ExposureCubeSun::load(const tip::Table * scData, bool verbose) {
	if (m_gtis.empty())
      throw std::runtime_error("ExposureCubeSun::ExposureCubeSun: GTIs are empty.\n"
                               "Cannot proceed with livetime calculation.");

   st_stream::StreamFormatter formatter("ExposureCubeSun", "load", 2);
   
   double ra, dec, rax, decx, ra_zenith, dec_zenith, start, stop, livetime;

   tip::Table::ConstIterator it(scData->end());
   tip::ConstTableRecord & row(*it);

// Count the rows within the user selected time interval (m_tmin, m_tmax)
// by counting inwards from the top and bottom of the scData.
   --it;
   tip::Index_t nrows(scData->getNumRecords());
   if (nrows == 0) {
      return;
   }
   for ( ; it != scData->begin(); --it, nrows--) {
      row["stop"].get(stop);
      if (stop < m_tmax) {
         break;
      }
   }

   it = scData->begin();
   for ( ; it != scData->end(); ++it, nrows--) {
      row["start"].get(start);
      if (start > m_tmin) {
         break;
      }
   }

// Reset to the FT2 interval to the one that preceeds the
// user-selected interval, if possible; and set the start time to that
// of the initial row.
   if (it != scData->begin()) {
      --it;
   }
   row["start"].get(start);

// Set the step size for the printing out the little progress dots.
   tip::Index_t istep(nrows/20);
   if (istep == 0) {
      istep = 1;
   }

   for (tip::Index_t irow = 0; it != scData->end() && start < m_tmax;
        ++it, ++irow) {
      if (verbose && (irow % istep) == 0 ) {
         formatter.warn() << "."; 
      }
      row["livetime"].get(livetime);
      row["start"].get(start);
      row["stop"].get(stop);
      double deltat = livetime;
      double fraction;
      if (acceptInterval(start, stop, m_timeCuts, m_gtis, fraction)) {
         row["ra_scz"].get(ra);
         row["dec_scz"].get(dec);
         row["ra_scx"].get(rax);
         row["dec_scx"].get(decx);
         row["ra_zenith"].get(ra_zenith);
         row["dec_zenith"].get(dec_zenith);
				//Usa astro to calculate the direction to the sun at center of bin
				const double mjd = (start+stop)/2./86400. + s_mjd_missionStart;
				astro::SkyDir scsun(m_solar_dir.direction(mjd));
         double weight(livetime/(stop - start));
         if (CosineBinner2D::nphibins() == 0) {
            m_exposure->fill(astro::SkyDir(ra, dec), scsun, astro::SkyDir(ra_zenith, dec_zenith), 
                 deltat*fraction);
            m_weightedExposure->fill(astro::SkyDir(ra, dec),
								                     scsun,
                                     astro::SkyDir(ra_zenith, dec_zenith), 
                                     deltat*fraction*weight);
         } else {
            m_exposure->fill_zenith(astro::SkyDir(ra, dec), astro::SkyDir(rax, decx),
								        scsun,
                        astro::SkyDir(ra_zenith, dec_zenith), 
                        deltat*fraction);
            m_weightedExposure->fill_zenith(astro::SkyDir(ra, dec),
                                            astro::SkyDir(rax, decx),
								                            scsun,
                                            astro::SkyDir(ra_zenith,dec_zenith),
                                            deltat*fraction*weight);
         }
         m_numIntervals++;
      }
   }
   if (verbose) {
      formatter.warn() << "!" << std::endl;
   }
}

void ExposureCubeSun::writeFile(const std::string & outfile) const {
   std::string dataPath = 
      facilities::commonUtilities::getDataPath("SolarSystemTools");
   std::string templateFile = 
      facilities::commonUtilities::joinPath(dataPath, "LivetimeCubeSunTemplate");
   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   fileSvc.createFile(outfile, templateFile);

   writeFilename(outfile);

	 m_exposure->write(outfile, "EXPOSURESUN");

   m_weightedExposure->write(outfile, "WEIGHTED_EXPOSURESUN");

   writeBins(outfile);
}

void ExposureCubeSun::writeFilename(const std::string & outfile) const {
   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   tip::Image * phdu(fileSvc.editImage(outfile, ""));
   phdu->getHeader()["FILENAME"].set(facilities::Util::basename(outfile));
   delete phdu;
}

void ExposureCubeSun::writeBins(const std::string & outfile) const {
   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   tip::Table * table = fileSvc.editTable(outfile, "CTHETABOUNDS");
   table->setNumRecords(CosineBinner2D::nbins());

   tip::Table::Iterator it(table->begin());
   tip::TableRecord & row(*it);
   
   std::vector<double> mubounds, thbounds;
   computeBins(mubounds,thbounds);

   for (size_t i(0); i < mubounds.size() -1; i++, ++it) {
      row["CTHETA_MIN"].set(mubounds.at(i+1));
      row["CTHETA_MAX"].set(mubounds.at(i));
   }
   delete table;

   table = fileSvc.editTable(outfile, "THETASUNBOUNDS");
   table->setNumRecords(CosineBinner2D::nbins2());

   tip::Table::Iterator itSun(table->begin());
   tip::TableRecord & rowSun(*itSun);
   
   for (size_t i(0); i < thbounds.size() -1; i++, ++itSun) {
      rowSun["THETA_MIN"].set(thbounds.at(i)*180/M_PI);
      rowSun["THETA_MAX"].set(thbounds.at(i+1)*180/M_PI);
   }
   delete table;
}
   
void ExposureCubeSun::setCosbinsFieldFormat(const std::string & outfile,
                                         const ExposureSun * self,
                                         const std::string & extname) const {
   int status(0);

   fitsfile * fptr(0);
   std::string extfilename(outfile + "[" + extname + "]");
   fits_open_file(&fptr, extfilename.c_str(), READWRITE, &status);
   fitsReportError(status, "ExposureCubeSun::setCosbinsFieldFormat");
   
   int colnum(1); // by assumption
   fits_modify_vector_len(fptr, colnum, CosineBinner2D::nbins()*CosineBinner2D::nbins2()*(1+CosineBinner2D::nphibins()), &status);
   fitsReportError(status, "ExposureCubeSun::setCosbinsFieldFormat");

   fits_close_file(fptr, &status);
   fitsReportError(status, "ExposureCubeSun::setCosbinsFieldFormat");
}

bool ExposureCubeSun::
acceptInterval(double start, double stop, 
               const std::vector< std::pair<double, double> > & timeCuts,
               const std::vector< std::pair<double, double> > & gtis,
               double & fraction) {
   if (stop < start) {
      std::ostringstream message;
      message << "FT2 files has an interval with START > STOP:\n START = "
              << std::setprecision(15) << start << "\n STOP = "
              << stop;
      throw std::runtime_error(message.str());
   }
   if (start == stop) {
      fraction = 0;
      return false;
   }

   std::pair<double, double> candidateInterval(start, stop);

   typedef std::vector< std::pair<double, double> > IntervalCont_t;
   IntervalCont_t::const_iterator it;

   for (it = timeCuts.begin(); it != timeCuts.end(); ++it) {
      if (!overlaps(*it, candidateInterval)) {
         fraction = 0;
         return false;
      }
   }
   
   double total(0);
   double maxTotal(candidateInterval.second - candidateInterval.first);
   
   IntervalCont_t::const_iterator gti_first = 
      std::upper_bound(gtis.begin(), gtis.end(), candidateInterval,
                       ::compareFirst);
   if (gti_first != gtis.begin()) {
      --gti_first;
   }

   IntervalCont_t::const_iterator gti_last = 
      std::upper_bound(gti_first, gtis.end(), candidateInterval,
                       ::compareSecond);

   if (gti_last != gtis.end()) {
      ++gti_last;
   }

   for (it = gti_first; it != gti_last; ++it) {
      double dt(overlap(*it, candidateInterval));
      total += dt;
   }
   if (total > maxTotal || gtis.size() == 0) {
      total = maxTotal;
   }
   fraction = total/(stop - start);
   return true;
}

bool ExposureCubeSun::
overlaps(const std::pair<double, double> & interval1,
         std::pair<double, double> & interval2) {
   double start = std::max(interval1.first, interval2.first);
   double stop = std::min(interval1.second, interval2.second);
   if (start < stop) {
      interval2.first = start;
      interval2.second = stop;
      return true;
   }
   return false;
}

double ExposureCubeSun::
overlap(const std::pair<double, double> & interval1,
        const std::pair<double, double> & interval2) {
   double start = std::max(interval1.first, interval2.first);
   double stop = std::min(interval1.second, interval2.second);
   if (start < stop) {
      return stop - start;
   }
   return 0;
}

void ExposureCubeSun::
fitsReportError(int status, const std::string & routine) const {
   if (status == 0) {
      return;
   }
   fits_report_error(stderr, status);
   std::ostringstream message;
   message << routine << ": CFITSIO error " << status;
   throw std::runtime_error(message.str());
}

void ExposureCubeSun::
computeBins(std::vector<double> & mubounds, std::vector<double> &thbounds) const {
	CosineBinner2D::cosThetaBins(mubounds);
	thetaBinsSun(thbounds);
}

double ExposureCubeSun::livetime(const astro::SkyDir & dir,
                              double costheta, double thetasun, double phi) const {
   return m_exposure->data()[dir](costheta, thetasun, phi);
}

ExposureCubeSun& ExposureCubeSun::operator += (const ExposureCubeSun &other) {
	*m_exposure += *(other.m_exposure);
	*m_weightedExposure += *(other.m_weightedExposure);
	return *this;
}

void ExposureCubeSun::
thetaBinsSun(std::vector<double> &thSunbounds) const {
	CosineBinner2D::cosTheta2Bins(thSunbounds);
	for (size_t i = 0; i < thSunbounds.size(); ++i) {
		thSunbounds[i] = acos(thSunbounds[i]);
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
