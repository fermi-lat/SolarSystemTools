/**
 * @file LikeExposureSun.cxx
 * @brief Implementation of ExposureSun class for use by the SolarSystemTools.
 * @author G. Johannesson
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/SolarSystemTools/src/LikeExposureSun.cxx,v 1.1.1.1 2012/02/11 02:26:40 gudlaugu Exp $
 */

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "facilities/commonUtilities.h"
#include "facilities/Util.h"

#include "st_stream/StreamFormatter.h"

#include "fitsio.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "SolarSystemTools/CosineBinner2D.h"
#include "healpix/HealpixArray.h"

#include "SolarSystemTools/LikeExposureSun.h"
#include "Likelihood/RoiCuts.h"

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

const double LikeExposureSun::s_mjd_missionStart(astro::JulianDate::missionStart());

LikeExposureSun::
LikeExposureSun(double skybin, double costhetabin, double costhetabinSun,
             const std::vector< std::pair<double, double> > & timeCuts,
             const std::vector< std::pair<double, double> > & gtis,
             double zenmax)
   : ExposureSun(skybin, costhetabin, costhetabinSun, std::cos(zenmax*M_PI/180.)),
     m_costhetabin(costhetabin), m_timeCuts(timeCuts), m_gtis(gtis),
     m_costhetabinSun(costhetabinSun), m_numIntervals(0), 
		 m_solar_dir(astro::SolarSystem::SUN),
     m_weightedExposure(new ExposureSun(skybin, costhetabin, costhetabinSun, 
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
      throw std::runtime_error("LikeExposureSun::LikeExposureSun: GTIs are empty.\n"
                               "Cannot proceed with livetime calculation.");
   }
}

void LikeExposureSun::load(const tip::Table * scData, bool verbose) {
   st_stream::StreamFormatter formatter("LikeExposureSun", "load", 2);
   
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
            fill(astro::SkyDir(ra, dec), scsun, astro::SkyDir(ra_zenith, dec_zenith), 
                 deltat*fraction);
            m_weightedExposure->fill(astro::SkyDir(ra, dec),
								                     scsun,
                                     astro::SkyDir(ra_zenith, dec_zenith), 
                                     deltat*fraction*weight);
         } else {
            fill_zenith(astro::SkyDir(ra, dec), astro::SkyDir(rax, decx),
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

void LikeExposureSun::writeFile(const std::string & outfile) const {
   std::string dataPath = 
      facilities::commonUtilities::getDataPath("SolarSystemTools");
   std::string templateFile = 
      facilities::commonUtilities::joinPath(dataPath, "LivetimeCubeSunTemplate");
   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   fileSvc.createFile(outfile, templateFile);

   writeFilename(outfile);

   writeLivetimes(outfile);

   writeLivetimes(outfile, m_weightedExposure, "WEIGHTED_EXPOSURESUN");

   writeCosbins(outfile);
}

void LikeExposureSun::writeFilename(const std::string & outfile) const {
   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   tip::Image * phdu(fileSvc.editImage(outfile, ""));
   phdu->getHeader()["FILENAME"].set(facilities::Util::basename(outfile));
   delete phdu;
}

void LikeExposureSun::writeLivetimes(const std::string & outfile,
                                  const ExposureSun * self,
                                  const std::string & extname) const {
   if (self == 0) {
      self = this;
   }
	 self->write(outfile,extname);
}
/*
   setCosbinsFieldFormat(outfile, self, extname);

   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   tip::Table * table = fileSvc.editTable(outfile, extname);
   table->setNumRecords(self->data().size());

   tip::Table::Iterator it(table->begin());
   tip::TableRecord & row(*it);

   healpix::HealpixArray<CosineBinner2D>::const_iterator 
      pixel(self->data().begin());

	 std::vector<double> vtmp(CosineBinner2D::nbins()*CosineBinner2D::nthbins()*(1+CosineBinner2D::nphibins()), 0.0);
   for ( ; pixel != self->data().end(); ++pixel, ++it) {
		 std::fill(vtmp.begin(), vtmp.end(), 0);
			for (CosineBinner2D::const_iterator it = (*pixel).begin(); it != (*pixel).end(); ++it) {
				vtmp[(*it).first] = (*it).second;
			}
      row["COSBINS"].set(vtmp);
      astro::SkyDir dir(self->data().dir(pixel));
      row["RA"].set(dir.ra());
      row["DEC"].set(dir.dec());
   }

   tip::Header & header(table->getHeader());
   header["PIXTYPE"].set("HEALPIX"); 
   header["ORDERING"].set("NESTED"); 
   header["COORDSYS"].set(self->data().healpix().galactic()? "GAL" : "EQU");
   header["NSIDE"].set(self->data().healpix().nside()); 
   header["FIRSTPIX"].set(0); 
   header["LASTPIX"].set(self->data().size() - 1); 
   header["THETABIN"].set(CosineBinner2D::thetaBinning());
   header["NBRBINS"].set(CosineBinner2D::nbins()*CosineBinner2D::nthbins());
   header["NBRBINS1"].set(CosineBinner2D::nbins());
   header["COSMIN1"].set(CosineBinner2D::cosmin());
   header["NBRBINS2"].set(CosineBinner2D::nthbins());
   header["COSMIN2"].set(CosineBinner2D::thmax());
   header["PHIBINS"].set(CosineBinner2D::nphibins());

   delete table;
}
	 */

void LikeExposureSun::writeCosbins(const std::string & outfile) const {
   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   tip::Table * table = fileSvc.editTable(outfile, "CTHETABOUNDS");
   table->setNumRecords(CosineBinner2D::nbins());

   tip::Table::Iterator it(table->begin());
   tip::TableRecord & row(*it);
   
   std::vector<double> mubounds, muSunbounds;
   computeCosbins(mubounds,muSunbounds);

   for (size_t i(0); i < mubounds.size() -1; i++, ++it) {
      row["CTHETA_MIN"].set(mubounds.at(i+1));
      row["CTHETA_MAX"].set(mubounds.at(i));
   }
   delete table;

   table = fileSvc.editTable(outfile, "CTHETASUNBOUNDS");
   table->setNumRecords(CosineBinner2D::nthbins());

   tip::Table::Iterator itSun(table->begin());
   tip::TableRecord & rowSun(*itSun);
   
   for (size_t i(0); i < muSunbounds.size() -1; i++, ++itSun) {
      rowSun["CTHETA_MIN"].set(muSunbounds.at(i+1));
      rowSun["CTHETA_MAX"].set(muSunbounds.at(i));
   }
   delete table;
}
   
void LikeExposureSun::setCosbinsFieldFormat(const std::string & outfile,
                                         const ExposureSun * self,
                                         const std::string & extname) const {
   int status(0);

   fitsfile * fptr(0);
   std::string extfilename(outfile + "[" + extname + "]");
   fits_open_file(&fptr, extfilename.c_str(), READWRITE, &status);
   fitsReportError(status, "LikeExposureSun::setCosbinsFieldFormat");
   
   int colnum(1); // by assumption
   fits_modify_vector_len(fptr, colnum, CosineBinner2D::nbins()*CosineBinner2D::nthbins()*(1+CosineBinner2D::nphibins()), &status);
   fitsReportError(status, "LikeExposureSun::setCosbinsFieldFormat");

   fits_close_file(fptr, &status);
   fitsReportError(status, "LikeExposureSun::setCosbinsFieldFormat");
}

bool LikeExposureSun::
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

bool LikeExposureSun::
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

double LikeExposureSun::
overlap(const std::pair<double, double> & interval1,
        const std::pair<double, double> & interval2) {
   double start = std::max(interval1.first, interval2.first);
   double stop = std::min(interval1.second, interval2.second);
   if (start < stop) {
      return stop - start;
   }
   return 0;
}

void LikeExposureSun::
fitsReportError(int status, const std::string & routine) const {
   if (status == 0) {
      return;
   }
   fits_report_error(stderr, status);
   std::ostringstream message;
   message << routine << ": CFITSIO error " << status;
   throw std::runtime_error(message.str());
}

void LikeExposureSun::
computeCosbins(std::vector<double> & mubounds, std::vector<double> &muSunbounds) const {
   bool sqrtbins(CosineBinner2D::thetaBinning() == "SQRT(1-COSTHETA)");
   double cosmin(CosineBinner2D::cosmin());
   size_t nbins(CosineBinner2D::nbins());
   double cosminSun(CosineBinner2D::thmax());
   size_t nbinsSun(CosineBinner2D::nthbins());
   mubounds.clear();
//   for (int i(nbins); i >= 0; i--) {
   for (size_t i(0); i < nbins+1; i++) {
      double factor(static_cast<double>(i)/nbins);
      if (sqrtbins) {
         factor *= factor;
      }
      mubounds.push_back(1. - factor*(1. - cosmin));
   }
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

} // namespace SolarSystemTools
