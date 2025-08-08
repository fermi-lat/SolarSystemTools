/** 
 * @file makeSolarExposureCubeSun.cxx
 * @brief Create an Exposure hypercube including distance from solar center
 * @author G. Johannesson
 *
 *  $Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/src/makeSolarExposureCube/makeSolarExposureCube.cxx,v 1.7 2012/08/24 13:48:48 gudlaugu Exp $
 */

#include <cstdlib>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <memory>
#include <stdexcept>
#include <cstdio>

#include "st_stream/StreamFormatter.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "st_facilities/Util.h"

#include "SolarSystemTools/CosineBinner2D.h"

#include "SolarSystemTools/ExposureCubeSun.h"
#include "Likelihood/RoiCuts.h"
#include "general_util/generic_utils.h"

namespace {
   void getTBounds(const std::vector< std::pair<double, double> > & gtis,
                   double & tmin, double & tmax) {
      if (gtis.size() > 0) {
         tmin = gtis.at(0).first;
         tmax = gtis.at(0).second;
      }
      for (size_t i(0); i < gtis.size(); i++) {
         if (gtis.at(i).first < tmin) {
            tmin = gtis.at(i).first;
         }
         if (gtis.at(i).second > tmax) {
            tmax = gtis.at(i).second;
         }
      }         
   }
   void getTimeBounds(const std::vector< std::pair<double, double> > & gtis,
                      const std::vector< std::pair<double, double> > & tranges,
                      double & tmin, double & tmax) {
      getTBounds(gtis, tmin, tmax);
      double t0(tmin), t1(tmax);
      getTBounds(tranges, t0, t1);
      if (t0 < tmin) {
         tmin = t0;
      }
      if (t1 > tmax) {
         tmax = t1;
      }
   }
}

/**
 * @class ExposureCubeSun
 * @brief Class to encapsulate methods for creating an exposure
 * hypercube in (ra, dec, cos_theta, cos_theta_sun) using the ExposureCubeSun class.
 *
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/src/makeSolarExposureCube/makeSolarExposureCube.cxx,v 1.7 2012/08/24 13:48:48 gudlaugu Exp $
 */
class ExposureCubeSun : public st_app::StApp {
public:
   ExposureCubeSun() : st_app::StApp(), 
                    m_pars(st_app::StApp::getParGroup("gtltcubesun")), 
                    m_exposure(0), m_roiCuts(0) {
      setVersion(s_cvs_id);
   }
   virtual ~ExposureCubeSun() throw() {
      try {
         delete m_exposure;
         delete m_roiCuts;
      } catch (std::exception &eObj) {
         std::cerr << eObj.what() << std::endl;
      } catch (...) {
      }
    }
   virtual void run();
   virtual void banner() const;
private:
   st_app::AppParGroup & m_pars;
   SolarSystemTools::ExposureCubeSun * m_exposure;
   Likelihood::RoiCuts * m_roiCuts;
   void readRoiCuts();
   void createDataCube();
   void writeTableKeywords(const std::string & outfile,
                           const std::string & tablename) const;
   void writeDateKeywords(const std::string & outfile, 
                          double tstart, double tstop) const;
   static std::string s_cvs_id;
};

st_app::StAppFactory<ExposureCubeSun> myAppFactory("gtltcubesun");

std::string ExposureCubeSun::s_cvs_id("$Name:  $");

void ExposureCubeSun::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

void ExposureCubeSun::run() {
   m_pars.Prompt();
   m_pars.Save();
   readRoiCuts();
   std::string output_file = m_pars["outfile"];
   if (st_facilities::Util::fileExists(output_file)) {
      if (m_pars["clobber"]) {
         std::remove(output_file.c_str());
      } else {
         st_stream::StreamFormatter formatter("gtltcubesun", "run", 2);
         formatter.err() << "Output file " << output_file 
                         << " already exists,\n"
                         << "and you have set 'clobber' to 'no'.\n"
                         << "Please provide a different output file name." 
                         << std::endl;
         std::exit(1);
      }
   }
   createDataCube();

   const double tstart(m_roiCuts->minTime());
   const double tstop(m_roiCuts->maxTime());

   m_exposure->writeFile(output_file, tstart, tstop, *m_roiCuts);

	 writeDateKeywords(output_file, tstart, tstop);
}

void ExposureCubeSun::writeTableKeywords(const std::string & outfile,
                                      const std::string & tablename) const {
   std::unique_ptr<tip::Table> 
      table(tip::IFileSvc::instance().editTable(outfile, tablename));
   m_roiCuts->writeDssTimeKeywords(table->getHeader());
   double tstart(m_roiCuts->minTime());
   double tstop(m_roiCuts->maxTime());
   tip::Header & header(table->getHeader());
   header["TSTART"].set(tstart);
   header["TSTOP"].set(tstop);
}

void ExposureCubeSun::writeDateKeywords(const std::string & outfile, 
                                     double tstart, double tstop) const {
   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   std::vector<std::string> extnames;
   extnames.push_back("");
   extnames.push_back("EXPOSURESUN");
   extnames.push_back("WEIGHTED_EXPOSURESUN");
   extnames.push_back("CTHETABOUNDS");
   extnames.push_back("THETASUNBOUNDS");
   extnames.push_back("GTI");
   for (std::vector<std::string>::const_iterator name(extnames.begin());
        name != extnames.end(); ++name) {
      tip::Extension * hdu(fileSvc.editExtension(outfile, *name));
      st_facilities::Util::writeDateKeywords(hdu, tstart, tstop, *name!="");
      if (*name == "") {
         // Update CREATOR value with: Tool name/Group/Version
         std::string creator_version =  GenericUtils::creator_banner("gtltcubesun");
         hdu->getHeader()["CREATOR"].set(creator_version);
         std::string file_version = m_pars["file_version"];
         hdu->getHeader()["VERSION"].set(file_version);
      }
      delete hdu;
   }
}

void ExposureCubeSun::readRoiCuts() {
   std::string event_file = m_pars["evfile"];
   m_roiCuts = new Likelihood::RoiCuts();
   if (event_file != "") {
      std::string evtable = m_pars["evtable"];
      std::vector<std::string> eventFiles;
      st_facilities::Util::resolve_fits_files(event_file, eventFiles);
      m_roiCuts->readCuts(eventFiles, evtable, false);
   }
	 //Override the maximum and minimum time, convenince for splitting time
	 //ranges into smaller with a monolithic ft1 file.
   double tmin = m_pars["tmin"];
   double tmax = m_pars["tmax"];
	 if ( tmin == 0 ) tmin = m_roiCuts->minTime();
	 if ( tmax == 0 ) tmax = m_roiCuts->maxTime();
   m_roiCuts->setCuts(0, 0, 180, 30, 3e5, tmin, tmax, -1, true);
}

void ExposureCubeSun::createDataCube() {
   st_stream::StreamFormatter formatter("gtltcubesun", 
                                        "createDataCube", 2);

   std::vector<std::pair<double, double> > timeCuts;
   std::vector<std::pair<double, double> > gtis;
   m_roiCuts->getTimeCuts(timeCuts);
   m_roiCuts->getGtis(gtis);

   double tmin, tmax;
   ::getTimeBounds(gtis, timeCuts, tmin, tmax);
   static double maxIntervalSize(30);
   tmin -= 2.*maxIntervalSize;
   tmax += 2.*maxIntervalSize;
   std::ostringstream filter;
   filter << std::setprecision(20);
   filter << "(START >= " << tmin << ") && (STOP <= " << tmax << ")";
   formatter.info(4) << "applying filter: " << filter.str() << std::endl;

   double zmax = m_pars["zmax"];
   if (zmax < 180.) {
      formatter.info(2) << "WARNING: You have chosen to apply a zenith angle cut of "
                        << zmax << " degrees." << std::endl
                        << "Applying such a cut for this tool is not equivalent to \n"
                        << "applying a zenith angle cut in gtselect." << std::endl
                        << "If you don't understand this comment, " << std::endl
                        << "then you probably shouldn't be applying this cut." 
                        << std::endl;
   }

   // Set the number of phibins using the static function interface
   // from healpix::CosineBinner (this is how
   // map_tools/exposure_cube.cxx does it.)
   long nphibins = m_pars["phibins"];
   if (nphibins > 0) {
      SolarSystemTools::CosineBinner2D::setPhiBins(nphibins);
   }

	 // Set the SolarSystem body
	 astro::SolarSystem::Body body = SolarSystemTools::ExposureCubeSun::stringToBody(m_pars["body"]);

	 // The first bin should always be 0.25 degrees to account for the sun and
	 // the moon.  We have to slightly adjust the powerbinsun parameter
	 const double thsunmax = m_pars["thetasunmax"];
	 const double costhsunmax = cos(thsunmax*M_PI/180.);
	 double powerbin = m_pars["powerbinsun"];
	 const double ratio = (1-costhsunmax)/(1-cos(0.25*M_PI/180.));
	 unsigned int cosbins2 = static_cast<unsigned int>( pow( ratio, 1./powerbin) );
	 // Need at least 2 bins for this to work
	 if (cosbins2 < 2)
		 ++cosbins2;
	 powerbin = log(ratio)/log(static_cast<double>(cosbins2)) + 1e-10; //The small number is to make sure the rounding does not go the wrong way

   m_exposure = new SolarSystemTools::ExposureCubeSun(m_pars["binsz"], 
                                             m_pars["dcostheta"],
																						 0.25,
                                             thsunmax,
																						 powerbin,
                                             timeCuts, gtis, body, zmax);
   std::string scFile = m_pars["scfile"];
   st_facilities::Util::file_ok(scFile);
   std::vector<std::string> scFiles;
   st_facilities::Util::resolve_fits_files(scFile, scFiles);
   std::vector<std::string>::const_iterator scIt = scFiles.begin();
   for ( ; scIt != scFiles.end(); scIt++) {
      st_facilities::Util::file_ok(*scIt);
      formatter.err() << "Working on file " << *scIt << std::endl;
      const tip::Table * scData = 
         tip::IFileSvc::instance().readTable(*scIt, m_pars["sctable"],
                                             filter.str());
      formatter.info(4) << "read " << scData->getNumRecords() 
                        << " rows" << std::endl;
      int chatter = m_pars["chatter"];
      bool print_output(true);
      if (chatter < 2) {
         print_output = false;
      }
      m_exposure->load(scData, print_output);
      delete scData;
   }

   if (m_exposure->numIntervals() == 0) {
      formatter.warn() << "WARNING: No intervals have been read in from "
                       << "the FT2 files that correspond to the FT1 data.\n"
                       << "All livetimes will be identically zero."
                       << std::endl;
   }
}
