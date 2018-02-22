/**
 * @file gtexpcubesun.cxx
 * @brief Application for creating binned exposure maps for moving sources
 * @author G. Johannesson
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/SolarSystemTools/src/gtexpcubesun/gtexpcubesun.cxx,v 1.1.1.1 2012/02/11 02:26:40 gudlaugu Exp $
 */

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <memory>
#include <stdexcept>

#include "st_stream/StreamFormatter.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "tip/IFileSvc.h"
#include "tip/Image.h"
#include "tip/Header.h"

#include "dataSubselector/Gti.h"

#include "Likelihood/AppHelpers.h"
#include "SolarSystemTools/BinnedExposureSun.h"
#include "Likelihood/CountsMap.h"
#include "SolarSystemTools/Observation.h"

using namespace SolarSystemTools;

class ExpCubeSun : public st_app::StApp {
public:
   ExpCubeSun();
   virtual ~ExpCubeSun() throw() {
      try {
         delete m_helper;
      } catch (std::exception &eObj) {
         std::cerr << eObj.what() << std::endl;
      } catch (...) {
      }
   }
   virtual void run();
   virtual void banner() const;
private:
	 Likelihood::AppHelpers * m_helper;
   st_app::AppParGroup & m_pars;
   double m_srRadius;
   void promptForParameters();
   void set_phi_status();
   void generateEnergies(std::vector<double> & energies) const;
   void copyGtis() const;
   void copyHeaderKeywords() const;
   void copyDssKeywords(tip::Header & header) const;
   static std::string s_cvs_id;
};

st_app::StAppFactory<ExpCubeSun> myAppFactory("gtexpcubesun");

std::string ExpCubeSun::s_cvs_id("$Name:  $");

ExpCubeSun::ExpCubeSun() : st_app::StApp(), m_helper(0), 
                     m_pars(st_app::StApp::getParGroup("gtexpcubesun")) {
   setVersion(s_cvs_id);
}

void ExpCubeSun::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

void ExpCubeSun::run() {
   promptForParameters();

   std::string cmap_file = m_pars["cmap"];
   std::string ltcube_file = m_pars["infile"];

   bool useEbounds(true);
   std::string bincalc = m_pars["bincalc"];
   if (bincalc == "CENTER") {
      useEbounds = false;
   }

   if (cmap_file != "none") {
      m_helper = new Likelihood::AppHelpers(&m_pars, "BINNED");
   } else {
      m_helper = new Likelihood::AppHelpers(&m_pars, "NONE");
   }
   set_phi_status();
   m_helper->checkOutputFile();
   ExposureCubeSun ltcube; 
	 SolarSystemTools::Observation observation(m_helper->observation(), &ltcube);
	 ltcube.setEfficiencyFactor(observation.respFuncs().efficiencyFactor());
   ltcube.readExposureCubeSun(ltcube_file);

   if (cmap_file != "none") {
// Create map to match counts map.
      m_helper->checkTimeCuts(cmap_file, "", ltcube_file, "Exposure");
			Likelihood::CountsMap cmap(cmap_file);
      BinnedExposureSun bexpmap(cmap, observation, useEbounds, 
                             &m_pars);
      bexpmap.writeOutput(m_pars["outfile"]);
      return;
   }

// Create map for user-defined geometry.
   std::vector<double> energies;
   generateEnergies(energies);
   if (!useEbounds) {
      for (size_t k(0); k < energies.size() - 1; k++) {
         energies[k] = std::sqrt(energies[k]*energies[k+1]);
      }
      energies.pop_back();
   }
   BinnedExposureSun bexpmap(energies, observation, &m_pars);
   bexpmap.writeOutput(m_pars["outfile"]);
   copyGtis();
}

void ExpCubeSun::generateEnergies(std::vector<double> & energies) const {
   energies.clear();
   std::string ebinalg = m_pars["ebinalg"];
   if (ebinalg == "FILE") {
      std::string ebinfile = m_pars["ebinfile"];
      const tip::Table * energybins = 
         tip::IFileSvc::instance().readTable(ebinfile, "ENERGYBINS");
      tip::Table::ConstIterator it = energybins->begin();
      tip::ConstTableRecord & row = *it;
      double energy;
      double emax;
      for ( ; it != energybins->end(); ++it) {
         row["E_MIN"].get(energy);
         // Note that energies in gtbindef output are in units of keV.
         energies.push_back(energy/1e3);
         row["E_MAX"].get(emax);
      }
      energies.push_back(emax/1e3);
      delete energybins;
   } else {
      double emin = m_pars["emin"];
      double emax = m_pars["emax"];
      int enumbins = m_pars["enumbins"];
      double estep = std::log(emax/emin)/enumbins;
      for (size_t k(0); k < enumbins + 1; k++) {
         energies.push_back(emin*std::exp(estep*k));
      }
   }
}

void ExpCubeSun::promptForParameters() {
   m_pars.Prompt("infile");
   m_pars.Prompt("cmap");
   m_pars.Prompt("outfile");
   m_pars.Prompt("irfs");
   std::string cmap = m_pars["cmap"];
   if (cmap == "none") {
      m_pars.Prompt("nxpix");
      m_pars.Prompt("nypix");
      m_pars.Prompt("binsz");
      m_pars.Prompt("coordsys");
      m_pars.Prompt("xref");
      m_pars.Prompt("yref");
      m_pars.Prompt("axisrot");
      m_pars.Prompt("proj");
      std::string ebinalg = m_pars["ebinalg"];
      if (ebinalg == "FILE") {
         m_pars.Prompt("ebinfile");
      } else {
         m_pars.Prompt("emin");
         m_pars.Prompt("emax");
         m_pars.Prompt("enumbins");
      }
   }
   m_pars.Save();
}

void ExpCubeSun::set_phi_status() {
   // If indicated, turn off phi-dependence for all IRFs.
   bool ignorephi = m_pars["ignorephi"];
   if (ignorephi) {
      const Likelihood::Observation & observation(m_helper->observation());
      std::map<unsigned int, irfInterface::Irfs *>::const_iterator respIt 
         = observation.respFuncs().begin();
      for ( ; respIt != observation.respFuncs().end(); ++respIt) {
         respIt->second->aeff()->setPhiDependence(false);
      }
   }
}

void ExpCubeSun::copyGtis() const {
   std::string infile = m_pars["infile"];
   dataSubselector::Gti gti(infile);
   std::string outfile = m_pars["outfile"];
   gti.writeExtension(outfile);
}

void ExpCubeSun::copyHeaderKeywords() {
   std::string infile = m_pars["infile"];
   std::string inext("EXPOSURESUN");
   const tip::Table * intab = 
     tip::IFileSvc::instance().readTable(infile, inext);
   const tip::Header & inheader(intab->getHeader());
   
   std::string outfile = m_pars["outfile"];
   std::string outext("EXPOSURE");
   tip::Table * outtab = 
     tip::IFileSvc::instance().readTable(outfile, outext);
   tip::Header & outheader(outtab->getHeader());

   // Unfortunately TIP does not provide access to header comments, so
   // we cannot copy them using it. Could use cfitsio directly but
   // prefer TIP as it is standard here.

#define COPYKEYWORD(type, name)			\
   try {					\
     type x;					\
     inheader[name].get(x);			\
     outheader[name].set(x);			\
   } catch(...) {				\
   }
   
   COPYKEYWORD(std::string, "DATE-OBS");
   COPYKEYWORD(std::string, "DATE-END");
   COPYKEYWORD(double,      "TSTART");
   COPYKEYWORD(double,      "TSTOP");
   COPYKEYWORD(double,      "MJDREFI");
   COPYKEYWORD(double,      "MJDREFF");
   COPYKEYWORD(std::string, "TIMEUNIT");
   COPYKEYWORD(double,      "TIMEZERO");
   COPYKEYWORD(std::string, "TIMESYS");
   COPYKEYWORD(std::string, "TIMEREF");
   COPYKEYWORD(bool,        "CLOCKAPP");
   COPYKEYWORD(bool,        "GPS_OUT");

#undef COPYKEYWORD

   copyDssKeywords(outheader);

   delete outtab;
   delete intab;
}

void ExpCubeSun::copyDssKeywords(tip::Header & header) const {
   dataSubselector::Cuts * irfs_cuts(0);
   irfs_cuts = new dataSubselector::Cuts();
   irfs_cuts->addVersionCut("IRF_VERSION", m_helper->irfsName());
   irfs_cuts->writeDssKeywords(header);
   delete irfs_cuts;
}
