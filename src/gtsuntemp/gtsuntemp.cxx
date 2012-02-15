/**
 * @file gtsuntemp.cxx
 * @brief Application for creating binned exposure maps for moving sources
 * @author G. Johannesson
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/SolarSystemTools/src/gtsuntemp/gtsuntemp.cxx,v 1.1.1.1 2012/02/11 02:26:40 gudlaugu Exp $
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
#include "SolarSystemTools/HealpixExposureSun.h"
#include "SolarSystemTools/SolarProfile.h"
#include "SolarSystemTools/FitsSolarProfile.h"
#include "SolarSystemTools/SolarTemplate.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/BinnedExposure.h"

using namespace SolarSystemTools;

class SunTemplate : public st_app::StApp {
public:
   SunTemplate();
   virtual ~SunTemplate() throw() {
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
   void generateEnergies(std::vector<double> & energies) const;
   void copyGtis() const;
   static std::string s_cvs_id;
};

st_app::StAppFactory<SunTemplate> myAppFactory("gtsuntemp");

std::string SunTemplate::s_cvs_id("$Name:  $");

SunTemplate::SunTemplate() : st_app::StApp(), m_helper(0), 
                     m_pars(st_app::StApp::getParGroup("gtsuntemp")) {
   setVersion(s_cvs_id);
}

void SunTemplate::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

void SunTemplate::run() {
   promptForParameters();

   std::string cmap_file = m_pars["cmap"];
   std::string expsun_file = m_pars["expsun"];
   std::string avgexp_file = m_pars["avgexp"];
   std::string solarprofile_file = m_pars["sunprof"];

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
   m_helper->checkOutputFile();
   HealpixExposureSun expsun(expsun_file); 
	 Likelihood::BinnedExposure avgexp(avgexp_file);
	 FitsSolarProfile prof(solarprofile_file);

   if (cmap_file != "none") {
// Create map to match counts map.
			Likelihood::CountsMap cmap(cmap_file);
      SolarTemplate temp(cmap, expsun, avgexp, prof, useEbounds, 
                             &m_pars);
      temp.writeOutput(m_pars["outfile"]);
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
   SolarTemplate temp(energies, expsun, avgexp, prof, &m_pars);
   temp.writeOutput(m_pars["outfile"]);
   copyGtis();
}

void SunTemplate::generateEnergies(std::vector<double> & energies) const {
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

void SunTemplate::promptForParameters() {
   m_pars.Prompt("expsun");
   m_pars.Prompt("avgexp");
	 m_pars.Prompt("sunprof");
   m_pars.Prompt("cmap");
   m_pars.Prompt("outfile");
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

void SunTemplate::copyGtis() const {
   std::string infile = m_pars["expsun"];
   dataSubselector::Gti gti(infile);
   std::string outfile = m_pars["outfile"];
   gti.writeExtension(outfile);
}
