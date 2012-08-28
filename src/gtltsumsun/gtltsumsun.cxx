/**
 * @file gtltsumsun.cxx
 * @brief Application for creating binned exposure maps for moving sources
 * @author G. Johannesson
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/src/gtltsumsun/gtltsumsun.cxx,v 1.1 2012/06/16 11:29:15 gudlaugu Exp $
 */

#include <cstdlib>
#include <string>
#include <cstdio>

#include <memory>
#include <stdexcept>

#include "st_facilities/Util.h"

#include "st_stream/StreamFormatter.h"

#include "st_app/AppParGroup.h"
#include "st_app/StApp.h"
#include "st_app/StAppFactory.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"
#include "tip/Header.h"

#include "dataSubselector/Cuts.h"
#include "Likelihood/AppHelpers.h"

#include "SolarSystemTools/ExposureCubeSun.h"

using namespace SolarSystemTools;

class SumSunLivetime : public st_app::StApp {
public:
   SumSunLivetime();
   virtual ~SumSunLivetime() throw() {
      try {
      } catch (std::exception &eObj) {
         std::cerr << eObj.what() << std::endl;
      } catch (...) {
      }
   }
   virtual void run();
   virtual void banner() const;
private:
   st_app::AppParGroup & m_pars;
   std::vector<std::string> m_fileList;
	 static std::string s_cvs_id;

   void promptForParameters();
};

st_app::StAppFactory<SumSunLivetime> myAppFactory("gtltsumsun");

std::string SumSunLivetime::s_cvs_id("$Name:  $");

SumSunLivetime::SumSunLivetime() : st_app::StApp(),  
                     m_pars(st_app::StApp::getParGroup("gtltsumsun")) {
   setVersion(s_cvs_id);
}

void SumSunLivetime::banner() const {
   int verbosity = m_pars["chatter"];
   if (verbosity > 2) {
      st_app::StApp::banner();
   }
}

void SumSunLivetime::run() {
   promptForParameters();

	 ExposureCubeSun ltcube;
	 ltcube.readExposureCubeSun(m_fileList.front());
	 
   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());

	 /*?
	 std::vector<dataSubselector::Cuts> my_cuts;
	 my_cuts.push_back(dataSubselector::Cuts(m_fileList.front(), "EXPOSURESUN",
				 false, true));
				 */

   const tip::Table * table(fileSvc.readTable(m_fileList.front(), "EXPOSURESUN"));
   double tstart, tstop;
	 table->getHeader()["TSTART"].get(tstart);
	 table->getHeader()["TSTOP"].get(tstop);
	 delete table;


	 //Store the number of bins for testing since it is not checked elsewhere
	 const size_t nbins = CosineBinner2D::nbins();
	 const double cosmin = CosineBinner2D::cosmin();
	 const size_t nbins2 = CosineBinner2D::nbins2();
	 const double cosmin2 = CosineBinner2D::cosmin2();
	 const double power2 = CosineBinner2D::power2();
	 const size_t nphibins = CosineBinner2D::nphibins();

   st_stream::StreamFormatter formatter("gtltsumsun", "run", 2);
   formatter.warn() << "Adding livetime cubes";
   for (size_t k=1; k < m_fileList.size(); ++k) {
		 formatter.warn() << ".";
		 ExposureCubeSun other;
		 other.readExposureCubeSun(m_fileList[k]);

		 if (nbins != CosineBinner2D::nbins() 
				 || nbins2 != CosineBinner2D::nbins2() 
				 || cosmin != CosineBinner2D::cosmin() 
				 || cosmin2 != CosineBinner2D::cosmin2() 
				 || power2 != CosineBinner2D::power2() 
				 || nphibins != CosineBinner2D::nphibins() ){
			 std::ostringstream message;
			 message << "The binning of the Exposure extension in " 
				 << m_fileList.at(k)
				 << " does not match the binning of the extension in  " 
				 << m_fileList.front()
				 << ": " << nbins <<" != "<< CosineBinner2D::nbins()
				 << ": " << nbins2 <<" != "<< CosineBinner2D::nbins2()
				 << ": " << cosmin <<" != "<< CosineBinner2D::cosmin()
				 << ": " << cosmin2 <<" != "<< CosineBinner2D::cosmin2()
				 << ": " << power2 <<" != "<< CosineBinner2D::power2()
				 << ": " << nphibins <<" != "<< CosineBinner2D::nphibins();
			 throw std::runtime_error(message.str());
		 }

		 try {
			 ltcube += other;
		 } catch (std::runtime_error &e) {
			 std::ostringstream message;
			 message << "Error while adding \""
				 << m_fileList.at(k)
				 << "\": "
				 << e.what();
			 throw std::runtime_error(message.str());
		 }

		 const double start = other.tstart();
		 const double stop = other.tstop();
		 if (start < tstart) {
			 tstart = start;
		 }
		 if (stop > tstop) {
			 tstop = stop;
		 }
		 //my_cuts.push_back(dataSubselector::Cuts(m_fileList.at(k), "EXPOSURESUN",
			//		 false, true));
	}

	 std::string outfile = m_pars["outfile"];
	 std::remove(outfile.c_str());

	 formatter.warn() << " Done" << std::endl << "Merging cuts";
	 Likelihood::RoiCuts cuts;
	 cuts.readCuts(m_fileList, "EXPOSURESUN", false);

	 formatter.warn() << " Done" << std::endl << "Writing file";
	 ltcube.writeFile(outfile, tstart, tstop, cuts);
	 formatter.warn() << " Done" << std::endl;

}


void SumSunLivetime::promptForParameters() {
   m_pars.Prompt("infile1");
   std::string infile1 = m_pars["infile1"];
   st_facilities::Util::resolve_fits_files(infile1, m_fileList);
   if (m_fileList.size() == 1) {
      m_pars.Prompt("infile2");
      std::string infile2 = m_pars["infile2"];
      m_fileList.push_back(infile2);
   }
   m_pars.Prompt("outfile");
   m_pars.Save();
}

