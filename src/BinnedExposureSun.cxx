/**
 * @file BinnedExposureSun.cxx
 * @brief Integral of effective area over time for the entire sky at
 * various energies and binned in angles from the sun.
 * @author G. Johannesson
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/src/BinnedExposureSun.cxx,v 1.3 2012/02/15 03:03:51 gudlaugu Exp $
 */

#include <cmath>
#include <cstdio>

#include <algorithm>
#include <memory>
#include <sstream>
#include <stdexcept>

#include "st_stream/StreamFormatter.h"

#include "st_app/AppParGroup.h"

#include "tip/Header.h"
#include "tip/IFileSvc.h"
#include "tip/Image.h"

#include "astro/SkyProj.h"

#include "st_facilities/Util.h"

#include "SolarSystemTools/BinnedExposureSun.h"
#include "Likelihood/CountsMap.h"
#include "SolarSystemTools/Observation.h"

namespace {
   inline double fracDiff(double target, double result) {
      return std::fabs((target - result)/target);
   }
   std::vector<double>::const_iterator 
   findNearest(const std::vector<double> & xx, double x, double tol=1e-5) {
      std::vector<double>::const_iterator ix = std::find(xx.begin(),
                                                         xx.end(), x);
      if (ix == xx.end()) { // no exact match, so look for nearest
         for (ix = xx.begin(); ix != xx.end(); ++ix) {
            if (fracDiff(x, *ix) < tol) {
               return ix;
            }
         }
         std::ostringstream what;
         what << "BinnedExposureSun::operator(): The energy " << x
              << " is not available.\nHere are the relevant energies:\n";
         for (size_t i(0); i < xx.size(); i++) {
            what << xx.at(i) << "\n";
         }
         throw std::runtime_error(what.str());
      }
      return ix;  // return the exact match
   }
}

namespace SolarSystemTools {

BinnedExposureSun::BinnedExposureSun() : m_observation(0), m_proj(0), 
                                   m_costhmin(-1), m_costhmax(1),
                                   m_enforce_boundaries(false) {}

BinnedExposureSun::BinnedExposureSun(const Likelihood::CountsMap & cmap,
                               const Observation & observation,
                               bool useEbounds,
                               const st_app::AppParGroup * pars)
   : m_observation(&observation), m_proj(0), m_costhmin(-1), m_costhmax(1),
     m_enforce_boundaries(false) {
   setMapGeometry(cmap);
   if (!useEbounds) {
      for (size_t k(0); k < m_energies.size()-1; k++) {
         m_energies[k] = std::sqrt(m_energies[k]*m_energies[k+1]);
      }
      m_energies.pop_back();
      m_naxes[2] -= 1;
   }
   if (pars) {
      setCosThetaBounds(*pars);
   }
   computeMap();
}

BinnedExposureSun::BinnedExposureSun(const std::vector<double> & energies,
                               const Observation & observation,
                               const st_app::AppParGroup * pars) 
   : m_energies(energies), m_observation(&observation), m_proj(0),
     m_costhmin(-1), m_costhmax(1),m_enforce_boundaries(false)  {
   if (pars) {
      setMapGeometry(*pars);
      setCosThetaBounds(*pars);
   } else {
      setMapGeometry();
   }
   computeMap();
}

BinnedExposureSun::BinnedExposureSun(const std::string & filename) 
   : m_observation(0), m_proj(0), m_costhmin(-1), m_costhmax(1),
     m_enforce_boundaries(false) {
   m_proj = new astro::SkyProj(filename);

   std::auto_ptr<const tip::Image> 
      image(tip::IFileSvc::instance().readImage(filename, ""));

   m_exposureMap.clear();
   image->get(m_exposureMap);

   m_naxes = image->getImageDimensions();

   std::auto_ptr<const tip::Table>
      energies(tip::IFileSvc::instance().readTable(filename, "Energies"));

   m_energies.clear();
   tip::Table::ConstIterator it = energies->begin();
   tip::ConstTableRecord & row = *it;
   for ( ; it != energies->end(); ++it) {
      double value;
      row["Energy"].get(value);
      m_energies.push_back(value);
   }
}

BinnedExposureSun::~BinnedExposureSun() {
   delete m_proj;
}

double BinnedExposureSun::operator()(double energy, double ra, double dec, double theta) const {
   std::vector<double>::const_iterator ie = ::findNearest(m_energies, energy);
   unsigned int k = ie - m_energies.begin();

   int l = 0;
   for ( ; l < m_naxes[3]; ++l){
     if ( m_thetasun[l] < theta && theta < m_thetasun[l+1] )
       break;
   }

   std::pair<double, double> pixel;
   st_facilities::Util::skyDir2pixel(*m_proj, astro::SkyDir(ra, dec),
                                     pixel.first, pixel.second);

   int i = static_cast<int>(pixel.first - 1);
   int j = static_cast<int>(pixel.second - 1);

	 const unsigned int indx = ((l*m_energies.size() + k)*m_naxes.at(1) + j)*m_naxes.at(0) + i;

   bool within_bounds = (i >= 0 && i < m_naxes[0] &&
                         j >= 0 && j < m_naxes[1] &&
                         k >= 0 && k < m_energies.size() &&
                         l >= 0 && l < m_naxes[3]);

   if (m_enforce_boundaries && !within_bounds) {
      throw std::runtime_error("Request for exposure at a sky position that "
                               "is outside of the map boundaries.");
   }

   try {
      return m_exposureMap.at(indx);
   } catch (std::out_of_range &) {
      // Range check performed already, so do nothing and return 0.
   }
   return 0;
}

void BinnedExposureSun::setMapGeometry(const Likelihood::CountsMap & cmap) {
   m_proj_name = cmap.proj_name();
   m_crpix[0] = cmap.crpix1();
   m_crpix[1] = cmap.crpix2();
   m_crval[0] = cmap.crval1();
   m_crval[1] = cmap.crval2();
   m_cdelt[0] = cmap.cdelt1();
   m_cdelt[1] = cmap.cdelt2();
   m_crota2 = cmap.crota2();
   m_naxes.resize(4, 0);
   m_naxes[0] = cmap.naxis1();
   m_naxes[1] = cmap.naxis2();
   m_naxes[2] = cmap.energies().size();
   m_naxes[3] = m_observation->expCubeSun().nthetaBinsSun();
   m_energies = cmap.energies();
   m_observation->expCubeSun().thetaBinsSun(m_thetasun);
   m_isGalactic = cmap.isGalactic();
}

void BinnedExposureSun::setMapGeometry(const st_app::AppParGroup & pars) {
   m_naxes.resize(4, 0);
   m_naxes[0] = pars["nxpix"];
   m_naxes[1] = pars["nypix"];
   m_naxes[2] = m_energies.size();
   m_naxes[3] = m_observation->expCubeSun().nthetaBinsSun();
   m_observation->expCubeSun().thetaBinsSun(m_thetasun);
   double binsz = pars["binsz"];
   m_cdelt[0] = -binsz;
   m_cdelt[1] = binsz;
   m_crval[0] = pars["xref"];
   m_crval[1] = pars["yref"];
   m_crota2 = pars["axisrot"];
   std::string proj_name = pars["proj"];
   m_proj_name = proj_name;
   m_crpix[0] = m_naxes[0]/2. + 0.5;
   m_crpix[1] = m_naxes[1]/2. + 0.5;
   std::string coordsys = pars["coordsys"];
   m_isGalactic = (coordsys == "GAL");
}

void BinnedExposureSun::setMapGeometry() {
   m_proj_name = "CAR";
   m_isGalactic = false;
   m_naxes.resize(4, 0);
   m_naxes[0] = 360;
   m_naxes[1] = 180;
   m_naxes[2] = m_energies.size();
   m_naxes[3] = m_observation->expCubeSun().nthetaBinsSun();
   m_observation->expCubeSun().thetaBinsSun(m_thetasun);
   m_crpix[0] = m_naxes[0]/2. + 0.5;
   m_crpix[1] = m_naxes[1]/2. + 0.5;
   m_crval[0] = 180.;
   m_crval[1] = 0;
   m_cdelt[0] = -1;
   m_cdelt[1] = 1;
   m_crota2 = 0;
}

void BinnedExposureSun::computeMap() {
   m_proj = new astro::SkyProj(m_proj_name, &m_crpix[0], &m_crval[0],
                               &m_cdelt[0], m_crota2, m_isGalactic);

   m_exposureMap.resize(m_naxes.at(0)*m_naxes.at(1)*m_energies.size()*m_naxes.at(3), 0);
   int iter(0);
   st_stream::StreamFormatter formatter("BinnedExposureSun", "computeMap", 2);
   formatter.warn() << "Computing binned exposure map";

   std::vector<double> thetasun(m_naxes.at(3));
   for (int j = 0; j < m_naxes.at(3); ++j) {
     thetasun[j] = (m_thetasun[j]+m_thetasun[j+1])/2.;
   }

	 //Create a cache for AEff calculations
	 std::vector<Aeff*> aeffs;
	 for (size_t i = 0; i < m_energies.size(); ++i)
     aeffs.push_back(new Aeff(m_energies[i], *m_observation, m_costhmin, m_costhmax));

   for (int j = 0; j < m_naxes.at(1); j++) {
      for (int i = 0; i < m_naxes.at(0); i++, iter++) {
         if ((iter % ((m_naxes.at(1)*m_naxes.at(0))/20)) == 0) {
            formatter.warn() << ".";
         }
         // std::pair<double, double> coord = m_proj->pix2sph(i + 1, j + 1);
         // astro::SkyDir dir(coord.first, coord.second, coordsys);
         astro::SkyDir dir;
         try {
            st_facilities::Util::pixel2SkyDir(*m_proj, i + 1, j + 1, dir);
         } catch (...) {
            // The astro::SkyProj class throws a SkyProjException
            // here, but SkyProjException is annoyingly defined in
            // SkyProj.cxx
            // http://www-glast.stanford.edu/cgi-bin/viewcvs/astro/src/SkyProj.cxx?revision=1.27&view=markup
            // so that client code cannot catch it directly. Amazing.
            continue;
         }
                                           
         for (unsigned int k = 0; k < m_energies.size(); k++) {
               const Aeff &aeff = *aeffs[k]; 
               for (int l = 0; l < m_naxes.at(3); ++l) {
								   const unsigned int indx = ((l*m_energies.size() + k)*m_naxes.at(1) + j)*m_naxes.at(0) + i;
                   m_exposureMap.at(indx)
                     +=m_observation->expCubeSun().value(dir, thetasun[l], aeff, m_energies.at(k));
               }
         }
      }
   }
   formatter.warn() << "!" << std::endl;
	 for (size_t i = 0; i < m_energies.size(); ++i)
     delete aeffs[i];
}

void BinnedExposureSun::writeOutput(const std::string & filename) const {
   std::remove(filename.c_str());

   std::string ext("PRIMARY");
   tip::IFileSvc::instance().appendImage(filename, ext, m_naxes);
   tip::Image * image = tip::IFileSvc::instance().editImage(filename, ext);

   image->set(m_exposureMap);

   tip::Header & header(image->getHeader());

   header["TELESCOP"].set("GLAST");
   header["INSTRUME"].set("LAT");
   header["DATE-OBS"].set("");
   header["DATE-END"].set("");

   header["CRVAL1"].set(m_crval[0]);
   double crpix1 = static_cast<double>(m_naxes[0] + 1)/2.;
   header["CRPIX1"].set(crpix1);
   header["CDELT1"].set(m_cdelt[0]);

   header["CRVAL2"].set(m_crval[1]);
   double crpix2 = static_cast<double>(m_naxes[1] + 1)/2.;
   header["CRPIX2"].set(crpix2);
   header["CDELT2"].set(m_cdelt[1]);

   if (m_isGalactic) {
      header["CTYPE1"].set("GLON-" + m_proj_name);
      header["CTYPE2"].set("GLAT-" + m_proj_name);
   } else {
      header["CTYPE1"].set("RA---" + m_proj_name);
      header["CTYPE2"].set("DEC--" + m_proj_name);
   }

   int nee = m_energies.size();
   header["CRVAL3"].set(log(m_energies.at(0)));
   header["CRPIX3"].set(1);
   if (nee == 1) {
      header["CDELT3"].set(0);
   } else {
      header["CDELT3"].set(log(m_energies.at(nee-1)/m_energies.at(0))/(nee-1));
   }
   header["CTYPE3"].set("log_Energy");

   header["CRVAL4"].set(m_thetasun[0]);
   header["CRPIX4"].set(1);
   if (m_naxes[3] == 1) {
     header["CDELT4"].set(0);
   } else {
     header["CDELT4"].set(m_thetasun[1]-m_thetasun[0]);
   }

   delete image;

   ext = "ENERGIES";
   tip::IFileSvc::instance().appendTable(filename, ext);
   tip::Table * table = tip::IFileSvc::instance().editTable(filename, ext);
   table->appendField("Energy", "1D");
   table->setNumRecords(m_energies.size());

   tip::Table::Iterator row = table->begin();
   tip::Table::Record & record = *row;

   std::vector<double>::const_iterator energy = m_energies.begin();
   for ( ; energy != m_energies.end(); ++energy, ++row) {
      record["Energy"].set(*energy);
   }

   delete table;

   ext = "COSTHETASUN";
   tip::IFileSvc::instance().appendTable(filename, ext);
   table = tip::IFileSvc::instance().editTable(filename, ext);
   table->appendField("CosTheta_Min", "1D");
   table->appendField("CosTheta_Max", "1D");
   table->setNumRecords(m_naxes[3]);

   row = table->begin();
   tip::Table::Record & record2 = *row;

   for (int i = 0; i < m_naxes[3]; ++i, ++row) {
      record2["CosTheta_Min"].set(m_thetasun[i]);
      record2["CosTheta_Max"].set(m_thetasun[i+1]);
   }

   delete table;
}

void BinnedExposureSun::setCosThetaBounds(const st_app::AppParGroup & pars) {
   double thmin = pars["thmin"];
   if (thmin > 0) {
      m_costhmax = std::cos(thmin*M_PI/180.);
   }
   double thmax = pars["thmax"];
   if (thmax < 180.) {
      m_costhmin = std::cos(thmax*M_PI/180.);
   }
}

double BinnedExposureSun::Aeff::operator()(double cosTheta, double phi) const {
   if (cosTheta < m_costhmin || cosTheta > m_costhmax) {
      return 0;
   }
	 //Check the cache for the value, add it if it does not exist
	 const std::pair<double,double> indx(std::make_pair(cosTheta,phi));
	 std::map<std::pair<double,double>,double>::iterator it = m_cache.find(indx);
	 if (it == m_cache.end()) {
		 it = m_cache.insert(std::make_pair(indx,ExposureCubeSun::Aeff::operator()(cosTheta, phi))).first;
	 }
	 //if ( (*it).second != ExposureCubeSun::Aeff::operator()(cosTheta,phi))
	//	 std::cout<<"Cache is not working properly: "<<cosTheta<<", "<<phi<<std::endl;
   return (*it).second;
}

} // namespace SolarSystemTools
