/**
 * @file SolarTemplate.cxx
 * @brief Methods to calculate a solar template for use in likelihood analysis
 * from solar binned exposure and solar profile
 * @author G. Johannesson
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/src/SolarTemplate.cxx,v 1.2 2012/03/21 22:50:20 gudlaugu Exp $
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

#include "SolarSystemTools/SolarTemplate.h"
#include "Likelihood/CountsMap.h"
#include "Likelihood/BinnedExposure.h"
#include "SolarSystemTools/HealpixExposureSun.h"
#include "SolarSystemTools/SolarProfile.h"


namespace SolarSystemTools {

SolarTemplate::SolarTemplate(const Likelihood::CountsMap & cmap, 
                               const HealpixExposureSun & hpexpsun, 
                               const Likelihood::BinnedExposure & avgexp, 
                               const SolarProfile & prof,
                               bool useEbounds,
                               const st_app::AppParGroup * pars)
   : m_expsun(hpexpsun), m_avgexp(avgexp), m_profile(prof), m_proj(0) {
   setMapGeometry(cmap);
   if (!useEbounds) {
      for (size_t k(0); k < m_energies.size()-1; k++) {
         m_energies[k] = std::sqrt(m_energies[k]*m_energies[k+1]);
      }
      m_energies.pop_back();
      m_naxes[2] -= 1;
   }
   computeMap();
}

SolarTemplate::SolarTemplate(const std::vector<double> & energies,
                               const HealpixExposureSun & hpexpsun, 
                               const Likelihood::BinnedExposure & avgexp, 
                               const SolarProfile & prof,
                               const st_app::AppParGroup * pars) 
   : m_energies(energies), m_expsun(hpexpsun), m_avgexp(avgexp), m_profile(prof), m_proj(0)  {
   if (pars) {
      setMapGeometry(*pars);
   } else {
      setMapGeometry();
   }
   computeMap();
}

SolarTemplate::~SolarTemplate() {
   delete m_proj;
}


void SolarTemplate::setMapGeometry(const Likelihood::CountsMap & cmap) {
   m_proj_name = cmap.proj_name();
   m_crpix[0] = cmap.crpix1();
   m_crpix[1] = cmap.crpix2();
   m_crval[0] = cmap.crval1();
   m_crval[1] = cmap.crval2();
   m_cdelt[0] = cmap.cdelt1();
   m_cdelt[1] = cmap.cdelt2();
   m_crota2 = cmap.crota2();
   m_naxes.resize(3, 0);
   m_naxes[0] = cmap.naxis1();
   m_naxes[1] = cmap.naxis2();
   m_naxes[2] = cmap.energies().size();
   m_energies = cmap.energies();
   m_isGalactic = cmap.isGalactic();
}

void SolarTemplate::setMapGeometry(const st_app::AppParGroup & pars) {
   m_naxes.resize(3, 0);
   m_naxes[0] = pars["nxpix"];
   m_naxes[1] = pars["nypix"];
   m_naxes[2] = m_energies.size();
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

void SolarTemplate::setMapGeometry() {
   m_proj_name = "CAR";
   m_isGalactic = false;
   m_naxes.resize(3, 0);
   m_naxes[0] = 360;
   m_naxes[1] = 180;
   m_naxes[2] = m_energies.size();
   m_crpix[0] = m_naxes[0]/2. + 0.5;
   m_crpix[1] = m_naxes[1]/2. + 0.5;
   m_crval[0] = 180.;
   m_crval[1] = 0;
   m_cdelt[0] = -1;
   m_cdelt[1] = 1;
   m_crota2 = 0;
}

SolarTemplate::ProfFun::ProfFun(const std::vector<double> & costhsun, 
					       const std::vector<double> & energies, 
					       const std::vector<double> & cache, size_t ie ) :
	m_costhsun(costhsun), m_energies(energies), m_cache(cache), m_ie(ie*(costhsun.size()-1)) {}

double SolarTemplate::ProfFun::operator() (double costhetasun) const {
	const std::vector<double>::const_iterator it = 
		std::upper_bound(m_costhsun.begin(), m_costhsun.end(), costhetasun, std::greater<double>());
	const size_t l = (it - m_costhsun.begin()) - 1;
	return m_cache[m_ie+l];
}


void SolarTemplate::computeMap() {
   m_proj = new astro::SkyProj(m_proj_name, &m_crpix[0], &m_crval[0],
                               &m_cdelt[0], m_crota2, m_isGalactic);

   m_template.resize(m_naxes.at(0)*m_naxes.at(1)*m_energies.size(), 0);
   int iter(0);
   st_stream::StreamFormatter formatter("SolarTemplate", "computeMap", 2);
   formatter.warn() << "Computing binned solarTemplate map";

   const std::vector<double> &thetasun = m_expsun.thetasun();
   std::vector<double> costhetasun(thetasun.size());
   for (size_t i = 0; i < costhetasun.size(); ++i){
	   costhetasun[i] = cos(thetasun[i]);
	 }

	 //Create the average intensity as a function of energy and angle
	 std::vector<double> intensityCache(m_energies.size()*(costhetasun.size()-1));
   for (unsigned int k = 0; k < m_energies.size(); k++) {
     for (int l = 0; l < costhetasun.size()-1; ++l) {
			   const size_t cInd = k*(costhetasun.size()-1) + l;
         intensityCache[cInd] = m_profile.average(costhetasun[l+1], costhetasun[l], m_energies[k]);
     }
	 }



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
					 //Get the average expsoure for the energy and pixel
					 const unsigned int indx = (k*m_naxes.at(1) + j)*m_naxes.at(0) + i;
					 const double avgExp = m_avgexp(m_energies[k], dir.ra(), dir.dec());
					 //Integrate the solar angles
	         const ProfFun f(costhetasun, m_energies, intensityCache, k);
					 m_template.at(indx) = m_expsun.integrate(m_energies[k], dir.ra(), dir.dec(), f);
					 m_template.at(indx) /= avgExp;  
         }
      }
   }
   formatter.warn() << "!" << std::endl;
}

void SolarTemplate::writeOutput(const std::string & filename) const {
   std::remove(filename.c_str());

   std::string ext("PRIMARY");
   tip::IFileSvc::instance().appendImage(filename, ext, m_naxes);
   tip::Image * image = tip::IFileSvc::instance().editImage(filename, ext);

   image->set(m_template);

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
}


} // namespace SolarSystemTools
