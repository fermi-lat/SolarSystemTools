/**
 * @file HealpixExposureSun.cxx
 * @brief Integral of effective area over time for the entire sky at
 * various energies.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/src/HealpixExposureSun.cxx,v 1.2 2012/02/13 22:24:37 gudlaugu Exp $
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

#include "SolarSystemTools/HealpixExposureSun.h"
#include "SolarSystemTools/BinnedExposureSun.h"
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
         what << "HealpixExposureSun::operator(): The energy " << x
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

HealpixExposureSun::HealpixExposureSun() : m_observation(0),
                                   m_costhmin(-1), m_costhmax(1),
                                   m_enforce_boundaries(false) {}

HealpixExposureSun::HealpixExposureSun(const std::vector<double> & energies,
                               const Observation & observation,
                               const st_app::AppParGroup * pars) 
   : m_energies(energies), m_observation(&observation),
     m_costhmin(-1), m_costhmax(1),m_enforce_boundaries(false)  {
   if (pars) {
      setMapGeometry(*pars);
      setCosThetaBounds(*pars);
   } else {
      setMapGeometry();
   }
   computeMap();
}

HealpixExposureSun::HealpixExposureSun(const std::string & filename) 
   : m_observation(0), m_costhmin(-1), m_costhmax(1),
     m_enforce_boundaries(false) {

    const tip::Table & table=*tip::IFileSvc::instance().readTable(filename, "EXPOSURE");
    const tip::Header& hdr = table.getHeader();
    int nside=0;
    hdr["NSIDE"].get(nside);
    std::string ordering;
    hdr["ORDERING"].get(ordering);
		healpix::Healpix::Ordering ord = (ordering == "NESTED")?
        healpix::Healpix::NEST: healpix::Healpix::RING;

    // Code for setting CoordSystem added 1/17/2008
    std::string check;
    try{ // new value
        hdr["COORDSYS"].get(check);
    }catch(const std::exception &){
        // allow old value
        hdr["COORDTYPE"].get(check);
    }

    astro::SkyDir::CoordSystem coordsys = (check == "GAL")?
        astro::SkyDir::GALACTIC: astro::SkyDir::EQUATORIAL;

    m_exposureMap.setHealpix(healpix::Healpix(nside, ord, coordsys));

    tip::Table::ConstIterator itor = table.begin();
		healpix::HealpixArray<pixel>::iterator haitor = m_exposureMap.begin();

    for( ; itor != table.end(); ++haitor, ++itor)
    {
			  std::vector<double> values;
			  std::vector<size_t> index;
        (*itor)["Values"].get(values);
        (*itor)["Index"].get(index);
				assert( values.size() == index.size() );
				
				for (size_t i = 0; i < index.size(); ++i) {
					(*haitor)[index[i]] = values[i];
				}

    }
    delete &table; 

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

   std::auto_ptr<const tip::Table>
      costhetasun(tip::IFileSvc::instance().readTable(filename, "Costhetasun"));

	 m_costhetasun.clear();
	 it = costhetasun->begin();
	 //Need to add the first value separately
	 double value;
	 row["CosTheta_Min"].get(value);
	 m_costhetasun.push_back(value);
   for ( ; it != costhetasun->end(); ++it) {
      row["CosTheta_Max"].get(value);
      m_costhetasun.push_back(value);
   }

}

double HealpixExposureSun::integrate(double energy, double ra, double dec, const Fun &f) const {
   std::vector<double>::const_iterator ie = ::findNearest(m_energies, energy);
   unsigned int k = ie - m_energies.begin();

	 //Get a reference to the pixel
	 const astro::SkyDir dir(ra,dec);
	 const pixel & pix = m_exposureMap[dir];

	 double integ(0);
	 pixel::const_iterator it = pix.lower_bound(k*(m_costhetasun.size()-1));
	 pixel::const_iterator end = pix.upper_bound((k+1)*(m_costhetasun.size()-1)-1);
	 for ( ; it != end; ++it ) {
		 const size_t cind = (*it).first-k*(m_costhetasun.size()-1);
		 integ += (*it).second * f(0.5*(m_costhetasun[cind]+m_costhetasun[cind+1]));
	 }

	 return integ;
}

double HealpixExposureSun::operator()(double energy, double ra, double dec, double costheta) const {
   std::vector<double>::const_iterator ie = ::findNearest(m_energies, energy);
   unsigned int k = ie - m_energies.begin();

   int l = 0;
   for ( ; l < m_costhetasun.size()-1; ++l){
     if ( m_costhetasun[l] < costheta && costheta < m_costhetasun[l+1] )
       break;
   }

	 const astro::SkyDir dir(ra,dec);

	 const unsigned int indx = l + k*(m_costhetasun.size()-1);

	 pixel::const_iterator it = m_exposureMap[dir].find(indx);

	 if (it != m_exposureMap[dir].end()) {
		 return (*it).second;
	 } else {
		 return 0;
	 }
}

/// return the closest power of 2 for the side parameter
/// 41252 square degrees for the sphere
/// nside=64 is about one degee: 49152 pixels
inline int side_from_degrees(double pixelsize){ 
    int n = 1;
    while( 12*n*n < 41252/(pixelsize*pixelsize) && n < 512){
        n *=2;
    }
    return n; 
} 

void HealpixExposureSun::setMapGeometry(const st_app::AppParGroup & pars) {
   m_observation->expCubeSun().costhetaBinsSun(m_costhetasun);
   double binsz = pars["binsz"];
	 if (binsz > 0) {
		m_exposureMap.setHealpix(healpix::Healpix(side_from_degrees(binsz), healpix::Healpix::NESTED, astro::SkyDir::EQUATORIAL));
	 } else {
	  m_exposureMap.setHealpix(m_observation->expCubeSun().getHealpix());
	 }
}

void HealpixExposureSun::setMapGeometry() {
   m_observation->expCubeSun().costhetaBinsSun(m_costhetasun);
	 m_exposureMap.setHealpix(m_observation->expCubeSun().getHealpix());
}

void HealpixExposureSun::computeMap() {
	int iter(0);
   st_stream::StreamFormatter formatter("HealpixExposureSun", "computeMap", 2);
   formatter.warn() << "Computing binned exposure map";

   std::vector<double> costhetasun(m_costhetasun.size()-1);
   for (int j = 0; j < costhetasun.size(); ++j) {
     costhetasun[j] = (m_costhetasun[j]+m_costhetasun[j+1])/2.;
   }

	 //Create a cache for AEff calculations
	 std::vector<BinnedExposureSun::Aeff*> aeffs;
	 for (size_t i = 0; i < m_energies.size(); ++i)
     aeffs.push_back(new BinnedExposureSun::Aeff(m_energies[i], *m_observation, m_costhmin, m_costhmax));

	 healpix::HealpixArray<pixel>::iterator haitor = m_exposureMap.begin();

	 for( ; haitor != m_exposureMap.end(); ++haitor, ++iter ) {

         if ((iter % ((m_exposureMap.size())/20)) == 0) {
            formatter.warn() << ".";
         }
         // std::pair<double, double> coord = m_proj->pix2sph(i + 1, j + 1);
         // astro::SkyDir dir(coord.first, coord.second, coordsys);
         astro::SkyDir dir = m_exposureMap.dir(haitor);
                                           
         for (unsigned int k = 0; k < m_energies.size(); k++) {
               const BinnedExposureSun::Aeff &aeff = *aeffs[k]; 
               for (int l = 0; l < costhetasun.size(); ++l) {
									 const double exposure = m_observation->expCubeSun().value(dir, costhetasun[l], aeff, m_energies.at(k));
									 if (exposure != 0) {
								      const unsigned int indx = l + k*costhetasun.size();
									    m_exposureMap[dir][indx] += exposure;
									 }
               }
         }
   }
   formatter.warn() << "!" << std::endl;
	 for (size_t i = 0; i < m_energies.size(); ++i)
     delete aeffs[i];
}

void HealpixExposureSun::writeOutput(const std::string & filename) const {
   std::remove(filename.c_str());

   std::string ext("EXPOSURE");
   tip::IFileSvc::instance().appendTable(filename, ext);
   tip::Table * table = tip::IFileSvc::instance().editTable(filename, ext);

	 table->appendField("Index", "1PJ");
	 table->appendField("Values", "1PE");
   tip::Index_t numrecs =  m_exposureMap.size() ;
	 table->setNumRecords(numrecs);

	 // get iterators for the Table and the HealpixArray
	 tip::Table::Iterator itor = table->begin();
	 healpix::HealpixArray<pixel>::const_iterator haitor = m_exposureMap.begin();

	 // now just copy
	 for( ; haitor != m_exposureMap.end(); ++haitor, ++itor)
	 {
		 std::vector<double> values((*haitor).size());
		 std::vector<size_t> index((*haitor).size());
		 pixel::const_iterator it = (*haitor).begin();
		 for (size_t i=0; it != (*haitor).end(); ++it, ++i){
			 index[i] = (*it).first;
			 values[i] = (*it).second;
		 }
		 (*itor)["Values"].set(values);
		 (*itor)["Index"].set(index);
	 }

	 tip::Header & header(table->getHeader());

	 header["PIXTYPE"].set("HEALPIX"); 
	 header["ORDERING"].set("NESTED"); 
	 header["COORDSYS"].set(m_exposureMap.healpix().galactic()? "GAL" : "EQU");
	 header["NSIDE"].set(m_exposureMap.healpix().nside()); 
	 header["FIRSTPIX"].set(0); 
	 header["LASTPIX"].set(m_exposureMap.size() - 1); 
	 header["NENERGIES"].set(m_energies.size());
	 header["NTHBINS"].set(m_costhetasun.size()-1);
   header["TELESCOP"].set("GLAST");
   header["INSTRUME"].set("LAT");
   header["DATE-OBS"].set("");
   header["DATE-END"].set("");

   delete table;

   ext = "ENERGIES";
   tip::IFileSvc::instance().appendTable(filename, ext);
   table = tip::IFileSvc::instance().editTable(filename, ext);
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
   table->setNumRecords(m_costhetasun.size()-1);

   row = table->begin();
   tip::Table::Record & record2 = *row;

   for (size_t i = 0; i < m_costhetasun.size()-1; ++i, ++row) {
      record2["CosTheta_Min"].set(m_costhetasun[i]);
      record2["CosTheta_Max"].set(m_costhetasun[i+1]);
   }

   delete table;
}

void HealpixExposureSun::setCosThetaBounds(const st_app::AppParGroup & pars) {
   double thmin = pars["thmin"];
   if (thmin > 0) {
      m_costhmax = std::cos(thmin*M_PI/180.);
   }
   double thmax = pars["thmax"];
   if (thmax < 180.) {
      m_costhmin = std::cos(thmax*M_PI/180.);
   }
}

} // namespace SolarSystemTools
