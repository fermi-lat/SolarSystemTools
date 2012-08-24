/*
 * @file HealpixExposureSun.cxx
 * @brief Integral of effective area over time for the entire sky at
 * various energies.
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/src/HealpixExposureSun.cxx,v 1.4 2012/03/27 22:39:32 gudlaugu Exp $
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
                                   m_costhmin(0), m_costhmax(1),
																	 m_distCosCut(-1), m_avgDist(0),
                                   m_enforce_boundaries(false) {}

HealpixExposureSun::HealpixExposureSun(const std::vector<double> & energies,
                               const Observation & observation,
                               const st_app::AppParGroup * pars) 
   : m_energies(energies), m_observation(&observation),
     m_costhmin(0), m_costhmax(1), m_distCosCut(-1), m_avgDist(0), m_enforce_boundaries(false)  {

	 m_distCosCut = m_observation->expCubeSun().distCosCut();
	 m_avgDist = m_observation->expCubeSun().avgDist();
	 pixel::setStride(energies.size());
   if (pars) {
      setMapGeometry(*pars);
      setCosThetaBounds(*pars);
   } else {
      setMapGeometry();
   }
   computeMap();
}

HealpixExposureSun::HealpixExposureSun(const std::string & filename) 
   : m_observation(0), m_costhmin(0), m_costhmax(1),
		 m_distCosCut(-1), m_avgDist(0),
     m_enforce_boundaries(false) {

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

	 pixel::setStride(m_energies.size());

   std::auto_ptr<const tip::Table>
      thetasun(tip::IFileSvc::instance().readTable(filename, "thetasun"));

	 m_thetasun.clear();
	 it = thetasun->begin();
	 //Need to add the first value separately
	 double value;
	 row["Theta_Min"].get(value);
	 m_thetasun.push_back(value*M_PI/180.);
   for ( ; it != thetasun->end(); ++it) {
      row["Theta_Max"].get(value);
      m_thetasun.push_back(value*M_PI/180.);
   }

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

		hdr["AVGDIST"].get(m_avgDist);
		double distThetaCut;
		hdr["DISTTCUT"].get(distThetaCut);
		m_distCosCut = cos(distThetaCut*M_PI/180.);

    astro::SkyDir::CoordSystem coordsys = (check == "GAL")?
        astro::SkyDir::GALACTIC: astro::SkyDir::EQUATORIAL;

    m_exposureMap.setHealpix(healpix::Healpix(nside, ord, coordsys));

    tip::Table::ConstIterator itor = table.begin();
		healpix::HealpixArray<pixel>::iterator haitor = m_exposureMap.begin();

    for( ; itor != table.end(); ++haitor, ++itor)
    {
			  std::vector<float> values;
			  std::vector<size_t> index;
        (*itor)["Values"].get(values);
        (*itor)["Index"].get(index);
				assert( values.size() == index.size() );
				
				for (size_t i = 0; i < index.size(); ++i) {
					(*haitor).value(index[i]) = values[i];
				}

    }
    delete &table; 

}

double HealpixExposureSun::integrate(double energy, double ra, double dec, const Fun &f) const {
   std::vector<double>::const_iterator ie = ::findNearest(m_energies, energy);
   unsigned int k = ie - m_energies.begin();

	 //Get a reference to the pixel
	 const astro::SkyDir dir(ra,dec);
	 const pixel & pix = m_exposureMap[dir];

	 double integ(0);
	 std::vector<size_t> ifakes = pix.sparseIndices();
	 for ( size_t i = 0; i < ifakes.size(); ++i ) {
		 const float value = pix[i*m_energies.size() + k];
		 if ( value != 0 )
		    integ += value * f(cos(0.5*(m_thetasun[ifakes[i]]+m_thetasun[ifakes[i]+1])));
	 }

	 return integ;
}

double HealpixExposureSun::operator()(double energy, double ra, double dec, double theta) const {

	// Check for theta out of bounds
	 if (m_thetasun[m_thetasun.size()-1] < theta) return 0;

   std::vector<double>::const_iterator ie = ::findNearest(m_energies, energy);
   unsigned int k = ie - m_energies.begin();

   int l = 0;
   for ( ; l < m_thetasun.size()-1; ++l){
     if ( m_thetasun[l] < theta && theta < m_thetasun[l+1] )
       break;
   }

	 const astro::SkyDir dir(ra,dec);

	 const size_t indx = l*m_energies.size() + k;

	 return m_exposureMap[dir].value(indx);
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
   m_observation->expCubeSun().thetaBinsSun(m_thetasun);
   double binsz = pars["binsz"];
	 if (binsz > 0) {
		m_exposureMap.setHealpix(healpix::Healpix(side_from_degrees(binsz), healpix::Healpix::NESTED, astro::SkyDir::EQUATORIAL));
	 } else {
	  m_exposureMap.setHealpix(m_observation->expCubeSun().getHealpix());
	 }
}

void HealpixExposureSun::setMapGeometry() {
   m_observation->expCubeSun().thetaBinsSun(m_thetasun);
	 m_exposureMap.setHealpix(m_observation->expCubeSun().getHealpix());
}

void HealpixExposureSun::computeMap() {
	int iter(0);
   st_stream::StreamFormatter formatter("HealpixExposureSun", "computeMap", 2);
   formatter.warn() << "Computing binned exposure map";

   std::vector<double> costhetasun(m_thetasun.size()-1);
   for (int j = 0; j < costhetasun.size(); ++j) {
     costhetasun[j] = cos((m_thetasun[j]+m_thetasun[j+1])/2.);
   }

	 //The distance scaling cutoff
	 const double distCosCut = m_observation->expCubeSun().distCosCut();
	 const double avgDist = m_observation->expCubeSun().avgDist();
	 const double avgDist2 = avgDist*avgDist;

	 //Create a cache for AEff calculations
	 std::vector<BinnedExposureSun::Aeff*> aeffs(m_energies.size());
	 for (size_t i = 0; i < m_energies.size(); ++i)
     aeffs[i] = new BinnedExposureSun::Aeff(m_energies[i], *m_observation, m_costhmin, m_costhmax);

	 healpix::HealpixArray<pixel>::iterator haitor = m_exposureMap.begin();

	 for( ; haitor != m_exposureMap.end(); ++haitor, ++iter ) {

         if ((iter % ((m_exposureMap.size())/20)) == 0) {
            formatter.warn() << ".";
         }
         // std::pair<double, double> coord = m_proj->pix2sph(i + 1, j + 1);
         // astro::SkyDir dir(coord.first, coord.second, coordsys);
         astro::SkyDir dir = m_exposureMap.dir(haitor);
                                           
				 for (int l = 0; l < costhetasun.size(); ++l) {
					 if (m_observation->expCubeSun().hasCosthetasun(dir,costhetasun[l])) {
						 for (unsigned int k = 0; k < m_energies.size(); k++) {
							 const BinnedExposureSun::Aeff &aeff = *aeffs[k]; 
					     double exposure = m_observation->expCubeSun().value(dir, costhetasun[l], aeff, m_energies.at(k));
							 //Here we assume the distCosCut aligns with one of the bins as
							 //it should
							 if (costhetasun[l] > distCosCut)
								 exposure *= avgDist2;
							 const unsigned int indx = l*m_energies.size() + k;
							 m_exposureMap[dir].value(indx) += exposure;
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

	 table->appendField("Index", "1PV");
	 table->appendField("Values", "1PE");
   tip::Index_t numrecs =  m_exposureMap.size() ;
	 table->setNumRecords(numrecs);

	 // get iterators for the Table and the HealpixArray
	 tip::Table::Iterator itor = table->begin();
	 healpix::HealpixArray<pixel>::const_iterator haitor = m_exposureMap.begin();

	 // now just copy
	 for( ; haitor != m_exposureMap.end(); ++haitor, ++itor)
	 {
		 std::vector<size_t> index = (*haitor).indices();
		 (*itor)["Index"].set(index);
		 (*itor)["Values"].set(*haitor);
	 }

	 tip::Header & header(table->getHeader());

	 header["PIXTYPE"].set("HEALPIX"); 
	 header["ORDERING"].set("NESTED"); 
	 header["COORDSYS"].set(m_exposureMap.healpix().galactic()? "GAL" : "EQU");
	 header["NSIDE"].set(m_exposureMap.healpix().nside()); 
	 header["FIRSTPIX"].set(0); 
	 header["LASTPIX"].set(m_exposureMap.size() - 1); 
	 header["NENERGIES"].set(m_energies.size());
	 header["NTHBINS"].set(m_thetasun.size()-1);
	 header["NBINS"].set((m_thetasun.size()-1)*m_energies.size());
	 header["DISTTCUT"].set(180./M_PI*acos(m_distCosCut));
	 header["AVGDIST"].set(m_avgDist);
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

   ext = "THETASUN";
   tip::IFileSvc::instance().appendTable(filename, ext);
   table = tip::IFileSvc::instance().editTable(filename, ext);
   table->appendField("Theta_Min", "1D");
   table->appendField("Theta_Max", "1D");
   table->setNumRecords(m_thetasun.size()-1);

   row = table->begin();
   tip::Table::Record & record2 = *row;

   for (size_t i = 0; i < m_thetasun.size()-1; ++i, ++row) {
      record2["Theta_Min"].set(m_thetasun[i]*180./M_PI);
      record2["Theta_Max"].set(m_thetasun[i+1]*180./M_PI);
   }

   delete table;
}

void HealpixExposureSun::setCosThetaBounds(const st_app::AppParGroup & pars) {
   double thmin = pars["thmin"];
   if (thmin > 0.) {
      m_costhmax = std::cos(thmin*M_PI/180.);
   }
   double thmax = pars["thmax"];
   if (thmax < 180.) {
      m_costhmin = std::cos(thmax*M_PI/180.);
   }
}

size_t HealpixExposureSun::pixel::m_stride(1);

void HealpixExposureSun::pixel::setStride(size_t stride) {m_stride = stride;}

float & HealpixExposureSun::pixel::value(size_t i) {
	const size_t ifake = i / m_stride;
	pixelIndex::iterator it = std::lower_bound(m_ifaketoreal.begin(), m_ifaketoreal.end(), std::pair<size_t,size_t>(ifake,0), less_than);

	//Add the bin if needed
	size_t ireal(0);
	if ( it == m_ifaketoreal.end() || it->first != ifake ) {
		ireal = m_ifaketoreal.size();
		m_ifaketoreal.insert(it, std::pair<size_t,size_t>(ifake,ireal));
		m_irealtofake.insert(m_irealtofake.end(), std::pair<size_t, size_t>(ireal,ifake));

		//Resize the storage
		resize(size()+m_stride, 0.0);
	} else {
		ireal = it->second;
	}
	return at(i + (ireal - ifake)*m_stride);
}

float HealpixExposureSun::pixel::value(size_t i) const {
	const size_t ifake = i / m_stride;
	pixelIndex::const_iterator it = std::lower_bound(m_ifaketoreal.begin(), m_ifaketoreal.end(), std::pair<size_t,size_t>(ifake,0), less_than);

	//Return 0 if bin not found
	if ( it == m_ifaketoreal.end() || it->first != ifake )  return 0;

  const size_t ireal = it->second;
	return at(i + (ireal - ifake)*m_stride);
}

std::vector<size_t> HealpixExposureSun::pixel::indices() const {
	std::vector<size_t> indices(size());
	for (size_t i = 0; i < m_irealtofake.size(); ++i) {
		const size_t irealstride = m_irealtofake[i].second*m_stride;
		for (size_t j = 0; j < m_stride; ++j) {
			indices[i*m_stride+j] = irealstride + j;
		}
	}
	return indices;
}

std::vector<size_t> HealpixExposureSun::pixel::sparseIndices() const {
	std::vector<size_t> indices(m_irealtofake.size());
	for (size_t i = 0; i < m_irealtofake.size(); ++i) {
		indices[i] = m_irealtofake[i].second;
	}
	return indices;
}

} // namespace SolarSystemTools
