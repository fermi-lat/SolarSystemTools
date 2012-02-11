/** @file ExposureSun.cxx
    @brief Implementation of class ExposureSun
		@author G. Johannesson
    
		$Header: $
*/
#include "SolarSystemTools/ExposureSun.h"
#include "healpix/HealpixArrayIO.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"
#include "astro/EarthCoordinate.h"
#include "astro/PointingTransform.h"
#include "astro/JulianDate.h"

#include <memory>
#include <algorithm>

using namespace SolarSystemTools;
using healpix::Healpix;

const double ExposureSun::s_mjd_missionStart(astro::JulianDate::missionStart());


ExposureSun::ExposureSun(const std::string& inputfile, const std::string& tablename)
: SkyExposure2D(SkyBinner2D(2))
{
	load(inputfile, tablename);
}

void ExposureSun::load(const std::string& inputFile, const std::string& tablename) 
{
    const tip::Table & table=*tip::IFileSvc::instance().readTable(inputFile, tablename);
    const tip::Header& hdr = table.getHeader();
    int nside=0;
    hdr["NSIDE"].get(nside);
    std::string ordering;
    hdr["ORDERING"].get(ordering);
    Healpix::Ordering ord = (ordering == "NESTED")?
        Healpix::NEST: Healpix::RING;
    std::string thetabinstring;
    hdr["THETABIN"].get(thetabinstring);
    bool thetabin = (thetabinstring == "COSTHETA")? false: true;
    int nbrbins;
    hdr["NBRBINS"].get(nbrbins);
    int nbrbins1;
    hdr["NBRBINS1"].get(nbrbins1);
    double cosmin1;
    hdr["COSMIN1"].get(cosmin1);
    int nbrbins2;
    hdr["NBRBINS2"].get(nbrbins2);
    double cosmin2;
    hdr["COSMIN2"].get(cosmin2);

    int nphibins(0);
    try{
        hdr["PHIBINS"].get(nphibins);
    }catch(...){}

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

    CosineBinner2D::setBinning(cosmin1, cosmin2, nbrbins1, nbrbins2, thetabin);
    if( nphibins>0){
        CosineBinner2D::setPhiBins(nphibins);
    }
    data().setHealpix(Healpix(nside, ord, coordsys));

    tip::Table::ConstIterator itor = table.begin();
		healpix::HealpixArray<CosineBinner2D>::iterator haitor = data().begin();

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
		create_cache();
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

ExposureSun::ExposureSun(double pixelsize, double cosbinsize, double cosbinsizesun, double zcut, bool weighted)
: SkyExposure2D(
    SkyBinner2D(Healpix(
      side_from_degrees(pixelsize),  // nside
      Healpix::NESTED, 
      astro::SkyDir::EQUATORIAL) )
  )
, m_zcut(zcut), m_lost(0)
, m_weighted(weighted)
, m_solar_dir(astro::SolarSystem::SUN)
{
    unsigned int cosbins = static_cast<unsigned int>(1./cosbinsize);
    unsigned int cosbinssun = static_cast<unsigned int>(2./cosbinsizesun);
    if( cosbins != CosineBinner2D::nbins1() || cosbinssun != CosineBinner2D::nbins2() ) { 
        unsigned int allbins(cosbins*cosbinssun);
        if( CosineBinner2D::nphibins() > 0 ){
            // add size for extra phi bins.
            allbins += allbins * CosineBinner2D::nphibins();
        }
        
        CosineBinner2D::setBinning(0, -1, cosbins, cosbinssun);
    }
    create_cache();
}

void ExposureSun::create_cache()
{
    size_t datasize(data().size());
    m_dir_cache.reserve(datasize);

    SkyBinner2D::iterator is = data().begin();
    for( ; is != data().end(); ++is){ // loop over all pixels
        Simple3Vector pixdir(data().dir(is)());
        m_dir_cache.push_back(std::make_pair(&*is, pixdir));
    }
}

/** @class Filler
    @brief private helper class used in for_each to fill a CosineBinner object
*/
class ExposureSun::Filler {
public:

    /** @brief ctor
        @param deltat time to add
        @param dirz direction to use to determine angle (presumably the spacecraft z-axis)
        @param dirx direction of x-axis
        @param zenith optional zenith direction for potential cut
        @param zcut optional cut: if -1, ignore
    */
    Filler( double deltat, const astro::SkyDir& dirz, const astro::SkyDir& dirx, const astro::SkyDir& dirsun, astro::SkyDir zenith=astro::SkyDir(), double zcut=-1)
        : m_dirz(dirz())
        , m_rot(astro::PointingTransform(dirz,dirx).localToCelestial().inverse())
				, m_dirsun(dirsun())
        , m_zenith(zenith())
        , m_deltat(deltat)
        , m_zcut(zcut)
        , m_total(0), m_lost(0)
        , m_use_phi(CosineBinner2D::nphibins()>0)
    {}
    Filler( double deltat, const astro::SkyDir& dirz, const astro::SkyDir& dirsun, astro::SkyDir zenith=astro::SkyDir(), double zcut=-1)
        : m_dirz(dirz())
				, m_dirsun(dirsun())
        , m_zenith(zenith())
        , m_deltat(deltat)
        , m_zcut(zcut)
        , m_total(0), m_lost(0)
        , m_use_phi(false)
    {}
    void operator()( const std::pair<CosineBinner2D*, Simple3Vector> & x)
    {
        // check if we are making a horizon cut:
        bool ok( m_zcut==-1);
        if( ! ok) {
            double z(x.second.dot(m_zenith));
            ok = z > m_zcut;
        }
        if( ok) {
            // if ok, add to the angle histogram
            const Simple3Vector& pixeldir(x.second);
            if( m_use_phi) {
                CLHEP::Hep3Vector instrument_dir( pixeldir.transform(m_rot) );
                 double costheta(instrument_dir.z()), phi(instrument_dir.phi());
                x.first->fill( costheta, pixeldir.dot(m_dirsun), phi , m_deltat);
						} else {
							  x.first->fill( pixeldir.dot(m_dirz), pixeldir.dot(m_dirsun), m_deltat);
						}

            m_total += m_deltat;
        }else{
            m_lost += m_deltat;
        }
    }
    double total()const{return m_total;}
    double lost()const{return m_lost;}
private:
    Simple3Vector m_dirz, m_dirsun;
    CLHEP::HepRotation m_rot;
    Simple3Vector m_zenith;
    double m_deltat, m_zcut;
    mutable double m_total, m_lost;
    bool m_use_phi;
};

void ExposureSun::fill(const astro::SkyDir& dirz, double deltat)
{
}

void ExposureSun::fill(const astro::SkyDir& dirz, const astro::SkyDir& dirsun, double deltat)
{
    Filler sum = for_each(m_dir_cache.begin(), m_dir_cache.end(), Filler(deltat, dirz, dirsun));
    double total(sum.total());
    addtotal(total);
    m_lost += sum.lost();
}


void ExposureSun::fill(const astro::SkyDir& dirz, const astro::SkyDir& dirsun, const astro::SkyDir& zenith, double deltat)
{
    Filler sum = for_each(m_dir_cache.begin(), m_dir_cache.end(), Filler(deltat, dirz, dirsun, zenith, m_zcut));
    double total(sum.total());
    addtotal(total);
    m_lost += sum.lost();
}

void ExposureSun::fill_zenith(const astro::SkyDir& dirz,const astro::SkyDir& dirx, 
		                       const astro::SkyDir& dirsun,
                           const astro::SkyDir& zenith, 
                           double deltat)
{
    Filler sum = for_each(m_dir_cache.begin(), m_dir_cache.end(), 
        Filler(deltat, dirz, dirx, dirsun, zenith, m_zcut));
    double total(sum.total());
    addtotal(total);
    m_lost += sum.lost();
}


void ExposureSun::write(const std::string& outputfile, const std::string& tablename)const
{
    // add a table to the file
    tip::IFileSvc::instance().appendTable(outputfile, tablename);
    tip::Table & table = *tip::IFileSvc::instance().editTable( outputfile, tablename);

    // this is a work-around for a bug in tip v2r1p1

    table.appendField("Index", "1PJ");
    table.appendField("Values", "1PE");
    tip::Index_t numrecs =  data().size() ;
    table.setNumRecords(numrecs);

    // get iterators for the Table and the HealpixArray
    tip::Table::Iterator itor = table.begin();
		healpix::HealpixArray<CosineBinner2D>::const_iterator haitor = data().begin();

    // now just copy
    for( ; haitor != data().end(); ++haitor, ++itor)
    {
			  std::vector<double> values((*haitor).size());
				std::vector<size_t> index((*haitor).size());
        CosineBinner2D::const_iterator it = (*haitor).begin();
				for (size_t i=0; it != (*haitor).end(); ++it, ++i){
					index[i] = (*it).first;
					values[i] = (*it).second;
				}
        (*itor)["Values"].set(values);
        (*itor)["Index"].set(index);
    }

    // set the headers (TODO: do the comments, too)
    tip::Header& hdr = table.getHeader();

    hdr["PIXTYPE"].set("HEALPIX"); 
    hdr["ORDERING"].set("NESTED"); 
    hdr["COORDSYS"].set(data().healpix().galactic()? "GAL" : "EQU");
    hdr["NSIDE"].set(data().healpix().nside()); 
    hdr["FIRSTPIX"].set(0); 
    hdr["LASTPIX"].set(data().size() - 1); 
    hdr["THETABIN"].set(CosineBinner2D::thetaBinning());
    hdr["NBRBINS"].set(CosineBinner2D::nbins1()*CosineBinner2D::nbins2());
    hdr["NBRBINS1"].set(CosineBinner2D::nbins1());
    hdr["COSMIN1"].set(CosineBinner2D::cosmin1());
    hdr["NBRBINS2"].set(CosineBinner2D::nbins2());
    hdr["COSMIN2"].set(CosineBinner2D::cosmin2());
    hdr["PHIBINS"].set(CosineBinner2D::nphibins());

    // need to do this to ensure file is closed when pointer goes out of scope
    delete &table;
}

void ExposureSun::load(const tip::Table * scData, 
                    const GTIvector& gti, 
                    bool verbose) {
   
   tip::Table::ConstIterator it = scData->begin();
   const tip::ConstTableRecord & row = *it;
   long nrows = scData->getNumRecords();

   for (long irow = 0; it != scData->end(); ++it, ++irow) {
      if (verbose && (irow % (nrows/20)) == 0 ) std::cerr << ".";
      if( processEntry( row, gti) )break;
   }
   if (verbose) std::cerr << "!" << std::endl;
}


bool ExposureSun::processEntry(const tip::ConstTableRecord & row, const GTIvector& gti)
{
    using astro::SkyDir;

    double  start, stop, livetime; 
    row["livetime"].get(livetime);
    if(livetime==0 ) return false; // assume this takes care of any entries during SAA
    row["start"].get(start);
    row["stop"].get(stop);
    double deltat = livetime; 


    double fraction(1); 
    bool  done(false);
    if( !gti.empty() ) {
        fraction = 0;

        GTIvector::const_iterator it  = gti.begin();
        for ( ; it != gti.end(); ++it) {
            double first = it->first,
                second=it->second;

            if( start < first ) {
                if( stop < first) continue; // history interval before gti
                if( stop <= second){
                    fraction = (stop-first)/(stop-start); // overlap start of gti
                    break;
                }
                fraction = (second-first)/(stop-start); // gti subset of history
                break;
            }else {
                if( start > second) continue; // interval after gti 
                if( stop <= second ) {
                    fraction = 1.0; break;  // fully contained
                }
                fraction = (second-start)/(stop-start);  // overlap end of gti
                break;
            }
 
        }
        done = fraction==0 && start > gti.back().second; 
    }
    if( fraction>0. ) {
        deltat *= fraction; // reduce if a boundary
        double ra, dec, razenith, deczenith;
        row["ra_scz"].get(ra);
        row["dec_scz"].get(dec);
        SkyDir scz(ra, dec);
        row["ra_scx"].get(ra);
        row["dec_scx"].get(dec);
        SkyDir scx(ra, dec);
	row["ra_zenith"].get(razenith);
	row["dec_zenith"].get(deczenith);
        SkyDir zenith(razenith, deczenith);
				//Usa astro to calculate the direction to the sun at center of bin
				const double mjd = (start+stop)/2./86400. + s_mjd_missionStart;
				SkyDir scsun(m_solar_dir.direction(mjd));
        if( m_weighted ){
            // adjust time by multiplying by livetime fraction
            deltat *= livetime/(stop-start);
        }
        fill_zenith(scz, scx, scsun,  zenith, deltat);
    }
    return done; 

}

