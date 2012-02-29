/** @file ExposureSun.h
    @brief definition of the class ExposureSun

    @author G. Johannesson

		$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/SolarSystemTools/SolarSystemTools/ExposureSun.h,v 1.1.1.1 2012/02/11 02:26:40 gudlaugu Exp $
*/
#ifndef SolarSystemTools_EXPOSURE_SUN_H
#define SolarSystemTools_EXPOSURE_SUN_H


#include "astro/SkyDir.h"
#include "astro/SolarSystem.h"
#include "healpix/HealpixArray.h"
#include "SolarSystemTools/CosineBinner2D.h"
#include "map_tools/Exposure.h"

#include <utility> // for std::pair


// define Exposure as specific instantiation of the above
typedef healpix::HealpixArray<SolarSystemTools::CosineBinner2D> SkyBinner2D;
typedef BasicExposure<SkyBinner2D, SolarSystemTools::CosineBinner2D> SkyExposure2D;

namespace SolarSystemTools {

/**
@class Exposure2D
@brief Manage a differential exposure database.

It is a pixelated using Healpix binning, and the CosineBinner2D class

$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/SolarSystemTools/SolarSystemTools/ExposureSun.h,v 1.1.1.1 2012/02/11 02:26:40 gudlaugu Exp $
*/

class ExposureSun : public SkyExposure2D {
public:
    //! create object with specified binning
    //! @param pixelsize (deg) Approximate pixel size, in degrees
    //! @param cosbinsize bin size in the cos(theta) binner
    //! @param weighted [false] set true to make a weighted table
    ExposureSun(double pixelsize=1., 
				double cosbinsize=1./CosineBinner2D::nbins1(), 
				double cosbinsizesun=1./CosineBinner2D::nbins2(), 
        double zcut=-1.0,
        bool   weighted=false
        );

    //! add a time interval at the given position, should never be used
    virtual void fill(const astro::SkyDir& dirz, double deltat);

    //! add a time interval at the given position
    virtual void fill(const astro::SkyDir& dirz, const astro::SkyDir& dirsun, double deltat);

    //! add a time interval at the given position
    virtual void fill(const astro::SkyDir& dirz, const astro::SkyDir& dirsun, const astro::SkyDir& dirzenith, double deltat);

    //! create object from the data file (FITS for now)
    ExposureSun(const std::string& inputfile, const std::string& tablename="ExposureSun");
    void load(const std::string& inputfile, const std::string& tablename="ExposureSun");

    //! Integral for costhetasun positions
    template<class F>
       double operator()(const astro::SkyDir&dir, double costhetasun, const F& fun) const
       {
          const SolarSystemTools::CosineBinner2D binner = data()[dir];
          return binner(fun,costhetasun);
       }

    template<class F>
       double integral(const astro::SkyDir&dir, double costhetasun, const F& fun) const
       {
          const SolarSystemTools::CosineBinner2D binner = data()[dir];
          return binner.integral(fun, costhetasun);
       }

    //! write out to a file.
    void write(const std::string& outputfile, const std::string& tablename="ExposureSun")const;

    typedef std::vector<std::pair<double, double> > GTIvector;

    //! load a set of history intervals from a table, qualified by a set of "good-time" intervals 
    void load(const tip::Table * scData, 
        const GTIvector & gti= GTIvector(), 
                    bool verbose=true);

    double lost()const{return m_lost;}

    /** @brief  allow horizon cut, possible if FOV includes horizon
        @param dirz direction of z-axis of instrument
        @param dirx direction of x-axis of instrument
        @param dirzenith direction of local zenith
        @param deltat time interval
        @param fraction livetime fractio
    */
    virtual void fill_zenith(const astro::SkyDir& dirz,const astro::SkyDir& dirx, const astro::SkyDir& dirsun,
        const astro::SkyDir& dirzenith, double deltat);

private:
    bool processEntry(const tip::ConstTableRecord & row, const GTIvector& gti);

		void create_cache();

        /** @class Simple3Vector 
    @brief replacement for Hep3Vector for speed of dot product

    */
    class Simple3Vector {
    public: 
        Simple3Vector(const CLHEP::Hep3Vector& v=CLHEP::Hep3Vector())
            : x(v.x()),y(v.y()),z(v.z()){};
        Simple3Vector(double a, double b, double c):x(a),y(b), z(c){}
        double dot(const Simple3Vector& u)const{return x*u.x+y*u.y+z*u.z;}
        CLHEP::Hep3Vector transform(const CLHEP::HepRotation& R)const{
            return CLHEP::Hep3Vector(
                R.xx()*x+R.xy()*y+R.xz()*z,
                R.yx()*x+R.yy()*y+R.yz()*z,
                R.zx()*x+R.zy()*y+R.zz()*z);
        }
        double x,y,z;
    };
    std::vector< std::pair<SolarSystemTools::CosineBinner2D* ,  Simple3Vector> > m_dir_cache;
		astro::SolarSystem m_solar_dir;
    class Filler ; ///< class used to fill a CosineBinner2D object with a value

    double m_zcut; ///< value for zenith angle cut
    double m_lost; ///< keep track of lost
    bool   m_weighted; ///< true if accumulating weighted livetime
		static const double s_mjd_missionStart;
};


} // namespace SolarSystemTools
#endif
