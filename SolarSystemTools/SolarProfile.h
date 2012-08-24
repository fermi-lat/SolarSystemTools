/**
 * @file SolarProfile.h
 * @brief Access to intensity as a function of angle and energy
 * @author G. Johannesson
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/SolarSystemTools/SolarProfile.h,v 1.1 2012/02/15 03:03:50 gudlaugu Exp $
 */

#ifndef SolarSystemTools_SolarProfile_h
#define SolarSystemTools_SolarProfile_h

#include <vector>

namespace SolarSystemTools {

/**
 * @class SolarProfile
 * @brief This class encapsulates the access to a solar profile, intensity as a
 * function of angle from the sun and energy.
 *
 * @author G. Johannesson
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/SolarSystemTools/SolarProfile.h,v 1.1 2012/02/15 03:03:50 gudlaugu Exp $
 */

class SolarProfile {

public:

	 virtual ~SolarProfile() {}

	 //! Energy is expected to be in MeV
	 double operator() (double costheta, double energy) const;

	 //! Return the average intensity over the theta range,  costhmin < costhmax
	 double average (double costhmin, double costhmax, double energy) const;

	 double avgDist () const { return m_avgDist; }


protected:

	 //! The index to this vector is i + j*ntheta where i is the energy index and
	 //! j is the theta index.  ntheta is the size of the theta array.
	 std::vector<double> m_profile;

	 //! These vectors should be sorted, m_costheta descending and m_energies
	 //! ascending
	 std::vector<double> m_costheta;
	 std::vector<double> m_energies;

	 //! The average distance to the source used in the profile squared in units
	 //! of lightseconds
	 double m_avgDist;

};

} //namespace SolarSystemTools

#endif // SolarSystemTools_SolarProfile_h
