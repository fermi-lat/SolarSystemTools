/**
 * @file SolarProfile.cxx
 * @brief Virtual access class to solar profiles for the solar template class
 * @author G. Johannesson
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/src/SolarProfile.cxx,v 1.1 2012/02/15 03:03:51 gudlaugu Exp $
 */

#include <vector>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <cmath>

#include "SolarSystemTools/SolarProfile.h"

namespace SolarSystemTools {

	double SolarProfile::operator() (double costheta, double energy) const {

		//Check to see that we are within range
		if (costheta > m_costheta.front() || costheta < m_costheta.back()) {
			throw std::runtime_error("Request for solar intensity at an angle that is outside the profile boundary");
			return 0;
		}
		if (energy < m_energies.front() || energy > m_energies.back()) {
			throw std::runtime_error("Request for solar intensity at an energy that is outside the profile range");
			return 0;
		}

		//Do interpolation in angle (linear in costheta) and energy (powerlaw)
		const std::vector<double>::const_iterator th_it = std::upper_bound(m_costheta.begin(), m_costheta.end(), costheta, std::greater<double>());
		const size_t i_th = (th_it == m_costheta.end()) ? m_costheta.size()-1 : th_it-m_costheta.begin();

		const std::vector<double>::const_iterator en_it = std::lower_bound(m_energies.begin(), m_energies.end(), energy);
		const size_t i_en = (en_it == m_energies.end()) ? m_energies.size()-1 : en_it-m_energies.begin();

		//The total value is frac_u*int_u + (1-frac_u)*int_l
		const double frac_u = (m_costheta[i_th-1] - costheta)/(m_costheta[i_th-1]-m_costheta[i_th]);

		const double int_enl = frac_u*m_profile[i_th + (i_en-1)*m_costheta.size()] + (1-frac_u)*m_profile[i_th-1 + (i_en-1)*m_costheta.size()];
		const double int_enu = frac_u*m_profile[i_th +  i_en   *m_costheta.size()] + (1-frac_u)*m_profile[i_th-1 +  i_en   *m_costheta.size()];

		//Do linear interpolation if either value is <= 0
		if ( int_enu > 0 && int_enl > 0 ) {
			const double ind = log(int_enu/int_enl)/log(m_energies[i_en]/m_energies[i_en-1]);
			const double out = int_enl * pow( energy/m_energies[i_en-1], ind );
			return out;
		} else {
			const double fr_l = (m_energies[i_en] - energy)/(m_energies[i_en]-m_energies[i_en-1]);
			const double out = int_enl*fr_l + (1-fr_l)*int_enu;
			return out;
		}
	}


	double SolarProfile::average(double costhmin, double costhmax, double energy) const {
		//Get the values at the boundaries first, throws an error if we are out of
		//range
		const double first = (*this)(costhmax, energy);
		const double last = (*this)(costhmin, energy);

		//Now loop over all the angles inbetween and do the integration
		std::vector<double>::const_iterator it = std::lower_bound(m_costheta.begin(), m_costheta.end(), costhmax, std::greater<double>());
		std::vector<double>::const_iterator end = std::upper_bound(m_costheta.begin(), m_costheta.end(), costhmin, std::greater<double>());
		
		double integ(0);
		double lastval(first);
		double lastcosth(costhmax);
		for ( ; it != end; ++it) {
			const double currval = (*this)(*it,energy);
			integ += 0.5*(lastval+currval)*(lastcosth - *it);
			lastval = currval;
			lastcosth = *it;
		}

		//Add the last step of the integral and divide with the total width
		integ += 0.5*(lastval + last)*(lastcosth - costhmin);
		integ /= costhmax - costhmin;

		return integ;
	}
	
}// namespace SolarSystemTools
