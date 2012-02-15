/**
 * @file FitsSolarProfile.cxx
 * @brief Access to solar intensity as a function of angle and energy from a
 * FITS file compatible with Andy's solar IC tool
 * @author G. Johannesson
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/SolarSystemTools/SolarSystemTools/FitsSolarProfile.h,v 1.1.1.1 2012/02/11 02:26:40 gudlaugu Exp $
 */

#include <string>
#include <cmath>

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "SolarSystemTools/FitsSolarProfile.h"

namespace SolarSystemTools {

	FitsSolarProfile::FitsSolarProfile(const std::string &filename) {

		const tip::Table * table = tip::IFileSvc::instance().readTable(filename, "ANGLES");

		tip::Table::ConstIterator angit = table->begin();
		for ( ; angit != table->end(); ++angit) {
			double theta;
			(*angit)["Angle"].get(theta);
			m_costheta.push_back(cos(theta*M_PI/180.));
		}

		delete table;

    table = tip::IFileSvc::instance().readTable(filename, "ENERGIES");

		tip::Table::ConstIterator enit = table->begin();
		for ( ; enit != table->end(); ++enit) {
			double energy;
			(*enit)["Energy"].get(energy);
			m_energies.push_back(energy);
		}

		delete table;

		table = tip::IFileSvc::instance().readTable(filename, "SPECTRUM_FOR_ANGLES");

		m_profile.resize(m_costheta.size()*m_energies.size());

		tip::Table::ConstIterator intit = table->begin();
		size_t j(0);
		for ( ; intit != table->end(); ++intit, ++j) {
			for ( size_t i = 0; i < m_costheta.size(); ++i) {
				const size_t index = j*m_costheta.size() + i;
				std::stringstream ss;
				ss << "intensity_" << i+1;
				double value;
			  (*intit)[ss.str()].get(value);
			  m_profile[index] = value;
			}
		}

		delete table;
	}

} //namespace SolarSystemTools


