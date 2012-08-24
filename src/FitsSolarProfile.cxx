/**
 * @file FitsSolarProfile.cxx
 * @brief Access to solar intensity as a function of angle and energy from a
 * FITS file.
 * @author G. Johannesson
 *
 * There should be 3 tables in the FITS file
 * *Angles with a column called Angle giving the angle from the moving source
 * in degrees
 * *Energies with a column called Energy giving the energy bins of the profile
 * in MeV
 * *SST_PROFILE with a vector column called Profile where each row gives the
 * angular distribution at a certain energy bin in units of photons/cm2/sr/s.  This 
 * table should have a header keyword
 * AVGDIST which gives the average distance to the source used for the profile
 * in units of lightseconds.
 * 
 * $Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/src/FitsSolarProfile.cxx,v 1.3 2012/06/10 13:15:52 gudlaugu Exp $
 */

#include <string>
#include <cmath>

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "SolarSystemTools/FitsSolarProfile.h"

namespace SolarSystemTools {

	FitsSolarProfile::FitsSolarProfile(const std::string &filename) 
	{

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

		table = tip::IFileSvc::instance().readTable(filename, "SST_PROFILE");

		m_profile.resize(m_costheta.size()*m_energies.size());

		tip::Table::ConstIterator intit = table->begin();
		size_t j(0);
		for ( ; intit != table->end(); ++intit, ++j) {

			std::vector<double> profile;
			(*intit)["Profile"].get(profile);
			
			for ( size_t i = 0; i < m_costheta.size(); ++i) {
				const size_t index = j*m_costheta.size() + i;
			  m_profile[index] = profile[i];
			}
		}

		const tip::Header& hdr = table->getHeader();

		double avgdist;
		hdr["AVGDIST"].get(avgdist);
		m_avgDist = avgdist;

		delete table;
	}

} //namespace SolarSystemTools


