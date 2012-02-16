/** @file CosineBinner.h
@brief Define the CosineBinner2D classixel 

@author G. Johannesson

$Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/SolarSystemTools/CosineBinner2D.h,v 1.1.1.1 2012/02/11 02:26:40 gudlaugu Exp $
*/

#ifndef SolarSystemTools_CosineBinner2D_h
#define SolarSystemTools_CosineBinner2D_h

#include <tr1/unordered_map>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

namespace SolarSystemTools {

    /** @class CosineBinner2D
        @brief manage two sets of bins in cos(theta)

        Note that it inherits from a vector of doubles, corresponding to the bins in cos(theta)

    */

class CosineBinner2D : std::tr1::unordered_map<size_t,double> {
public:
    CosineBinner2D();
    
		typedef std::tr1::unordered_map<size_t,double>::const_iterator const_iterator;
		typedef std::tr1::unordered_map<size_t,double>::iterator iterator;
    
    /// the binning function: add value to the selected bin, if costheta1 and costheta2 in range
    void fill(double costheta1, double costheta2, double value);
    
    /// the binning function: add value to the selected bin, if costheta1 and costheta2 in range
    //! alternate version for phi binning.
    //! @param costheta1 cos(theta1)
    //! @param costheta2 cos(theta2)
    //! @param phi value of phi, radians
    //! @return the index of the bin that was selected
    size_t fill(double costheta1, double costheta2, double phi, double value);

    /// reference to the contents of the bin containing the cos(theta) value
    //! version that has phi as well (operator() allows multiple args)
    //! @param costheta1 cos(theta1)
    //! @param costheta2 cos(theta2)
    double& operator()(double costheta1, double costheta2, double phi=-3*M_PI);
    double operator()(double costheta1, double costheta2, double phi=-3*M_PI)const;

		double& operator[](size_t i);
		double operator[](size_t i) const;

		size_t size() const {return std::tr1::unordered_map<size_t,double>::size();}
		iterator begin() {return std::tr1::unordered_map<size_t,double>::begin();}
		const_iterator begin() const {return std::tr1::unordered_map<size_t,double>::begin();}
		iterator end() {return std::tr1::unordered_map<size_t,double>::end();}
		const_iterator end() const {return std::tr1::unordered_map<size_t,double>::end();}

    /// cos(theta1) for the iterator
    double costheta1(const const_iterator &i)const;

    /// cos(theta2) for the iterator
    double costheta2(const const_iterator &i) const;

    /// phi for the iterator, returns negative numbers for iterators not in the phi range.
    double phi(const const_iterator &i) const;

	//	const_iterator end_costh()const{return lower_bound(s_nbins1*s_nbins2);}

    /// integral over the range with functor accepting costheta1 and costheta2 as an arg. 
    template<class F>
    double operator()(const F& f)const
    {   
			  sort_index();
        double sum=0;
				std::vector<size_t>::const_iterator it = m_sortedIndex.begin();
				std::vector<size_t>::const_iterator end = std::lower_bound(m_sortedIndex.begin(), m_sortedIndex.end(), s_nbins1*s_nbins2);
        for( ; it!=end; ++it){
					  const const_iterator mit = find(*it);
            sum += (*mit).second*f(costheta1(mit), costheta2(mit));
        }
        return sum; 

    }
    /// integral over the costheta1 range for a given costheta2 with functor accepting costheta1 as an arg. 
    template<class F>
    double operator()(const F& f, double costheta2)const
    {   
			  sort_index();
        double sum=0;
				const size_t i2 = cosine_index2(costheta2);
				std::vector<size_t>::const_iterator it = std::lower_bound(m_sortedIndex.begin(), m_sortedIndex.end(), i2*s_nbins1);
				std::vector<size_t>::const_iterator end = std::lower_bound(m_sortedIndex.begin(), m_sortedIndex.end(), (i2+1)*s_nbins1);
        for( ; it != end; ++it){
					  const const_iterator mit = find(*it);
            sum += (*mit).second*f(costheta1(mit));
        }
        return sum; 

    }
                           
    /// integral over the costheta1,phi range for a given costheta2 with functor accepting costheta1,phi as an arg. 
    template<class G>
    double integral(const G& f, double costheta2)const
    {   
			  sort_index();
        double sum=0;
				const size_t i2 = cosine_index2(costheta2);
				const size_t nbins = s_nbins1*s_nbins2;
				for(size_t iphi(0); iphi < nphibins(); ++iphi) {
					  const size_t phistart(iphi*nbins);
						std::vector<size_t>::const_iterator it = std::lower_bound(m_sortedIndex.begin(), m_sortedIndex.end(), s_nbins1*s_nbins2+i2*s_nbins1+phistart);
						std::vector<size_t>::const_iterator end = std::lower_bound(m_sortedIndex.begin(), m_sortedIndex.end(), s_nbins1*s_nbins2+(i2+1)*s_nbins1+phistart);
            for( ; it != end; ++it){
							  const const_iterator mit = find(*it);
                sum += (*mit).second*f.integral(costheta1(mit), phi(mit));
						}
        }
        return sum; 

    }
    /// integral over the range with functor accepting costheta, phi as args. 
    template<class G>
    double integral(const G& f)const
    {   
			  sort_index();
        double sum=0;
				std::vector<size_t>::const_iterator it = std::lower_bound(m_sortedIndex.begin(), m_sortedIndex.end(), s_nbins1*s_nbins2);
				std::vector<size_t>::const_iterator end = m_sortedIndex.end();
        for( ; it!=end; ++it){
					  const const_iterator mit = find(*it);
            sum += (*mit)*f.integral(costheta1(mit), costheta2(mit), phi(mit) );
        }
        return sum; 

    }


    /// define the binning scheme with class (static) variables
    static void setBinning(double cosmin1=0., double cosmin2=-1., size_t nbins1=40, size_t nbins2=80, bool sqrt_weight=true);

    //! set phibins to enable binning in phi
    //! phibins [15] note default is 3 degrees per bin, since folded about 0-45 degrees.
    static void setPhiBins(size_t phibins=15);

    static std::string thetaBinning();
    static double cosmin1();
    static double cosmin2();
    static size_t nbins1();
    static size_t nbins2();
    static size_t nphibins();

    //! access to fundamental binning
    static size_t cosine_index1(double costheta1);
    static size_t cosine_index2(double costheta2);
    static size_t phi_index(double phi);
		static size_t index(double costheta1, double costheta2, double phi);

    //! the translation from 0-2pi to 0-pi/4
    static double folded_phi(double phi);
private:

		mutable std::vector<size_t> m_sortedIndex;
		mutable bool m_sorted;

		void sort_index() const;

    static double s_cosmin1, s_cosmin2; ///< minimum value of cos(theta)
    static size_t s_nbins1, s_nbins2;  ///< number of costheta bins
    static bool  s_sqrt_weight; ///< true to use sqrt function, otherwise linear
    static size_t s_phibins; ///< number of phi bins
};

}
#endif
