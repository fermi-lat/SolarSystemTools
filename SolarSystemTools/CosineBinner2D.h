/** @file CosineBinner.h
@brief Define the CosineBinner2D classixel 

@author G. Johannesson

$Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/SolarSystemTools/CosineBinner2D.h,v 1.2 2012/02/16 23:19:52 gudlaugu Exp $
*/

#ifndef SolarSystemTools_CosineBinner2D_h
#define SolarSystemTools_CosineBinner2D_h

#include <vector>
#include <string>
#include <cmath>
#include <iostream>

namespace SolarSystemTools {

    /** @class CosineBinner2D
        @brief manage two sets of bins in cos(theta)

        Note that it inherits from a vector of doubles, corresponding to the bins in cos(theta)

    */

class CosineBinner2D : public std::vector<double> {
public:
    CosineBinner2D();
    
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
    const double& operator()(double costheta1, double costheta2, double phi=-3*M_PI)const;

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
        double sum=0;
				const_iterator it = begin();
				const_iterator end = begin()+ s_nbins1*s_nbins2;
        for( ; it!=end; ++it){
				  if ( (*it) != 0 )
            sum += (*it)*f(costheta1(it), costheta2(it));
        }
        return sum; 

    }
    /// integral over the costheta1 range for a given costheta2 with functor accepting costheta1 as an arg. 
    template<class F>
    double operator()(const F& f, double costheta2)const
    {   
        double sum=0;
				const size_t i2 = cosine_index2(costheta2);
				const_iterator it = begin() + i2*s_nbins1;
				const_iterator end = begin() + (i2+1)*s_nbins1;
        for( ; it != end; ++it){
					if ( (*it) != 0 ) 
            sum += (*it)*f(costheta1(it));
        }
        return sum; 

    }
                           
    /// integral over the costheta1,phi range for a given costheta2 with functor accepting costheta1,phi as an arg. 
    template<class G>
    double integral(const G& f, double costheta2)const
    {   
        double sum=0;
				const size_t i2 = cosine_index2(costheta2);
				const size_t nbins = s_nbins1*s_nbins2;
				for(size_t iphi(0); iphi < nphibins(); ++iphi) {
					  const size_t phistart(iphi*nbins);
						const_iterator it = begin()+s_nbins1*s_nbins2+i2*s_nbins1+phistart;
						const_iterator end = begin()+s_nbins1*s_nbins2+(i2+1)*s_nbins1+phistart;
            for( ; it != end; ++it){
							if ( (*it) != 0 )
                sum += (*it)*f.integral(costheta1(it), phi(it));
						}
        }
        return sum; 

    }
    /// integral over the range with functor accepting costheta, phi as args. 
    template<class G>
    double integral(const G& f)const
    {   
        double sum=0;
				const_iterator it = begin() + s_nbins1*s_nbins2;
				const_iterator en = end();
        for( ; it!=en; ++it){
					if ( (*it) != 0 )
            sum += (*it)*f.integral(costheta1(it), costheta2(it), phi(it) );
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

    static double s_cosmin1, s_cosmin2; ///< minimum value of cos(theta)
    static size_t s_nbins1, s_nbins2;  ///< number of costheta bins
    static bool  s_sqrt_weight; ///< true to use sqrt function, otherwise linear
    static size_t s_phibins; ///< number of phi bins
};

}
#endif
