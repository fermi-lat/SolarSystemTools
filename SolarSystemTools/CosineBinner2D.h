/** @file CosineBinner.h
@brief Define the CosineBinner2D classixel 

@author G. Johannesson

$Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/SolarSystemTools/CosineBinner2D.h,v 1.3 2012/02/17 01:49:50 gudlaugu Exp $
*/

#ifndef SolarSystemTools_CosineBinner2D_h
#define SolarSystemTools_CosineBinner2D_h

#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

namespace SolarSystemTools {

    /** @class CosineBinner2D
        @brief manage two sets of bins in cos(theta)

        Note that it inherits from a vector of doubles, corresponding to the bins in cos(theta)

    */

class CosineBinner2D : std::vector<double> {
public:
    CosineBinner2D();
    
    typedef std::vector<double>::const_iterator const_iterator;
    typedef std::vector<double>::iterator iterator;

    /// the binning function: add value to the selected bin, if costheta1 and costheta2 in range
    void fill(double costheta1, double costheta2, double value);
    
    /// the binning function: add value to the selected bin, if costheta1 and costheta2 in range
    //! alternate version for phi binning.
    //! @param costheta1 cos(theta1)
    //! @param costheta2 cos(theta2)
    //! @param phi value of phi, radians
    void fill(double costheta1, double costheta2, double phi, double value);

    /// reference to the contents of the bin containing the cos(theta) value
    //! version that has phi as well (operator() allows multiple args)
    //! @param costheta1 cos(theta1)
    //! @param costheta2 cos(theta2)
    double& operator()(double costheta1, double costheta2, double phi=-3*M_PI);
    double operator()(double costheta1, double costheta2, double phi=-3*M_PI)const;

    /// Provide access through the real index
    double& operator[](size_t i);
    double operator[](size_t i) const;

    size_t size() const {return std::vector<double>::size();}
    iterator begin() {return std::vector<double>::begin();}
    const_iterator begin() const {return std::vector<double>::begin();}
    iterator end() {return std::vector<double>::end();}
    const_iterator end() const {return std::vector<double>::end();}

    /// cos(theta1) for the iterator
    double costheta1(const const_iterator &i)const;

    /// cos(theta2) for the iterator
    double costheta2(const const_iterator &i) const;

    /// phi for the iterator, returns negative numbers for iterators not in the phi range.
    double phi(const const_iterator &i) const;

    /// True index from an iterator
    size_t index(const const_iterator &i) const;

    /// integral over the range with functor accepting costheta1 and costheta2 as an arg. 
    template<class F>
    double operator()(const F& f)const
    {   
       double sum=0;
       for ( size_t i = 0; i < m_icostheta2toi.size(); ++i ) {
          const_iterator it = begin() + i*s_nbins1*(s_phibins+1);
          const_iterator endit = it + s_nbins1;
          for( ; it!=endit; ++it){
             if ( (*it) != 0 )
                sum += (*it)*f(costheta1(it), costheta2(it));
          }
       }
       return sum; 

    }
    /// integral over the costheta1 range for a given costheta2 with functor accepting costheta1 as an arg. 
    template<class F>
    double operator()(const F& f, double costheta2)const
    {   
       double sum=0;
			 const size_t icostheta2 = cosine_index2(costheta2);
       std::vector<std::pair<size_t,size_t> >::const_iterator iit = icostheta2toi(icostheta2);
       if ( iit != m_icostheta2toi.end() && iit->first == icostheta2 ) {
          const size_t i2 = iit->second;
          const_iterator it = begin() + i2*s_nbins1*(s_phibins+1);
          const_iterator endit = it + s_nbins1;
          for( ; it != endit; ++it){
             if ( (*it) != 0 ) 
                sum += (*it)*f(costheta1(it));
          }
       }
       return sum; 

    }
                           
    /// integral over the costheta1,phi range for a given costheta2 with functor accepting costheta1,phi as an arg. 
    template<class G>
    double integral(const G& f, double costheta2)const
    {   
       double sum=0;
			 const size_t icostheta2 = cosine_index2(costheta2);
       std::vector<std::pair<size_t,size_t> >::const_iterator iit = icostheta2toi(icostheta2);
       if ( iit != m_icostheta2toi.end() && iit->first == icostheta2 ) {
          const size_t i2 = iit->second;
          const_iterator it = begin()+i2*s_nbins1*(s_phibins+1)+s_nbins1;
          const_iterator endit = it+s_nbins1*s_phibins;
          for( ; it != endit; ++it){
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
       for ( size_t i = 0; i < m_icostheta2toi.size(); ++i ) {
          const_iterator it = begin() + i*s_nbins1*(s_phibins+1)+s_nbins1;
          const_iterator endit = it + s_nbins1*s_phibins;
          for( ; it!=endit; ++it){
             if ( (*it) != 0 )
                sum += (*it)*f(costheta1(it), costheta2(it), phi(it));
          }
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

		//Store sparse values of costheta2
    std::vector<std::pair<size_t,size_t> > m_icostheta2toi;
    std::vector<std::pair<size_t,size_t> > m_itoicostheta2;

		static bool less_than(const std::pair<size_t,size_t> &a, const std::pair<size_t,size_t> &b) { return a.first < b.first; }
		std::vector<std::pair<size_t,size_t> >::const_iterator icostheta2toi(size_t icostheta) const;
		std::vector<std::pair<size_t,size_t> >::iterator icostheta2toi(size_t icostheta);
		std::vector<std::pair<size_t,size_t> >::const_iterator itoicostheta2(size_t i) const;

    size_t icostheta2(const CosineBinner2D::const_iterator &i) const;
    size_t aindex(double costheta1, double costheta2, double phi) const;

    static double s_cosmin1, s_cosmin2; ///< minimum value of cos(theta)
    static size_t s_nbins1, s_nbins2;  ///< number of costheta bins
    static bool  s_sqrt_weight; ///< true to use sqrt function, otherwise linear
    static size_t s_phibins; ///< number of phi bins
};

}
#endif
