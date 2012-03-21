/** @file CosineBinner.h
@brief Define the CosineBinner2D classixel 

@author G. Johannesson

$Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/SolarSystemTools/CosineBinner2D.h,v 1.4 2012/02/29 11:28:21 gudlaugu Exp $
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

    /// the binning function: add value to the selected bin, if costheta and theta in range
    void fill(double costheta, double theta, double value);
    
    /// the binning function: add value to the selected bin, if costheta and theta in range
    //! alternate version for phi binning.
    //! @param costheta cos(theta1)
    //! @param theta cos(theta2)
    //! @param phi value of phi, radians
    void fill(double costheta, double theta, double phi, double value);

    /// reference to the contents of the bin containing the cos(theta) value
    //! version that has phi as well (operator() allows multiple args)
    //! @param costheta cos(theta1)
    //! @param theta cos(theta2)
    double& operator()(double costheta, double theta, double phi=-3*M_PI);
    double operator()(double costheta, double theta, double phi=-3*M_PI)const;

    /// Provide access through the real index
    double& operator[](size_t i);
    double operator[](size_t i) const;

    size_t size() const {return std::vector<double>::size();}
    iterator begin() {return std::vector<double>::begin();}
    const_iterator begin() const {return std::vector<double>::begin();}
    iterator end() {return std::vector<double>::end();}
    const_iterator end() const {return std::vector<double>::end();}

		CosineBinner2D & operator += (const CosineBinner2D &other);

    /// cos(theta1) for the iterator
    double costheta(const const_iterator &i)const;

    /// cos(theta2) for the iterator
    double theta(const const_iterator &i) const;

    /// phi for the iterator, returns negative numbers for iterators not in the phi range.
    double phi(const const_iterator &i) const;

    /// True index from an iterator
    size_t index(const const_iterator &i) const;

    /// integral over the range with functor accepting costheta and theta as an arg. 
    template<class F>
    double operator()(const F& f)const
    {   
       double sum=0;
       for ( size_t i = 0; i < m_ithetatoi.size(); ++i ) {
          const_iterator it = begin() + i*s_nbins*(s_phibins+1);
          const_iterator endit = it + s_nbins;
          for( ; it!=endit; ++it){
             if ( (*it) != 0 )
                sum += (*it)*f(costheta(it), theta(it));
          }
       }
       return sum; 

    }
    /// integral over the costheta range for a given theta with functor accepting costheta as an arg. 
    template<class F>
    double operator()(const F& f, double theta)const
    {   
       double sum=0;
			 const size_t itheta = theta_index(theta);
       std::vector<std::pair<size_t,size_t> >::const_iterator iit = ithetatoi(itheta);
       if ( iit != m_ithetatoi.end() && iit->first == itheta ) {
          const size_t i2 = iit->second;
          const_iterator it = begin() + i2*s_nbins*(s_phibins+1);
          const_iterator endit = it + s_nbins;
          for( ; it != endit; ++it){
             if ( (*it) != 0 ) 
                sum += (*it)*f(costheta(it));
          }
       }
       return sum; 

    }
                           
    /// integral over the costheta,phi range for a given theta with functor accepting costheta,phi as an arg. 
    template<class G>
    double integral(const G& f, double theta)const
    {   
       double sum=0;
			 const size_t itheta = theta_index(theta);
       std::vector<std::pair<size_t,size_t> >::const_iterator iit = ithetatoi(itheta);
       if ( iit != m_ithetatoi.end() && iit->first == itheta ) {
          const size_t i2 = iit->second;
          const_iterator it = begin()+i2*s_nbins*(s_phibins+1)+s_nbins;
          const_iterator endit = it+s_nbins*s_phibins;
          for( ; it != endit; ++it){
             if ( (*it) != 0 )
                sum += (*it)*f.integral(costheta(it), phi(it));
          }
       }
       return sum; 

    }
    /// integral over the range with functor accepting costheta, phi as args. 
    template<class G>
    double integral(const G& f)const
    {   
       double sum=0;
       for ( size_t i = 0; i < m_ithetatoi.size(); ++i ) {
          const_iterator it = begin() + i*s_nbins*(s_phibins+1)+s_nbins;
          const_iterator endit = it + s_nbins*s_phibins;
          for( ; it!=endit; ++it){
             if ( (*it) != 0 )
                sum += (*it)*f(costheta(it), theta(it), phi(it));
          }
       }
       return sum; 

    }


    /// define the binning scheme with class (static) variables
    static void setBinning(double cosmin=0., double thmax=-1., size_t nbins=40, size_t nthbins=80, bool sqrt_weight=true);

    //! set phibins to enable binning in phi
    //! phibins [15] note default is 3 degrees per bin, since folded about 0-45 degrees.
    static void setPhiBins(size_t phibins=15);

    static std::string thetaBinning();
    static double cosmin();
    static double thmax();
    static size_t nbins();
    static size_t nthbins();
    static size_t nphibins();

    //! access to fundamental binning
    static size_t cosine_index(double costheta);
    static size_t theta_index(double theta);
    static size_t phi_index(double phi);
    static size_t index(double costheta, double theta, double phi);
		static void cosThetaBins(std::vector<double> &costheta);
		static void thetaBins(std::vector<double> &theta);
		static void phiBins(std::vector<double> &phi);

    //! the translation from 0-2pi to 0-pi/4
    static double folded_phi(double phi);
private:

		//Store sparse values of theta
    std::vector<std::pair<size_t,size_t> > m_ithetatoi;
    std::vector<std::pair<size_t,size_t> > m_itoitheta;

		static bool less_than(const std::pair<size_t,size_t> &a, const std::pair<size_t,size_t> &b) { return a.first < b.first; }
		std::vector<std::pair<size_t,size_t> >::const_iterator ithetatoi(size_t itheta) const;
		std::vector<std::pair<size_t,size_t> >::iterator ithetatoi(size_t itheta);
		std::vector<std::pair<size_t,size_t> >::const_iterator itoitheta(size_t i) const;

    size_t itheta(const CosineBinner2D::const_iterator &i) const;
    size_t aindex(double costheta, double theta, double phi) const;

    static double s_cosmin, s_thmax; ///< minimum value of cos(theta)
    static size_t s_nbins, s_nthbins;  ///< number of costheta bins
    static bool  s_sqrt_weight; ///< true to use sqrt function, otherwise linear
    static size_t s_phibins; ///< number of phi bins
};

}
#endif
