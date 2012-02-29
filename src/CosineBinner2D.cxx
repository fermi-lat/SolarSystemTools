/** @file CosineBinner2D.cxx
@brief Define the CosineBinner2D class 

@author G. Johannesson

$Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/src/CosineBinner2D.cxx,v 1.3 2012/02/17 01:49:51 gudlaugu Exp $
*/


#include "SolarSystemTools/CosineBinner2D.h"
#include <map>
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace{
    const double pi4( M_PI/4.0), pi2(M_PI/2.);

}
using namespace SolarSystemTools;

double CosineBinner2D::s_cosmin1=0.0;
double CosineBinner2D::s_cosmin2=-1.0;
size_t CosineBinner2D::s_nbins1=40;
size_t CosineBinner2D::s_nbins2=80;
bool   CosineBinner2D::s_sqrt_weight=true;
size_t CosineBinner2D::s_phibins=0;

size_t CosineBinner2D::nphibins(){return s_phibins;}
void   CosineBinner2D::setPhiBins(size_t n){s_phibins=n;}
double CosineBinner2D::cosmin1() { return s_cosmin1;}
double CosineBinner2D::cosmin2() { return s_cosmin2;}
size_t CosineBinner2D::nbins1() {return s_nbins1;}
size_t CosineBinner2D::nbins2() {return s_nbins2;}

size_t CosineBinner2D::cosine_index1(double costheta)
{
    double f ( (1.-costheta)/(1-s_cosmin1) );
    if(s_sqrt_weight) f=sqrt(f);
    const size_t i = static_cast<size_t>(f*s_nbins1);
    return i<s_nbins1-1 ? i : s_nbins1-1;
}
size_t CosineBinner2D::cosine_index2(double costheta)
{
    double f ( (1.-costheta)/(1-s_cosmin2) );
    if(s_sqrt_weight) f=sqrt(f);
    const size_t i = static_cast<size_t>(f*s_nbins2);
    return i<s_nbins2-1 ? i : s_nbins2-1;
}
double CosineBinner2D::folded_phi(double phi)
{
    if( phi<0) phi+= 2*M_PI;
    return  1.- fabs(pi4-fmod(phi, pi2) )/pi4 ;
}
size_t CosineBinner2D::phi_index(double phi)
{
   if ( phi < -2*M_PI ) return 0;
   return static_cast<size_t>(folded_phi(phi)*s_phibins+1);
}
size_t CosineBinner2D::index(double costheta1, double costheta2, double phi)
{
	const size_t i1 = cosine_index1(costheta1);
	const size_t i2 = cosine_index2(costheta2);
	const size_t n = phi_index(phi);
	return (i2*(s_phibins+1) + n) * s_nbins1 + i1;
}

//----------------------- Instance functions ------------------------------
CosineBinner2D::CosineBinner2D() 
{
    //resize(s_phibins==0? s_nbins1*s_nbins2: s_nbins1*s_nbins2*(1+s_phibins) );
}

std::vector<std::pair<size_t,size_t> >::iterator CosineBinner2D::icostheta2toi(size_t icostheta2){
	return std::lower_bound(m_icostheta2toi.begin(), m_icostheta2toi.end(), std::pair<size_t,size_t>(icostheta2,size_t(0)), less_than);
}

std::vector<std::pair<size_t,size_t> >::const_iterator CosineBinner2D::icostheta2toi(size_t icostheta2) const{
	return std::lower_bound(m_icostheta2toi.begin(), m_icostheta2toi.end(), std::pair<size_t,size_t>(icostheta2,size_t(0)), less_than);
}

std::vector<std::pair<size_t,size_t> >::const_iterator CosineBinner2D::itoicostheta2(size_t i) const{
	return std::lower_bound(m_itoicostheta2.begin(), m_itoicostheta2.end(), std::pair<size_t,size_t>(i,size_t(0)), less_than);
}

/// the binning function: add value to the selected bin
void CosineBinner2D::fill(double costheta1, double costheta2, double value)
{
    if( costheta1<=s_cosmin1 || costheta2<=s_cosmin2) return;
    (*this)(costheta1,costheta2) += value;
}
/// the binning function: add value to the selected bin
void CosineBinner2D::fill(double costheta1, double costheta2, double phi, double value)
{
    if( costheta1<=s_cosmin1 || costheta2<=s_cosmin2) return;
    (*this)(costheta1,costheta2) += value;
    // get the phi index and index into the array
    (*this)(costheta1,costheta2,phi) += value;
}

double& CosineBinner2D::operator()(double costheta1, double costheta2, double phi)
{
	const size_t i = index(costheta1, costheta2, phi);
	(*this)[i];
}

double CosineBinner2D::operator()(double costheta1, double costheta2, double phi)const
{
	const size_t i = index(costheta1, costheta2, phi);
	(*this)[i];
}

double & CosineBinner2D::operator[] (size_t i) {
   const size_t icostheta2 = i / ( s_nbins1*(s_phibins+1) );
   std::vector<std::pair<size_t,size_t> >::iterator it2 = icostheta2toi(icostheta2);

   //Add the bin if needed
   size_t i2(0);
   if ( it2 == m_icostheta2toi.end() || it2->first != icostheta2 ) {
      //Insert the mapping between indexes
      i2 = m_icostheta2toi.size();
      m_icostheta2toi.insert(it2, std::pair<size_t,size_t>(icostheta2,i2));
      m_itoicostheta2.insert(m_itoicostheta2.end(), std::pair<size_t,size_t>(i2,icostheta2));

      //Resize the storage, ~5% at a time
			const size_t alloc(s_nbins2/20+1);
      if (i2 % alloc == 0){
         reserve( std::min( size()+alloc*s_nbins1*(s_phibins+1), s_nbins1*s_nbins2*(s_phibins+1) ) );
      }
      resize(size()+s_nbins1*(s_phibins+1), 0.0);
   } else {
      i2 = it2->second;
   }
   return std::vector<double>::operator[](i + (i2 - icostheta2)*s_nbins1*(s_phibins+1));
   //return at(i + (i2 - icostheta2)*s_nbins1*(s_phibins+1));
}

double CosineBinner2D::operator[] (size_t i) const{
   const size_t icostheta2 = i / ( s_nbins1*(s_phibins+1) );
   std::vector<std::pair<size_t,size_t> >::const_iterator it = icostheta2toi(icostheta2);

   //Return 0 if not found
   if ( it == m_icostheta2toi.end() ) return 0;

   return std::vector<double>::operator[](i + (it->second - icostheta2)*s_nbins1*(s_phibins+1));
   //return at(i + (it->second - icostheta2)*s_nbins1*(s_phibins+1));
}

size_t CosineBinner2D::index(const CosineBinner2D::const_iterator &i) const
{
	const size_t ind = i - begin();
	const size_t icostheta2array = ind / (s_nbins1*(s_phibins+1));
	const size_t icostheta2 = itoicostheta2(icostheta2array)->second;
	return ind + (icostheta2 - icostheta2array)*s_nbins1*(s_phibins+1);
}

/// cos(theta1) for the iterator
double CosineBinner2D::costheta1(const CosineBinner2D::const_iterator &i)const
{
	  const size_t ind = i - begin();
    const size_t bin = ( ind % ( s_nbins1*(s_phibins+1) ) ) % s_nbins1 ;
    double f = (bin+0.5)/s_nbins1;
    if( s_sqrt_weight) f=f*f;
    return 1. - f*(1-s_cosmin1); 
}
/// cos(theta2) for the iterator
double CosineBinner2D::costheta2(const CosineBinner2D::const_iterator &i)const
{
	  const size_t ind = i - begin();
	  const size_t icostheta2array = ind / (s_nbins1*(s_phibins+1));
	  const size_t bin =  itoicostheta2(icostheta2array)->second;
    double f = (bin+0.5)/s_nbins2;
    if( s_sqrt_weight) f=f*f;
    return 1. - f*(1-s_cosmin2); 
}
/// phi for the iterator: note only 0-pi/4 (0-45 degrees)
double CosineBinner2D::phi(const CosineBinner2D::const_iterator &i)const
{
	  const size_t ind = i - begin();
    const size_t bin = ( ind % ( s_nbins1*(s_phibins+1) ) ) / s_nbins1;
		if (bin == 0) return -3*M_PI;
    double f = (bin-0.5)/s_phibins; // minus because we have to subtract 1 from bin
    return M_PI/4. * f;
}


std::string CosineBinner2D::thetaBinning(){ 
    if( s_sqrt_weight) {
        return "SQRT(1-COSTHETA)";
    }else{
        return "COSTHETA";
    }
}

void CosineBinner2D::setBinning(double cosmin1, double cosmin2, size_t nbins1, size_t nbins2, bool sqrt_weight)
{
    s_cosmin1=cosmin1, s_cosmin2=cosmin2, s_nbins1=nbins1, s_nbins2=nbins2, s_sqrt_weight=sqrt_weight;
}

