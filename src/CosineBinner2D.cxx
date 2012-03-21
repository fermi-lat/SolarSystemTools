/** @file CosineBinner2D.cxx
@brief Define the CosineBinner2D class 

@author G. Johannesson

$Header: /nfs/slac/g/glast/ground/cvs/SolarSystemTools/src/CosineBinner2D.cxx,v 1.4 2012/02/29 11:28:22 gudlaugu Exp $
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

double CosineBinner2D::s_cosmin=0.0;
double CosineBinner2D::s_thmax=M_PI;
size_t CosineBinner2D::s_nbins=40;
size_t CosineBinner2D::s_nthbins=40;
bool   CosineBinner2D::s_sqrt_weight=true;
size_t CosineBinner2D::s_phibins=0;

size_t CosineBinner2D::nphibins(){return s_phibins;}
void   CosineBinner2D::setPhiBins(size_t n){s_phibins=n;}
double CosineBinner2D::cosmin() { return s_cosmin;}
double CosineBinner2D::thmax() { return s_thmax;}
size_t CosineBinner2D::nbins() {return s_nbins;}
size_t CosineBinner2D::nthbins() {return s_nthbins;}

void CosineBinner2D::cosThetaBins(std::vector<double> &costheta)
{
	costheta.resize(s_nbins+1);
	for (size_t i = 0; i < s_nbins+1; ++i) {
		double f (static_cast<double>(i)/s_nbins);
    if(s_sqrt_weight) f*=f;
		costheta[i] = 1. - f*(1. - s_cosmin);
	}
}

void CosineBinner2D::thetaBins(std::vector<double> &theta)
{
	theta.resize(s_nthbins+1);
	for (size_t i = 0; i < s_nthbins+1; ++i) {
		double f (static_cast<double>(i)/s_nthbins);
    f*=f;
		theta[i] = f*s_thmax;
	}
}

void CosineBinner2D::phiBins(std::vector<double> &phi)
{
	if (s_phibins == 0) return; 
	phi.resize(s_phibins+1);
	for (size_t i = 0; i < s_phibins+1; ++i) {
		double f (static_cast<double>(i)/s_phibins);
		phi[i] = f*pi4;
	}
}

size_t CosineBinner2D::cosine_index(double costheta)
{
    double f ( (1.-costheta)/(1-s_cosmin) );
    if(s_sqrt_weight) f=sqrt(f);
    const size_t i = static_cast<size_t>(f*s_nbins);
    return i<s_nbins-1 ? i : s_nbins-1;
}
size_t CosineBinner2D::theta_index(double theta)
{
		double f ( theta/s_thmax );
		f = sqrt(f);
		const size_t i = static_cast<size_t>(f*s_nthbins);
		return i<s_nthbins-1 ? i : s_nthbins-1;
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
size_t CosineBinner2D::index(double costheta, double theta, double phi)
{
	const size_t i1 = cosine_index(costheta);
	const size_t i2 = theta_index(theta);
	const size_t n = phi_index(phi);
	return (i2*(s_phibins+1) + n) * s_nbins + i1;
}

//----------------------- Instance functions ------------------------------
CosineBinner2D::CosineBinner2D() 
{
    //resize(s_phibins==0? s_nbins*s_nthbins: s_nbins*s_nthbins*(1+s_phibins) );
}

std::vector<std::pair<size_t,size_t> >::iterator CosineBinner2D::ithetatoi(size_t itheta){
	return std::lower_bound(m_ithetatoi.begin(), m_ithetatoi.end(), std::pair<size_t,size_t>(itheta,size_t(0)), less_than);
}

std::vector<std::pair<size_t,size_t> >::const_iterator CosineBinner2D::ithetatoi(size_t itheta) const{
	return std::lower_bound(m_ithetatoi.begin(), m_ithetatoi.end(), std::pair<size_t,size_t>(itheta,size_t(0)), less_than);
}

std::vector<std::pair<size_t,size_t> >::const_iterator CosineBinner2D::itoitheta(size_t i) const{
	return std::lower_bound(m_itoitheta.begin(), m_itoitheta.end(), std::pair<size_t,size_t>(i,size_t(0)), less_than);
}

/// the binning function: add value to the selected bin
void CosineBinner2D::fill(double costheta, double theta, double value)
{
    if( costheta<=s_cosmin || theta>s_thmax) return;
    (*this)(costheta,theta) += value;
}
/// the binning function: add value to the selected bin
void CosineBinner2D::fill(double costheta, double theta, double phi, double value)
{
    if( costheta<=s_cosmin || theta>s_thmax) return;
    (*this)(costheta,theta) += value;
    // get the phi index and index into the array
    (*this)(costheta,theta,phi) += value;
}

double& CosineBinner2D::operator()(double costheta, double theta, double phi)
{
	const size_t i = index(costheta, theta, phi);
	(*this)[i];
}

double CosineBinner2D::operator()(double costheta, double theta, double phi)const
{
	const size_t i = index(costheta, theta, phi);
	(*this)[i];
}

double & CosineBinner2D::operator[] (size_t i) {
   const size_t itheta = i / ( s_nbins*(s_phibins+1) );
   std::vector<std::pair<size_t,size_t> >::iterator it2 = ithetatoi(itheta);

   //Add the bin if needed
   size_t i2(0);
   if ( it2 == m_ithetatoi.end() || it2->first != itheta ) {
      //Insert the mapping between indexes
      i2 = m_ithetatoi.size();
      m_ithetatoi.insert(it2, std::pair<size_t,size_t>(itheta,i2));
      m_itoitheta.insert(m_itoitheta.end(), std::pair<size_t,size_t>(i2,itheta));

      //Resize the storage, ~5% at a time
			const size_t alloc((s_nthbins)/20+1);
      if (i2 % alloc == 0){
         reserve( std::min( size()+alloc*s_nbins*(s_phibins+1), s_nbins*(s_nthbins)*(s_phibins+1) ) );
      }
      resize(size()+s_nbins*(s_phibins+1), 0.0);
   } else {
      i2 = it2->second;
   }
   return std::vector<double>::operator[](i + (i2 - itheta)*s_nbins*(s_phibins+1));
   //return at(i + (i2 - itheta)*s_nbins*(s_phibins+1));
}

double CosineBinner2D::operator[] (size_t i) const{
   const size_t itheta = i / ( s_nbins*(s_phibins+1) );
   std::vector<std::pair<size_t,size_t> >::const_iterator it = ithetatoi(itheta);

   //Return 0 if not found
   if ( it == m_ithetatoi.end() || it->first != itheta ) return 0;

   return std::vector<double>::operator[](i + (it->second - itheta)*s_nbins*(s_phibins+1));
   //return at(i + (it->second - itheta)*s_nbins*(s_phibins+1));
}

CosineBinner2D& CosineBinner2D::operator += (const CosineBinner2D &other) {
	//Check that their sizes match (currently have no change to check cosbin2)
	size_t sother = other.size()/other.m_ithetatoi.size();
	size_t s = size()/m_ithetatoi.size();
	if (sother != s)
		throw std::runtime_error("CosineBinner2D::operator+=: sizes don't match");

	//Loop over icostheta in other and add to current map
	std::vector<std::pair<size_t,size_t> >::const_iterator it = other.m_itoitheta.begin();
	for ( ; it != other.m_itoitheta.end(); ++it) {
		const size_t itheta = it->second;
		const size_t iadd = itheta*s_nbins*(s_phibins+1);
		const size_t otheradd = it->first*s_nbins*(s_phibins+1);
		for (size_t i = 0; i < s_nbins*(s_phibins+1); ++i) {
			(*this)[iadd+i] += other.at(otheradd+i);
		}
	}

	return *this;
	
}

size_t CosineBinner2D::index(const CosineBinner2D::const_iterator &i) const
{
	const size_t ind = i - begin();
	const size_t ithetaarray = ind / (s_nbins*(s_phibins+1));
	const size_t itheta = itoitheta(ithetaarray)->second;
	return ind + (itheta - ithetaarray)*s_nbins*(s_phibins+1);
}

/// cos(theta1) for the iterator
double CosineBinner2D::costheta(const CosineBinner2D::const_iterator &i)const
{
	  const size_t ind = i - begin();
    const size_t bin = ( ind % ( s_nbins*(s_phibins+1) ) ) % s_nbins ;
    double f = (bin+0.5)/s_nbins;
    if( s_sqrt_weight) f=f*f;
    return 1. - f*(1-s_cosmin); 
}
/// cos(theta2) for the iterator
double CosineBinner2D::theta(const CosineBinner2D::const_iterator &i)const
{
	  const size_t ind = i - begin();
	  const size_t ithetaarray = ind / (s_nbins*(s_phibins+1));
	  const size_t bin =  itoitheta(ithetaarray)->second;
		double f = (bin+0.5)/s_nthbins;
		f *= f;
		return f*s_thmax; 
}
/// phi for the iterator: note only 0-pi/4 (0-45 degrees)
double CosineBinner2D::phi(const CosineBinner2D::const_iterator &i)const
{
	  const size_t ind = i - begin();
    const size_t bin = ( ind % ( s_nbins*(s_phibins+1) ) ) / s_nbins;
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

void CosineBinner2D::setBinning(double cosmin, double thmax, size_t nbins, size_t nthbins, bool sqrt_weight)
{
    s_cosmin=cosmin, s_thmax=thmax, s_nbins=nbins, s_nthbins=nthbins, s_sqrt_weight=sqrt_weight;
}

