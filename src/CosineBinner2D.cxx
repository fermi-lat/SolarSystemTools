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
		const double ln2( log(2.) );
}
using namespace SolarSystemTools;

double CosineBinner2D::s_cosmin=0.0;
double CosineBinner2D::s_cosmin2=-1.0;
size_t CosineBinner2D::s_nbins=40;
size_t CosineBinner2D::s_nbins2=40;
bool   CosineBinner2D::s_sqrt_weight=true;
size_t CosineBinner2D::s_phibins=0;
double CosineBinner2D::s_power2=3.;

size_t CosineBinner2D::nphibins(){return s_phibins;}
void   CosineBinner2D::setPhiBins(size_t n){s_phibins=n;}
double CosineBinner2D::cosmin() { return s_cosmin;}
double CosineBinner2D::cosmin2() { return s_cosmin2;}
double CosineBinner2D::power2() { return s_power2;}
size_t CosineBinner2D::nbins() {return s_nbins;}
size_t CosineBinner2D::nbins2() {return s_nbins2;}

void CosineBinner2D::cosThetaBins(std::vector<double> &costheta)
{
	costheta.resize(s_nbins+1);
	for (size_t i = 0; i < s_nbins+1; ++i) {
		double f (static_cast<double>(i)/s_nbins);
    if(s_sqrt_weight) f*=f;
		costheta[i] = 1. - f*(1. - s_cosmin);
	}
}

void CosineBinner2D::cosTheta2Bins(std::vector<double> &costheta2)
{
	costheta2.resize(s_nbins2+1);
	for (size_t i = 0; i < s_nbins2+1; ++i) {
		double f (static_cast<double>(i)/s_nbins2);
		f = pow(f,s_power2);
		costheta2[i] = 1. - f*(1. - s_cosmin2);
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
size_t CosineBinner2D::costheta2_index(double costheta2)
{
		double f ( (1.-costheta2)/(1-s_cosmin2) );
		f = pow(f,1./s_power2);
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
size_t CosineBinner2D::index(double costheta, double costheta2, double phi)
{
	const size_t i1 = cosine_index(costheta);
	const size_t i2 = costheta2_index(costheta2);
	const size_t n = phi_index(phi);
	return (i2*(s_phibins+1) + n) * s_nbins + i1;
}

//----------------------- Instance functions ------------------------------
CosineBinner2D::CosineBinner2D() 
{
    //resize(s_phibins==0? s_nbins*s_nbins2: s_nbins*s_nbins2*(1+s_phibins) );
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
void CosineBinner2D::fill(double costheta, double costheta2, double value)
{
    if( costheta<=s_cosmin || costheta2<s_cosmin2) return;
    (*this)(costheta,costheta2) += value;
}
/// the binning function: add value to the selected bin
void CosineBinner2D::fill(double costheta, double costheta2, double phi, double value)
{
    if( costheta<=s_cosmin || costheta2<s_cosmin2) return;
    (*this)(costheta,costheta2) += value;
    // get the phi index and index into the array
    (*this)(costheta,costheta2,phi) += value;
}

double& CosineBinner2D::operator()(double costheta, double costheta2, double phi)
{
	const size_t i = index(costheta, costheta2, phi);
	(*this)[i];
}

double CosineBinner2D::operator()(double costheta, double costheta2, double phi)const
{
	const size_t i = index(costheta, costheta2, phi);
	(*this)[i];
}

double & CosineBinner2D::operator[] (size_t i) {
   const size_t icostheta2 = i / ( s_nbins*(s_phibins+1) );
   std::vector<std::pair<size_t,size_t> >::iterator it2 = icostheta2toi(icostheta2);

   //Add the bin if needed
   size_t i2(0);
   if ( it2 == m_icostheta2toi.end() || it2->first != icostheta2 ) {
      //Insert the mapping between indexes
      i2 = m_icostheta2toi.size();
      m_icostheta2toi.insert(it2, std::pair<size_t,size_t>(icostheta2,i2));
      m_itoicostheta2.insert(m_itoicostheta2.end(), std::pair<size_t,size_t>(i2,icostheta2));

      //Resize the storage, ~5% at a time
			const size_t alloc((s_nbins2)/20+1);
      if (i2 % alloc == 0){
         reserve( std::min( size()+alloc*s_nbins*(s_phibins+1), s_nbins*(s_nbins2)*(s_phibins+1) ) );
      }
      resize(size()+s_nbins*(s_phibins+1), 0.0);
   } else {
      i2 = it2->second;
   }
   return std::vector<double>::operator[](i + (i2 - icostheta2)*s_nbins*(s_phibins+1));
   //return at(i + (i2 - icostheta2)*s_nbins*(s_phibins+1));
}

double CosineBinner2D::operator[] (size_t i) const{
   const size_t icostheta2 = i / ( s_nbins*(s_phibins+1) );
   std::vector<std::pair<size_t,size_t> >::const_iterator it = icostheta2toi(icostheta2);

   //Return 0 if not found
   if ( it == m_icostheta2toi.end() || it->first != icostheta2 ) return 0;

   return std::vector<double>::operator[](i + (it->second - icostheta2)*s_nbins*(s_phibins+1));
   //return at(i + (it->second - icostheta2)*s_nbins*(s_phibins+1));
}

bool CosineBinner2D::hasCostheta2 ( double costheta ) const {
	const size_t icostheta2 = costheta2_index(costheta);
   std::vector<std::pair<size_t,size_t> >::const_iterator it = icostheta2toi(icostheta2);

   if ( it == m_icostheta2toi.end() || it->first != icostheta2 ) return false;
	 return true;
}

CosineBinner2D& CosineBinner2D::operator += (const CosineBinner2D &other) {
	//Check that their sizes match (currently have no change to check cosbin2)
	size_t sother = other.size()/other.m_icostheta2toi.size();
	size_t s = size()/m_icostheta2toi.size();
	if (sother != s)
		throw std::runtime_error("CosineBinner2D::operator+=: sizes don't match");

	//Loop over icostheta in other and add to current map
	std::vector<std::pair<size_t,size_t> >::const_iterator it = other.m_itoicostheta2.begin();
	for ( ; it != other.m_itoicostheta2.end(); ++it) {
		const size_t icostheta2 = it->second;
		const size_t iadd = icostheta2*s_nbins*(s_phibins+1);
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
	const size_t icostheta2array = ind / (s_nbins*(s_phibins+1));
	const size_t icostheta2 = itoicostheta2(icostheta2array)->second;
	return ind + (icostheta2 - icostheta2array)*s_nbins*(s_phibins+1);
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
double CosineBinner2D::costheta2(const CosineBinner2D::const_iterator &i)const
{
	  const size_t ind = i - begin();
	  const size_t icostheta2array = ind / (s_nbins*(s_phibins+1));
	  const size_t bin =  itoicostheta2(icostheta2array)->second;
		double f = (bin+0.5)/s_nbins2;
    f = pow(f,s_power2);
    return 1. - f*(1-s_cosmin2); 
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

void CosineBinner2D::setPower2(double power2) {  s_power2 = power2; }

void CosineBinner2D::setBinning(double cosmin, double cosmin2, size_t nbins, size_t nbins2, bool sqrt_weight)
{
    s_cosmin=cosmin, s_cosmin2=cosmin2, s_nbins=nbins, s_nbins2=nbins2, s_sqrt_weight=sqrt_weight;
}

