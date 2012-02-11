/** @file CosineBinner2D.cxx
@brief Define the CosineBinner2D class 

@author G. Johannesson

$Header: $
*/


#include "SolarSystemTools/CosineBinner2D.h"
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
    return static_cast<size_t>(folded_phi(phi)*s_phibins);
}
size_t CosineBinner2D::index(double costheta1, double costheta2, double phi)
{
	const size_t i1 = cosine_index1(costheta1);
	const size_t i2 = cosine_index2(costheta2);
	const size_t ci(i1 + i2*s_nbins1);
	if (phi < -2*M_PI) return ci;
	const size_t n = phi_index(phi);
	return ci + (n+1)*s_nbins1*s_nbins2;
}

//----------------------- Instance functions ------------------------------
CosineBinner2D::CosineBinner2D()
{
    //resize(s_phibins==0? s_nbins1*s_nbins2: s_nbins1*s_nbins2*(1+s_phibins) );
}

/// the binning function: add value to the selected bin
void CosineBinner2D::fill(double costheta1, double costheta2, double value)
{
    if( costheta1<=s_cosmin1 || costheta2<=s_cosmin2) return;
		(*this)(costheta1,costheta2) += value;
}
/// the binning function: add value to the selected bin
size_t CosineBinner2D::fill(double costheta1, double costheta2, double phi, double value)
{
    if( costheta1<=s_cosmin1 || costheta2<=s_cosmin2) return std::string::npos;
    (*this)(costheta1,costheta2) += value;
    // get the phi index and index into the array
    const size_t pi = index(costheta1,costheta2,phi);
    (*this)[pi] += value;
    return pi;
}

double& CosineBinner2D::operator()(double costheta1, double costheta2, double phi)
{
	const size_t i=index(costheta1,costheta2,phi);
  return (*this)[i];
}

const double& CosineBinner2D::operator()(double costheta1, double costheta2, double phi)const
{
	const size_t i=index(costheta1,costheta2,phi);
  return (*this)[i];
}


/// cos(theta1) for the iterator
double CosineBinner2D::costheta1(std::map<size_t,double>::const_iterator i)const
{
    const size_t bin = ( (*i).first % ( s_nbins1*s_nbins2 ) ) % s_nbins1 ;
    double f = (bin+0.5)/s_nbins1;
    if( s_sqrt_weight) f=f*f;
    return 1. - f*(1-s_cosmin1); 
}
/// cos(theta2) for the iterator
double CosineBinner2D::costheta2(std::map<size_t,double>::const_iterator i)const
{
    const size_t bin = ( (*i).first % ( s_nbins1*s_nbins2 ) ) / s_nbins1 ;
    double f = (bin+0.5)/s_nbins2;
    if( s_sqrt_weight) f=f*f;
    return 1. - f*(1-s_cosmin2); 
}
/// phi for the iterator: note only 0-pi/4 (0-45 degrees)
double CosineBinner2D::phi(std::map<size_t,double>::const_iterator i)const
{
    const size_t bin = (*i).first / ( s_nbins1*s_nbins2 ) - 1;
    double f = (bin+0.5)/s_phibins;
    if( s_sqrt_weight) f=f*f;
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

