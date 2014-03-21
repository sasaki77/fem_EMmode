#ifndef TRI_1_ELEMENT
#define TRI_1_ELEMENT

#include "common.h"

void Tri1Element::setParams(){
    area = 0.5*fabs( (p[2]->z - p[0]->z)*(p[1]->r - p[0]->r)
		 - (p[2]->r - p[0]->r)*(p[1]->z - p[0]->z));
    length = ( Node::len(*p[0],*p[1]) + Node::len(*p[1],*p[2]) + Node::len(*p[2],*p[0]))/3;
}

double Tri1Element::N(int i,double xi,double eta)
{
  double L1 = xi;
  double L2 = eta;
  double L3 = 1-xi-eta;
  switch(i){
  case 1: return L1;
  case 2: return L2;
  case 3: return L3;
  }
  return -1;
}

double Tri1Element::N_xi(int i,double xi,double eta)
{
  double L1 = xi;
  double L2 = eta;
  double L3 = 1-xi-eta;
  switch(i){
  case 1: return 1;
  case 2: return 0;
  case 3: return -1;
  }
  return -1;
}

double Tri1Element::N_eta(int i,double xi,double eta)
{
  double L1 = xi;
  double L2 = eta;
  double L3 = 1-xi-eta;
  switch(i){
  case 1: return 0;
  case 2: return 1;
  case 3: return -1;
  }
  return -1;
}

double Tri1Element::z_xi(double xi,double eta)
{
  double sum = 0;
  for(int i=0;i<3;i++)
    sum += p[i]->z*N_xi(i+1,xi,eta);
  return sum;
}
double Tri1Element::z_eta(double xi,double eta)
{
  double sum = 0;
  for(int i=0;i<3;i++)
    sum += p[i]->z*N_eta(i+1,xi,eta);
  return sum;
}
double Tri1Element::r_xi(double xi,double eta)
{
  double sum = 0;
  for(int i=0;i<3;i++)
    sum += p[i]->r*N_xi(i+1,xi,eta);
  return sum;
}
double Tri1Element::r_eta(double xi,double eta)
{
  double sum = 0;
  for(int i=0;i<3;i++)
    sum += p[i]->r*N_eta(i+1,xi,eta);
  return sum;
}

double Tri1Element::jcb(double xi,double eta)
{
  return ( z_xi(xi,eta)*r_eta(xi,eta) - z_eta(xi,eta)*r_xi(xi,eta) );
}

double Tri1Element::N_z(int i,double xi,double eta)
{
  return (( r_eta(xi,eta)*N_xi(i,xi,eta) - r_xi(xi,eta)*N_eta(i,xi,eta))/jcb(xi,eta));
}
double Tri1Element::N_r(int i,double xi,double eta)
{
  return ((-z_eta(xi,eta)*N_xi(i,xi,eta) + z_xi(xi,eta)*N_eta(i,xi,eta))/jcb(xi,eta));
}

void Tri1Element::setMatrix()
{
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      double eAsum=0, eBsum=0;
      for(int k=0;k<3;k++){
	double xi,eta,w;
	double f,g;
	double rr=0,r0=0;

	// 積分点数3 辺上0,要素内3 precision 2
	switch(k){
	case 0: xi = 2./3; eta = 1./6; w = 1./6; break;
	case 1: xi = 1./6; eta = 2./3; w = 1./6; break;
	case 2: xi = 1./6; eta = 1./6; w = 1./6; break;
	}

	for(int l=0;l<3;l++){
	  rr += N(l+1,xi,eta)*p[l]->r;
	  r0 += p[l]->r;
	}
	r0 /= 6;
 
 	f =  ( N_r(i+1,xi,eta)  * N_r(j+1,xi,eta)*rr 
	       + N_z(i+1,xi,eta)  * N_z(j+1,xi,eta)*rr
	       +   N(i+1,xi,eta)  * N_r(j+1,xi,eta)
	       + N_r(i+1,xi,eta)  *   N(j+1,xi,eta)
	       +   N(i+1,xi,eta)  *   N(j+1,xi,eta)/((fabs(rr)<CALC_EPS)?r0:rr)
	       ) * jcb(xi,eta);
 	g = N(i+1,xi,eta) * N(j+1,xi,eta) * rr * jcb(xi,eta);
	eAsum += f * w ;
	eBsum += g * w;
      }
      eA[i][j] = eAsum;
      eB[i][j] = eBsum;
    }
  }
}

void Tri1Element::calcEfield(vector<double> *u,double omega){
  for(int j=0;j<3;j++){
    if( p[j]->r < CALC_EPS &&
	p[(j+1)%3]->r > CALC_EPS &&
	p[(j+2)%3]->r > CALC_EPS )
      continue;
    int id = p[j]->id;
    p[j]->Ez = p[j]->Er = 0;
      
    double xi,eta;
    switch(j+1){
    case 1: xi = 1.0; eta = 0.0; break;
    case 2: xi = 0.0; eta = 1.0; break;
    case 3: xi = 0.0; eta = 0.0; break;
    }
    for(int k=0;k<3;k++)
      p[j]->Er -= (*u)[p[k]->id]*N_z(k+1,xi,eta) / (EPS_0*omega);
	
    if(fabs(p[j]->r)>CALC_EPS){
      for(int k=0;k<3;k++)
	p[j]->Ez += (*u)[p[k]->id]*N(k+1,xi,eta)/p[j]->r/(EPS_0*omega);
    }else{
      int vnum;
      if(p[(j+1)%3]->r < CALC_EPS){
	vnum = (j+2)%3;
      }else{
	vnum = (j+1)%3;
      }
      p[j]->Ez = (*u)[p[vnum]->id] / p[vnum]->r /(EPS_0*omega);
    }
    for(int k=0;k<3;k++)
      p[j]->Ez += (*u)[p[k]->id]*N_r(k+1,xi,eta)/(EPS_0*omega);
  }
}

double Tri1Element::calcPsum(vector<double> *u){
  double sum = 0;
  for(int m=0;m<3;m++){
    double xi,eta,w;
    switch(m){
    case 0: xi = 2./3; eta = 1./6; w = 1./6; break;
    case 1: xi = 1./6; eta = 2./3; w = 1./6; break;
    case 2: xi = 1./6; eta = 1./6; w = 1./6; break;
    }
    double NH=0, Nr=0;
    for(int j=0;j<3;j++){
      NH += N(j+1,xi,eta) * (*u)[p[j]->id] ;
      Nr += N(j+1,xi,eta) * p[j]->r;
    }
    sum += NH * NH * Nr * jcb(xi,eta) * w;
  }
  return sum;
}

#endif
