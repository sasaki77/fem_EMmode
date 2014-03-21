#ifndef TRI_2_ELEMENT
#define TRI_2_ELEMENT

#include "common.h"

void Tri2Element::setParams(){
    area = 0.5*fabs( (p[2]->z - p[0]->z)*(p[1]->r - p[0]->r)
		 - (p[2]->r - p[0]->r)*(p[1]->z - p[0]->z));
    length = ( Node::len(*p[0],*p[1]) + Node::len(*p[1],*p[2]) + Node::len(*p[2],*p[0]))/3;
}

double Tri2Element::N(int i,double xi,double eta)
{
  double L1 = xi;
  double L2 = eta;
  double L3 = 1-xi-eta;
  switch(i){
  case 1: return L1*(2*L1 - 1);
  case 2: return L2*(2*L2 - 1);
  case 3: return L3*(2*L3 - 1);
  case 4: return 4*L1*L2;
  case 5: return 4*L2*L3;
  case 6: return 4*L3*L1;
  }
  return -1;
}

double Tri2Element::N_xi(int i,double xi,double eta)
{
  double L1 = xi;
  double L2 = eta;
  double L3 = 1-xi-eta;
  switch(i){
  case 1: return 4*L1 - 1;
  case 2: return 0;
  case 3: return -4*L3 + 1;
  case 4: return 4*L2;
  case 5: return -4*L2;
  case 6: return 4*L3 - 4*L1;
  }
  return -1;
}

double Tri2Element::N_eta(int i,double xi,double eta)
{
  double L1 = xi;
  double L2 = eta;
  double L3 = 1-xi-eta;
  switch(i){
  case 1: return 0;
  case 2: return 4*L2-1;
  case 3: return -4*L3 + 1;
  case 4: return 4*L1;
  case 5: return 4*L3 - 4*L2;
  case 6: return -4*L1;
  }
  return -1;
}

double Tri2Element::z_xi(double xi,double eta)
{
  double sum = 0;
  for(int i=0;i<6;i++)
    sum += p[i]->z*N_xi(i+1,xi,eta);
  return sum;
}
double Tri2Element::z_eta(double xi,double eta)
{
  double sum = 0;
  for(int i=0;i<6;i++)
    sum += p[i]->z*N_eta(i+1,xi,eta);
  return sum;
}
double Tri2Element::r_xi(double xi,double eta)
{
  double sum = 0;
  for(int i=0;i<6;i++)
    sum += p[i]->r*N_xi(i+1,xi,eta);
  return sum;
}
double Tri2Element::r_eta(double xi,double eta)
{
  double sum = 0;
  for(int i=0;i<6;i++)
    sum += p[i]->r*N_eta(i+1,xi,eta);
  return sum;
}

double Tri2Element::jcb(double xi,double eta)
{
  return ( z_xi(xi,eta)*r_eta(xi,eta) - z_eta(xi,eta)*r_xi(xi,eta) );
}

double Tri2Element::N_z(int i,double xi,double eta)
{
  return (( r_eta(xi,eta)*N_xi(i,xi,eta) - r_xi(xi,eta)*N_eta(i,xi,eta))/jcb(xi,eta));
}
double Tri2Element::N_r(int i,double xi,double eta)
{
  return ((-z_eta(xi,eta)*N_xi(i,xi,eta) + z_xi(xi,eta)*N_eta(i,xi,eta))/jcb(xi,eta));
}

void Tri2Element::setMatrix()
{
  for(int i=0;i<6;i++){
    for(int j=0;j<6;j++){
      double eAsum=0, eBsum=0;
      for(int k=0;k<3;k++){
	double xi,eta,w;
	double f,g;
	double rr=0,r0=0;
	// 積分点数3 辺上3,要素内0 precision 2
// 	switch(k){
// 	case 0: xi = 0.5; eta = 0.5; w = 1./3; break;
// 	case 1: xi = 0;   eta = 0.5; w = 1./3; break;
// 	case 2: xi = 0.5; eta = 0;   w = 1./3; break;
// 	}

	// 積分点数3 辺上0,要素内3 precision 2
	switch(k){
	case 0: xi = 2./3; eta = 1./6; w = 1./6; break;
	case 1: xi = 1./6; eta = 2./3; w = 1./6; break;
	case 2: xi = 1./6; eta = 1./6; w = 1./6; break;
	}

	// 積分点数7 辺上6,要素内1 precision 3
// 	switch(k){
// 	case 0: xi = 1./3; eta = 1./3; w = 27./60; break;
// 	case 1: xi = 0.5;  eta = 0.5;  w = 8./60 ; break;
// 	case 2: xi = 0;    eta = 0.5;  w = 8./60 ; break;
// 	case 3: xi = 0.5;  eta = 0  ;  w = 8./60 ; break;
// 	case 4: xi = 1;    eta = 0  ;  w = 3./60 ; break;
// 	case 5: xi = 0;    eta = 1  ;  w = 3./60 ; break;
// 	case 6: xi = 0;    eta = 0  ;  w = 3./60 ; break;
// 	}

	// 積分点数6 辺上0,要素内6  precision 4
// 	switch(k){
// 	case 0: xi = 0.81684757; eta = 0.09157621; w = 0.10995174; break;
// 	case 1: xi = 0.09157621; eta = 0.81684757; w = 0.10995174; break;
// 	case 2: xi = 0.09157621; eta = 0.09157621; w = 0.10995174; break;
// 	case 3: xi = 0.44594849; eta = 0.44594849; w = 0.22338159; break;
// 	case 4: xi = 0.10810302; eta = 0.44594849; w = 0.22338159; break;
// 	case 5: xi = 0.44594849; eta = 0.10810302; w = 0.22338159; break;
// 	}

	// 積分点数7 辺上0,要素内7 precision 5
// 	double a1 = 0.05971587, b1 = 0.47014206;
// 	double a2 = 0.79742699, b2 = 0.10128651;
// 	switch(k){
// 	case 0: xi = 1./3; eta = 1./3; w = 0.225; break;
// 	case 1: xi = a1;   eta = b1;   w = 0.13239415; break;
// 	case 2: xi = b1;   eta = a1;   w = 0.13239415; break;
// 	case 3: xi = b1;   eta = b1;   w = 0.13239415; break;
// 	case 4: xi = a2;   eta = b2;   w = 0.12593918; break; 
// 	case 5: xi = b2;   eta = a2;   w = 0.12593918; break; 
// 	case 6: xi = b2;   eta = b2;   w = 0.12593918; break; 
// 	}

	// 積分点数19 辺上0,要素内19 precision 9
// 	double a1 = 0.91054097, b1 = 0.04472951, w1 = 0.02557768;
// 	double a2 = 0.62359293, b2 = 0.18820354, w2 = 0.07964774;
// 	double a3 = 0.12582082, b3 = 0.43708959, w3 = 0.07782754;
// 	double a4 = 0.02063496, b4 = 0.48968252, w4 = 0.03133470;
// 	double a5 = 0.74119860, b5 = 0.22196299, c5 = 0.03683841, w5 = 0.04328354;
// 	switch(k){
// 	case 0:  xi = 0.33333333; eta = 0.33333333; w = 0.09713580; break;

// 	case 1:  xi = a1; eta = b1; w = w1; break;
// 	case 2:  xi = b1; eta = a1; w = w1; break;
// 	case 3:  xi = b1; eta = b1; w = w1; break;

// 	case 4:  xi = a2; eta = b2; w = w2; break;
// 	case 5:  xi = b2; eta = a2; w = w2; break;
// 	case 6:  xi = b2; eta = b2; w = w2; break;

// 	case 7:  xi = a3; eta = b3; w = w3; break; 
// 	case 8:  xi = b3; eta = a3; w = w3; break;
// 	case 9:  xi = b3; eta = b3; w = w3; break;

// 	case 10: xi = a4; eta = b4; w = w4; break;
// 	case 11: xi = b4; eta = a4; w = w4; break;
// 	case 12: xi = b4; eta = b4; w = w4; break;

// 	case 13: xi = a5; eta = b5; w = w5; break;
// 	case 14: xi = a5; eta = c5; w = w5; break;
// 	case 15: xi = b5; eta = a5; w = w5; break;
// 	case 16: xi = b5; eta = c5; w = w5; break;
// 	case 17: xi = c5; eta = a5; w = w5; break;
// 	case 18: xi = c5; eta = b5; w = w5; break;
// 	}

	for(int l=0;l<6;l++){
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

void Tri2Element::calcEfield(vector<double> *u,double omega){
  for(int j=0;j<6;j++){
      if( j < 3 ){
	if( p[j]->r < CALC_EPS &&
	    p[(j+1)%3]->r > CALC_EPS &&
	    p[(j+2)%3]->r > CALC_EPS )
	  continue;
      }
      int id = p[j]->id;
      p[j]->Ez = p[j]->Er = 0;

      double xi,eta;
      switch(j+1){
      case 1: xi = 1.0; eta = 0.0; break;
      case 2: xi = 0.0; eta = 1.0; break;
      case 3: xi = 0.0; eta = 0.0; break;
      case 4: xi = 0.5; eta = 0.5; break;
      case 5: xi = 0.0; eta = 0.5; break;
      case 6: xi = 0.5; eta = 0.0; break;
      }
      for(int k=0;k<6;k++)
	p[j]->Er -= (*u)[p[k]->id]*N_z(k+1,xi,eta) / (EPS_0*omega);
	
      if(fabs(p[j]->r)>CALC_EPS){
	for(int k=0;k<6;k++)
	  p[j]->Ez += (*u)[p[k]->id]*N(k+1,xi,eta)/p[j]->r/(EPS_0*omega);
      }else{
	if( j >= 3 ){
	  // 中点を見る場合
	  int mnum = j-3;
	  p[j]->Ez =
	    (-(*u)[p[(mnum+2)%3]->id] + 2*(*u)[p[3+(mnum+1)%3]->id] + 
	     2*(*u)[p[3+(mnum+2)%3]->id]) /
	    (-p[(mnum+2)%3]->r + 2*p[3+(mnum+1)%3]->r + 
	     2*p[3+(mnum+2)%3]->r)/(EPS_0*omega);
	}else{
	  // 頂点を見る場合
	  int vnum;
	  if(p[(j+1)%3]->r < CALC_EPS){
	    vnum = (j+2)%3;
	    p[j]->Ez = (-(*u)[p[vnum]->id] + 4*(*u)[p[3+vnum]->id])
	      /(-p[vnum]->r + 4*p[3+vnum]->r)/(EPS_0*omega);
	  }else{
	    vnum = (j+1)%3;
	    p[j]->Ez = (-(*u)[p[vnum]->id] + 4*(*u)[p[3+j]->id])
	      /(-p[vnum]->r + 4*p[3+j]->r)/(EPS_0*omega);
	  }
	}
      }
      for(int k=0;k<6;k++)
	p[j]->Ez += (*u)[p[k]->id]*N_r(k+1,xi,eta)/(EPS_0*omega);
    }
}

double Tri2Element::calcPsum(vector<double>* u){
  double sum = 0;
  for(int m=0;m<3;m++){
    double xi,eta,w;
    switch(m){
    case 0: xi = 2./3; eta = 1./6; w = 1./6; break;
    case 1: xi = 1./6; eta = 2./3; w = 1./6; break;
    case 2: xi = 1./6; eta = 1./6; w = 1./6; break;
    }
    double NH=0, Nr=0;
    for(int j=0;j<6;j++){
      NH += N(j+1,xi,eta) * (*u)[p[j]->id] ;
      Nr += N(j+1,xi,eta) * p[j]->r;
    }
    sum += NH * NH * Nr * jcb(xi,eta) * w;
  }
  
  return sum;
}

#endif
