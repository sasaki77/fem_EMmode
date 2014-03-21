#ifndef QUAD_2_ELEMENT
#define QUAD_2_ELEMENT

#include "common.h"

void Quad1Element::setParams(){
  Node vnode[4];
  vnode[0].setNode(p[0]->z - p[1]->z ,p[0]->r - p[1]->r);
  vnode[1].setNode(p[2]->z - p[1]->z ,p[2]->r - p[1]->r);
  vnode[2].setNode(p[3]->z - p[2]->z ,p[3]->r - p[2]->r);
  vnode[3].setNode(p[1]->z - p[2]->z ,p[1]->r - p[2]->r);
  
  area = fabs((vnode[0].z * vnode[1].r - vnode[0].r * vnode[1].z))/2 
    + fabs((vnode[2].z * vnode[3].r - vnode[2].r * vnode[3].z))/2; 

  double sum=0;
  for(int i=0;i<4;i++) sum = Node::len(*p[i],*p[(i+1)%4]);
  length = sum/4.0;
}

double Quad1Element::N(int i,double xi,double eta)
{
  switch(i){
  case 1: return 1.0/4.0*(1.0+xi)*(1.0+eta);
  case 2: return 1.0/4.0*(1.0-xi)*(1.0+eta);
  case 3: return 1.0/4.0*(1.0-xi)*(1.0-eta);
  case 4: return 1.0/4.0*(1.0+xi)*(1.0-eta);
  }
  return -1;
}

double Quad1Element::N_xi(int i,double xi,double eta)
{
  switch(i){
  case 1: return 1.0/4.0*(1.0+eta);
  case 2: return -1.0/4.0*(1.0+eta);
  case 3: return -1.0/4.0*(1.0-eta);
  case 4: return 1.0/4.0*(1.0-eta);
  }
  return -1;
}

double Quad1Element::N_eta(int i,double xi,double eta)
{
  switch(i){
  case 1: return 1.0/4.0*(1.0+xi);
  case 2: return 1.0/4.0*(1.0-xi);
  case 3: return -1.0/4.0*(1.0-xi);
  case 4: return -1.0/4.0*(1.0+xi);
  }
  return -1;
}

double Quad1Element::z_xi(double xi,double eta)
{
  double sum = 0;
  for(int i=0;i<6;i++)
    sum += p[i]->z*N_xi(i+1,xi,eta);
  return sum;
}
double Quad1Element::z_eta(double xi,double eta)
{
  double sum = 0;
  for(int i=0;i<6;i++)
    sum += p[i]->z*N_eta(i+1,xi,eta);
  return sum;
}
double Quad1Element::r_xi(double xi,double eta)
{
  double sum = 0;
  for(int i=0;i<6;i++)
    sum += p[i]->r*N_xi(i+1,xi,eta);
  return sum;
}
double Quad1Element::r_eta(double xi,double eta)
{
  double sum = 0;
  for(int i=0;i<6;i++)
    sum += p[i]->r*N_eta(i+1,xi,eta);
  return sum;
}

double Quad1Element::jcb(double xi,double eta)
{
  return (1.0/8.0) * (((p[0]->z - p[1]->z) * (p[2]->r - p[3]->r)- (p[2]->z - p[3]->z) *
		       (p[0]->r - p[1]->r)) * xi + 
		      ((p[0]->z - p[3]->z) * (p[1]->r - p[2]->r) - (p[1]->z - p[2]->z) * 
		       (p[0]->r - p[3]->r)) * eta
		      + ((p[0]->z - p[2]->z) * (p[1]->r - p[3]->r) - (p[1]->z - p[3]->z) *
			 (p[0]->r - p[2]->r))
		      );
}

double Quad1Element::N_z(int i,double xi,double eta)
{
  switch(i){
  case 1:
    return 1.0/(8.0*jcb(xi,eta)) * ((p[2]->r - p[3]->r) * xi + 
			       (p[1]->r - p[2]->r) * eta 
			       + p[1]->r - p[3]->r);
  case 2:
    return 1.0/(8.0*jcb(xi,eta)) * (-(p[2]->r - p[3]->r) * xi -
			       (p[0]->r - p[3]->r) * eta 
			       - p[0]->r + p[2]->r);
  case 3:
    return  1.0/(8.0*jcb(xi,eta)) * (-(p[0]->r - p[1]->r) * xi +
				(p[0]->r - p[3]->r) * eta 
				- p[1]->r + p[3]->r);
  case 4:
    return 1.0/(8.0*jcb(xi,eta)) * ((p[0]->r - p[1]->r) * xi -
				(p[1]->r - p[2]->r) * eta 
				+ p[0]->r - p[2]->r);
  }
}
double Quad1Element::N_r(int i,double xi,double eta)
{
  switch(i){
  case 1:
    return 1.0/(8.0*jcb(xi,eta)) * (-(p[2]->z - p[3]->z) * xi - 
			       (p[1]->z - p[2]->z) * eta 
			       - p[1]->z + p[3]->z);
  case 2:
    return 1.0/(8.0*jcb(xi,eta)) * ((p[2]->z - p[3]->z) * xi +
			       (p[0]->z - p[3]->z) * eta 
			       + p[0]->z - p[2]->z);
  case 3:
    return  1.0/(8.0*jcb(xi,eta)) * ((p[0]->z - p[1]->z) * xi -
				(p[0]->z - p[3]->z) * eta 
				+ p[1]->z - p[3]->z);
  case 4:
    return 1.0/(8.0*jcb(xi,eta)) * (-(p[0]->z - p[1]->z) * xi +
			       (p[1]->z - p[2]->z) * eta 
			       - p[0]->z + p[2]->z);
  }
} 

void Quad1Element::setMatrix()
{
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      double eAsum=0, eBsum=0;
      for(int k=0;k<3;k++){
	for(int l=0;l<3;l++){
	  double xi,eta,wi,wj;
	  double f,g;
	  double rr=0;

	  switch(k){
	  case 0: xi = -sqrt(3.0/5.0);  wi = 5.0/9.0; break;
	  case 1: xi = 0;               wi = 8.0/9.0; break;
	  case 2: xi = sqrt(3.0/5.0);   wi = 5.0/9.0; break;
	  }
	  switch(l){
	  case 0: eta = -sqrt(3.0/5.0); wj = 5.0/9.0; break;
	  case 1: eta = 0;              wj = 8.0/9.0; break;
	  case 2: eta = sqrt(3.0/5.0);  wj = 5.0/9.0; break;
	  }

	  for(int m=0;m<4;m++){
	    // 積分点が四角形内部にあるためrrが0になることはない
	    rr += N(m+1,xi,eta)*p[m]->r;
	  }
	  
	  f = ((N_z(j+1,xi,eta) * N_z(i+1,xi,eta) * rr)
	       + (N_r(j+1,xi,eta) * N_r(i+1,xi,eta) * rr)
	       + (N(i+1,xi,eta) * N_r(j+1,xi,eta))
	       + (N(j+1,xi,eta) * N_r(i+1,xi,eta))
	       + (N(j+1,xi,eta) * N(i+1,xi,eta) / rr)
	       ) * jcb(xi,eta);

	  g = N(j+1,xi,eta) * N(i+1,xi,eta) * rr * jcb(xi,eta);
	  
	  eAsum += f * wi * wj ;
	  eBsum += g * wi * wj;
	}
	eA[i][j] = eAsum;
	eB[i][j] = eBsum;
      }
    }
  }
}

void Quad1Element::calcEfield(vector<double> *u,double omega){
  for(int j=0;j<4;j++){
    // 一点だけが軸上であるような四角形の軸上での計算を無視
    if( p[j]->r < CALC_EPS &&
	p[(j+1)%4]->r > CALC_EPS &&
	p[(j+2)%4]->r > CALC_EPS &&
	p[(j+3)%4]->r > CALC_EPS )
      continue;

    //int id = p[j]->id;
    p[j]->Ez = p[j]->Er = 0;

    // ここから下をElementの関数で
    double xi,eta;
    switch(j+1){
    case 1: xi = 1.0;  eta = 1.0;  break;
    case 2: xi = -1.0; eta = 1.0;  break;
    case 3: xi = -1.0; eta = -1.0; break;
    case 4: xi = 1.0;  eta = -1.0; break;
    }
    for(int k=0;k<4;k++)
      p[j]->Er -= (*u)[p[k]->id]*N_z(k+1,xi,eta) / (EPS_0*omega);
	
    // Ezの第一項を計算
    if(fabs(p[j]->r)>CALC_EPS){
      for(int k=0;k<4;k++)
	// 電場計算する節点がr>0の場合
	p[j]->Ez += (*u)[p[k]->id]*N(k+1,xi,eta)/p[j]->r/(EPS_0*omega);
    }else{
      int vnum;
      if(p[(j+1)%4]->r < CALC_EPS){
	// 電場計算する節点の次の節点がr=0の場合
	vnum = (j+3)%4;
      }else{
	// 電場計算する節点の前の節点がr=0の場合
	vnum = (j+1)%4;
      }
      p[j]->Ez = (*u)[p[vnum]->id] / p[vnum]->r /(EPS_0*omega);
    }
    for(int k=0;k<4;k++)
      // Ezの第二項を計算
      p[j]->Ez += (*u)[p[k]->id]*N_r(k+1,xi,eta)/(EPS_0*omega);
  }
}


double Quad1Element::calcPsum(vector<double> *u){
  double sum = 0;
  for(int k=0;k<3;k++){
    for(int l=0;l<3;l++){
      double xi,eta,wi,wj;
      switch(k){
      case 0: xi = -sqrt(3.0/5.0);  wi = 5.0/9.0; break;
      case 1: xi = 0;               wi = 8.0/9.0; break;
      case 2: xi = sqrt(3.0/5.0);   wi = 5.0/9.0; break;
      }
      switch(l){
      case 0: eta = -sqrt(3.0/5.0); wj = 5.0/9.0; break;
      case 1: eta = 0;              wj = 8.0/9.0; break;
      case 2: eta = sqrt(3.0/5.0);  wj = 5.0/9.0; break;
      }
      double NH=0, Nr=0;
      for(int j=0;j<4;j++){
	NH += N(j+1,xi,eta) * (*u)[p[j]->id];
	Nr += N(j+1,xi,eta) * p[j]->r;
      }
      sum += NH * NH * Nr * jcb(xi,eta) * wi * wj;
    }
  }
  return sum;
}

#endif
