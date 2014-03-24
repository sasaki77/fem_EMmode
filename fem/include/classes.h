#ifndef CLASSES_H
#define CLASSES_H

#include "common.h"

class Node{
 public:
  double z,r;
  bool   known;
  double val;
  double Ez,Er;
  int    id;
  
  void setNode( double _z, double _r ){ z=_z;r=_r; }

  static double len(Node p1,Node p2){
    double dx = fabs(p1.z - p2.z);
    double dy = fabs(p1.r - p2.r);
    
    return sqrt( dx*dx + dy*dy );
  }
};

class Edge{
 public:
  string material;

  virtual double calcPowerLoss(vector<double>* u,double) = 0;
  virtual double N(int,double)                           = 0;
  virtual double N_xi(int,double)                        = 0;
  
  double getSurfResist(double omega){
    double sigma;
    if( material == "$Cu" ){
      sigma = 58.8e6;
    }
    return (sqrt(MU0*omega/(2*sigma)));
  }
};

class P2Edge: public Edge{
 public:
  Node *p[2];
  
  P2Edge(){
    p[0] = NULL;
    p[1] = NULL;
  }

  P2Edge( Node *_p0, Node *_p1 ){
    p[0] = _p0;
    p[1] = _p1;
  }
  void set( Node *_p0, Node *_p1 ){
    p[0] = _p0;
    p[1] = _p1;
  }

  double N(int i,double xi)
  {
    switch(i){
    case 1: return -1./2 * xi + 1./2;
    case 2: return  1./2 * xi + 1./2;
    }
    return -1;
  }

  double N_xi(int i,double xi)
  {
    switch(i){
    case 1: return -1./2;
    case 2: return  1./2;
    }
    return -1;
  }

  double calcPowerLoss(vector<double>* u,double omega){
    double f=0;
    double P_wall = 0;
    
    for(int m=0;m<2;m++){
      double xi,w;
      
      switch(m+1){
      case 1: xi = -1/sqrt(3); w = 1; break;
      case 2: xi =  1/sqrt(3); w = 1; break;
      }
      
      double NH=0, Nz=0, Nr=0, NR=0;
      
      for(int j=0;j<2;j++){
	NH += N(j+1,xi) * (*u)[p[j]->id];
	NR += N(j+1,xi) * p[j]->r;
	Nz += N_xi(j+1,xi) * p[j]->z;
	Nr += N_xi(j+1,xi) * p[j]->r;
      }
      //      cout << "NR = " << NR << endl;
      if( p[0]->z == p[1]->z && 
	  p[1]->z == p[2]->z ){
	f += (NH*NH*NR*Nr)*w;
      }else{
	f += (NR*NH*NH*sqrt(Nz*Nz+Nr*Nr))*w;
      }
      double rs = getSurfResist(omega);
      //    cout << "f = " << f <<endl;
      P_wall += M_PI*rs*f;
    }
    return P_wall;
  };

  bool operator==( const P2Edge& obj ){
    return ( (p[1] == obj.p[1] && p[0] == obj.p[0])  ||
	     (p[0] == obj.p[1] && p[1] == obj.p[0]));
  }
  bool operator!=( const P2Edge& obj ){
    return ( !(*this == obj) );
  }
};

class P3Edge: public Edge{
 public:
  Node *p[3];
  
  P3Edge(){
    p[0] = NULL;
    p[1] = NULL;
    p[2] = NULL;
  }
  
  // p[0]:始点(xi=-1) p[1]:中点(xi=0) p[2]:終点(xi=1)
  P3Edge( Node *_p0, Node *_p1, Node *_p2){
    p[0] = _p0;
    p[1] = _p1;
    p[2] = _p2;
  }
  
  void set( Node *_p0, Node *_p1, Node *_p2){
    p[0] = _p0;
    p[1] = _p1;
    p[2] = _p2;
  }

  double N(int i,double xi)
  {
    switch(i){
    case 1: return -0.5*xi + 0.5*xi*xi;
    case 2: return 1 - xi*xi;
    case 3: return  0.5*xi + 0.5*xi*xi;
    }
    return -1;
  }

  double N_xi(int i,double xi)
  {
    switch(i){
    case 1: return -0.5 + xi;
    case 2: return -2*xi;
    case 3: return  0.5 + xi;
    }
    return -1;
  }

  double calcPowerLoss(vector<double>* u,double omega){
    double f=0;
    double P_wall = 0;
    
    for(int m=0;m<3;m++){
      double xi,w;
      
      switch(m+1){
      case 1: xi = -sqrt(3./5); w = 5./9; break;
      case 2: xi =  0;          w = 8./9; break;
      case 3: xi =  sqrt(3./5); w = 5./9; break;
      }
      
      double NH=0, Nz=0, Nr=0, NR=0;
      
      for(int j=0;j<3;j++){
	NH += N(j+1,xi) * (*u)[p[j]->id];
	NR += N(j+1,xi) * p[j]->r;
	Nz += N_xi(j+1,xi) * p[j]->z;
	Nr += N_xi(j+1,xi) * p[j]->r;
      }
      
      //      cout << "NR = " << NR << endl;
      if( p[0]->z == p[1]->z && 
	  p[1]->z == p[2]->z &&
	  p[2]->z == p[0]->z ){
	f += (NH*NH*NR*Nr)*w;
      }else{
	f += (NR*NH*NH*sqrt(Nz*Nz+Nr*Nr))*w;
      }
      
      double rs = getSurfResist(omega);
      //    cout << "f = " << f <<endl;
      P_wall += M_PI*rs*f;
    }
    return P_wall;
  }

  bool operator==( const P3Edge& obj ){
    return ( p[1] == obj.p[1] &&
	     ((p[0] == obj.p[0] && p[2] == obj.p[2]) ||
	      (p[0] == obj.p[2] && p[2] == obj.p[0])) );
  }
  bool operator!=( const P3Edge& obj ){
    return ( !(*this == obj) );
  }
};

class Element{
 public:
  int    id;
  double area;
  double length;

  virtual void setParams() = 0;
  virtual void setMatrix() = 0;
  
  virtual int getVertexNum()        = 0;
  virtual int getVertexId(int)      = 0;
  virtual void setNode(int,Node*)   = 0;
  virtual double getAvalue(int,int) = 0;
  virtual double getBvalue(int,int) = 0;

  virtual void calcEfield(vector<double>*,double) = 0;
  virtual double calcPsum(vector<double>*)        = 0;

  // ===== 以下のメソッドは setMatrix() で用いる =====
 public:
  virtual double N(int,double,double)     = 0;
  virtual double N_xi(int,double,double)  = 0;
  virtual double N_eta(int,double,double) = 0;
  virtual double N_z(int,double,double)   = 0;
  virtual double N_r(int,double,double)   = 0;
  virtual double z_xi(double,double)      = 0;
  virtual double z_eta(double,double)     = 0;
  virtual double r_xi(double,double)      = 0;
  virtual double r_eta(double,double)     = 0;
  virtual double jcb(double,double)       = 0; 
};

class Tri1Element: public Element{
 public:
  Node   *p[3];
  double eA[3][3];
  double eB[3][3];  
  
  Tri1Element(){
    p[0] = NULL; p[1] = NULL; p[2] = NULL;
  }

  void   setParams();
  void   setMatrix();
  
  int    getVertexNum(){ return 3; };
  void   setNode(int i,Node* _p){ p[i]=_p; };
  int    getVertexId(int i){ return p[i]->id; };
  double getAvalue(int i,int j){ return eA[i][j];  };
  double getBvalue(int i,int j){ return eB[i][j];  };

  void   calcEfield(vector<double>*,double);
  double calcPsum(vector<double>*);

  // ===== 以下のメソッドは setMatrix() で用いる =====
 public:
  double N(int,double,double);
  double N_xi(int,double,double);
  double N_eta(int,double,double);
  double N_z(int,double,double);
  double N_r(int,double,double);
  double z_xi(double,double);
  double z_eta(double,double);
  double r_xi(double,double);
  double r_eta(double,double);
  double jcb(double,double); 
};

class Tri2Element: public Element{
 public:
  Node   *p[6];
  double eA[6][6];
  double eB[6][6];  
  
  Tri2Element(){
    p[0] = NULL; p[1] = NULL; p[2] = NULL;
    p[3] = NULL; p[4] = NULL; p[5] = NULL;
  }

  void   setParams();
  void   setMatrix();
  int    getVertexNum(){ return 6; };
  void   setNode(int i,Node* _p){ p[i]=_p; };
  int    getVertexId(int i){ return p[i]->id; };
  double getAvalue(int i,int j){ return eA[i][j];  };
  double getBvalue(int i,int j){ return eB[i][j];  };

  void   calcEfield(vector<double>*,double);
  double calcPsum(vector<double>*);

  // ===== 以下のメソッドは setMatrix() で用いる =====
 public:
  double N(int,double,double);
  double N_xi(int,double,double);
  double N_eta(int,double,double);
  double N_z(int,double,double);
  double N_r(int,double,double);
  double z_xi(double,double);
  double z_eta(double,double);
  double r_xi(double,double);
  double r_eta(double,double);
  double jcb(double,double); 
};

class Quad1Element: public Element{
 public:
  Node   *p[4];
  double eA[4][4];
  double eB[4][4];  
  
  Quad1Element(){
    p[0] = NULL; p[1] = NULL; p[2] = NULL; p[3] = NULL;
  }

  void   setParams();
  void   setMatrix();
  int    getVertexNum(){ return 4; };
  void   setNode(int i,Node* _p){ p[i]=_p; };
  int    getVertexId(int i){ return p[i]->id; };
  double getAvalue(int i,int j){ return eA[i][j];  };
  double getBvalue(int i,int j){ return eB[i][j];  };

  void calcEfield(vector<double>*,double);
  double calcPsum(vector<double>*);

  // ===== 以下のメソッドは setMatrix() で用いる =====
 public:
  double N(int,double,double);
  double N_xi(int,double,double);
  double N_eta(int,double,double);
  double N_z(int,double,double);
  double N_r(int,double,double);
  double z_xi(double,double);
  double z_eta(double,double);
  double r_xi(double,double);
  double r_eta(double,double);
  double jcb(double,double); 
};

class Quad2Element: public Element{
 public:
  Node   *p[8];
  double eA[8][8];
  double eB[8][8];  
  
  Quad2Element(){
    p[0] = NULL; p[1] = NULL; p[2] = NULL; p[3] = NULL;
    p[4] = NULL; p[5] = NULL; p[6] = NULL; p[7] = NULL;
  }

  void   setParams();
  void   setMatrix();
  int    getVertexNum(){ return 8; };
  void   setNode(int i,Node* _p){ p[i]=_p; };
  int    getVertexId(int i){ return p[i]->id; };
  double getAvalue(int i,int j){ return eA[i][j];  };
  double getBvalue(int i,int j){ return eB[i][j];  };

  void   calcEfield(vector<double>*,double);
  double calcPsum(vector<double>*);

  // ===== 以下のメソッドは setMatrix() で用いる =====
 public:
  double N(int,double,double);
  double N_xi(int,double,double);
  double N_eta(int,double,double);
  double N_z(int,double,double);
  double N_r(int,double,double);
  double z_xi(double,double);
  double z_eta(double,double);
  double r_xi(double,double);
  double r_eta(double,double);
  double jcb(double,double); 
};

class FEM{
 public:
  int form;
  vector< Element* > elem;
  vector< Node >     node;
  vector< Edge* >    mt_edge;	 // edges for material of surface of cavity
  double lambda;		  
  vector< vector< double > > A;	  
  vector< vector< double > > B;	  
  vector< double > u;		 //H_theta, u[i]はnode[i]における磁場 [A/m]
 				  
  double P_in;		// input power [j]
  double P_wall;	// power dissipation on metal wall [W]
  double omega;		// 2*pi*frequency [1/s]
  double Q;		// Q-value        [(dimensionless)]
  double E0;		// average of electric field on Z axis.
  double LenZ;		// length of Z axis.
  double MAX_H;			  
  double R;		// shunt impedance

  void input(string);
  void input_re(string); // input already calculated file
  void output(string,double);
  void outputEonRAxis( string );
  void disp();
  void memo(double);

  void setLargeMatrix();

  void cgm();
  void setbc();
  void calcMField();
  void calcEField();
  void calcPowerLoss();
  
  double getAveragedArea()
  {
    double sum=0;
    for(unsigned int i=0;i<elem.size();i++)
      sum += elem[i]->area;
    sum /= elem.size();
    return sum;
  }
  
  double getAveragedLength()
  {
    double sum = 0;
    for(unsigned int i=0;i<elem.size();i++)
      sum += elem[i]->length;
    sum /= elem.size();
    return sum;
  }
};

#endif
