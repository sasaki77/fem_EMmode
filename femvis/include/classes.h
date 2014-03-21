#ifndef CLASSES
#define CLASSES

#include "common.h"

class Node{
 public:
  double r,z;
  double u;
  double Er,Ez;
  int id;
  Node(){
    r=z=u=Er=Ez=0;
  }
  void set( double _z, double _r ){
    z = _z;
    r = _r;
  }
  void set( double _z, double _r, int _id ){
    z = _z;
    r = _r;
    id = _id;
  }
};

class Element{
 public:
  Node *p[8];

  // cf. calcEF.cpp
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

class Params
{
 public:
  int form;
  string ifname;
  int vertexNum;
  double max_u,min_u;
  double max_e,max_ev;
  double omega;

  bool isDispEV1,isDispEV2, isDispEL;

  vector< Node > node;
  vector< Element > elem;

  // cf. calcEF.cpp
  vector< Node > ev;
  double ev_dev;
  void getElectroVector();

  Params(){
    ifname = "";
    max_u = max_e = max_ev = DBL_MIN;
    min_u = DBL_MAX;
    isDispEV1 = true;
    isDispEV2 = true;
    isDispEL = true;
    ev_dev = 30;
  }

  void input(string fname){
    double dummy;
    int nodenum,elemnum;
    if(fname == ""){
      cout << "input filename?" << flush << ">>" ;
      cin >> ifname;
    }else{
      ifname = fname;
    }

    ifname += ".out" ;

    ifstream ifs(ifname.c_str());
    if( !ifs ){
      cerr << "cannot found \"" << ifname << "\"" << endl;
      exit(2);
    }
    
    // $で囲まれている始めの部分を読み飛ばす
    for(int i=0;i<2;i++)
      ifs.ignore(numeric_limits<std::streamsize>::max() , '$');

    ifs >> form;
    switch(form){
    case 1: vertexNum = 3; break;
    case 2: vertexNum = 6; break;
    case 3: vertexNum = 4; break;
    case 4: vertexNum = 8; break;
    }
    ifs >> dummy; //lambda value
    omega = LIGHT_SPEED * sqrt(dummy);

    ifs >> nodenum;
    ifs >> elemnum;
    node.resize(nodenum);
    elem.resize(elemnum);

    for(int i=0;i<node.size();i++){
      ifs >> node[i].id;
      ifs >> node[i].z;
      ifs >> node[i].r;
      ifs >> node[i].u;
      ifs >> node[i].Ez;
      ifs >> node[i].Er;
    }

    int t;
    for(int i=0;i<elem.size();i++){
      ifs >> dummy;
      for(int j=0;j<vertexNum;j++){
	ifs >> t;
	elem[i].p[j] = &node[t];
      }
    }
    //disp();
  }

  void disp(){
    cout << "node size : " << node.size() << endl;
    cout << "elem size : " << elem.size() << endl;
    cout << "max_e : " << max_e << endl; 
    for(int i=0;i<node.size();i++){
      cout << "node[" << i << "] : " << node[i] << endl;
      cout << " (Ez,Er) = (" << node[i].Ez << "," << node[i].Er << ")" << endl;
    }
    for(int i=0;i<elem.size();i++)
      cout << "elem[" << i << "] : " << elem[i] << endl;
  }

  void normalize(){
    double maxz,maxr;
    double minz,minr;
    double rz,rr,rp;
    double square;
    double erms;
    maxz = maxr = DBL_MIN;
    minz = minr = DBL_MAX;

    for(int i=0;i<node.size();i++){
      if(minz > node[i].z) minz = node[i].z;
      if(minr > node[i].r) minr = node[i].r;
      if(maxz < node[i].z) maxz = node[i].z;
      if(maxr < node[i].r) maxr = node[i].r;
      if(max_u < node[i].u) max_u = node[i].u;
      if(min_u > node[i].u) min_u = node[i].u;
      erms = sqrt(node[i].Er*node[i].Er + node[i].Ez*node[i].Ez);
      if(max_e < erms )	max_e = erms;
    }
    rz = maxz-minz;
    rr = maxr-minr;
    rp = max_u-min_u;
    square = (rz>rr)?rz:rr;
    for(int i=0;i<node.size();i++){
      node[i].z = ( (node[i].z - minz) / square )*1.8 - 0.9;
      node[i].r = ( (node[i].r - minr) / square )*1.8 - 0.9;
      node[i].Er /= max_e;
      node[i].Ez /= max_e;
      node[i].u = (node[i].u - min_u)/max_u;
    }
    for(int i=0;i<ev.size();i++){
      ev[i].z = ( (ev[i].z - minz) / square )*1.8 - 0.9;
      ev[i].r = ( (ev[i].r - minr) / square )*1.8 - 0.9;
      ev[i].Er /= max_e;
      ev[i].Ez /= max_e;
    }
  }
};

#endif
