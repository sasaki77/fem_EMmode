#ifndef IO_CPP
#define IO_CPP

#include "common.h"

void FEM::input(string ifname)
{
  if( ifname.empty() ){
    cout << " input >> " ;
    cin  >> ifname;
  }

  string   fname_in = ifname + ".elem";
  ifstream ifs_in(fname_in.c_str());

  if( !ifs_in ){
    cout << " cannot open \"" << fname_in << "\" in INPUT()" << endl;
    exit(EXIT_FAILURE);
  }else{
    cout << " open \"" << fname_in << "\"" << endl;
  }

  double dummy;
  int    nodenum,elemnum;
  double mind=DBL_MAX,maxd=DBL_MIN;
  
  ifs_in >> form;
  ifs_in >> nodenum;
  ifs_in >> elemnum;
  node.resize( nodenum );
  elem.resize( elemnum );

  for(int i=0;i<node.size();i++){
    ifs_in >> dummy;
    node[i].id = i;
    ifs_in >> node[i].z;
    ifs_in >> node[i].r;
    if(node[i].r<CALC_EPS){
      if(mind>node[i].z) mind = node[i].z;
      if(maxd<node[i].z) maxd = node[i].z;
    }
  }
  LenZ = fabs(maxd-mind);

  int t;
  for(int i=0;i<elem.size();i++){
    switch(form){
    case 1:elem[i] = new Tri1Element();break;
    case 2:elem[i] = new Tri2Element();break;
    case 3:elem[i] = new Quad1Element();break;
    case 4:elem[i] = new Quad2Element();break;
    }
    
    ifs_in >> dummy;

    elem[i]->id = i;
    for(int j=0;j<elem[i]->getVertexNum();j++){
      ifs_in >> t;
      t--;
      elem[i]->setNode(j,&node[t]);
    }
    elem[i]->setMatrix();
    elem[i]->setParams();
  }

  // 全体節点行列の大きさを決定
  A.resize( node.size() );
  for(int i=0;i<A.size();i++) A[i].resize( node.size() );
  B.resize( node.size() );
  for(int i=0;i<B.size();i++) B[i].resize( node.size() );
  u.resize( node.size() );


  string   fname_bc = ifname + ".bc";
  ifstream ifs_bc(fname_bc.c_str());

  if( !ifs_bc ){
    cout << " cannot open \"" << fname_bc << "\"" << endl;
    exit(EXIT_FAILURE);
  }else{
    cout << " open \"" << fname_bc << "\"" << endl;
  }

  char buff[BUFF_SIZE];
  int  read_state = 0;

  while( ifs_bc.getline( buff,BUFF_SIZE  ) ){
    if( buff[0] == '\n' || buff[0] == '\r' || buff[0] == '\0') continue;
    
    if( (string)buff == "$begin_bc" ) read_state = 1;
    else if ( (string)buff == "$begin_material" ) read_state = 2;
    
    if(read_state == 0) continue;

    switch(read_state){
    case 1:
      // 境界条件の読み込み
      while( ifs_bc.getline( buff,BUFF_SIZE ) ){
	if( (string)buff == "$end" ) break;
	
	stringstream     ss;
	vector< string > vs;
	char             minibuff[BUFF_SIZE];

	ss << (string)buff;

	while( ss.getline(minibuff,BUFF_SIZE,' '))
	  if( minibuff[0] != '\0')
	    vs.push_back((string)minibuff);

	int id;
	id = atoi( vs[0].c_str()  ) - 1;
	cout << "node[" << id << "] =" << node[id] << endl;
	node[id].val = atof( vs[1].c_str()  );
	node[id].known = true;
      }
      break;

    case 2:
      // 材料特性の読み込み
      // materialはedgeの数が示されているためfor文で読み出す
      int edgenum;
      int p1,p2,p3;
      ifs_bc >> edgenum;
      mt_edge.resize(edgenum);

      for(int i=0;i<edgenum;i++){
	if(form ==1 || form ==3){
	  int P1,P2;
	  ifs_bc >> p1;
	  ifs_bc >> p2;
	  mt_edge[i] = new P2Edge(&node[p1-1],&node[p2-1]);
	}else{
	  int P1,P2,P3;
	  ifs_bc >> p1;
	  ifs_bc >> p2;
	  ifs_bc >> p3;
	  mt_edge[i] = new P3Edge(&node[p1-1],&node[p2-1],&node[p3-1]);
	}
	ifs_bc >> mt_edge[i]->material;
	read_state = 0;
      }
      break;
    }
  }
  disp();
}

void FEM::input_re(string ifname){
  string ifname_re = ifname + ".out";
  ifstream ifs(ifname_re.c_str());
  if( !ifs ){
    cout << " cannot open \"" << ifname_re << "\" in INPUT()" << endl;
    exit(EXIT_FAILURE);
  }else{
    cout << " open \"" << ifname_re << "\"" << endl;
  }
  
  // $で囲まれている始めの部分を読み飛ばす
  for(int i=0;i<2;i++){
    ifs.ignore(numeric_limits<std::streamsize>::max() , '$');
  }
  double dummy;
  ifs >> dummy; // ignore information of order
  ifs >> lambda;
  ifs >> dummy; // ignore node size
  ifs >> dummy; // ignore elem size
  for(int i=0;i<node.size();i++){
    ifs >> dummy; // ignore output number
    ifs >> dummy; // ignore node[i].z
    ifs >> dummy; // ignore node[i].r
    ifs >> u[i];
    ifs >> dummy;
    ifs >> dummy;
  }
  omega = LIGHT_SPEED* sqrt(lambda);
}

void FEM::output(string ofname,double timer)
{
  ofstream ofs((ofname+".out").c_str());
  
  // 始めに$で囲まれている部分は"femvis"で読み飛ばされる
  ofs << "$----------------------------------------" << endl;
  // 実行した日時の表示
  time_t     current;
  time(&current);                    
  struct tm  *local = localtime(&current);
  ofs << " " << asctime(local);
  ofs << " time [s]        = " << timer << endl;
  ofs << " nodenum         = " << node.size() << endl;
  ofs << " elemnum         = " << elem.size() << endl;
  ofs << " edge size       = " << getAveragedLength() << endl;
  ofs << " area size       = " << getAveragedArea() << endl;
  ofs << " lambda          = " << setprecision(16) << lambda << setprecision(6)<< endl;
  ofs << " frequency [Hz]  = " << omega/(2*M_PI) << endl;
  ofs << " input power [J] = " << P_in << endl;
  ofs << " power loss [W]  = " << P_wall << endl;
  ofs << " Q               = " << Q << endl;
  ofs << " E0 [V/m]        = " << E0 << endl;
  ofs << " R [ohm/m]       = " << R << endl;
  ofs << "----------------------------------------$" << endl;

  ofs << form << endl;
  ofs << lambda << endl;

  ofs << node.size() << endl;
  ofs << elem.size() << endl;

  for(int i=0;i<node.size();i++)
    ofs << i+1 << " " << node[i].z << " " << node[i].r << " " << u[i] << " "
	<< node[i].Ez << " " << node[i].Er << " " << endl;
  
  for(int i=0;i<elem.size();i++){
    ofs << i+1;
    for(int j=0;j<elem[i]->getVertexNum();j++){
      ofs << " " << elem[i]->getVertexId(j);
    }
    ofs << endl;
  }
}

void FEM::outputEonRAxis( string ofname )
{
  ofstream ofs((ofname + "_Ez_z.txt").c_str());
  ofs << "node num = "    << node.size() << endl;
  ofs << "element num = " << elem.size() << endl;

  for(int i=0;i<node.size();i++){
	  if( fabs(node[i].r)<CALC_EPS ){
		  ofs << node[i].z << "\t" << node[i].Ez << endl;
	  }
  }

}

void FEM::disp()
{
  cout << "  number of nodes    :" << node.size() << endl;
  cout << "  number of elements :" << elem.size() << endl;
  for(int i=0;i<node.size();i++) cout << node[i] << endl;
  for(int i=0;i<mt_edge.size();i++) cout << "  edge[" << i << "] = " << *mt_edge[i] << endl;
  for(int i=0;i<elem.size();i++){
    cout << *elem[i] << endl;
  }
//   for(int i=0;i<elem.size();i++){
//     cout << "  elem-" << elem[i]->id << " " << elem[i] << endl;
//     cout << "  eA" << endl;
//     for(int j=0;j<6;j++){
//       cout << "  | " ;
//       for(int k=0;k<6;k++)
// 	cout << setw(10) << elem[i]->eA[j][k] << " ";
//       cout << " |" << endl;
//     }
//     cout << "  eB" << endl;
//     for(int j=0;j<6;j++){
//       cout << "  | " ;
//       for(int k=0;k<6;k++)
// 	cout << setw(10) << elem[i]->eB[j][k] << " ";
//       cout << " |" << endl;
//     }
//     cout << endl;
//   }
}

void FEM::memo(double time)
{
  ofstream ofs("precision.txt",std::ios::out | std::ios::app);
  ofs << time << '\t' << setprecision(16) << lambda << endl;
}


void FEM::setLargeMatrix()
{
  for(int i=0;i<elem.size();i++){
    for(int j=0;j<elem[i]->getVertexNum();j++){
      for(int k=0;k<elem[i]->getVertexNum();k++){
	A[elem[i]->getVertexId(j)][elem[i]->getVertexId(k)] += elem[i]->getAvalue(j,k);
	B[elem[i]->getVertexId(j)][elem[i]->getVertexId(k)] += elem[i]->getBvalue(j,k);
      }
    }
  }
}

int wid(int n)
{
  int cnt=1;
  while( n>=10 ){
    n/=10;
    cnt++;
  }
  return cnt;
}

extern FEM fem;

ostream& operator<<(ostream& os,const Node& obj)
{
  return ( os << "  node[" << setw(wid(fem.node.size())) << obj.id << "] : (" << setw(10) << obj.z << " , " << setw(10) << obj.r << " ) , known="
	   << obj.known << " , value=" << obj.val );
} 

ostream& operator<<(ostream& os,/*const*/ Element& obj)
{
  Element *e = &obj;
  if(typeid(obj) == typeid(Tri1Element)){
    Tri1Element *e1 = dynamic_cast<Tri1Element*>(e);
    os << "  Tri1elem[" << setw(wid(fem.elem.size())) << obj.id << "] : ( " ;
    for(int i=0;i<3;i++) os << e1->p[i]->id << " ";
    os << ")";
  }
  else if(typeid(obj) == typeid(Tri2Element)){
    Tri2Element *e2 = dynamic_cast<Tri2Element*>(e);
    os << "  Tri2elem[" << setw(wid(fem.elem.size())) << obj.id << "] : ( " ;
    for(int i=0;i<6;i++) os << e2->p[i]->id << " ";
    os << ")";
  }
  else if(typeid(obj) == typeid(Quad1Element)){
    Quad1Element *e3 = dynamic_cast<Quad1Element*>(e);
    os << "  Quad1elem[" << setw(wid(fem.elem.size())) << obj.id << "] : ( " ;
    for(int i=0;i<4;i++) os << e3->p[i]->id << " ";
    os << ")";
  }
  else if(typeid(obj) == typeid(Quad2Element)){
    Quad2Element *e4 = dynamic_cast<Quad2Element*>(e);
    os << "  Quad2elem[" << setw(wid(fem.elem.size())) << obj.id << "] : ( " ;
    for(int i=0;i<8;i++) os << e4->p[i]->id << " ";
    os << ")";
  }
  return ( os );
}

ostream& operator<<(ostream& os,/*const*/ Edge& obj)
{
  Edge *e = &obj;
  if(typeid(obj) == typeid(P2Edge)){
    P2Edge *e2 = dynamic_cast<P2Edge*>(e);
    os << "  ( node[" << e2->p[0]->id << "], node[" << e2->p[1]->id 
       << "]) material = " << e2->material ;
  }
  else if(typeid(obj) == typeid(P3Edge)){
    P3Edge *e3 = dynamic_cast<P3Edge*>(e);
    os << "  ( node[" << e3->p[0]->id << "], node[" << e3->p[1]->id 
       << "], node[" << e3->p[2]->id << "]) material = " << e3->material ;
  }
  return os; 
}

#endif
