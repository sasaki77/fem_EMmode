#ifndef MAIN
#define MAIN

#include <GLUT/glut.h>
#include "common.h"
#include "gl2ps.h"

Params data;

void display(void)
{
  glClear(GL_COLOR_BUFFER_BIT);
  
  // 三角形一次要素の描画
  if(data.form == 1){
    for(int i=0;i<data.elem.size();i++){
      Element e = data.elem[i];
      double h[3];
      glBegin(GL_TRIANGLES);
      for(int k=0;k<3;k++) h[k] = e.p[k]->u;
      glBegin(GL_TRIANGLES);
      for(int k=0;k<3;k++){
	double red = (h[k] - 0.5)*2;
	double green = 1 - fabs(h[k]-0.5)*2;
	double blue = 1 - (h[k])*2;
	if( red<0 ) red = 0;
	if( green<0 ) green = 0;
	if( blue<0 ) blue = 0;
	glColor3d(red,green,blue);
	glVertex2d(e.p[k]->z,e.p[k]->r);
      }
      glEnd();
    }
  }  

  // 三角形二次要素の描画
  else if(data.form == 2){
    for(int i=0;i<data.elem.size();i++){
      Element e = data.elem[i];
      Node *t[3];
      double h[3];
      for(int j=0;j<4;j++){
	switch(j){
	case 0:t[0] = e.p[0]; t[1] = e.p[3]; t[2] = e.p[5]; break;
	case 1:t[0] = e.p[1]; t[1] = e.p[4]; t[2] = e.p[3]; break;
	case 2:t[0] = e.p[3]; t[1] = e.p[4]; t[2] = e.p[5]; break;
	case 3:t[0] = e.p[2]; t[1] = e.p[5]; t[2] = e.p[4]; break;
	}
	for(int k=0;k<3;k++) h[k] = t[k]->u;
	glBegin(GL_TRIANGLES);
	for(int k=0;k<3;k++){
	  double red = (h[k] - 0.5)*2;
	  double green = 1 - fabs(h[k]-0.5)*2;
	  double blue = 1 - (h[k])*2;
	  if( red<0 ) red = 0;
	  if( green<0 ) green = 0;
	  if( blue<0 ) blue = 0;
	  glColor3d(red,green,blue);
	  glVertex2d(t[k]->z,t[k]->r);
	}
	glEnd();
      }
    }
  }

  // 四角形一次要素の描画
  else if(data.form == 3){
    glBegin(GL_QUADS);
    for(int i=0;i<data.elem.size();i++){
      Element e = data.elem[i];
      glBegin(GL_QUADS);
      for(int j=0;j<4;j++){
	double red = (e.p[j]->u - 0.5)*2;
	double green = 1 - fabs(e.p[j]->u-0.5)*2;
	double blue = 1 - (e.p[j]->u)*2;
	if( red<0 ) red = 0;
	if( green<0 ) green = 0;
	if( blue<0 ) blue = 0;
	glColor3d(red,green,blue);
	glVertex2d(e.p[j]->z,e.p[j]->r);
      }
      glEnd();
    }
  }

  // 四角形要素の描画(４つの三角形と１つの四角形で描画)
  else if(data.form == 4){
    for(int i=0;i<data.elem.size();i++){
      Element e = data.elem[i];
    
      // -----四角形の描画------
      glBegin(GL_QUADS);
      for(int j=4;j<8;j++){
	double red = (e.p[j]->u - 0.5)*2;
	double green = 1 - fabs(e.p[j]->u-0.5)*2;
	double blue = 1 - (e.p[j]->u)*2;
	if( red<0 ) red = 0;
	if( green<0 ) green = 0;
	if( blue<0 ) blue = 0;
	glColor3d(red,green,blue);
	glVertex2d(e.p[j]->z,e.p[j]->r);
      }
      glEnd();
    
      // 4つの三角形の描画
      Node *t[3];
      for(int j=0;j<4;j++){
	switch(j){
	case 0:t[0] = e.p[0]; t[1] = e.p[4]; t[2] = e.p[7]; break;
	case 1:t[0] = e.p[1]; t[1] = e.p[5]; t[2] = e.p[4]; break;
	case 2:t[0] = e.p[2]; t[1] = e.p[6]; t[2] = e.p[5]; break;
	case 3:t[0] = e.p[3]; t[1] = e.p[7]; t[2] = e.p[6]; break;
	}
	glBegin(GL_TRIANGLES);
	for(int k=0;k<3;k++){
	  double red = (t[k]->u - 0.5)*2;
	  double green = 1 - fabs(t[k]->u-0.5)*2;
	  double blue = 1 - (t[k]->u)*2;
	  if( red<0 ) red = 0;
	  if( green<0 ) green = 0;
	  if( blue<0 ) blue = 0;
	  glColor3d(red,green,blue);
	  glVertex2d(t[k]->z,t[k]->r);
	}
	glEnd();
      }
    }
  }

  // 要素間の線
  // 三角形要素
  if(data.form == 1 || data.form == 2){
    glLineWidth(1.0f);
    glColor3d(0,0,0);
    for(int i=0;i<data.elem.size();i++){
      Element e = data.elem[i];
      glBegin(GL_LINE_LOOP);
      glVertex2d(e.p[0]->z,e.p[0]->r);
      glVertex2d(e.p[1]->z,e.p[1]->r);
      glVertex2d(e.p[2]->z,e.p[2]->r);
      glEnd();
    }
  }

  // 四角形要素
  else if(data.form == 3 || data.form == 4){
    glLineWidth(1.0f);
    glBegin(GL_LINES);
    glColor3d(0,0,0);
    for(int i=0;i<data.elem.size();i++){
      Element e = data.elem[i];
      for(int j=0;j<4;j++){
	glVertex2d( e.p[j]->z , e.p[j]->r );
	glVertex2d( e.p[(j+1)%4]->z , e.p[(j+1)%4]->r );
      }
    }
    glEnd();
  }

  glColor3d(0.7,0.0,0.7);
  glColor3d(0,0,0);
  // 電場のベクトル1
  if( data.isDispEV1 ){
    glLineWidth(2.0f);

    for(int i=0;i<data.node.size();i++){
      Node p = data.node[i];
      double dz = data.node[i].Ez/5;
      double dr = data.node[i].Er/5;
      double a = sqrt(dz*dz+dr*dr)/3;
      double theta = atan2(dr,dz), phi = M_PI/12;

      glBegin(GL_LINES);
      // ベクトルの長さ部分
      glVertex2d(p.z,p.r);
      glVertex2d(p.z+dz,p.r+dr);
      // ベクトルの矢印部分
      glVertex2d(p.z+dz,p.r+dr);
      glVertex2d(p.z+dz+a*cos(theta + M_PI - phi),p.r+dr+a*sin(theta + M_PI - phi));
      glVertex2d(p.z+dz,p.r+dr);
      glVertex2d(p.z+dz+a*cos(theta + M_PI + phi),p.r+dr+a*sin(theta + M_PI + phi));

      glEnd();
    }
  }

  // 電場のベクトル2
  //   if( data.isDispEV2 ){
  //     glLineWidth(2.0f);
  //     for(int i=0;i<data.ev.size();i++){
  //       Node p = data.ev[i];
  //       double dz = data.ev[i].Ez/data.ev_dev/2;
  //       double dr = data.ev[i].Er/data.ev_dev/2;
  //       double a = sqrt(dz*dz+dr*dr)/3;
  //       double theta = atan2(dr,dz), phi = M_PI/12;

  //       glBegin(GL_LINES);
  //       // ベクトルの長さ部分
  //       glVertex2d(p.z,p.r);
  //       glVertex2d(p.z+dz,p.r+dr);
  //       // ベクトルの矢印部分
  //       glVertex2d(p.z+dz,p.r+dr);
  //       glVertex2d(p.z+dz+a*cos(theta + M_PI - phi),p.r+dr+a*sin(theta + M_PI - phi));
  //       glVertex2d(p.z+dz,p.r+dr);
  //       glVertex2d(p.z+dz+a*cos(theta + M_PI + phi),p.r+dr+a*sin(theta + M_PI + phi));

  //       glEnd();
  //     }
  //   }


  // 電界強度の同じ点を線で結ぶ
  //   if( data.isDispEL ){
  //     glLineWidth(1.0f);
  //     glColor3d(0.0,0.7,0.7);
  //     double de = 0.05;
  //     for(int i=0;i<data.elem.size();i++){
  //       Element elem = data.elem[i];
  //       for(int t_itr=0;t_itr<4;t_itr++){
  // 	Node p[5];
  // 	double e[3];
  // 	switch(t_itr){
  // 	case 0:p[0] = *(elem.p[0]); p[1] = *(elem.p[3]); p[2] = *(elem.p[5]); break;
  // 	case 1:p[0] = *(elem.p[1]); p[1] = *(elem.p[4]); p[2] = *(elem.p[3]); break;
  // 	case 2:p[0] = *(elem.p[3]); p[1] = *(elem.p[4]); p[2] = *(elem.p[5]); break;
  // 	case 3:p[0] = *(elem.p[2]); p[1] = *(elem.p[5]); p[2] = *(elem.p[4]); break;
  // 	}
  // 	for(int j=0;j<3;j++)
  // 	  e[j] = sqrt( p[j].Ez * p[j].Ez + p[j].Er * p[j].Er ) / data.max_e;
	
  // 	int max_id,mid_id,min_id;
  // 	max_id = (e[2]>e[1])?2:1;
  // 	max_id = (e[max_id]>e[0])?max_id:0;
  // 	min_id = (e[2]<e[1])?2:1;
  // 	min_id = (e[min_id]<e[0])?min_id:0;
  // 	if(e[2]>e[1]) mid_id = (e[2]<e[0])?2:(e[1]>e[0])?1:0;
  // 	else          mid_id = (e[1]<e[0])?1:(e[2]>e[0])?2:0;
  // 	//cout << e[max_id] << " " << e[mid_id] << " " << e[min_id] << endl;

  // 	double dk = 0.05;
  // 	for(double k=0;k<1.0;k+=dk){
  // 	  //ここで等高線毎にひく
  // 	  if( e[max_id] > k && k > e[mid_id] ){
  // 	    //このkの値ならば要素内に線を引くことができる．
  // 	    //cout << "k = " << k << endl;
  // 	    p[3].z = (k-e[mid_id])/(e[max_id]-e[mid_id])*(p[max_id].z-p[mid_id].z) + p[mid_id].z;
  // 	    p[3].r = (k-e[mid_id])/(e[max_id]-e[mid_id])*(p[max_id].r-p[mid_id].r) + p[mid_id].r;
  // 	    p[4].z = (k-e[min_id])/(e[max_id]-e[min_id])*(p[max_id].z-p[min_id].z) + p[min_id].z;
  // 	    p[4].r = (k-e[min_id])/(e[max_id]-e[min_id])*(p[max_id].r-p[min_id].r) + p[min_id].r;
  // 	    glBegin( GL_LINES );
  // 	    glVertex2d( p[3].z, p[3].r );
  // 	    glVertex2d( p[4].z, p[4].r );
  // 	    glEnd();
  // 	  }else if( e[mid_id] > k && k > e[min_id] ){
  // 	    //このkの値ならば要素内に線を引くことができる．
  // 	    p[3].z = (k-e[min_id])/(e[max_id]-e[min_id])*(p[max_id].z-p[min_id].z) + p[min_id].z;
  // 	    p[3].r = (k-e[min_id])/(e[max_id]-e[min_id])*(p[max_id].r-p[min_id].r) + p[min_id].r;
  // 	    p[4].z = (k-e[min_id])/(e[mid_id]-e[min_id])*(p[mid_id].z-p[min_id].z) + p[min_id].z;
  // 	    p[4].r = (k-e[min_id])/(e[mid_id]-e[min_id])*(p[mid_id].r-p[min_id].r) + p[min_id].r;
  // 	    glBegin( GL_LINES );
  // 	    glVertex2d( p[3].z, p[3].r );
  // 	    glVertex2d( p[4].z, p[4].r );
  // 	    glEnd();	
  // 	  }
  // 	}
  //       }
  //     }
  //   }

  glFlush();
}

  void init(string fname)
  {
    glClearColor(0.8, 0.8, 0.8, 1.0);
    data.input(fname);
    //data.getElectroVector();
    data.normalize();

    //   cout << data.node.size() << endl;
    //   for(int i=0;i<data.node.size();i++)
    //     cout << data.node[i] << endl;
    //   cout << data.elem.size() << endl;
    //   for(int i=0;i<data.elem.size();i++)
    //     cout << data.elem[i] << endl;
  }

  void saveAsEps(){
    FILE *fp;
    char file[]="visualized.eps";
    int state = GL2PS_OVERFLOW, buffsize = 0;
    GLint viewport[4];
    int options=GL2PS_BEST_ROOT | GL2PS_OCCLUSION_CULL | GL2PS_USE_CURRENT_VIEWPORT;
    if( (fp = fopen(file, "wb")) == NULL ){
      cout << "cant open \"" << file << "\"" << endl;
      exit(0);
    }

    cout << "start to output \"" << file << "\" ..." << endl; 

    while(state == GL2PS_OVERFLOW){
      buffsize += 1024*1024;
      gl2psBeginPage(file, "gl2psTest", viewport, GL2PS_EPS, GL2PS_NO_SORT, options,
		     GL_RGBA, 0, NULL, 0, 0, 0,
		     buffsize, fp, file);
      display();
      state = gl2psEndPage();
    }
    fclose(fp);

    cout << "Done." << endl;
  }

  void keyboard( unsigned char key, int x, int y ){
    switch(key){
    case 'q': exit(0);
    case 's': saveAsEps(); break;
    case 'v':
      data.isDispEV1 = !data.isDispEV1;
      glutPostRedisplay();
      break;
    case 'b':
      data.isDispEV2 = !data.isDispEV2;
      glutPostRedisplay();
      break;
    case 'e':
      data.isDispEL = !data.isDispEL;
      glutPostRedisplay();
      break;
    default: break;
    }
  }

  void glMain(int argc,char *argv[])
  {
    string fname;
    if( argc == 1 ){
      cout << "input filename?\n>>" << endl;
      cin >> fname;
    }else if( argc == 2 ){
      fname = argv[1];
    }else{
      cout << "check number of arguments" << endl;
      return;
    }
    glutInitWindowSize(600,600);
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA);
    glutCreateWindow("fem visualizer");
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    init(fname);
    glutMainLoop();
  }

  int main(int argc, char *argv[])
  {
    srand((unsigned)time(NULL));
    glMain(argc,argv);
    return 0;
  }

  std::ostream& operator<<(std::ostream& os,const Node& obj)
  {
    return ( os << "(" << obj.z << " , " << obj.r << " ) "  << obj.id );
  } 
  std::ostream& operator<<(std::ostream& os,const Element& obj)
  {
    return ( os << "( " << obj.p[0]->id << " , " << obj.p[1]->id << " , " << obj.p[2]->id << " , " << obj.p[3]->id << " , " << obj.p[4]->id << " , " << obj.p[5]->id << " )" );
  }


#endif
