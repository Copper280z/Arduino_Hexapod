/* Hexapod
Setup a serial port controller Hexapod

*/

#include <Servo.h>
#include <BasicLinearAlgebra.h>
#include <Geometry.h>
#include <OtherRotations.h>
#include "math.h"
using namespace BLA;
using namespace Geometry;



class Hexapod{

  Matrix<6,3> B;
  Matrix<6,3> P;
  Matrix<6,3> P_Home;
  Matrix<6,3> Shift;
  Matrix<6,3> leg_vectors;
  Matrix<1,6> a;
  Matrix<1,6> Beta;
  Matrix<1,6> CosBeta;
  Matrix<1,6> SinBeta;

  Matrix<1,6> L;

  Matrix<1,3> P_origin;
  
  
  RefMatrix<1,3,Array<6,3>> leg1_vec = leg_vectors.Submatrix<1,3>(0,0);
  RefMatrix<1,3,Array<6,3>> leg2_vec = leg_vectors.Submatrix<1,3>(1,0);
  RefMatrix<1,3,Array<6,3>> leg3_vec = leg_vectors.Submatrix<1,3>(2,0);
  RefMatrix<1,3,Array<6,3>> leg4_vec = leg_vectors.Submatrix<1,3>(3,0);
  RefMatrix<1,3,Array<6,3>> leg5_vec = leg_vectors.Submatrix<1,3>(4,0);
  RefMatrix<1,3,Array<6,3>> leg6_vec = leg_vectors.Submatrix<1,3>(5,0);
  
  Matrix<1,6> d;
  float h;
    
  public:
    Servo leg1;
    Servo leg2;
    Servo leg3;
    Servo leg4;
    Servo leg5;
    Servo leg6;
    
    Matrix<1,3> rot_point;
    
    Hexapod();
    void begin();
    void write_all(int);
    void write(int, int, int, int, int, int);
    void write(Matrix<1,6>);
    void relMove(float, float, float, float, float, float);
    float calc_angle(int);
    Matrix<6,3> shift_by(Matrix<6,3>,Matrix<1,3>);
    Matrix<6,3> shift_by(Matrix<6,3>,float,float,float);
    Matrix<1,3> shift_by(Matrix<1,3>,float,float,float);
    Matrix<1,3> shift_by(Matrix<1,3>,Matrix<1,3>);
  };

Hexapod::Hexapod () {
//    B.Fill(0); //B = {6.54, 3.66, 2.95, 3.22, 7.54, 5.12, 8.98, 9.99, 1.56} is also valid

    B={32.077,-31.799, 11.95,
         43.577, -11.88, 11.95,
         11.50, 43.679, 11.95,
         -11.50, 43.679, 11.95,
         -43.577, -11.88, 11.95,
         -32.077, -31.799, 11.95};

    P={21.396,-22.169, 70.002,
         29.898, -7.446, 70.002,
         8.50, 29.615, 70.002,
         -8.5, 29.615, 70.002,
         -29.818, -7.446, 70.002,
         -21.398, -22.169, 70.002};

     P_Home = {21.396,-22.169, 70.002,
         29.898, -7.446, 70.002,
         8.50, 29.615, 70.002,
         -8.5, 29.615, 70.002,
         -29.818, -7.446, 70.002,
         -21.398, -22.169, 70.002};

     Beta = {240,60,0,180,120,300}; //{0,0,-30,120,60,240}
     Beta = Beta * DEG_TO_RAD;
     CosBeta = {cos(Beta(0,0)),
                cos(Beta(0,1)),
                cos(Beta(0,2)),
                cos(Beta(0,3)),
                cos(Beta(0,4)),
                cos(Beta(0,5))};
     SinBeta = {sin(Beta(0,0)),
                sin(Beta(0,1)),
                sin(Beta(0,2)),
                sin(Beta(0,3)),
                sin(Beta(0,4)),
                sin(Beta(0,5))};
     
     d = {65.4,64.6,65.1,63.9,64.3,63.1};
     h = 14.0;
     P_origin = {0,0,75};
     rot_point = P_origin;
     Shift.Fill(0);
};

void Hexapod::begin() {
    leg1.attach(8,600,2500);
    leg2.attach(9);
    leg3.attach(4,544,2300);
    leg4.attach(5);
    leg5.attach(6);//,600,2450);
    leg6.attach(7,500,2200);
};
void Hexapod::write_all(int angle){
  leg1.write(angle);
  leg2.write(180-angle);
  leg3.write(angle);
  leg4.write(180-angle);
  leg5.write(angle);
  leg6.write(180-angle);
};

void Hexapod::write(int a, int b, int c, int d, int e, int f){
  leg1.write(a);
  leg2.write(180-b);
  leg3.write(c);
  leg4.write(180-d);
  leg5.write(e);
  leg6.write(180-f);
};
void Hexapod::write(Matrix<1,6> A){
  leg1.write( A(0,0)+90);
  leg2.write(180- (A(0,1)+90) );
  leg3.write(A(0,2)+90);
  leg4.write(180- (A(0,3)+90));
  leg5.write(A(0,4)+90);
  leg6.write(180- (A(0,5)+90));
};
float Hexapod::calc_angle(int idx){
  float e = 2*h* leg_vectors(idx,2);
//  Serial << "l(0,2): " << leg_vectors(idx,2) << "\n";
//  Serial << "e: " << e << "\n";

  float f = 2*h*(CosBeta(0,idx)*leg_vectors(idx,0)+SinBeta(0,idx)*leg_vectors(idx,1));
//  Serial << "f: " << f << "\n";

  float g = L(0,idx)*L(0,idx) - (d(0,idx)*d(0,idx) - h*h);
//  Serial << "g: " << g << "\n";

  float ak = asin(g/sqrt(e*e+f*f))-atan2(f,e);
  ak = RAD_TO_DEG*ak;
//  Serial << "ak: " << ak << "\n";

  return ak;

};

Matrix<6,3> Hexapod::shift_by(Matrix<6,3> points,Matrix<1,3> vec){
    Shift={vec(0,0),vec(0,1), vec(0,2),
         vec(0,0),vec(0,1), vec(0,2),
         vec(0,0),vec(0,1), vec(0,2),
         vec(0,0),vec(0,1), vec(0,2),
         vec(0,0),vec(0,1), vec(0,2),
         vec(0,0),vec(0,1), vec(0,2)};

      return points+Shift; 
};

Matrix<6,3> Hexapod::shift_by(Matrix<6,3> points,float x, float y, float z){
    Shift={x,y, z,
         x,y, z,
         x,y, z,
         x,y, z,
         x,y, z,
         x,y, z};
     Serial << "Shift: " << Shift << "\n";

      return points+Shift; 
};
Matrix<1,3> Hexapod::shift_by(Matrix<1,3> points,Matrix<1,3> vec){
      return points+vec; 
};

Matrix<1,3> Hexapod::shift_by(Matrix<1,3> points,float x, float y, float z){
    Matrix<1,3> vec = {x,y,z};
    return points+vec; 
};
void Hexapod::relMove(float x, float y, float z, float u, float v , float w){
      Serial << "P: " << P<<"\n";
      EulerAngles rot(u*DEG_TO_RAD,v*DEG_TO_RAD,w*DEG_TO_RAD);//, EulerAngles::RotationFrame::Static, EulerAngles::RotationOrder::XYZ);
       
      P=shift_by(P,x,y,z); 
      rot_point = shift_by(rot_point,x,y,z);
      P = shift_by(shift_by(P,-rot_point) * rot.to_rotation_matrix(),rot_point) ;
      
      leg_vectors = P-B;
      Serial << "Leg vectors: " << leg_vectors<<"\n";
      
      L = {Norm(leg1_vec),
      Norm(leg2_vec),
      Norm(leg3_vec),
      Norm(leg4_vec),
      Norm(leg5_vec),
      Norm(leg6_vec)};
      Serial << "L vector: " << L <<"\n";

      a(0,0) = calc_angle(0);
      a(0,1) = calc_angle(1);
      a(0,2) = calc_angle(2);
      a(0,3) = calc_angle(3);
      a(0,4) = calc_angle(4);
      a(0,5) = calc_angle(5);
      Serial << "a: " << a << "\n\n";
      write(-a);
  
};

int pos = 90;    // variable to store the servo position
int i=90;
int j=90;
Hexapod hex;


void setup() {
//  Serial.begin(115200);
  hex.begin();
  hex.rot_point(0,2) +=20; 
//  hex.write_all(pos);
  
//  hex.write(90,90,90,90,90,90);
  hex.relMove(0,0,0,0,0,0);
  delay(500);
  hex.relMove(0,0,0,0,-10,0);
  delay(1000);
}

void loop() {

//  delay(2000);  


    hex.relMove(0,0,0,0,20,0);
    delay(1000);
    hex.relMove(0,0,0,0,-20,0);
    delay(1000);

//  for (pos = 0; pos <= 5; pos += 1) { // goes from 0 degrees to 180 degrees
//    // in steps of 1 degree
//    hex.relMove(0,0,0,1,0,0);
//    delay(50);                       // waits 15ms for the servo to reach the position
//  }
//  for (pos = 10; pos >= 0; pos -= 1) { // goes from 180 degrees to 0 degrees
//    hex.relMove(0,0,0,-1,0,0);
//    delay(50);                       // waits 15ms for the servo to reach the position
//  }
//  for (pos = 0; pos <= 5; pos += 1) { // goes from 0 degrees to 180 degrees
//    // in steps of 1 degree
//    hex.relMove(0,0,0,1,0,0);
//    delay(50);                       // waits 15ms for the servo to reach the position
//  }
}
