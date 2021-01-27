#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <random>
static std::default_random_engine engine; // random seed=10
static std::uniform_real_distribution<double> uniform(1.,0);

int count;
float total, inBox;

// user defined function below
float f (float x){
  return std::pow(cos(x),10);
}
//

//function to calculate a definite integral given bounds of integration (xmin/max) & bounds of function (ymin/ymax)
float integral (float (*f)(float), float xmin, float xmax, float ymin, float ymax){
  for (count=0; count < 1000000; count++){
    float u1 = uniform(engine);
    float u2 = uniform(engine);

    float xcoord = ((xmax - xmin)*u1) + xmin;
    float ycoord = ((ymax - ymin)*u2) + ymin;
    float val = f(xcoord);

    total++;

    if (val > ycoord){
      inBox++;
    }
  }

  float density = inBox/total;

  std::cout<<(xmax - xmin)*(ymax - ymin)*density<<std::endl;
}

int main(){
  std::cout<< "RESULT: " <<std::endl;
  integral(f,-M_PI/2,M_PI/2,0,1);
}