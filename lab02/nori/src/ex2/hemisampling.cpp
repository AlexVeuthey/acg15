/*******************************************************************************
 *  hemisampling.cpp
 *******************************************************************************
 *  Copyright (c) 2013 Alexandre Kaspar <alexandre.kaspar@a3.epfl.ch>
 *  For Advanced Computer Graphics, at the LGG / EPFL
 * 
 *        DO NOT REDISTRIBUTE
 ***********************************************/

#include <nori/vector.h>

NORI_NAMESPACE_BEGIN

inline float length_squared(const Vector3f &vec){
   return vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2];
}

inline float length(const Vector3f &vec){
   return sqrt(length_squared(vec));
}

/**
 * \brief Cosine hemisphere sampling
 * 
 * \param sample a 2d uniform sample
 */
inline Vector3f squareToCosineHemisphere(const Point2f &sample) {
   
   // TODO implement cosine hemisphere sampling using a given uniform [0;1]^2 sample
   
   /*Vector3f v;
   do{
      v[0] = 1-2*drand48();
      v[1] = 1-2*drand48();
      v[2] = 1-2*drand48();
   } while(length_squared(v) > 1 || v[2] < 0);
   
   v /= length(v);
   
   return Vector3f(v[0], v[1], v[2]); // this is wrong, replace it with your solution!*/
   
   float x = (sample[0]-0.5)*2;
   float y = (sample[1]-0.5)*2;
   
   if(sqrt(x*x+y*y) > 1){
      return Vector3f(0, 0, 0);
   }
   else{
      float z = sqrt(1-x*x-y*y);
      
      return Vector3f(x, y, z);
   }
}

NORI_NAMESPACE_END
