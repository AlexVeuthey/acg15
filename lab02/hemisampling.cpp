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

/**
 * \brief Cosine hemisphere sampling
 * 
 * \param sample a 2d uniform sample
 */
inline Vector3f squareToCosineHemisphere(const Point2f &sample) {

        // TODO implement cosine hemisphere sampling using a given uniform [0;1]^2 sample
        float r = sqrt(sample.x());
        float theta = 2 * M_PI * sample.y();
        float x = r * cos(theta);
        float y = r * sin(theta);
        float z = sin(acos(sqrt(x*x + y*y)));
        return Vector3f(x,y,z);

}

NORI_NAMESPACE_END
