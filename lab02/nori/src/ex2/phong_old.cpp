/*******************************************************************************
 *  phong.cpp
 *******************************************************************************
 *  Copyright (c) 2013 Alexandre Kaspar <alexandre.kaspar@a3.epfl.ch>
 *  For Advanced Computer Graphics, at the LGG / EPFL
 * 
 *        DO NOT REDISTRIBUTE
 ***********************************************/

#include <nori/acg.h>
#include <nori/bsdf.h>
#include <nori/frame.h>
#include "hemisampling.cpp"
#include <stdlib.h>
#include <time.h>

NORI_NAMESPACE_BEGIN

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// put your group number here!
#define GROUP_NUMBER 0
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define PI 3.141592

GROUP_NAMESPACE_BEGIN()

class uhdas{
   
};

/**
 * \brief Phong BRDF model
 */
class Phong : public BSDF {
public:

    Phong(const PropertyList &propList) {
        m_Kd = propList.getColor("kd", Color3f(0.5f));
        m_Ks = propList.getColor("ks", Color3f(0.5f));
        m_exp = propList.getFloat("n", 20.0f);
        srand(time(NULL));
    }

    /// Reflection in local coordinates
    inline Vector3f reflect(const Vector3f &wi) const {
        return Vector3f(-wi.x(), -wi.y(), wi.z());
    }

    /// Evaluate the BRDF model

    Color3f eval(const BSDFQueryRecord &bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
                || Frame::cosTheta(bRec.wi) <= 0
                || Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);

        // Based on http://www.cs.virginia.edu/~jdl/importance.doc‎
        float dot_result = bRec.wi[0]*bRec.wo[0]+bRec.wi[1]*bRec.wo[1]+bRec.wi[2]*bRec.wo[2];
        float length_wi = sqrt(bRec.wi[0]*bRec.wi[0]+bRec.wi[1]*bRec.wi[1]+bRec.wi[2]*bRec.wi[2]);
        float length_wo = sqrt(bRec.wo[0]*bRec.wo[0]+bRec.wo[1]*bRec.wo[1]+bRec.wo[2]*bRec.wo[2]);
        float alpha = acos(dot_result/(length_wi*length_wo));
        
        if(alpha > PI/2) alpha = PI/2;
              
        // TODO evaluate the phong brdf for the pair of given incoming
        //      and outgoing directions bRec.wi & bRec.wo
        
                        
        return m_Kd*1/PI + m_Ks*(m_exp+2)/(2*PI)*pow(cos(alpha), m_exp);
    }

    /// Compute the density of \ref sample() wrt. solid angles

    float pdf(const BSDFQueryRecord &bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
                || Frame::cosTheta(bRec.wi) <= 0
                || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

        // TODO implement the pdf of phong importance sampling
        float dot_result = bRec.wi[0]*bRec.wo[0]+bRec.wi[1]*bRec.wo[1]+bRec.wi[2]*bRec.wo[2];
        float length_wi = sqrt(bRec.wi[0]*bRec.wi[0]+bRec.wi[1]*bRec.wi[1]+bRec.wi[2]*bRec.wi[2]);
        float length_wo = sqrt(bRec.wo[0]*bRec.wo[0]+bRec.wo[1]*bRec.wo[1]+bRec.wo[2]*bRec.wo[2]);
        float alpha = acos(dot_result/(length_wi*length_wo));
        if(alpha > PI/2) alpha = PI/2;
        
        return (m_exp+1)/(2*PI)*pow(cos(alpha), m_exp);
    }

    /// Draw a a sample from the BRDF model

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample_) const {
        Point2f sample(sample_);
        if (Frame::cosTheta(bRec.wi) <= 0)
            return Color3f(0.0f);

        float r = float(rand())/RAND_MAX;
        // 1. Select diffuse or specular
        float kd = m_Kd.getLuminance();
        float ks = m_Ks.getLuminance();
        float specSamplingWeight = ks / (ks + kd);
        bool useSpecular = true;
        if (r > specSamplingWeight) {
            useSpecular = false;
        }

        bRec.measure = ESolidAngle;
        bRec.eta = 1.0f; // no change in relative index
        
        // TODO implement phong importance sampling
        
        //bRec.wo = Vector3f(0.0f, 0.0f, 0.0f); // this is utterly wrong!
        
        return eval(bRec) * Frame::cosTheta(bRec.wo)/pdf(bRec);
    }

    /// Return a human-readable summary

    QString toString() const {
        return QString(
                "Phong[\n"
                "  Kd = %1\n"
                "  Ks = %2\n"
                "  n  = %3\n"
                "]").arg(m_Kd.toString()).arg(m_Ks.toString()).arg(m_exp);
    }

    Color3f getColor() const {
        return m_Kd;
    }

    EClassType getClassType() const {
        return EBSDF;
    }
private:
    Color3f m_Kd, m_Ks;
    float m_exp;
};

GROUP_NAMESPACE_END

NORI_REGISTER_GROUP_CLASS(Phong, "phong");
NORI_NAMESPACE_END
