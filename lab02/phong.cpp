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
#define GROUP_NUMBER 8
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

GROUP_NAMESPACE_BEGIN()

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

        // Compute perfect specular reflection
        Vector3f perfectSpecular = this->reflect(bRec.wi);

        // Compute cos(alpha) (angle between perfect specular reflection and current outgoing light ray (incoming tracer ray)
        float alphaCos = perfectSpecular[0] * bRec.wo[0] + perfectSpecular[1] * bRec.wo[1] + perfectSpecular[2] * bRec.wo[2];
        if (alphaCos < 0)
            alphaCos = 0;

        // Evaluate Phong BRDF
        return this->m_Kd * (1/M_PI) + this->m_Ks * ((this->m_exp + 2) / M_PI) * pow(alphaCos, this->m_exp);
    }

    /// Compute the density of \ref sample() wrt. solid angles

    float pdf(const BSDFQueryRecord &bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
                || Frame::cosTheta(bRec.wi) <= 0
                || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

        // Compute perfect specular reflection
        Vector3f perfectSpecular = this->reflect(bRec.wi);

        // Compute cos(alpha) (angle between perfect specular reflection and current outgoing light ray (incoming tracer ray)
        float alphaCos = perfectSpecular[0] * bRec.wo[0] + perfectSpecular[1] * bRec.wo[1] + perfectSpecular[2] * bRec.wo[2];
        if (alphaCos < 0)
            alphaCos = 0;

        // Compute ks
        float ksum = this->m_Kd.getLuminance() + this->m_Ks.getLuminance();
        float ks = this->m_Ks.getLuminance() / ksum;
        float kd = this->m_Kd.getLuminance() / ksum;
        float exp = this->m_exp;

        // Combined PDF is unbiased (see report)
        return kd * ((1.0f/M_PI) * Frame::cosTheta(bRec.wo)) + ks * (((exp + 1)/(float)(2*M_PI)) * pow(alphaCos, exp));
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

        // Generate a specular sample
        if (useSpecular)
        {
            // Compute perfect specular direction
            Vector3f perfectSpecular = reflect(bRec.wi);

            // Construct frame
            Frame frameSpecular(perfectSpecular);

            // Create random vector
            float alpha = acos(pow(sample.x(),1/(this->m_exp + 1)));
            float phi = 2 * M_PI * sample.y();
            Vector3f randomVector(sin(alpha) * cos(phi), sin(alpha) * sin(phi), cos(alpha));

            // Transform random vector to coordinates of the frame
            bRec.wo = frameSpecular.toWorld(randomVector);
        }

        // Generate a diffuse sample
        else
        {
            float theta = acos(sqrt(sample.x()));
            float phi = 2 * M_PI * sample.y();
            bRec.wo = Vector3f(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
        }

        // Return 0 if generated backside out vector
        if (Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);

        // Return evaluation of the sample
        return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
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
