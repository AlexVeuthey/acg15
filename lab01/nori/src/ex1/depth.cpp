/*******************************************************************************
 *  depth.cpp
 *******************************************************************************
 *  Copyright (c) 2013 Alexandre Kaspar <alexandre.kaspar@a3.epfl.ch>
 *  For Advanced Computer Graphics, at the LGG / EPFL
 * 
 *        DO NOT REDISTRIBUTE
 ***********************************************/

#include <nori/acg.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// put your group number here!
#define GROUP_NUMBER 38736
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

GROUP_NAMESPACE_BEGIN()

class sadads{
   
};

/**
 * \brief Depth mapping: a simple rendering technique that 
 * displays the depth of an object.
 */
class Depth : public Integrator{
public:
	Depth(const PropertyList &propList) {
                /* Depth near and far plane distance */
                m_near = propList.getFloat("near", 1e-4);
                m_far = propList.getFloat("far", 1e2);
                /* Min intensity */
                m_Ka = propList.getFloat("ambiant", 0.1);
                m_gamma = propList.getFloat("gamma", 5.0);
	}

	Color3f Li(const Scene *scene, Sampler *samp, const Ray3f &ray) const {
                
		/* Find the surface that is visible in the requested direction */
		Intersection its;
		if (!scene->rayIntersect(ray, its)){
			return Color3f(0.0f);
      }
      
      
      //TODO implement the depth color integration here
      //float dist = 1.0;
      //Color3f colour;
      Point3f intersectPoint = its.p;
      Point3f origin = ray.o; //using ray for origin, not the best
      
      float dist = sqrt((intersectPoint[0]-origin[0])*(intersectPoint[0]-origin[0]) + (intersectPoint[1]-origin[1])*(intersectPoint[1]-origin[1]) + (intersectPoint[2]-origin[2])*(intersectPoint[2]-origin[2]) );
      
      //qDebug() << intersectPoint.toString();
      //qDebug() << origin.toString();
      
      //Vector3f vec;
      //its.toWorld(vec);
                      
      if(dist >= m_far) dist = 0;
      if(dist <= m_near) dist = 1.0;
      else dist = 1.0 - ( (dist-m_near) / (m_far-m_near) );
      
      
      dist = m_Ka + (1.0-m_Ka)* pow(dist, m_gamma);
      //apply m_gamma and find ka
      
		return Color3f(dist, dist, dist);
	}

	QString toString() const {
		return QString("Depth[near=%1, far=%2, Ka=%3]").arg(m_near).arg(m_far).arg(m_Ka);
	}
private:
	float m_near, m_far, m_Ka, m_gamma;
};

GROUP_NAMESPACE_END

NORI_REGISTER_GROUP_CLASS(Depth, "depth");
NORI_NAMESPACE_END
