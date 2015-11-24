/*******************************************************************************
 *  diffuse.cpp
 *******************************************************************************
 *  Template taken from Nori, without its implementation
 *  Copyright (c) 2012 by Wenzel Jakob and Steve Marschner.
 ***********************************************/

#include <nori/bsdf.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Diffuse / Lambertian BRDF model
 */
class Diffuse : public BSDF {
public:
	Diffuse(const PropertyList &propList) {
		m_albedo = propList.getColor("albedo", Color3f(0.5f));
	}

	/// Evaluate the BRDF model
	Color3f eval(const BSDFQueryRecord &) const {
            throw NoriException("Not implemented!");
	}

	/// Compute the density of \ref sample() wrt. solid angles
	float pdf(const BSDFQueryRecord &) const {
		throw NoriException("Not implemented!");
	}

	/// Draw a a sample from the BRDF model
	Color3f sample(BSDFQueryRecord &, const Point2f &) const {
		throw NoriException("Not implemented!");
	}

	/// Return a human-readable summary
	QString toString() const {
		return QString(
			"Diffuse[\n"
			"  albedo = %1\n"
			"]").arg(m_albedo.toString());
	}
        
        Color3f getColor() const { return m_albedo; }

	EClassType getClassType() const { return EBSDF; }
private:
	Color3f m_albedo;
};

NORI_REGISTER_CLASS(Diffuse, "diffuse");
NORI_NAMESPACE_END