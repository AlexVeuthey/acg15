/*******************************************************************************
 *  path.cpp
 *******************************************************************************
 *  Copyright (c) 2013 Alexandre Kaspar <alexandre.kaspar@a3.epfl.ch>
 *  For Advanced Computer Graphics, at the LGG / EPFL
 *
 *        DO NOT REDISTRIBUTE
 ***********************************************/

#include <nori/acg.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/luminaire.h>
#include <nori/bsdf.h>
#include <nori/medium.h>
#include <nori/phase.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// put your group number here!
#define GROUP_NUMBER 8
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

GROUP_NAMESPACE_BEGIN()

float dot(Normal3f &vec1, const Vector3f &vec2){
   return vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
}

/**
 * Simple path tracer implementation
 */
class PathTracer : public Integrator {
public:

        PathTracer(const PropertyList &) {
        }

        /// Return the mesh corresponding to a given luminaire
        inline const Mesh *getMesh(const Luminaire *lum) const {
                const Mesh *mesh = dynamic_cast<const Mesh *> (lum->getParent());
                if (!mesh) throw NoriException("Unhandled type of luminaire!");
                return mesh;
        }

        /**
         * \brief Directly sample the lights, providing a sample weighted by 1/pdf
         * where pdf is the probability of sampling that given sample
         * 
         * \param scene
         * the scene to work with
         * 
         * \param lRec
         * the luminaire information storage
         * 
         * \param _sample
         * the 2d uniform sample
         * 
         * \return the sampled light radiance including its geometric, visibility and pdf weights
         */
        inline Color3f sampleLights(const Scene *scene, LuminaireQueryRecord &lRec, const Point2f &_sample) const {
                Point2f sample(_sample);
                const std::vector<Luminaire *> &luminaires = scene->getLuminaires();

                if (luminaires.size() == 0)
                        throw NoriException("LightIntegrator::sampleLights(): No luminaires were defined!");

                // 1. Choose one luminaire at random
                int index = std::min((int) (luminaires.size() * sample.x()), (int) luminaires.size() - 1);
                sample.x() = luminaires.size() * sample.x() - index; // process sample to be Unif[0;1] again

                // 2. Sample the position on the luminaire mesh
                // using Mesh::samplePosition(const Point2d &sample, Point3f &p, Normal3f &n)
                lRec.luminaire = luminaires[index];
                const Mesh *mesh = getMesh(lRec.luminaire);
                mesh->samplePosition(sample, lRec.p, lRec.n);
                lRec.d = lRec.p - lRec.ref;

                // 3. Compute distance between the two points (from first mesh, to luminaire mesh)
                float dist2 = lRec.d.squaredNorm();
                lRec.dist = std::sqrt(dist2);
                lRec.d /= lRec.dist;

                // 4. Correct side of luminaire
                // /!\ if on the wrong side, then we get no contribution!
                float dp = -lRec.n.dot(lRec.d);
                lRec.pdf = dp > 0 ? mesh->pdf() * dist2 / dp : 0.0f;

                if (dp > 0) {
                        // 5. Check the visibility
                        if (scene->rayIntersect(Ray3f(lRec.ref, lRec.d, Epsilon, lRec.dist * (1 - 1e-4f))))
                                return Color3f(0.0f);
                        // 6. Geometry term on luminaire's side
                        // Visiblity + Geometric term on the luminaire's side 
                        //      G(x, x', w, w') = ( cos(w) cos(w') ) / ||x - x'||^2
                        float G_lum = dp / dist2;

                        // 7. Radiance from luminaire
                        Color3f value = lRec.luminaire->getColor();

                        return value * G_lum * luminaires.size() / mesh->pdf();
                } else {
                        // wrong side of luminaire!
                        return Color3f(0.0f);
                }
        }

        Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &_ray) const {
                Ray3f ray(_ray);
                Intersection its;
                Color3f result(0.0f), throughput(1.0f);
                int bounces = 0;
                const int MAX_BOUNCES = 10000; //put a limit to the max number of bounces
                bool cutOff = false;
                double q = 1.0;

                //for (bounces = 0; bounces < 2; bounces++)
                while (!cutOff && bounces < MAX_BOUNCES)
                {
                    // Step 1: Intersect the ray with the scene. Return environment
                    // luminaire if no hit.
                    // This means that we have to return the current result we have,
                    // which will never be augmented again because ray is going "nowhere"
                    // but can already be non-zero if bounces > 0.
                    if (!scene->rayIntersect(ray, its))
                    {
                        return result;
                    }

                    // Step 2: Check if the ray hit a light source.
                    if (its.mesh->isLuminaire())
                    {
                        // TODO Don't know what to do here exactly...
                        //LuminaireQueryRecord lRecDir(its.mesh->getLuminaire(), ray.o, its.p, its.shFrame.n);
                        //result += throughput * its.mesh->getLuminaire()->eval(lRecDir);
                        if (bounces == 0)
                        {
                            return its.mesh->getLuminaire()->getColor();
                        }
                        else
                        {
                            return result;
                        }
                    }

                    // Step 3: Direct illumination sampling.
                    // Sample the light intensity from a random light,
                    // then apply the BRDF for sampled direction and take current throughput into account.
                    LuminaireQueryRecord lRec(its.p);
                    Color3f lightIntensity = sampleLights(scene, lRec, sampler->next2D());
                    BSDFQueryRecord bRec(its.toLocal(-ray.d/sqrt(ray.d[0]*ray.d[0]+ray.d[1]*ray.d[1]+ray.d[2]*ray.d[2])), its.toLocal(lRec.d), ESolidAngle);
                    result += throughput * lightIntensity * its.mesh->getBSDF()->eval(bRec) * dot(its.shFrame.n, lRec.d);

                    // Step 4: Recursively sample indirect illumination
                    // Sample a new direction according to BRDF, modify throughput accordingly,
                    // then prepare the corresponding ray for the next iterative recursion step.
                    BSDFQueryRecord bRecInd(its.toLocal(-ray.d/sqrt(ray.d[0]*ray.d[0]+ray.d[1]*ray.d[1]+ray.d[2]*ray.d[2])));
                    throughput *= its.mesh->getBSDF()->sample(bRecInd, sampler->next2D());
                    ray = Ray3f(its.p, its.toWorld(bRecInd.wo));

                    // Step 5. Apply Russioan Roullette after 2 main bounces.
                    if (bounces >= 1)
                    {
                        double qx = sampler->next1D();
                        if (qx > q)
                        {
                            throughput *= (double) 1.0 / (double)(1-q);
                        }
                        else
                        {
                            cutOff = true;
                        }
                    }

                    // Increment bounces number
                    bounces++;
                }

                return result;
        }

        QString toString() const {
                return "PathTracer[]";
        }
};

GROUP_NAMESPACE_END

NORI_REGISTER_GROUP_CLASS(PathTracer, "path");
NORI_NAMESPACE_END
