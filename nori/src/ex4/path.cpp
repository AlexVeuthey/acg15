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

float length(Vector3f &vec){
   return sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
}

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
         const int NB_BOUNCES = 10;
         float q = 0.1;
         
         // TODO implement a path tracer
         
         for(int i = 0; i <= NB_BOUNCES; i++){
         
            // Step 1: Intersect the ray with the scene. Return environment
            // luminaire if no hit.
            if (!scene->rayIntersect(ray, its))
               return result;
            
            // Step 2: Check if the ray hit a light source.
            const Mesh *mesh = its.mesh;
            const BSDF *bsdf = mesh->getBSDF();
            
            if( mesh->isLuminaire() ){
               if(i != 0){
                  return result;
               }
               const Luminaire *luminaire = its.mesh->getLuminaire();
               LuminaireQueryRecord lqr(luminaire, ray.o, its.p, its.shFrame.n);
               return luminaire->eval(lqr);
            }
            
            // Step 3: Direct illumination sampling.
            LuminaireQueryRecord lqr(its.p);
            
            //lqr will be modified in this function and will have the wo at x' in lqr.d
            Color3f light_sample = sampleLights(scene, lqr, sampler->next2D());
            
            BSDFQueryRecord bRec(its.toLocal(-ray.d/length(ray.d)), its.toLocal(lqr.d), ESolidAngle);
            
            result += light_sample*bsdf->eval(bRec) * throughput * dot(its.shFrame.n, lqr.d);
            
            // Step 4: Recursively sample indirect illumination
            throughput *= bsdf->sample(bRec, sampler->next2D());
            
            //next ray, going from the current point to the direction of W_o
            ray = Ray3f(its.p, its.toWorld(bRec.wo));
            
            // Step 5. Apply Russion Roullette after 2 main bounces.
            if(i >= 2 && sampler->next1D() < q){
               return result;
            }
            else if(i >= 2){
               throughput *= 1.f/(1-q);
            }
         }
         
         
         ///TODO
         
         return result;
        }

        QString toString() const {
                return "PathTracer[]";
        }
};

GROUP_NAMESPACE_END

NORI_REGISTER_GROUP_CLASS(PathTracer, "path");
NORI_NAMESPACE_END
