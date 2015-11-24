/*******************************************************************************
 *  light_integrator.cpp
 *******************************************************************************
 *  Copyright (c) 2013 Alexandre Kaspar <alexandre.kaspar@a3.epfl.ch>
 *  For Advanced Computer Graphics, at the LGG / EPFL
 * 
 *        DO NOT REDISTRIBUTE
 ***********************************************/

#include <nori/acg.h>
#include <nori/bsdf.h>
#include <nori/common.h>
#include <nori/integrator.h>
#include <nori/luminaire.h>
#include <nori/sampler.h>
#include <nori/scene.h>
#include <vector>

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
 * \brief Simple local illumination integrator
 * using light area sampling
 */
class LightIntegrator : public Integrator {
public:

        LightIntegrator(const PropertyList &propList) {
                Q_UNUSED(propList);
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

                // TODO Implement the following steps
                // and take care of using the good G, V terms to work with the Li method below
                
                // 1. Choose one luminaire at random
                const Luminaire *rand_lum = luminaires.at(sample[0]+sample[1]);
                
                                  
                // 2. Sample the position on the luminaire mesh
                // using Mesh::samplePosition(const Point2d &sample, Point3f &p, Normal3f &n)
                const Mesh *lum_mesh = getMesh(rand_lum);
                Point3f p_lum;
                Normal3f n_lum;
                //find a random point and normal on the luminaire mesh
                lum_mesh->samplePosition(sample, p_lum, n_lum);
                

                // 3. Compute geometry term G and visibility term on the luminaire's side (no cos(w) of the mesh side)
                // as well as the pdf of that point being found
                // use Mesh::pdf to get the probability of choosing the point in Mesh::samplePosition
                
                
                Vector3f vec_mesh_to_light(p_lum[0]-lRec.p[0], p_lum[1]-lRec.p[1], p_lum[2]-lRec.p[2]);
                Vector3f vec_mesh_to_light_norm = vec_mesh_to_light/length(vec_mesh_to_light);
                
                Ray3f ray(lRec.p, vec_mesh_to_light_norm);
                
                float G = 1.0;
                
                //if sampled point don't see the sampled light point, we put the G parameters at 0
                Intersection its;
                if (scene->rayIntersect(ray, its)){
                   if( !its.mesh->isLuminaire() ){
                      G = 0.0;
                   }
                }
                
                //cos with normals at point x' and x''
                float cos_pp = dot(n_lum, -vec_mesh_to_light_norm);
                if(cos_pp < 0.0) cos_pp = 0.0;
                if(cos_pp > 1.0) cos_pp = 1.0;
                
                float cos_p = dot(lRec.n, vec_mesh_to_light_norm);
                if(cos_p < 0.0) cos_p = 0.0;
                if(cos_p > 1.0) cos_p = 1.0;
                
                G *= cos_p*cos_pp/pow(length(vec_mesh_to_light), 2);
                
                //weight the pdf of the light by all lights 
                //(since we have a 1/#lights chance to choose each light)
                float pdf = getMesh(rand_lum)->pdf();
                pdf = pdf/luminaires.size();
                
                lRec.d = vec_mesh_to_light_norm;
                lRec.luminaire = rand_lum;
                
                //luminaire record to evaluate the luminaire we are sampling
                LuminaireQueryRecord lqr(
                         rand_lum,
                         lRec.p, p_lum, n_lum);
                
                // 4. Return radiance emitted from luminaire multiplied by the appropriate terms G, V ...
                return rand_lum->eval(lqr)*G*pdf;
        }

        /**
         * \brief Simple local illumination integration:
         * We cast a ray from the camera, intersects it with the first element
         * in the scene it goes through and from there we directly sample the
         * light's in the scene to compute direct lighting.
         */ 
        Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray_) const {
                Ray3f ray(ray_);

                /* Find the surface that is visible in the requested direction */
                Intersection its;
                if (!scene->rayIntersect(ray, its))
                        return Color3f(0.0f);
                

                const Mesh *mesh = its.mesh;
                const BSDF *bsdf = mesh->getBSDF();

                /// TODO implement direct lighting using light sampling using
                //      sampleLights(const Scene *, LuminaireQueryRecord &, const Point2d &)
                // which you also have to implement
                
                Point2f sample = sampler->next2D();
                
                if (mesh->isLuminaire()) {
                        const Luminaire *luminaire = its.mesh->getLuminaire();
                        LuminaireQueryRecord lRec(luminaire, ray.o, its.p, its.shFrame.n);
                        return luminaire->eval(lRec);
                }
                
                LuminaireQueryRecord lqr(NULL, ray.o, its.p, its.shFrame.n);
                
                //lqr will be modified in this function and will have the wo at x' in lqr.d
                Color3f radiance = sampleLights(scene, lqr, sample);
                
                BSDFQueryRecord bRec(its.toLocal(-ray.d/length(ray.d)), its.toLocal(lqr.d), ESolidAngle);
                
                //returns the phong calculation with the radiance
                return bsdf->eval(bRec) * radiance;
        }

        QString toString() const {
                return QString("LightIntegrator[]");
        }
};

GROUP_NAMESPACE_END

NORI_REGISTER_GROUP_CLASS(LightIntegrator, "light");
NORI_NAMESPACE_END
