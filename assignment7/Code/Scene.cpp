//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"

void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}


Vector3f Scene::shade(const Ray &ray) const
{
    // return  Vector3f(1.f, 1.f, 1.f);
    Vector3f direct_energy(0.f, 0.f, 0.f);
    Vector3f bounced_energy(0.f, 0.f, 0.f);
    Intersection pos_on_obj = intersect(ray);
    if(!pos_on_obj.happened) return direct_energy;
    Intersection pos_on_ls;
    float pdf_ls;    
    sampleLight(pos_on_ls, pdf_ls);
    // here we need to use hasEmit
    if(pos_on_obj.obj->hasEmit()) return pos_on_ls.emit;            
    // shoot to the light
    Ray direct_ray = Ray(pos_on_obj.coords, normalize(pos_on_ls.coords - pos_on_obj.coords));
    // check if the light is blocked by any other objects
    Intersection blocked = intersect(direct_ray);


    if(!blocked.happened || blocked.obj->hasEmit())
    {
        // brdf of light source

        //here the cos_theta is zero.
        Vector3f f_r_ls = pos_on_obj.m->eval( - ray.direction, direct_ray.direction, pos_on_obj.normal);
        float cos_theta = dotProduct(direct_ray.direction, pos_on_obj.normal);
        float cos_theta_ = dotProduct(-direct_ray.direction, pos_on_ls.normal);
        direct_energy = pos_on_ls.emit * f_r_ls * cos_theta * cos_theta_ / squaredNorm(pos_on_ls.coords - pos_on_obj.coords) / pdf_ls;
    }
    // Use russian Roulette to determine whether we should terminate the process
    float p = get_random_float();
    if(p < RussianRoulette && !pos_on_obj.obj->hasEmit())
    {
        Vector3f w_o = pos_on_obj.m->sample(- ray.direction, pos_on_obj.normal);
        Vector3f f_r_obj = pos_on_obj.m->eval(- ray.direction, w_o, pos_on_obj.normal);
        float pdf_obj = pos_on_obj.m->pdf(- ray.direction, w_o, pos_on_obj.normal);
        Ray next_ray(pos_on_obj.coords, w_o);
        Vector3f radiance = shade(next_ray);
        float cos_theta = dotProduct(next_ray.direction, pos_on_obj.normal);
        bounced_energy += f_r_obj * radiance * cos_theta / pdf_obj / RussianRoulette;
    }
    return  direct_energy + bounced_energy;
}
// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    // TO DO Implement Path Tracing Algorithm here
    // firstly, we compute the intersection, if the ray doesn't hit any object, return zero.
    // secondly, we sample the light and compute the pdf
    return shade(ray);
}