#ifndef SQUARE_H
#define SQUARE_H
#include "Vec3.h"
#include <vector>
#include "Mesh.h"
#include <cmath>

struct RaySquareIntersection{
    bool intersectionExists;
    float t;
    float u,v;
    Vec3 intersection;
    Vec3 normal;
};


class Square : public Mesh {
public:
    Vec3 m_normal;
    Vec3 m_bottom_left;
    Vec3 m_right_vector;
    Vec3 m_up_vector;

    Square() : Mesh() {}
    Square(Vec3 const & bottomLeft , Vec3 const & rightVector , Vec3 const & upVector , float width=1. , float height=1. ,
           float uMin = 0.f , float uMax = 1.f , float vMin = 0.f , float vMax = 1.f) : Mesh() {
        setQuad(bottomLeft, rightVector, upVector, width, height, uMin, uMax, vMin, vMax);
    }

    void setQuad( Vec3 const & bottomLeft , Vec3 const & rightVector , Vec3 const & upVector , float width=1. , float height=1. ,
                  float uMin = 0.f , float uMax = 1.f , float vMin = 0.f , float vMax = 1.f) {
        m_right_vector = rightVector;
        m_up_vector = upVector;
        m_normal = Vec3::cross(rightVector , upVector);
        m_bottom_left = bottomLeft;

        m_normal.normalize();
        m_right_vector.normalize();
        m_up_vector.normalize();

        m_right_vector = m_right_vector*width;
        m_up_vector = m_up_vector*height;

        vertices.clear();
        vertices.resize(4);
        vertices[0].position = bottomLeft;                                      vertices[0].u = uMin; vertices[0].v = vMin;
        vertices[1].position = bottomLeft + m_right_vector;                     vertices[1].u = uMax; vertices[1].v = vMin;
        vertices[2].position = bottomLeft + m_right_vector + m_up_vector;       vertices[2].u = uMax; vertices[2].v = vMax;
        vertices[3].position = bottomLeft + m_up_vector;                        vertices[3].u = uMin; vertices[3].v = vMax;
        vertices[0].normal = vertices[1].normal = vertices[2].normal = vertices[3].normal = m_normal;
        triangles.clear();
        triangles.resize(2);
        triangles[0][0] = 0;
        triangles[0][1] = 1;
        triangles[0][2] = 2;
        triangles[1][0] = 0;
        triangles[1][1] = 2;
        triangles[1][2] = 3;


    }

    virtual RayIntersection intersect(const Ray &ray) const {
        RayIntersection intersection;
        // bornes du quad
        Vec3 high;
        Vec3 low;

        // calcul des axes du quad
        Vec3 dep = vertices[2].position - vertices[0].position;
        Vec3 to_process;
        for(int i = 0; i < 3; i ++){
            high[i] = std::max(vertices[2].position[i], vertices[0].position[i]);   // up_right, bottom_left
            low[i] = std::min(vertices[2].position[i], vertices[0].position[i]);

            // axes a prendre en compte dans to_process (square: 2D)
            // pour restreindre le plan au carré.
            if(fabs(dep[i]) > 0.0001f)
                to_process[i] = 1;
        }

        // Intersection droite/plan
        float denom = Vec3::dot(ray.direction(), vertices[0].normal);
        Vec3 ray_origin_plane = vertices[0].position - ray.origin();
        float D = Vec3::dot(ray_origin_plane, vertices[0].normal);

        if(fabs(denom) > 0.0001f) {

            // t n'est pas infini (car dénominateur pas proche de 0) -> il y a une intersection avec le plan
            float t = D/denom;

            if(t >= 0.0001){
                //intersection devant la caméra
                Vec3 inter = ray.origin() + t * ray.direction();

                // on détermine si le rayon intersecte bien le quad défini (et non plus le plan)
                bool ray_in_plane = true;
                for(int i = 0; i < 3; i ++){
                    if(to_process[i]){
                        if(inter[i] > high[i] || inter[i] < low[i])
                            ray_in_plane = false;
                    }
                }
               
                if(ray_in_plane){
                    intersection.t = t;
                    intersection.intersectionExists = true;
                    intersection.intersection = inter;
                    intersection.normal = vertices[0].normal;
                    intersection.objType = 2;
                }else{

                    intersection.intersectionExists = false;
                }

            }
        } 
    
        
        return intersection;
    }
};
#endif // SQUARE_H