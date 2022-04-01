#ifndef Sphere_H
#define Sphere_H
#include "Vec3.h"
#include <vector>
#include "Mesh.h"
#include <cmath>

struct RaySphereIntersection {
    bool intersectionExists;
    float t;
    float theta,phi; //pour les coordonnées sphériques (utiles pour plaquer les textures)
    Vec3 intersection;
    Vec3 secondintersection;
    Vec3 normal;
};

static
Vec3 SphericalCoordinatesToEuclidean( Vec3 ThetaPhiR ) {
    return ThetaPhiR[2] * Vec3( cos(ThetaPhiR[0]) * cos(ThetaPhiR[1]) , sin(ThetaPhiR[0]) * cos(ThetaPhiR[1]) , sin(ThetaPhiR[1]) );
}
static
Vec3 SphericalCoordinatesToEuclidean( float theta , float phi ) {
    return Vec3( cos(theta) * cos(phi) , sin(theta) * cos(phi) , sin(phi) );
}

static
Vec3 EuclideanCoordinatesToSpherical( Vec3 xyz ) {
    float R = xyz.length();
    float phi = asin( xyz[2] / R );
    float theta = atan2( xyz[1] , xyz[0] );
    return Vec3( theta , phi , R );
}



class Sphere : public Mesh {
public:
    Vec3 m_center;
    float m_radius;

    Sphere() : Mesh() {}
    Sphere(Vec3 c , float r) : Mesh() , m_center(c) , m_radius(r) {}

    void build_arrays(){
        unsigned int nTheta = 20 , nPhi = 20;
        positions_array.resize(3 * nTheta * nPhi );
        normalsArray.resize(3 * nTheta * nPhi );
        uvs_array.resize(2 * nTheta * nPhi );
        for( unsigned int thetaIt = 0 ; thetaIt < nTheta ; ++thetaIt ) {
            float u = (float)(thetaIt) / (float)(nTheta-1);
            float theta = u * 2 * M_PI;
            for( unsigned int phiIt = 0 ; phiIt < nPhi ; ++phiIt ) {
                unsigned int vertexIndex = thetaIt + phiIt * nTheta;
                float v = (float)(phiIt) / (float)(nPhi-1);
                float phi = - M_PI/2.0 + v * M_PI;
                Vec3 xyz = SphericalCoordinatesToEuclidean( theta , phi );
                positions_array[ 3 * vertexIndex + 0 ] = m_center[0] + m_radius * xyz[0];
                positions_array[ 3 * vertexIndex + 1 ] = m_center[1] + m_radius * xyz[1];
                positions_array[ 3 * vertexIndex + 2 ] = m_center[2] + m_radius * xyz[2];
                normalsArray[ 3 * vertexIndex + 0 ] = xyz[0];
                normalsArray[ 3 * vertexIndex + 1 ] = xyz[1];
                normalsArray[ 3 * vertexIndex + 2 ] = xyz[2];
                uvs_array[ 2 * vertexIndex + 0 ] = u;
                uvs_array[ 2 * vertexIndex + 1 ] = v;
            }
        }
        triangles_array.clear();
        for( unsigned int thetaIt = 0 ; thetaIt < nTheta - 1 ; ++thetaIt ) {
            for( unsigned int phiIt = 0 ; phiIt < nPhi - 1 ; ++phiIt ) {
                unsigned int vertexuv = thetaIt + phiIt * nTheta;
                unsigned int vertexUv = thetaIt + 1 + phiIt * nTheta;
                unsigned int vertexuV = thetaIt + (phiIt+1) * nTheta;
                unsigned int vertexUV = thetaIt + 1 + (phiIt+1) * nTheta;
                triangles_array.push_back( vertexuv );
                triangles_array.push_back( vertexUv );
                triangles_array.push_back( vertexUV );
                triangles_array.push_back( vertexuv );
                triangles_array.push_back( vertexUV );
                triangles_array.push_back( vertexuV );
            }
        }
    }


    virtual RayIntersection intersect(const Ray &ray) const {
        RayIntersection intersection;
        float a = (ray.direction()).squareNorm();
        float b = 2 * (Vec3::dot(ray.direction(), (ray.origin() - m_center) ));
        float c = (ray.origin() - m_center).squareNorm() - pow(m_radius, 2);

        float discr = (b*b) - 4 *(a * c);
        float s1, s2;
        if(discr > 0.0001f){

            s1 = (-b - sqrt(discr)) / (2*a);
            s2 = (-b + sqrt(discr)) / (2*a);

            float lower_t, higher_t;
            
            if(s1 < s2) {
                lower_t = s1;
                higher_t = s2;
            }else{

                lower_t = s2;
                higher_t = s1;
            }


            if(lower_t >= 0.0001f){

                intersection.intersectionExists = true;
                intersection.t = lower_t;
                //remplacement du t dans l'expression parametrique du rayon
                Vec3 inter1 = ray.origin() + intersection.t*(ray.direction());
                Vec3 inter2 = ray.origin() + higher_t*(ray.direction());

                // calcul de la normale de la première intersection
                Vec3 inter_normal = inter1  - m_center;
                inter_normal.normalize();

                intersection.intersection = inter1;
                intersection.intersection_2 = inter2;
                intersection.normal = inter_normal;
                intersection.objType = 1;
            }

        }else {

            intersection.intersectionExists = false;
        }

        return intersection;
    }
};
#endif
