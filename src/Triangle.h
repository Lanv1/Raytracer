#ifndef TRIANGLE_H
#define TRIANGLE_H
#include "Vec3.h"
#include "Ray.h"
#include "Plane.h"

struct RayIntersection{
    bool intersectionExists;
    int objType;                // 0: mesh, 1:sphere, 2:square
    float t;
    float w0,w1,w2;         //coordonées barycenrtiques de l'intersection
    unsigned int tIndex =  -1;
    Vec3 intersection;
    Vec3 intersection_2;
    Vec3 normal;
};

class Triangle {
private:
    Vec3 m_c[3] , m_normal;
    float area;
public:
    Triangle() {}
    Triangle( Vec3 const & c0 , Vec3 const & c1 , Vec3 const & c2 ) {
        m_c[0] = c0;
        m_c[1] = c1;
        m_c[2] = c2;
        updateAreaAndNormal();
        if(area == 0)
            std::cout<<"aire nulle: "<<c0<<", "<<c1<<", "<<c2<<std::endl;
    }
    void updateAreaAndNormal() {
        Vec3 nNotNormalized = Vec3::cross( m_c[1] - m_c[0] , m_c[2] - m_c[0] );
        float norm = nNotNormalized.length();
        m_normal = nNotNormalized / norm;
        area = norm / 2.f;
    }
    void setC0( Vec3 const & c0 ) { m_c[0] = c0; } // remember to update the area and normal afterwards!
    void setC1( Vec3 const & c1 ) { m_c[1] = c1; } // remember to update the area and normal afterwards!
    void setC2( Vec3 const & c2 ) { m_c[2] = c2; } // remember to update the area and normal afterwards!
    Vec3 const & normal() const { return m_normal; }
    Vec3 projectOnSupportPlane( Vec3 const & p ) const {
        Vec3 result;
        return result;
    }

    float squareDistanceToSupportPlane( Vec3 const & p ) const {
        float result;
        
        return result;
    }

    float distanceToSupportPlane( Vec3 const & p ) const {
        float result;
        Vec3 PA = m_c[0] - p;

        result = Vec3::dot(PA, m_normal) / m_normal.length();

        return result;
    }

    bool isParallelTo( Line const & L ) const {
        if(fabs(Vec3::dot(L.direction(), m_normal)) <= 0.0001f)
            return true;
        return false;
    }

    float getIntersectionTWithSupportPlane( Line const & L ) const {
        float result = -1;
        if(!isParallelTo(L)){
            float denom = Vec3::dot(L.direction(), m_normal);

            float D = Vec3::dot(m_c[0] - L.origin(), m_normal);
            float t = D / denom;
            
            // check si l'intersection est bien devant "la caméra".
            if(t > 0.0) 
                result = t;

        }

        return result;
    }

    void computeBarycentricCoordinates( Vec3 const & p , float & u0 , float & u1 , float & u2 ) const {
        //TODO Complete
        Vec3 e1 = m_c[1] - m_c[0];
        Vec3 e2 = m_c[2] - m_c[1];
        Vec3 e3 = m_c[0] - m_c[2];

        Vec3 p1 = p - m_c[0];
        Vec3 p2 = p - m_c[1];
        Vec3 p3 = p - m_c[2];

        Vec3 cross_0 = Vec3::cross(e1, p1);
        float area_0 = cross_0.length() / 2.0f;

        Vec3 cross_1 = Vec3::cross(e2, p2);
        float area_1 = cross_1.length() / 2.0f; 
        
        Vec3 cross_2 = Vec3::cross(e3, p3);
        float area_2 = cross_2.length() / 2.0f; 

        u2 = area_0 / this->area;
        u0 = area_1 / this->area;
        u1 = area_2 / this->area;
    }

    RayIntersection getIntersection( Ray const & ray ) const {
        RayIntersection result;
        float t = getIntersectionTWithSupportPlane(ray);
        float eps = 0.000001;

        if(t > eps) {

            Vec3 inter = ray.origin() + t * ray.direction();
            float a, b, c;
            computeBarycentricCoordinates(inter, a, b, c);
 
            if(b >= 0 && c >= 0 && a >= 0 && (a + b + c <= 1)){    

                result.t = t;
                result.intersection = inter;
                result.objType = 0;
                result.w0 = a;
                result.w1 = b;
                result.w2 = c;
                result.intersectionExists = true;
                result.normal = m_normal;

            }else {
                result.intersectionExists = false;
            }
        }else {

            result.intersectionExists = false;
        }


        return result;
    }
};
#endif
