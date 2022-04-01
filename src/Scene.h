#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <string>
#include "Mesh.h"
#include "Sphere.h"
#include "Square.h"


#include <GL/glut.h>

using namespace std;


enum LightType {
    LightType_Spherical,
    LightType_Quad
};


struct Light {
    Vec3 material;
    bool isInCamSpace;
    LightType type;

    Vec3 pos;
    float radius;

    Mesh quad;

    float powerCorrection;

    Light() : powerCorrection(1.0) {}

};

struct RaySceneIntersection{
    bool intersectionExists;
    unsigned int typeOfIntersectedObject;
    unsigned int objectIndex;
    float t;
    RayIntersection rayIntersection;
    // RayIntersection rayMeshIntersection;
    // RayIntersection raySphereIntersection;
    // RayIntersection raySquareIntersection;

    RaySceneIntersection() : intersectionExists(false) , t(FLT_MAX) {
    }
    //si t == FLT_MAX après le calcul -> pas d'intersection
};



class Scene {
    std::vector< Mesh > meshes;         //typeOfIntersectedObject: 0
    std::vector< Sphere > spheres;      //1
    std::vector< Square > squares;      //2
    std::vector< Light > lights;

    // SceneObj contiendra tous les objets de la scène, besoin d'utiliser
    // des références sur ceux-ci pour appeler la méthode "intersect" la plus spécialisée (en fonction du type de l'objet).
    std::vector<Mesh*> sceneObj;


public:


    Scene() {
        
    }

    void updateSceneObj(){
        int total_size = meshes.size() + spheres.size() +squares.size();
        sceneObj.clear();
        sceneObj.resize(total_size);

        int size = meshes.size();
        int k = 0;
        for(int i = 0; i < size; i ++) {
            sceneObj[k] =  &meshes[i];
            k++;
        }

        size = spheres.size();
        for(int i = 0; i < size; i ++) {
            sceneObj[k] =  &spheres[i];
            k++;
        }

        size = squares.size();
        for(int i = 0; i < size; i ++) {
            sceneObj[k] =  &squares[i];
            k++;
        }

        cout<<"SCENEOBJ.SZ= "<<sceneObj.size()<<", meshes: "<<meshes.size()<<endl;

    }

    void draw() {
        // iterer sur l'ensemble des objets, et faire leur rendu :
        for( unsigned int It = 0 ; It < meshes.size() ; ++It ) {
            Mesh const & mesh = meshes[It];
            mesh.draw();
        }
        for( unsigned int It = 0 ; It < spheres.size() ; ++It ) {
            Mesh const & sphere = spheres[It];
            sphere.draw();
        }
        for( unsigned int It = 0 ; It < squares.size() ; ++It ) {
            Mesh const & square = squares[It];
            square.draw();
        }
    }


    RaySceneIntersection computeIntersection(Ray const & ray) {
        RaySceneIntersection result;

        int sceneSize = sceneObj.size();
        for(int i = 0; i < sceneSize; i ++) {
            RayIntersection inter = sceneObj[i]->intersect(ray);

            if(inter.intersectionExists) {
                if(inter.t < result.t){
                    result.t = inter.t;
                    result.objectIndex = i;
                    result.typeOfIntersectedObject = inter.objType;
                    result.intersectionExists = true;
                    result.rayIntersection = inter;
                }
            }
        
        }

        return result;

    }

    // intersection avec un objet quelconque -> ombre
    bool isInShadow(Ray const & shadowRay, Vec3 lightPos, int interObjIndex) {
        Vec3 interVec, toLight;
        int objSize = sceneObj.size();

        for(int i = 0; i < objSize; i ++) {

            RayIntersection inter = sceneObj[i]->intersect(shadowRay);
            // pas de self-intersect ni de shadow-acne.
            if(inter.intersectionExists && inter.t > 0.0001f && i != interObjIndex ) {
                
                // comparer les distances [origin, light] et [origin, intersection]
                // pour vérifier qu'il s'agit d'une intersection entre la position de la lumière
                // et le point d'origine du shadow-ray.

                toLight = lightPos - shadowRay.origin();
                interVec = inter.intersection - shadowRay.origin();

                if(interVec.squareLength() < toLight.squareLength() && !(inter.normal == Vec3(0., 0., 0.)))
                    return true;

            }
        }
        return false;


    }


    Vec3 computeBlinnPhong(float cosTheta, std::vector<Vec3> K, std::vector<Vec3> I, Vec3 R, Vec3 V, float alpha) {
        Vec3 color = Vec3(0., 0., 0.);

        float current_chan;         // R G B
        for(int i = 0; i < 3; i ++) {

            current_chan =  (K[0][i] * I[0][i]) * cosTheta +        
                            (K[1][i] * I[1][i]) * pow(max(Vec3::dot(R, V), (float)0.), alpha);
                       

            if(current_chan < 0)
                current_chan = 0;
            
            color[i] = current_chan;
        }
        
        return color;
    }

    Vec3 getRandPoint(float x_min, float x_max, float y_min, float y_max){
        // srand(time(NULL));
        float x = x_min + float(rand()) / float((RAND_MAX)) * (x_max-x_min);    // % x_max;
        float y = y_min + float(rand()) / float((RAND_MAX)) * (y_max-y_min);    // % x_max;
        return Vec3(x, 0, y);
    }

    // Echantillonnage aléatoire de la grille de lumières (grille xz). 
    // u -> z, v -> x
    std::vector<Vec3> getAreaLight(Vec3 center,Vec3 u, Vec3 v, int resolution) {
        std::vector<Vec3> grid;
        Vec3 topLeft = center - v/2. - u/2.;

        grid.resize(resolution * resolution);   //grille 2D applatie
        Vec3 yStep = u/float(resolution);
        Vec3 xStep = v/float(resolution);
        Vec3 minPt, maxPt, randXZ;
        for(size_t i = 0; i < resolution; i ++) {
            for(size_t j = 0; j < resolution; j ++) {

                minPt = topLeft + (j*xStep) + (i*yStep); // point haut-gauche de la cellule courante.
                maxPt = minPt + xStep + yStep;           // point bas-droite de la cellule actuelle.

                randXZ = getRandPoint(minPt[0], maxPt[0], minPt[2], maxPt[2]);
                randXZ[1] = center[1];
                grid[i * resolution + j] = randXZ;

            }
        }

        return grid;
    }

    

    Vec3 rayTraceRecursive( Ray ray , int NRemainingBounces ) {

        int areaLightRes = 1; // resolution de la grille a échantillonner aléatoirement.
        int gdSize = areaLightRes * areaLightRes;   // nombre de shadowrays envoyés.
        float soften = 0;
        Ray shadowRay;
        Vec3 color = Vec3(0., 0., 0.);
        Vec3 radiance = Vec3(0., 0., 0.);
        RaySceneIntersection raySceneIntersection = computeIntersection(ray);
        Vec3 L, interN, R, V, inter;       // 
        float cosTheta;
        int type, index;
        double shininess;
        bool mirror_surf = false;
        std::vector<Vec3> K;    // constantes ambiantes, diffuses et speculaires de l'objet.
        std::vector<Vec3> I;    // constante de la lumière (de base une seule constante mais on peut en ajouter).

        K.resize(3);
        I.resize(3);

        if(raySceneIntersection.intersectionExists){
            type = raySceneIntersection.typeOfIntersectedObject;
            inter = raySceneIntersection.rayIntersection.intersection;
            index = raySceneIntersection.objectIndex;
            
            K[0] = sceneObj[index]->material.diffuse_material;
            K[1] = sceneObj[index]->material.specular_material;
            K[2] = sceneObj[index]->material.ambient_material; 

            I[0] = lights[0].material;
            I[1] = lights[0].material;

            //Normale de l'intersection trouvée.
            interN = raySceneIntersection.rayIntersection.normal;

            L = lights[0].pos - inter;
            L.normalize();

            cosTheta = Vec3::dot(L, interN);
            R = (2*cosTheta * interN) - L;
            R.normalize();                  

            V = ray.origin() - inter;
            V.normalize();
            
            float ct = Vec3::dot(interN, ray.direction());

            Vec3 ambient(0., 0., 0.);   // la composante ambiante ne doit pas être atténuée.
            for(int i = 0; i < 3; i ++) {
                ambient[i] += I[0][i] * K[2][i];
            }

            //réflexion
            if(sceneObj[index]->material.type == MaterialType::Material_Mirror && NRemainingBounces > 0){
                
                Vec3 R1 = ray.direction() - 2*ct * interN;
                R1.normalize();
                Ray reflect(inter, R1 * 0.0001);
                color += 0.8 * rayTraceRecursive(reflect, NRemainingBounces - 1); 
            
            //réfraction
            }else if(sceneObj[index]->material.type == MaterialType::Material_Glass && NRemainingBounces > 0) {
                
                float ior = 1. / sceneObj[index]->material.index_medium;    //index de l'air / index du matériau.
                float sinT2 = ior * ior * (1.0 - ct * ct);
                if(sinT2 <= 1.0){
                    
                    Vec3 TR = ior * ray.direction() - (ior + sqrt(1.0 - sinT2)) * interN;
                    Ray refract(inter, TR);
                    color += 0.8 * rayTraceRecursive(refract, NRemainingBounces - 1);
                
                }

            }else {

                shininess = sceneObj[index]->material.shininess;
                color = computeBlinnPhong(cosTheta, K, I, R, V, shininess);

                //SOFT SHADOWS
                std::vector<Vec3> areaLight = getAreaLight(lights[0].pos, Vec3(0,0,0.2), Vec3(0.2,0,0), areaLightRes);
                for(size_t i = 0; i < gdSize; i ++){

                    L = areaLight[i] - inter;
                    L.normalize();
                    Ray shadowRay(inter, L);

                    if(!isInShadow(shadowRay, lights[0].pos, index))
                        soften ++;

                }
                soften /= (float) gdSize;  // Proportion de shadowRays n'intersectant aucun objets.
                
                color *= soften;  
                color += ambient;
                
            }
            

        }
        return color;
    }


    Vec3 rayTrace( Ray const & rayStart ) {
        //appeler la fonction recursive
        Vec3 color = rayTraceRecursive(rayStart, 3) ;
        return color;
    }

    void setup_single_sphere(Vec3 centerPos) {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5,5,5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }
        {
            spheres.resize( spheres.size() + 1 );
            Sphere & s =  spheres[spheres.size() - 1];
            s.m_center = centerPos;
            s.m_radius = 1.f;
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.7,0.3,1 );
            s.material.specular_material = Vec3( 0.2,0.2,0.2 );
            s.material.shininess = 20;
        }

        updateSceneObj();
    }

    void add_sphere(Vec3 sphere_pos) {
        spheres.resize(spheres.size() + 1);
        Sphere &s = spheres[spheres.size() - 1];
        s.m_center = sphere_pos;
        s.m_radius = 1.;
        s.build_arrays();

        s.material.type = Material_Mirror;
        s.material.diffuse_material = Vec3( 0.7,0.3,1 );
        s.material.specular_material = Vec3( 0.2,0.2,0.2 );
        s.material.shininess = 20;

        updateSceneObj();
    }

    void setup_single_square() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5,5,5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }

        {
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1,0.8,0.8 );
            s.material.specular_material = Vec3( 0.8,0.8,0.8 );
            s.material.ambient_material = Vec3( 0.2,0.2,0.1 );
            s.material.shininess = 20;
        }

        updateSceneObj();
    }

    void setup_cornell_box(){
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();
        sceneObj.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3( 0.0, 1.5, 0.0 );
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }

        { //Back Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.6,0.6,0.6 );
            s.material.specular_material = Vec3( 0.6,0.6,0.6 );
            s.material.shininess = 16;
        }

        { //Left Wall

            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));            
            s.rotate_y(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.,0.7,0.8 );
            s.material.specular_material = Vec3( 0.,0.4,0.3 );
            s.material.shininess = 16;
        }

        { //Right Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.0,0.3,0.4 );
            s.material.specular_material = Vec3( 0.0,0.2,0.1 );
            s.material.shininess = 16;
        }

        { //Floor
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.6, 0.6 ,0.6 );
            s.material.specular_material = Vec3( 0.6,0.6,0.6 );
            s.material.shininess = 16;
        }

        { //Ceiling
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(0.,0.5,0.6 );
            s.material.specular_material = Vec3( 0.,0.4,0.5 );
            s.material.shininess = 16;
        }

        { //MIRRORED Sphere
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,1.,1. );
            s.material.specular_material = Vec3(  1.,1.,1. );
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }


        { //GLASS Sphere

            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0, -1.25, 0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
            s.material.transparency = 1.0;
            s.material.index_medium = 1.4;
        }
        { //GLASS Sphere

            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-0.3, -1.9, 0.8);
            s.m_radius = 0.10f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
            s.material.transparency = 1.0;
            s.material.index_medium = 1.4;
        }

        {
            meshes.resize(meshes.size() + 1);
            Mesh & mesh = meshes[meshes.size() - 1];
            mesh.loadOFF("star.off");
            mesh.build_arrays();
            mesh.material.ambient_material = Vec3(0.0, 0.02, 0.05);
            mesh.material.diffuse_material = Vec3(0.8, 0.2, 0.);
            mesh.material.specular_material = Vec3(0.1, 0.4, 0.);
            mesh.material.shininess = 16;
            mesh.translate(Vec3(-1.8, -2.0, 0.8));
            mesh.scale(Vec3(0.8, 0.8, 0.8));
        }
        {
            meshes.resize(meshes.size() + 1);
            Mesh & mesh = meshes[meshes.size() - 1];
            mesh.loadOFF("star.off");
            mesh.build_arrays();
            mesh.material.ambient_material = Vec3(0.0, 0.02, 0.05);
            mesh.material.diffuse_material = Vec3(0.8, 0.2, 0.);
            mesh.material.specular_material = Vec3(0.1, 0.4, 0.);
            mesh.material.shininess = 16;
            mesh.translate(Vec3(0.2, -1.4, -1.5));
            mesh.scale(Vec3(0.8, 0.8, 0.8));
        }

        updateSceneObj();

    }

    void setup_mesh_box(){
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();
        sceneObj.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3( 0.0, 1.5, 0.0 );
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }

        { //Back Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.,0.6,0.6 );
            s.material.specular_material = Vec3( 0.,0.6,0.6 );
            s.material.shininess = 16;
        }

        { //Left Wall

            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));            
            s.rotate_y(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
        }

        { //Right Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.0,1.0,0.0 );
            s.material.specular_material = Vec3( 0.0,1.0,0.0 );
            s.material.shininess = 16;
        }

        { //Floor
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.5, 0.6 ,1.0 );
            s.material.specular_material = Vec3( 0.5,0.6,1.0 );
            s.material.shininess = 16;
        }

        { //Ceiling
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(0.3,0.,0.2 );
            s.material.specular_material = Vec3( 0.3,0.,0.2 );
            s.material.shininess = 16;
        }


        {
            meshes.resize(meshes.size() + 1);
            Mesh & mesh = meshes[meshes.size() - 1];
            mesh.loadOFF("star.off");
            mesh.build_arrays();
            mesh.material.ambient_material = Vec3(0.0, 0.02, 0.05);
            mesh.material.diffuse_material = Vec3(0.8, 0.2, 0.);
            mesh.material.specular_material = Vec3(0.1, 0.4, 0.);
            mesh.material.shininess = 16;
            mesh.centerAndScaleToUnit();
            mesh.translate(Vec3(0., -0.7, 0.));
            mesh.scale(Vec3(0.5, 0.5, 0.5));
        }
        updateSceneObj();

    }

    void add_square(float z_pos) {
        squares.resize( squares.size() + 1 );
        Square & s = squares[squares.size() - 1];
        s.setQuad(Vec3(-1., -1., z_pos), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
        s.build_arrays();
        s.material.diffuse_material = Vec3( 1,0.8,0.8 );
        s.material.specular_material = Vec3( 0.8,0.8,0.8 );
        s.material.shininess = 20;
    }

    void setup_mesh(std::string meshName) {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();
        sceneObj.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(5., 5., 5.);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }

        {
            meshes.resize(meshes.size() + 1);
            Mesh & mesh = meshes[meshes.size() - 1];
            mesh.loadOFF(meshName);
            mesh.build_arrays();
            mesh.material.ambient_material = Vec3(0.07, 0.02, 0.);
            mesh.material.diffuse_material = Vec3(0.8, 0.5, 0.3);
            mesh.material.specular_material = Vec3(0.1, 0.4, 0.2);
            mesh.material.shininess = 16;
        }

        updateSceneObj();



    }
};



#endif
