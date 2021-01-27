#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <algorithm>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <random>
static std::default_random_engine engine(10); // random seed=10
static std::uniform_real_distribution<double> uniform(1.,0);

class Vector{
    public:
        explicit Vector(double x=0, double y=0, double z=0) { 
            coords[0] = x;
            coords[1] = y;
            coords[2] = z;
        };
        double operator[](int i) const {return coords[i];};
        double &operator[](int i) {return coords[i];};
        double sqrNorm() const {
            return( coords[0] * coords[0]+ coords[1] * coords[1]+ coords[2] * coords[2] );
        };
        const Vector& operator+=(const Vector& a){
            coords[0] += a[0];
            coords[1] += a[1];
            coords[2] += a[2];
            return *this;
        }
        Vector get_normalized(){
            double norm = sqrt(sqrNorm());
            return(Vector(coords[0]/norm, coords[1]/norm, coords[2]/norm));
        };

    private:
        double coords [3];
};
// Surcharges d'Opérateur 
Vector operator+ (const Vector &a, const Vector &b){
    return(Vector(a[0]+b[0],a[1]+b[1],a[2]+b[2]));
};
Vector operator- (const Vector &a, const Vector &b){
    return(Vector(a[0]-b[0],a[1]-b[1],a[2]-b[2]));
};
Vector operator- (const Vector &a){
    return(Vector(-a[0],-a[1],-a[2]));
};
Vector operator* (const double &a, const Vector &b){
    return(Vector(a*b[0],a*b[1],a*b[2]));
};
Vector operator* (const Vector &a, const double &b){
    return(Vector(a[0]*b,a[1]*b,a[2]*b));
};
Vector operator/ (const Vector &a, const double &b){
    return(Vector(a[0]/b,a[1]/b,a[2]/b));
};
Vector cross(const Vector &a, const Vector &b){ // Mat product
    return Vector(
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    );
};
Vector operator* (const Vector &a, const Vector &b){ // term-to-term multiplication
    return(Vector(a[0]*b[0],a[1]*b[1],a[2]*b[2]));
};

double dot(const Vector &a, const Vector &b){ // Scalar product
    return(a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
};
double sqr(double x){ // Compute square
    return x*x;
}

Vector random_cos(const Vector &N){
    double u1 = uniform(engine);
    double u2 = uniform(engine);

    double x = cos(2 * M_PI*u1) * sqrt(1 - u2);
    double y = sin(2 * M_PI*u1) * sqrt(1 - u2);
    double z = sqrt(u2);
    Vector T1;

    if(N[0] < N[1] && N[0] < N[2]){ // Composante sur x est la plus petite
        T1 = Vector(0,N[2], -N[1]);
    }
    else{
        if(N[1] < N[0] && N[1] < N[2]){ // Composante sur y est la plus petite
            T1 = Vector(N[2],0,-N[0]);
        }
        else{ // Composante sur z est la plus petite
            T1 = Vector(N[1], -N[0],0);
        }
        

    }
    T1 = T1.get_normalized();
    Vector T2 = cross(N,T1);
    return (x*T1 + y*T2 + z*N);
}
//////////////// Github pour le rendu /////////////////////////
// Image des différents trucs qui marchent avec explications
// Section feedback sur le cours à la fin du rapport
class Ray{
    public:
        Ray(const Vector& C, const Vector&u) : C(C), u(u) {}
        Vector C, u;
};


class Sphere{
public:
    Sphere(const Vector& O, double R, const Vector &albedo, bool isMirror = false, bool isTransparent = false): O(O), R(R), albedo(albedo), isMirror(isMirror), isTransparent(isTransparent) {
    }
    bool intersect(const Ray& r, Vector& P, Vector& N, double &t){ 
        // Solves the intersection equation
        double a = 1;
        double b = 2 * dot(r.u, r.C - O);
        double c = (r.C-O).sqrNorm() - R*R;

        double delta = b*b - 4*a*c;
        if (delta < 0) {
            t = 1E10;
            return false;
        }
        else{
            double sqrtD(sqrt(delta));
            double t2 = (-b + sqrtD) / (2*a);

            if(t2 < 0)return false;

            double t1 = (-b - sqrtD) / (2*a);

            if(t1 > 0){
                t = t1;
            }
            else{
                t = t2;
            }

            P = r.C + (r.u * t);
            N = (P - O).get_normalized();
            return true;
        };
    }
    Vector O;
    double R;  
    Vector albedo;
    bool isMirror, isTransparent;     
};

class Scene{
public:
    Scene(){};
    bool intersect(const Ray& r, Vector& P, Vector& N, Vector &albedo, bool &mirror, bool &transp, double &t){
        t = 1E10;
        bool is_inter = false;

        for(int i=0; i < objects.size(); i++){
            Vector loc_P, loc_N;
            double loc_T;
            if(objects[i].intersect(r, loc_P, loc_N, loc_T) && loc_T < t){
                t = loc_T;
                is_inter = true;
                albedo = objects[i].albedo;
                P = loc_P;
                N = loc_N;
                mirror = objects[i].isMirror;
                transp = objects[i].isTransparent;
            };
        };
        return is_inter;
    };

    Vector getColor(const Ray& r, int rebond){
        Vector P, N, albedo;
        double t;
        bool mirror, transp;

        bool inter = intersect(r, P, N, albedo, mirror, transp, t); 

        Vector p_color(0,0,0); 

        if(rebond >5){
            return Vector(0.,0.,0.);
        }
        if(inter){ 
            if (mirror){ // Miroir
                Vector reflectedDir = r.u - 2* dot(r.u, N)*N; // Reflected Ray Direction
                Ray reflectedRay(P + 0.001*N, reflectedDir);

                return getColor(reflectedRay, rebond + 1);
            }
            else{
                if (transp){ // Transp
                    // Refractive index (n1,n2) & Normal Vector (N2)
                    double n1 = 1, n2 = 1.4;  
                    Vector N2 = N;

                    if(dot(r.u,N) > 0){ // Leaving sphere, change normal vector and swap refractive indexes
                        std::swap(n1, n2);
                        N2 = -N;
                    };
                                   
                    // Tangential (Tt) & normal (Tn) vector, incidence index (rad)
                    Vector Tt = n1/n2 * (r.u - dot(r.u, N2)*N2);
                    double rad = 1 - sqr(n1/n2) * (1 - sqr(dot(r.u,N2))); 

                    if(rad < 0){ // No transmission, TOTAL REFLECTION
                        Vector reflectedDir = r.u - 2* dot(r.u, N)*N;
                        Ray reflectedRay(P + 0.0001*N, reflectedDir);
                        return getColor(reflectedRay, rebond + 1);
                    }

                    Vector Tn = -sqrt(rad)*N2; // /!\ Signe
                    Vector refraDir = Tt + Tn;
                    Ray refraRay(P - 0.0001*N2, refraDir);
                    return getColor(refraRay, rebond + 1); // Enlever +1
                }
                else{
                    // Eclairage direct
                    Vector PL = (L - P);
                    double norm = sqrt(PL.sqrNorm()); 
                    
                    // Shadow
                    Vector shadowP, shadowN, shadowAlbedo;
                    double shadowt;
                    bool shadowMirror, shadowTransp;

                    Ray shadowRay(P + 0.0001 * N, PL/norm);
                    
                    bool shadowInter = intersect(shadowRay, shadowP, shadowN, shadowAlbedo, shadowMirror, shadowTransp, shadowt);

                    if(shadowInter && shadowt < norm){
                        p_color = Vector(0.,0.,0.);
                    } 
                    else{
                        p_color =  I / (4 * M_PI * norm*norm) * albedo/M_PI * std::max(0., dot(N, PL/norm));
                    }

                    // Eclairage indirect
                    Vector wi = random_cos(N); // Reflected Ray Direction
                    Ray wiRay(P + 0.0001*N, wi);
                    p_color += albedo*getColor(wiRay, rebond + 1);

                }     
            }  
        }
        return p_color;
    };

    std::vector<Sphere> objects;
    Vector L;
    double I;
};

int main() {
    // Define image size
    int W = 512;
    int H = 512;
    
    // Create Scene and light. Setup camera position
    Vector C(0,0,55);
    Scene scene;
    scene.I = 5E9;
    scene.L = Vector(-10,20,40);

    // Add elemnts to the scene (Walls approximated with spheres)
    Sphere S(Vector(0,0,0),10,Vector(1,1.,1.),false,true); // White
    Sphere SmurBot(Vector(0, -1000, 0),990,Vector(1, 1., 1.)); // White
    Sphere SmurTop(Vector(0, 1000, 0),940,Vector(1, 1., 1.)); // White
    Sphere SmurLeft(Vector(-1000, 0 ,0),940,Vector(1, 0.5, 0.5)); // Red
    Sphere SmurRight(Vector(1000, 0, 0),940,Vector(0.5, 1., 0.5));// Green
    Sphere SmurBack(Vector(0, 0, -1000),940,Vector(0.2, 1., 1.));// Cyan
    Sphere SmurFront(Vector(0, 0, 1000),940,Vector(1, 1., 0.5)); // Yellow
    
    scene.objects.push_back(S);
    scene.objects.push_back(SmurLeft);
    scene.objects.push_back(SmurRight);
    scene.objects.push_back(SmurBack);
    scene.objects.push_back(SmurFront);
    scene.objects.push_back(SmurTop);
    scene.objects.push_back(SmurBot);

    double angle = 60;
    double f = angle * M_PI / 180;
    int nb_rays = 71;

    std::vector<unsigned char> image(W*H * 3, 0);
    #pragma omp parallel for schedule(dynamic,1) //+ regarder prise en charge openMP avec vscode ?
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector u(j-W/2, i-H/2, -W/(2.*tan(f/2)));
            u = u.get_normalized();
            Ray r(C, u);

            // Complete here
            Vector p_color(0,0,0);

            for(int k=0; k < nb_rays; k++){
                p_color += scene.getColor(r,0);
            };
            p_color = p_color / nb_rays;
            // Till here
            image[((H-i-1)*W + j) * 3 + 0] = std::min(255.,std::pow(p_color[0],0.45));
            image[((H-i-1)*W + j) * 3 + 1] = std::min(255.,std::pow(p_color[1],0.45));
            image[((H-i-1)*W + j) * 3 + 2] = std::min(255.,std::pow(p_color[2],0.45));
        }
    }
    stbi_write_png("Images/image_71Ray.png", W, H, 3, &image[0], 0);
 
    return 0;
}
// cd "d:\Documents\CentralePougne\3A\MOS\Informatique_Graphique\" ; if ($?) { g++ -fopenmp test.cpp -o test } ; if ($?) { .\test }
