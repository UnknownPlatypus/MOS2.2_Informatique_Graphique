#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <algorithm>

#include <string> 
#include <iostream>
#include <stdio.h>
#include <list>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <random>
static std::default_random_engine engine(10); // random seed=10
static std::uniform_real_distribution<double> uniform(1.,0);

//////////////// Github pour le rendu /////////////////////////
// Image des différents trucs qui marchent avec explications
// Section feedback sur le cours à la fin du rapport

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
            return(coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2]);
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
    return(Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]));
};
Vector operator- (const Vector &a, const Vector &b){
    return(Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]));
};
Vector operator- (const Vector &a){
    return(Vector(-a[0], -a[1], -a[2]));
};
Vector operator* (const double &a, const Vector &b){
    return(Vector(a * b[0], a * b[1], a * b[2]));
};
Vector operator* (const Vector &a, const double &b){
    return(Vector(a[0] * b, a[1] * b, a[2] * b));
};
Vector operator/ (const Vector &a, const double &b){
    return(Vector(a[0] / b, a[1] / b, a[2] / b));
};
Vector cross(const Vector &a, const Vector &b){ // Mat product
    return Vector(
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    );
};
Vector operator* (const Vector &a, const Vector &b){ // term-to-term multiplication
    return(Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]));
};

double dot(const Vector &a, const Vector &b){ // Scalar product
    return(a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
};
double sqr(double x){ // Square
    return x*x;
}

Vector random_cos(const Vector &N){
    double u1 = uniform(engine);
    double u2 = uniform(engine);

    double x = cos(2 * M_PI * u1) * sqrt(1 - u2);
    double y = sin(2 * M_PI * u1) * sqrt(1 - u2);
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

class Ray{
public:
        Ray(const Vector& C, const Vector&u) : C(C), u(u) {}
        Vector C, u;
};

class Object{
public:
    Object(){};
    virtual bool intersect(const Ray& r, Vector& P, Vector& normale, double &t) = 0;
    Vector albedo;
    bool isMirror, isTransparent;    
};

class Sphere : public Object {
public:
    Sphere(const Vector& O, double R, const Vector &albedo, bool isMirror = false, bool isTransparent = false): O(O), R(R){
        this->albedo = albedo;
        this->isMirror = isMirror;
        this->isTransparent = isTransparent;
    }
    bool intersect(const Ray& r, Vector& P, Vector& N, double &t){ 
        // Solves the intersection equation
        double a = 1;
        double b = 2 * dot(r.u, r.C - O);
        double c = (r.C - O).sqrNorm() - R * R;

        double delta = b * b - 4 * a * c;
        if (delta < 0) {
            t = 1E10;
            return false;
        }
        else{
            double sqrtD = sqrt(delta);
            double t2 = (-b + sqrtD) / (2 * a);

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
 
};

class Boundingbox{
public:
    bool intersect(const Ray& r){
        double tx1 = (mini[0] - r.C[0]) / r.u[0];
        double tx2 = (maxi[0] - r.C[0]) / r.u[0];
        double txmin = std::min(tx1, tx2);
        double txmax = std::max(tx1, tx2);

        double ty1 = (mini[1] - r.C[1]) / r.u[1];
        double ty2 = (maxi[1] - r.C[1]) / r.u[1];
        double tymin = std::min(ty1, ty2);
        double tymax = std::max(ty1, ty2);

        double tz1 = (mini[2] - r.C[2]) / r.u[2];
        double tz2 = (maxi[2] - r.C[2]) / r.u[2];
        double tzmin = std::min(tz1, tz2);
        double tzmax = std::max(tz1, tz2);

        double tMax = std::min(txmax,std::min(tymax,tzmax));
        double tMin = std::max(txmin,std::min(tymin,tzmin));

        if (tMax <0) return false;
        return tMax > tMin;
    }
    Vector mini, maxi;
};


class Scene{
public:
    Scene(){};
    bool intersect(const Ray& r, Vector& P, Vector& N, Vector &albedo, bool &mirror, bool &transp, double &t, int& objectid){
        t = 1E10; // Infini
        bool is_inter = false; // Init

        for(int i=0; i < objects.size(); i++){
            Vector loc_P, loc_N;
            double loc_T;

            if(objects[i]->intersect(r, loc_P, loc_N, loc_T) && loc_T < t){
                t = loc_T;
                is_inter = true;
                albedo = objects[i]->albedo;
                P = loc_P;
                N = loc_N;
                mirror = objects[i]->isMirror;
                transp = objects[i]->isTransparent;
                objectid = 1;
            };
        };
        return is_inter;
    };

    Vector getColor(const Ray& r, int rebond, bool lastdefuse){

        double epsi = 0.00001; // Deal with intersection issues
        if(rebond >5){
            return Vector(0.,0.,0.);
        }
        // Initialize
        Vector P, N, albedo;
        double t;
        bool mirror, transp;
        int objectid;
        bool inter = intersect(r, P, N, albedo, mirror, transp, t, objectid); 
        Vector p_color(0,0,0); 
        
        if(inter){ 
            if(objectid == 0){
                if (rebond == 0 || !lastdefuse) {
                    return Vector(I,I,I) / (4 * M_PI * M_PI * dynamic_cast<Sphere*>(objects[0])->R * dynamic_cast<Sphere*>(objects[0])->R);
                }
                else{
                    return Vector(0, 0, 0);
                }
            }

            else{            
                if (mirror){ // Mirror effect
                    Vector reflectedDir = r.u - 2 * dot(r.u, N) * N; // Reflected Ray Direction
                    Ray reflectedRay(P + epsi * N, reflectedDir);

                    return getColor(reflectedRay, rebond + 1, false);
                }
                else{
                    if (transp){ // Transparent effect
                        // Refractive index (n1,n2) & Normal Vector (N2)
                        double n1 = 1, n2 = 1.4;  //n2 ~= 1.4 for glass
                        Vector N2 = N;

                        if(dot(r.u, N) > 0){ // Leaving sphere, invert normal vector and swap refractive indexes
                            std::swap(n1, n2);
                            N2 = -N;
                        };
                                    
                        // Tangential (Tt) & Normal (Tn) vectors, incidence index (rad)
                        Vector Tt = n1/n2 * (r.u - dot(r.u, N2) * N2);
                        double rad = 1 - sqr(n1 / n2) * (1 - sqr(dot(r.u,N2))); 

                        if(rad < 0){ // TOTAL REFLECTION (No transmission)
                            Vector reflectedDir = r.u - 2 * dot(r.u, N) * N;
                            Ray reflectedRay(P + epsi * N, reflectedDir);
                            return getColor(reflectedRay, rebond + 1, false);
                        }

                        Vector Tn = -sqrt(rad) * N2; // /!\ Signe
                        Vector refraDir = Tt + Tn;
                        Ray refraRay(P - epsi*N2, refraDir);
                        return getColor(refraRay, rebond, false); // Enlever +1 dans rebond ?
                    }
                    else{
                        // Vector PL = (L - P);
                        // double norm = sqrt(PL.sqrNorm()); 
                        
                        // // Shadow
                        // Vector shadowP, shadowN, shadowAlbedo;
                        // double shadowt;
                        // bool shadowMirror, shadowTransp;
                        // int shadowid;

                        // Ray shadowRay(P + 0.0001 * N, PL/norm);
                        
                        // bool shadowInter = intersect(shadowRay, shadowP, shadowN, shadowAlbedo, shadowMirror, shadowTransp, shadowt, shadowid);

                        // if(shadowInter && shadowt < norm){
                        //     p_color = Vector(0.,0.,0.);
                        // } 
                        // else{
                        //     p_color =  I / (4 * M_PI * norm*norm) * albedo/M_PI * std::max(0., dot(N, PL/norm));
                        // }

                        // Eclairage direct
                        Vector PL = (L - P);
                        PL = PL.get_normalized();
                        Vector w = random_cos(-PL);
                        Vector xprim = w * dynamic_cast<Sphere*>(objects[0])->R + dynamic_cast<Sphere*>(objects[0])->O;
                        Vector Pxprime = xprim - P;
                        double d = sqrt(Pxprime.sqrNorm());
                        Pxprime = Pxprime/d;

                        Vector shadowP, shadowN, shadowAlbedo;
                        double shadowt;
                        bool shadowMirror, shadowTransp;
                        int shadowid;

                        Ray shadowRay(P + epsi * N, Pxprime);
                        
                        bool shadowInter = intersect(shadowRay, shadowP, shadowN, shadowAlbedo, shadowMirror, shadowTransp, shadowt, shadowid);

                        if(shadowInter && shadowt < d - 3*epsi){ // Might change to 3-4*epsi
                            p_color = Vector(0.,0.,0.);
                        } 
                        else{
                            double SR = sqr(dynamic_cast<Sphere*>(objects[0])->R);
                            double proba = std::max(0.,dot(-PL,w)) / (M_PI * SR * SR);
                            double J = std::max(0.,dot(w, -Pxprime) / (d*d));
                            p_color =  I / (4 * M_PI * SR * SR) * albedo  * std::max(0., dot(N, Pxprime)) * J / proba; 
                        }

                        // Eclairage indirect
                        Vector dir_ir = random_cos(N); // Reflected Ray Direction
                        Ray irRay(P + epsi * N, dir_ir.get_normalized()); // normalize ?
                        p_color += albedo / M_PI * getColor(irRay, rebond + 1, true); //  /M_PI ?

                    }     
                }  
            }
        }
        return p_color;
    };

    std::vector<Object*> objects;
    Vector L;
    double I;
};


class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
	};
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};

class Noeud {
public:
        Noeud *fg, *fd;
        Boundingbox b;
        int start, end;
};

class TriangleMesh : public Object {
public:
  ~TriangleMesh() {}
	TriangleMesh(const Vector& albedo, bool mirror = false, bool transp = false) {
        this->albedo = albedo;
        isMirror = mirror;
        transp = transp;
        BVH = new Noeud;
    };

    Boundingbox buildBB(int start, int end){// triangle indexes
        bool inf = 1E9;
        bb.mini = Vector(inf,inf,inf);
        bb.maxi = Vector(-inf,-inf,-inf);

        for(int i=start; i < end; i++){
            for(int j=0; j < 3; j++){
                bb.mini[j] = std::min(bb.mini[j], vertices[indices[i].vtxi][j]); 
                bb.maxi[j] = std::max(bb.maxi[j], vertices[indices[i].vtxi][j]); 
                bb.mini[j] = std::min(bb.mini[j], vertices[indices[i].vtxj][j]); 
                bb.maxi[j] = std::max(bb.maxi[j], vertices[indices[i].vtxj][j]); 
                bb.mini[j] = std::min(bb.mini[j], vertices[indices[i].vtxk][j]); 
                bb.maxi[j] = std::max(bb.maxi[j], vertices[indices[i].vtxk][j]); 
            }            
        }
    }

void buildBVH(Noeud* n, int start, int end){
       
        n->start = start;
        n->end = end;
        n->b = buildBB(n->start, n->end);
        Vector diag = n->b.maxi - n->b.mini;

        int d;
        if (diag[0]>=diag[1] && diag[0]>=diag[2]){
            d = 0;
        }
        else{
            if (diag[1]>=diag[0] && diag[1]>=diag[2]){
                d = 1;
            }
            else{
                d = 2;
            }
        }

        double mid = (n->b.mini[d] + n->b.maxi[d]) / 2;
        int id_pivot = n-> start;

        for (int i = n->start; i < n->end; i++){
            double mid_tri = (vertices[indices[i].vtxi][d] + vertices[indices[i].vtxj][d] + vertices[indices[i].vtxk][d]) / 3;
            if(mid_tri < mid){
                std::swap(indices[i], indices[id_pivot]);
                id_pivot++;
            }
        }

        n->fg = NULL;
        n->fd = NULL;
        if(id_pivot == start || id_pivot == end || (end-start < 5)) return;

        n->fg = new Noeud;
        n->fd = new Noeud;
        buildBVH(n->fg, n->start, id_pivot);
        buildBVH(n->fd, id_pivot,n->end);
    }
	
	void readOBJ(const char* obj) {

		char matfile[255];
		char grp[255];

		FILE* f;
		f = fopen(obj, "r");
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				curGroup++;
			}

			if (line[0] == 'v' && line[1] == ' ') {
				Vector vec;

				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
					col[0] = std::min(1., std::max(0., col[0]));
					col[1] = std::min(1., std::max(0., col[1]));
					col[2] = std::min(1., std::max(0., col[2]));

					vertices.push_back(vec);
					vertexcolors.push_back(col);

				} else {
					sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector vec;
				sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;
				t.group = curGroup;

				char* consumedline = line + 1;
				int offset;

				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
					indices.push_back(t);
				} else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
						indices.push_back(t);
					} else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							indices.push_back(t);
						} else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}

				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.group = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					} else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						} else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;								
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							} else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								} else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}

			}

		}
		fclose(f);

	}

    bool intersect(const Ray& r, Vector& P, Vector& normale, double &t){ 
        if(!BVH->b.intersect(r)) return false;

        t = 1E10; // double t
        bool is_inter = false;

        std::list<Noeud*> l;
        l.push_back(BVH);

        while(!l.empty()){
            Noeud* c = l.front();
            l.pop_front();
            if(c->fg){
                if(!c->fg->b.intersect(r)){
                    l.push_front(c->fg);
                }
                if(!c->fd->b.intersect(r)){
                    l.push_front(c->fd);
                }
            }
            else{
                for(int i = c->start; i < c->end; i++){
                    const Vector &A = vertices[indices[i].vtxi]; // First vertice
                    const Vector &B = vertices[indices[i].vtxj]; // Second vertice
                    const Vector &C = vertices[indices[i].vtxk]; // Third vertice

                    Vector e1 = B - A;
                    Vector e2 = C - A;
                    Vector N = cross(e1,e2);
                    Vector AO = r.C - A;
                    Vector AOu = cross(AO, r.u);
                    double invUN = 1. / dot(r.u, N);

                    double beta = -dot(e2,AOu) * invUN;
                    double gamma = dot(e1,AOu) * invUN;
                    double alpha = 1 - beta - gamma;
                    double localt = -dot(AO,N) * invUN;

                    if(beta >= 0 && gamma >= 0 && beta <=1 && gamma <= 1 && alpha >= 0 && localt > 0){ // intesection
                        is_inter = true;
                        if(localt < t){
                            t = localt;
                            // normale = N.get_normalized(); // Modifier pour stocker index utile et calculer a la fin seulement
                            normale = alpha*normals[indices[i].ni] + beta*normals[indices[i].nj] + gamma*normals[indices[i].nk];
                            normale = normale.get_normalized();
                            P = r.C  + r.u *t;
                        }
                    }
                }

            }
            
        }       
        return is_inter;        
    }

	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
    Boundingbox bb; 

    Noeud* BVH;
};



int main() {
    // Define image size & nb of Rays
    int W = 128; // 512
    int H = 128; // 512
    int nb_rays = 2; // 100

    // Create Scene and light. Setup camera position
    Vector C(0,0,55);
    Scene scene;
    scene.I = 5E9;
    scene.L = Vector(-10, 20, 40);
    double angle = 60;
    double alpha = angle * M_PI / 180;

    // Add elemnts to the scene (Walls approximated with spheres)
    Sphere Lumiere(scene.L, 5, Vector(1.,1.,1.)); // White
    //Sphere S(Vector(0,0,0),10,Vector(1.,1.,1.),false,false); // White
    Sphere S1(Vector(0, 0, 0), 10, Vector(1, 1, 1), false, false);// White
    Sphere S2(Vector(-10, 0, -20), 10, Vector(1, 1, 1), false, false);// White
    Sphere S3(Vector(10, 0, 20), 10, Vector(1, 1, 1), false, false);// White
    Sphere SmurBot(Vector(0, -1000, 0),990,Vector(0., 0., 1.)); // Blue
    Sphere SmurTop(Vector(0, 1000, 0),940,Vector(1., 0., 0.)); // Red
    Sphere SmurLeft(Vector(-1000, 0 ,0),960,Vector(1., 0.5, 0.5)); // Pale Red
    Sphere SmurRight(Vector(1000, 0, 0),960,Vector(0.5, 1., 0.5));// Green
    Sphere SmurBack(Vector(0, 0, -1000),940,Vector(0.2, 1., 1.));// Cyan
    Sphere SmurFront(Vector(0, 0, 1000),940,Vector(1., 1., 0.5)); // Yellow
    TriangleMesh m(Vector(1.,1.,1.));

    m.readOBJ("Images/Dog/13463_Australian_Cattle_Dog_v3.obj");
    for(int i=0; i < m.vertices.size(); i++){
        std::swap(m.vertices[i][1], m.vertices[i][2]);
        m.vertices[i][2] = -m.vertices[i][2];
        m.vertices[i][2] += 10;
        m.vertices[i][1] -= 10;
        
    }
    for(int i=0; i < m.normals.size(); i++){
        std::swap(m.normals[i][1], m.normals[i][2]);
    }
    //m.buildBB(); // complete ?
    m.buildBVH(m.BVH, 0, m.indices.size()); 

    scene.objects.push_back(&Lumiere); 
    //scene.objects.push_back(S);
    // scene.objects.push_back(&S1);
    // scene.objects.push_back(&S2);
    // scene.objects.push_back(&S3);
    scene.objects.push_back(&m);
    scene.objects.push_back(&SmurLeft);
    scene.objects.push_back(&SmurRight);
    scene.objects.push_back(&SmurBack);
    scene.objects.push_back(&SmurFront);
    scene.objects.push_back(&SmurTop);
    scene.objects.push_back(&SmurBot);

    // Depth of field
    double capSize = 0.00001;
    double focale = 55;

    std::vector<unsigned char> image(W * H * 3, 0);
    #pragma omp parallel for schedule(dynamic,1) 
    for (int i = 0; i < H; i++) { // Height
        for (int j = 0; j < W; j++) { // Width

            Vector p_color(0,0,0);

            for(int k=0; k < nb_rays; k++){ // Rays
                
                double u1 = uniform(engine);
                double u2 = uniform(engine);
                double dx = 0.25 * cos(2 * M_PI * u1) * sqrt(-2 * log(u2));
                double dy = 0.25 * sin(2 * M_PI * u1) * sqrt(-2 * log(u2));

                u1 = uniform(engine);
                u2 = uniform(engine);
                double dx2 = capSize * cos(2 * M_PI * u1) * sqrt(-2 * log(u2));
                double dy2 = capSize * sin(2 * M_PI * u1) * sqrt(-2 * log(u2));

                Vector dir(j - W/2 + dx +0.5, i - H/2 + dy + 0.5, -W / (2. * tan(alpha/2)));  
                dir = dir.get_normalized();

                Vector target = C + focale * dir;
                Vector C_prime = C + Vector(dx2, dy2, 0);
                Vector dir_prime = (target - C_prime).get_normalized();

                Ray r(C_prime,dir_prime);
                
                p_color += scene.getColor(r, 0, false);
            };
            p_color = p_color / nb_rays;
            // Till here
            image[((H - i - 1) * W + j) * 3 + 0] = std::min(255., std::pow(p_color[0],0.45));
            image[((H - i - 1) * W + j) * 3 + 1] = std::min(255., std::pow(p_color[1],0.45));
            image[((H - i - 1) * W + j) * 3 + 2] = std::min(255., std::pow(p_color[2],0.45));
        }
    }
    stbi_write_png("Images/image_Mesh2.png", W, H, 3, &image[0], 0);
 
    return 0;
}
// cd "d:\Documents\CentralePougne\3A\MOS\Informatique_Graphique\" ; if ($?) { g++ -fopenmp test.cpp -o test } ; if ($?) { .\test }