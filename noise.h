#ifndef NOISE
#define NOISE

#ifndef GLEW_STATIC
#define GLEW_STATIC
#include "GL/glew.h"
#include <glm/glm.hpp>
#include <glm/ext.hpp>
#endif

struct Noise3D{
	~Noise3D(){}
	virtual float getValueAt(float xin, float yin, float zin){ return 0; }
	protected:
		glm::ivec3 size;
		glm::vec3 scale;
};

struct FractalNoise3D{
	
	Noise3D **layer;
	int layers; // octaves
	float persistence;
	
	FractalNoise3D(int n, float p){ // fetch layers of noise
		layers = n;
		layer = (Noise3D**)calloc(sizeof(Noise3D*), layers);
		persistence = p;
	}
	~FractalNoise3D(){
		for(int i = 0; i < layers; i++){
			if(layer[i] != NULL)
				delete layer[i];
		}
		free(layer);
	}
	void addLayer(Noise3D *l, int i){
		layer[i] = l;
	}
	float getValueAt(float xin, float yin, float zin){
		float total = 0;
		float maxValue = 0;
		float frequency = 1;
		float amplitude = 1;
		for(int i = 0; i < layers; i++){
			total += amplitude * layer[i]->getValueAt(xin * frequency, yin * frequency, zin * frequency);
			maxValue += amplitude;
			amplitude *= persistence;
			frequency *= 2;
		}
		return total / maxValue;
	}
};

// Perlin Noise in 3 dimensions, using 3 different methods

struct PerlinNoise3D : Noise3D{
	
	// external usage
	PerlinNoise3D(glm::ivec3 pn, glm::vec3 s);
	virtual float getValueAt(float xin, float yin, float zin);
	
	// general internal functionality
	protected:
		void getCoordinates(glm::vec3 p, glm::ivec3 *i0, glm::ivec3 *i1, glm::vec3 *t, glm::vec3 *d0, glm::vec3 *d1);
		void getCornerCoordinates(glm::ivec3 *c, glm::ivec3 i0, glm::ivec3 i1);
		void getDotProducts(float *n, glm::vec3 d0, glm::vec3 d1, glm::vec3 *g);
		float interpolateValues(float a, float b, float w);
		int getPermutation(glm::ivec3 c);
		
		// later-updated functionality
		virtual void getInterpolation(glm::vec3 t, glm::vec3* s);
		virtual void getGradientProducts(float *n, glm::ivec3 i0, glm::ivec3 i1, glm::vec3 d0, glm::vec3 d1);
		virtual void getGradientVectors(glm::vec3 *g, glm::ivec3 i0, glm::ivec3 i1){}
};

struct PerlinNoise3DOriginal : PerlinNoise3D{
	
	// initialise random gradient lookup table for grid intersections
	PerlinNoise3DOriginal(glm::ivec3 pn, glm::vec3 s);
	~PerlinNoise3DOriginal();
	protected:
		void getGradientVectors(glm::vec3 *g, glm::ivec3 i0, glm::ivec3 i1);
	private:
		glm::vec3 *gradient;
		int gradientn;
};

struct PerlinNoise3DImproved : PerlinNoise3D{
	PerlinNoise3DImproved(glm::ivec3 pn, glm::vec3 s) : PerlinNoise3D(pn, s){}
	protected:
		
		// replace random gradient distribution with predetermined axis-aligned 12 gradient vectors
		virtual void getGradientVectors(glm::vec3 *g, glm::ivec3 i0, glm::ivec3 i1);
		
		// smootherstep, or "quintic", weights
		virtual void getInterpolation(glm::vec3 t, glm::vec3* s);
};

struct PerlinNoise3DImplicit : PerlinNoise3DImproved{
	PerlinNoise3DImplicit(glm::ivec3 pn, glm::vec3 s) : PerlinNoise3DImproved(pn, s){}
	
	// replace gradient vector lookup and dot product with optimised implicit calculation
	protected:
		virtual void getGradientProducts(float *n, glm::ivec3 i0, glm::ivec3 i1, glm::vec3 d0, glm::vec3 d1);
	private:
		float getGradientProductImplicit(int p, float x, float y, float z);
};

// main functions

PerlinNoise3D::PerlinNoise3D(glm::ivec3 pn, glm::vec3 s){
	size = pn;
	scale = s;
}

float PerlinNoise3D::getValueAt(float xin, float yin, float zin){ // retrieve single noise positional value
	
	// coordinates
	glm::vec3 p = glm::vec3(xin / scale.x, yin / scale.y, zin / scale.z);
	glm::ivec3 i0, i1;
	glm::vec3 t, d0, d1;
	getCoordinates(p, &i0, &i1, &t, &d0, &d1);
	
	// interpolation
	glm::vec3 s;
	getInterpolation(t, &s);
	
	// gradient dot products
	float n[0b1000];
	getGradientProducts((float*)&n, i0, i1, d0, d1);
	
	// interpolation
	float v00, v10, v01, v11, v0, v1;
	v00 = interpolateValues(n[0b000], n[0b100], s.x);
	v10 = interpolateValues(n[0b010], n[0b110], s.x);
	v01 = interpolateValues(n[0b001], n[0b101], s.x);
	v11 = interpolateValues(n[0b011], n[0b111], s.x);
	v0 = interpolateValues(v00, v10, s.y);
	v1 = interpolateValues(v01, v11, s.y);
	return interpolateValues(v0, v1, s.z);
}

// internal functions

void PerlinNoise3D::getCoordinates(glm::vec3 p, glm::ivec3 *i0, glm::ivec3 *i1, glm::vec3 *t, glm::vec3 *d0, glm::vec3 *d1){
	
	// grid-cell coordinates
	*i0 = glm::ivec3((int)p.x & 255, (int)p.y & 255, (int)p.z & 255);
	*i1 = glm::ivec3((i0->x + 1) & 255, (i0->y + 1) & 255, (i0->z + 1) & 255);
	
	// grid-cell distances
	*t = glm::vec3(p.x - (float)i0->x, p.y - (float)i0->y, p.z - (float)i0->z);
	
	// grid-cell distance vectors
	*d0 = *t;
	*d1 = *t - glm::vec3(1);
}

void PerlinNoise3D::getCornerCoordinates(glm::ivec3 *c, glm::ivec3 i0, glm::ivec3 i1){
	c[0b000] = glm::ivec3(i0.x, i0.y, i0.z);
	c[0b001] = glm::ivec3(i1.x, i0.y, i0.z);
	c[0b010] = glm::ivec3(i0.x, i1.y, i0.z);
	c[0b011] = glm::ivec3(i1.x, i1.y, i0.z);
	c[0b100] = glm::ivec3(i0.x, i0.y, i1.z);
	c[0b101] = glm::ivec3(i1.x, i0.y, i1.z);
	c[0b110] = glm::ivec3(i0.x, i1.y, i1.z);
	c[0b111] = glm::ivec3(i1.x, i1.y, i1.z);
}

void PerlinNoise3D::getDotProducts(float *n, glm::vec3 d0, glm::vec3 d1, glm::vec3 *g){
	glm::vec3 d[0b1000];
	d[0b000] = glm::vec3(d0.x, d0.y, d0.z);
	d[0b100] = glm::vec3(d1.x, d0.y, d0.z);
	d[0b010] = glm::vec3(d0.x, d1.y, d0.z);
	d[0b110] = glm::vec3(d1.x, d1.y, d0.z);
	d[0b001] = glm::vec3(d0.x, d0.y, d1.z);
	d[0b101] = glm::vec3(d1.x, d0.y, d1.z);
	d[0b011] = glm::vec3(d0.x, d1.y, d1.z);
	d[0b111] = glm::vec3(d1.x, d1.y, d1.z);
	n[0b000] = glm::dot(d[0b000], g[0b000]);
	n[0b100] = glm::dot(d[0b100], g[0b100]);
	n[0b010] = glm::dot(d[0b010], g[0b010]);
	n[0b110] = glm::dot(d[0b110], g[0b110]);
	n[0b001] = glm::dot(d[0b001], g[0b001]);
	n[0b101] = glm::dot(d[0b101], g[0b101]);
	n[0b011] = glm::dot(d[0b011], g[0b011]);
	n[0b111] = glm::dot(d[0b111], g[0b111]);
}

float PerlinNoise3D::interpolateValues(float a, float b, float w){
	return (b - a) * w + a;
}

int PerlinNoise3D::getPermutation(glm::ivec3 c){ // fetch coherent pseudo-random permutation value
	int permutation[512] = {
		151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 
        103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 247, 120, 234, 75, 0, 
        26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33, 88, 237, 149, 56, 
        87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166, 
        77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 
        46, 245, 40, 244, 102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 
        187, 208, 89, 18, 169, 200, 196, 135, 130, 116, 188, 159, 86, 164, 100, 109, 
        198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123, 5, 202, 38, 147, 118, 126, 
        255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 
        170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 
        172, 9, 129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 
        104, 218, 246, 97, 228, 251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 
        241, 81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31, 181, 199, 106, 
        157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 
        93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180, 
		// ... repeat ...
		151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 
        103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 247, 120, 234, 75, 0, 
        26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33, 88, 237, 149, 56, 
        87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166, 
        77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 
        46, 245, 40, 244, 102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 
        187, 208, 89, 18, 169, 200, 196, 135, 130, 116, 188, 159, 86, 164, 100, 109, 
        198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123, 5, 202, 38, 147, 118, 126, 
        255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 
        170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 
        172, 9, 129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 
        104, 218, 246, 97, 228, 251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 
        241, 81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31, 181, 199, 106, 
        157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 
        93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180};
	return permutation[permutation[permutation[c.x] + c.y] + c.z];
}

// updated functions

PerlinNoise3DOriginal::PerlinNoise3DOriginal(glm::ivec3 pn, glm::vec3 s) : PerlinNoise3D(pn, s){
	
	// generate random coherent gradient vector hash
	srand(time(NULL));
	gradient = (glm::vec3*)malloc(sizeof(glm::vec3) * 256);
	for(int i = 0; i < 256; i++){
		float theta = acos(2.f * ((float)rand() / RAND_MAX) - 1.f);
		float phi = 2.f * ((float)rand() / RAND_MAX) * PI;
		gradient[i] = glm::normalize(glm::vec3(
			cos(phi) * sin(theta), 
			sin(phi) * sin(theta), 
			cos(theta)));
	}
}

PerlinNoise3DOriginal::PerlinNoise3DOriginal::~PerlinNoise3DOriginal(){
	
	// cleanup gradient vector hash
	free(gradient);
}

void PerlinNoise3D::getInterpolation(glm::vec3 t, glm::vec3* s){ // smoothstep weights
	*s = (-2.f * t + 3.f) * t*t;
}

void PerlinNoise3D::getGradientProducts(float *n, glm::ivec3 i0, glm::ivec3 i1, glm::vec3 d0, glm::vec3 d1){ // calculate dot products from direction vectors
	
	// corner coordinates
	glm::ivec3 c[0b1000];
	getCornerCoordinates((glm::ivec3*)c, i0, i1);
	
	// gradient vectors dot products
	glm::vec3 g[0b1000];
	getGradientVectors((glm::vec3*)&g, i0, i1);
	getDotProducts(n, d0, d1, (glm::vec3*)&g);
}

void PerlinNoise3DOriginal::getGradientVectors(glm::vec3 *g, glm::ivec3 i0, glm::ivec3 i1){ // fetch hashed gradient vectors
	
	// corner coordinates
	glm::ivec3 c[0b1000];
	getCornerCoordinates((glm::ivec3*)&c, i0, i1);
	
	// random directions
	g[0b000] = gradient[getPermutation(c[0b000])];
	g[0b100] = gradient[getPermutation(c[0b001])];
	g[0b010] = gradient[getPermutation(c[0b010])];
	g[0b110] = gradient[getPermutation(c[0b011])];
	g[0b001] = gradient[getPermutation(c[0b100])];
	g[0b101] = gradient[getPermutation(c[0b101])];
	g[0b011] = gradient[getPermutation(c[0b110])];
	g[0b111] = gradient[getPermutation(c[0b111])];
}

void PerlinNoise3DImproved::getGradientVectors(glm::vec3 *g, glm::ivec3 i0, glm::ivec3 i1){
	
	// corner coordinates
	glm::ivec3 c[0b1000];
	getCornerCoordinates((glm::ivec3*)&c, i0, i1);

	// twelve directions
	glm::vec3 G[16] = { // lookup index taken from lo-value 4 bits sample
		glm::vec3(1,1,0), glm::vec3(-1,1,0), glm::vec3(1,-1,0), glm::vec3(-1,-1,0), 
		glm::vec3(1,0,1), glm::vec3(-1,0,1), glm::vec3(1,0,-1), glm::vec3(-1,0,-1), 
		glm::vec3(0,1,1), glm::vec3(0,-1,1), glm::vec3(0,1,-1), glm::vec3(0,-1,-1), 
		glm::vec3(1,1,0), glm::vec3(-1,1,0), glm::vec3(0,-1,1), glm::vec3(0,-1,-1)};
	g[0b000] = G[getPermutation(c[0b000]) & 15];
	g[0b100] = G[getPermutation(c[0b001]) & 15];
	g[0b010] = G[getPermutation(c[0b010]) & 15];
	g[0b110] = G[getPermutation(c[0b011]) & 15];
	g[0b001] = G[getPermutation(c[0b100]) & 15];
	g[0b101] = G[getPermutation(c[0b101]) & 15];
	g[0b011] = G[getPermutation(c[0b110]) & 15];
	g[0b111] = G[getPermutation(c[0b111]) & 15];
}

void PerlinNoise3DImproved::getInterpolation(glm::vec3 t, glm::vec3* s){
	*s = ((6.f * t - 15.f) * t + 10.f) * t*t*t;
}

void PerlinNoise3DImplicit::getGradientProducts(float *n, glm::ivec3 i0, glm::ivec3 i1, glm::vec3 d0, glm::vec3 d1){
	
	// corner coordinates
	glm::ivec3 c[0b1000];
	getCornerCoordinates((glm::ivec3*)c, i0, i1);
	
	// implicit values
	n[0b000] = getGradientProductImplicit(getPermutation(c[0b000]), d0.x, d0.y, d0.z);
	n[0b100] = getGradientProductImplicit(getPermutation(c[0b001]), d1.x, d0.y, d0.z);
	n[0b010] = getGradientProductImplicit(getPermutation(c[0b010]), d0.x, d1.y, d0.z);
	n[0b110] = getGradientProductImplicit(getPermutation(c[0b011]), d1.x, d1.y, d0.z);
	n[0b001] = getGradientProductImplicit(getPermutation(c[0b100]), d0.x, d0.y, d1.z);
	n[0b101] = getGradientProductImplicit(getPermutation(c[0b101]), d1.x, d0.y, d1.z);
	n[0b011] = getGradientProductImplicit(getPermutation(c[0b110]), d0.x, d1.y, d1.z);
	n[0b111] = getGradientProductImplicit(getPermutation(c[0b111]), d1.x, d1.y, d1.z);
}

// private functions

float PerlinNoise3DImplicit::getGradientProductImplicit(int p, float x, float y, float z){ // since the chosen gradients are along the dimensional axes, dot products are sums of the position's components
	int h = p & 15; // taken lo-value 4 bits as sample from hash
	float u = h<8 ? x : y; // distribute as x/y/z based on bit-2/3/4 of sample
	float v = h<4 ? y : h==12||h==14 ? x : z;
	return ((h&1) == 0 ? u : -u) + ((h&2) == 0 ? v : -v); // sign as +/- based on bit-1 of sample
}

// Perlin Noise in 2 dimensions

struct PerlinNoise2D{
	
	// external usage
	PerlinNoise2D(glm::ivec2 pn, glm::vec2 s);
	float getValueAt(float xin, float yin);
	
	// internal functionality
	protected:
		glm::ivec2 size;
		glm::ivec2 scale;
		glm::vec2 randomGradient2D(int ix, int iy);
		float dotProductGradient2D(int ix, int iy, float x, float y);
		float interpolate(float a, float b, float w);
};

PerlinNoise2D::PerlinNoise2D(glm::ivec2 pn, glm::vec2 s){
	size = pn;
	scale = s;
}

float PerlinNoise2D::getValueAt(float xin, float yin){
	
	float x = (float)xin / scale.x;
	float y = (float)yin / scale.y;
	
	int x0 = (int)x;
	int x1 = x0 + 1;
	int y0 = (int)y;
	int y1 = y0 + 1;
	
	float sx = x - (float)x0;
	float sy = y - (float)y0;
	
	float n0, n1, ix0, ix1;
	
	n0 = dotProductGradient2D(x0, y0, x, y);
	n1 = dotProductGradient2D(x1, y0, x, y);
	ix0 = interpolate(n0, n1, sx);
	
	n0 = dotProductGradient2D(x0, y1, x, y);
	n1 = dotProductGradient2D(x1, y1, x, y);
	ix1 = interpolate(n0, n1, sx);
	
	return interpolate(ix0, ix1, sy);
}

glm::vec2 PerlinNoise2D::randomGradient2D(int ix, int iy){
	const unsigned w = 8 * sizeof(unsigned);
	const unsigned s = w / 2;
	unsigned a = ix, b = iy;
	a *= 3284157443; b ^= a << s | a >> w-s;
	b *= 1911520717; a ^= b << s | b >> w-s;
	a *= 2048419325;
	float random = a * (3.14159265 / ~(~0u >> 1)); // in 0 <-> 2PI
	return glm::vec2(sin(random), cos(random));
}

float PerlinNoise2D::dotProductGradient2D(int ix, int iy, float x, float y){
	glm::vec2 gradient = randomGradient2D(ix, iy); // get random direction vector
	float dx = x - (float)ix;
	float dy = y - (float)iy;
	return (dx * gradient.x + dy * gradient.y);
}

float PerlinNoise2D::interpolate(float a, float b, float w){
	return (b - a) * w + a;
}

#endif