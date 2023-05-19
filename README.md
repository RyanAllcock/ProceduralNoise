~~~ Perlin Noise Collection ~~~

~ How To Use:
	- compile with C++ using GLEW header files
	- attach to appropriate program then run (e.g. point data for Marching Cubes/Squares)

~ In This Project:
	- 3 versions of 3D Perlin Noise
	- 2D Perlin Noise
	- Fractal Noise interface for layered noise

~ Perlin Noise:
	- Used widespread in early CG industry
	- Initially used academically for RAM-limited complex image synthesis
	- Defined early instances of computer-generated noise, and is usually just called "noise" by early papers

~ General Algorithm of Perlin Noise: (https://en.wikipedia.org/wiki/Perlin_noise#cite_ref-ken_powerpoint_1-0)
	- define n-dimensional grid of cells, with a gradient vector at each grid-intersection/cell-corner
	- for each point...
	- -> calculate an offset vector to the 2^n surrounding cell corners
	- -> calculate the dot product of this offset vector with each surrounding cell corner's gradient vector
	- -> interpolate the dot product values using a smoothing function to find the final noise value

~ Types of Perlin Noise: https://stackoverflow.com/questions/64781549/perlin-noise-different-implementations
	- Original (smoothstep fade, 256 random direction gradient vectors lookup) https://web.archive.org/web/20180616011332/http://www.scratchapixel.com/lessons/procedural-generation-virtual-worlds%20/perlin-noise-part-2
	- Improved (smootherstep fade, 12 direction gradient vectors lookup OR direct dot product calculation) https://mrl.cs.nyu.edu/~perlin/paper445.pdf
	- Implicit (implicitly calculates the dot products without using gradient vector lookup; gradient vectors are axis-aligned) https://mrl.cs.nyu.edu/~perlin/noise/

~ What to add next:
	- add more comments
	- do perlin.h <- noise.h & fractal.h <- noise.h superclass stucture
		> separate noise struct for noise selections & layering (fractal noise) (is this possible in c++, with polymorphism; only the most general superclass's functionality can be used for a "Noise[]" list?)
	- do other types of noise (e.g. value noise, worley noise, etc.)
	- add test environment (start marching squares/cubes repository and use barebones implementation in a test-env. folder)
	- add precomputing bulk lattices later; for now, just compute the surrounding 2^n corner gradient vectors per point done
	- split noiseND into noise2D & noise3D, then those into Perlin2D and Perlin3D
	- compute vertex normals from the derivatives used in Perlin Noise
	- add options to clamp/scale-with-maximum-value/etc. the Fractal Noise point value output
	- find possibly more versions of Perlin Noise to add (e.g. Simplex Noise https://en.wikipedia.org/wiki/Simplex_noise)

~ More Interesting/Extra Links:
	- (Perlin Noise for dynamic textures, with interesting code to inspect) https://stackoverflow.com/questions/70985360/how-to-create-a-3d-random-gradient-out-of-3-passed-values-in-a-fragment-shader
	- (source finder) https://www.bing.com/search?q=random+3d+direction+vector&qs=n&form=QBRE&sp=-1&ghc=1&lq=0&pq=random+3d+direction+vecto&sc=10-25&sk=&cvid=ADC8A6DD38EA453BA83775BBA3300EEE&ghsh=0&ghacc=0&ghpl=
	- (source finder) https://www.youtube.com/results?search_query=3d+perlin+noise
	- (source finder) https://www.bing.com/search?q=3d+perlin+noise&PC=U316&FORM=CHROMN
	- (Perlin Noise potential parameters) https://substance3d.adobe.com/documentation/sddoc/3d-perlin-noise-168199157.html
	- (noise distortion techniques) https://www.bing.com/search?q=perlin+noise+distortion&PC=U316&FORM=CHROMN
	- (source finder) https://www.bing.com/search?q=versions+of+perlin+noise+implementation&qs=n&form=QBRE&sp=-1&ghc=1&lq=0&pq=versions+of+perlin+noise+implementa&sc=10-35&sk=&cvid=E4D1C90C68284C4F973689A775FAA6DC&ghsh=0&ghacc=0&ghpl=
	- (Perlin Noise creator insight and noise-generated bump maps, an interesting NVidia book to investigate on GPU programming) https://developer.nvidia.com/gpugems/gpugems/part-i-natural-effects/chapter-5-implementing-improved-perlin-noises
	- (Perlin Noise function explanations and "octave"s) https://adrianb.io/2014/08/09/perlinnoise.html
	- (Fractal Brownian Motion introduction, possible Fractal Noise with specific parameters and persistence? octave and frequency elaborations) https://rtouti.github.io/graphics/perlin-noise-algorithm
	- (Perlin Noise limitations) https://www.reddit.com/r/proceduralgeneration/comments/qq7ia5/looking_for_alternatives_to_perlinsimplex_noise/
	- (procedural noise analysis) https://inria.hal.science/hal-00920177/document
	- (Simplex Noise implementation concepts) https://stackoverflow.com/questions/18885440/why-does-simplex-noise-seem-to-have-more-artifacts-than-classic-perlin-noise
	- (original Perlin Noise code and sources) https://web.archive.org/web/20070807065855/http://mrl.nyu.edu/~perlin/doc/oscar.html
	- (academic book incorporating Perlin Noise for texturing and modelling) https://www.bing.com/search?q=EBERT%2C+D.+ET+AL.+1998.+Texturing+and+Modeling%3B+A+Procedural+Approach%2C+Second+Edition.+AP+Professional%2C+Cambridge.&PC=U316&FORM=CHROMN