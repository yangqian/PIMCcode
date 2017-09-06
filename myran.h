#include "pcg64.h"
struct Ran {
    pcg64_random_t rng;
	Ran(uint64_t v=42u,uint64_t w=52u){
    pcg64_srandom_r(&rng, PCG_128BIT_CONSTANT(0ULL, v), PCG_128BIT_CONSTANT(0ULL, w));
	}
	inline uint64_t int64() {
    return pcg64_random_r(&rng);
	}
	inline double doub() { return 5.42101086242752217E-20 * int64(); }
};
