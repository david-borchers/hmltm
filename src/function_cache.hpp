#ifndef _FUNCACHE_HPP_FDASHUIFDAS
#define _FUNCACHE_HPP_FDASHUIFDAS

union jsdpun {
struct {
   unsigned int   lo;
   unsigned int   hi;
} s  ;
double 	d;
} ;

typedef double(*dblfun)(double);

template<dblfun F>
class FunctionCache
{
	    static const unsigned SZ = 0x1000;
	    struct Entry { double in, out; };
	    Entry table[SZ];
	  public:

      FunctionCache()
      {
        for(unsigned i = 0; i < SZ; ++i)
        {
          table[i].in = 0;
          table[i].out = F(0);
        }
      }

	    double operator()(double x)
	    {
	        jsdpun pun;
	        pun.d = x;
	        unsigned int hash32 = pun.s.lo ^ pun.s.hi;
	        unsigned short hash16 = hash32 ^ (hash32 >> 16);
	        unsigned short index = hash16 & (SZ - 1);
	
	        // FIXME: static init is only valid when F(0) == 0
	        // FIXME: need to put data in JSThreadData
	
	        Entry &e = table[index];
	        if (e.in == x)
	            return e.out;

          e.in = x;
          e.out = F(x);
	        return e.out;
	   }
};

#endif