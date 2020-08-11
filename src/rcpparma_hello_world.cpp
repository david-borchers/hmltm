#include "rcpparma_hello_world.h"

using namespace std;
using namespace Rcpp ;

typedef double floattype;
#include <assert.h>
#include <RcppArmadillo.h>
#include <vector>
#include <string>

// This slightly nasty looking macro provides us a general method
// of switching between the string names, and types, of the different
// measurement functions

#define TYPE_DISPATCH(type_string, funct) \
  if(type_string == "h.EP1") \
  { funct(EP1); } \
  else if(type_string == "h.EP1.0") \
  { funct(EP1_0); } \
  else if(type_string == "h.EP1x.0") \
  { funct(EP1x_0); } \
  else if(type_string == "h.EP2") \
  { funct(EP2); } \
  else if(type_string == "h.EP2.0") \
  { funct(EP2_0); } \
  else if(type_string == "h.EP2x.0") \
  { funct(EP2x_0); } \
  else if(type_string == "h.IP") \
  { funct(IP); } \
  else if(type_string == "h.IP.0") \
  { funct(IP_0); } \
  else \
  { \
    throw std::runtime_error("Unrecognised hazard function name"); \
  }

// This is just a nice macro we can reuse.
#define RVAR(type, x) type x = as<type>(I ## x); 


//#define DOTIME

#ifdef DOTIME
#include <sys/time.h>
#include <sys/resource.h>


inline long double get_cpu_time()
{
    rusage r;
    getrusage(RUSAGE_SELF, &r);
    long double cpu_time = r.ru_utime.tv_sec;
    cpu_time += static_cast<long double>(r.ru_utime.tv_usec) / 1000000.0;
    return cpu_time;
}

double cpu_time = 0;
int calls = 0;
#endif


//#define P(x) cout << x << endl;
#define P(x) ;


#include "function_cache.hpp"

////////////////////////////////
/// Some boring tiny functions
////////////////////////////////

inline double deg2rad(double theta)
{ return theta*PI/180; }

inline double rad2deg(double theta)
{ return theta*180/PI; }

inline double tan_deg2rad(double x)
{ return tan(deg2rad(x)); }

FunctionCache<tan_deg2rad> tan_deg2radCache;


inline double rfrom_xtheta(double x, double theta)
{ return x/sin(deg2rad(theta)); }

inline double yfrom_xtheta(double x, double theta)
{ return x/tan_deg2radCache(theta); }

inline double rfrom_xy(double x, double y)
{ return(sqrt(x*x+y*y)); }

inline double yfrom_xr(double x, double r)
{ return(sqrt(r*r-x*x)); }


inline double logit(double p)
{ return log(p/1-p); }

inline double inv_logit(double x)
{ return 1/(1+exp(-x)); }

inline double myexp(double x)
{ return exp(x); }

inline double mylog(double x)
{ return log(x); }

FunctionCache<inv_logit> inv_logitCache;
FunctionCache<logit> logitCache;
FunctionCache<myexp> expCache;
FunctionCache<mylog> logCache;


////////////////////////////////////////////////
// The different estimation(?) functions
////////////////////////////////////////////////

struct EP1
{
  static double h(double x, double y, const vector<double>& b)
  {
    // vector<double> par = invtfm(b);
    double par1 = inv_logitCache(b[1-1]);
    double par2 = expCache(b[2-1]);
    double par3 = expCache(b[3-1]);

    return (par1*exp(-( pow(abs(x),par2) + pow(y,par2) )/( pow(par3,par2) )));
  }

  static vector<double> tfm(const vector<double>& par)
  {
    vector<double> ret;
    ret.push_back(logitCache(par[1-1]));
    for(int i = 1; i < (int)par.size(); ++i)
      ret.push_back(log(par[i]));
//    vector<double> ret;
//    ret.push_back(logitCache(par[1-1]));
//    ret.push_back(logCache(par[2-1]));
//    ret.push_back(logCache(par[3-1]));
    return ret;
  }

  static vector<double> invtfm(const vector<double>& b)
  {
    vector<double> ret;
    ret.push_back(inv_logitCache(b[1-1]));
    for(int i = 1; i < (int)b.size(); ++i)
      ret.push_back(expCache(b[i]));
//    ret.push_back(inv_logitCache(b[1-1]));
//    ret.push_back(expCache(b[2-1]));
//    ret.push_back(expCache(b[3-1]));
    return ret;
  }
};



struct EP1_0
{
 static double h(double x, double y, const vector<double>& b)
  {
    double par1 = expCache(b[0]);
    double par2 = expCache(b[1]);
    double par3 = expCache(b[2]);
    return(par1*exp(-( pow(abs(x),par2) + pow(y,par2) )/( pow(par3,par2) )));
  }

  static vector<double> tfm(const vector<double>& par)
  {
    vector<double> ret(par.size());
    for(int i = 0; i < (int)par.size(); ++i)
      ret[i] = logCache(par[i]);
    return ret;
  }

  static vector<double> invtfm(const vector<double>& par)
  {
    vector<double> ret(par.size());
    for(int i = 0; i < (int)par.size(); ++i)
      ret[i] = expCache(par[i]);
    return ret;
  }
};

struct EP1x_0
{
    static double h(double x, double y, const vector<double>& b)
    {
        //vector<double> par = invtfm(b);
        // note: inlined invtfm straight in
        return (exp(-(  pow((abs(x)/expCache(b[3-1])),expCache(b[1-1]) ) + pow((abs(y)/expCache(b[2-1])),expCache(b[1-1])) )));
    }
    
    static vector<double> tfm(const vector<double>& par)
    {
        vector<double> ret;
        for(int i = 0; i < (int)par.size(); ++i)
            ret.push_back(log(par[i]));
        return ret;
    }
    
    static vector<double> invtfm(const vector<double>& par)
    {
        vector<double> ret;
        for(int i = 0; i < (int)par.size(); ++i)
            ret.push_back(expCache(par[i]));
        return ret;
    }
};

struct EP2
{
  static double h(double x, double y, const vector<double>& b)
  {
    double par1 = inv_logitCache(b[1-1]);
    double par2 = expCache(b[2-1]);
    double par3 = expCache(b[3-1]);
    double par4 = expCache(b[4-1]);

    return (par1*exp(-( pow((abs(x)/par4),par2) + pow((abs(y)/par4),par3) )));
  }

  static vector<double> tfm(const vector<double>& par)
  {
    vector<double> ret(par.size());
    ret.push_back(logitCache(par[1-1]));
    for(int i = 1; i < (int)par.size(); ++i)
      ret.push_back(log(par[i]));
//    vector<double> ret(4);
//    ret[0] = logitCache(par[0]);
//    ret[1] = logCache(par[1]);
//    ret[2] = logCache(par[2]);
//    ret[3] = logCache(par[3]);
    return ret;
  }

  static vector<double> invtfm(const vector<double>& par)
  {
    vector<double> ret(par.size());
    ret.push_back(inv_logitCache(par[1-1]));
    for(int i = 1; i < (int)par.size(); ++i)
      ret.push_back(expCache(par[i]));
//    vector<double> ret(4);
//    ret[0] = inv_logitCache(par[0]);
//    ret[1] = expCache(par[1]);
//    ret[2] = expCache(par[2]);
//    ret[3] = expCache(par[3]);
    return ret;
  }
};

struct EP2x_0
{
  static double h(double x, double y, const vector<double>& b)
  {
    //vector<double> par = invtfm(b);
    // note: inlined invtfm straight in
    return (exp(-(  pow((abs(x)/expCache(b[4-1])),expCache(b[1-1]) ) + pow((abs(y)/expCache(b[3-1])),expCache(b[2-1])) )));
  }

  static vector<double> tfm(const vector<double>& par)
  {
    vector<double> ret;
    for(int i = 0; i < (int)par.size(); ++i)
      ret.push_back(log(par[i]));
    return ret;
  }

  static vector<double> invtfm(const vector<double>& par)
  {
    vector<double> ret;
    for(int i = 0; i < (int)par.size(); ++i)
      ret.push_back(expCache(par[i]));
    return ret;
  }
};

struct EP2_0
{
  static double h(double x, double y, const vector<double>& b)
  {
    //vector<double> par = invtfm(b);
    // note: inlined invtfm straight in
    return (exp(-(  pow((abs(x)/expCache(b[3-1])),expCache(b[1-1]) ) + pow((abs(y)/expCache(b[3-1])),expCache(b[2-1])) )));
  }

  static vector<double> tfm(const vector<double>& par)
  {
    vector<double> ret;
    for(int i = 0; i < (int)par.size(); ++i)
      ret.push_back(log(par[i]));
    return ret;
  }

  static vector<double> invtfm(const vector<double>& par)
  {
    vector<double> ret;
    for(int i = 0; i < (int)par.size(); ++i)
      ret.push_back(expCache(par[i]));
    return ret;
  }
};

struct IP
{
  static double h(double x, double y, const vector<double>& b)
  {
    double par1 = inv_logitCache(b[1-1]);
    double par2 = expCache(b[2-1]);
    double par3 = expCache(b[3-1]);

    return (par1*exp(par2*log(par3)-(par2/2)*log(par3*par3+x*x+y*y)));
  }

  static vector<double> tfm(const vector<double>& par)
  {
    vector<double> ret;
    ret.push_back(logCache(par[0]));
    for(int i = 1; i < (int)par.size(); ++i)
      ret.push_back(log(par[i]));
//    vector<double> ret(3);
//    ret[0] = logitCache(par[0]);
//    ret[1] = logCache(par[1]);
//   ret[2] = logCache(par[2]);
    return ret;
  }

  static vector<double> invtfm(const vector<double>& par)
  {
    vector<double> ret;
    ret.push_back(inv_logitCache(par[0]));
    for(int i = 1; i < (int)par.size(); ++i)
      ret.push_back(expCache(par[i]));
//    vector<double> ret(3);
//    ret[0] = inv_logitCache(par[0]);
//    ret[1] = expCache(par[1]);
//    ret[2] = expCache(par[2]);
    return ret;
  }
};

struct IP_0
{
  static double h(double x, double y, const vector<double>& b)
  {
    double par1 = expCache(b[1-1]);
    double par2 = expCache(b[2-1]);

    return(exp(par1*logCache(par2)-(par1/2)*log(par2*par2+x*x+y*y)));
  }

  static vector<double> tfm(const vector<double>& par)
  {
    vector<double> ret;
    for(int i = 0; i < (int)par.size(); ++i)
      ret.push_back(log(par[i]));
//    vector<double> ret(2);
//    ret[0] = logCache(par[0]);
//    ret[1] = logCache(par[1]);
    return ret;
  }

  static vector<double> invtfm(const vector<double>& par)
  {
    vector<double> ret;
    for(int i = 0; i < (int)par.size(); ++i)
      ret.push_back(expCache(par[i]));
//    vector<double> ret(2);
//    ret[0] = expCache(par[0]);
//    ret[1] = expCache(par[1]);
    return ret;
  }
};

#define GET_TFM(T) v = T::tfm(v)
#define GET_INVTFM(T) v = T::invtfm(v)

SEXP hmltm_get_tfm(SEXP Iv, SEXP Ifun)
{
BEGIN_RCPP
  RVAR(std::string, fun);
  RVAR(vector<double>, v);
  TYPE_DISPATCH(fun, GET_TFM);
  return wrap(v);
END_RCPP
}

SEXP hmltm_get_invtfm(SEXP Iv, SEXP Ifun)
{
BEGIN_RCPP
  RVAR(std::string, fun);
  RVAR(vector<double>, v);
  TYPE_DISPATCH(fun, GET_INVTFM);
  return wrap(v);
END_RCPP
}



vector<double> gety_obs(double x, bool null_yobs, double yobs, double theta_f, double theta_b, double ymax, double dy)
{
  double maxy=-999, miny=-999;
  if(theta_f==0) { maxy=ymax; } else { maxy=yfrom_xtheta(x,theta_f); }  //# max y is as enter dectectable pie slice
  
  if(theta_b==180) 
  { miny=-ymax; }
  else if(theta_b==0)
  { miny = 0; }
  else 
  { miny=yfrom_xtheta(x,theta_b); } //  # min y is as leave dectectable pie slice
  
  if(!null_yobs) { miny = max(miny, yobs); }
              // # ensure that ys dont fall short of maxy:
  if(maxy<=miny) {
    vector<double> yi;
    yi.push_back(miny);
    return yi; 
  }
  double yrange=maxy-miny;
  double ndiv=floor(yrange/dy);
              // This loop was risky, now re-written
              //  if((yrange-ndiv*dy)>0) maxy=miny+(ndiv+1)*dy;
  vector<double> yi;
  yi.reserve(ndiv+1);
  for(int i = 0; i <= ndiv + 1; ++i)
  {
    yi.push_back(miny + i * dy);
    P(yi.back() << " ");
  }
  P("loop done");
  return yi;
}

template<typename h>
void do_pxy(vector<double>& x, vector<double>& y, vector<double>& b, vector<double>& pcu,
            arma::mat& Pi, vector<double>& delta, double ymax, double dy, double theta_f,
            double theta_b, bool ally, bool cdf, vector<double>& p)
{
  int m = pcu.size();

  // stuff for manipulating b:
  int n = x.size();
  int nb = b.size()/n;

  vector<double> bb (nb);
  vector<double> yi;
  vector<double> one(m, 1);

  // Only declare these once for efficency
  arma::vec prodB(m);
  arma::vec lambda(m);
  arma::mat pi_mult(m,m);

  for(unsigned i = 0; i < n; ++i)
  {
    prodB = one;

    for(unsigned j = 0; j <nb; j++) bb[j] = b[i*nb+j]; // select appropriate row of b

    if(ally) 
      yi=gety_obs(x[i],true,0,theta_f,theta_b,ymax,dy); // # get ys from min y in view to ymax
    else
      yi=gety_obs(x[i],false,y[i],theta_f,theta_b,ymax,dy); // # get y's from y[i] to ymax
    
    int nT = yi.size();

    if(nT > 1)
    {

      for(int u = nT - 1; u >= 1; --u) // # do calculation for nT-1 non-detection events
      {
        double hval = h::h(x[i], yi[u], bb);
        
        P("," << hval);
        for(int i = 0; i < m; ++i)
          lambda.at(i) = 1 - pcu[i] * hval;

        for(int i = 0; i < m; ++i)
        {
          for(int j = 0; j < m; ++j)
          {
            pi_mult.at(i,j) = Pi.at(i,j) * lambda(j);
          }
        }
        prodB = (pi_mult * prodB);
      }
    }

    if(cdf || ally) // # if not seen, calculation for last non-detection event
    {
      for(int j = 0; j < m; ++j)
      { lambda.at(j,j) = 1 - pcu[j] * h::h(x[i], yi[0], bb); }
    }
    else // # if seen, calculation for last detection event
    {
      for(int j = 0; j < m; ++j)
      { lambda.at(j,j) = pcu[j] * h::h(x[i], yi[0], bb); }
    }

    for(int t = 0; t < m; ++t)
    {
      for(int j = 0; j < m; ++j)
      {
        pi_mult.at(t,j) = Pi.at(t,j) * lambda(j);
      }
    }
    prodB = (pi_mult * prodB); // # do last update

    p[i] = as_scalar(arma::rowvec(&*delta.begin(), delta.size(), false)*prodB); // # calculate final prob

    if(abs(p[i]-1) < 1e-10) p[i]=1; //# very clunky way of dealing with numerical underflow
  }

  if(cdf || ally)
  {
    for(int i = 0; i < p.size(); ++i)
      p[i] = 1 - p[i]; // # since in this case p calculated so far is prob(NOT seen)

  }
}

#ifdef DOTIME
double cpu_start = get_cpu_time();
#endif

using namespace std;
using namespace Rcpp;



SEXP bias_p_xy1(SEXP Ix, SEXP Iy, SEXP Ihfun, SEXP Ib, SEXP Ipcu, SEXP IPi, SEXP Idelta, SEXP Iymax, 
				SEXP Idy, SEXP Itheta_f, SEXP Itheta_b, SEXP Ially, SEXP Icdf)
{
BEGIN_RCPP
  RVAR(vector<double>, x);

  vector<double> y;
  if(!RObject(Iy).isNULL())
  { y = as<vector<double> >(Iy); }

  RVAR(string, hfun);
  RVAR(vector<double>, b);
  RVAR(vector<double>, pcu);

  NumericMatrix NMPi = IPi;

  arma::mat Pi(NMPi.begin(), NMPi.nrow(), NMPi.ncol(), true);

  RVAR(vector<double>, delta);
  RVAR(double, ymax);
  RVAR(double, dy);
  RVAR(double, theta_f);
  RVAR(double, theta_b);
  RVAR(bool, ally);
  RVAR(bool, cdf);

  vector<double> p(x.size());
  #define DOFUN(Type) do_pxy< Type >(x, y, b, pcu, Pi, delta, ymax, dy, theta_f, theta_b, ally, cdf, p);



  TYPE_DISPATCH(hfun, DOFUN)


#ifdef DOTIME
  double cpu_end = get_cpu_time();
  double cpu_cache = cpu_time;
  cpu_time += cpu_end - cpu_start;
  calls++;
 
  if((int)cpu_cache < (int)cpu_time)
  cout << "CPU: " << cpu_time << ":" << cpu_time/calls << endl;
#endif

  return wrap(p);
END_RCPP
}

SEXP bias_gety_obs(SEXP Ix, SEXP Inull_yobs, SEXP Iyobs, SEXP Itheta_f,
                   SEXP Itheta_b, SEXP Iymax, SEXP Idy)
{
BEGIN_RCPP
  RVAR(double, x);
  RVAR(bool, null_yobs);
  RVAR(double, yobs);
  RVAR(double, theta_f);
  RVAR(double, theta_b);
  RVAR(double, ymax);
  RVAR(double, dy);

  vector<double> p = gety_obs(x,null_yobs,yobs,theta_f,theta_b,ymax,dy);
  return wrap(p);
END_RCPP
}
