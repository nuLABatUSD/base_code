#include "base_arrays.hh"
#include "gl_vals.hh"
#include "gel_vals.hh"

#include <iostream>
#include <cmath>

using std::cout;
using std::endl;
using std::ostream;

array::array(int N_values)
{   N = N_values;}

double array::get_value(int i)
{   return values[i];}

double array::get_value_from_end(int i)
{   
    if(i < 1){
        cout << "ERROR: in get_value_from_end attempting to access a non-positive element, i = " << i << endl;
        return 0.;
    }
    return values[N-i];}

int array::get_length()
{    return N;}

int array::get_len()
{   return get_length();}

int array::length()
{   return get_length();}

void array::set_value(int i, double v)
{   values[i] = v;}

void array::print_all(){
    for(int i =0; i<N; i++){
        cout << values[i] << endl;
    }
}

void array::print(int N_top, int N_bot) // default values N_top = 3; N_bot = 1
{
    if (N <= N_top + N_bot)
        print_all();
    else if (N_top < 0 || N_bot < 0)
        print_all();
    else
    {
        for (int i = 0; i < N_top; i++)
            cout << values[i] << endl;
        cout << "..." << endl;
        for (int i = 0; i < N_bot; i++)
            cout << values[N - N_bot + i] << endl;
    }
}


dummy_vars::dummy_vars() : array(-1)
{   ; }

dummy_vars::dummy_vars(int num) : array(num){
    values = new double[N]();
    weights = new double[N]();
    max_linspace = 0.;
}

dummy_vars::dummy_vars(dummy_vars* copy_me) : array(copy_me->get_length())
{
    values = new double[N]();
    weights = new double[N]();
    max_linspace = copy_me->get_max_linspace();

    for(int i = 0; i<N; i++)
        {
            values[i] = copy_me->get_value(i);
            weights[i] = copy_me->get_weight(i);
        }
}

dummy_vars::~dummy_vars(){
    if(N != -1){
        delete[] values;
        delete[] weights;
    }
}

double dummy_vars::get_weight(int i)
{ return weights[i]; }

double dummy_vars::get_weight_from_end(int i)
{   return weights[N-i];}

double dummy_vars::get_max_linspace()
{   return max_linspace;  }

double dummy_vars::get_max_value()
{   return values[N-1];}

void dummy_vars::set_weight(int i, double w)
{weights[i] = w;}

void dummy_vars::set_trap_weights(){
   weights[0] = 0.5 * (values[1] - values[0]);
   weights[N-1] = 0.5 * (values[N-1] - values[N-2]);
    for(int i=1; i<N-1; i++){
       weights[i] = 0.5 * (values[i+1] - values[i-1]);
   }
}

int dummy_vars::bin_below(double x){
    if(N < 2)
        return N;
        
    if(values[N-1] <= x)
        return N-1;
        
    int guess = (int) ( (x - values[0]) / ((values[N-1] - values[0]) / (N-1)) );
        
    if (guess >= N-1){
        if (values[N-1] <= x)
            return N-1;
        guess = N-2;
    }
    else if( values[guess] <= x && x < values[guess+1])
        return guess;
        
    int pm = 1;
    if(values[guess] > x)
        pm = -1;
        
    int g;
    for(int i = 1; i < N; i++){
        g = guess + pm * i;
        if (g == 0 || g == N-1)
            return g;
            
        if (values[g] <= x && x < values[g+1])
            return g;
    }
    return g;
}

int dummy_vars::index_below_for_interpolation(double e_val)
{
    if(e_val > values[N-1])
        return N-1;
   
    if(e_val < values[0])
    {
        cout << "WARNING: energy val to interpolate=" << e_val << ", least dummy_vars val=" << values[0] << endl;
        return -999;
    }
    
    if(e_val == values[0])
    {
        return 0;
    }
       
    return bin_below(e_val);
/*    for(int i = 0; i < N; i++){
        if(e_val <= values[i]){
            return i-1;
        }
    }
    
    
    return -99; */
}


double dummy_vars::integrate(dep_vars* fvals){
    double result = 0;
    for (int i = 0; i<N; i++){
       result += fvals->get_value(i) * weights[i]; 
    }
    return result;    
}

double dummy_vars::partial_integrate_end(int start_bin, dep_vars* fvals){
    double result = 0.;
    for(int i = start_bin; i < N; i++)
       result += fvals->get_value(i) * weights[i]; 
    return result;
}

double dummy_vars::partial_integrate_pow_end(int start_bin, dep_vars* fvals, double n_exp){
    double result = 0.;
    for(int i = start_bin; i < N; i++)
       result += fvals->get_value(i) * weights[i] * pow(values[i], n_exp); 
    return result;    
}

double dummy_vars::integrate_pow(dep_vars* fvals, double n_exp)
{   return partial_integrate_pow_end(0, fvals, n_exp);  }

void dummy_vars::print_csv(ostream& os)
{
    for (int i = 0; i < N-1; i++)
        os << values[i] << ", ";
    os << values[N-1] << endl;
    
    for(int i = 0; i < N-1; i++)
        os << weights[i] << ", ";
    os << weights[N-1];
}

gl_dummy_vars::gl_dummy_vars(int num_gl, double start) : dummy_vars(num_gl)
{
    const double* val;
    const double* w;
        
    switch(num_gl){
        case 2:
            val = xvals_2;
            w = wvals_2;
            break;
        case 5:
            val = xvals_5;
            w = wvals_5;
            break;
        case 10:
            val = xvals_10;
            w = wvals_10;
            break;
        case 50:
            val = xvals_50;
            w = wvals_50;
            break;
        default:
            cout << "Error: This Gauss Legendre number is not supported" << endl;
            return;    
    }
    
    for(int i = 0; i < num_gl; i++){
        values[i] = val[i] + start;
        weights[i] = w[i] * exp(val[i]);
    }
}

gl_dummy_vars::gl_dummy_vars(gl_dummy_vars* copy_me) : dummy_vars(copy_me)
{   ;}

gel_dummy_vars::gel_dummy_vars(int num_gel, double start, double end) : dummy_vars(num_gel)
{
    double half_width = (end - start) / 2.;
    double slope_shift = half_width;
    double shift = (end + start) / 2.;
    
    const double* val;
    const double* w;
        
    switch(num_gel){
        case 2:
            val = gel_vals_2;
            w = gel_weights_2;
            break;
        case 5:
            val = gel_vals_5;
            w = gel_weights_5;
            break;
        case 10:
            val = gel_vals_10;
            w = gel_weights_10;
            break;
        case 50:
            val = gel_vals_50;
            w = gel_weights_50;
            break;
        case 100:
            val = gel_vals;
            w = gel_weights;
            break;
        default:
            cout << "Error: This Gauss Legendre number is not supported" << endl;
            return;
    }
    
    for(int i = 0; i < num_gel; i++){
        values[i] = slope_shift * val[i] + shift;
        weights[i] = half_width * w[i];
    }
    
}

gel_dummy_vars::gel_dummy_vars(gel_dummy_vars* copy_me) : dummy_vars(copy_me)
{   ;}



dep_vars::dep_vars(int size) : array(size)
{
    values = new double[N]();
}
    
dep_vars::dep_vars(double* copy_me, int size) : array(size)
{
    values = new double[N]();
    for (int i = 0; i < N; i++)
        values[i] = copy_me[i];
}
    
dep_vars::dep_vars(dep_vars* copy_me) : array(copy_me->get_length())
{
    values = new double[N]();
    for (int i = 0; i < N; i++)
        values[i] = copy_me->get_value(i);
}
    
dep_vars::~dep_vars()
{delete[] values;}

bool dep_vars::isnan(){
    for(int i = 0; i < N; i++)
        if(std::isnan(values[i])){
            cout << "ERROR: dep_vars object is nan at index " << i << endl;
            return true;
        }
    return false;
}

void dep_vars::zeros()
{
    for (int i = 0; i < N; i++)
        values[i] = 0.0;
}

void dep_vars::copy(dep_vars* z)
{
    for (int i = 0; i < N; i++) 
        values[i] = z -> get_value(i);
}

void dep_vars::multiply_by(double scalar)
{
    for (int i = 0; i < N; ++i) 
        values[i] *= scalar;
}

void dep_vars::multiply_by(dep_vars* vector)
{
    if(vector->get_length() != N){
        cout << "Error using dep_vars::multiply_by(dep_vars*). Two dep_vars objects have different lengths." << endl;
        return;
    }
    
    for(int i = 0; i < N; i++)
        values[i] *= vector->get_value(i);
}

void dep_vars::add_to(double c, dep_vars* z)
{
    for (int i = 0; i < N; ++i)
        values[i] += c * z -> get_value(i);
}


void dep_vars::print_csv(ostream& os)
{
    for (int i = 0; i < N-1; i++)
        os << values[i] << ", ";
    os << values[N-1];
}


