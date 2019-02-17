//============================================================================
// Name        : polynomial.cpp
// Author      : Juliette Rocco
// Date        : November 20th 2018
//============================================================================
//change ints to unsigned ints
#include <iostream>
#include "Polynomial.h"

#ifndef MARMOSET_TESTING
int main() {

        poly_t p_polly{nullptr,0};
        poly_t q_polly{nullptr,0};
        double a_array[5] = {24,10,-15, 0,1};
        unsigned int degree = 4;
//double x = 1.1;
//double r = 3;
        init_poly(p_polly,a_array, degree);
//double a = 2.4;
//double b = 3.5;
//double n = 2;
//poly_multiply (p_polly,p_polly);
        poly_add(p_polly, p_polly);
        return 0;
}
#endif

void init_poly( poly_t &p, double const init_coeffs[],unsigned int const init_degree )
{
        //if nullptr then delete contents
        if (p.a_coeffs != nullptr)
        {
                delete[] p.a_coeffs;
                p.a_coeffs = nullptr;

        }
        //create a new poly with sufficent space and change its degree
        p.a_coeffs = new double[init_degree + 1];
        p.degree = init_degree;
//assign it values
        for (unsigned int i{0}; i < (init_degree + 1); ++i)
        {
                p.a_coeffs[i] = init_coeffs[i];
        }
}
//destroy the polynomial
void destroy_poly( poly_t &p )
{
        delete[] p.a_coeffs;
        p.a_coeffs = nullptr;
}
//return the degree of the polynomial
unsigned int poly_degree( poly_t const &p )
{
        if (p.a_coeffs == nullptr)
        {
                throw 0;
        }
        return p.degree;
}
//return the coefficent at n
double poly_coeff ( poly_t const &p, unsigned int n )
{
        if (p.a_coeffs == nullptr)
        {
                throw 0;
        }

        return p.a_coeffs[n];
}
double poly_val( poly_t const &p, double x )
{
        //if nullptr then throw 0
        if (p.a_coeffs == nullptr)
        {
                throw 0;
        }

        double value = p.a_coeffs[p.degree];

        int new_degree {p.degree - 1};

        //plug in the x and compute the value
        for (int i{new_degree}; i >= 0; --i)
        {
                value = value * x + p.a_coeffs[i];
        }
        return value;
}

void poly_add( poly_t &p, poly_t const &q )
{
        if (p.a_coeffs == nullptr || q.a_coeffs == nullptr)
        {
                throw 0;
        }
//create temp space to store coes
        double *temp_coe{new double[std::max(p.degree,q.degree) + 1]{}};
//add them together
        for (int i{std::max(p.degree,q.degree)}; i >= 0; --i )
        {
                if(p.degree >= i)
                {
                        temp_coe[i] += p.a_coeffs[i];
                }
                if (q.degree >= i)
                {
                        temp_coe[i] += q.a_coeffs[i];
                }
                std::cout << temp_coe[i] <<std::endl;
        }

        p.degree = std::max(p.degree,q.degree);

        int temp_degree =std::max(p.degree,q.degree);

        while (temp_coe[temp_degree]== 0)
        {
                p.degree -=1;
        }

        if (p.degree < 0)
        {
                p.degree = 0;
        }

//assign new coe values  back to p poly
        double *new_coes{ new double[p.degree + 1 ]{}};
        for (unsigned int i{0}; i <= p.degree; ++i)
        {
                new_coes[i] = temp_coe[i];
        }

        delete [] temp_coe;
        temp_coe = nullptr;

        delete [] p.a_coeffs;
        p.a_coeffs = new_coes;

}
void poly_subtract( poly_t &p, poly_t const &q )
{
        //if empty than throw 0
        if (p.a_coeffs == nullptr || q.a_coeffs == nullptr)
        {
                throw 0;
        }
//create new temp coe
        double *temp_coe{new double[std::max(p.degree,q.degree) + 1]{}};

        //subtract them
        for (int i{std::max(p.degree,q.degree)}; i >= 0; --i )
        {
                if(p.degree >= i)
                {
                        temp_coe[i] += p.a_coeffs[i];
                }
                if (q.degree >= i)
                {
                        temp_coe[i] -= q.a_coeffs[i];
                }
        }
//degree is the highest degree
        p.degree = std::max(p.degree,q.degree);

        int temp_degree =std::max(p.degree,q.degree);

        while (temp_coe[temp_degree]== 0)
        {
                p.degree -=1;
        }

        if (p.degree < 0)
        {
                p.degree = 0;
        }


        double *new_coes{ new double[p.degree + 1 ]{}};
//reassign the space to p coe
        for (unsigned int i{0}; i <= p.degree; ++i)
        {
                new_coes[i] = temp_coe[i];
        }

        delete [] temp_coe;
        temp_coe = nullptr;

        delete [] p.a_coeffs;
        p.a_coeffs = new_coes;

}

void poly_multiply( poly_t &p, poly_t const &q )
{
        //if its empty
        if (p.a_coeffs == nullptr || q.a_coeffs == nullptr)
        {
                throw 0;
        }
        double *temp_coe{new double[p.degree + q.degree + 1]{}};
        //mulitply the first term by the second term
        for (int i=0; i<= p.degree; i++)
        {

                // multiply first term by all the second ones
                for (int j=0; j<=q.degree; j++)
                        temp_coe[i+j] += p.a_coeffs[i]*q.a_coeffs[j];
        }
//updates the degree
        p.degree += q.degree;

        delete [] p.a_coeffs;
        p.a_coeffs = temp_coe;
        temp_coe = nullptr;

}
double poly_divide( poly_t &p, double r )
{
        double value{0};
//if there is no coes throw 0
        if (p.a_coeffs == nullptr)
        {
                throw 0;
        }
        //if the degree of the function is 0 then return 0
        if (p.degree == 0)
        {
                double temp_coe = p.a_coeffs[0];
                p.a_coeffs[0] = 0;
                return temp_coe;
        }
        else
        {

                double *new_coes{ new double[p.degree]{}};

//perform synthetic division
                for (unsigned int i{p.degree}; i > 0; --i)
                {
                        value *= r;
                        value += p.a_coeffs[i];
                        new_coes[i-1] = value;
                }
                double remainder = (new_coes[0]*r + p.a_coeffs[0]);

                p.degree -= 1;
//delete the coe
                delete [] p.a_coeffs;
                p.a_coeffs = new_coes;
                new_coes = nullptr;

                return remainder;
        }
}



void poly_diff( poly_t &p )
{
        //check if empty
        if (p.a_coeffs == nullptr)
        {
                throw 0;
        }

        //check if the degree is 0
        if (p.degree == 0)
        {
                p.a_coeffs[0] = 0;
        }
        else
        {
                double *new_coes{new double[p.degree]{}};
//perform differentiation
                for (unsigned int i{0}; i < p.degree; ++i)
                {
                        new_coes[i] = (i + 1) * p.a_coeffs[i + 1];
                        std::cout <<new_coes[i] <<std::endl;
                }
                p.degree -= 1;
//delete the polys and reassign
                delete [] p.a_coeffs;
                p.a_coeffs = new_coes;
                new_coes = nullptr;

        }

}
double poly_approx_int( poly_t const &p, double a, double b, unsigned int n )
{
        if (p.a_coeffs == nullptr)
        {
                throw 0;
        }
//perform sum math

        double h_value = (b-a)/n;
        double sum{0};
        double x_eval{0};
        sum = poly_val(p,b) + poly_val(p,a);

        for (unsigned int k{1}; k < n; ++k )
        {
                //loop to perform the calculation
                x_eval = a + ( k* h_value);
                sum += 2 * poly_val(p, x_eval);
        }

        sum *= h_value /2;

        std::cout << sum;

        return sum;
}




