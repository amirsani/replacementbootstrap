#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <random>
#include <cmath>
#include <ctime>
#include <unordered_map>
#include <vector>
#include <chrono>
#include "functions.h"
#include "bootstrap_methods.h"

using namespace std;
using namespace chrono;

#define output_to_file 0
#define print_option 0
#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

/********************************************************************
 *
 * 						R-Bootstrap
 *
 *
 ********************************************************************/

inline int get_block (
        const int Kn_1,
        const int A,
        const int idx0,
        const int idx1,
        const int idx2,
        const int idx3,
        const int idx4)
{
    return idx0*Kn_1*A*A*Kn_1 + idx1*A*A*Kn_1 + idx2*A*Kn_1 + idx3*Kn_1 + idx4;
}

inline int get_original (
        const int Kn_1,
        const int n,
        const int idx0,
        const int idx1,
        const int idx2)
{
    return idx0*Kn_1*n + idx1*n + idx2;
}

inline int get_done (
        const int  Kn_1,
        const int  A,
        const int  idx0,
        const int  idx1,
        const int  idx2,
        const int  idx3)
{
    return idx0*Kn_1*A*A + idx1*A*A + idx2*A + idx3;
}

inline int get_nu (
        const int  A,
        const int idx0,
        const int  idx1)
{
    return idx0*A + idx1;
}

inline int get_K (
        const int  n,
        const int  idx0,
        const int  idx1)
{
    return idx0*n + idx1;
}

inline int find_idx (
        const int pattern_length,
        const int t,
        const int x[],
        const int Alphabet_exp[])
{
    int	idx = x[t-pattern_length]+1;
    for (int counter = 1; counter<pattern_length; ++counter)
        idx += (x[t-(pattern_length-counter)]+1) * Alphabet_exp[counter];

    return idx;
}

template <typename K, typename V>
V Get_nusum (
        const unordered_map<K,V> & m,
        const K & key,
        double A2)
{
    typename unordered_map<K,V>::const_iterator it = m.find(key);
    if(unlikely(it != m.end()))
        return it->second;
    else
        return A2;
}

template <typename K, typename V>
V Get_nucount (
        const unordered_map<K,V> & m,
        const K & key)
{
    typename unordered_map<K,V>::const_iterator it = m.find(key);
    if(unlikely(it != m.end()))
        return it->second;
    else
        return 0.5;
}

double compute_syn_prod (
        const int Kn_1,
        const int n,
        const int order_i,
        const int order_j,
        const int from,
        const int to,
        const double R_SYNTHETIC[])
{
    double res = 1.0;
    for (int t = from; t<to; ++t) // Take the product from t_prime_ij to n
        res *= R_SYNTHETIC[get_original(Kn_1,n,order_i,order_j,t)];
    return res;
}

double compute_kirchevsky (
        const int order,
        int const sequence[],
        const int pos,
        const unordered_map<int,double>& nu_COUNT_MAP,
        const unordered_map<int,double>& nu_SUM_MAP,
        const double nu_0,
        const int Alphabet_exp[],
        const int A,
        const double A2)
{
    if (likely(pos > order))
    {
        int idx = find_idx(order, pos, sequence, Alphabet_exp);
        return Get_nucount(nu_COUNT_MAP,get_nu(A,idx,sequence[pos])) / Get_nusum(nu_SUM_MAP, idx, A2); //nu_sum[idx];
    }
    else if (order != 0)
        return nu_0;
    else
    {
        int zero = 0;
        return Get_nucount(nu_COUNT_MAP,get_nu(A,zero,sequence[pos])) / Get_nusum(nu_SUM_MAP, zero, A2); //nu_sum[0];
    }

}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Notice that all the indexes like the insertion point and t_prime_ij are zero-based, while
// in the notes are all 1-based. This explains some potential inconsistencies between
// the code and the notes in the for-loops. Furthermore, note the in the code
// t_prime_ij is defined with a +1 included, unlike the definition used in the notes.
//////////////////////////////////////////////////////////////////////////////////////////////////

void RBootOneBootstrap(
        const int Kn_1,
        const int A,
        const int n,
        const int x[],
        const int replacements,
        const unordered_map<int,double> &nu_COUNT_MAP,
        const unordered_map<int,double> &nu_SUM_MAP,
        const double R_ORIGINAL[],
        const double w[],
        int xboot[],
        const int seed,
        const int Alphabet_exp[],
        const double A2,
        const int change_option)
{
    //Set Random Number Generators
    mt19937 rand_number(seed);
    uniform_int_distribution<> uni_distribution(0, n-1);

    int     end_point = 0,
            t_prime_ij = 0,
            min_n_t_prime_ij = 0,
            max_ij = 0,
            y[n],
            yprime[n],
            original_value = 0,
            insertion_value = 0,
            insertion_positions[3],
            insert_pos = 0,
            next_t_prime_ij,
            min_n_next_t_prime_ij,
            min_n_past_t_prime_ij,
            past_t_prime_ij;

    double  i_sum = 0.0,
            j_sum = 0.0,
            b_sum = 0.0,
            s_prod = 1.0,
            kirchevsky_i = 0.0,
            kirchevsky_j = 0.0,
            pi_prev = 1.0,
            pi_post = 1.0,
            pi_cur = 1.0,
            nu_0 = 1.0/A,
            KIAS[Kn_1],
            *R_BLOCK = new double[Kn_1*Kn_1*A*A*Kn_1],
            R_SYNTHETIC[Kn_1*Kn_1*n],
            **PI_PROD = new double*[Kn_1],
            **PI_NEXT = new double*[Kn_1];

    fill(KIAS, KIAS+Kn_1, 0.0);

    for (int i=0; i<Kn_1; ++i)
    {
        PI_PROD[i] = new double[Kn_1];
        PI_NEXT[i] = new double[Kn_1];
    }

    for (int i=0; i<Kn_1; ++i)
    {
        for (int j=0; j<Kn_1; ++j)
        {
            PI_PROD[i][j] = 0.0;
            PI_NEXT[i][j] = 0.0;
        }
    }

    int done_test[Kn_1*Kn_1*A*A];

    vector<double> probabilities(A);

    copy(x, x+n, xboot);
    copy(R_ORIGINAL, R_ORIGINAL+(Kn_1*Kn_1*n), R_SYNTHETIC);

    ///////////////////////////////////////////////////////
    // Pre-computation of the insertion points
    ///////////////////////////////////////////////////////
    for (int r = 0; r<3; ++r)
        insertion_positions[r] = uni_distribution(rand_number);
    copy(xboot, xboot+n, y);
    copy(xboot, xboot+n, yprime);

    ///////////////////////////////////////////////////////
    // Generation of the new bootstrapped sequence
    ///////////////////////////////////////////////////////
    int r = 0;
    while (r<replacements)
    {
        insert_pos = insertion_positions[1];

        ///////////////////////////////////////////////////////
        // Reset all the structures
        ///////////////////////////////////////////////////////
        for (int a = 0; a<A; ++a)
            probabilities[a] = 0;

        end_point = Kn_1*Kn_1*A*A*Kn_1;
        fill(R_BLOCK, R_BLOCK+end_point, 1.0);

        end_point = Kn_1*Kn_1*A*A;
        fill(done_test, done_test+end_point, 0);

        ///////////////////////////////////////////////////////
        // Computation of the PI values
        ///////////////////////////////////////////////////////
        int zero = 0;
        for (int i = 0; i < Kn_1; ++i)
        {
            for (int j = 0; j < Kn_1; ++j)
            {
                max_ij = max(i,j);
                t_prime_ij = insert_pos + max_ij + 1;
                min_n_t_prime_ij = min(n, t_prime_ij);

                if (likely(r > 0))
                {
                    past_t_prime_ij = insertion_positions[0] + max_ij + 1;
                    min_n_past_t_prime_ij = min(n, past_t_prime_ij);

                    pi_cur = compute_syn_prod(Kn_1, n, i, j, insertion_positions[0], min_n_past_t_prime_ij, R_SYNTHETIC);

                    PI_PROD[i][j] = PI_PROD[i][j]*pi_cur / PI_NEXT[i][j];

                    if (r < replacements-1)
                    {
                        next_t_prime_ij = insertion_positions[2] + max_ij + 1;
                        min_n_next_t_prime_ij = min(n, next_t_prime_ij);
                        PI_NEXT[i][j] = compute_syn_prod(Kn_1, n, i, j, insertion_positions[2], min_n_next_t_prime_ij, R_SYNTHETIC);
                    }
                }
                else
                {
                    pi_prev = compute_syn_prod(Kn_1, n, i, j, zero, insert_pos, R_SYNTHETIC);
                    pi_post = compute_syn_prod(Kn_1, n, i, j, t_prime_ij, n, R_SYNTHETIC);
                    PI_PROD[i][j] = pi_prev*pi_post;

                    if (replacements > 1)
                    {
                        next_t_prime_ij = insertion_positions[2] + max_ij + 1;
                        min_n_next_t_prime_ij = min(n, next_t_prime_ij);
                        PI_NEXT[i][j] = compute_syn_prod(Kn_1, n, i, j, insertion_positions[2], min_n_next_t_prime_ij, R_SYNTHETIC);
                    }
                }
            }
        }

        insertion_positions[0] = insertion_positions[1];
        insertion_positions[1] = insertion_positions[2];
        insertion_positions[2] = uni_distribution(rand_number);

        for (int a = 0; a < A-1; ++a)
        {
            y[insert_pos] = a;
            i_sum = 0.0;

            for (int i = 0; i < Kn_1; ++i)
            {
                // precompute all the K_{i,a,s} terms that could be used later on
                end_point = insert_pos+Kn_1;
                for (int s = insert_pos; s < end_point; ++s)
                    KIAS[s-insert_pos] = compute_kirchevsky(i, y, s, nu_COUNT_MAP, nu_SUM_MAP, nu_0, Alphabet_exp, A, A2);

                j_sum = 0.0;

                for (int j = 0; j < Kn_1; ++j)
                {
                    b_sum = 0.0;

                    max_ij = max(i,j);
                    t_prime_ij = insert_pos + max_ij + 1;
                    min_n_t_prime_ij = min(n, t_prime_ij);

                    for (int b = 0; b < A; ++b)
                    {
                        s_prod = 1.0;
                        if (i != j || a != b)
                        {
                            yprime[insert_pos] = b;

                            bool R_to_do = false;
                            if( !done_test[get_done(Kn_1, A,i,j,a,b)] )
                            {
                                done_test[get_done(Kn_1, A,i,j,a,b)] = 1;
                                done_test[get_done(Kn_1, A,j,i,b,a)] = 1;
                                R_to_do = true;
                            }

                            // Calculate R Measure
                            for (int s = insert_pos; s < min_n_t_prime_ij; ++s)
                            {
                                int shift = s-insert_pos;
                                if (R_to_do)
                                {
                                    kirchevsky_i = KIAS[shift];
                                    kirchevsky_j = compute_kirchevsky(j, yprime, s, nu_COUNT_MAP, nu_SUM_MAP, nu_0, Alphabet_exp, A, A2);

                                    R_BLOCK[get_block(Kn_1,A,i,j,a,b,shift)] = kirchevsky_j/kirchevsky_i;
                                    R_BLOCK[get_block(Kn_1,A,j,i,b,a,shift)] = kirchevsky_i/kirchevsky_j;
                                }

                                s_prod *= R_BLOCK[get_block(Kn_1,A,i,j,a,b,shift)];
                            }
                        }
                        b_sum += s_prod;
//                        cout << "a: " << a << ", i: " << i << ", j: " << j << ", b: " << b << endl;
                    } // END LOOP B
                    j_sum += w[j]*PI_PROD[i][j]*b_sum;
                } // END LOOP J
                i_sum += w[i]/j_sum;
            } // END LOOP I
            probabilities[a] = i_sum;
//            cout << "P(" << a << ") = " << i_sum << endl;
        } // END LOOP A
        probabilities[A-1] = 1.0;
        for (int a=0; a<A-1; ++a)
            probabilities[A-1] -= probabilities[a];

        //original_value = xboot[insert_pos];
	original_value = x[insert_pos];
        discrete_distribution<> insertion_dist (probabilities.begin(), probabilities.end());
        insertion_value = insertion_dist(rand_number);

        if (insertion_value != original_value)
        {
            if (change_option==1)
                r++; // Actual replacement, so increment

            xboot[insert_pos] = insertion_value;
            for (int i = 0; i<Kn_1; ++i)
            {
                for (int j = 0; j<Kn_1; ++j)
                {
                    max_ij = max(i,j);
                    t_prime_ij = insert_pos + max_ij + 1;
                    min_n_t_prime_ij = min(t_prime_ij,n);

                    for (int s = insert_pos; s<min_n_t_prime_ij; ++s)
                    {
                        int shift = s-insert_pos;
                        R_SYNTHETIC[get_original(Kn_1,n,i,j,s)] = R_BLOCK[get_block(Kn_1,A,i,j,insertion_value,insertion_value,shift)];
                    }
                }
            }
        }

        if (change_option==0)
            r++;

        y[insert_pos] = xboot[insert_pos];
        yprime[insert_pos] = xboot[insert_pos];

    } // END REPLACEMENT LOOP

    for (int i=0; i<Kn_1; i++) {
        delete PI_PROD[i];
        delete PI_NEXT[i];
    }
    delete [] PI_PROD;
    delete [] PI_NEXT;
    delete [] R_BLOCK;
}

/********************************************************************
 *
 *                      RBOOT ORIGINAL
 *
 *
 *******************************************************************/
void RBoot (
        const int n,
        const int A,
        const int x[],
        const int replacements,
        const int bootstraps,
        const int Kn,
        const int change_option)
{
    int zero = 0,
        end_point = 0,
        idx = 0,
        dictionary_size = int(2*(pow(2.0,Kn))-1),
        Alphabet_exp[Kn],
        Kn_1 = Kn+1;

    double  nu_0 = 1.0/A,
            A2 = A/2.0,
            r = 0.0,
            w[Kn_1];

    int     **xboots = new int*[bootstraps];

    for(int i=0; i<bootstraps; ++i)
        xboots[i] = new int[n];

    unordered_map<int,double> 	nu_SUM_MAP,
                                nu_COUNT_MAP;

    double *R_ORIGINAL = new double[Kn_1*Kn_1*n],
           *K_VALUES = new double[Kn_1*n];

    end_point = Kn_1*Kn_1*n;
    fill(R_ORIGINAL, R_ORIGINAL+end_point, 1.0);

    end_point = Kn_1*n;
    fill(K_VALUES, K_VALUES+end_point, nu_0);

    for (int i=0; i<Kn; ++i)
        Alphabet_exp[i] = pow(1.0*A,i);

    for (int i=0; i<Kn_1; ++i)
        w[i] = 1/log(i+2)-1/log(i+3);

    /* PRE-COMPUTATION OF COUNTERS */
    for (int i=0; i<A; ++i)
        for (int s=0; s<n; ++s)
            if (x[s] == i)
                nu_COUNT_MAP[get_nu(A,zero,i)]++;

    for (int i=1; i<Kn_1; ++i)
        for (int t=i+1; t<n; ++t)
        {
            idx = find_idx(i, t, x, Alphabet_exp);
            nu_COUNT_MAP[get_nu(A,idx,x[t])]++;
        }

    for (int i=0; i<dictionary_size; ++i)
    {
        double temp_value = 0.0;
        for (int j=0; j<A; ++j)
            temp_value += Get_nucount(nu_COUNT_MAP,get_nu(A,i,j));

        if (temp_value > 0.0)
            nu_SUM_MAP[i] = temp_value+A2;
    }

    typedef unordered_map<int,double>::iterator it_type;
    for(it_type iterator = nu_COUNT_MAP.begin(); iterator != nu_COUNT_MAP.end(); iterator++)
    {
        nu_COUNT_MAP[iterator->first] += 0.5;
    }

    /* PRE-COMPUTATION OF KRICHEVSKY */
    for (int s=0; s<n; ++s)
        K_VALUES[get_K(n,zero,s)] = Get_nucount(nu_COUNT_MAP,get_nu(A,zero,x[s])) / Get_nusum(nu_SUM_MAP, zero, A2); //nu_sum[0];

    for (int i=1; i<Kn_1; ++i)
    {
        for (int s=i; s<n; ++s)
        {
            idx = find_idx(i, s, x, Alphabet_exp);
            K_VALUES[get_K(n,i,s)] = Get_nucount(nu_COUNT_MAP,get_nu(A,idx,x[s])) / Get_nusum(nu_SUM_MAP, idx, A2); //nu_sum[idx];
        }
    }

    /* PRE-COMPUTATION OF R-MEASURES */
    for (int i=0; i<Kn_1; ++i)
        for (int j=i+1; j<Kn_1; ++j)
            for (int s=0; s<n; ++s)
            {
                r = K_VALUES[get_K(n,j,s)]/ K_VALUES[get_K(n,i,s)];
                R_ORIGINAL[get_original(Kn_1, n,i,j,s)] = r;
                R_ORIGINAL[get_original(Kn_1, n,j,i,s)] = 1.0/r;
//                cout << "K_VALUES[get_K(seq," << j << "," << s << ")]/" << "K_VALUES[get_K(seq," << i << "," << s << ")]" << " = " << r << endl;
            }

    for(int seq_idx=0; seq_idx< bootstraps; ++seq_idx)
    {
        RBootOneBootstrap(Kn_1, A, n, x, replacements, nu_COUNT_MAP, nu_SUM_MAP, R_ORIGINAL, w, xboots[seq_idx], seq_idx, Alphabet_exp, A2, change_option);
    }

    stringstream out_name;
    out_name << "boot_n" << n << "_b" << bootstraps << "_r" << replacements << ".tar.gz";
    cout << "Writing bootstrap sequences to " << out_name.str() << endl;
    writeAllSequencesToFile(n, bootstraps, xboots, out_name.str());

    delete [] K_VALUES;
    delete [] R_ORIGINAL;
}
