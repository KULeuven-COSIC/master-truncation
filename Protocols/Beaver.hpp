/*
 * Beaver.cpp
 *
 */

#ifndef PROTOCOLS_BEAVER_HPP_
#define PROTOCOLS_BEAVER_HPP_

#include "Beaver.h"

#include "Replicated.hpp"

#include <array>

template<class T>
inline Beaver<T>::Beaver(Player& P) : prep(0), MC(0), P(P)
{
    assert(P.num_players() == 3);
	if (not P.is_encrypted())
		insecure("unencrypted communication", false);

  #ifdef OUR_TRUNC
  if (P.my_num() == 0) {
    shared_prngs[1].SeedPairwise(P, 2); // 1
    shared_prngs[0].SeedPairwise(P, 1); // 2
  }else if(P.my_num() == 1) {
    shared_prngs[1].SeedPairwise(P, 0); // 2
    shared_prngs[0].SeedPairwise(P, 2); // 3
  }else {
    shared_prngs[0].SeedPairwise(P, 0); // 1
    shared_prngs[1].SeedPairwise(P, 1); // 3
  }
  #else
	shared_prngs[0].ReSeed();
	octetStream os;
	os.append(shared_prngs[0].get_seed(), SEED_SIZE);
	P.send_relative(1, os);
	P.receive_relative(-1, os);
	shared_prngs[1].SetSeed(os.get_data());
  #endif
}

#ifdef OUR_TRUNC
template<class T>
inline Beaver<T>::Beaver(Player& P, array<PairwisePRNG, 2>& prngs) :
#else
template<class T>
inline Beaver<T>::Beaver(Player& P, array<PRNG, 2>& prngs) :
#endif
        P(P)
{
    for (int i = 0; i < 2; i++)
        shared_prngs[i].SetSeed(prngs[i]);
}

template<class T>
typename T::Protocol Beaver<T>::branch()
{
    typename T::Protocol res(P);
    res.prep = prep;
    res.MC = MC;
    res.init_mul();
    return res;
}

template<class T>
void Beaver<T>::init(Preprocessing<T>& prep, typename T::MAC_Check& MC)
{
    this->prep = &prep;
    this->MC = &MC;
}

template<class T>
void Beaver<T>::init_mul()
{
    assert(this->prep);
    assert(this->MC);
    shares.clear();
    opened.clear();
    triples.clear();
    lengths.clear();
}

template<class T>
void Beaver<T>::prepare_mul(const T& x, const T& y, int n)
{
    (void) n;
    triples.push_back({{}});
    auto& triple = triples.back();
    triple = prep->get_triple(n);
    shares.push_back(x - triple[0]);
    shares.push_back(y - triple[1]);
    lengths.push_back(n);
}

template<class T>
void Beaver<T>::exchange()
{
    assert(shares.size() == 2 * lengths.size());
    MC->init_open(P, shares.size());
    for (size_t i = 0; i < shares.size(); i++)
        MC->prepare_open(shares[i], lengths[i / 2]);
    MC->exchange(P);
    for (size_t i = 0; i < shares.size(); i++)
        opened.push_back(MC->finalize_raw());
    it = opened.begin();
    triple = triples.begin();
}

template<class T>
void Beaver<T>::start_exchange()
{
    MC->POpen_Begin(opened, shares, P);
}

template<class T>
void Beaver<T>::stop_exchange()
{
    MC->POpen_End(opened, shares, P);
    it = opened.begin();
    triple = triples.begin();
}

template<class T>
T Beaver<T>::finalize_mul(int n)
{
    (void) n;
    typename T::open_type masked[2];
    T& tmp = (*triple)[2];
    for (int k = 0; k < 2; k++)
    {
        masked[k] = *it++;
    }
    tmp += (masked[0] * (*triple)[1]);
    tmp += ((*triple)[0] * masked[1]);
    tmp += T::constant(masked[0] * masked[1], P.my_num(), MC->get_alphai());
    triple++;
    return tmp;
}

template<class T>
void Beaver<T>::check()
{
    assert(MC);
    MC->Check(P);
}

template<class T>
void Beaver<T>::randoms(T& res, int n_bits)
{
    for (int i = 0; i < 2; i++)
        res[i].randomize_part(shared_prngs[i], n_bits);
}

template<class T>
template<class U>
void Beaver<T>::trunc_pr(const vector<int>& regs, int size, U& proc)
{
    // std::cout << "truncpr" << std::endl;
    assert(regs.size() % 4 == 0);
    assert(proc.P.num_players() == 3);
    assert(proc.Proc != 0);
    typedef typename T::clear value_type;
    ArgList<TruncPrTupleWithGap<value_type>> infos(regs);
    auto& S = proc.get_S();

    // use https://eprint.iacr.org/2019/131
    bool have_small_gap = false;
    // use https://eprint.iacr.org/2018/403
    bool have_big_gap = false;

    for (auto info : infos)
        if (info.small_gap())
            have_small_gap = true;
        else
            have_big_gap = true;

    #ifdef OUR_TRUNC
    if (have_big_gap){
        using Z2 = typename T::T;
        int comp_player = 1;
        int check_player1 = 0;
        int check_player2 = 2;
        bool compute = P.my_num() == comp_player;
        bool checker1 = P.my_num() == check_player1;
        bool checker2 = P.my_num() == check_player2;

        std::vector<Z2> Gammas1, Gammas2, Gammas3;

        octetStream cs;
        octetStream rs;
    
        // ROUND 1
        if (compute)
        {
            for (auto info : infos)
                for (int i = 0; i < size; i++)
                {
                    auto& x = S[info.source_base + i];
                    auto& y = S[info.dest_base + i];
                    // y[0] = 0;
                    y[0] = this->shared_prngs[0].template get<value_type>();
                    y[1] = (x.sum() >> info.m) - y[0];
                    y[1].pack(cs);
                }
            P.send_to(check_player1, cs);
        }

        else if (checker1)
        {
            P.receive_player(comp_player, cs); 
            for (auto info : infos)
                for (int i = 0; i < size; i++)
                {
                    auto& x = S[info.source_base + i];
                    auto& y = S[info.dest_base + i]; 
                    y[1] = -value_type(-value_type(x[1]) >> info.m);;
                    y[0] = cs.get<value_type>();

                    // CHECK PREP
                    Z2 s2 = -value_type(-value_type(x[0]) >> info.m);
                    auto r_pp = this->shared_prngs[1].template get<value_type>();
                    // auto r_pp = 0;
                    Z2 s3 = r_pp;
                    Z2 gamma2 = y[0] - s2;
                    Gammas2.push_back(gamma2);
                    Z2 gamma3 = y[1] - s3;
                    Gammas3.push_back(gamma3);
                    gamma2.pack(rs);
                }

            P.send_to(check_player2, rs);
            P.receive_player(check_player2, rs);

            // CHECK
            for (int i = 0; i < size; i++)
            {
                Z2 gamma1 = rs.get<Z2>();
                Z2 check = gamma1 + Gammas2[i] + Gammas3[i];
                if (check != Z2(0) and check != Z2(1) and check != Z2(-1))
                {
                    std::cout << "Abort needed, gamma value: " << check << std::endl;
                    throw std::runtime_error(std::string("Failed gamma check"));    //abort
                }
            }
        }

        else if (checker2) 
        {   
            for (auto info : infos)
                for (int i = 0; i < size; i++)
                {
                    auto& x = S[info.source_base + i];
                    auto& y = S[info.dest_base + i];
                    auto r_prime = this->shared_prngs[1].template get<value_type>();
                    // auto r_prime = 0;
                    y[0] = -value_type(-value_type(x[0]) >> info.m);
                    y[1] = r_prime;

                    //CHECK PREP
                    Z2 r_pp = this->shared_prngs[0].template get<value_type>();
                    // Z2 r_pp = 0;
                    Z2 s1 = (x.sum() >> info.m) - r_pp;
                    Z2 s3 = r_pp;
                    Z2 gamma1 = y[1] - s1;
                    Gammas1.push_back(gamma1);
                    Z2 gamma3 = y[0] - s3; 
                    Gammas3.push_back(gamma3);
                    gamma1.pack(rs);
                }

            P.send_to(check_player1, rs);
            P.receive_player(check_player1, rs);

            // CHECK
            for (int i = 0; i < size; i++)
            {
                Z2 gamma2 = rs.get<Z2>(); 
                Z2 check = Gammas1[i] + gamma2 + Gammas3[i];
                if (check != Z2(0) and check != Z2(1) and check != Z2(-1))
                {
                    std::cout << "Abort needed, gamma value: " << check << std::endl;
                    throw std::runtime_error(std::string("Failed gamma check"));    //abort
                }
            }
        }
    }
    else if (have_small_gap) {
        int gen_player = 2;
        int comp_player = 1;
        bool generate = P.my_num() == gen_player;
        bool compute = P.my_num() == comp_player;
        ArgList<TruncPrTupleWithGap<value_type>> infos(regs);
        auto& S = proc.get_S();

        octetStream cs;
        ReplicatedInput<T> input(P);

        if (generate)
        {
            SeededPRNG G;
            for (auto info : infos)
                for (int i = 0; i < size; i++)
                {
                    auto r = G.get<value_type>();
                    input.add_mine(info.upper(r));
                    if (info.small_gap())
                        input.add_mine(info.msb(r));
                    (r + S[info.source_base + i][0]).pack(cs);
                }
            P.send_to(comp_player, cs);
        }
        else
            input.add_other(gen_player);

        if (compute)
        {
            P.receive_player(gen_player, cs);
            for (auto info : infos)
                for (int i = 0; i < size; i++)
                {
                    auto c = cs.get<value_type>() + S[info.source_base + i].sum();
                    input.add_mine(info.upper(c));
                    if (info.small_gap())
                        input.add_mine(info.msb(c));
                }
        }

        input.add_other(comp_player);
        input.exchange();
        init_mul();

        for (auto info : infos)
            for (int i = 0; i < size; i++)
            {
                this->trunc_pr_counter++;
                auto c_prime = input.finalize(comp_player);
                auto r_prime = input.finalize(gen_player);
                S[info.dest_base + i] = c_prime - r_prime;

                if (info.small_gap())
                {
                    auto c_dprime = input.finalize(comp_player);
                    auto r_msb = input.finalize(gen_player);
                    S[info.dest_base + i] += ((r_msb + c_dprime)
                            << (info.k - info.m));
                    prepare_mul(r_msb, c_dprime);
                }
            }

        exchange();

        for (auto info : infos)
            for (int i = 0; i < size; i++)
                if (info.small_gap())
                    S[info.dest_base + i] -= finalize_mul()
                            << (info.k - info.m + 1);
        
    }
    

    #elif defined(ABY3_MAL_TRUNC)
    if (have_big_gap){
        using Z2 = typename T::T;
        
        octetStream cs;
        std::vector<Z2> masked_s_all;
        value_type r_prime = 0;
        value_type r = 0;

        for (auto info : infos)
            for (int i = 0; i < size; i++)
            {
                auto& x = S[info.source_base + i];
                value_type masked_s = x[0] - r_prime;
                masked_s_all.push_back(masked_s);
                masked_s.pack(cs);
            }

        P.send_relative(-1, cs);
        P.send_relative(1, cs);
        octetStream os;
        P.receive_relative(1, cs);
        P.receive_relative(-1, os);

        for (auto info : infos)
            for (int i = 0; i < size; i++)
            {
                Z2 masked_x = (masked_s_all[i] + cs.get<value_type>() + os.get<value_type>()) >> info.m;
                auto& y = S[info.dest_base + i];
                
                if (P.my_num() == 0)
                {
                    y[0] = r + masked_x;
                    y[1] = r;  
                }
                if (P.my_num() == 1)
                {
                    y[0] = r;
                    y[1] = r + masked_x;  
                }
                if (P.my_num() == 2)
                {
                    y[0] = r;
                    y[1] = r;  
                }    
            }
    }
    else if (have_small_gap) {
        int gen_player = 2;
        int comp_player = 1;
        bool generate = P.my_num() == gen_player;
        bool compute = P.my_num() == comp_player;
        ArgList<TruncPrTupleWithGap<value_type>> infos(regs);
        auto& S = proc.get_S();

        octetStream cs;
        ReplicatedInput<T> input(P);

        if (generate)
        {
            SeededPRNG G;
            for (auto info : infos)
                for (int i = 0; i < size; i++)
                {
                    auto r = G.get<value_type>();
                    input.add_mine(info.upper(r));
                    if (info.small_gap())
                        input.add_mine(info.msb(r));
                    (r + S[info.source_base + i][0]).pack(cs);
                }
            P.send_to(comp_player, cs);
        }
        else
            input.add_other(gen_player);

        if (compute)
        {
            P.receive_player(gen_player, cs);
            for (auto info : infos)
                for (int i = 0; i < size; i++)
                {
                    auto c = cs.get<value_type>() + S[info.source_base + i].sum();
                    input.add_mine(info.upper(c));
                    if (info.small_gap())
                        input.add_mine(info.msb(c));
                }
        }

        input.add_other(comp_player);
        input.exchange();
        init_mul();

        for (auto info : infos)
            for (int i = 0; i < size; i++)
            {
                this->trunc_pr_counter++;
                auto c_prime = input.finalize(comp_player);
                auto r_prime = input.finalize(gen_player);
                S[info.dest_base + i] = c_prime - r_prime;

                if (info.small_gap())
                {
                    auto c_dprime = input.finalize(comp_player);
                    auto r_msb = input.finalize(gen_player);
                    S[info.dest_base + i] += ((r_msb + c_dprime)
                            << (info.k - info.m));
                    prepare_mul(r_msb, c_dprime);
                }
            }

        exchange();

        for (auto info : infos)
            for (int i = 0; i < size; i++)
                if (info.small_gap())
                    S[info.dest_base + i] -= finalize_mul()
                            << (info.k - info.m + 1);
    }
    
    #else 
    int gen_player = 2;
    int comp_player = 1;
    bool generate = P.my_num() == gen_player;
    bool compute = P.my_num() == comp_player;

    octetStream cs;
    ReplicatedInput<T> input(0, *this);

    if (generate)
    {
        SeededPRNG G;
        for (auto info : infos)
            for (int i = 0; i < size; i++)
            {
                auto& x = S[info.source_base + i];
                if (info.small_gap())
                {
                    auto r = G.get<value_type>();
                    input.add_mine(info.upper(r));
                    input.add_mine(info.msb(r));
                    (r + x[0]).pack(cs);
                }
                else
                {
                    
                    auto& y = S[info.dest_base + i];
                    auto r = this->shared_prngs[0].template get<value_type>();
                    y[1] = -value_type(-value_type(x.sum()) >> info.m) - r;
                    y[1].pack(cs);
                    y[0] = r;
                }
            }

        P.send_to(comp_player, cs);
    }
    else if (have_small_gap)
        input.add_other(gen_player);

    if (compute)
    {
        P.receive_player(gen_player, cs);
        for (auto info : infos)
            for (int i = 0; i < size; i++)
            {
                auto& x = S[info.source_base + i];
                if (info.small_gap())
                {
                    auto c = cs.get<value_type>() + x.sum();
                    input.add_mine(info.upper(c));
                    input.add_mine(info.msb(c));
                }
                else
                {
                    auto& y = S[info.dest_base + i];
                    y[0] = cs.get<value_type>();
                    y[1] = x[1] >> info.m;
                }
            }
    }

    if (have_big_gap and not (compute or generate))
    {
        for (auto info : infos)
            if (info.big_gap())
                for (int i = 0; i < size; i++)
                {
                    auto& x = S[info.source_base + i];
                    auto& y = S[info.dest_base + i];
                    y[0] = x[0] >> info.m;
                    y[1] = this->shared_prngs[1].template get<value_type>();
                }
    }

    if (have_small_gap)
    {
        input.add_other(comp_player);
        input.exchange();
        init_mul();

        for (auto info : infos)
            for (int i = 0; i < size; i++)
            {
                if (info.small_gap())
                {
                    this->trunc_pr_counter++;
                    auto c_prime = input.finalize(comp_player);
                    auto r_prime = input.finalize(gen_player);
                    S[info.dest_base + i] = c_prime - r_prime;

                    auto c_dprime = input.finalize(comp_player);
                    auto r_msb = input.finalize(gen_player);
                    S[info.dest_base + i] += ((r_msb + c_dprime)
                            << (info.k - info.m));
                    prepare_mul(r_msb, c_dprime);
                }
            }

        exchange();

        for (auto info : infos)
            for (int i = 0; i < size; i++)
                if (info.small_gap())
                    S[info.dest_base + i] -= finalize_mul()
                            << (info.k - info.m + 1);
    }
#endif
}

#endif
