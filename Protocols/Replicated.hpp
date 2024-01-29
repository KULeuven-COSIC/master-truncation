/*
 * Replicated.cpp
 *
 */

#ifndef PROTOCOLS_REPLICATED_HPP_
#define PROTOCOLS_REPLICATED_HPP_

#include <cmath>

#include "Replicated.h"
#include "Processor/Processor.h"
#include "Processor/TruncPrTuple.h"
#include "Tools/benchmarking.h"
#include "Tools/Bundle.h"

#include "ReplicatedInput.h"
#include "Rep3Share2k.h"

#include "ReplicatedPO.hpp"
#include "Math/Z2k.hpp"

template<class T>
ProtocolBase<T>::ProtocolBase() :
        trunc_pr_counter(0), rounds(0), trunc_rounds(0), dot_counter(0),
        bit_counter(0), counter(0)
{
}

template<class T>
Replicated<T>::Replicated(Player& P) : ReplicatedBase(P)
{
    assert(T::vector_length == 2);
}

template<class T>
Replicated<T>::Replicated(const ReplicatedBase& other) :
        ReplicatedBase(other)
{
}

inline ReplicatedBase::ReplicatedBase(Player& P) : P(P)
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
inline ReplicatedBase::ReplicatedBase(Player& P, array<PairwisePRNG, 2>& prngs) :
#else
inline ReplicatedBase::ReplicatedBase(Player& P, array<PRNG, 2>& prngs) :
#endif
        P(P)
{
    for (int i = 0; i < 2; i++)
        shared_prngs[i].SetSeed(prngs[i]);
}

inline ReplicatedBase ReplicatedBase::branch() const
{
    return {P, shared_prngs};
}

template<class T>
ProtocolBase<T>::~ProtocolBase()
{
#ifdef VERBOSE_COUNT
    if (counter or rounds)
        cerr << "Number of " << T::type_string() << " multiplications: "
                << counter << " (" << bit_counter << " bits) in " << rounds
                << " rounds" << endl;
    if (counter or rounds)
        cerr << "Number of " << T::type_string() << " dot products: " << dot_counter << endl;
    if (trunc_pr_counter or trunc_rounds)
        cerr << "Number of probabilistic truncations: " << trunc_pr_counter << " in " << trunc_rounds << " rounds" << endl;
#endif
}

template<class T>
void ProtocolBase<T>::mulrs(const vector<int>& reg,
        SubProcessor<T>& proc)
{
    proc.mulrs(reg);
}

template<class T>
void ProtocolBase<T>::multiply(vector<T>& products,
        vector<pair<T, T> >& multiplicands, int begin, int end,
        SubProcessor<T>& proc)
{
#ifdef VERBOSE_CENTRAL
    fprintf(stderr, "multiply from %d to %d in %d\n", begin, end,
            BaseMachine::thread_num);
#endif

    init(proc.DataF, proc.MC);
    init_mul();
    for (int i = begin; i < end; i++)
        prepare_mul(multiplicands[i].first, multiplicands[i].second);
    exchange();
    for (int i = begin; i < end; i++)
        products[i] = finalize_mul();
}

template<class T>
T ProtocolBase<T>::mul(const T& x, const T& y)
{
    init_mul();
    prepare_mul(x, y);
    exchange();
    return finalize_mul();
}

template<class T>
void ProtocolBase<T>::prepare_mult(const T& x, const T& y, int n,
		bool)
{
    prepare_mul(x, y, n);
}

template<class T>
void ProtocolBase<T>::finalize_mult(T& res, int n)
{
    res = finalize_mul(n);
}

template<class T>
T ProtocolBase<T>::finalize_dotprod(int length)
{
    counter += length;
    dot_counter++;
    T res;
    for (int i = 0; i < length; i++)
        res += finalize_mul();
    return res;
}

template<class T>
T ProtocolBase<T>::get_random()
{
    if (random.empty())
    {
        buffer_random();
        assert(not random.empty());
    }

    auto res = random.back();
    random.pop_back();
    return res;
}

template<class T>
vector<int> ProtocolBase<T>::get_relevant_players()
{
    vector<int> res;
    int n = dynamic_cast<typename T::Protocol&>(*this).P.num_players();
    for (int i = 0; i < T::threshold(n) + 1; i++)
        res.push_back(i);
    return res;
}

template<class T>
void Replicated<T>::init_mul()
{
    for (auto& o : os)
        o.reset_write_head();
    add_shares.clear();
}

template<class T>
void Replicated<T>::prepare_mul(const T& x,
        const T& y, int n)
{
    typename T::value_type add_share = x.local_mul(y);
    prepare_reshare(add_share, n);
}

template<class T>
void Replicated<T>::prepare_reshare(const typename T::clear& share,
        int n)
{
    typename T::value_type tmp[2];
    for (int i = 0; i < 2; i++)
        tmp[i].randomize(shared_prngs[i], n);
    auto add_share = share + tmp[0] - tmp[1];
    add_share.pack(os[0], n);
    add_shares.push_back(add_share);
}

template<class T>
void Replicated<T>::exchange()
{
    os[0].append(0);
    if (os[0].get_length() > 0)
        P.pass_around(os[0], os[1], 1);
    this->rounds++;
}

template<class T>
void Replicated<T>::start_exchange()
{
    os[0].append(0);
    P.send_relative(1, os[0]);
    this->rounds++;
}

template<class T>
void Replicated<T>::stop_exchange()
{
    P.receive_relative(-1, os[1]);
}

template<class T>
inline T Replicated<T>::finalize_mul(int n)
{
    this->counter++;
    this->bit_counter += n;
    T result;
    result[0] = add_shares.next();
    result[1].unpack(os[1], n);
    return result;
}

template<class T>
inline void Replicated<T>::init_dotprod()
{
    init_mul();
    dotprod_share.assign_zero();
}

template<class T>
inline void Replicated<T>::prepare_dotprod(const T& x, const T& y)
{
    dotprod_share = dotprod_share.lazy_add(x.local_mul(y));
}

template<class T>
inline void Replicated<T>::next_dotprod()
{
    dotprod_share.normalize();
    prepare_reshare(dotprod_share);
    dotprod_share.assign_zero();
}

template<class T>
inline T Replicated<T>::finalize_dotprod(int length)
{
    (void) length;
    this->dot_counter++;
    return finalize_mul();
}

template<class T>
T Replicated<T>::get_random()
{
    T res;
    for (int i = 0; i < 2; i++)
        res[i].randomize(shared_prngs[i]);
    return res;
}

template<class T>
void ProtocolBase<T>::randoms_inst(vector<T>& S,
		const Instruction& instruction)
{
    for (int j = 0; j < instruction.get_size(); j++)
    {
        auto& res = S[instruction.get_r(0) + j];
        randoms(res, instruction.get_n());
    }
}

template<class T>
void Replicated<T>::randoms(T& res, int n_bits)
{
    for (int i = 0; i < 2; i++)
        res[i].randomize_part(shared_prngs[i], n_bits);
}

// void print_os(const octetStream &os) {
//   std::string content = os.str();
//   for(char c : content) {
//     std::cout << std::hex << (int)(c);
//   }
// }

template<class T>
template<class U>
void Replicated<T>::trunc_pr(const vector<int>& regs, int size, U& proc,
        false_type)
{
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

    // std::cout<<"truncpr"<<std::endl;

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
                    y[0] = -value_type(-value_type(x[0]) >> info.m);
                    y[1] = r_prime;

                    //CHECK PREP
                    Z2 r_pp = this->shared_prngs[0].template get<value_type>();
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

template<class T>
template<class U>
void Replicated<T>::trunc_pr(const vector<int>& regs, int size, U& proc,
        true_type)
{
    // (void) regs, (void) size, (void) proc;
    // throw runtime_error("trunc_pr not implemented");
    this->trunc_rounds++;
    trunc_pr(regs, size, proc, T::clear::characteristic_two);
}

template<class T>
template<class U>
void Replicated<T>::trunc_pr(const vector<int>& regs, int size,
        U& proc)
{
    this->trunc_rounds++;
    trunc_pr(regs, size, proc, T::clear::characteristic_two);
}

template<class T> // this is Rep3Share2<BITSIZE>> where BITSIZE in {64, 72, 128, 192}
template<class U> // this is SubProcessor<Rep3Share2<BITSIZE>
void Replicated<T>::relu(U &proc, int reg0, int reg1, int k) {
  assert(proc.P.num_players() == 3);
  assert(proc.Proc != 0);
  using Z2 = typename T::T;
  // storage of secret ints S
  auto& S = proc.get_S();
  // address of result: reg0
  // address of input: reg1
  int f = 16;
  const T &input = S[reg1];
  T &output = S[reg0];
  // placeholder: copy input to output
  //std::cout << "my P number: " << P.my_num() << "\n" ;
  //std::cout << "Input: " << input << " " << input.v[0] << "\n";

  int eval_player = 2;
  bool evaluator = P.my_num() == eval_player; 
  int other_player;
  if (P.my_num() == 0) 
      { other_player = 1; } 
  else 
      { other_player = 0; }
  int Drelu_prime;
  T Drelu;
  GlobalPRNG G(P);

  if (!evaluator){
      octetStream os;
      //GlobalPRNG G(P);
      PairwisePRNG PR(proc.P, other_player);
      // Generate random bit t and mask input
      bool t = PR.get_bit();
      //std::cout << "Bit t: " << t << "\n";
      T x;
      if(t) {x.v[0] = Z2(0)-input.v[0]; x.v[1] = Z2(0)-input.v[1];} else { x = input; }
      
         
      // Generate random vector r
      std::vector<int> r;
      for (size_t i = 0; i < size_t(k+2); i++) {
          r.push_back(PR.get_uint()); 
      }

      if (P.my_num() == 1)
	  { x.v[0] += x.v[1]; }

      // Compute u
      std::vector<Z2> u(r.size());
      u[0] = (t) ? Z2(-1) : Z2(1);
      u[1] = x.v[0];
      for (size_t i = 2; i < r.size(); i++) {
          u[i] = x.v[0] >> (i-1);
	  //std::cout << "u[i]: " << u[i] << "\n";
      }

      // Compute v
      std::vector<Z2> v(r.size());
      v[0] = u[0] + Z2(3)*u[1] - Z2(1); //-T::constant(1,P.my_num()); 
      for (size_t i = 1; i < r.size(); i++) {
 	  Z2 u_sum;
	  for (size_t j = i; j < u.size(); j++) {
	      u_sum += u[j]; 
          }
	  //v[i] = u_sum - T::constant(1,P.my_num());
	  v[i] = u_sum - Z2(1);
	  //std::cout << "v[i]: " << v[i] << "\n";
      }

      // Compute w  
      std::vector<T> w(r.size());
      for (size_t i = 0; i < r.size(); i++) {
          w[i].v[0] = Z2(r[i])*v[i];
	  //std::cout << "w[i]: " << w[i] << "\n";	
      }
       
      //uint g = G.get_uint();
      //std::shuffle(std::begin(w), std::end(w), int(g));
      //std::random_shuffle(w.begin(), w.end(), G.get_uint()%i); //shuffle the output

      // Send the w vector to the evaluator
      os.store(w);
      P.send_to(eval_player, os); 
      this->rounds++; 

      // EVALUATOR DOES COMPUTATIONS //
	  
      if (P.my_num() == 0){
	      // Generate triples with P2
	      PairwisePRNG PR0(proc.P, eval_player);
	      int a0 = PR0.get_uint();
	      int b0 = PR0.get_uint();
	      int c0 = PR0.get_uint();
	      
	      Z2 e0 = input[0] - Z2(a0);
	      octetStream os;
	      os.store(e0);
	      P.send_to(1, os); 
              //this->rounds++; 

	      // Open e
	      octetStream rs;
              P.receive_player(1, rs);
	      Z2 e1;
	      rs.get(e1);
	      Z2 e = e0 + e1;

	      // Get d
	      octetStream rs1;
	      int d;
	      P.receive_player(2, rs1);
	      rs1.get(d);

	      // Compute share of relu
	      Z2 relu = Z2(1-2*int(t))*(Z2(d)*e + Z2(d*a0) + e*Z2(b0) + Z2(c0)) + Z2(int(t))*input.v[0];

	      octetStream os1;
	      os1.store(relu);
	      P.send_to(eval_player, os1); 
              //this->rounds++; 

	      octetStream rs2;
              P.receive_player(1, rs2);
	      Z2 x2;
	      rs2.get(x2);

	      Drelu.v[0] = relu;
	      Drelu.v[1] = x2;

      } else {
	      //Generate triples with P2
	      PairwisePRNG PR1(proc.P, eval_player);
	      int a1 = PR1.get_uint();
	      int b1 = PR1.get_uint();
	      
              Z2 e1 = input[0] + input[1] - Z2(a1);
	      octetStream os;
	      os.store(e1);
	      P.send_to(0, os); 
              //this->rounds++; 

	      // Open e
	      octetStream rs;
              P.receive_player(0, rs);
	      Z2 e0;
	      rs.get(e0);
	      Z2 e = e0 +e1;

	      // Get c1 and d
	      octetStream rs1;
	      std::vector<int> vals(2);
	      P.receive_player(2, rs1);
	      rs1.get(vals);
	      int c1 = vals[0];
	      int d = vals[1];

	      // Compute share of relu
	      Z2 relu = Z2(1-2*int(t))*(Z2(d*a1) + e*Z2(b1) + Z2(c1)) + Z2(int(t))*input.v[0];

	      SeededPRNG PP;
	      int r = PP.get_uint() * pow(2, f);
	      Z2 x2 = relu - Z2(r);
	      Z2 x3 = r;
	      octetStream os1;
	      os1.store(x2);
	      P.send_to(0, os1); 
              //this->rounds++; 

              octetStream os2;
	      os2.store(x3);
	      P.send_to(2, os1); 
              //this->rounds++; 

	      Drelu.v[0] = x2;
	      Drelu.v[1] = x3;
      } 
      octetStream rss;
      P.receive_player(eval_player, rss);
      rss.get(Drelu_prime);
      //std::cout << "CompRec_Drelu: " << Drelu << int(t) << " " << (Drelu+int(t))%2 << "\n";
      Drelu_prime = (Drelu_prime+int(t))%2;
  }

  if (evaluator){
      octetStream rs;
      octetStream rs2;
      P.receive_player(0, rs);
      P.receive_player(1, rs2);
      std::vector<T> w1, w2;
      rs.get(w1);
      rs2.get(w2);
          
      std::vector<Z2> w(w1.size());
      for (size_t i = 0; i < w.size(); i++) {
          //w[i] = w1[i] + w2[i];
	  w[i] = (w1[i].v[0] + w2[i].v[0]);	      
      }

      if (std::count(w.begin(), w.end(), Z2(0))) {
          Drelu_prime = 1;
      }
      else {
          Drelu_prime = 0;
      }

      
      //std::cout << "Eval_Drelu: " << Drelu << "\n";

      PairwisePRNG PR0(proc.P, 0);
      int a0 = PR0.get_uint();
      int b0 = PR0.get_uint();
      int c0 = PR0.get_uint();

      PairwisePRNG PR1(proc.P, 1);
      int a1 = PR1.get_uint();
      int b1 = PR1.get_uint();

      // Compute c1 and d and send them o comp. parties
      int c1 = (a0 + a1)*(b0 + b1) - c0;
      int d = Drelu_prime - (b0 + b1);

      octetStream os;
      os.store(d);
      P.send_to(0, os); 
      //this->rounds++; 

      octetStream os1;
      std::vector<int> vals(2);
      vals[0] = c1;
      vals[1] = d;
      os1.store(vals);
      P.send_to(1, os1); 
      //this->rounds++; 

//
      octetStream os2;
      os2.store(Drelu_prime);
      P.send_all(os2);
      //this->rounds++; 


      octetStream rs3;
      P.receive_player(0, rs3);
      Z2 x1;
      rs3.get(x1);

      octetStream rs4;
      P.receive_player(1, rs4);
      Z2 x3;
      rs4.get(x3);

      Drelu.v[0] = x3;
      Drelu.v[1] = x1;

  }

  init_mul();
  prepare_mul(input, Drelu);
  exchange();
  T a = finalize_mul();
  std::cout << "fin-mul" << a << std::endl;
  output.v[0] = input.v[0]*Z2(Drelu_prime);
  output.v[1] = input.v[1]*Z2(Drelu_prime);

}


#endif
