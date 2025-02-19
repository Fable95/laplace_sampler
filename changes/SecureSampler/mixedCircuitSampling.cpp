/*
 * mixed-example.cpp
 *
 */

#include "Protocols/ProtocolSet.h"
#include <numeric>
#include <chrono>
#include "Machines/maximal.hpp"

#define VERBOSE 1
int dabits_count = 0;
template<class T>
void run(char** argv);

int main(int argc, char** argv){
    // need player number and number of players
    if (argc != 5){
        cerr << "Usage: " << argv[0] << "<my number: 0/1/...> <total number of players> <sec param> <num samples>" << endl;
        exit(1);
    }    
    // cout << "----------------    SEMI   SHARE    ----------------" << endl;
    // run<SemiShare<gfp_<0,2>>>(argv);
    // cout << "----------------    SEMI2k SHARE    ----------------" << endl;
    // run<Semi2kShare<64>>(argv);
    cout << "----------------    SPDZ2k SHARE    ----------------" << endl;
    run<Spdz2kShare<64, 64>>(argv);
    // cout << "----------------    MASCOT SHARE    ----------------" << endl;
    // run<Share<gfp_<0, 2>>>(argv);
}

template<class T>
vector<T> get_binary_representation(double alpha, int d, const typename T::mac_key_type &mac_key, int my_num){
    if(alpha >= 1.0 || alpha <= 0.0){
        throw runtime_error("wrong bernoulli probability");
    }
    vector<T> representation(d);
    int n_bits = 1;
    // Numer is in 0 < x < 1
    for (int i = 0; i < d; i++){
        alpha *= 2;
        if(alpha >= 1.0){
            representation.at(i) = T::constant(1, my_num, mac_key, n_bits);
            alpha -= 1;
        } else {
            representation.at(i) = T::constant(0, my_num, mac_key, n_bits);
        }
    }
    return representation;
}

template<class T, class V>
void open_print_vector(const vector<T>& vec, Player& P, string name, V& output){
    int size = vec.size();
    output.init_open(P, size);
    for (int i = 0; i < size; i++){
        output.prepare_open(vec.at(i));
    }
    output.exchange(P);
    cout << name  << " = [";
    for (int i = 0; i < size; i++){
        if(i!=0) std::cout << ", ";
        cout << output.finalize_open();
    }
    cout << "]" << endl;
}

// First implement naive, later implement optimized prefix-or with the methor of Chandra Fortune and Lipton
template<class T>
void PRE_OR(const vector<typename T::bit_type>& op1, vector<typename T::bit_type>& res, int d, Player& P, MixedProtocolSet<T>& set){
    (void)P;
    int size = op1.size();
    if(size % d != 0){
        throw runtime_error("Op1 size has to be multiple of batch size d");
    }
    int batches = size / d;
    auto& bit_protocol = set.binary.protocol;
    res.resize(size);
    if(op1.empty()){
        throw runtime_error("Cannot run Pre_or with empty input");
    }
    for (int i = 0; i < batches; i++)
    {
        res[i * d] = op1[i * d];
    }
    
    for (int i = 1; i < d; i++){
        bit_protocol.init_mul();    
        for (int j = 0; j < batches; j++){
            bit_protocol.prepare_mul(op1[j*d + i], res[j*d + i - 1], 1);
            res[j*d + i] = op1[j*d + i] + res[j*d + i -1];
        }
        bit_protocol.exchange();
        for (int j = 0; j < batches; j++){
            res[j*d + i] += bit_protocol.finalize_mul(1);    
        }
    }
    // open_print_vector<T, typename T::bit_type::MAC_Check>(op1, P, "a       ", set.binary.output);
    // open_print_vector<T, typename T::bit_type::MAC_Check>(res, P, "PRE(a)  ", set.binary.output);
    
}

template<class T>
void insert_samples(vector<typename T::bit_type>& full_vector, vector<typename T::bit_type>& to_insert){
    if(to_insert.size() > full_vector.size())
        throw runtime_error("Error cannot insert larger vector without free places");
    if(full_vector.size()%to_insert.size() != 0)
        throw runtime_error("Small vector has to divide large vector");
    size_t index_multiplier = (full_vector.size() / to_insert.size());
    for (size_t i = 0; i < to_insert.size(); i++){
        size_t index = i * index_multiplier;
        full_vector.at(index) = to_insert.at(i);
    }
}

template<class T>
vector<typename T::bit_type> approx_bernoulli(int n_parties, int N, int d, vector<typename T::bit_type>& alpha_l,
                      MixedProtocolSet<T>& set, MixedProtocolSetup<T>& setup, Player& P, bool first_empty=false, int partition=1){
    auto& bit_input = set.binary.input;
    auto& bit_protocol = set.binary.protocol;
    
    const auto& binary_mac = setup.binary.get_mac_key();
    vector<typename T::bit_type> u_l(d * N);
    vector<typename T::bit_type> c_l(d * N);
    vector<typename T::bit_type> e_l(d * N);
    typename T::bit_type f_l;
    
    vector<typename T::bit_type> prod(N);
    if(N % partition != 0){
        throw runtime_error("Error Cannot partition if partition count does not divide number of samples");
    }
    int len_partition = 1 + N/partition;
    if(first_empty)
        prod.resize(N+partition);
    SeededPRNG g;
    vector<vector<int>> s(N);

    for (size_t i = 0; i < prod.size(); i++){
        prod.at(i) = T::bit_type::constant(1, P.my_num(), binary_mac, 1);
    }

    // inputs in binary domain
    bit_input.reset_all(P);
    for (int num = 0; num < N; num++){
        for (int i = 0; i < d; i++){
            bit_input.add_from_all(g.get_bit(), 1);
        }
    }
    bit_input.exchange();
    for (int num = 0; num < N; num++){
        for (int i = 0; i < d; i++){
            u_l.at((num * d) + i) = bit_input.finalize(0, 1);
            for (int j = 1; j < n_parties; j++){
                u_l.at((num * d) + i) += bit_input.finalize(j, 1);    
            }
            c_l.at((num * d) + i) = u_l.at((num * d) + i) + alpha_l.at(i);
        }
    }
    PRE_OR<T>(c_l, e_l, d, P, set);
    bit_protocol.init_mul();
    for (int num = 0; num < N; num++)
    {
        f_l = e_l.at(num * d);
        bit_protocol.prepare_mul(f_l, u_l.at(num * d), 1);
        for (int i = 1; i < d; i++){    
            f_l = e_l.at(num * d + i) + e_l.at(num * d + i - 1);
            bit_protocol.prepare_mul(f_l, u_l.at(num * d + i), 1);    
        }
    }
    bit_protocol.exchange();
    
    for (size_t num = 0; num < prod.size(); num++){
        // This creates an empty slot at the beginning of each partition
        if(first_empty && (num % len_partition == 0))
            num++;
        for (int i = 0; i < d; i++){
            prod.at(num) += bit_protocol.finalize_mul(1);
        }
    }
    bit_protocol.check();
    // string result_name = string("bernoulli (") + to_string(alpha) + string(")");
    // open_print_vector<typename T::bit_type, decltype(bit_output)>(prod, P, result_name, set.binary.output);
    return prod;
    
}


template<class T>
vector<T> transform_share(vector<typename T::bit_type>& bit_vector, MixedProtocolSet<T>& set, MixedProtocolSetup<T>& setup, Player& P, vector<pair<T, typename T::bit_type>>& dabits){
    auto& bit_output = set.binary.output;
    auto& prep = set.preprocessing;
    bool dabit_necessary = dabits.empty();
    int size = bit_vector.size();
    bit_output.init_open(P, size);
    pair<T, typename T::bit_type> dabit;

    if(dabit_necessary){
        for (int i = 0; i < size; i++){
            dabits.push_back({});
            dabit = dabits.back();
            prep.get_dabit(dabit.first, dabit.second);
            dabits_count++;
        }
    }

    for (int i = 0; i < size; i++){
        
        dabit = dabits.at(dabits.size() - i - 1);
        bit_output.prepare_open(typename T::bit_type::part_type(
                                dabit.second.get_bit(0) + bit_vector.at(i).get_bit(0)));    
    }
    bit_output.exchange(P);
    vector<T> res(size);
    for (int i = 0; i < size; i++){ 
        typename T::clear masked = bit_output.finalize_open().get_bit(0);
        auto mask = dabits.back().first;
        res.at(i) = (mask - mask * masked * 2 + T::constant(masked, P.my_num(), setup.get_mac_key()));
        dabits.pop_back();
    }   
    // open_print_vector<typename T::bit_type, decltype(bit_output)>(bit_vector, P, "bit vector", bit_output);
    // open_print_vector<T, decltype(set.output)>(res, P, "Zp vector ", set.output);
    return res;
}

template<class T>
void update_time(T& start, chrono::microseconds& total_time, string name){
    auto end_time = std::chrono::high_resolution_clock::now();
    auto diff = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start);
    (void)name;
    // cout << name << diff.count() << " ms" << std::endl;
    total_time += diff;
    start = std::chrono::high_resolution_clock::now();
}

template<class T>
vector<T> FDL(int n_parties, int N, int d, double p, int num_iterations,
    MixedProtocolSet<T>& set, MixedProtocolSetup<T>& setup, Player& P, vector<pair<T, typename T::bit_type>>& da_bits){
    SeededPRNG g;
    double p0 = (1.0 - p) / (1.0 + p);
    double p1 = 1.0 - p;
    auto& bit_input = set.binary.input;
    auto& protocol = set.protocol;
    vector<T> k_vec(num_iterations);
    vector<T> s_vec(num_iterations);
    
    vector<typename T::bit_type> c_vector(N * num_iterations);
    vector<typename T::bit_type> sigma_vector(num_iterations);
        
    auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = start_time;
    chrono::microseconds total_time = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    vector<typename T::bit_type> alpha_0 = get_binary_representation<typename T::bit_type>(p0, d, setup.binary.get_mac_key(), P.my_num());
    vector<typename T::bit_type> alpha_1 = get_binary_representation<typename T::bit_type>(p1, d, setup.binary.get_mac_key(), P.my_num());
    vector<typename T::bit_type> b0 = approx_bernoulli(n_parties, num_iterations, d, alpha_0, set, setup, P);
    update_time(start_time, total_time, string("Bernoulli(") + to_string(num_iterations) + string("): "));
    vector<typename T::bit_type> bi = approx_bernoulli(n_parties, num_iterations * (N-1), d, alpha_1, set, setup, P, true, num_iterations);
    update_time(start_time, total_time, string("Bernoulli(") + to_string(num_iterations * (N-1)) + string("): "));
    insert_samples<T>(bi, b0);
    PRE_OR(bi, c_vector, N, P, set);
    update_time(start_time, total_time,"PREFIX_OR: ");
    bit_input.reset_all(P);
    for (int num = 0; num < num_iterations; num++)
        bit_input.add_from_all(g.get_bit(), 1);
    bit_input.exchange();
    for (int num = 0; num < num_iterations; num++){
        sigma_vector.at(num) =  bit_input.finalize(0, 1);
        for (int j = 1; j < n_parties; j++){
            sigma_vector.at(num) += bit_input.finalize(j, 1);    
        }
    }
    update_time(start_time, total_time, "input sigma: ");
    auto c_Zp_vector = transform_share(c_vector, set, setup, P, da_bits);
    auto sigma_Zp_vector = transform_share(sigma_vector, set, setup, P, da_bits);
    std::cout << "after usage " << da_bits.size() << " dabits remain" << std::endl; 
    update_time(start_time, total_time,"Transform Share: ");
    protocol.init_mul();
    for (int i = 0; i < num_iterations; i++){
        auto l = T::constant(N, P.my_num(), setup.get_mac_key());
        for (int j = 0; j < N; j++){
            l -= c_Zp_vector.at(i*N + j);
        }
        auto& sigma_Zp = sigma_Zp_vector.at(i);
        sigma_Zp *= -2;
        sigma_Zp += T::constant(1, P.my_num(), setup.get_mac_key());
        protocol.prepare_mul(sigma_Zp, l);
    }
    protocol.exchange();
    for (int i = 0; i < num_iterations; i++){
        k_vec.at(i) = protocol.finalize_mul();
    }
    update_time(start_time, total_time,"s, l and k vec: ");    
    cout << "Online time taken:   " << total_time.count() * 1e-3 << " ms" << std::endl;
    // cout << "average time taken: " << static_cast<double>(total_time.count())/num_iterations << "ms" << std::endl;
    // cout << "total time taken:   " << total_time.count() / 1000 << " s" << std::endl;

    // cout << "multiplication time: " << total_time.count() << "ms" << std::endl;
    // open_print_vector<T, decltype(set.output)>(s_vec, P, "s", set.output);
    // open_print_vector<T, decltype(set.output)>(l_vec, P, "l", set.output);
    open_print_vector<T, decltype(set.output)>(k_vec, P, "k", set.output);
    return k_vec;
}

template<class T>
void run(char** argv)
{
    
    double p = 0.5;
    // int d = 32;
    // int LaplaceWidth = 32;
    // double p0 = (1.0 - p) / (1.0 + p);
    // double p1 = 1.0 - p0;
    // reduce batch size
    // vector<int> bin_rep = get_binary_representation(0.75, 5);

    // set up networking on localhost
    int my_number = atoi(argv[1]);
    int n_parties = atoi(argv[2]);
    int security_parameter   = atoi(argv[3]);
    int NumberSamples = atoi(argv[4]);
    int port_base = 9999;
    
    OnlineOptions::singleton.bucket_size = 5;
    OnlineOptions::singleton.batch_size = (NumberSamples < 64) ? 64 : NumberSamples;

    Names N(my_number, n_parties, "localhost", port_base);
    CryptoPlayer P(N);

    // protocol setup (domain, MAC key if needed etc)
    MixedProtocolSetup<T> setup(P);
    
    // set of protocols (bit_input, multiplication, output)
    MixedProtocolSet<T> set(P, setup);
    
    auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = start_time;
    // Preprocessing
    set.preprocessing.buffer_triples();
    set.preprocessing.buffer_bits();
    vector<pair<T, typename T::bit_type>> dabits(NumberSamples*(security_parameter+1));
    // vector<pair<T, typename T::bit_type>> dabits(0);
    for (size_t i = 0; i < dabits.size(); i++){
        auto& dabit = dabits.at(i);
        set.preprocessing.get_dabit(dabit.first, dabit.second);
    }
    end_time = std::chrono::high_resolution_clock::now();
    chrono::microseconds total_time = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time); 
    auto network = P.total_comm().sent;
    cout << "Preprocessing time taken: " << total_time.count() * 1e-3 << " ms" << std::endl;
    cout << "Preprocessing data transmission: "  << (network) * 1e-6
                << " MB (this party only)" << endl;
    cout << endl << "Security Parameter " << security_parameter << ":" << endl;

    vector<T> k = FDL(n_parties, security_parameter, security_parameter, p, NumberSamples, set, setup, P, dabits);
    
    std::cout << "dabit count: " << dabits_count << std::endl;
    cout << "Online data transmission: "  << (P.total_comm().sent - network) * 1e-6
                << " MB (this party only)" << endl;
    set.check();
    
}
