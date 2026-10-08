#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
using namespace std;

typedef unsigned long long ull;

// Subcadena común más larga: inicio en a, inicio en b, longitud (base 0)
struct Comun { int ia; int ib; int largo; };

const ull MOD = (1ULL << 61) - 1;
const ull BASE = 911382323ULL;

inline ull mulmod(ull a, ull b) {
    __uint128_t r = (__uint128_t)a * b;
    ull x = (ull)(r & MOD) + (ull)(r >> 61);
    return x >= MOD ? x - MOD : x;
}

// Hash de todas las ventanas de largo L
inline vector<ull> hashes_ventana(const string& s, int L, ull potencia) {
    vector<ull> h;
    ull act = 0;
    for(int i=0; i<(int)s.size(); i++){
        act = (mulmod(act, BASE) + (ull)s[i]) % MOD;
        if(i >= L){
            ull quitar = mulmod((ull)s[i-L], potencia);
            act = (act + MOD - quitar) % MOD;
        }
        if(i >= L-1) h.push_back(act);
    }
    return h;
}

// Busca una subcadena común de largo L, verifica caracter a caracter
inline bool existe_comun(const string& a, const string& b, int L, int& ia, int& ib) {
    ull potencia = 1;
    for(int i=0; i<L; i++) potencia = mulmod(potencia, BASE);

    vector<ull> ha = hashes_ventana(a, L, potencia);
    vector<ull> hb = hashes_ventana(b, L, potencia);

    unordered_map<ull,int> pos;
    for(int i=0; i<(int)ha.size(); i++) pos[ha[i]] = i;

    for(int j=0; j<(int)hb.size(); j++){
        auto it = pos.find(hb[j]);
        if(it == pos.end()) continue;
        int i = it->second;
        bool igual = true;
        for(int k=0; k<L && igual; k++){
            if(a[i+k] != b[j+k]) igual = false;
        }
        if(igual){
            ia = i;
            ib = j;
            return true;
        }
    }
    return false;
}

// Búsqueda binaria sobre el largo + hash rodante
inline Comun subcadena_comun(const string& a, const string& b) {
    int lo = 1, hi = min(a.size(), b.size());
    Comun mejor = {-1, -1, 0};
    while(lo <= hi){
        int mid = (lo + hi) / 2;
        int ia, ib;
        if(existe_comun(a, b, mid, ia, ib)){
            mejor = {ia, ib, mid};
            lo = mid + 1;
        } else {
            hi = mid - 1;
        }
    }
    return mejor;
}