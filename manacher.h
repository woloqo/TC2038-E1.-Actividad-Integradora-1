#pragma once
#include <string>
#include <vector>
#include <utility>
using namespace std;

// Devuelve (inicio, longitud) del palíndromo más largo, inicio en base 0
inline pair<int,int> manacher(const string& o) {
    string s;
    s.push_back('@');
    s.push_back('$');
    for(auto a: o){
        s.push_back(a);
        s.push_back('$');
    }
    s.push_back('#');

    int n = s.length();
    vector<int> p(n);

    int centro = 0, limite = 0, pmax = 0, ipmax = 0;
    for(int i=1; i<n-1; i++){
        if(i < limite){
            int simetrica = 2 * centro - i;
            p[i] = min(limite-i, p[simetrica]);
        }
        int gap = p[i] + 1;
        while(s[i-gap] == s[i+gap]){
            p[i]++;
            gap++;
        }
        if(i+p[i] > limite){
            limite = i + p[i];
            centro = i;
        }
        if(p[i] > pmax){
            pmax = p[i];
            ipmax = i;
        }
    }
    return {(ipmax - pmax - 1) / 2, pmax};
}