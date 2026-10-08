#pragma once
#include <string>
#include <vector>
#include <utility>
#include <unordered_map>
#include <fstream>
#include <iostream>
using namespace std;

// Codón -> aminoácido. STOP = '*'
inline unordered_map<string,char> construir_tabla_codones() {
    const string bases = "TCAG";
    const string aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

    unordered_map<string,char> tabla;
    int idx = 0;
    for (char a : bases)
        for (char b : bases)
            for (char c : bases) {
                string codon;
                codon += a; codon += b; codon += c;
                tabla[codon] = aa[idx++];
            }
    return tabla;
}

// Agrupa de 3 en 3 desde la posición 'marco' (0, 1 o 2) y traduce
inline string traducir(const string& seq, int marco) {
    static const unordered_map<string,char> tabla = construir_tabla_codones();
    string resultado;
    for(size_t i = marco; i +3 <= seq.size(); i += 3){
        string codon = seq.substr(i,3);
        auto it = tabla.find(codon);
        resultado += (it != tabla.end()) ? it->second : 'X';
    }
    return resultado;
}

// Devuelve (nombre, aminoácidos) por cada proteína de seq-proteins.txt
inline vector<pair<string,string>> leer_proteinas(const string& ruta) {
    vector<pair<string,string>> proteinas;
    ifstream archivo(ruta);
    if (!archivo.is_open()){
        cerr << "No se pudo abrir" <<ruta<<"\n";
        return proteinas;
    }
    string linea;
    while (getline(archivo, linea)) {
        while (!linea.empty() && (linea.back() == '\r' || linea.back() == ' '))
            linea.pop_back();
        if (linea.empty()) continue;
        if (linea[0] == '>') {
            size_t k = 0;
            while (k < linea.size() && linea[k] == '>') k++;
            proteinas.push_back({linea.substr(k), ""});
        } else if (!proteinas.empty()) {
            proteinas.back().second += linea;
        }
    }
    return proteinas;
}