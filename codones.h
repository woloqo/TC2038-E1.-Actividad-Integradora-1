#pragma once
#include <string>
#include <vector>
#include <utility>
#include <unordered_map>
using namespace std;

// Codón -> aminoácido. STOP = '*'
inline unordered_map<string,char> construir_tabla_codones() {
    // TODO Persona C
    return {};
}

// Agrupa de 3 en 3 desde la posición 'marco' (0, 1 o 2) y traduce
inline string traducir(const string& seq, int marco) {
    // TODO Persona C
    return "";
}

// Devuelve (nombre, aminoácidos) por cada proteína de seq-proteins.txt
inline vector<pair<string,string>> leer_proteinas(const string& ruta) {
    // TODO Persona C
    return {};
}