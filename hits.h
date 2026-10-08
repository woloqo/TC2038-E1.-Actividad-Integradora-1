// hits.h
#pragma once
#include <string>
#include <vector>
#include <utility>
using namespace std;

struct Hit { string nombre; int inicio; int fin; int marco; };  // inicio y fin en base 0

// TODO Persona C: implementación real (esta versión vacía es temporal)
inline vector<Hit> buscar_proteinas(const string& genoma,
                                    const vector<pair<string,string>>& proteinas) {
    return {};
}