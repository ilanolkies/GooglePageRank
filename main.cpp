#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <cmath>

using namespace std;

#define tolerance 0.000001

int main(int argc, char *argv[]) {
    if (argc != 3) {
        cout << "El programa debe cumplir el siguiente formato: ./tp1 archivo p" << std::endl;
        return 1;
    }

    string input_file = argv[1];
    float p = atof(argv[2]);


    fstream input(input_file, ios_base::in);

    // W ∈ Rn×n -> wij = j->i = 1 : 0
    // l = sum ij wij != 0 ? 1 : 0

    uint n, l;
    input >> n >> l;

    vector<map<uint, double>> W(n, map<uint, double>());

    int j, i;

    for (int k = 0; k < l; k++) {
        input >> j >> i;
        W[i - 1][j - 1] = 1;
    }

    input.close();

    // C ∈ Rn -> cj = sum i wij
    vector<double> c(n , 0);

    for (i = 0; i < n; ++i)
        for (map<uint, double>::iterator col_it = W[i].begin(); col_it != W[i].end(); col_it++)
            c[i] += col_it->second;

    // A = I−pWD
    vector<map<uint, double>> A(n, map<uint, double>());

    // WD ∈ Rn×n -> wdij cj != 0 ? wij/cj : 0
    // pWD ∈ Rn×n -> wdij cj != 0 ? p wij/cj : 0
    // I-pWD ∈ Rn×n
    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j) {
            if (c[j] != 0) A[i][j] = p * W[i][j] / c[j];
            if (i == j) A[i][j] = 1 - A[i][j];
        }

    vector<double> e(n, 1);

    vector<double> ranking(n, 0);

    double A_jj, A_ij;
    short status = 0;

    // Iteracion sobre las columnas de A, excepto la ultima porque ahi no hay que hacer nada
    for (j = 0; j < n - 1; ++j) {
        map<uint, double>::iterator it_fila_j = A[j].begin();

        i = j;
        //Mientras el i-esimo elemento no nulo de la fila i no está en la columna j...
        while (i < n && (it_fila_j == A[j].end() || it_fila_j->first != j))
            it_fila_j = A[++i].begin();     //avanzo a la siguiente fila.
            //Si encontre una fila con elemento no nulo en la columna j...
            if (i < n) {
                //Cambio de lugar las filas (para que no haya un 0 en la diagonal).
                A[j].swap(A[i]);
                //En consecuencia debo cambiar también el orden de e.
                double bb_m = e[j];
                e[j] = e[i];
                e[i] = bb_m;
                //Debido al swap, it_fila_j es iterador de la fila j.
                A_jj = it_fila_j->second;
                ++i;
                //Reviso las siguientes filas.
                while (i < n) {
                    map<uint, double>::iterator it1 = A[i].begin();
                    map<uint, double>::iterator fin1 = A[i].end();
                    if (it1 != fin1 && it1->first == j) {  //Si el elemento en la columna j no es nulo, resto filas.
                        A_ij = it1->second;
                        map<uint, double>::iterator it2 = it_fila_j;
                        ++it2;  //Voy a la siguiente columna relevante de la fila j.
                        map<uint, double>::iterator fin2 = A[j].end();
                        A[i].erase(it1++); //El elemento en la columna j debe quedar en 0. Voy a la siguiente columna relevante de la fila i (con it1).
                    while (it2 != fin2) { //Mientras no haya acabado la fila para restar:
                        double resultado_de_la_resta =-(A_ij / A_jj) * (it2->second); //empiezo restando lo que hay que restar
                        while (it1 != fin1 && it1->first < it2->first) ++it1;
                        if (it1 == fin1 || it1->first > it2->first)
                            if (abs(resultado_de_la_resta) > tolerance)
                                A[i].insert(it1, make_pair(it2->first, resultado_de_la_resta));
                        else {
                            it1->second += resultado_de_la_resta;
                            if (abs(it1->second) < tolerance)
                                A[i].erase(it1++);
                            else
                                ++it1;
                        }
                        ++it2;
                    }
                    e[i] -= (A_ij / A_jj) * e[j]; // Actualizo ranking
                }
                ++i;
            }
        }
    }
    for (long int i = n -1; i >= 0; --i) {
        auto it = A[i].begin();
        A_jj = ((it == A[i].end() || it->first != i) ? 0 : it->second);
        if (A_jj != 0) ++it;  //Si A_jj no es 0 avanzo "it", pués no forma parte de la siguiente resta.
        while (it != A[i].end()) {
            e[i] -= (it->second) * ranking[it->first];   //b_i - sum_j(A_ij*x_j)
            ++it;
        }
        if (A_jj == 0 && e[i] != 0) {
            status = -1; //el sistema es incompatible
            break;
        } else if (A_jj == 0 && e[i] == 0) {
            status = 1; //hay infinitos resultados
            ranking[i] = 0;
        } else
            ranking[i] = e[i] / A_jj;
    }

    // Normalize ranking
    vector<double> v = ranking;
    sort(v.begin(), v.end()); // Order to get a more precise sum calculation

    double s = 0;
    for (int i=0; i<v.size(); i++) s += v[i];

    for (int i=0; i<v.size(); i++) ranking[i] = ranking[i] / s;

    // Print result
    for (int i = 0; i < ranking.size(); ++i) cout << ranking[i] << "\n";

    return 0;
}
