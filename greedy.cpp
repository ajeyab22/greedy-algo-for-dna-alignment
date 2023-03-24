#include <iostream>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstring>
#include <fstream>

using namespace std;
using namespace std::chrono;

const int INF = -1e9;
double SPrime(int i, int j, int d, int mat, int mis) {
    return ((i + j) * mat / 2) - d * (mat - mis);
}

int main() {
    string a,b,line;
    ifstream in_file("g1.fasta");
    if (in_file.is_open())
    {
        getline(in_file, line);
        while (getline(in_file, line) )
        {
            a += line;
        }
        in_file.close();
    }
    ifstream in_file_2("g2.fasta");
    if (in_file_2.good())
    {
        getline(in_file_2, line);
        while (getline(in_file_2, line) )
        {
            b += line;
        }

        in_file_2.close();
    }

    int i = 0;

    int M = a.length();
    int N = b.length();

    int mat = 2;
    int mis = -1;
    int ind = -2;

    double X = 5;

    double TPrime = 0;
    double T = 0;
    auto start =steady_clock::now();

    double** R = (double**)malloc((max(M, N) + 1) * sizeof(double*));
    for (int i = 0; i < (max(M, N) + 1); i++) {
        R[i] = (double*)malloc((2*max(M, N) + 1) * sizeof(double));     
    }
    double T_arr[max(M, N) + 1];
    i=0;
    while ((i < min(M, N)) && (a[i] == b[i])) {
        i++;
    }
    R[0][max(M, N)] = i;
    TPrime = SPrime(i, i, 0,mat, mis);
    T_arr[0] = TPrime;
    int d = 0;
    int L = 0;
    int U = 0;
    cout<<"Greedy algorithm"<<endl;
    do {
        
        d++;
        int dprime = max((int) (d - floor((X + mat / 2) / (mat - mis)) - 1), 0);
        for (int k = L - 1; k <= U + 1; k++) {
            int firstI = INF;
            int secondI = INF;
            int thirdI = INF;

            if (L < k) {
                firstI = R[d - 1][k - 1+max(M, N)] + 1;
            }
            if ((L <= k) && (k <= U)) {
                secondI = R[d - 1][k+max(M, N)] + 1;
            }
            if (k < U) {
                thirdI = R[d - 1][k + 1+max(M, N)];
            }

            i = max(firstI, max(secondI, thirdI));

            int j = i - k;
            if ((i > INF) && (SPrime(i, j, d,mat, mis) >= (T_arr[dprime] - X))) {
                while ((i < M - 1) && (j < N - 1) && (a[i] == b[j])) {
                    i++;
                    j++;
                }
                R[d][k+max(M, N)] = i;
                TPrime = max(TPrime, SPrime(i, j, d,mat, mis));
            } else {
                R[d][k+max(M, N)] = INF;
            }
        }
        T_arr[d] = TPrime;
        for (int k = 0; k <= max(M, N); k++) {
            if (R[d][k+max(M, N)] > INF) {
                L = k;
                break;
            }
        }
        for (int k = max(M, N); k >= 0; k--) {
            if (R[d][k+max(M, N)] > INF) {
                U = k;
                break;
            }
        }

        for(int k = max(M, N); k >= 0; k--){
            if(R[d][k+max(M, N)] == N + k){
                L = max(L, k + 2);
                break;
            }
        }
        for(int k = 0; k <= max(M, N); k++){
            if(R[d][k+max(M, N)] == M + N - k){
                U = min(U, k - 2);
                break;
            }
        }

    } while ((L > U + 2));
    auto stop =steady_clock::now();
    auto totime=duration_cast<milliseconds>(stop - start).count();
    cout<<"Optimal score:"<<TPrime << endl;
    cout<<"Time taken in millisec:"<<totime<<endl;
}
