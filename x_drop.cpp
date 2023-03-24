#include <iostream>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstring>
#include <fstream>

using namespace std;
using namespace std::chrono;
const int INF = -1e9; 

int main() {
    double T = 0;
    double T_prime = T;
    int k = 0,j;
    int L = 0;
    int U = 0;
    double X = 5;
    double mat = 2;
    double mis = -1;
    double ind = -2;
    string sequence_1,sequence_2,line;
    ifstream in_file("g1.fasta");
    if (in_file.is_open())
    {
        getline(in_file, line);
        int count=0;

        while (getline(in_file, line)  )
        {
            sequence_1 += line;
        }
        in_file.close();
    }
    ifstream in_file_2("g2.fasta");
    if (in_file_2.good())
    {
        getline(in_file_2, line);
        int count=0;
        while (getline(in_file_2, line) )
        {
            sequence_2 += line;
        }

        in_file_2.close();
    }
    cout<<"X drop algorithm"<<endl;
    auto start =steady_clock::now();
    int N = sequence_1.size(); 
    int M = sequence_2.size(); 
    double** S = (double**)malloc((2*N+1) * sizeof(double*));
    for (int i = 0; i < (2*M+1); i++) {
        S[i] = (double*)malloc((2*M+1) * sizeof(double));     
    }
    S[0][0] = T;

    while (L <= U + 2) {
        k+=2;
        for (int i = L; i <= U + 2; i += 1) {
            j = k - i;
            if (i%2) {
                
                S[i][j] = (sequence_1[i/2] == sequence_2[j/2]) ? (S[i-1][j-1] + mat/2) : (S[i-1][j-1] + mis/2);
                
            } else if(i>0 && j>0){
                
                double score1 = (i >= L && i - 2 <= U && sequence_1[(i/2)-1] == sequence_2[(j/2)-1]) ? (S[i-1][j-1] + mat/2) : INF;
                double score2 = (i >= L && i - 2 <= U && sequence_1[(i/2)-1] != sequence_2[(j/2)-1]) ? (S[i-1][j-1] + mis/2) : INF;
                double score3 = (i <= U) ? S[i][j-2] + ind : INF;
                double score4 = (i - 2 >= L) ? S[i-2][j] + ind : INF;
                S[i][j] = max(max(score1, score2), max(score3, score4));

            }
            else{
                double score3 = (i <= U) ? S[i][j-2] + ind : INF;
                double score4 = (i - 2 >= L) ? S[i-2][j] + ind : INF;
                S[i][j] = max(score3, score4);

            }
            T_prime = max(T_prime, S[i][j]);
            if (S[i][j] < T - X) {
                S[i][j] = INF;
            }
        }
        L = -INF;
        U = INF;
        for (int i = 0; i <= k; i++) {
            if (S[i][k-i] > INF) {
                L = min(L, i);
                U = max(U, i);
            }
        }
        L = max(L, 2+k-2*N);
        U = min(U, 2*M-2);
        T = T_prime;
    }
    free(S);
    auto stop =steady_clock::now();
    auto totime=duration_cast<milliseconds>(stop - start).count();
    cout<<"Optimal score:"<<T_prime << endl;
    cout<<"Time taken in millisec:"<<totime<<endl;
    return 0;
}
