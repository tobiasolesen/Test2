//Bruk &k nar du vil endre variabelen k inni funksjonen. Brukes i prototypen av funksjonen og i def av funksjonen.
//Ikke i kallet pa funksjonen. Trenger da ikke bruke *k noe sted tror jeg.
//Matriser skal vaere dobbeltpekere.
//Vektorer skal vaere enkeltpekere.
//Vet at de 3 laveste egenverdiene skal vaere: 3, 7, 11

//i d) skal jeg bruke koden fra b), endrer bare pa potensialet V(rho)
//rho er avstanden mellom partiklene (?)
//Vil finne den minste egenverdien (tilsvarer energien til ground state) og tilhorende egenvektor (bolgefunksjonen i ground state)
//apne tekstfil i terminal: less "filnavn"
//Komme ut av fil og tilbake i terminal: shift + q
//Med w_r=0.25 skal de minste egenverdiene bli: 1.25, 2.19, 3.18 ca

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include <iostream>
#include <cmath>
#include <array>
#include <vector>
#include <fstream>
//#include "jacobi.h"

using namespace std;

//Function prototypes:
void printMatrixA(double** A, int N);
double maxoffdiag(double** A, int& k, int& l, int N);
vector<int> jacobi_method(double** A, double** R, int N);  //returneren en vektor med doubles
void rotate ( double ** A, double ** R, int k, int l, int N );

int main(){

    //Definerer var som trengs for a lage matrisen A
    //N = 200 er passe?
    int N = 200;
    double rho_0 = 0;
    double rho_max = 8;
    double h = (rho_max - rho_0)/((double)N);
    double hh = h*h;
    double e = -1/hh;
    //double w_r = 0.25;
    double w_r = 5;
    double w_r2 = w_r*w_r;

    //Lager vektorer som trengs:
    double* rho = new double[N];
    double* V = new double[N];
    double* d = new double[N];
    double* eigenvec0 = new double[N];

    //Fyller vektorene
    for (int i = 0; i < N; i++){
        rho[i] = rho_0 + (i+1)*h;
        V[i] = w_r2*rho[i]*rho[i] + 1.0/rho[i];
        d[i] = 2./hh + V[i];
    }

    //Deklarerer/lager matrisene A og R:
    double** A = new double*[N];
    double** R = new double*[N];
    for (int i = 0; i < N; i++){
        A[i] = new double[N];
        R[i] = new double[N];
    }

    //Fyller inn offdiag-elementene i A:
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (i == j-1){
                A[i][j] = e;
            } else if (i == j+1){
                A[i][j] = e;
            } else {
                A[i][j] = 0.0;
            }
        }
    }
    //Fyller diagonalelementene i A:
    for (int i = 0; i < N; i++){
        A[i][i] = d[i];
    }
    //printMatrixA(A, N);  //Skriver ut matrisen A som vi starter med

    vector <int> ind = jacobi_method(A, R, N); //Kaller jacobi_method. ind inneholder

    //Skriver til filen data.dat
    ofstream outFile;
    outFile.open("egenverdier.dat", std::ios::out);

    for (int j = 0; j < N; j++){
        eigenvec0[j] = R[ind[0]][j];  //eigenvec0 blir egenvektoren til grunntilstanden (dvs til laveste egenverdi)
        outFile << eigenvec0[j] << endl;  //Skriver egenvektoren til fil
    }
    outFile.close();

    return 0;
}  //Slutt pa main

//Funksjon som skriver ut matrise
void printMatrixA(double** A, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%6.2f ", A[i][j]);
        }
        cout << endl;
    }
    cout << endl;
    fflush(stdout);
}

//Setting up the eigenvector matrix (som identitetsmatrisen?):
vector<int> jacobi_method(double** A, double** R, int N) {

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (i == j){
                R[i][j] = 1.0;
            }
            else{
                R[i][j] = 0.0;
            }
        }
    }

    int k, l;
    double epsilon = 1.0e-8;
    double max_number_iterations = (double) N * (double) N * (double) N;
    int iterations = 0;
    double max_offdiag = maxoffdiag(A, k, l, N);

    while ( fabs(max_offdiag) > epsilon && (double) iterations < max_number_iterations ){
        max_offdiag = maxoffdiag ( A, k, l, N );
        rotate ( A, R, k, l, N );
        iterations++;
        //Hopper ut av denne lokken nar offdiag elementene til A er essentially zero
    }

    cout << "Number of iterations: " << iterations << "\n";

    vector<double>diag(N);  //Deklarerer vektor diag med lengde N?

    for (int i = 0; i < N; i++){
        diag[i] = A[i][i];  //Fyller diag med A sine diagonalelementer
    }

    //finner indeks til minste egenverdiene og lagrer indeksen i lowestIndex
    int ll = 3;  //Antall egenverdier jeg bryr meg om
    vector<double> lowest(ll);
    vector<int> lowestIndex(ll);  //vektor som skal inneholde indeksene til de minste egenverdiene
    for (int i = 0; i < ll; i++) {
        int index = 0;
        double low = 200;

        for (unsigned int j = 0; j < diag.size(); j++) {
            if (diag[j] < low) {
                low = diag[j];
                index = j;
            }
        }
        lowest[i] = low;
        lowestIndex[i] = index;
        diag[index] = 91385938659368;
    }

    //printMatrixA(A,N);
    for (int i = 0; i < ll; i++) {
        cout <<  "Laveste egenverdier: " << lowest[i] << " tilhorende index: " << lowestIndex[i] << endl;
    }

    sort(diag.begin(), diag.end());  //Sorterer vektoren fra minst til storste verdier
    //Skriver ut diag sine 3 forste elementer (3 forste egenverdiene):
    //cout << "3 forste egenverdier: " << diag[0] << " " << diag[1] << " " << diag[2] << endl;

    return lowestIndex;
}

//Function to find the maximum element of A:
double maxoffdiag ( double** A, int& k, int& l, int N) {

    double max = 0.0;
    for ( int i = 0; i < N; i++ ) {
        for ( int j = i + 1; j < N; j++ ) {
            if ( fabs(A[i][j]) > max ){
                max = fabs(A[i][j]);
                l = i;
                k = j;  //Maks-elementet girs indeksene [k][l], og [l][k] senere?
            }
        }
    }
    return max;
}

// Function to find the values of cos and sin
void rotate ( double ** A, double ** R, int k, int l, int N ){

    double s, c;
    //cout << k << " " << l << endl;

    if ( A[k][l] != 0.0 ) {

        double t, tau;
        tau = (A[l][l] - A[k][k])/(2*A[k][l]);

        if ( tau > 0 ) {
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        }
        else {
            t = -1.0/( -tau + sqrt(1.0 + tau*tau));
        }
        c = 1.0/sqrt(1.0+t*t);
        s = c*t;
    } else{
        c = 1.0;
        s = 0.0;
    }

    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A[k][k];
    a_ll = A[l][l];

    // changing the matrix elements with indices k and l A[k][k] = c*c*a_kk - 2.0*c*s*A[k][l] + s*s*a_ll; A[l][l] = s*s*a_kk + 2.0*c*s*A[k][l] + c*c*a_ll; A[k][l] = 0.0; // hard-coding of the zeros

    A[k][k] = a_kk*c*c - 2*A[k][l]*c*s + a_ll*s*s;
    A[l][l] = a_ll*c*c + 2*A[k][l]*c*s + a_kk*s*s;
    A[l][k] = 0.0;
    A[k][l] = 0.0;
    // and then we change the remaining elements
    for ( int i = 0; i < N; i++ ) {
        if ( i != k && i != l ){
            a_ik = A[i][k];
            a_il = A[i][l];
            A[i][k] = c*a_ik - s*a_il;
            A[k][i] = A[i][k];
            A[i][l] = c*a_il + s*a_ik;
            A[l][i] = A[i][l];
        }
        // Finally, we compute the new eigenvectors
        r_ik = R[i][k];
        r_il = R[i][l];
        R[i][k] = c*r_ik - s*r_il;
        R[i][l] = c*r_il + s*r_ik;
    }

    return;

}
