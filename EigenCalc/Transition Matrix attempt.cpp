#include <iostream>
#include <Eigen/Dense>

#include <vector>
#include <iterator>
#include <Eigen/Eigenvalues>

using namespace std;
using namespace Eigen;


void InitTransMatrix(Matrix <double, 6, 6> (&arr), double Ps, double Pe, double Pi, double Pr, double Pv)
    {
        arr(0, 0) = Ps;
        arr(1, 0) = 1 - Ps;
        for(int n = 0;n < 2; n++)
        {
            arr(n+2, 0) = 0;
        }

        arr(0, 1) = Pe;
        arr(1, 1) = 0;
        for(int n = 0;n < 2; n++)
        {
            arr(n+2, 1) = (1-Pe)/2;
        }

        for(int n = 0;n < 3; n++)
        {
            arr(n, 2) = 0;
        }
        arr(3, 2) = Pi;

        arr(0, 3) = 1 - Pr;
        for(int n = 0;n < 2; n++)
        {
            arr(n+1, 3) = 0;
        }
        arr(3, 3) = Pr;

        for(int n = 0; n < 4; n++)
        {
            arr(4, n) = 0;
        }

        //change to make it so arr(4, 3) is set equal to 1 - (Pr/2) along with arr(0, 3))
        
        arr(5, 3) = 0, arr(0, 3) = 1 - Pr;

        for(int n = 0; n < 4; n++)
        {
            arr(n, 4) = 0;
        }

        arr(4, 4) = Pv;


        for(int n = 0; n < 4; n++)
        {
            arr(5, n) = 0;
        }

        arr(5, 4) = 1 - Pv;

        for(int n = 0; n < 5; n++)
        {
            arr(n, 5) = 0;
        }

        arr(5, 5) = 1;        

    }

void IterateStepSEIR(Matrix <double, 6, 6> (&arr1), double (&arr2)[6])
    {

        double Vec[4] = {arr2[0], arr2[1], arr2[2], arr2[3]};

        for(int i = 0; i < 4; i++)
        {

            double sum = 0.0;

            for(int j = 0; j < 4; j++)
            {
                sum += arr1(i, j) * Vec[j];

            }

            arr2[i] = sum;

        }

        double normal_sum = arr2[0] + arr2[1] + arr2[2] + arr2[3];

        for(int i = 0; i < 4; i++)
        {
            arr2[i] /= normal_sum;
        }

        arr2[4] = 0, arr2[5] = 0;

    }

void IterateStepSEIRIm(Matrix <double, 6, 6> (&arr1), double (&arr2)[6])
    {

        double Vec[6] = {arr2[0], arr2[1], arr2[2], arr2[3], 0, arr2[5]};

        for(int i = 0; i < 6; i++)
        {

            double sum = 0.0;
            if(i != 4)
            {
                for(int j = 0; j < 6; j++)
                {
                    if(j != 4)
                    {
                        sum += arr1(i, j) * Vec[j];
                    }

                }
            }

            arr2[i] = sum;

        }

        double normal_sum = arr2[0] + arr2[1] + arr2[2] + arr2[3] + arr2[5];

        for(int i = 0; i < 6; i++)
        {
            if(i != 4)
            arr2[i] /= normal_sum;
        }

    }

void IterateStepSEIRVIm(Matrix <double, 6, 6> (&arr1), double (&arr2)[6])
    {

        double Vec[6] = {arr2[0], arr2[1], arr2[2], arr2[3], arr2[4], arr2[5]};

        for(int i = 0; i < 6; i++)
        {

            double sum = 0.0;

            for(int j = 0; j < 6; j++)
            {
                sum += arr1(i, j) * Vec[j];

            }

            arr2[i] = sum;

        }

        double normal_sum = arr2[0] + arr2[1] + arr2[2] + arr2[3] + arr2[4] + arr2[5];

        for(int i = 0; i < 6; i++)
        {
            arr2[i] /= normal_sum;
        }

    }

void InitRef(Matrix <double, 6, 6> (&arr1), Matrix <double, 6, 6> (&arr3))
    {
        for(int i = 0; i < 6; i++)
    {
        for(int j = 0; j < 6; j++)
        {
            arr3(i, j) = arr1(i, j);
        }
    }
    }

void InitRef_matSEIR(Matrix <double, 6, 6> (&arr1), MatrixXd (&ref_mat))
    {
        ref_mat.resize(4, 4);

        for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            ref_mat(i, j) = arr1(i, j);
        }
    }
    }

void InitRef_matSEIRIm(Matrix <double, 6, 6> (&arr1), MatrixXd (&ref_mat))
    {

        ref_mat = arr1;
        ref_mat.resize(5, 5);

        int l = 0;

        for(int i = 0; i < 6; i++)
        {

            if(i != 4)
            {

                int k = 0;

                for(int j = 0; j < 6; j++)
                {
                    if(j != 4)
                    {
                        arr1(i, j) = ref_mat(i - l, j - k);
                    }
                    else
                    {
                        k += 1;
                    }

                }
            }
            else
            {
                l += 1;
            }

        }
    
    }


void InitRef_matSEIRVIm(Matrix <double, 6, 6> (&arr1), MatrixXd (&ref_mat))
    {

        ref_mat.resize(6, 6);

        for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            ref_mat(i, j) = arr1(i, j);
        }
    }
    }


int main(){

    Matrix <double, 6, 6> transmatrix;
    Matrix <double, 6, 6> reftransmatrix;

    double InitialVec[6];
    double RefVec[6];

    double Ps, Pe, Pi, Pr, Pv;

    int iterations;

    
    cout << "What is Ps?" << endl;
    cin >> Ps;
    cout << "What is Pe?" << endl;
    cin >> Pe;
    cout << "What is Pi?" << endl;
    cin >> Pi;
    cout << "What is Pr?" << endl;
    cin >> Pr;
    cout << "What is Pv?" << endl;
    cin >> Pv;

    InitTransMatrix(transmatrix, Ps, Pe, Pi, Pr, Pv);
    InitRef(transmatrix, reftransmatrix);

    cout << "What is the Initial S value?" << endl;
    cin >> InitialVec[0];
    RefVec[0] = InitialVec[0];
    cout << "What is the Initial E value?" << endl;
    cin >> InitialVec[1];
    RefVec[1] = InitialVec[1];
    cout << "What is the Initial I value?" << endl;
    cin >> InitialVec[2];
    RefVec[2] = InitialVec[2];
    cout << "What is the Initial R value?" << endl;
    cin >> InitialVec[3];
    RefVec[3] = InitialVec[3];
    cout << "What is the Initial V value?" << endl;
    cin >> InitialVec[4];
    RefVec[4] = InitialVec[4];
    cout << "What is the Initial Im value?" << endl;
    cin >> InitialVec[5];
    RefVec[5] = InitialVec[5];

    int exitcon = 1;

    MatrixXcd eigenvalueMatrix;
    MatrixXcd eigenvectorMatrix;
    MatrixXd ref_mat(6, 6);

    while(exitcon != 0)
    {

        InitRef(transmatrix, reftransmatrix);

        for(int i = 0; i < 6; i++)
        {
            RefVec[i] = InitialVec[i];
        }

        cout << "What would you like to do?\n0: End Program\n 1: Iterate SEIR model\n 2: Iterate SIERIm model\n 3: Iterate SEIRVIm model\n 4: Find Eigenvalues and Eigenvectors for Porbability Matrix\n";
        cin >> exitcon;

        switch(exitcon)
        {
            

            case 0:
            break;

            case 1:

                cout << "How many iterations do you want to calculate?" << endl;
                cin >> iterations;

                for(int i = 0; i < iterations; i++)
                {
                    IterateStepSEIR(reftransmatrix, RefVec);
                }

                cout << "Your vector is " << endl;

                for(int i = 0; i < 4; i++)
                {
                    cout << RefVec[i] << endl;
                }
                cout << endl;
                cin >> exitcon;

            break;

            case 2:

                cout << "How many iterations do you want to calculate?" << endl;
                cin >> iterations;
                
                reftransmatrix(5, 3) = reftransmatrix(0, 3) / 2;
                reftransmatrix(0, 3) *= 0.5;

                for(int i = 0; i < iterations; i++)
                {
                    IterateStepSEIRIm(reftransmatrix, RefVec);
                }

                for(int i = 0; i < 6; i++)
                {
                    cout << RefVec[i] << endl;
                }
                cout << endl;
                cin >> exitcon;

            break;

            case 3:
            
                cout << "How many iterations do you want to calculate?" << endl;
                cin >> iterations;

                for(int i = 0; i < iterations; i++)
                {
                    IterateStepSEIRVIm(reftransmatrix, RefVec);
                }

                for(int i = 0; i < 6; i++)
                {
                    cout << RefVec[i] << endl;
                }
                cout << endl;
                cin >> exitcon;

            break;

            case 4:

            int l;
            cout << "Which model are you finding the eigenvectors and/or eigenvalues for?\n 1. SEIR\n 2. SEIRIm\n 3. SEIRVIm\n";
            cin >> l;
            switch(l)
            {
                case 1:
                {

                    InitRef_matSEIR(reftransmatrix, ref_mat);

                    EigenSolver<MatrixXd> eigenValueSolver(ref_mat);

                    eigenvalueMatrix = eigenValueSolver.eigenvalues();
                    eigenvectorMatrix = eigenValueSolver.eigenvectors();

                    cout << "\nComputed Eigenvalues:\n\n" << eigenvalueMatrix << endl;
                    cout << "\nComputed EigenVectors\n\n" << eigenvectorMatrix << endl;
                    cout << endl;

                    cin >> exitcon;
                }
                break;

                case 2:
                {

                    reftransmatrix(5, 3) = reftransmatrix(0, 3) / 2;
                    reftransmatrix(0, 3) *= 0.5;

                    InitRef_matSEIRIm(reftransmatrix, ref_mat);

                    EigenSolver<MatrixXd> eigenValueSolver(ref_mat);

                    eigenvalueMatrix = eigenValueSolver.eigenvalues();
                    eigenvectorMatrix = eigenValueSolver.eigenvectors();

                    cout << "Computed Eigenvalues:\n" << eigenvalueMatrix << endl;
                    cout << "\nComputed EigenVectors\n" << eigenvectorMatrix << endl;
                    cout << endl;

                    cin >> exitcon;
                }

                break;

                case 3:

                {

                    InitRef_matSEIRVIm(reftransmatrix, ref_mat);

                    EigenSolver<MatrixXd> eigenValueSolver(ref_mat);

                    eigenvalueMatrix = eigenValueSolver.eigenvalues();
                    eigenvectorMatrix = eigenValueSolver.eigenvectors();

                    cout << "Computed Eigenvalues:\n" << eigenvalueMatrix << endl;
                    cout << "\nComputed EigenVectors\n" << eigenvectorMatrix << endl;
                    cout << endl;

                    cin >> exitcon;
                }
                
                break;

                // normalize eigen vectors 
            }

            

        }
    }


    return 0;
}