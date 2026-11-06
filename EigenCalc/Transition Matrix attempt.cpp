#include <iostream>
#include <vector>
#include <iterator>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "pbPlots.hpp"
#include "supportLib.hpp"


using namespace std;
using namespace Eigen;

void InitTransMatrix(Matrix<double, 6, 6> (&arr), double Ps, double Pe, double Pi, double Pr, double Pv)
{
    arr(0, 0) = Ps;
    arr(1, 0) = 1 - Ps;
    for(int n = 0; n < 2; n++)
    {
        arr(n+2, 0) = 0;
    }

    arr(0, 1) = Pe;
    arr(1, 1) = 0;
    for(int n = 0; n < 2; n++)
    {
        arr(n+2, 1) = (1-Pe)/2;
    }

    for(int n = 0; n < 3; n++)
    {
        arr(n, 2) = 0;
    }
    arr(3, 2) = Pi;

    arr(0, 3) = 1 - Pr;
    for(int n = 0; n < 2; n++)
    {
        arr(n+1, 3) = 0;
    }
    arr(3, 3) = Pr;

    for(int n = 0; n < 4; n++)
    {
        arr(4, n) = 0;
    }

    arr(5, 3) = 0; 
    arr(0, 3) = 1 - Pr;

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

void IterateStepSEIR(Matrix<double, 6, 6> (&arr1), double (&arr2)[6], int iter, vector<double> (&x), vector<double> (&y1), vector<double> (&y2), vector<double> (&y3), vector<double> (&y4), vector<double> (&y5), vector<double> (&y6))
{
    if(iter == 1)
    {
    x[0] = 0;
    y1[0] = arr2[0];
    y2[0] = arr2[1];
    y3[0] = arr2[2];
    y4[0] = arr2[3];
    y5[0] = arr2[4];
    y6[0] = arr2[5];
    }
    
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

    if (normal_sum > 0)
    {
        for(int i = 0; i < 4; i++)
        {
            arr2[i] /= normal_sum;
        }
    }

    arr2[4] = 0; 
    arr2[5] = 0;

    x[iter] = iter;
    y1[iter] = arr2[0];
    y2[iter] = arr2[1];
    y3[iter] = arr2[2];
    y4[iter] = arr2[3];
    y5[iter] = arr2[4];
    y6[iter] = arr2[5];
}

void IterateStepSEIRIm(Matrix <double, 6, 6> (&arr1), double (&arr2)[6], int iter, vector<double> (&x), vector<double> (&y1), vector<double> (&y2), vector<double> (&y3), vector<double> (&y4), vector<double> (&y5), vector<double> (&y6))
{
    if(iter == 1)
    {
    x[0] = 0;
    y1[0] = arr2[0];
    y2[0] = arr2[1];
    y3[0] = arr2[2];
    y4[0] = arr2[3];
    y5[0] = arr2[4];
    y6[0] = arr2[5];
    }

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

    if (normal_sum > 0)
    {
        for(int i = 0; i < 6; i++)
        {
            if(i != 4)
                arr2[i] /= normal_sum;
        }
    }

    x[iter] = iter;
    y1[iter] = arr2[0];
    y2[iter] = arr2[1];
    y3[iter] = arr2[2];
    y4[iter] = arr2[3];
    y5[iter] = arr2[4];
    y6[iter] = arr2[5];
}

void IterateStepSEIRVIm(Matrix <double, 6, 6> (&arr1), double (&arr2)[6], int iter, vector<double> (&x), vector<double> (&y1), vector<double> (&y2), vector<double> (&y3), vector<double> (&y4), vector<double> (&y5), vector<double> (&y6))
{

    if(iter == 1)
    {
    x[0] = 0;
    y1[0] = arr2[0];
    y2[0] = arr2[1];
    y3[0] = arr2[2];
    y4[0] = arr2[3];
    y5[0] = arr2[4];
    y6[0] = arr2[5];
    }

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

    if (normal_sum > 0)
    {
        for(int i = 0; i < 6; i++)
        {
            arr2[i] /= normal_sum;
        }
    }

    x[iter] = iter;
    y1[iter] = arr2[0];
    y2[iter] = arr2[1];
    y3[iter] = arr2[2];
    y4[iter] = arr2[3];
    y5[iter] = arr2[4];
    y6[iter] = arr2[5];
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
    int remove4[5] = {0, 1, 2, 3, 5};

    for(int i = 0; i < 5; i++)
    {
        for(int j = 0; j < 5; j++)
        {
            arr1(i, j) = ref_mat(remove4[i], remove4[j]);
        }
    }
}

void InitRef_matSEIRVIm(Matrix <double, 6, 6> (&arr1), MatrixXd (&ref_mat))
{
    ref_mat.resize(6, 6);
    for(int i = 0; i < 6; i++)
    {
        for(int j = 0; j < 6; j++)
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
    int iterations = 0;

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
    cin >> InitialVec[0]; RefVec[0] = InitialVec[0];
    cout << "What is the Initial E value?" << endl;
    cin >> InitialVec[1]; RefVec[1] = InitialVec[1];
    cout << "What is the Initial I value?" << endl;
    cin >> InitialVec[2]; RefVec[2] = InitialVec[2];
    cout << "What is the Initial R value?" << endl;
    cin >> InitialVec[3]; RefVec[3] = InitialVec[3];
    cout << "What is the Initial V value?" << endl;
    cin >> InitialVec[4]; RefVec[4] = InitialVec[4];
    cout << "What is the Initial Im value?" << endl;
    cin >> InitialVec[5]; RefVec[5] = InitialVec[5];

    int exitcon = 1;

    MatrixXcd eigenvalueMatrix;
    MatrixXcd eigenvectorMatrix;
    MatrixXd ref_mat(6, 6);

    vector<double> xall;
    vector<double> y1, y2, y3, y4, y5, y6;

    ScatterPlotSeries *susc = GetDefaultScatterPlotSeriesSettings();
    susc->xs = &xall; susc->ys = &y1;
    //susc->linearInterpolation = true; susc->lineThickness = 2;
    susc->linearInterpolation = false;
    susc->color = CreateRGBColor(1.0, 0.0, 0.0);
    susc->pointType = toVector(L"dots");

    ScatterPlotSeries *expo = GetDefaultScatterPlotSeriesSettings();
    expo->xs = &xall; expo->ys = &y2;
    //expo->linearInterpolation = true; expo->lineThickness = 2;
    expo->linearInterpolation = false;
    expo->color = CreateRGBColor(0.5, 0.5, 0.0);
    expo->pointType = toVector(L"dots");

    ScatterPlotSeries *infe = GetDefaultScatterPlotSeriesSettings();
    infe->xs = &xall; infe->ys = &y3;
    //infe->linearInterpolation = true; infe->lineThickness = 2;
    infe->linearInterpolation = false;
    infe->color = CreateRGBColor(0.0, 1.0, 0.0);
    infe->pointType = toVector(L"dots");

    ScatterPlotSeries *reco = GetDefaultScatterPlotSeriesSettings();
    reco->xs = &xall; reco->ys = &y4;
    //reco->linearInterpolation = true; reco->lineThickness = 2;
    reco->linearInterpolation = false;
    reco->color = CreateRGBColor(0.0, 0.75, 0.25);
    reco->pointType = toVector(L"dots");

    ScatterPlotSeries *vacc = GetDefaultScatterPlotSeriesSettings();
    vacc->xs = &xall; vacc->ys = &y5;
    //vacc->linearInterpolation = true; vacc->lineThickness = 2;
    vacc->linearInterpolation = false;
    vacc->color = CreateRGBColor(0.0, 0.25, 0.75);
    vacc->pointType = toVector(L"dots");

    ScatterPlotSeries *immu = GetDefaultScatterPlotSeriesSettings();
    immu->xs = &xall; immu->ys = &y6;
    //immu->linearInterpolation = true; immu->lineThickness = 2;
    immu->linearInterpolation = false;
    immu->color = CreateRGBColor(0.0, 0.0, 1.0);
    immu->pointType = toVector(L"dots");

    ScatterPlotSettings *settings = GetDefaultScatterPlotSettings();
    settings->width = 600;
    settings->height = 400;
    settings->autoBoundaries = false;
    settings->xMin = 0;
    settings->yMin = 0;
    settings->yMax = 1;
    settings->autoPadding = true;
    settings->xLabel = toVector(L"Iterations");
    settings->yLabel = toVector(L"Probability");

    vector<ScatterPlotSeries*> GraphVec = {susc, expo, infe, reco, vacc, immu};

    while(exitcon != 0)
    {
        InitRef(transmatrix, reftransmatrix);

        for(int i = 0; i < 6; i++)
        {
            RefVec[i] = InitialVec[i];
        }

        cout << "What would you like to do?\n0: End Program\n1: Iterate SEIR model\n2: Iterate SEIRIm model\n3: Iterate SEIRVIm model\n4: Find Eigenvalues and Eigenvectors for Probability Matrix\n";
        cin >> exitcon;

        switch(exitcon)
        {
            case 0:
                break;

            case 1:

                xall.clear();
                y1.clear();
                y2.clear();
                y3.clear();
                y4.clear();
                y5.clear();
                y6.clear();

                cout << "How many iterations do you want to calculate?" << endl;
                cin >> iterations;

                settings->title = toVector(L"SEIR Model Simulation");
                settings->xMin = 0;
                settings->xMax = iterations;

                xall.resize(iterations + 1);
                y1.resize(iterations + 1);
                y2.resize(iterations + 1);
                y3.resize(iterations + 1);
                y4.resize(iterations + 1);
                y5.resize(iterations + 1);
                y6.resize(iterations + 1);

                for(int i = 0; i < iterations; i++)
                {
                    IterateStepSEIR(reftransmatrix, RefVec, i, xall, y1, y2, y3, y4, y5, y6);
                }

                {
                    settings->scatterPlotSeries = &GraphVec;
                    RGBABitmapImageReference *imageReference = CreateRGBABitmapImageReference();
                    StringReference *errorMessage = CreateStringReferenceLengthValue(0, L' ');
                    
                    bool success = DrawScatterPlotFromSettings(imageReference, settings, errorMessage);
                    if (success) 
                    {
                        vector<double> *pngdata = ConvertToPNG(imageReference->image);
                        WriteToFile(pngdata, "seir_plot.png");
                        cout << "Plot saved to seir_plot.png\n";
                    }
                }

                cout << "Your vector is:\n";
                for(int i = 0; i < 4; i++)
                {
                    cout << RefVec[i] << endl;
                }
                cout << endl;
                break;

            case 2:

                xall.clear();
                y1.clear();
                y2.clear();
                y3.clear();
                y4.clear();
                y5.clear();
                y6.clear();

                cout << "How many iterations do you want to calculate?" << endl;
                cin >> iterations;

                settings->title = toVector(L"SEIRIm Model Simulation");
                settings->xMin = 0;
                settings->xMax = iterations;
                
                reftransmatrix(5, 3) = reftransmatrix(0, 3) / 2;
                reftransmatrix(0, 3) *= 0.5;

                xall.resize(iterations + 1);
                y1.resize(iterations + 1);
                y2.resize(iterations + 1);
                y3.resize(iterations + 1);
                y4.resize(iterations + 1);
                y5.resize(iterations + 1);
                y6.resize(iterations + 1);

                for(int i = 0; i < iterations; i++)
                {
                    IterateStepSEIRIm(reftransmatrix, RefVec, i, xall, y1, y2, y3, y4, y5, y6);
                }

                {
                    settings->scatterPlotSeries = &GraphVec;
                    RGBABitmapImageReference *imageReference = CreateRGBABitmapImageReference();
                    StringReference *errorMessage = CreateStringReferenceLengthValue(0, L' ');
                    
                    bool success = DrawScatterPlotFromSettings(imageReference, settings, errorMessage);
                    if (success) 
                    {
                        vector<double> *pngdata = ConvertToPNG(imageReference->image);
                        WriteToFile(pngdata, "seirim_plot.png");
                        cout << "Plot saved to seirim_plot.png\n";
                    }
                }

                for(int i = 0; i < 6; i++)
                {
                    cout << RefVec[i] << endl;
                }
                cout << endl;
                break;

            case 3:

                xall.clear();
                y1.clear();
                y2.clear();
                y3.clear();
                y4.clear();
                y5.clear();
                y6.clear();

                cout << "How many iterations do you want to calculate?" << endl;
                cin >> iterations;

                settings->title = toVector(L"SEIRVIm Model Simulation");
                settings->xMin = 0;
                settings->xMax = iterations;

                xall.resize(iterations + 1);
                y1.resize(iterations + 1);
                y2.resize(iterations + 1);
                y3.resize(iterations + 1);
                y4.resize(iterations + 1);
                y5.resize(iterations + 1);
                y6.resize(iterations + 1);

                for(int i = 0; i < iterations; i++)
                {
                    IterateStepSEIRVIm(reftransmatrix, RefVec, i, xall, y1, y2, y3, y4, y5, y6);
                }

                {
                    settings->scatterPlotSeries = &GraphVec;
                    RGBABitmapImageReference *imageReference = CreateRGBABitmapImageReference();
                    StringReference *errorMessage = CreateStringReferenceLengthValue(0, L' ');
                    
                    bool success = DrawScatterPlotFromSettings(imageReference, settings, errorMessage);
                    if (success) 
                    {
                        vector<double> *pngdata = ConvertToPNG(imageReference->image);
                        WriteToFile(pngdata, "seirvim_plot.png");
                        cout << "Plot saved to seirvim_plot.png\n";
                    }
                }

                for(int i = 0; i < 6; i++)
                {
                    cout << RefVec[i] << endl;
                }
                cout << endl;
                break;

            case 4:
            {
                int l;
                cout << "Which model are you finding the eigenvectors and/or eigenvalues for?\n1. SEIR\n2. SEIRIm\n3. SEIRVIm\n";
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
                        cout << "\nComputed EigenVectors:\n\n" << eigenvectorMatrix << endl;
                        cout << endl;
                        break;
                    }
                    case 2:
                    {
                        reftransmatrix(5, 3) = reftransmatrix(0, 3) / 2;
                        reftransmatrix(0, 3) *= 0.5;
                        InitRef_matSEIRIm(reftransmatrix, ref_mat);
                        EigenSolver<MatrixXd> eigenValueSolver(ref_mat);
                        eigenvalueMatrix = eigenValueSolver.eigenvalues();
                        eigenvectorMatrix = eigenValueSolver.eigenvectors();

                        cout << "Computed Eigenvalues:\n" << eigenvalueMatrix << endl;
                        cout << "\nComputed EigenVectors:\n" << eigenvectorMatrix << endl;
                        cout << endl;
                        break;
                    }
                    case 3:
                    {
                        InitRef_matSEIRVIm(reftransmatrix, ref_mat);
                        EigenSolver<MatrixXd> eigenValueSolver(ref_mat);
                        eigenvalueMatrix = eigenValueSolver.eigenvalues();
                        eigenvectorMatrix = eigenValueSolver.eigenvectors();

                        cout << "Computed Eigenvalues:\n" << eigenvalueMatrix << endl;
                        cout << "\nComputed EigenVectors:\n" << eigenvectorMatrix << endl;
                        cout << endl;
                        break;
                    }
                }
                break;
            }
        }
    }

    return 0;
}

