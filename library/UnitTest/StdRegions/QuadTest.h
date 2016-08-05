//
// Created by li12242 on 16/6/30.
//

#ifndef DGOM_QUADTEST_H
#define DGOM_QUADTEST_H

#include "StdRegions/StdRegions.h"

/* set parameters for standard quadrilateral element */
#define Deg 2                   /* degree */
#define Dfp Deg+1               /* number of points at each face */
#define Dnp (Deg+1)*(Deg+1)     /* number of points */
#define Dnfaces 4               /* number of faces */

double TestQuad_r[Dnp] = {-1,  0,  1, -1, 0, 1, -1, 0, 1};
double TestQuad_s[Dnp] = {-1, -1, -1,  0, 0, 0,  1, 1, 1};

int    TestQuad_fmask[Dnfaces][Dfp] = {{0,1,2},{2,5,8},{8,7,6},{6,3,0}};
double TestQuad_V[Dnp][Dnp] =
        {{5.0000000000000011e-01,-8.6602540378443860e-01,1.1180339887498949e+00,-8.6602540378443860e-01,1.4999999999999998e+00,-1.9364916731037083e+00,1.1180339887498949e+00,-1.9364916731037083e+00,2.5000000000000004e+00,},
        {5.0000000000000011e-01,-8.6602540378443860e-01,1.1180339887498949e+00,0.0000000000000000e+00,-0.0000000000000000e+00,0.0000000000000000e+00,-5.5901699437494745e-01,9.6824583655185414e-01,-1.2500000000000002e+00,},
        {5.0000000000000011e-01,-8.6602540378443860e-01,1.1180339887498949e+00,8.6602540378443860e-01,-1.4999999999999998e+00,1.9364916731037083e+00,1.1180339887498949e+00,-1.9364916731037083e+00,2.5000000000000004e+00,},
        {5.0000000000000011e-01,0.0000000000000000e+00,-5.5901699437494745e-01,-8.6602540378443860e-01,-0.0000000000000000e+00,9.6824583655185414e-01,1.1180339887498949e+00,0.0000000000000000e+00,-1.2500000000000002e+00,},
        {5.0000000000000011e-01,0.0000000000000000e+00,-5.5901699437494745e-01,0.0000000000000000e+00,0.0000000000000000e+00,-0.0000000000000000e+00,-5.5901699437494745e-01,-0.0000000000000000e+00,6.2500000000000011e-01,},
        {5.0000000000000011e-01,0.0000000000000000e+00,-5.5901699437494745e-01,8.6602540378443860e-01,0.0000000000000000e+00,-9.6824583655185414e-01,1.1180339887498949e+00,0.0000000000000000e+00,-1.2500000000000002e+00,},
        {5.0000000000000011e-01,8.6602540378443860e-01,1.1180339887498949e+00,-8.6602540378443860e-01,-1.4999999999999998e+00,-1.9364916731037083e+00,1.1180339887498949e+00,1.9364916731037083e+00,2.5000000000000004e+00,},
        {5.0000000000000011e-01,8.6602540378443860e-01,1.1180339887498949e+00,0.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,-5.5901699437494745e-01,-9.6824583655185414e-01,-1.2500000000000002e+00,},
        {5.0000000000000011e-01,8.6602540378443860e-01,1.1180339887498949e+00,8.6602540378443860e-01,1.4999999999999998e+00,1.9364916731037083e+00,1.1180339887498949e+00,1.9364916731037083e+00,2.5000000000000004e+00,}};

double TestQuad_M[Dnp][Dnp] =
        {{7.11111111111111110e-02,3.55555555555555902e-02,-1.77777777777777951e-02,3.55555555555555416e-02,1.77777777777777118e-02,-8.88888888888895133e-03,-1.77777777777778020e-02,-8.88888888888893051e-03,4.44444444444443056e-03},
        {3.55555555555555902e-02,2.84444444444444444e-01,3.55555555555555208e-02,1.77777777777777500e-02,1.42222222222222222e-01,1.77777777777777569e-02,-8.88888888888894960e-03,-7.11111111111111666e-02,-8.88888888888890276e-03,},
        {-1.77777777777777951e-02,3.55555555555555208e-02,7.11111111111111249e-02,-8.88888888888892878e-03,1.77777777777776633e-02,3.55555555555555486e-02,4.44444444444441061e-03,-8.88888888888894266e-03,-1.77777777777777951e-02,},
        {3.55555555555555416e-02,1.77777777777777500e-02,-8.88888888888892878e-03,2.84444444444444444e-01,1.42222222222222250e-01,-7.11111111111111666e-02,3.55555555555555208e-02,1.77777777777777569e-02,-8.88888888888890276e-03,},
        {1.77777777777777118e-02,1.42222222222222222e-01,1.77777777777776633e-02,1.42222222222222250e-01,1.13777777777777778e+00,1.42222222222222167e-01,1.77777777777776008e-02,1.42222222222222167e-01,1.77777777777777812e-02,},
        {-8.88888888888895133e-03,1.77777777777777569e-02,3.55555555555555486e-02,-7.11111111111111666e-02,1.42222222222222167e-01,2.84444444444444444e-01,-8.88888888888894960e-03,1.77777777777777639e-02,3.55555555555555625e-02,},
        {-1.77777777777778020e-02,-8.88888888888894960e-03,4.44444444444441061e-03,3.55555555555555208e-02,1.77777777777776008e-02,-8.88888888888894960e-03,7.11111111111111249e-02,3.55555555555555347e-02,-1.77777777777777986e-02,},
        {-8.88888888888893051e-03,-7.11111111111111666e-02,-8.88888888888894266e-03,1.77777777777777569e-02,1.42222222222222167e-01,1.77777777777777639e-02,3.55555555555555347e-02,2.84444444444444389e-01,3.55555555555555694e-02,},
        {4.44444444444443056e-03,-8.88888888888890276e-03,-1.77777777777777951e-02,-8.88888888888890276e-03,1.77777777777777812e-02,3.55555555555555625e-02,-1.77777777777777986e-02,3.55555555555555694e-02,7.11111111111111249e-02,}};

double TestQuad_Dr[Dnp][Dnp] =
        {{-1.50000000000000000e+00,1.99999999999999978e+00,-4.99999999999999889e-01,-7.89491928622333323e-17,1.57898385724466665e-16,-7.89491928622333323e-17,3.94745964311166662e-17,-7.89491928622333323e-17,3.94745964311166662e-17},
        {-5.00000000000000111e-01,0.00000000000000000e+00,5.00000000000000111e-01,0.00000000000000000e+00,0.00000000000000000e+00,0.00000000000000000e+00,0.00000000000000000e+00,-0.00000000000000000e+00,0.00000000000000000e+00},
        {4.99999999999999889e-01,-1.99999999999999978e+00,1.50000000000000000e+00,7.89491928622333323e-17,-1.57898385724466665e-16,7.89491928622333323e-17,-3.94745964311166662e-17,7.89491928622333323e-17,-3.94745964311166662e-17,},
        {2.22044604925031308e-16,1.11022302462515654e-16,5.55111512312578270e-17,-1.50000000000000000e+00,1.99999999999999978e+00,-4.99999999999999889e-01,0.00000000000000000e+00,1.11022302462515654e-16,-5.55111512312578270e-17,},
        {-1.11022302462515642e-16,9.86864910777916654e-18,1.11022302462515654e-16,-5.00000000000000000e-01,-1.97372982155583331e-17,5.00000000000000000e-01,0.00000000000000000e+00,9.86864910777916654e-18,0.00000000000000000e+00,},
        {-2.22044604925031308e-16,-1.11022302462515654e-16,0.00000000000000000e+00,4.99999999999999778e-01,-1.99999999999999978e+00,1.50000000000000000e+00,1.11022302462515654e-16,-1.11022302462515654e-16,0.00000000000000000e+00,},
        {-1.57898385724466665e-16,-2.22044604925031308e-16,1.11022302462515654e-16,-7.89491928622333323e-17,1.57898385724466665e-16,-7.89491928622333323e-17,-1.50000000000000000e+00,2.00000000000000000e+00,-4.99999999999999944e-01,},
        {-1.11022302462515654e-16,0.00000000000000000e+00,1.11022302462515654e-16,0.00000000000000000e+00,0.00000000000000000e+00,0.00000000000000000e+00,-5.00000000000000000e-01,-0.00000000000000000e+00,5.00000000000000000e-01,},
        {7.89491928622333323e-17,7.89491928622333323e-17,1.50496898893632314e-16,3.94745964311166662e-17,-7.89491928622333323e-17,3.94745964311166662e-17,5.00000000000000000e-01,-2.00000000000000000e+00,1.50000000000000000e+00,}};

double TestQuad_Ds[Dnp][Dnp] = {
        {-1.5000000000000002e+00,-7.8949192862233332e-17,3.9474596431116666e-17,1.9999999999999998e+00,1.5789838572446666e-16,-7.8949192862233332e-17,-4.9999999999999978e-01,-7.8949192862233332e-17,3.9474596431116666e-17,},
        {-2.2204460492503131e-16,-1.5000000000000000e+00,-5.5511151231257827e-17,1.1102230246251565e-16,1.9999999999999998e+00,1.1102230246251565e-16,5.5511151231257827e-17,-4.9999999999999989e-01,-5.5511151231257827e-17,},
        {0.0000000000000000e+00,-7.8949192862233332e-17,-1.5000000000000000e+00,-2.2204460492503131e-16,1.5789838572446666e-16,2.0000000000000000e+00,2.2204460492503131e-16,-7.8949192862233332e-17,-5.0000000000000000e-01,},
        {-5.0000000000000011e-01,0.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,5.0000000000000011e-01,-0.0000000000000000e+00,0.0000000000000000e+00,},
        {-9.1285004246957324e-17,-5.0000000000000000e-01,-4.9343245538895833e-18,9.8686491077791665e-18,-1.9737298215558333e-17,9.8686491077791665e-18,1.0608797790862607e-16,5.0000000000000000e-01,-4.9343245538895833e-18,},
        {-1.1102230246251565e-16,0.0000000000000000e+00,-5.0000000000000000e-01,0.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,1.1102230246251565e-16,-0.0000000000000000e+00,5.0000000000000000e-01,},
        {4.9999999999999989e-01,7.8949192862233332e-17,-3.9474596431116666e-17,-1.9999999999999998e+00,-1.5789838572446666e-16,7.8949192862233332e-17,1.5000000000000000e+00,7.8949192862233332e-17,-3.9474596431116666e-17,},
        {0.0000000000000000e+00,4.9999999999999978e-01,5.5511151231257827e-17,-1.1102230246251565e-16,-1.9999999999999998e+00,-1.1102230246251565e-16,5.5511151231257827e-17,1.5000000000000000e+00,5.5511151231257827e-17,},
        {0.0000000000000000e+00,3.9474596431116666e-17,5.0000000000000000e-01,7.8949192862233332e-17,-7.8949192862233332e-17,-2.0000000000000000e+00,1.1102230246251565e-16,3.9474596431116666e-17,1.5000000000000000e+00,}};

double TestQuad_LIFT[Dnp][12] =
        {{4.5000000000000018e+00,0.0000000000000000e+00,2.2204460492503131e-16,1.5000000000000009e+00,-4.4408920985006262e-16,6.6613381477509392e-16,5.5511151231257827e-16,-9.4368957093138306e-16,1.5000000000000004e+00,2.2204460492503131e-16,-1.4432899320127035e-15,4.5000000000000009e+00},
        {0.0000000000000000e+00,4.5000000000000009e+00,-3.3306690738754696e-16,-7.5000000000000022e-01,4.4408920985006262e-16,-1.1102230246251565e-16,-5.5511151231257827e-17,1.5000000000000011e+00,-2.0816681711721685e-16,-1.1102230246251565e-16,6.9388939039072284e-16,-7.5000000000000000e-01},
        {2.2204460492503131e-16,-1.3322676295501878e-15,4.5000000000000009e+00,4.5000000000000018e+00,0.0000000000000000e+00,2.2204460492503131e-16,1.5000000000000009e+00,-4.4408920985006262e-16,6.1062266354383610e-16,5.5511151231257827e-16,-9.4368957093138306e-16,1.5000000000000004e+00},
        {-7.5000000000000022e-01,4.4408920985006262e-16,-1.1102230246251565e-16,-6.9388939039072284e-17,1.5000000000000011e+00,-2.2204460492503131e-16,-1.1102230246251565e-16,6.9388939039072284e-16,-7.5000000000000000e-01,0.0000000000000000e+00,4.5000000000000009e+00,-3.3306690738754696e-16},
        {9.0205620750793969e-17,-7.5000000000000011e-01,1.6653345369377348e-16,9.0205620750793969e-17,-7.5000000000000011e-01,1.6653345369377348e-16,8.3266726846886741e-17,-7.5000000000000011e-01,1.5959455978986625e-16,8.3266726846886741e-17,-7.5000000000000011e-01,1.5959455978986625e-16},
        {-1.1102230246251565e-16,6.6613381477509392e-16,-7.5000000000000000e-01,0.0000000000000000e+00,4.5000000000000009e+00,-3.3306690738754696e-16,-7.5000000000000022e-01,4.4408920985006262e-16,-1.3877787807814457e-16,-5.5511151231257827e-17,1.5000000000000011e+00,-2.0816681711721685e-16},
        {1.5000000000000009e+00,-4.4408920985006262e-16,6.6613381477509392e-16,4.9960036108132044e-16,-7.7715611723760958e-16,1.5000000000000007e+00,4.4408920985006262e-16,-1.7763568394002505e-15,4.5000000000000009e+00,4.5000000000000018e+00,0.0000000000000000e+00,0.0000000000000000e+00},
        {-6.9388939039072284e-17,1.5000000000000011e+00,-2.2204460492503131e-16,-8.3266726846886741e-17,6.1062266354383610e-16,-7.5000000000000022e-01,-1.1102230246251565e-16,4.5000000000000009e+00,-4.1633363423443370e-16,-7.5000000000000033e-01,3.8857805861880479e-16,-1.1102230246251565e-16},
        {5.5511151231257827e-16,-8.8817841970012523e-16,1.5000000000000004e+00,2.2204460492503131e-16,-1.3322676295501878e-15,4.5000000000000009e+00,4.5000000000000018e+00,-4.4408920985006262e-16,2.2204460492503131e-16,1.5000000000000009e+00,-3.3306690738754696e-16,5.5511151231257827e-16}};

#endif //DGOM_QUADTEST_H