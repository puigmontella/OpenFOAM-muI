/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2312                                  |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       surfaceScalarField;
    location    "150";
    object      phi.b;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 3 -1 0 0 0 0];

oriented        oriented;

internalField   nonuniform List<scalar> 
199
(
2.79632e-11
5.62313e-11
8.4689e-11
1.13253e-10
1.42151e-10
1.71117e-10
2.00416e-10
2.29805e-10
2.5951e-10
2.89328e-10
3.19444e-10
3.49701e-10
3.80234e-10
4.10935e-10
4.41894e-10
4.73046e-10
5.04442e-10
5.3605e-10
5.67892e-10
5.99961e-10
6.32262e-10
6.64797e-10
6.97567e-10
7.30575e-10
7.63822e-10
7.97313e-10
8.31044e-10
8.65025e-10
8.99253e-10
9.3373e-10
9.68463e-10
1.00345e-09
1.03869e-09
1.07419e-09
1.10996e-09
1.14599e-09
1.18228e-09
1.21884e-09
1.25568e-09
1.29278e-09
1.33017e-09
1.36783e-09
1.40577e-09
1.444e-09
1.48251e-09
1.5213e-09
1.56039e-09
1.59978e-09
1.63945e-09
1.67943e-09
1.7197e-09
1.76028e-09
1.80116e-09
1.84235e-09
1.88386e-09
1.92567e-09
1.9678e-09
2.01025e-09
2.05301e-09
2.09611e-09
2.13952e-09
2.18326e-09
2.22734e-09
2.27174e-09
2.31649e-09
2.36156e-09
2.40698e-09
2.45274e-09
2.49884e-09
2.5453e-09
2.59209e-09
2.63925e-09
2.68675e-09
2.73461e-09
2.78282e-09
2.83139e-09
2.88033e-09
2.92963e-09
2.97929e-09
3.02931e-09
3.07971e-09
3.13047e-09
3.1816e-09
3.23311e-09
3.28498e-09
3.33723e-09
3.38985e-09
3.44284e-09
3.49621e-09
3.54995e-09
3.60406e-09
3.65855e-09
3.71341e-09
3.76865e-09
3.82426e-09
3.88023e-09
3.93658e-09
3.99329e-09
4.05036e-09
4.10781e-09
4.1656e-09
4.22376e-09
4.28227e-09
4.34112e-09
4.40032e-09
4.45986e-09
4.51973e-09
4.57994e-09
4.64045e-09
4.70129e-09
4.76244e-09
4.82387e-09
4.88561e-09
4.94762e-09
5.00989e-09
5.07244e-09
5.13522e-09
5.19825e-09
5.2615e-09
5.32494e-09
5.38859e-09
5.45241e-09
5.51639e-09
5.58051e-09
5.64475e-09
5.7091e-09
5.77352e-09
5.83801e-09
5.90253e-09
5.96706e-09
6.03157e-09
6.09605e-09
6.16045e-09
6.22475e-09
6.28892e-09
6.35293e-09
6.41673e-09
6.48031e-09
6.5436e-09
6.60661e-09
6.66925e-09
6.73151e-09
6.79334e-09
6.85469e-09
6.91553e-09
6.97579e-09
7.03545e-09
7.09445e-09
7.15273e-09
7.21026e-09
7.26697e-09
7.32281e-09
7.37774e-09
7.43169e-09
7.48461e-09
7.53646e-09
7.58717e-09
7.63668e-09
7.68495e-09
7.73192e-09
7.77753e-09
7.82175e-09
7.86451e-09
7.90576e-09
7.94545e-09
7.98352e-09
8.01997e-09
8.05471e-09
8.0877e-09
8.11892e-09
8.14831e-09
8.17582e-09
8.20143e-09
8.22509e-09
8.24678e-09
8.2664e-09
8.28394e-09
8.29935e-09
8.31255e-09
8.32349e-09
8.33207e-09
8.33822e-09
8.34181e-09
8.34276e-09
8.34084e-09
8.33574e-09
8.32745e-09
8.31214e-09
8.29259e-09
8.25759e-09
5.62348e-09
1.67896e-09
1.51524e-10
1.28093e-11
1.08019e-12
5.26291e-14
-4.68576e-14
-5.77412e-14
-5.91226e-14
)
;

boundaryField
{
    inlet
    {
        type            cyclic;
        value           nonuniform List<scalar> 
200
(
1.58808e-28
4.77322e-28
7.97807e-28
1.12208e-27
1.45072e-27
1.78561e-27
2.12839e-27
2.48024e-27
2.83866e-27
3.20214e-27
3.56658e-27
3.93304e-27
4.291e-27
4.65826e-27
5.02876e-27
5.40155e-27
5.77525e-27
6.14264e-27
6.511e-27
6.88347e-27
7.25642e-27
7.63433e-27
7.98572e-27
8.34621e-27
8.69742e-27
9.05167e-27
9.44544e-27
9.83729e-27
1.02278e-26
1.06194e-26
1.10134e-26
1.14008e-26
1.17849e-26
1.21839e-26
1.25744e-26
1.2969e-26
1.33633e-26
1.375e-26
1.41364e-26
1.45211e-26
1.4898e-26
1.52793e-26
1.56504e-26
1.60157e-26
1.63722e-26
1.67191e-26
1.70531e-26
1.7372e-26
1.76715e-26
1.79475e-26
1.81938e-26
1.855e-26
1.90184e-26
1.94541e-26
1.98604e-26
2.04007e-26
2.09234e-26
2.1268e-26
2.16025e-26
2.21013e-26
2.25965e-26
2.30951e-26
2.32291e-26
2.29778e-26
2.30833e-26
2.33395e-26
2.35409e-26
2.36871e-26
2.41851e-26
2.46205e-26
2.49961e-26
2.53107e-26
2.60287e-26
2.67125e-26
2.69068e-26
2.70775e-26
2.72197e-26
2.73201e-26
2.73699e-26
2.78732e-26
2.83224e-26
2.87237e-26
2.90843e-26
2.94069e-26
2.96935e-26
2.99435e-26
3.01552e-26
3.03229e-26
3.04429e-26
3.05023e-26
2.92355e-26
2.90832e-26
3.00482e-26
3.08715e-26
3.15758e-26
3.2213e-26
3.27502e-26
3.3239e-26
3.37125e-26
3.41452e-26
3.46038e-26
3.50635e-26
3.55267e-26
3.60211e-26
3.57797e-26
3.55571e-26
3.61302e-26
3.67181e-26
3.73407e-26
3.71839e-26
3.70625e-26
3.69568e-26
3.68491e-26
3.67194e-26
3.65524e-26
3.63209e-26
3.5084e-26
3.46083e-26
3.39268e-26
3.49061e-26
3.56746e-26
3.62616e-26
3.66969e-26
3.6005e-26
3.29997e-26
3.37708e-26
3.42764e-26
3.44177e-26
3.53192e-26
3.59382e-26
3.63045e-26
3.64283e-26
3.63094e-26
3.70817e-26
3.87993e-26
4.03941e-26
4.19534e-26
4.3566e-26
4.52482e-26
4.59323e-26
4.69547e-26
4.82949e-26
5.00348e-26
5.21921e-26
5.23676e-26
5.32221e-26
5.47504e-26
5.70235e-26
6.01141e-26
5.87771e-26
5.84943e-26
5.91799e-26
6.07161e-26
6.0484e-26
6.14255e-26
6.34871e-26
6.38555e-26
6.53364e-26
6.51634e-26
6.64605e-26
6.91153e-26
6.72301e-26
6.66292e-26
6.74355e-26
6.66905e-26
6.72855e-26
6.92404e-26
6.62911e-26
6.4736e-26
6.44803e-26
6.54503e-26
6.76713e-26
6.80082e-26
6.98671e-26
7.33377e-26
7.18295e-26
7.18971e-26
7.02256e-26
7.03147e-26
7.20209e-26
6.82912e-26
6.95905e-26
6.9057e-26
6.64392e-26
6.51848e-26
6.51607e-26
6.27897e-26
6.15876e-26
6.15585e-26
6.24537e-26
6.47225e-26
6.33689e-26
5.91982e-26
4.78852e-26
4.5136e-26
3.71258e-26
2.89624e-26
2.07208e-26
1.24429e-26
4.14894e-27
)
;
    }
    outlet
    {
        type            cyclic;
        value           nonuniform List<scalar> 
200
(
-1.58808e-28
-4.77322e-28
-7.97807e-28
-1.12208e-27
-1.45072e-27
-1.78561e-27
-2.12839e-27
-2.48024e-27
-2.83866e-27
-3.20214e-27
-3.56658e-27
-3.93304e-27
-4.291e-27
-4.65826e-27
-5.02876e-27
-5.40155e-27
-5.77525e-27
-6.14264e-27
-6.511e-27
-6.88347e-27
-7.25642e-27
-7.63433e-27
-7.98572e-27
-8.34621e-27
-8.69742e-27
-9.05167e-27
-9.44544e-27
-9.83729e-27
-1.02278e-26
-1.06194e-26
-1.10134e-26
-1.14008e-26
-1.17849e-26
-1.21839e-26
-1.25744e-26
-1.2969e-26
-1.33633e-26
-1.375e-26
-1.41364e-26
-1.45211e-26
-1.4898e-26
-1.52793e-26
-1.56504e-26
-1.60157e-26
-1.63722e-26
-1.67191e-26
-1.70531e-26
-1.7372e-26
-1.76715e-26
-1.79475e-26
-1.81938e-26
-1.855e-26
-1.90184e-26
-1.94541e-26
-1.98604e-26
-2.04007e-26
-2.09234e-26
-2.1268e-26
-2.16025e-26
-2.21013e-26
-2.25965e-26
-2.30951e-26
-2.32291e-26
-2.29778e-26
-2.30833e-26
-2.33395e-26
-2.35409e-26
-2.36871e-26
-2.41851e-26
-2.46205e-26
-2.49961e-26
-2.53107e-26
-2.60287e-26
-2.67125e-26
-2.69068e-26
-2.70775e-26
-2.72197e-26
-2.73201e-26
-2.73699e-26
-2.78732e-26
-2.83224e-26
-2.87237e-26
-2.90843e-26
-2.94069e-26
-2.96935e-26
-2.99435e-26
-3.01552e-26
-3.03229e-26
-3.04429e-26
-3.05023e-26
-2.92355e-26
-2.90832e-26
-3.00482e-26
-3.08715e-26
-3.15758e-26
-3.2213e-26
-3.27502e-26
-3.3239e-26
-3.37125e-26
-3.41452e-26
-3.46038e-26
-3.50635e-26
-3.55267e-26
-3.60211e-26
-3.57797e-26
-3.55571e-26
-3.61302e-26
-3.67181e-26
-3.73407e-26
-3.71839e-26
-3.70625e-26
-3.69568e-26
-3.68491e-26
-3.67194e-26
-3.65524e-26
-3.63209e-26
-3.5084e-26
-3.46083e-26
-3.39268e-26
-3.49061e-26
-3.56746e-26
-3.62616e-26
-3.66969e-26
-3.6005e-26
-3.29997e-26
-3.37708e-26
-3.42764e-26
-3.44177e-26
-3.53192e-26
-3.59382e-26
-3.63045e-26
-3.64283e-26
-3.63094e-26
-3.70817e-26
-3.87993e-26
-4.03941e-26
-4.19534e-26
-4.3566e-26
-4.52482e-26
-4.59323e-26
-4.69547e-26
-4.82949e-26
-5.00348e-26
-5.21921e-26
-5.23676e-26
-5.32221e-26
-5.47504e-26
-5.70235e-26
-6.01141e-26
-5.87771e-26
-5.84943e-26
-5.91799e-26
-6.07161e-26
-6.0484e-26
-6.14255e-26
-6.34871e-26
-6.38555e-26
-6.53364e-26
-6.51634e-26
-6.64605e-26
-6.91153e-26
-6.72301e-26
-6.66292e-26
-6.74355e-26
-6.66905e-26
-6.72855e-26
-6.92404e-26
-6.62911e-26
-6.4736e-26
-6.44803e-26
-6.54503e-26
-6.76713e-26
-6.80082e-26
-6.98671e-26
-7.33377e-26
-7.18295e-26
-7.18971e-26
-7.02256e-26
-7.03147e-26
-7.20209e-26
-6.82912e-26
-6.95905e-26
-6.9057e-26
-6.64392e-26
-6.51848e-26
-6.51607e-26
-6.27897e-26
-6.15876e-26
-6.15585e-26
-6.24537e-26
-6.47225e-26
-6.33689e-26
-5.91982e-26
-4.78852e-26
-4.5136e-26
-3.71258e-26
-2.89624e-26
-2.07208e-26
-1.24429e-26
-4.14894e-27
)
;
    }
    top
    {
        type            fixedValue;
        value           uniform 0;
    }
    bottom
    {
        type            fixedValue;
        value           uniform 0;
    }
    backPlane
    {
        type            empty;
        value           nonuniform List<scalar> 0();
    }
    frontPlane
    {
        type            empty;
        value           nonuniform List<scalar> 0();
    }
}


// ************************************************************************* //