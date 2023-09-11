{

dimensionedScalar shear_mod ("shear_mod", dimensionSet(1,-1,-2,0,0,0,0),80e9);

volTensorField TgradU = (fvc::grad(U)().T());

omega_dot = 0.5*(fvc::grad(U)-TgradU);

epsilon_dot = 0.5*(fvc::grad(U)-TgradU) - (1.0/3.0)*(fvc::div(phi))*I;

dsdt = 2.0*shear_mod*(epsilon_dot/* - epsilon_P_dot*/) - ((s & omega_dot) - (omega_dot & s));

s=s + runTime.deltaT()*(dsdt);

    // solve
    // (
    //     fvm::ddt(s) ==
    //     dsdt
    // );




/*
tensor T1(1, 2, 3, 4, 5, 6, 7, 8, 9);

tensor T2(10, 11, 12, 13, 14, 15, 16, 17, 18);

tensor T3 = (T1 & T2);

tensor T4 = (T2 & T1);

Info<<"T3 ="<<T3<<endl;

Info<<"T4 ="<<T4<<endl;

exit(0);

*/
}