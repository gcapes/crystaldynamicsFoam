

//((fvc::ddt(rho)+fvc::div(mixture.rhoPhi())));



 fvScalarMatrix  rhoEqn
(


fvm::ddt(rho)
+fvm::div(fvc::flux(U),rho)

);

        rhoEqn.relax();
        fvConstraints.constrain(rhoEqn);

        rhoEqn.solve();