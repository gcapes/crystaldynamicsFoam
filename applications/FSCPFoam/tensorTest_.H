// Add this as a new file: tensorTest.H
// Include it once at the start of the time loop to verify operations

{
    Info << "\n===== TENSOR OPERATION VERIFICATION =====" << endl;
    
    // Test case: uniaxial compression F_yy = 0.99
    tensor F_test(1, 0, 0, 0, 0.99, 0, 0, 0, 1);
    
    Info << "Test F = " << F_test << endl;
    Info << "F_yy (component 4) = " << F_test.component(4) << endl;
    
    // C = F^T * F
    tensor C_test = F_test.T() & F_test;
    Info << "C = F^T & F = " << C_test << endl;
    Info << "C_yy = " << C_test.component(4) << " (should be " << 0.99*0.99 << ")" << endl;
    
    // E = 0.5*(C - I)
    tensor E_test = 0.5 * (C_test - tensor::I);
    Info << "E = 0.5*(C-I) = " << E_test << endl;
    Info << "E_yy = " << E_test.component(4) << " (should be " << 0.5*(0.99*0.99 - 1) << ")" << endl;
    
    // S = lambda*tr(E)*I + 2*mu*E
    scalar lambda_val = 39645522388.05971;
    scalar mu_val = 18656716417.910446;
    scalar trE_test = tr(E_test);
    Info << "tr(E) = " << trE_test << endl;
    
    tensor S_test = lambda_val * trE_test * tensor::I + 2.0 * mu_val * E_test;
    Info << "S = " << S_test << endl;
    Info << "S_yy = " << S_test.component(4) << " (should be NEGATIVE for compression)" << endl;
    
    // sigma = (1/J) * F * S * F^T
    scalar J_test = det(F_test);
    Info << "J = det(F) = " << J_test << " (should be ~0.99)" << endl;
    
    tensor sigma_test = (1.0/J_test) * (F_test & S_test & F_test.T());
    Info << "sigma = (1/J)*F*S*F^T = " << sigma_test << endl;
    Info << "sigma_yy = " << sigma_test.component(4) << " (should be NEGATIVE for compression)" << endl;
    
    // Now check what OpenFOAM's & operator actually does
    Info << "\n--- Verifying & operator ---" << endl;
    tensor A(1, 2, 3, 4, 5, 6, 7, 8, 9);
    tensor B(9, 8, 7, 6, 5, 4, 3, 2, 1);
    tensor AB = A & B;
    Info << "A = " << A << endl;
    Info << "B = " << B << endl;
    Info << "A & B = " << AB << endl;
    
    // Manual calculation: (A&B)_ij = A_ik * B_kj
    // (A&B)_00 = A_00*B_00 + A_01*B_10 + A_02*B_20 = 1*9 + 2*6 + 3*3 = 9+12+9 = 30
    Info << "(A&B)_00 should be 30, got " << AB.component(0) << endl;
    
    // (A&B)_11 = A_10*B_01 + A_11*B_11 + A_12*B_21 = 4*8 + 5*5 + 6*2 = 32+25+12 = 69
    Info << "(A&B)_11 should be 69, got " << AB.component(4) << endl;
    
    Info << "=====================================\n" << endl;
}
