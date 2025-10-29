
gridMass = dimensionedScalar("zero", dimMass, 0.0);
gridMomentum = dimensionedVector("zero", dimMass*dimVelocity, vector::zero);
gridStress = dimensionedTensor("zero", dimPressure, tensor::zero);
gridForce = dimensionedVector("zero", dimForce/dimVolume, vector::zero);
