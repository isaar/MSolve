using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti
{
    internal interface IInterfaceProjection
    {
        Vector CalculateRigidBodyModesCoefficients(Vector flexibilityTimeslagrangeMultipliers,
            Vector boundaryDisplacements);

        void InitializeLagrangeMultipliers(Vector rigidBodyModesWork, Vector lagrange);

        void InvertCoarseProblemMatrix();

        void ProjectVector(Vector original, Vector projected);
    }
}
