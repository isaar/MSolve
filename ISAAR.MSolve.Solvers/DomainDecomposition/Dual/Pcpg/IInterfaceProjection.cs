using System;
using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcpg
{
    internal interface IInterfaceProjection
    {
        Vector CalcParticularLagrangeMultipliers(Vector rigidBodyModesWork);

        Vector CalcRigidBodyModesCoefficients(Vector flexibilityTimeslagrangeMultipliers,
            Vector boundaryDisplacements);

        void InvertCoarseProblemMatrix();

        void ProjectVector(Vector original, Vector projected, bool transposeProjector);
    }
}
