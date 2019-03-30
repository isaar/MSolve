using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcpg
{
    internal interface IInterfaceFlexibilityMatrix
    {
        int Order { get; }
        void Multiply(Vector lhs, Vector rhs);
        Vector Multiply(Vector lhs); //TODO: This should be an extension
    }
}
