using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Not sure if there needs to be an interface. In FETI-1 there is only one flexibility matrix. However in FETI-DP there is
//      FIrr and FIrc.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg
{
    internal interface IInterfaceFlexibilityMatrix 
    {
        int Order { get; }
        void Multiply(Vector lhs, Vector rhs);
        Vector Multiply(Vector lhs); //TODO: This should be an extension
    }
}
