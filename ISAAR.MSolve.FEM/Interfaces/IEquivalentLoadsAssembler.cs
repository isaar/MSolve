using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IEquivalentLoadsAssembler
    {
        IVector GetEquivalentNodalLoads(IVector solution, double constraintScalingFactor);
    }
}
