using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: this and its concrete implementations should be in Analyzers project
namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IEquivalentLoadsAssembler_v2
    {
        IVector GetEquivalentNodalLoads(IVectorView solution, double constraintScalingFactor);
    }
}
