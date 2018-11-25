using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: this and its concrete implementations should be in Analyzers project
//TODO: delete the original one in FEM.Interfaces
namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IEquivalentLoadsAssembler_v2
    {
        IVector GetEquivalentNodalLoads(IVectorView solution, double constraintScalingFactor);
    }
}
