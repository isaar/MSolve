using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: The linear system classes need to be redesigned taking into account a) 1 linear system per subdomain, b) 1 IVector that
//      hides many concrete Vectors for the subdomains.
//TODO: I would like to make this internal, but it does not work.
namespace ISAAR.MSolve.Solvers.LinearSystems
{
    public interface ISingleSubdomainLinearSystem : ILinearSystem
    {
        Vector RhsConcrete { get; set; }
        Vector SolutionConcrete { get; set; }
        Vector CreateZeroVectorConcrete();
    }
}
