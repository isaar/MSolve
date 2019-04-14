using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: Add Gray scale filter, as presented in "An efficient 3D topology optimization code written in Matlab (2014)", which
//      requires calculations during the OC optimizer
//TODO: Add filters presented in other publications with Matlab codes.
namespace ISAAR.MSolve.Optimization.Structural.Topology.SIMP.Filtering
{
    public interface IDensityFilter
    {
        void FilterSensitivities(Vector densities, ref Vector sensitivities); //TODO: perhaps use IVectorView?
    }
}
