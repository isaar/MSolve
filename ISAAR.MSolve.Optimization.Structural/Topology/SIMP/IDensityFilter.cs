using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Optimization.Structural.Topology.SIMP
{
    public interface IDensityFilter
    {
        void FilterSensitivities(Vector densities, ref Vector sensitivities); //TODO: perhaps use IVectorView?
    }
}
