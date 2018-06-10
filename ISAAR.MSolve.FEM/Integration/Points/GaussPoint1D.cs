using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Integration.Points
{
    /// <summary>
    /// Integration point (coordinates + weight) defined in the natural coordinate system of a finite element. Immutable.
    /// </summary>
    public class GaussPoint1D
    {
        public GaussPoint1D(double xi, double weight)
        {
            this.Xi = xi;
            this.Weight = weight;
        }

        public double Xi { get; }
        public double Weight { get; }
    }
}
