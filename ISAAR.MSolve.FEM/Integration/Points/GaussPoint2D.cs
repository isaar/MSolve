using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Integration.Points
{
    /// <summary>
    /// Integration point (coordinates + weight) defined in the natural coordinate system of a finite element. Immutable.
    /// </summary>
    public class GaussPoint2D
    {
        public GaussPoint2D(double xi, double eta, double weight)
        {
            this.Xi = xi;
            this.Eta = eta;
            this.Weight = weight;
        }

        public double Xi { get; }
        public double Eta { get; }
        public double Weight { get; }
    }
}
