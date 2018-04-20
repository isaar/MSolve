using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.XFEM.Geometry.CoordinateSystems
{
    class PolarPoint2D
    {
        public double R { get; }
        public double Theta { get; }

        public PolarPoint2D (double r, double theta)
        {
            this.R = r;
            this.Theta = theta;
        }

        public PolarPoint2D(Vector2 polarCoordinates)
        {
            this.R = polarCoordinates[0];
            this.Theta = polarCoordinates[1];
        }

        /// <summary>
        /// Returns an array of length = 2 containing the coordinates. The point is not mutated if the array is changed.
        /// </summary>
        public Vector2 Coordinates { get { return Vector2.Create(R, Theta); } }
    }
}
