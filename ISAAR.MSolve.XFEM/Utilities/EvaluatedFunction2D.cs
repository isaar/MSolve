using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.XFEM.Utilities
{
    /// <summary>
    /// Data transfer object to store and pass around the value and derivatives of a 2D function, evaluated at some 
    /// It mainly serves to avoid obscure Tuple<double, Tuple<double, double>> objects.
    /// </summary>
    class EvaluatedFunction2D
    {
        public EvaluatedFunction2D(double value, Vector2 cartesianDerivatives)
        {
            this.Value = value;
            this.CartesianDerivatives = cartesianDerivatives;
        }

        public double Value { get; }
        public Vector2 CartesianDerivatives { get; }
    }
}
