using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Interpolation
{
    class ShapeFunctionDerivatives2D
    {
        private readonly double[,] naturalDerivatives; // An immutable collection might be better

        public Jacobian2D Jacobian { get; }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="naturalDerivatives">They are not deep copied. Each row corresponds to a node. 
        ///     Column 1 stores dN/dXi, column 2 stores dN/dEta.</param>
        /// <param name="jacobian"></param>
        public ShapeFunctionDerivatives2D(double[,] naturalDerivatives, Jacobian2D jacobian) // No deep copy
        {
            this.naturalDerivatives = naturalDerivatives;
            this.Jacobian = jacobian;
        }

        public Tuple<double, double> NaturalDerivativesOfNode(int nodeIndex) // TODO: use sth safer than indices
        {
            return new Tuple<double, double>(naturalDerivatives[nodeIndex, 0], naturalDerivatives[nodeIndex, 1]);
        }

        // Perhaps I should only store cartesian derivatives to avoid recalculating them.
        public Tuple<double, double> CartesianDerivativesOfNode(int nodeIndex) 
        {
            return Jacobian.TransformNaturalDerivativesToCartesian(
                naturalDerivatives[nodeIndex, 0], naturalDerivatives[nodeIndex, 1]);
        }


    }
}
