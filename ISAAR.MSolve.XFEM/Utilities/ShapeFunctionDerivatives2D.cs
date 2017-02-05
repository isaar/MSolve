using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Utilities
{
    //TODO: add methods to facilitate building without breaking immutability
    class ShapeFunctionDerivatives2D
    {
        private readonly double[,] data;
        public int NodesCount { get; }

        // Deep copies of input arrays
        public ShapeFunctionDerivatives2D(double[] xiDerivatives, double[] etaDerivatives)
        {
            Debug.Assert(xiDerivatives.Length == etaDerivatives.Length);
            this.NodesCount = xiDerivatives.Length;
            this.data = new double[NodesCount, 2];
            for (int nodeIndex = 0; nodeIndex < NodesCount; ++nodeIndex)
            {
                data[nodeIndex, 0] = xiDerivatives[nodeIndex];
                data[nodeIndex, 1] = etaDerivatives[nodeIndex];
            }
        }

        public Tuple<double, double> DerivativesOfNode(int nodeIndex)
        {
            return new Tuple<double, double>(data[nodeIndex, 0], data[nodeIndex, 1]);
        }

        public double XiDerivativeOfNode(int nodeIndex)
        {
            return data[nodeIndex, 0];
        }

        public double EtaDerivativeOfNode(int nodeIndex)
        {
            return data[nodeIndex, 1];
        }

        public double[] XiDerivatives()
        {
            double[] result = new double[4];
            for (int nodeIndex = 0; nodeIndex < NodesCount; ++nodeIndex)
            {
                result[nodeIndex] = data[nodeIndex, 0];
            }
            return result;
        }

        public double[] EtaDerivatives()
        {
            double[] result = new double[4];
            for (int nodeIndex = 0; nodeIndex < NodesCount; ++nodeIndex)
            {
                result[nodeIndex] = data[nodeIndex, 1];
            }
            return result;
        }
    }
}
