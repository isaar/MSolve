using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Tensors;

namespace ISAAR.MSolve.XFEM.Interpolation.GaussPointSystems
{
    // TODO: Use the existing interoplation classes. 
    // Would be ideal to use the existing quadrature too, but the order is different.
    class Quad4GPSystem: IGaussPointSystem
    {
        // The order is important
        private static readonly INaturalPoint2D[] gaussPoins = new INaturalPoint2D[]
        {
            new NaturalPoint2D(-1/Math.Sqrt(3), -1/Math.Sqrt(3)),
            new NaturalPoint2D(-1/Math.Sqrt(3), 1/Math.Sqrt(3)),
            new NaturalPoint2D(1/Math.Sqrt(3), 1/Math.Sqrt(3)),
            new NaturalPoint2D(-1/Math.Sqrt(3), 1/Math.Sqrt(3))
        };

        private static readonly double[,] shapeFunctionsAtNodes;

        static Quad4GPSystem()
        {
            Func<double, double, double>[] shapeFunctions = new Func<double, double, double>[]
            {
                (r, s) => { return 0.25 * (1 - r) * (1 - s); },
                (r, s) => { return 0.25 * (1 + r) * (1 - s); },
                (r, s) => { return 0.25 * (1 + r) * (1 + s); },
                (r, s) => { return 0.25 * (1 - r) * (1 + s); }
            };

            // Coordinates of nodes in the GP system: r=sqrt(3)*xi, s=sqrt(3)*eta
            double[,] nodeCoordinates = new double[,]
            {
                { -Math.Sqrt(3), -Math.Sqrt(3) },
                { Math.Sqrt(3), -Math.Sqrt(3) },
                { Math.Sqrt(3), Math.Sqrt(3) },
                { -Math.Sqrt(3), Math.Sqrt(3) }
            };

            shapeFunctionsAtNodes = new double[4, 4];
            for (int node = 0; node < 4; ++node)
            {
                double r = nodeCoordinates[node, 0];
                double s = nodeCoordinates[node, 1];
                for (int gp = 0; gp < 4; ++gp)
                {
                    shapeFunctionsAtNodes[node, gp] = shapeFunctions[gp](r, s);
                }
            }
        }

        public IReadOnlyList<INaturalPoint2D> GaussPoints { get { return gaussPoins; } }

        public IReadOnlyList<Tensor2D> ExtrapolateTensorFromGaussPointsToNodes(IReadOnlyList<Tensor2D> tensorsAtGPs)
        {
            var nodalTensors = new Tensor2D[4];
            for (int node = 0; node < 4; ++node)
            {
                double tensorXX = 0, tensorYY = 0, tensorXY = 0;
                for (int gp = 0; gp < 4; ++gp)
                {
                    double shapeValue = shapeFunctionsAtNodes[node, gp];
                    tensorXX += shapeValue * tensorsAtGPs[gp].XX;
                    tensorYY += shapeValue * tensorsAtGPs[gp].YY;
                    tensorXY += shapeValue * tensorsAtGPs[gp].XY;
                }
                nodalTensors[node] = new Tensor2D(tensorXX, tensorYY, tensorXY);
            }
            return nodalTensors;
        }
    }
}
