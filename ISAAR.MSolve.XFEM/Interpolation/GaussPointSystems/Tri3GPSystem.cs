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
    class Tri3GPSystem: IGaussPointSystem
    {
        // The order is important
        private static readonly INaturalPoint2D[] gaussPoins = new INaturalPoint2D[]
        {
            new NaturalPoint2D(1.0/6.0, 1.0/6.0),
            new NaturalPoint2D(2.0/3.0, 1.0/6.0),
            new NaturalPoint2D(1.0/6.0, 2.0/3.0)
        };

        private static readonly double[,] shapeFunctionsAtNodes;

        static Tri3GPSystem()
        {
            Func<double, double, double>[] shapeFunctions = new Func<double, double, double>[]
            {
                (r, s) => { return 1 - r - s; },
                (r, s) => { return r; },
                (r, s) => { return s; }
            };

            // Coordinates of nodes in the GP system: r=2*xi-1/3, s=2*eta-1/3
            double[,] nodeCoordinates = new double[,]
            {
                { -1.0/3.0, -1.0/3.0 },
                { 5.0/3.0, -1.0/3.0 },
                { -1.0/3.0, 5.0/3.0 }
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
            for (int node = 0; node < 3; ++node)
            {
                double tensorXX = 0, tensorYY = 0, tensorXY = 0;
                for (int gp = 0; gp < 3; ++gp)
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
