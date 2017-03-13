using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Integration.Rules;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.XFEM.Elements
{
    class IsoparametricQuad9: ContinuumElement2D
    {
        public IsoparametricQuad9(IReadOnlyList<Node2D> nodes, IIntegrationStrategyFactory2D integrationFactory):
            base (nodes, IsoparametricInterpolation2D.Quad9, GaussLegendre2D.Order3x3, integrationFactory)
        {
            if (nodes.Count != 9) throw new ArgumentException("A Quad4 finite element has 4 nodes, but "
                + nodes.Count + " nodes were provided.");
        }
    }
}
