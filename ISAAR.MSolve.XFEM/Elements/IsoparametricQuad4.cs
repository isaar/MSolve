using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Integration.Rules;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.XFEM.Elements
{
    // I could very easily make ContinuumElement2D generic on the Interpolation type. However, these concrete classes
    // are more readable (and expected in FEM software) and may be useful for geometric operations later on.
    class IsoparametricQuad4: ContinuumElement2D
    {
        public IsoparametricQuad4(IReadOnlyList<Node2D> nodes, IIntegrationStrategyFactory2D integrationFactory):
            base (nodes, IsoparametricInterpolation2D.Quad4, GaussLegendre2D.Order2x2, integrationFactory)
        {
            if (nodes.Count != 4) throw new ArgumentException("A Quad4 finite element has 4 nodes, but " 
                + nodes.Count + " nodes were provided.");
        }
    }
}
