using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.XFEM.Integration.Strategies
{
    /// <summary>
    /// Integration with standard quadratures. Lightweight and fast since all integration points are cached.
    /// TODO: Needs a better name
    /// </summary>
    class SimpleIntegration2D : IIntegrationStrategy2D<ContinuumElement2D>
    {
        public SimpleIntegration2D()
        {
        }

        public IReadOnlyList<GaussPoint2D> GenerateIntegrationPoints(ContinuumElement2D element)
        {
            return element.StandardQuadrature.IntegrationPoints;
        }
    }
}