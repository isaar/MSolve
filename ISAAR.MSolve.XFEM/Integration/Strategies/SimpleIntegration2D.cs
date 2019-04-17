using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Integration;
using ISAAR.MSolve.XFEM.Elements;

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

        public IReadOnlyList<GaussPoint> GenerateIntegrationPoints(ContinuumElement2D element)
        {
            return element.StandardQuadrature.IntegrationPoints;
        }
    }
}