using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Integration;
using ISAAR.MSolve.XFEM.Elements;

namespace ISAAR.MSolve.XFEM.Integration
{
    public class SimpleIntegration2D : IIntegrationStrategy2D<XContinuumElement2D>
    {
        public SimpleIntegration2D()
        {
        }

        public IReadOnlyList<GaussPoint> GenerateIntegrationPoints(XContinuumElement2D element)
        {
            return element.StandardQuadrature.IntegrationPoints;
        }
    }
}
