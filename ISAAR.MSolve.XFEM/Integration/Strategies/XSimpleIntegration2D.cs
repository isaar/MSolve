using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Integration;
using ISAAR.MSolve.XFEM.Elements;

namespace ISAAR.MSolve.XFEM.Integration.Strategies
{
    // TODO: Remove the plain ContinuumElement, SimpleIntegration2D
    class XSimpleIntegration2D : IIntegrationStrategy2D<XContinuumElement2D>
    {
        public XSimpleIntegration2D()
        {
        }

        public IReadOnlyList<GaussPoint> GenerateIntegrationPoints(XContinuumElement2D element)
        {
            return element.StandardQuadrature.IntegrationPoints;
        }
    }
}
