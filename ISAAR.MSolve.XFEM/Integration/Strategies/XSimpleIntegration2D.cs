using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Integration.Points;

namespace ISAAR.MSolve.XFEM.Integration.Strategies
{
    // TODO: Remove the plain ContinuumElement, SimpleIntegration2D
    class XSimpleIntegration2D : IIntegrationStrategy2D<XContinuumElement2D>
    {
        public XSimpleIntegration2D()
        {
        }

        public IReadOnlyList<GaussPoint2D> GenerateIntegrationPoints(XContinuumElement2D element)
        {
            return element.StandardQuadrature.IntegrationPoints;
        }
    }
}
