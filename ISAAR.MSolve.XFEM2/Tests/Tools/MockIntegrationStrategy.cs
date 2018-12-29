using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.XFEM.Tests.Tools

{
    class MockIntegrationStrategy: IIntegrationStrategy2D<XContinuumElement2D>
    {
        private readonly GaussPoint2D[] knownIntegrationPoints;

        public MockIntegrationStrategy(GaussPoint2D[] knownIntegrationPoints)
        {
            this.knownIntegrationPoints = knownIntegrationPoints;
        }

        public IReadOnlyList<GaussPoint2D> GenerateIntegrationPoints(XContinuumElement2D element)
        {
            return knownIntegrationPoints;
        }
    }
}
