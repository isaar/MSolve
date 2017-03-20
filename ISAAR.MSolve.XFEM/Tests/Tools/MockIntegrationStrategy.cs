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
        private readonly IFiniteElementMaterial2D material;

        public MockIntegrationStrategy(GaussPoint2D[] knownIntegrationPoints, IFiniteElementMaterial2D material)
        {
            this.knownIntegrationPoints = knownIntegrationPoints;
            this.material = material;
        }

        public IReadOnlyDictionary<GaussPoint2D, IFiniteElementMaterial2D> GetIntegrationPointsAndMaterials(XContinuumElement2D element)
        {
            var pointsAndMaterials = new Dictionary<GaussPoint2D, IFiniteElementMaterial2D>();
            foreach (GaussPoint2D gaussPoint in knownIntegrationPoints)
            {
                pointsAndMaterials.Add(gaussPoint, material.Clone());
            }
            return pointsAndMaterials;
        }

        public void Update(XContinuumElement2D element)
        {

        }
    }
}
