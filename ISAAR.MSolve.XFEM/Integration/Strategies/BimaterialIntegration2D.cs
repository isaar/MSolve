using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Integration.Rules;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.XFEM.Integration.Strategies
{
    class BimaterialIntegration2D: IIntegrationStrategy2D
    {
        public class Factory : IIntegrationStrategyFactory2D
        {
            private readonly IIntegrationRule2D<ContinuumElement2D> integrationRule;
            private readonly MaterialInterface2D bimaterialInterface;

            public Factory(IIntegrationRule2D<ContinuumElement2D> integrationRule, MaterialInterface2D bimaterialInterface)
            {
                this.integrationRule = integrationRule;
                this.bimaterialInterface = bimaterialInterface;
            }

            public IIntegrationStrategy2D CreateStrategy(ContinuumElement2D element)
            {
                return new BimaterialIntegration2D(element, integrationRule, bimaterialInterface);
            }
        }

        private readonly IEnumerable<Tuple<GaussPoint2D, IFiniteElementMaterial2D>> pointsAndMaterials;

        private BimaterialIntegration2D(ContinuumElement2D element, 
            IIntegrationRule2D<ContinuumElement2D> integrationRule, MaterialInterface2D bimaterialInterface)
        {
            IReadOnlyList<GaussPoint2D> points = integrationRule.GenerateIntegrationPoints(element);

            var pairs = new Tuple<GaussPoint2D, IFiniteElementMaterial2D>[points.Count];
            for (int i = 0; i < points.Count; ++i)
            {
                ICartesianPoint2D cartesianPoint = 
                    element.Interpolation.TransformNaturalToCartesian(element.Nodes, points[i]); 

                IFiniteElementMaterial2D material;
                bool onInterface = bimaterialInterface.LocatePoint(cartesianPoint, out material);
                if (onInterface)
                {
                    throw new ArgumentException("This integration rule generates a point " + 
                        "that lies directly on the bimaterial interface. Use another one.");
                }

                pairs[i] = new Tuple<GaussPoint2D, IFiniteElementMaterial2D>(points[i], material.Clone());
            }
            pointsAndMaterials = pairs;
        }

        public IEnumerable<Tuple<GaussPoint2D, IFiniteElementMaterial2D>> GetIntegrationPointsAndMaterials()
        {
            return pointsAndMaterials;
        }

        public void Update()
        {
            throw new NotImplementedException();
        }
    }
}
