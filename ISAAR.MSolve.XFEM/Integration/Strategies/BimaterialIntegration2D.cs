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
    class BimaterialIntegration2D: IIntegrationStrategy2D<XContinuumElement2D>
    {
        private readonly IIntegrationRule2D<XContinuumElement2D> integrationRule;
        private readonly MaterialInterface2D bimaterialInterface;

        /// <summary>
        /// This may be updated throughout the analysis.
        /// </summary>
        private Dictionary<GaussPoint2D, IFiniteElementMaterial2D> pointsAndMaterials;

        public BimaterialIntegration2D(IIntegrationRule2D<XContinuumElement2D> integrationRule, 
            MaterialInterface2D bimaterialInterface)
        {
            this.integrationRule = integrationRule;
            this.bimaterialInterface = bimaterialInterface;
        }

        public IReadOnlyDictionary<GaussPoint2D, IFiniteElementMaterial2D> GetIntegrationPointsAndMaterials(
            XContinuumElement2D element)
        {
            if (pointsAndMaterials == null) Initialize(element); // TODO: This should be thread safe.
            return pointsAndMaterials;
        }

        private void Initialize(XContinuumElement2D element)
        {
            pointsAndMaterials = new Dictionary<GaussPoint2D, IFiniteElementMaterial2D>();
            foreach (GaussPoint2D gaussPoint in integrationRule.GenerateIntegrationPoints(element))
            {
                // TODO: remove the train wreck code. Actually this is everywhere since element types expose most of 
                // their composed classes, instead of wrapping their functionality. (Should they?)
                ICartesianPoint2D cartesianPoint =
                    element.Interpolation.TransformNaturalToCartesian(element.Nodes, gaussPoint);

                IFiniteElementMaterial2D material;
                bool onInterface = bimaterialInterface.LocatePoint(cartesianPoint, out material);
                if (onInterface)
                {
                    throw new ArgumentException("The provided integration rule generates a point " +
                        "that lies directly on the bimaterial interface. Use another one.");
                }

                pointsAndMaterials.Add(gaussPoint, material.Clone());
            }
        }

        public void Update(XContinuumElement2D element)
        {
            throw new NotImplementedException();
        }
    }
}
