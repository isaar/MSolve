using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.Enrichments.Items.CrackTip;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Integration.Rules;
using ISAAR.MSolve.XFEM.Integration.Strategies;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.XFEM.Tests.Khoei
{
    class JintegralStrategy: IIntegrationStrategy2D<XContinuumElement2D>
    {
        /// <summary>
        /// Might need it later when updating. TODO: If I had a MaterialView, I could store the same material object 
        /// for all elements and save up memory.
        /// </summary>
        private readonly IFiniteElementMaterial2D commonMaterial;
        private readonly IStandardQuadrature2D standardIntegrationRule;
        private readonly IIntegrationRule2D<XContinuumElement2D> enrichedIntegrationRule;

        /// <summary>
        /// This may be updated throughout the analysis.
        /// </summary>
        private Dictionary<GaussPoint2D, IFiniteElementMaterial2D> pointsAndMaterials;

        public JintegralStrategy(IStandardQuadrature2D stdIntegrationRule, 
            RectangularSubgridIntegration2D<XContinuumElement2D> enrichedIntegrationRule,
            IFiniteElementMaterial2D commonMaterial)
        {
            this.standardIntegrationRule = stdIntegrationRule;
            this.enrichedIntegrationRule = enrichedIntegrationRule;
            this.commonMaterial = commonMaterial.Clone(); // The object passed might be mutated later on 
        }

        public IReadOnlyDictionary<GaussPoint2D, IFiniteElementMaterial2D> GetIntegrationPointsAndMaterials(
            XContinuumElement2D element)
        {
            if (pointsAndMaterials == null) Initialize(element); // TODO: This should be thread safe.
            return pointsAndMaterials;
        }

        private void Initialize(XContinuumElement2D element)
        {
            IReadOnlyList<GaussPoint2D> integrationPoints;
            if (element.EnrichmentItems.Count == 0) // Cases 1,2 : standard and blending elements
            {
                integrationPoints = standardIntegrationRule.IntegrationPoints;
            }
            else // Case 3: Enriched element
            {
                integrationPoints = enrichedIntegrationRule.GenerateIntegrationPoints(element);
            }

            pointsAndMaterials = new Dictionary<GaussPoint2D, IFiniteElementMaterial2D>();
            foreach (GaussPoint2D gaussPoint in integrationPoints)
            {
                pointsAndMaterials.Add(gaussPoint, commonMaterial.Clone());
            }
        }

        public void Update(XContinuumElement2D element)
        {
            throw new NotImplementedException();
        }
    }
}
