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
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.XFEM.Integration.Strategies
{
    class HomogeneousIntegration2D: IIntegrationStrategy2D<XContinuumElement2D>
    {
        /// <summary>
        /// Might need it later when updating. TODO: If I had a MaterialView, I could store the same material object 
        /// for all elements and save up memory.
        /// </summary>
        private readonly IFiniteElementMaterial2D commonMaterial;
        private readonly IIntegrationRule2D<XContinuumElement2D> enrichedIntegrationRule;

        /// <summary>
        /// This may be updated throughout the analysis.
        /// </summary>
        private Dictionary<GaussPoint2D, IFiniteElementMaterial2D> pointsAndMaterials;

        public HomogeneousIntegration2D(IIntegrationRule2D<XContinuumElement2D> enrichedIntegrationRule, 
            IFiniteElementMaterial2D commonMaterial)
        {
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
            if (element.EnrichmentItems.Count == 0)
            {
                // Case 1: Blending element with at least one tip enriched node. TODO: abstract the handling of blending elements
                if (HasTipEnrichedNodes(element))
                {
                    // TODO: verify the need to use higher order integration for tip blending elements
                    if (element.StandardQuadrature == GaussLegendre2D.Order2x2)
                    {
                        integrationPoints = GaussLegendre2D.Order4x4.IntegrationPoints;
                    }
                    else
                    {
                        throw new NotImplementedException("TODO: this is hardcoded for Quad4. " +
                            "Have the standard quadrature give its next(and then that ones next) order.");
                    }
                }
                else // Case 2: Standard element or blending element without tip enriched nodes
                {
                    integrationPoints = element.StandardQuadrature.IntegrationPoints;
                }
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

        // TODO: determining the state of the element is the responsibility of XContinuumElement2D.
        private bool HasTipEnrichedNodes(XContinuumElement2D element)
        {
            foreach (XNode2D node in element.Nodes)
            {
                foreach (IEnrichmentItem2D enrichment in node.EnrichmentItems)
                {
                    if (enrichment is CrackTip2D) return true;
                }
            }
            return false;
        }
    }
}
