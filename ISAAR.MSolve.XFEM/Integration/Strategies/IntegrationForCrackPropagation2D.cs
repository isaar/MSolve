using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Integration.Points;

namespace ISAAR.MSolve.XFEM.Integration.Strategies
{
    class IntegrationForCrackPropagation2D: IIntegrationStrategy2D<XContinuumElement2D>
    {
        // TODO: verify the need to use higher order integration for tip blending elements
        private readonly IIntegrationStrategy2D<XContinuumElement2D> integrationForTipBlendingElements;
        private readonly IIntegrationStrategy2D<XContinuumElement2D> integrationForEnrichedElements;
        
        public IntegrationForCrackPropagation2D(
            IIntegrationStrategy2D<XContinuumElement2D> integrationForEnrichedElements,
            IIntegrationStrategy2D<XContinuumElement2D> integrationForTipBlendingElements)
        {
            this.integrationForEnrichedElements = integrationForEnrichedElements;
            this.integrationForTipBlendingElements = integrationForTipBlendingElements;
        }

        public IReadOnlyList<GaussPoint2D> GenerateIntegrationPoints(XContinuumElement2D element)
        {
            if (element.EnrichmentItems.Count == 0)
            {
                // Case 1: Blending element with at least one tip enriched node. TODO: abstract the handling of blending elements
                if (HasTipEnrichedNodes(element))
                {
                    return integrationForTipBlendingElements.GenerateIntegrationPoints(element);
                    #region legacy hard-coded integration for Quad 4 blending elements
                    //if (element.StandardQuadrature == GaussLegendre2D.Order2x2)
                    //{
                    //    return GaussLegendre2D.Order4x4.IntegrationPoints;
                    //}
                    //else
                    //{
                    //    throw new NotImplementedException("TODO: this is hardcoded for Quad4. " +
                    //        "Have the standard quadrature give its next (and then that one's next) order.");
                    //}
                    #endregion
                }
                else // Case 2: Standard element or blending element without tip enriched nodes
                {
                    return element.StandardQuadrature.IntegrationPoints;
                }
            }
            else // Case 3: Enriched element
            {
                return integrationForEnrichedElements.GenerateIntegrationPoints(element);
            }
        }

        // TODO: determining the state of the element might be the responsibility of XContinuumElement2D.
        private bool HasTipEnrichedNodes(XContinuumElement2D element)
        {
            foreach (XNode2D node in element.Nodes)
            {
                foreach (IEnrichmentItem2D enrichment in node.EnrichmentItems.Keys)
                {
                    if (enrichment is CrackTipEnrichments2D) return true;
                }
            }
            return false;
        }
    }
}
