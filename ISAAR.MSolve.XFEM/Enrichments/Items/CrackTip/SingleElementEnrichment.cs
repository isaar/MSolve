using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.Shapes;

namespace ISAAR.MSolve.XFEM.Enrichments.Items.CrackTip
{
    /// <summary>
    /// TODO: Handle the case where the crack tip is located on an edge or node belonging to multiple elements.
    /// </summary>
    class SingleElementEnrichment: ITipEnrichmentAreaStrategy
    {
        private readonly double magnificationOfJintegralRadius;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="magnificationOfJintegralRadius">The outer countour of the J-integral domain is defined as:
        ///     radius = magnification * sqrt(areaOfElementContainingTip). This parameter is the magnification.</param>
        public SingleElementEnrichment(double magnificationOfJintegralRadius)
        {
            // TODO: Add checks for valid values of the magnitude.
            this.magnificationOfJintegralRadius = magnificationOfJintegralRadius;
        }

        public double ComputeRadiusOfJintegralOuterContour(CrackTip2D tipItem)
        {
            // TODO: what happens if the tip is on an element's edge/node?
            var outline = ConvexPolygon2D.CreateUnsafe(tipItem.TipElement.Nodes);
            return magnificationOfJintegralRadius * Math.Sqrt(outline.ComputeArea()); 
        }

        public IReadOnlyList<XNode2D> SelectNodesForEnrichment(CrackTip2D tipItem)
        {
            return tipItem.TipElement.Nodes;
        }
    }
}
