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
        ///     radius = magnification * sqrt(areaOfElementContainingTip). This parameter is the magnification. 
        ///     It should be at least 1.5 (see "Modeling quasi-static crack growth with the extended finite element 
        ///     method Part II: Numerical applications, Huang et al, 2003" page 7546). Usually values 2-3 are selected 
        ///     (see Ahmed thesis, 2009).</param>
        public SingleElementEnrichment(double magnificationOfJintegralRadius)
        {
            // TODO: Add checks for valid values of the magnitude.
            this.magnificationOfJintegralRadius = magnificationOfJintegralRadius;
        }

        public double ComputeRadiusOfJintegralOuterContour(CrackTip2D tipItem)
        {
            double maxTipElementArea = -1.0;
            foreach (var element in tipItem.TipElements)
            {
                var outline = ConvexPolygon2D.CreateUnsafe(element.Nodes);
                double elementArea = outline.ComputeArea();
                if (elementArea > maxTipElementArea) maxTipElementArea = elementArea;
            }
            return magnificationOfJintegralRadius * Math.Sqrt(maxTipElementArea); 
        }

        public IReadOnlyList<XNode2D> SelectNodesForEnrichment(CrackTip2D tipItem)
        {
            var nodes = new List<XNode2D>();
            foreach (var element in tipItem.TipElements) nodes.AddRange(element.Nodes);
            return nodes;
        }
    }
}
