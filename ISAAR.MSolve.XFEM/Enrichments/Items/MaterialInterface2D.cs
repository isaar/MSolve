using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.XFEM.Enrichments.Items
{
    class MaterialInterface2D : AbstractEnrichmentItem2D
    {
        private readonly IFiniteElementMaterial2D materialInNegativeDomain;
        private readonly IFiniteElementMaterial2D materialInPositiveDomain;
        
        public ICurve2D Discontinuity { get; }

        public MaterialInterface2D(ICurve2D geometry, 
            IFiniteElementMaterial2D materialInNegativeDomain, IFiniteElementMaterial2D materialInPositiveDomain)
        {
            this.Discontinuity = geometry;
            this.EnrichmentFunctions = new IEnrichmentFunction2D[] { new RampFunction2D(this) };

            // I am cloning these to be safe. If I had immutable materials or material factories or material views I 
            // could avoid redundant copies in many classes that use materials
            this.materialInNegativeDomain = materialInNegativeDomain.Clone();
            this.materialInPositiveDomain = materialInPositiveDomain.Clone();
        }

        // TODO: add some tolerance when checking around 0. Perhaps all this is not needed though and I could even 
        // ignore points on the interface. It certainly needs a better name
        /// <summary>
        /// Finds the material at the requested cartesian point and returns false. If the point lies directly on the 
        /// interface, a null reference is set as the material (out param) and true is returned.
        /// </summary>
        /// <param name="point"></param>
        /// <param name="materialAtThisPoint"></param>
        /// <returns></returns>
        public bool LocatePoint(ICartesianPoint2D point, out IFiniteElementMaterial2D materialAtThisPoint)
        {
            int sign = Math.Sign(Discontinuity.SignedDistanceOf(point));
            if (sign < 0)
            {
                materialAtThisPoint = materialInNegativeDomain;
                return false;
            }
            else if (sign > 0)
            {
                materialAtThisPoint = materialInPositiveDomain;
                return false;
            }
            else
            {
                materialAtThisPoint = null;
                return false;
            }
        }

        public override IReadOnlyList<ICartesianPoint2D> IntersectionPointsForIntegration(XElement2D element)
        {
            return Discontinuity.IntersectionWith(element);
        }
    }
}
