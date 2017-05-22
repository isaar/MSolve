using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Functions;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Enrichments.Items
{
    // Connects the geometry, model and enrichment function entities.
    // TODO: At this point it does most of the work in 1 class. Appropriate decomposition is needed.
    interface IEnrichmentItem2D
    {
        // Perhaps the nodal dof types should be decided by the element type (structural, continuum) in combination with the EnrichmentItem2D and drawn from XContinuumElement2D
        IReadOnlyList<ArtificialDOFType> DOFs { get; } 

        IReadOnlyList<XContinuumElement2D> AffectedElements { get; }

        /// <summary>
        /// Assigns enrichment functions and their nodal values to each enriched node.
        /// </summary>
        void EnrichNodes();

        void EnrichNode(XNode2D node); // TODO: delete after debugging

        void EnrichElement(XContinuumElement2D element);

        IReadOnlyList<ICartesianPoint2D> IntersectionPointsForIntegration(XContinuumElement2D element);

        double[] EvaluateFunctionsAt(ICartesianPoint2D point);

        EvaluatedFunction2D[] EvaluateAllAt(INaturalPoint2D point, IReadOnlyList<XNode2D> elementNodes,
             EvaluatedInterpolation2D interpolation);
    }
}
