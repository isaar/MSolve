using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.XFEM.Interpolation;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Enrichments.Items
{
    // Connects the geometry, model and enrichment function entities.
    // TODO: At this point it does most of the work in 1 class. Appropriate decomposition is needed.
    public interface IEnrichmentItem2D
    {
        // Perhaps the nodal dof types should be decided by the element type (structural, continuum) in combination with the EnrichmentItem2D and drawn from XContinuumElement2D
        IReadOnlyList<EnrichedDof> Dofs { get; } 

        //IReadOnlyList<XContinuumElement2D> AffectedElements { get; }

        ///// <summary>
        ///// Assigns enrichment functions and their nodal values to each enriched node.
        ///// </summary>
        //void EnrichNodes();

        //void EnrichElement(XContinuumElement2D element);

        IReadOnlyList<CartesianPoint2D> IntersectionPointsForIntegration(XContinuumElement2D element);

        double[] EvaluateFunctionsAt(XNode2D node);

        EvaluatedFunction2D[] EvaluateAllAt(NaturalPoint2D point, XContinuumElement2D element,
             EvaluatedInterpolation2D interpolation);
    }
}
